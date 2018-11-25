% =========================================================================
% Part of the Matlab code to accompany the paper
% 'Solving heterogeneous agent models in discrete time with
% many idiosyncratic states by perturbation methods', CEPR Discussion Paper
% DP13071
% by Christian Bayer and Ralph Luetticke
% http://www.erc-hump-uni-bonn.net/
% =========================================================================

%% Initialize workspace and load directories
addpath(genpath('functions'))
addpath(genpath('latex'))

% Switch options
casename='SS_BASELINE_HANC_KS';
printIRFs     = false;

%% Select aggregate parameters
grid.K          = linspace(0.975*targets.K,1.025*targets.K,mpar.nK);
% Guess for Law of Motion
par.beta_K      = zeros(2,mpar.ns);
par.beta_K(2,:) = 1;

mpar.nblocks    = 10;
mpar.weight     = 0.25;
%% Use Tauchen method to approximate TFP process

% Define indeces for converting joint states to individual ones
h_index=repmat((1:mpar.nh),mpar.ns,1);
h_index=h_index(:); % maps joint states to h dimension
s_index=repmat((1:mpar.ns)',1,mpar.nh);
s_index=s_index(:); % maps joint states to s dimension

% Calculate transition probability matrix for joint process
P=NaN(mpar.ns*mpar.nh);
for hh=1:mpar.nh
    for s=1:mpar.ns
        P(s + (hh-1)*mpar.ns,:)=kron(P_H(hh,:),P_S(s,:));
    end
end

P=P./repmat(sum(P,2),[1 size(P,2)]);

hmesh=repmat(grid.h,[mpar.ns 1]);
hmesh=hmesh(:)';

smesh=repmat(grid.s',[1 mpar.nh]);
smesh=smesh(:)';

[meshes.m,meshes.h] = ndgrid(grid.m,hmesh,grid.K);
[TFP] = ndgrid(grid.s,grid.K);
[mesh.m,mesh.h] = ndgrid(grid.m,grid.h);

%%
distKS=9999;
countOL=0;
while distKS > mpar.crit
    countOL=countOL+1;
    
    %% Prepare items for EGM
    % Layout of matrices:
    %   Dimension 1: capital k
    %   Dimension 2: stochastic sigma s x stochastic human capital h
    %   Dimension 3: aggregate capital K
    
    % Construct relevant returns and interpolation matrices from KS LoMs    
    K_aux  = repmat(grid.K, [mpar.ns 1]);
    K_prime=NaN(mpar.ns,mpar.nK);
    for ss=1:mpar.ns %PI=f(s,M(-1),M,K) <=> PI'=f(s',M,M',K')
        K_prime(ss,:)     = exp(par.beta_K(1,ss)   + par.beta_K(2,ss).*log(K_aux(ss,:)));
    end
    
    % Implied price of capital
    Q=par.phi*(K_prime./K_aux-1) + 1;
    
    aux=permute(K_prime,[2,1]); %f(K,s)
    
    [~,indexK]=histc(aux,grid.K);
    
    % Replace the index of out of bounds forecasts (=1) by lower and upper bound for
    % the end of interval gridpoints, i.e. 2 for too small forecasts and
    % mpar.nM for too large forecasts.
    indexK(aux<=grid.K(1))=1;
    indexK(aux>=grid.K(end))=mpar.nK-1;
    weightsK = max(min(1,(aux-grid.K(indexK))./(grid.K(indexK+1)-grid.K(indexK))),0);
    
    P_K=zeros(mpar.nK, mpar.nK, mpar.ns);
    for ss=1:mpar.ns
        for KK=1:mpar.nK
            P_K(KK,indexK(KK,ss)+1,ss) = weightsK(KK,ss)   ;
            P_K(KK,indexK(KK,ss),ss) = (1-weightsK(KK,ss));
        end
    end
    
    
    idm=repmat((1:mpar.nm)',[1 mpar.ns*mpar.nh mpar.nK]);
    idm=idm(:);
    idsh=repmat((1:mpar.ns*mpar.nh)',[1 mpar.nm mpar.nK]);
    idsh=permute(idsh,[2 1 3]);
    idsh=idsh(:);
    
    idK = repmat(reshape(repmat(permute(indexK,[2 1]),...
        [mpar.nh 1]),[1 mpar.ns*mpar.nh*mpar.nK]),[mpar.nm 1]);
    idK=idK(:);
    
    W_K = repmat(reshape(repmat(permute(weightsK,[2 1]),...
        [mpar.nh 1]),[1 mpar.ns*mpar.nh*mpar.nK]),[mpar.nm 1]);
    W_K=W_K(:);
    
    index1=sub2ind([mpar.nm mpar.ns*mpar.nh mpar.nK],idm,idsh,idK);
    index2=sub2ind([mpar.nm mpar.ns*mpar.nh mpar.nK],idm,idsh,idK+1);
    
    weight1 = (1-W_K);
    weight2 = (W_K);
    NN=mpar.nm*mpar.ns*mpar.nh*mpar.nK;
    rowindex=repmat(1:NN,[1 2]);
    TT = sparse(rowindex,[index1(:); index2(:)],[weight1(:); weight2(:)],NN,NN);  
    
    % Returns
    N    =  (TFP.*par.alpha.*K_aux.^(1-par.alpha)).^(1/(1-par.alpha+par.gamma));
    W_fc =  TFP.*par.alpha.* (K_aux./N).^(1-par.alpha);
    R_fc = TFP.*(1-par.alpha).* (N./K_aux).^(par.alpha)- par.delta;
    Y = TFP.*N.^(par.alpha).*K_aux.^(1-par.alpha);
    Profits=  1/2*par.phi*((K_prime-K_aux).^2)./K_aux;
    WW=par.gamma/(1+par.gamma).*(N./par.H).*W_fc;
    
    WW=reshape(permute(repmat(WW,[1 1 mpar.nh mpar.nm])...
        ,[4 1 3 2]),[mpar.nm mpar.ns*mpar.nh mpar.nK ]); %Wages
    RR=reshape(permute(repmat(R_fc,[1 1 mpar.nh mpar.nm])...
        ,[4 1 3 2]),[mpar.nm mpar.ns*mpar.nh mpar.nK ]); %Dividend
    QQ=reshape(permute(repmat(Q,[1 1 mpar.nh mpar.nm])...
        ,[4 1 3 2]),[mpar.nm mpar.ns*mpar.nh mpar.nK ]); %Capital
    PP=reshape(permute(repmat(Profits,[1 1 mpar.nh mpar.nm])...
        ,[4 1 3 2]),[mpar.nm mpar.ns*mpar.nh mpar.nK ]); %Profits
    
    inc.labor   = WW.*meshes.h;
    inc.rent    = RR.*meshes.m;
    inc.capital = QQ.*meshes.m;
    inc.profits = PP;
    
    % Consumption guess 
    c_guess = inc.labor + inc.rent + inc.capital + inc.profits;
       
    %% Use EGM to update household policies until convergence
    distC=9999;
    
    disp('Solving household problem by EGM')
    
    while distC>mpar.crit
        mutil_c = 1./(c_guess.^par.xi); % marginal utility at consumption policy no adjustment
        
        mutil_c=(1+RR).*mutil_c;
        aux=reshape(permute(mutil_c,[2 1 3]),[mpar.ns*mpar.nh mpar.nm*mpar.nK]);
        % form expectations
        aux = permute(reshape(P*aux,[mpar.ns*mpar.nh mpar.nm mpar.nK]),[2 1 3]);
        EMU_aux = par.beta* TT*aux(:);
        EMU = reshape(EMU_aux,[mpar.nm mpar.ns*mpar.nh mpar.nK]);
        
        c_aux = 1./(EMU.^(1/par.xi));
        
        k_aux = (c_aux + inc.capital - inc.labor - inc.profits);
        k_aux = k_aux./(RR+1);
        
        % Identify binding constraints
        binding_constraints = meshes.m < repmat(k_aux(1,:,:),[mpar.nm 1 1]);
        
        % Consumption when drawing assets m' to zero: Eat all Resources
        Resource = inc.labor + inc.capital + inc.profits;
        
        k_aux = reshape(k_aux,[mpar.nm mpar.ns*mpar.nh*mpar.nK]);
        c_aux= reshape(c_aux,[mpar.nm mpar.ns*mpar.nh*mpar.nK]);
        
        %Interpolate grid.m and c_n_aux defined on k_aux over grid.m
        c_update=zeros(mpar.nm,mpar.ns*mpar.nh*mpar.nK);
        k_update=zeros(mpar.nm,mpar.ns*mpar.nh*mpar.nK);

        for hh=1:mpar.ns*mpar.nh*mpar.nK
            Savings=griddedInterpolant(k_aux(:,hh),grid.m); % generate savings function a(s,a*)=a'
            k_update(:,hh)=Savings(grid.m); % Obtain m'(m,h) by Interpolation
            Consumption=griddedInterpolant(k_aux(:,hh),c_aux(:,hh)); % generate consumption function c(s,a*(s,a'))
            c_update(:,hh)=Consumption(grid.m);  % Obtain c(m,h) by interpolation (notice this is out of grid, used linear interpolation)
        end
        
        c_update = reshape(c_update,[mpar.nm, mpar.ns*mpar.nh, mpar.nK]);
        k_update = reshape(k_update,[mpar.nm, mpar.ns*mpar.nh, mpar.nK]);
        
        c_update(binding_constraints) = Resource(binding_constraints)-grid.m(1);
        k_update(binding_constraints) = min(grid.m);
        
        distC = max((abs(c_guess(:)-c_update(:))));
        
        c_guess=c_update;
        
    end
    
    c_star=c_guess;
    k_star=k_update;
    k_star=reshape(k_star,[mpar.nm mpar.ns mpar.nh mpar.nK]);
    
    
    %% Simulation Krusell Smith
    %   Initialize time series data
    disp('KS simulation')
    joint_distr=joint_distr(:)';
    
    pr_s=rand(RandStream('mcg16807', 'Seed',20180621),mpar.nblocks,mpar.T(1,1));
    PS_S=cumsum(P_S,2);
    
    K_sim  = zeros(mpar.T(1,1),mpar.nblocks);
    ss_sim     = zeros(mpar.T(1,1),mpar.nblocks);
        
    parfor bb = 1:mpar.nblocks
                   
        [K_sim(:,bb), ss_sim(:,bb)]=KSsimulation(k_star,P_H,PS_S,pr_s(bb,:),joint_distr,grid,mesh,mpar);
        
    end
        
    %% Update coefficients of LOM
    
    % Regression inputs
    ssaux=ss_sim(mpar.T(2,1):mpar.T(1,1)-2,:);
    ssaux=ssaux(:);
    
    K=log((K_sim(mpar.T(2,1):mpar.T(1,1)-2,:)));
    K=K(:);
    K_prime=log((K_sim(mpar.T(2,1)+1:mpar.T(1,1)-1,:)));
    K_prime=K_prime(:);
    
    ss_dummyaux=zeros(size(ss_sim(mpar.T(2,1)+1:mpar.T(1,1)-1,:)));
    ss_dummies=NaN(mpar.ns,size(ss_dummyaux(:),1));
    ss_aux=ss_sim(mpar.T(2,1)+1:mpar.T(1,1)-1,:);
    for j=1:mpar.ns
        ss_dummyaux(ss_aux==j)=1;
        
        ss_dummies(j,:)=ss_dummyaux(:);
        
        ss_dummyaux=zeros(size(ss_sim(mpar.T(2,1)+1:mpar.T(1,1)-1,:)));
    end
    
    % Regressors
    X=[ss_dummies(2:mpar.ns,:)' repmat(K,[1 mpar.ns]).*ss_dummies'];
  
    % Regressions
    statsK = regstats(K_prime,X,'linear',{'beta' 'covb' 'rsquare' 'mse','r','tstat'});
    % Statistics
    R2_K  = statsK.rsquare;
    
    clear statsbetaK_aux
    statsK.beta(2:mpar.ns)=statsK.beta(2:mpar.ns)+statsK.beta(1);
    statsbetaK_aux(1,:)=statsK.beta(1:mpar.ns)';
    statsbetaK_aux(2,:)=statsK.beta(mpar.ns+1:2*mpar.ns)';    
    
    % Check convergence of LoM coefficients
    distKS = max(max(abs(statsbetaK_aux-par.beta_K)))
    
    par.beta_K   = (1-mpar.weight)*real(statsbetaK_aux)+(mpar.weight)*par.beta_K;
    
    
end
% toc
