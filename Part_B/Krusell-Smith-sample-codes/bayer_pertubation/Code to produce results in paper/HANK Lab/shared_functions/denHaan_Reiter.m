function [JD,RB,EqQuant]  = denHaan_Reiter(JDminus,Sminus,RBminus,...
    Control_F,Control_GUESS,ControlSS,indexMUdct,indexVKdct,par,mpar,grid,targets,meshes,P,aggrshock)

%% Initializations
mutil = @(c)(1./(c.^par.xi));

% Number of states, controls
Ny = length(indexMUdct)+length(indexVKdct);
NN   = mpar.nm*mpar.nk*mpar.nh; % Number of points in the full grid

%% Indexes for LHS/RHS
% Indexes for Controls
mutil_cind = 1:length(indexMUdct);
Vkind = length(indexMUdct)+(1:length(indexVKdct));

Qind  = Ny+1;
PIind = Ny+2;
Yind  = Ny+3;
Kind  = Ny+10;


%% Control Variables (Change Value functions according to sparse polynomial)
Control      = Control_F;
Controlminus = Control_GUESS;

Control(end-mpar.oc+1:end)       = ControlSS(end-mpar.oc+1:end) + Control_F(end-mpar.oc+1:end,:);
Controlminus(end-mpar.oc+1:end)  = ControlSS(end-mpar.oc+1:end) + (Control_GUESS(end-mpar.oc+1:end,:));

% Controls
XX             = zeros(NN,1);
XX(indexMUdct) = Control(mutil_cind);
aux = reshape(XX,[mpar.nm, mpar.nk, mpar.nh]);
aux = myidct(aux,1); % do dct-transformation
aux = myidct(aux,2); % do dct-transformation
mutil_c_dev = myidct(aux,3); % do dct-transformation
mutil_c        = (mutil(mutil_c_dev(:)+ControlSS((1:NN))));

XX             = zeros(NN,1);
XX(indexVKdct) = Control(Vkind);
aux = reshape(XX,[mpar.nm, mpar.nk, mpar.nh]);
aux = myidct(aux,1); % do dct-transformation
aux = myidct(aux,2); % do dct-transformation
Vk_dev = myidct(aux,3); % do dct-transformation
Vk             = (mutil(Vk_dev(:)+ControlSS((1:NN)+NN)));

% Aggregate Controls (t+1)
PI = exp(Control(PIind));
Y  = exp(Control(Yind ));
K  = exp(Control(Kind ));

% Aggregate Controls (t)
PIminus = exp(Controlminus(PIind));
Qminus  = exp(Controlminus(Qind ));
Yminus  = exp(Controlminus(Yind ));

%% States
Hminus  = par.H;%meshes.h(:)'*JDminus(:); 
Kminus  = meshes.k(:)'*JDminus(:);
Bminus  = meshes.m(:)'*JDminus(:);
RBminus = exp(RBminus);
%% Initialize
DD = @(Qg)PriceD(Qg,Vk,PI,Y/Yminus,Bminus,Kminus,Hminus,RBminus,Sminus,mutil_c,mpar,par,grid,aggrshock,P,JDminus,targets,meshes);
QuantGuess=[K,PIminus];

options=optimset('Display','off','TolFun',1e-10,'TolX',1e-10);

[EqQuant]=fsolve(@(p)DD(p),QuantGuess,options);
[~,JD,RB]= DD(EqQuant);
fprintf('%2.2f %8.2f \n',(abs(EqQuant-QuantGuess)./QuantGuess)*100);

end

function [QuantDiff,JD_new,RB] = PriceD(QuantGuess,Vk,PI,YG,Bminus,Kminus,Hminus,RBminus,Sminus,mutil_c,mpar,par,grid,aggrshock,P,JDminus,targets,meshes)
K=QuantGuess(1);
PIminus=QuantGuess(2);

switch(aggrshock)
    case('TFP')
        TFP=exp(Sminus);
        MPshock=0;
    case('Uncertainty')
         % Tauchen style for Probability distribution next period
        [P,~,~] = ExTransitions(exp(Sminus),grid,mpar,par);
        TFP=1;%(Hminus/par.H).^(-par.alpha);
        MPshock=0;
    case('MP')
        TFP=1;
        MPshock=Sminus;
end

Qminus=((K/Kminus)-1)*par.phi+1;

RB = log(par.RB) + par.rho_R* log(RBminus/par.RB)  + log(PIminus/par.PI).*((1-par.rho_R)*par.theta_pi) + MPshock;

Tauminus  = par.tau;

%% Aggregate Output
mc        =  par.mu- (par.beta * log(PI)*YG - log(PIminus))/par.kappa;

Nminus    =  (Tauminus*TFP*par.alpha*Kminus.^(1-par.alpha).*mc).^(1/(1-par.alpha+par.gamma));

YminusNEW = (TFP*(Nminus).^(par.alpha).*Kminus.^(1-par.alpha));

%% Prices that are not part of Control Vector

% Wage Rate
Wminus =  TFP *par.alpha       *mc.* (Kminus./(Nminus)).^(1-par.alpha);
% Return on Capital
Rminus =  TFP *(1-par.alpha)   *mc.* ((Nminus)./Kminus).^(par.alpha)  - par.delta;

% Profits for Entrepreneurs
Profitminus = (1-mc)*YminusNEW - YminusNEW.*(1/(1-par.mu))./par.kappa./2 .*log(PIminus).^2 + 1/2*par.phi*((K-Kminus).^2)./Kminus;

% Wages net of leisure services
WW=par.gamma/(1+par.gamma)*(Nminus./Hminus)*Wminus*ones(mpar.nm,mpar.nk,mpar.nh);
WW(:,:,end)=Profitminus*par.profitshare;

%% Incomes (grids)
inc.labor   = par.tau*WW.*(meshes.h);
inc.rent    = meshes.k*Rminus;
inc.capital = meshes.k*Qminus;
inc.money   = (RBminus/PIminus).*meshes.m...
    + (meshes.m<0).*(par.borrwedge/PIminus).*meshes.m;

%% First Set: Value Functions, Marginal Values
%% Update policies

EVk = reshape(reshape(Vk,[mpar.nm*mpar.nk mpar.nh])*P',[mpar.nm,mpar.nk mpar.nh]);
RBaux = exp(RB)/PI + (meshes.m<0).*(par.borrwedge/PI);
EVm = reshape(reshape(RBaux(:).*mutil_c,[mpar.nm*mpar.nk mpar.nh])*P',[mpar.nm,mpar.nk mpar.nh]);

[c_a_star,m_a_star,k_a_star,c_n_star,m_n_star] = EGM_policyupdate(EVm,EVk,Qminus,PIminus,RBminus,inc,meshes,grid,par,mpar);

%% Differences for distributions
% Initialize matrices
weight11  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
weight12  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
weight21  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
weight22  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);

% find next smallest on-grid value for money and capital choices
[Dist_m_a,idm_a] = genweight(m_a_star,grid.m);
[Dist_m_n,idm_n] = genweight(m_n_star,grid.m);
[Dist_k,idk_a]   = genweight(k_a_star,grid.k);
idk_n = repmat(ones(mpar.nm,1)*(1:mpar.nk),[1 1 mpar.nh]); %This is the actual point on grid


% Transition matrix adjustment case
idm_a=repmat(idm_a(:),[1 mpar.nh]);
idk_a=repmat(idk_a(:),[1 mpar.nh]);
idh=kron(1:mpar.nh,ones(1,mpar.nm*mpar.nk*mpar.nh));

index11 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_a(:),idk_a(:),idh(:));
index12 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_a(:),idk_a(:)+1,idh(:));
index21 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_a(:)+1,idk_a(:),idh(:));
index22 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_a(:)+1,idk_a(:)+1,idh(:));


for hh=1:mpar.nh
    
    %Corresponding weights
    weight21_aux =  Dist_m_a(:,:,hh).*(1-Dist_k(:,:,hh));
    weight11_aux = (1-Dist_m_a(:,:,hh)).*(1-Dist_k(:,:,hh));
    weight22_aux =  (Dist_m_a(:,:,hh)).*(Dist_k(:,:,hh));
    weight12_aux =  (1-Dist_m_a(:,:,hh)).*(Dist_k(:,:,hh));
    
    % Dimensions (mxk,h',h)
    weight11(:,:,hh)=weight11_aux(:)*P(hh,:);
    weight12(:,:,hh)=weight12_aux(:)*P(hh,:);
    weight21(:,:,hh)=weight21_aux(:)*P(hh,:);
    weight22(:,:,hh)=weight22_aux(:)*P(hh,:);
end
%Dimensions (mxk,h,h')
weight11=permute(weight11,[1 3 2]);
weight22=permute(weight22,[1 3 2]);
weight12=permute(weight12,[1 3 2]);
weight21=permute(weight21,[1 3 2]);
rowindex=repmat(1:mpar.nm*mpar.nk*mpar.nh,[1 4*mpar.nh]);

H_a=sparse(rowindex(:),[index11(:); index21(:); index12(:); index22(:)],...
    [weight11(:); weight21(:); weight12(:); weight22(:)],mpar.nm*mpar.nk*mpar.nh,mpar.nm*mpar.nk*mpar.nh); % mu'(h',k'), a without interest

% Policy Transition Matrix for no-adjustment case
weight11  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
weight21  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
idm_n=repmat(idm_n(:),[1 mpar.nh]);
idk_n=repmat(idk_n(:),[1 mpar.nh]);
index11 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_n(:),idk_n(:),idh(:));
index21 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_n(:)+1,idk_n(:),idh(:));


for hh=1:mpar.nh
    
    %Corresponding weights
    weight21_aux =  Dist_m_n(:,:,hh);
    weight11_aux = (1-Dist_m_n(:,:,hh));
    
    weight21(:,:,hh)=weight21_aux(:)*P(hh,:);
    weight11(:,:,hh)=weight11_aux(:)*P(hh,:);
    
end
weight11=permute(weight11,[1 3 2]);
weight21=permute(weight21,[1 3 2]);

rowindex=repmat(1:mpar.nm*mpar.nk*mpar.nh,[1 2*mpar.nh]);

H_n=sparse(rowindex,[index11(:); index21(:)],...
    [weight11(:); weight21(:)],mpar.nm*mpar.nk*mpar.nh,mpar.nm*mpar.nk*mpar.nh); % mu'(h',k'), a without interest

% Joint transition matrix and transitions

H=par.nu*H_a + (1-par.nu)*H_n;

JD_new=JDminus(:)'*H;
JD_new = reshape(JD_new(:),[mpar.nm,mpar.nk,mpar.nh]);

%% Third Set: Government Budget constraint
% Return on bonds (Taylor Rule)
taxrevenue =(1-Tauminus).*Wminus.*Nminus +(1-Tauminus).*Profitminus;
      
B = grid.B* exp(par.rho_B * log((Bminus)/(targets.B)) ...
    + par.rho_B * log(RBminus/par.RB)...
    - (par.rho_B+par.gamma_pi) * log(PIminus/par.PI) ...
    - par.gamma_T * log((taxrevenue)/(targets.T)));

BNEW=meshes.m(:)'*JD_new(:);
KNEW=meshes.k(:)'*JD_new(:);

QuantDiff=[ KNEW-K,BNEW-B];

end
function [ weight,index ] = genweight( x,xgrid )
% function: GENWEIGHT generates weights and indexes used for linear interpolation
%
% X: Points at which function is to be interpolated.
% xgrid: grid points at which function is measured
% no extrapolation allowed
[~,index] = histc(x,xgrid);
index(x<=xgrid(1))=1;
index(x>=xgrid(end))=length(xgrid)-1;

weight = (x-xgrid(index))./(xgrid(index+1)-xgrid(index)); % weight xm of higher gridpoint
weight(weight<=0) = 1.e-16; % no extrapolation
weight(weight>=1) = 1-1.e-16; % no extrapolation

end  % function
