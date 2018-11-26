function [Difference,LHS,RHS,JD_new,c_a_star,m_a_star,k_a_star,c_n_star,m_n_star,P]  = Fsys(State,Stateminus,...
    Control_sparse,Controlminus_sparse,StateSS,...
    ControlSS,Gamma_state,indexMUdct,indexVKdct,DC1,DC2,DC3,...
    par,mpar,grid,targets,Copula,P,aggrshock)
% System of equations written in Schmitt-Grohé-Uribe generic form with states and controls
% STATE: Vector of state variables t+1 (only marginal distributions for histogram)
% STATEMINUS: Vector of state variables t (only marginal distributions for histogram)
% CONTROL: Vector of state variables t+1 (only coefficients of sparse polynomial)
% CONTROLMINUS: Vector of state variables t (only coefficients of sparse polynomial)
% STATESS and CONTROLSS: Value of the state and control variables in steady
% state. For the Value functions these are at full grids.
% GAMMA_STATE: Mapping such that perturbationof marginals are still
% distributions (sum to 1).
% PAR, MPAR: Model and numerical parameters (structure)
% GRID: Liquid, illiquid and productivity grid
% TARGETS: Stores targets for government policy
% COPULA: Interpolant that allows to map marginals back to full-grid
% distribuitions
% P: steady state transition matrix
% aggrshock: sets whether the Aggregate shock is TFP or uncertainty
%
% =========================================================================
% Part of the Matlab code to accompany the paper
% 'Solving heterogeneous agent models in discrete time with
% many idiosyncratic states by perturbation methods', CEPR Discussion Paper
% DP13071
% by Christian Bayer and Ralph Luetticke
% http://www.erc-hump-uni-bonn.net/
% =========================================================================


%% Initializations
util  = @(c)(c.^(1-par.xi))./(1-par.xi);
mutil = @(c)(1./(c.^par.xi));
invutil = @(u)(((1-par.xi).*u).^(1/(1-par.xi)));
invmutil = @(mu)((1./mu).^(1/par.xi));

% Number of states, controls
nx   = mpar.numstates; % Number of states
ny   = mpar.numcontrols; % number of Controls
NxNx = nx-mpar.os; % Number of states without aggregate
Ny = length(indexMUdct)+length(indexVKdct);
NN   = mpar.nm*mpar.nk*mpar.nh; % Number of points in the full grid

% Initialize LHS and RHS
LHS  = zeros(nx+Ny+mpar.oc,1);
RHS  = zeros(nx+Ny+mpar.oc,1);

%% Indexes for LHS/RHS
% Indexes for Controls
mutil_cind = 1:length(indexMUdct);
Vkind = length(indexMUdct)+(1:length(indexVKdct));

Qind  = Ny+1;
PIind = Ny+2;
Yind  = Ny+3;
Gind  = Ny+4;
Wind  = Ny+5;
Rind  = Ny+6;
Profitind  = Ny+7;
Nind  = Ny+8;
Tind  = Ny+9;
Kind  = Ny+10;
Bind  = Ny+11;

% Indexes for States
marginal_mind = (1:mpar.nm-1);
marginal_kind = (mpar.nm-1+(1:mpar.nk-1));
marginal_hind = (mpar.nm+mpar.nk-2 + (1:(mpar.nh-2)));

RBind = NxNx+1;
Sind  = NxNx+2;

%% Control Variables (Change Value functions according to sparse polynomial)
Control      = Control_sparse;
Controlminus = Controlminus_sparse;

Control(end-mpar.oc+1:end)       = ControlSS(end-mpar.oc+1:end) + Control_sparse(end-mpar.oc+1:end,:);
Controlminus(end-mpar.oc+1:end)  = ControlSS(end-mpar.oc+1:end) + (Controlminus_sparse(end-mpar.oc+1:end,:));

%% State Variables
% read out marginal histogramm in t+1, t
Distribution      = StateSS(1:end-2) + Gamma_state * State(1:NxNx);
Distributionminus = StateSS(1:end-2) + Gamma_state * Stateminus(1:NxNx);

% Aggregate Endogenous States
RB      = StateSS(end-1) + (State(end-1));
RBminus = StateSS(end-1) + (Stateminus(end-1));

% Aggregate Exogenous States
S       = StateSS(end) + (State(end));
Sminus  = StateSS(end) + (Stateminus(end));

%% Split the Control vector into items with names
% Controls
XX             = zeros(NN,1);
XX(indexMUdct) = Control(mutil_cind);
aux = reshape(XX,[mpar.nm, mpar.nk, mpar.nh]);
aux = myidct(aux,1,DC1); % do dct-transformation
aux = myidct(aux,2,DC2); % do dct-transformation
mutil_c_dev = myidct(aux,3,DC3); % do dct-transformation
mutil_c     = (mutil(mutil_c_dev(:)+ControlSS((1:NN))));

XX             = zeros(NN,1);
XX(indexVKdct) = Control(Vkind);
aux = reshape(XX,[mpar.nm, mpar.nk, mpar.nh]);
aux = myidct(aux,1,DC1); % do dct-transformation
aux = myidct(aux,2,DC2); % do dct-transformation
Vk_dev = myidct(aux,3,DC3); % do dct-transformation
Vk     = (mutil(Vk_dev(:)+ControlSS((1:NN)+NN)));


% Aggregate Controls (t+1)
PI = exp(Control(PIind));
Y  = exp(Control(Yind ));
K  = exp(Control(Kind ));
B  = exp(Control(Bind ));

% Aggregate Controls (t)
PIminus = exp(Controlminus(PIind));
Qminus  = exp(Controlminus(Qind ));
Yminus  = exp(Controlminus(Yind ));
Gminus  = exp(Controlminus(Gind ));
Wminus  = exp(Controlminus(Wind ));
Rminus  = exp(Controlminus(Rind ));
Profitminus  = exp(Controlminus(Profitind ));
Nminus  = exp(Controlminus(Nind ));
Tminus  = exp(Controlminus(Tind ));
Kminus  = exp(Controlminus(Kind ));
Bminus  = exp(Controlminus(Bind ));

%% Write LHS values
% Controls
LHS(nx+Vkind)      = Controlminus(Vkind);
LHS(nx+mutil_cind) = Controlminus(mutil_cind);
LHS(nx+Qind)       = (Qminus);
LHS(nx+Yind)       = (Yminus);

LHS(nx+Gind)       = (Gminus);
LHS(nx+Wind)       = (Wminus);
LHS(nx+Rind)       = (Rminus);
LHS(nx+Profitind)  = (Profitminus);
LHS(nx+Nind)       = (Nminus);
LHS(nx+Tind)       = (Tminus);
LHS(nx+Kind)       = (Kminus);
LHS(nx+Bind)       = (Bminus);

% States
% Marginal Distributions (Marginal histograms)
LHS(marginal_mind) = Distribution(1:mpar.nm-1);
ba=mpar.nm;
LHS(marginal_kind) = Distribution(ba+(1:mpar.nk-1));
ba=mpar.nm+mpar.nk;
LHS(marginal_hind) = Distribution(ba+(1:mpar.nh-2));
LHS(RBind)         = (RB);
LHS(Sind)          = (S);

% take into account that RB is in logs
RB=exp(RB); RBminus=exp(RBminus);


%% Set of Differences for exogenous process
RHS(Sind) = (par.rhoS * (Sminus));

switch(aggrshock)
    case('MP')
        EPS_TAYLOR=Sminus;
        TFP=1;
    case('TFP')
        TFP=exp(Sminus);
        EPS_TAYLOR=0;
    case('Uncertainty')
        TFP=1;
        EPS_TAYLOR=0;
        % Tauchen style for Probability distribution next period
        [P,~,~] = ExTransitions(exp(Sminus),grid,mpar,par);
end

% Calculate aggregate Capital, Bonds and Human Capital Supply in t
marginal_mminus = Distributionminus(1:mpar.nm)';
marginal_kminus = Distributionminus((1:mpar.nk)+mpar.nm)';
marginal_hminus = Distributionminus(mpar.nm+mpar.nk+(1:mpar.nh))';

Hminus  = sum(grid.h(1:end-1).*marginal_hminus(1:end-1)); %Last column is entrepreneurs.
Lminus  = sum(grid.m.*marginal_mminus);

RHS(nx+Bind) = Lminus ;
RHS(nx+Kind) = (sum(grid.k.*marginal_kminus));

% Calculate joint distributions
cumdist = zeros(mpar.nm+1,mpar.nk+1,mpar.nh+1);
cumdist(2:end,2:end,2:end) = Copula({cumsum(marginal_mminus),cumsum(marginal_kminus)',cumsum(marginal_hminus)});
JDminus = diff(diff(diff(cumdist,1,1),1,2),1,3);

% Generate meshes for b,k,h
[meshes.m,meshes.k,meshes.h] = ndgrid(grid.m,grid.k,grid.h);

%% Aggregate Output
mc              =  par.mu- (par.beta * log(PI)*Y/Yminus - log(PIminus))/par.kappa;

RHS(nx+Nind)    =  (par.tau*TFP*par.alpha*Kminus.^(1-par.alpha).*mc).^(1/(1-par.alpha+par.gamma));

RHS(nx+Yind)    = (TFP*(Nminus).^(par.alpha).*Kminus.^(1-par.alpha));

%% Prices that are not part of Control Vector

% Wage Rate
RHS(nx+Wind) =  TFP *par.alpha       *mc.* (Kminus./(Nminus)).^(1-par.alpha);
% Return on Capital
RHS(nx+Rind) =  TFP *(1-par.alpha)   *mc.* ((Nminus)./Kminus).^(par.alpha)  - par.delta;

% Profits for Entrepreneurs
RHS(nx+Profitind) = (1-mc)*Yminus - Yminus.*(1/(1-par.mu))./par.kappa./2 .*log(PIminus).^2 + 1/2*par.phi*((K-Kminus).^2)./Kminus;

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
RBaux = RB/PI + (meshes.m<0).*(par.borrwedge/PI);
EVm = reshape(reshape(RBaux(:).*mutil_c,[mpar.nm*mpar.nk mpar.nh])*P',[mpar.nm,mpar.nk mpar.nh]);

[c_a_star,m_a_star,k_a_star,c_n_star,m_n_star] = EGM_policyupdate(EVm,EVk,Qminus,PIminus,RBminus,inc,meshes,grid,par,mpar);

meshaux=meshes;
meshaux.h(:,:,end)=1000;

%% Update Marginal Value Bonds
mutil_c_n = mutil(c_n_star); % marginal utility at consumption policy no adjustment
mutil_c_a = mutil(c_a_star); % marginal utility at consumption policy adjustment
mutil_c_aux    = par.nu.*mutil_c_a + (1-par.nu).*mutil_c_n; % Expected marginal utility at consumption policy (w &w/o adjustment)
aux=((invmutil(mutil_c_aux(:)))-ControlSS((1:NN)));
aux = reshape(aux,[mpar.nm, mpar.nk, mpar.nh]);
aux = mydct(aux,1,DC1); % do dct-transformation
aux = mydct(aux,2,DC2); % do dct-transformation
aux = mydct(aux,3,DC3); % do dct-transformation
DC=aux(:);

RHS(nx+mutil_cind) = (DC(indexMUdct)); % Write Marginal Utility to RHS of F

%% Update marginal Value of Capital
EVk     = reshape(Vk,[mpar.nm*mpar.nk, mpar.nh])*P';
Vk_next = griddedInterpolant(meshaux.m,meshaux.k,meshaux.h,reshape(EVk,[mpar.nm mpar.nk mpar.nh]));
Vk_aux  = par.nu.*(Rminus+Qminus).*mutil_c_a + (1-par.nu).*Rminus.*mutil_c_n+ par.beta.*(1-par.nu).*Vk_next(m_n_star,meshaux.k,meshaux.h); % Expected marginal utility at consumption policy (w &w/o adjustment)

aux=((invmutil(Vk_aux(:)))-ControlSS((1:NN)+NN));
aux = reshape(aux,[mpar.nm, mpar.nk, mpar.nh]);
aux = mydct(aux,1,DC1); % do dct-transformation
aux = mydct(aux,2,DC2); % do dct-transformation
aux = mydct(aux,3,DC3); % do dct-transformation
DC=aux(:);
RHS(nx+Vkind) = (DC(indexVKdct));  % Write Marginal Value of K to RHS of F

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

% Next period marginal histograms
% liquid assets
aux_m = squeeze(sum(sum(JD_new,2),3));
RHS(marginal_mind) = aux_m(1:end-1); %Leave out last state
% illiquid assets
aux_k =  squeeze(sum(sum(JD_new,1),3))' ;
RHS(marginal_kind) = aux_k(1:end-1); %Leave out last state
% human capital
aux_h = squeeze(sum(sum(JD_new,1),2));
RHS(marginal_hind) = aux_h(1:end-2); %Leave out last state & entrepreneurs

%% Third Set: Government Budget constraint
% Return on bonds (Taylor Rule)
RHS(RBind) = log(par.RB) + par.rho_R* log(RBminus/par.RB)  + log(PIminus/par.PI).*((1-par.rho_R)*par.theta_pi)+EPS_TAYLOR;

taxrevenue =(1-par.tau).*Wminus.*Nminus +(1-par.tau).*Profitminus;

% Inflation jumps to equilibrate real bond supply and demand
RHS(nx+PIind) = par.rho_B * log((Bminus)/(targets.B)) ...
    + par.rho_B * log(RBminus/par.RB)...
    - (par.rho_B+par.gamma_pi) * log(PIminus/par.PI) ...
    - par.gamma_T * log((Tminus)/(targets.T));
LHS(nx+PIind) = log((B)/(targets.B));

% Government expenditures
RHS(nx+Gind) =  B - Bminus*RBminus/PIminus + Tminus;

RHS(nx+Tind) = taxrevenue;

% Resulting Price of Capital
RHS(nx+Qind)=(par.phi*(K./Kminus-1)+1)-par.ABS;

%% Difference
Difference=((LHS-RHS));

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
