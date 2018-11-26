function [Difference,LHS,RHS,JD_new,c_star,m_star,P]  = Fsys(State,Stateminus,...
    Control_sparse,Controlminus_sparse,StateSS,...
    ControlSS,Gamma_state,indexMUdct,Copula,...
    par,mpar,grid,meshes,P,aggrshock,oc,os)
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
% GRID: Asset and productivity grid
% COPULA: Interpolant that allows to map marginals back to full-grid
% distribuitions
% P: steady state productivity transition matrix
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
mutil = @(c)(1./(c.^par.xi));
invmutil = @(mu)((1./mu).^(1/par.xi));

% Number of states, controls
nx   = mpar.numstates; % Number of states
ny   = mpar.numcontrols; % number of Controls
NxNx = nx-os; % Number of states without aggregate
Ny = length(indexMUdct);
NN   = mpar.nm*mpar.nh; % Number of points in the full grid

% Initialize LHS and RHS
LHS  = zeros(nx+Ny+oc,1);
RHS  = zeros(nx+Ny+oc,1);

%% Indexes for LHS/RHS
% Indexes for Controls
mutil_cind = 1:Ny;
Qind  = Ny+1;
Yind  = Ny+2;
Wind  = Ny+3;
Rind  = Ny+4;
Nind  = Ny+5;
Kind  = Ny+6;

% Indexes for States
distr_ind = 1:mpar.nm*mpar.nh-1;

Sind  = NxNx+1;

%% Control Variables (Change Value functions according to sparse polynomial)
Control      = Control_sparse;
Controlminus = Controlminus_sparse;

Control(end-oc+1:end)       = ControlSS(end-oc+1:end) + Control_sparse(end-oc+1:end,:);
Controlminus(end-oc+1:end)  = ControlSS(end-oc+1:end) + (Controlminus_sparse(end-oc+1:end,:));


%% State Variables
% read out marginal histogramm in t+1, t
Distribution      = StateSS(1:end-os) + Gamma_state * State(1:NxNx);
Distributionminus = StateSS(1:end-os) + Gamma_state * Stateminus(1:NxNx);

% Aggregate Exogenous States
S       = StateSS(end) + (State(end));
Sminus  = StateSS(end) + (Stateminus(end));

%% Split the Control vector into items with names
% Controls
XX             = zeros(NN,1);
XX(indexMUdct) = Control(mutil_cind);
mutil_c        = (reshape(XX,[mpar.nm,mpar.nh]));
mutil_c        = (mutil(mutil_c(:)+ControlSS((1:NN))));

% Aggregate Controls (t+1)
K  = exp(Control(Kind ));
Q  = exp(Control(Qind ));
R  = exp(Control(Rind ));

% Aggregate Controls (t)
Qminus  = exp(Controlminus(Qind ));
Yminus  = exp(Controlminus(Yind ));
Wminus  = exp(Controlminus(Wind ));
Rminus  = exp(Controlminus(Rind ));
Nminus  = exp(Controlminus(Nind ));
Kminus  = exp(Controlminus(Kind ));

%% Write LHS values
% Controls
LHS(nx+mutil_cind) = Controlminus(mutil_cind);
LHS(nx+Qind)       = (Qminus);
LHS(nx+Yind)       = (Yminus);
LHS(nx+Wind)       = (Wminus);
LHS(nx+Rind)       = (Rminus);
LHS(nx+Nind)       = (Nminus);
LHS(nx+Kind)       = (Kminus);

% States
% Joint Distribution
LHS(distr_ind) = Distribution(1:mpar.nm*mpar.nh-1);

LHS(Sind)          = (S);

%% Set of Differences for exogenous process
RHS(Sind) = (par.rhoS * (Sminus));

switch(aggrshock)
    case('TFP')
        TFP=exp(Sminus);
    case('Uncertainty')
        TFP=1;
        % Tauchen style for Probability distribution next period
        [P,~,~] = ExTransitions(exp(Sminus),grid,mpar,par);
    
end

marginal_mminus = sum(reshape(Distributionminus,[mpar.nm mpar.nh]),2);
marginal_hminus = squeeze(sum(reshape(Distributionminus,[mpar.nm mpar.nh]),1));

Hminus  = grid.h(:)'*marginal_hminus(:); %Last column is entrepreneurs.

RHS(nx+Kind)= grid.m(:)'*marginal_mminus(:);

%% Update controls
RHS(nx+Nind)    =  (TFP*par.alpha*Kminus.^(1-par.alpha)).^(1/(1-par.alpha+par.gamma));

RHS(nx+Yind)    = (TFP*(Nminus).^(par.alpha).*Kminus.^(1-par.alpha));

% Wage Rate
RHS(nx+Wind) = TFP *par.alpha.* (Kminus./(Nminus)).^(1-par.alpha);
% Return on Capital
RHS(nx+Rind) = TFP *(1-par.alpha).* ((Nminus)./Kminus).^(par.alpha)  - par.delta;

RHS(nx+Qind)=(par.phi*(K./Kminus-1)+1);

%% Wages net of leisure services
WW=par.gamma/(1+par.gamma)*(Nminus./Hminus)*Wminus*ones(mpar.nm,mpar.nh);

%% Incomes (grids)
inc.labor   = WW.*(meshes.h);
inc.money   = (Rminus+Qminus)*meshes.m;
inc.profits  = 1/2*par.phi*((K-Kminus).^2)./Kminus;

%% Update policies
Raux = (R+Q)/(Qminus); % Expected marginal utility at consumption policy (w &w/o adjustment)
EVm = reshape(reshape(Raux(:).*mutil_c,[mpar.nm mpar.nh])*P',[mpar.nm, mpar.nh]);

[c_star,m_star] = EGM_policyupdate(EVm,Rminus,Qminus,inc,meshes,grid,par,mpar);

%% Update Marginal Value Bonds
mutil_c_aux = mutil(c_star); % marginal utility at consumption policy no adjustment
aux=((invmutil(mutil_c_aux(:)))-ControlSS((1:NN)));
DC=(reshape(aux,[mpar.nm,mpar.nh]));

DC=DC(:);

RHS(nx+mutil_cind) = (DC(indexMUdct)); % Write Marginal Utility to RHS of F
%% Update distribution
% find next smallest on-grid value for money choices
weight11  = zeros(mpar.nm, mpar.nh,mpar.nh);
weight12  = zeros(mpar.nm, mpar.nh,mpar.nh);

% Adjustment case
[Dist_m,idm] = genweight(m_star,grid.m);

idm=repmat(idm(:),[1 mpar.nh]);
idh=kron(1:mpar.nh,ones(1,mpar.nm*mpar.nh));

index11 = sub2ind([mpar.nm mpar.nh],idm(:),idh(:));
index12 = sub2ind([mpar.nm mpar.nh],idm(:)+1,idh(:));

for hh=1:mpar.nh
    
    %Corresponding weights
    weight11_aux = (1-Dist_m(:,hh));
    weight12_aux =  (Dist_m(:,hh));
    
    % Dimensions (mxk,h',h)
    weight11(:,:,hh)=weight11_aux(:)*P(hh,:);
    weight12(:,:,hh)=weight12_aux(:)*P(hh,:);
end

weight11=permute(weight11,[1 3 2]);
weight12=permute(weight12,[1 3 2]);

rowindex=repmat(1:mpar.nm*mpar.nh,[1 2*mpar.nh]);

H=sparse(rowindex,[index11(:); index12(:)],...
    [weight11(:); weight12(:)],mpar.nm*mpar.nh,mpar.nm*mpar.nh); % mu'(h',k'), a without interest


JD_new=Distributionminus(:)'*H;

RHS(distr_ind) = JD_new(1:mpar.nm*mpar.nh-1); %Leave out last state

JD_new = reshape(JD_new(:),[mpar.nm,mpar.nh]);


%% Difference from SS
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
