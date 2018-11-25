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


%% Select options
% Search for steady state and name it
FindNewSS           = false;
casename            = 'SS_BASELINE_TWOASSET_REDUCTION';

% Options
mpar.overrideEigen  = false; % Warning appears, but critical Eigenvalue shifted
print_IRFs          = false; % Plot and print as PDF

%% Produce matrices to reduce state-space
mainskript_statereduc

%% Initialize System
% disp('Computing system for SS.');
F = @(a,b,c,d)Fsys(a,b,c,d,Xss,Yss,Gamma_state,indexMUdct,indexVKdct,DC1,DC2,DC3,par,mpar,grid,targets,Copula,P_H,aggrshock);

[Fss,LHS,RHS,Distr] = F(State,State_m,Contr,Contr_m);


%% Solve RE via Schmitt-Grohe-Uribe Form
[hx,gx,F1,F2,F3,F4,par] = SGU_solver(F,mpar,par,p);

