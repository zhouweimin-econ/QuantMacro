%% Parameters
% Household Parameters
par.beta        = 0.99;     % Discount factor
par.xi          = 1;          % CRRA
par.gamma       = 1;          % Inverse Frisch elasticity
par.nu          = 0.065;          % Prob. of trade given adj. decision

% Income Process
par.rhoH        = 0.98;    % Persistence of productivity
par.sigmaH      = 0.06;    % STD of productivity shocks
mpar.in         = 0.00075;  % Prob. to become entrepreneur
mpar.out        = 0.0625;   % Prob. to become worker again

% Firm Side Parameters
par.eta         = 20;
par.mu          = (par.eta-1)/par.eta;       % Markup 5%
par.alpha       = 0.64/par.mu;  % Labor share 2/3
par.delta       = 0.1/4;     % Depreciation rate
par.phi         = 0;        % Capital adj costs

% Phillips Curve
par.prob_priceadj = 1/1000; % average price duration of 4 quarters = 1/(1-par.prob_priceadj)
par.kappa         = (1-par.prob_priceadj)*(1-par.prob_priceadj*par.beta)/par.prob_priceadj;% Phillips-curve parameter (from Calvo prob.)

% Central Bank Policy
par.theta_pi    = 1.25; % Reaction to inflation
par.rho_R       = 0.8;  % Inertia

% Tax Schedule
par.tau         = 0.7;   % Proportional tax on labor and profit income 

% Debt rule
par.gamma_pi    = 0;   % Reaction to inflation
par.gamma_T     = 0; % Reaction to tax revenue
par.rho_B       = 0.86;  % Autocorrelation

%% Returns
par.PI  = 1.00^.25;     % Gross inflation
par.RB  = (par.PI*1.01)^0.25;    % Real return times inflation

par.ABS = 0;                    % Loan to value ratio max.
par.borrwedge = par.PI*(1.125^0.25-1); % Wedge on borrowing beyond secured borrowing

par.Q  = 1;

%% Grids
% Idiosyncratic States
mpar.nm         = 100;
mpar.nk         = 100;
mpar.nh         = 4;
mpar.tauchen    ='importance';


%% Numerical Parameters
mpar.crit    = 1e-10;

