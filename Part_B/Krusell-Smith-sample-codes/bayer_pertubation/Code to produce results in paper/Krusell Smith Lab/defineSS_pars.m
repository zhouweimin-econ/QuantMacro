%% Parameters
% Household Parameters
par.beta        = 0.99;     % Discount factor
par.xi          = 1;          % CRRA
par.gamma       = 1;          % Inverse Frisch elasticity

% Income Process
par.rhoH        = 0.95;    % Persistence of productivity
par.sigmaH      = 0.06;    % STD of productivity shocks

% Firm Side Parameters
par.alpha       = 0.64;%2/3;       % Labor share 2/3
par.delta       = 0.1/4;     % Depreciation rate
par.phi         = 0;        % Capital adj costs

% Price of Capital in SS
par.Q  = 1;

%% Grids
% Idiosyncratic States
mpar.nm         = 100;
mpar.nk         = 0;
mpar.nh         = 2;
mpar.tauchen    ='importance';


%% Numerical Parameters
mpar.crit    = 1e-10;
