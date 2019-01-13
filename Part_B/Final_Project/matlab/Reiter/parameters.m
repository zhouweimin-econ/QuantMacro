%% KS parameters setup: Reiter version
% require matlab 2016+ 

global n indx nodes Grid
global betta alf sig delta Zbar rhoz sigma_epsz TransMat_IdShock amin amax

n.I      = 100;  % points for asset distribution
n.S      = 2;    % idiosyncratic shock
n.Aggr   = 6;    % aggregate states 
n.Shocks = 1;    % aggregate shocks 

n.ExpErr = n.I*n.S; % number of equations with expectations 

%betta = .99;   % Discount Factor
alf     = 0.3;      % Capital Income Share
sig     = 1;        % Risk Aversion
delta   = 0.025;  % Depreciation rate

% AggrShock = (1-rhoz)*Zbar + rhoz*AggrShock_lag + epsz
Zbar    = (1/.99 - (1-delta))/alf; 
rhoz    = 0.9;    
sigma_epsz=1;  % std.dev 

% idShock = rho_idShock*idShock_lag + eps_idShock
sigma_idShock = 0.1; % std.dev
rho_idShock  = 0.9;   

% AMPLITUDE OF THE GRID ON ASSETS
amin  = 0;            % lower bound on individual assets
amax  = 5;            % upper bound on individual assets 


%% Labels
vlab    = ...      
char('Consumption ',...
'Investment',...
'Capital',...
'Interest Rate',...
'Output',...
'Productivity Shock');

%% BUILD THE VECTOR with the POSITION INDEX of all elements
n.Agents  = n.S*n.I;
n.Var     = 2*n.Agents+n.Aggr; % number of variables
indx.current  = 1:n.Var;
indx.lag      = n.Var+1:2*n.Var;
indx.iidShock = max(indx.lag)+1:max(indx.lag)+n.Shocks;
indx.ExpErr   = max(indx.iidShock)+1:max(indx.iidShock)+n.ExpErr;
indx.total    = max(indx.ExpErr);

indx.AggVar = [(indx.current(2*n.Agents+1):indx.current(2*n.Agents+n.Aggr))';(indx.lag(2*n.Agents+1):indx.lag(2*n.Agents+n.Aggr))'];

%% BUILD GRID for individual state variables 
% assets
nodes.assets    = linspace(amin,amax,n.I)';

% idiosyncratic shock
[TransMat_IdShock, log_idShock] = tauchen(rho_idShock,0, (sigma_idShock)^2, n.S,1); 
nodes.idShock  = exp(log_idShock);

% calculate stationary distribution of idiosyncratic shocks
crit=1;
distr_idshock_old = ones(n.S,1)/n.S;
while abs(crit)>1e-8
    distr_idshock=TransMat_IdShock*distr_idshock_old;
    crit = norm(distr_idshock-distr_idshock_old);
    distr_idshock_old=distr_idshock;
end

% make sure the mean equals to one
idshock_mean  = distr_idshock'*nodes.idShock;
nodes.idShock = nodes.idShock/idshock_mean;

% build combined Grid (first column: assets, second column: with (nI*nS) rows)  
Grid  = [kron(ones(n.S,1),nodes.assets) kron(nodes.idShock,ones(n.I,1))];           

