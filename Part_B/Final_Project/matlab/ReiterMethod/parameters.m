%% KS parameters setup: Reiter version
% require matlab 2016+ 

global n indx nodes Grid
global betta alf sig delta Zbar rhoz sigma_epsz TransMat_IdShock amin amax

n.I      = 100;  % points for asset distribution
n.S      = 2;    % idiosyncratic shock
n.Aggr   = 6;    % aggregate states (c, inv, k, y, z, R)
n.Shocks = 1;    % aggregate shocks 

n.ExpErr = n.I*n.S; % number of equations with expectations 

%betta = .99;   % Discount Factor
alf = 0.3;sig = 1;delta = 0.025;

% AggrShock = (1-rhoz)*Zbar + rhoz*AggrShock_lag + epsz
Zbar = (1/.99 - (1-delta))/alf; 
rhoz = 0.9;sigma_epsz = 1;  % std.dev 

% idShock = rho_idShock*idShock_lag + eps_idShock
rho_idShock  = 0.9;sigma_idShock = 0.1; % std.dev
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
%% build grids for individual states

% assets: evenly spaced
amin = 0;amax = 5;
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

idshock_mean  = distr_idshock'*nodes.idShock;
nodes.idShock = nodes.idShock/idshock_mean;
% kron build combined Grid (assets ; shocks with (nI*nS) rows)  
Grid  = [kron(ones(n.S,1),nodes.assets) kron(nodes.idShock,ones(n.I,1))];           