function resid = model_equations(x)
% KS setup

global n indx Grid
global betta alf sig delta Zbar rhoz sigma_epsz TransMat_IdShock amin amax

x(indx.AggVar) =exp(x(indx.AggVar));          % these values of x are intended to be the log of endog. variables

% VARIABLES
coeff = x(indx.current(1:n.Agents));               % coefficients for individual choices
distr   = x(indx.current(n.Agents+1:2*n.Agents));   % distribution of agents
c   = x(indx.current(2*n.Agents+1));               % Aggregate Consumption
inv = x(indx.current(2*n.Agents+2));               % Investment
k   = x(indx.current(2*n.Agents+3));               % Capital
R   = x(indx.current(2*n.Agents+4));               % Interest Rate
y   = x(indx.current(2*n.Agents+5));               % Output
z   = x(indx.current(2*n.Agents+6));               % Productivity Shock

coeff_lag = x(indx.lag(1:n.Agents));               % coefficients of individual choices
distr_lag = x(indx.lag(n.Agents+1:2*n.Agents));     % distribution of agents
c_lag   = x(indx.lag(2*n.Agents+1));               % Consumption
inv_lag = x(indx.lag(2*n.Agents+2));               % Investment
k_lag   = x(indx.lag(2*n.Agents+3));               % Capital
R_lag   = x(indx.lag(2*n.Agents+4));               % Interest Rate
y_lag   = x(indx.lag(2*n.Agents+5));               % Output
z_lag   = x(indx.lag(2*n.Agents+6));               % Productivity Shock

% EXOGENOUS SHOCKS
epsz = x(indx.iidShock);       % Innovation of Productivity Shock

% EXPECTATION ERRORS
eta_EE  = x(indx.ExpErr);   % for EULER EQUATION


%% recover residuals and transition matrix from individual problem
W = (1-alf)*y;
W_lag = (1-alf)*y_lag;
AggVar = [R,W];
AggVar_lag = [R_lag, W_lag];
[EulerResid,~,Pmat] = EulerResiduals(coeff_lag, coeff, AggVar_lag, AggVar);   

%AggVar = [k,z];
%AggVar_lag = [k_lag, z_lag];
%[EulerResid,~,asset_next,Pmat] = EulerResiduals(coeff_lag, coeff, AggVar_lag, AggVar);   

%% EQUATIONS OF THE MODEL
resid(1:n.Agents) = EulerResid + eta_EE;                         % Consumption Euler Equation
resid(n.Agents+1:2*n.Agents) = - distr + Pmat'*distr_lag;         % Law of Motion of Distribution
resid(2*n.Agents+1) = - c + sum(distr.*exp(coeff));               % Aggregate Consumption
resid(2*n.Agents+2) = - inv + k - (1-delta)*k_lag;                % Aggregate Investment
resid(2*n.Agents+3) = - k_lag + sum(distr.*Grid(:,1));            % Market Clearing of Capital
resid(2*n.Agents+4) = - R + (1+alf*y/k_lag - delta);              % Interest Rate
resid(2*n.Agents+5) = - y + Zbar*z*k_lag^alf;                     % Production Function
resid(2*n.Agents+6) = - log(z) + rhoz*log(z_lag) + epsz;          % Exogenous Process for z

resid = resid(:);
