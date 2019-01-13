%{
 Final Project of Part B Quantitative Macroeconomics, UAB
 Weimin Zhou
 Jan, 2019
 codes modification based on Github Resources (mentioned in my readme.md):
 Descriptions:
  For Step 2. approximte model dynamics by perturbation
%}
function resid = model_equations(x)

global n indx Grid
global betta alf sig delta Zbar rhoz sigma_epsz TransMat_IdShock amin amax

x(indx.AggVar) =exp(x(indx.AggVar));          % log of endog. variables


coeff = x(indx.current(1:n.Agents));               % coefficients for individual choices (log(c_now))
distr = x(indx.current(n.Agents+1:2*n.Agents));    % distribution of agents
c     = x(indx.current(2*n.Agents+1));              
inv   = x(indx.current(2*n.Agents+2));              
k     = x(indx.current(2*n.Agents+3));              
R     = x(indx.current(2*n.Agents+4));               
y     = x(indx.current(2*n.Agents+5));               
z     = x(indx.current(2*n.Agents+6));              

coeff_lag = x(indx.lag(1:n.Agents));           
distr_lag = x(indx.lag(n.Agents+1:2*n.Agents));     
c_lag   = x(indx.lag(2*n.Agents+1));               
inv_lag = x(indx.lag(2*n.Agents+2));               
k_lag   = x(indx.lag(2*n.Agents+3));               
R_lag   = x(indx.lag(2*n.Agents+4));               
y_lag   = x(indx.lag(2*n.Agents+5));               
z_lag   = x(indx.lag(2*n.Agents+6));               

epsz = x(indx.iidShock);
eta_EE  = x(indx.ExpErr);   % for euler equation
%% recover residuals and transition matrix from HH
W = (1-alf)*y; % no endogenous labor
W_lag = (1-alf)*y_lag;
AggVar = [R,W];
AggVar_lag = [R_lag, W_lag];
[EulerResid,~,Pmat] = EulerResiduals(coeff_lag, coeff, AggVar_lag, AggVar);   
%% equations dynamics
resid(1:n.Agents) = EulerResid + eta_EE;                          % Consumption Euler Equation
resid(n.Agents+1:2*n.Agents) = - distr + Pmat'*distr_lag;         % law of motion
resid(2*n.Agents+1) = - c + sum(distr.*exp(coeff));               % Consumption
resid(2*n.Agents+2) = - inv + k - (1-delta)*k_lag;                % Investment
resid(2*n.Agents+3) = - k_lag + sum(distr.*Grid(:,1));            % Market Clearing of Capital
resid(2*n.Agents+4) = - R + (1+alf*y/k_lag - delta);              % Interest Rate from Firm
resid(2*n.Agents+5) = - y + Zbar*z*k_lag^alf;                     % Production Function from Firm
resid(2*n.Agents+6) = - log(z) + rhoz*log(z_lag) + epsz;          % DGP z 

resid = resid(:);
