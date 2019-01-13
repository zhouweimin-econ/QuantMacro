%{
 Final Project of Part B Quantitative Macroeconomics, UAB
 Weimin Zhou
 Jan, 2019
 codes modification based on Github Resources (mentioned in my readme.md):
 Descriptions:
  compute residual of Market Clearing condition
%}

function [resid,coeff_new,Pmat] = EulerResiduals(coeff_now, coeff_next, AggVar_now, AggVar_next)
global n indx nodes Grid
global betta alf sig delta Zbar rhoz sigma_epsz TransMat_IdShock amin amax

R_now = AggVar_now(1);
W_now = AggVar_now(2);
R_next = AggVar_next(1);
W_next = AggVar_next(2);
%% current variables
asset_now = Grid(:,1);
shock_now = Grid(:,2);

c_now = exp(coeff_now);
asset_next  = W_now.*shock_now + R_now.*asset_now - c_now; % budget constraint
%% future variables
% obtain consumption for all possible realizations of future shocks

c_val = reshape(exp(coeff_next),n.I,n.S);

for iS=1:n.S
    c_next(:,iS)=interp1(nodes.assets,c_val(:,iS),asset_next,'linear','extrap');
    RHS(:,iS) = kron(TransMat_IdShock(:,iS),ones(n.I,1)).*(betta*c_next(:,iS).^(-sig).*R_next);     
end
ExpRHS = sum(RHS,2);

%% Euler residual (take into account borrowing limit)
resid = min((c_now.^-sig - ExpRHS)./c_now.^(-sig),asset_next-amin);
Pmat = BuildTrans(asset_next,TransMat_IdShock,nodes.assets,n.S);

% update coefficients (if we use time iteration)
coeff_new = (ExpRHS.^(-1/sig));
is_constrained = ((asset_next-amin)<(c_now.^-sig - ExpRHS));
coeff_new(is_constrained) = (W_now.*shock_now(is_constrained) + R_now.*asset_now(is_constrained) - amin);
coeff_new = log(coeff_new);
