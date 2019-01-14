%{
 Final Project of Part B Quantitative Macroeconomics, UAB
 Weimin Zhou
 Jan, 2019
 Descriptions:

function required in Step 1: compute SS

given guess K 
 1) solve HH, get coefficients
 2) compute stationary distribution of agents on (assets, idShock)
 3) check if residual of Market Clearing condition
%}
function [res,coeff,distr] = model_stst(bettaguess)
global n indx nodes Grid
global betta alf sig delta Zbar rhoz sigma_epsz TransMat_IdShock amin amax
global coeffguess distrguess
%% 1) solve HH, get coefficients

% initialize 
coeff = coeffguess; 
betta = exp(bettaguess);
K_stst = 1;
R = 1+alf*Zbar*K_stst^(alf-1)-delta;
W = (1-alf)*Zbar*K_stst^alf; 
AggVar=[R;W]; 


iter = 1;maxiter = 1000; tol = 1e-10;     
smooth = 0.5;       % weight on new coeff when updating guess

while iter<maxiter 
    [resid,coeff_new,Pmat] = EulerResiduals(coeff, coeff, AggVar, AggVar);
    diff = max(abs(resid));
    
    if diff<tol
        coeffguess = coeff_new;
        break; 
    end
    
    fprintf('Iter: %d, Max. diff.: %g \n',iter, diff);
    iter = iter + 1;
    coeff  = smooth*coeff_new + (1-smooth)*coeff; 
end
%% 2) compute stationary distribution of agents on (assets, idShock)
crit = 1; distr= distrguess;

while crit>1e-10
    distr_new = Pmat'*distr;
    distr_new = distr_new/sum(distr_new);
    crit=norm(distr_new - distr);
    distr = distr_new;
end

distrguess=distr;
%% 3) check if residual of Market Clearing condition

Knew = sum(distr.*Grid(:,1)); % newasset = distribution .* assets
res = (Knew - K_stst);       

