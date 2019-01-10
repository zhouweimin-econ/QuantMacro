function [res,coeff,distr] = model_stst(bettaguess)
% Function to compute the steady state of the Krusell-Smith model
% Steps: for a given guess of aggregate capital K 
% 1) solve the individual problem 
% 2) calculate steady state distribution (s) of agents on (assets, idShock)
% 3) check if the market clearing condition K=sum(s.*assets) is satisfied

% Input: the guess of aggregate capital (K)
% Output: 
%       - resid: residual in the market clearing condition
%       - coeff: coefficients solving the individual problem
%       - s: stationary distribution
global n indx nodes Grid
global betta alf sig delta Zbar rhoz sigma_epsz TransMat_IdShock amin amax
global coeffguess distrguess

%% Step 1: Solve the individual problem (i.e. find coeff), for given aggregate variables (K)

% initialize variabeles
coeff = coeffguess; 
betta = exp(bettaguess);
K_stst = 1;
R = 1+alf*Zbar*K_stst^(alf-1)-delta;
W = (1-alf)*Zbar*K_stst^alf; 
AggVar=[R;W]; % aggregate capital and productuvity in steady state

% EulerResiduals_stst = @(x) EulerResiduals(x, x, AggVar, AggVar);
% opt_options = [];%optimset('Display','iter');
% coeff = fsolve(EulerResiduals_stst, coeffguess,opt_options);
% coeffguess= coeff;
% [~,~,Pmat] = EulerResiduals(coeff, coeff, AggVar, AggVar);

% options for iterative scheme
iter = 1;           % initialize counter;
maxiter = 10000;     % maximum number of iterations
tol = 1e-12;         % convergence criterion
smooth = 0.5;       % weight on new coeff when updating guess

while iter<maxiter  
    [resid,coeff_new,Pmat] = EulerResiduals(coeff, coeff, AggVar, AggVar);
    diff = max(abs(resid));
    
    if diff<tol
        coeffguess = coeff_new; %store coefficients (coeffguess is global: speed up code for next iteration on K) 
        break; 
    end
    
    %fprintf('Iter: %d, Max. diff.: %g \n',iter, diff);
    iter = iter + 1;
    coeff  = smooth*coeff_new + (1-smooth)*coeff; 
end

if iter == maxiter; error('Individual Problem not solved'); end

%% Step 2: Find steady state distribution of agents, using iteration
% ... alternative: solve system (I-P)*s=0, 
% ... i.e. find eigenvector of P associated with unit eigenvalues (very slow with large matrices)

crit = 1; % initialize stopping criterion (arbitrary large number)
distr= distrguess;
while crit>1e-12
    distr_new = Pmat'*distr;
    distr_new = distr_new/sum(distr_new);
    crit=norm(distr_new - distr);
    distr = distr_new;
end

distrguess=distr;
%% Step 3: verify guess
Knew = sum(distr.*Grid(:,1));
res = (Knew - K_stst); % check if guess was correct

