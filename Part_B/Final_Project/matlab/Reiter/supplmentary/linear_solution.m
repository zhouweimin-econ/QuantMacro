% Function to solve the model using a first-order perturbation 
% using the gensys algorithm (by C. Sims)
%
% Inputs:
%    - funcname: name of the function containing the model equations
%    - stst:     vector with steady state values of endogeneous variables
%                (must be in the same order as in the model equations file)
%    - indx:     structure containing indices of position of variables
%                 (current, lag, iid shocks, expectation errors)
%    - doAutoDiff: 1 = use Automatic Differentiation (fast and accurate, require AutoDiff toolbox) 
%                  0 = use Numerical Differentiation (slow and less accurate, require file jacob.m)
%    - varargin: possibly other inputs for the model equations
% 
% Outputs: the State-Space Form  x(t) = C + G1*x(t-1) + impact*eps(t); 
%    - C:       vector of constants
%    - G1:      response to endogenous variables 
%    - impact:  response to exogenous shocks (iid innovation)

function [Gx,CONSTANT,IMPACT]=linear_solution(funcname,stst,indx, doAutoDiff)

% 1) Check that equations are satisfied at steady state
stst2 = zeros(indx.total,1);  
stst2([indx.current indx.lag]) = repmat(stst,2,1);

check = max(abs(funcname(stst2)));

if check>1e-5
   fprintf(1,'Maximum residual at stst is %e;',check);
   warning('Possible error in steady state calcuation or in model equations file');
end

% 2) Linearize the model around the steady state and obtain the linear system
%    Gamma0 X{t} = Gamma1 X_{t-1} + PSI*eps{t} + PI*exp_err{t};  (see Sims, eq. 1)
if doAutoDiff==1   % using Automatic Differentation
    
    stst2_ad = amatinit(stst2);          % initialize automatic differentiation at stst2
    resid_ad = feval(funcname,stst2_ad); % calculate level and 1st order derivatives of model equations
    JJ = ajac(resid_ad,0);               % Jacobian of vector field  
    
else              % using Numerical Differentiation
    JJ = jacob(funcname,stst2,1e-8);
end

Gamma0  = -JJ(:,indx.current);
Gamma1  = JJ(:,indx.lag);
PSI = JJ(:,indx.iidShock);
PI  = JJ(:,indx.ExpErr);
c   = zeros(size(Gamma0,1),1);

% 3) solve the model using the Sims' approach (as done e.g. in DYNARE)
div=1;
[Gx,CONSTANT,IMPACT]=gensys(Gamma0,Gamma1,c,PSI,PI,div); %see eq. 44 in Sims' paper for the meaning.

