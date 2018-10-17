% ====
% Weimin Zhou
% Value Function Approximation
function v = value(k,param,theta)
kmin    = param(1);
kmax    = param(2);
n       = param(3);
k       = 2*(k-kmin)/(kmax-kmin)-1;
v       = chebyshev(k,n)*theta;
end