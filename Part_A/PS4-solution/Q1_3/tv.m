function res=tv(kp,param,theta)

alpha   = param(1);
beta    = param(2);
delta   = param(3);
kmin    = param(4);
kmax    = param(5);
n       = param(6);
k       = param(7);
kp      = sqrt(kp.^2); % restrict to be positive
v       = value(kp,[kmin kmax n],theta);
c       = k.^alpha+(1-delta)*k-kp;
d       = find(c<=0);
c(d)    = NaN;
util    = log(c);
util(d) = -1000;
res     = -(util+beta*v);
% insures positivity of k?
% computes the value function
% computes consumption
% find negative consumption
% get rid off negative c
% computes utility
% utility = low number for c<0
% compute -TV (we are minimizing)
end
                                