var y I k a c w R r;
varexo e;
parameters alpha beta delta rho sigma sigmae;

load param_nc;
set_param_value('alpha',alpha);
set_param_value('beta',beta);
set_param_value('sigma',sigma);
set_param_value('delta',delta);
set_param_value('rho',rho);
set_param_value('sigmae',sigmae);

model;

exp(c)^(-sigma) = beta*exp(c(+1))^(-sigma)*(alpha*exp(a(+1))*exp(k)^(alpha-1) + (1-delta));
exp(y) = exp(a)*exp(k(-1))^(alpha);
exp(k) = exp(a)*exp(k(-1))^(alpha) - exp(c) + (1-delta)*exp(k(-1));
a = rho*a(-1) + e;
exp(y) = exp(c) + exp(I);
exp(c)^(-sigma) = beta*exp(c(+1))^(-sigma)*(1+r);
exp(R) = alpha*exp(a)*exp(k(-1))^(alpha-1);
exp(w) = (1-alpha)*exp(a)*exp(k(-1))^(alpha);


end;

initval;
k = log(30);
y = log(3);
c = log(2.5);
I = log(0.5);
a=0;
r = (1/beta) - 1;
R = log((1/beta) - (1-delta));
w = log(1);
end;


shocks;
var e = sigmae^2;
end;

steady;

stoch_simul(order=1,irf=20);
