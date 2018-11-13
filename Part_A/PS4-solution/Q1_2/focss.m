function  F = focss(x)
% x =[kstar hstar  css]';
beta  = .988;
delta = 0.013;
alfa  = 0.321; 
kappa = 5.24; 
nu    = 2.0;

F(1) = x(1) - ((1/beta-(1-delta))/(alfa*x(2)^(1-alfa)))^(1/(alfa-1)); % 1 = beta*(f'k + 1-delta)
F(2) = (1-alfa)*x(1)^alfa *  x(2)^(-alfa) - kappa*x(2)^(1/nu)*x(3);   % f'h * u'c = - u'h 
F(3) = x(1)^alfa *x(2)^(1-alfa) - x(3) -  delta*x(1) ; %Budge Constraint
end
