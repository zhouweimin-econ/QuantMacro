function  F = focss(x)
% x =[kstar hstar  css]';
global alpha beta delta kappa nu

F(1) = x(1) - ((1/beta-(1-delta))/(alpha*x(2)^(1-alpha)))^(1/(alpha-1)); % 1 = beta*(f'k + 1-delta)
F(2) = (1-alpha)*x(1)^alpha *  x(2)^(-alpha) - kappa*x(2)^(1/nu)*x(3);   % f'h * u'c = - u'h 
F(3) = x(1)^alpha *x(2)^(1-alpha) - x(3) -  delta*x(1) ; %Budge Constraint
end
