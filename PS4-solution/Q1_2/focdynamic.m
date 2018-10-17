function  F = focdynamic(x,kt,invest)
% x =[ht  ct]';
beta  = .988;
delta = 0.013;
alfa  = 0.321; 
kappa = 5.24; 
nu    = 2.0;

F(1) = (1-alfa)*kt^alfa * x(1)^(-alfa) - kappa*x(1)^(1/nu)*x(2);   % f'h * u'c = - u'h the intertemporal euler equation
F(2) = kt^alfa*x(1)^(1-alfa) - invest - x(2);  %Budge Constraint
end