function  F = focdynamic1(x,kt,invest)
% x =[ht  ct]';
alfa  = 0.321; 
kappa = 5.24; 
nu    = 2.0;

F = (1-alfa)*kt^alfa * x^(-alfa) - kappa*x^(1/nu)*(invest - kt^alfa * x^(1-alfa)); % f'h * u'c = - u'h 
end