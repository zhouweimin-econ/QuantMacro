function  F = focs(x,y0,w,R)
% x =[c a h hph hpl]';
beta = .99;
sigma = 3;
kappa = 4; 
nu = 4;

% c - c'
F(1) = x(1,1) - (R*beta*(.5*(x(4,1)*(w+.05) + R*x(2,1))^(-sigma)+.5*(x(5,1)*(w-.05) + R*x(2,1))^(-sigma)))^(-1/sigma);
% Budge Constraint 1
F(2) = x(1,1) + x(2,1) - w*x(3,1) - y0;
% h - c
F(3) = x(3,1) - ((w/kappa)*(x(1,1)^(-sigma)))^(nu);

% c'-h'
F(4) = .5*(x(4,1)^(1/nu)-((w+.05)/kappa)*(x(4,1)*(w+.05) + R*x(2,1))^(-sigma))+.5*(x(5,1)^(1/nu)-((w - .05)/kappa)*(x(5,1)*(w-.05) + R*x(2,1))^(-sigma));
% h - h'
% F(5) = x(3,1) - ((w/kappa)*R*beta*( .5*(kappa*x(4,1)^(1/nu)/(w+.05)) + .5*(kappa*x(5,1)^(1/nu)/(w-.05)) ))^(nu);

%c - h'
F(6) = x(1,1) - (R*beta*( .5*(kappa*x(4,1)^(1/nu)/(w+.05)) + .5*(kappa*x(5,1)^(1/nu)/(w-.05))))^(-1/sigma);
% c' - h
F(7) = x(3,1) - ((w/kappa) * (R*beta*(.5*(x(4,1)*(w+.05) + R*x(2,1))^(-sigma)+.5*(x(5,1)*(w-.05) + R*x(2,1))^(-sigma))))^(nu); 
end
