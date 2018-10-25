function F = root1(x,y0,eta,R) %Can do other way like Sarah: 3 equs 3 unknowns
nu = 4;
sigma = 3;
kappa = 4; 
beta = 0.99;

c = x(1);
hph = x(2);
hpl = x(3);
h = x(4);
a = x(5);
cph = x(6);
cpl = x(7);


% c - c'
F(1) = (1/( R*beta*( .5*(1/cph)^(sigma) + .5*(1/(cpl)^(sigma)) )))^(1/sigma) - c;  % for c<a always hold

% BC_1
F(2) = eta*h + y0 - c - a;

% BC_2
F(3) = (eta+.05)*hph + R*a - cph;
F(4) = (eta-.05)*hpl + R*a - cpl;

% F(3) = .5*(kappa*hp)^(1/nu)/(eta+0.05) + .5*(kappa*hp)^(1/nu)/(eta-0.05) - (.5*(1/(hp*(eta+0.05)+R*a))^sigma+ .5*(1/(hp*(eta-0.05)+R*a))^sigma);
% F(3) =  .5*((kappa*hp^(1/nu))/((1/((eta+0.05)*hp + R*a))^sigma * (eta+0.05))) + .5*((kappa*hp^(1/nu))/((1/((eta-0.05)*hp + R*a))^sigma * (eta-0.05))) -1;
%F(3) = .5*((kappa*hp^(1/nu)) -(eta+0.05)*(1/((eta+0.05)*hp + R*a))^(sigma)) + .5*((kappa*hp^(1/nu)) -(eta-0.05)*(1/((eta-0.05)*hp + R*a))^(sigma))             ; 
% h' - c'
F(5) = .5*(hph^(1/nu)-((eta+.05)/kappa)*(cph)^(-sigma))+.5*(hpl^(1/nu)-((eta - .05)/kappa)*(cpl)^(-sigma));

% F(4) = eta*c^(-sigma)-kappa*h^(1/nu);
% c-h
F(6) = c - (eta/(kappa*h^(1/nu)))^(1/sigma);

% h - h'
F(7) = h - ((eta/kappa)*R*beta*( .5*(kappa*hph^(1/nu)/(eta+.05)) + .5*(kappa*hpl^(1/nu)/(eta-.05)) ))^(nu);

% c - h'
%F(8) = c - (R*beta*( .5*(kappa*hph^(1/nu)/(eta+.05)) + .5*(kappa*hpl^(1/nu)/(eta-.05))))^(-1/sigma);
% c' - h
%F(9) = h - ((eta/kappa) * (R*beta*(.5*(cph)^(-sigma)+.5*(cpl)^(-sigma))))^(nu); 
end