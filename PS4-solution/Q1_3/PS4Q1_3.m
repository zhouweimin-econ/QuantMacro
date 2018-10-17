% =======================================================================
% Quant Macro PS-4  
% Weimin Zhou
% Due: 17, Oct, 2018
% This file is to use Chebyshev Approximation, and redo item 1.
% Notice that Chebyshev polynomials donot put any assumption 
% on the shape of the value function (i.e. concave and strictly increasing)
% =======================================================================
clear;clf;close all;
cd '~/Desktop/PS4/Q1_3'  % in order to save png
disp('Solving deterministic (no labor) by Chebyshev Approximation of value function')
%%
tic
delta   = 0.031;
beta    = 0.988;
alpha   = 0.321;
nbk     = 100;
p       = 10;
crit    = 1;
iter    = 1;
epsi    = 1e-4;
maxits  = 100;
ks      = ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));
dev     = 0.8;
kmin    = (1-dev)*ks;
kmax    = (1+dev)*ks;

% Interpolating nodes
rk      = -cos((2*[1:nbk]'-1)*pi/(2*nbk));  
kgrid   = kmin+(rk+1)*(kmax-kmin)/2;   %transform     

% Initial guess for the approximation
v = zeros(nbk,1);
% v       = (log(kgrid.^alpha))/(1-beta);
X       = chebyshev(rk,nbk);
th0     = X\v;
Tv      = zeros(nbk,1);
kp      = zeros(nbk,1);

% options = optimset('GradObj','on');
while (crit>epsi && iter < maxits)
    k0   = kgrid(1);
    for i=1:nbk
        param = [alpha beta delta kmin kmax nbk kgrid(i)];
        kp(i) = fminunc(@tv,k0,[],param,th0);
        k0    = kp(i);
        Tv(i) = -tv(kp(i),param,th0);
    end;
    
    theta= X\Tv;
    crit = norm(Tv-v); 
    v =Tv;
    th0 = theta;
    iter= iter+1;
end
disp('Value Function Iterations for the deterministic growth model are completed.')
toc
% takes over 117 seconds
plot(kp,v);
