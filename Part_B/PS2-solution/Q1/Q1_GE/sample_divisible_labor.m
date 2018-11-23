%__________________________________________________________________________
%
% Quant Macro Part-2
% PS2-Q1: Solution of the model without aggregate risk
% Author: Weimin Zhou
% Date: Nov, 2018
%__________________________________________________________________________
clear;clf;clc;close all;
cd '~/Desktop/PS2'  % in order to save png
%% 0. preamble set-up
% The system of equations 
A= [ 1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0 ; ...
     0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1 ; ...
     0  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  1  0  0  0  1  0  0 ; ...
     0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0 ; ...
     1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0 5.6 0 -1  0  0  0  0  0 ; ...
    -1 0 28/3 0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
  .02 .48 .05 .45 0 0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0 0 .02 .48 .05 .45 0 0 0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0 ];

b= [7/8; 7/8; 7/8; 7/8; 1/8; 1/8; 1/8; 1/8; 7/24; 21/40; 0; 0; 0.02; 0.005; 0.05; 0.02];

pize = reshape(A^-1*b,4,4); 
% fix z=z_g, normalize the transition of good times:
pie_g0 = pize(1:2,1:2);
pie_b0 = pize(3:4,3:4);   
for i=1:2
    for j=1:2
        pie_g(i,j)=pie_g0(i,j)/sum(pie_g0(i,:));
    end
end
for i=1:2
    for j=1:2
        pie_b(i,j)=pie_b0(i,j)/sum(pie_b0(i,:));
    end
end

% transtion matrix aggregate state
piZ = [ 7/8  1/8;...
        1/8  7/8];
% or equilvalently:
% pie_g2 =pie_g/piZ(1,1); pie_b2 = pie_b/piZ(2,2);

% above for Question.1-3 part(a): compare the policy function for varying
% the transition matrix of the employment state e.
%% 1. parameters and functional forms

% parameters
alpha  = 0.36;                   % capital income share
beta   = 95/100;                % subjective discount factor
delta  = 2.5/1000;                 % depreciation rate of physical capital
gamma  = 3/2;                   % inverse elasticity of intertemporal substitution
varphi = 2/3;                   % Frisch elasticity of labor supply
z  = [0.7, 1.3]';
N      = 2;                     % number of possible productivity realizations
y1     = 0;                     % low productivity realization
y2     = 1;                     % high productivity realization
Gamma  = 0.35;

% transition probability matrix (productivity)
pi      = zeros(N,N);

question =1; %for good time transition

% probabilities for 1.1.3 - (a)
if question==1 % i
    pi = pie_g;
else           % ii
    pi = pie_b;
end

% (inverse) marginal utility functions
up    = @(c) c.^(-gamma);        % marginal utility of consumption
invup = @(x) x.^(-1/gamma);      % inverse of marginal utility of consumption
vp    = @(h) Gamma*h.^(1/varphi);      % marginal disutility of labor supply
invvp = @(x) (x/Gamma).^(varphi);        % inverse marginal disutility of labor supply    
%% 2. discretization

% set up asset grid
b  = 0;                          % borrowing constraint
M  = 250;                        % number of asset grid points
aM = 45;                         % maximum asset level
A  = linspace(b,aM,M)';          % equally-spaced asset grid from a_1=b to a_M

% set up productivity grid (idiosyncratic risk)
Y  = [y1,y2]';                   % grid for employment state

% vectorize the grid in two dimensions
Amat = repmat(A,1,N);            % values of A change vertically
Ymat = repmat(Y',M,1);           % values of Y change horizontally

% % this is the built-in alternative
% [Amat,Ymat] = ndgrid(A,Y);
%% 3. endogenous functions

% optimal labor supply
H  = @(c,y,w) invvp(up(c).*y.*w);

% current consumption level, cp0(anext,ynext) is the guess
C0 = @(cp0,r) invup(beta*(1+r)*up(cp0)*pi');
                
% current asset level, c0 = C0(cp0(anext,ynext))
A0 = @(anext,y,c0,r,w) 1/(1+r)...
                    *(c0+anext-H(c0,y,w).*y.*w);
%% 4. solve for the stationary equilibrium

% convergence criterion for consumption iteration
crit = 10^(-6);

% parameters of the simulation
% note: 
% it takes quite some simulation periods to get to the stationary
% distribution. Choose a high T >= 10^(-4) once the algorithm is running.
I = 10^(4);             % number of individuals
T = 10^(4);             % number of periods

%% General Equilibrium set-up 
% choose interval where to search for the stationary interest rate
% note: 
% the staionary distribution is very sensitive to the interst rate. 
% make use of the theoretical result that the stationary rate is slightly 
% below 1/beta-1

r0  = (1/beta-1)-[10^(-12),10^(-4)];
r0  = (1/beta-1)-10^(-4);
w0 = (1-alpha)*(alpha./(r0+delta)).^(alpha/(1-alpha));
% set up an anonymous function
fprintf('Start solving the Aiyagari model... \n');
tic;
myfun   = @(r) stationary_equilibrium(r,crit,I,T,Amat,Ymat,alpha,b,delta,pi,varphi,A0,C0,H);
options = optimset('display','iter','TolX',1e-8,'MaxIter',20);
rstar   = fzero(myfun,r0,options);
options1 = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
rstar   = fminunc(myfun, r0, options1);
fprintf('Done with the Aiyagari model in %f sec. \n',toc);

[residual1,at] = stationary_equilibrium(r,crit,I,T,Amat,Ymat,alpha,b,delta,pi,varphi,A0,C0,H);
% get the simulated asset levels
fprintf('Fetching the wealth distribution... \n');
[r,at,cp0] = stationary_equilibrium(rstar,crit,I,T,Amat,Ymat,alpha,b,delta,rho,varphi,A0,C0,H);
%% 5. plot the wealth distribution

% use the last 100 periods
[n,xout] = hist(at(:,T-100:T),M);     % choose M bins
bar(xout,n/sum(n));                   % relative frequency is n/sum(n)
title(sprintf('Wealth distribution, rho = %2.1f',rho));
xlabel('Asset level')
ylabel('Relative frequency')
%%
[M,N] = size(Amat);

% get productivity realizations from the first row
y1 = Ymat(1,1);
y2 = Ymat(1,2);
%%
r0  = (1/beta-1)-10^(-4);
%%
question =1; %for good time transition

% probabilities for 1.1.3 - (a)
if question==1 % i
    pi = pie_g;
else           % ii
    pi = pie_b;
end

zed = 1.8;
%zed = 0.6;
residual = 1
while abs(residual)>0.01
    
    
    % compute the wage according to marginal pricing and a Cobb-Douglas 
% production function
    w0 = zed*(1-alpha)*(alpha./(r0+delta)).^(alpha/(1-alpha));

% initial guess on future consumption (consume asset income plus labor
% income from working h=1.
    cp0 = r0*Amat+Ymat*w0;
    cp0(cp0==0)=0.01;
    cp0(cp0==NaN)=0.01;
    cp0(cp0==Inf)=0.01;
% distance between two successive functions start with a distance above the
% convergence criterion
    dist    = crit+1;
% maximum iterations (avoids infinite loops)
    maxiter = 10^(3);
% counter-
    iter    = 1;

    fprintf('Inner loop, running... \n');
  
    while (dist>crit) && (iter<maxiter)

% derive current consumption
        c0 = C0(cp0,r0);
%c0(c0==0)=0.01;
%c0(c0==NaN)=0.01;
%c0(c0==Inf)=0.01;
% derive current assets
        a0 = A0(Amat,Ymat,c0,r0,w0);
%a0 = real(a0);
%a0(a0==0)=0.01;
%a0(a0==NaN)=0.01;
%a0(a0==Inf)=0.01;
%%% update the guess for the consumption decision

% consumption decision rule for a binding borrowing constraint can be
% solved as a quadratic equation c^(2)-((1+r)a+b)c-(yw)^(1+2/3) = 0 general
% formula and notation: ax^(2)+bx+c = 0 x = (-b+sqrt(b^2-4ac))/2a
        cpbind = ((1+r0)*Amat+b)/2+sqrt(((1+r0)*Amat+b).^(2)+4*(Ymat*w0).^(1+varphi))/2;
%cpbind(cpbind==0)=0.01;
%cpbind(cpbind==NaN)=0.01;
%cpbind(cpbind==Inf)=0.01;
% % slow alternative: rootfinding gives the same result cpbind  =
% zeros(M,N); options = optimset('Display','off'); cpbind(:,1) =
% fsolve(@(c)
% (1+r0)*Amat(:,1)+H(c,Y(1),w0).*Y(1).*w0-c+b,cp0(:,1),options);
% cpbind(:,2) = fsolve(@(c)
% (1+r0)*Amat(:,2)+H(c,Y(2),w0).*Y(2).*w0-c+b,cp0(:,2),options);

% consumption for nonbinding borrowing constraint
        cpnon = zeros(M,N);
% interpolation conditional on productivity realization instead of
% extrapolation use the highest value in a0


        cpnon(:,1)  = interp1(a0(:,1),c0(:,1),Amat(:,1),'spline');
        cpnon(:,2)  = interp1(a0(:,2),c0(:,2),Amat(:,2),'spline');
%cpnon(cpnon==0) = 0.01;
%cpnon(cpnon==NaN)=0.01;
%cpnon(cpnon==Inf)=0.01;
% merge the two, separate grid points that induce a binding borrowing
% constraint for the future asset level (the first observation of the
% endogenous current asset grid is the threshold where the borrowing
% constraint starts binding, for all lower values it will also bind as the
% future asset holdings are monotonically increasing in the current asset
% holdings).
        cpnext(:,1) = (Amat(:,1)>a0(1,1)).*cpnon(:,1)+(Amat(:,1)<=a0(1,1)).*cpbind(:,1);
        cpnext(:,2) = (Amat(:,2)>a0(1,2)).*cpnon(:,2)+(Amat(:,2)<=a0(1,2)).*cpbind(:,2);
%cpnext(cpnext==0) = 0.01;
%cpnext(cpnext==NaN) = 0.01;
%cpnext(cpnext==Inf) = 0.01;
% distance measure
        dist = norm((cpnext-cp0)./cp0);

% display every 100th iteration
%if mod(iter,100) == 1
%    fprintf('Inner loop, iteration: %3i, Norm: %2.6f \n',[iter,dist]);
%end

% increase counter
        iter = iter+1;

% update the guess on the consumption function
%cp0 = real(cpnext);

    end

    fprintf('Inner loop, done. \n');

    fprintf('Starting simulation... \n');

%%% simulate the stationary wealth distribution

% initialize variables
    at      = zeros(I,T+1);
    yt      = zeros(I,T);
    ct      = zeros(I,T);
    ht      = zeros(I,T);

    at(:,1) = 1;            % initial asset level

    for t=1:T
   % draw uniform random numbers across individuals
        s = unifrnd(0,1,I,1);
        if t>1   % persistence for t>1
            yt(:,t) = ((s<=pi(1,1))&(yt(:,t-1)==y1)).*y1+((s>pi(2,1))&(yt(:,t-1)==y2)).*y1...
                +((s<=pi(2,2))&(yt(:,t-1)==y2)).*y2+((s>pi(1,2))&(yt(:,t-1)==y1)).*y2;
        else     % random allocation in t=0
            yt(:,t) = (s<=1/2).*y1+(s>1/2).*y2;
    end

   % consumption
        ct(:,t)   = (yt(:,t)==y1).*interp1(Amat(:,1),cp0(:,1),at(:,t),'spline')+(yt(:,t)==y2).*interp1(Amat(:,2),cp0(:,2),at(:,t),'spline');

   % labor supply
        ht(:,t)   = H(ct(:,t),yt(:,t),w0);

   % future assets
        at(:,t+1) = (1+r0)*at(:,t)+ht(:,t).*yt(:,t).*w0-ct(:,t); 

    end

    fprintf('simulation done... \n');

% compute aggregates from market clearing
    K = mean(mean(at(:,T-100:T)));
    L = mean(mean(yt(:,T-100:T).*ht(:,T-100:T)));
    r = alpha*(K/L)^(alpha-1)-delta;

% compute the distance between the two
    residual = (r-r0)/r0
    r0 = r

end

% for z_H = 1.6:
% r1 = 0.0525
% r2 = 0.0466, residual = -0.1132
% r3 = 0.0436, residual = -0.0650 done..
% ...
% rn = 0.0409, residual = -0.0098
% for z_L = 0.8：
% r1 = 0.0525
% r2 = 0.0693, residual = -0.0620
% r3 = 0.0469, residual = -0.0490
% ...
% rn = 0.0404, residual = -0.0085