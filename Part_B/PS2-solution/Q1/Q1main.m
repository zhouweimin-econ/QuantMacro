%__________________________________________________________________________
%
% Quant Macro Part-2
% PS2-Q1: Solution of the model without aggregate risk
% Author: Weimin Zhou
% Date: Nov, 2018
%__________________________________________________________________________
clear;clc;clf;close all; cd '~/Desktop/PS2';  % in order to save png

disp('Solving Aiyagari model with endogenous labor choice by VFI by generating labor grid and consumption grid separately around their own SS values')

%% 1. Partial Equilibrium: guess r0 and w0 
disp('Paritial Equilibrium Part')
%% 1.0. parameters and preamble set-up
%{
% otherwise, if we know each value in transition matrix, calculate pize as
follows:

Zvec = [ 0.80 1.20 ]; nz=2;
PZ = [ 0.90 0.10
       0.10 0.90 ];

Avec =  [ 0.95 1.05 ]; na=2;
PA = [ 0.90 0.10
       0.10 0.90 ];

% Lump-together idiosyncratic and aggregate shock in a single P

[ AME, ZME ]=ndgrid(Avec,Zvec);
pize = kron(PA,PZ);
%}
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
%or equilvalently:
% pie_g2 =pie_g/piZ(1,1); pie_b2 = pie_b/piZ(2,2);

% above for Question.1-3 part(a): compare the policy function for varying
% the transition matrix of the employment state e.


alpha = .36;
beta = .95;
delta = .0025;
kappa = .975;
psi = 1.78;
sigma = 1.5;
nu = 2.0;

z =[0,1];
Xi = (1/beta - 1 + delta)/alpha;
pi    = zeros(2,2);


question =2; %for good time transition

% probabilities for 1.1.3 - (a)
if question==1 % i
    pi = pie_g;
else           % ii
    pi = pie_b;
end

r0  = (1/beta-1)-[10^(-12),10^(-4)];
w0 = (1-alpha)*(alpha./(r0+delta)).^(alpha/(1-alpha));

dim = 2;    % Number of values the exogenous state variable might take
hg  = 11;   % Number of grid points for the control variable h
kg  = 101; % Number of grid points for both the control var and state var k
khratio = ((1/beta - 1 + delta)/alpha)^(1/(alpha-1));
chratio = khratio^alpha - delta*khratio;

hstar = 1/(1+psi*chratio/((1-alpha)*khratio^alpha));
kstar = khratio*hstar;
cstar = chratio*hstar;

clear khratio chratio;

kgrid =linspace(.15*kstar ,1.85*kstar , kg)';
hgrid = linspace(.1*hstar, 1.5*hstar, hg)';

c = zeros(hg,kg,dim,kg);
u = zeros(hg,kg,dim,kg);
%% 1.1 VFI
%i is a counter for the control variable h
%j is a counter for the control variable k'
%m is a counter for the control variable z
%n is a counter for the control variable k
warning off;
for i=1:hg
 for j = 1:kg
     for m=1:dim
         for n=1:kg
           c(i,j,m,n)= (1+r0(m))*kgrid(n) +  w0(m)*hgrid(i)*z(m) + (1-delta)*kgrid(n) - kgrid(j);
           if c(i,j,m,n)<0
               c(i,j,m,n)=0;
           end
         u(i,j,m,n) = c(i,j,m,n)^(1-sigma)/(1-sigma) - psi*hgrid(i)^(1+1/nu)/(1+1/nu);
         end
     end
 end
end


tic
warning on;
clear c 
clear i j m n

v = zeros(kg,dim);
convcrit = 1e-4; diff = 1; iter = 0;
while diff>convcrit && iter<1000
    diff = 0;
    for m= 1:dim
        for n=1:kg
            objfn(:,:,m,n) = u(:,:,m,n) + beta *(pi(m,1)*(v(:,1)*ones(1,hg))' + pi(m,2)*(v(:,2)*ones(1,hg))');
            Tv(n,m)=max(max(objfn(:,:,m,n)));
        end
    end
    diff = norm(v-Tv);
    v = Tv;
    iter = iter + 1;
    if (mod(iter,10)==0 || iter ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, diff); 
    end
end

for m=1:dim
    for n = 1:kg
        objfn(:,:,m,n)= u(:,:,m,n) + beta *(pi(m,1)*(v(:,1)*ones(1,hg))' + pi(m,2)*(v(:,2)*ones(1,hg))');
        [tmp1,x1] = max(objfn(:,:,m,n),[],1);
        [tmp2,x2] = max(tmp1,[],2);
        kgridrule(m,n)=x2;
        hgridrule(m,n)=x1(x2);
        
        kdecrule(m,n) = kgrid(kgridrule(m,n));
        hdecrule(m,n) = hgrid(hgridrule(m,n));
        
        cdecrule(m,n) = (1+r0(m))*kgrid(n) +  w0(m)*hdecrule(m,n)*z(m) + (1-delta)*kgrid(n) - kdecrule(m,n);
    end
end

clear tmp1 tmp2 x1 x2
clear diff convcrit m n 
clear objfn Tv
clear u

disp('Value Function Iterations for the growth model are completed.')
toc
disp('.')
disp('Number of iterations:'); 
disp(iter-1)

% plot part: 
%{ 
%% 1.2a good transition
figure;
hold on;
subplot(1,2,1)
plot(kgrid,kdecrule)
title('decision rule of assets');

subplot(1,2,2)
plot(kgrid, hdecrule)
title('decision rule of labor supply');
hold off
saveas(gcf,'Q1.png')
%% 1.2b bad transition
% change to 
% question = 2;
% do again above and plot figure:
figure;
hold on;
subplot(1,2,1)
plot(kgrid,kdecrule)
title('decision rule of assets');

subplot(1,2,2)
plot(kgrid, hdecrule)
title('decision rule of labor supply');
hold off
saveas(gcf,'Q1-bad.png')
%}
%% 1.2.1,2,3 General Equilibrium
disp('General Equilibrium Part')
disp('refer to sample_divisible_labor.m')

%% 1.2.4 Calibration based on  [sample_divisible_labor.m]

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
Gamma  = 0.35;
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
%% A calibration for parameters value of Gamma
% Since my simulation is 10^4 * 10^4, so it will be quiet a bit long

question =1; %for good time transition

% probabilities for 1.1.3 - (a)
if question==1 % i
    pi = pie_g;
else           % ii
    pi = pie_b;
end

zed = [1.8, 0.6];
itera  = 0.01;
maxima = 0.6; 
Xcrit = 1;
gamma = 2;
while (abs(Xcrit -0.97)>0.01) && (itera < maxima)
    itera = itera + 0.01;
    disp('updating value of \Gamma')
    Gamma = 0.5 + itera
    % (inverse) marginal utility functions
    up    = @(c) c.^(-gamma);        % marginal utility of consumption
    invup = @(x) x.^(-1/gamma);      % inverse of marginal utility of consumption
    vp    = @(h) Gamma*h.^(1/varphi);      % marginal disutility of labor supply
    invvp = @(x) (x/Gamma).^(varphi);        % inverse marginal disutility of labor supply 
    r0  = (1/beta-1)-10^(-4);
    residual = 1;
    clear w0 cp0 c0 a0 cpbind cpnon cpnext
    for q = 1:2 
        while abs(residual)>0.01   
            w0 = zed(q)*(1-alpha)*(alpha./(r0+delta)).^(alpha/(1-alpha));
            cp0 = r0*Amat+Ymat*w0;
            cp0(cp0==0)=0.01;
            cp0(cp0==NaN)=0.01;
            cp0(cp0==Inf)=0.01;

            dist    = crit+1;
            % maximum iterations (avoids infinite loops)
            maxiter = 10^(3);
            % counter-
            iter    = 1;

            fprintf('Inner loop, running... \n');
  
            while (dist>crit) && (iter<maxiter)

                c0 = C0(cp0,r0);

                a0 = A0(Amat,Ymat,c0,r0,w0);

                cpbind = ((1+r0)*Amat+b)/2+sqrt(((1+r0)*Amat+b).^(2)+4*(Ymat*w0).^(1+varphi))/2;

                cpnon = zeros(M,N);


                cpnon(:,1)  = interp1(a0(:,1),c0(:,1),Amat(:,1),'spline');
                cpnon(:,2)  = interp1(a0(:,2),c0(:,2),Amat(:,2),'spline');

                cpnext(:,1) = (Amat(:,1)>a0(1,1)).*cpnon(:,1)+(Amat(:,1)<=a0(1,1)).*cpbind(:,1);
                cpnext(:,2) = (Amat(:,2)>a0(1,2)).*cpnon(:,2)+(Amat(:,2)<=a0(1,2)).*cpbind(:,2);

                dist = norm((cpnext-cp0)./cp0);
                iter = iter+1;
            end

            fprintf('Inner loop, done. \n');

            fprintf('Starting simulation... \n');

            clear at yt ct ht 
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
            clear K L
% compute aggregates from market clearing
            K = mean(mean(at(:,T-100:T)));
            L(q) = mean(mean(yt(:,T-100:T).*ht(:,T-100:T)));
            r = alpha*(K/L(q))^(alpha-1)-delta;
            residual = (r-r0)/r0;
            r0 = r;
        end
    end
    disp('(N_h + N_l)/2')
    Xcrit = mean(L)
end
