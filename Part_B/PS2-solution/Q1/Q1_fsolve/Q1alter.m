%% Alternative 1: revise Luis's sample code
%__________________________________________________________________________
% Quant Macro Part-2
% PS2-Q1: Solution of the model without aggregate risk
% Author: Weimin Zhou
% Date: Nov, 2018
%__________________________________________________________________________

clear;clf;clc;close all;
cd '~/Desktop/PS2'  % in order to save png
%% 0. parameters and preamble set-up
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


%_____________________________________________
%
% Part 1: Aiyagari with endogenous labor supply
% ____________________________________________
%% 1.1 parameters and functional forms

% parameters
z=[1 1];                        % no aggregate shock
alpha  = 0.36;                   % capital income share
b      = 0;                     % borrowing constraint
beta   = 95/100;                % subjective discount factor
delta  = 2.5/1000;                 % depreciation rate of physical capital
gamma  = 3/2; %2.5                   % inverse elasticity of intertemporal substitution
varphi = 2/3;                   % Frisch elasticity of labor supply

rho    = 5/10;                  % persistence parameter prodictivity
N      = 2;                     % number of possible productivity realizations
y1     = 95/100;                % low productivity realization
y2     = 105/100;               % high productivity realization
Gamma  = 0.3;

% transition probability matrix (productivity)
pi      = zeros(N,N);

question = 1;

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
%% 3. endogenous functions

% optimal labor supply
H  = @(c,y,w) invvp(up(c).*y.*w);

% current consumption level, cp0(anext,ynext) is the guess
C0 = @(cp0,r) invup(beta*(1+r)*up(cp0)*pi');
                
% current asset level, c0 = C0(cp0(anext,ynext))
A0 = @(anext,y,c0,r,w) 1/(1+r)...
                    *(c0+anext-H(c0,y,w).*y.*w);
 
%% 1.2 Equilibrium: previously given a price, now to converge on two markets
% c + k = w*n
% find the invariant distribution mu(e,a):
% mu0 = TRAN*mu0 for stationary 
% mu(e1,a1)  =sum( pi_e1e * sum( mu0(a,e)) for all a1 = a'(a,e))

%% 1.1 Household problem 

%=========== Starting values for V ===========%     
% v0_{1g} = log(c) + beta* v0_{1g} not changed 
v1 = @(k) log( z(1)*(r0).*k + (1-alfa)*z(1)*(K/L(1))^(alfa) -delta.*k )/(1-betta);
v0 = @(k) log( alfa*z(1)*(K/L(1))^(alfa-1)*k -delta*k )/(1-betta);

k_grid=[0:0.1:5,5.3:0.3:50]; % grid for individual capital: {0,0.1,0.2,...,4.9,5}U{5.3,5.6,5.9,...,50}

[m n] = size(k_grid);

V0 = zeros(n,n);
V1 = zeros(n,n);

e = [0;1]';

% guess prices
r0  = (1/beta-1);
w0 = (1-alpha)*(alpha/(r0+delta))^(alpha/(1-alpha));

%disp('solving the Aiyagari model under partial equilibrium');
%[consumption]   = stationary(r0,crit,Amat,Ymat,alpha,b,delta,varphi,A0,C0);



x0_em = ((1+r0)*Amat + w0*1)*ones(2,M);
x0_un = ((1+r0)*Amat)*ones(2,M);
% test for initial guess and algorithm for fsovle:
% consumfunc = @(consum)  focdynamic(consum,Amat(1),Amat(1),Gamma,varphi,1,gamma,r0,w0);
% consump(1,1) = max(fsolve(consumfunc,3,optimoptions('fsolve','Display','iter')),0);


% fsolve obtain: c(a,e,a'), n(a,e,c)
for i = 1:M 
    for j = 1:M 
        % employed    
        consumfunc = @(consum) focdynamic(consum,Amat(j),Amat(i),Gamma,varphi,1,gamma,r0,w0);
        consump(i,j,1) = max(fsolve(consumfunc,3,optimoptions('fsolve','Display','iter')),0);
        % unemployed
        consumfunc2 = @(consum)  focdynamic(consum,Amat(j),Amat(i),Gamma,varphi,0,gamma,r0,w0);
        consump(i,j,2) = max(fsolve(consumfunc2,3,optimoptions('fsolve','Display','iter')),0);
    end
end

for i = 1:M
    for j=1:M
        labor(i,j,1) = H(consump(i,j,1),1,w0(1));
        labor(i,j,2) = H(consump(i,j,2),0,w0(1));
    end
end

% plot n and c (not policy func)
% u(c) - v(n)
nu = 2.0;
utilem = log(consump(:,:,1)) - Gamma.* labor(:,:,1).^(1+1/nu)/(1+1/nu);
utilun = log(consump(:,:,2)) - Gamma.* labor(:,:,2).^(1+1/nu)/(1+1/nu);

for iter=1:100  % VFI         
for i=1:size(k_grid,2)
    for I=1:size(k_grid,2)        
        V0t(i,I)= max(log(consump(i,I,1)) - Gamma.* labor(i,I,1).^(1+1/nu)/(1+1/nu))' + beta * ([pi(1,:)]*([V0g(:,I),V1g(:,I)]'))');
        V1t(i,I)= max(log(consump(i,I,2))' - Gamma.* labor(i,I,2).^(1+1/nu)/(1+1/nu))'+ beta * ([pi(2,:)]*([V0g(:,I),V1g(:,I)]'))');       
    end     
end

dev= max(max(abs( [V0t-V0,V1t-V1])));

if dev<0.001
    break
else
    V0=V0t;
    V1=V1t;
end 
end % VFI         
        
% Recover the policy function only after converging for saving time and memory 
for i=1:size(k_grid,2)
    for I=1:size(k_grid,2)
        [V0t(i,I),a(i,I,2)]= max(log(c(i,I,0,2))' + betta * ([pize(3,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
        [V1t(i,I),a(i,I,1)]= max(log(c(i,I,1,2))' + betta * ([pize(4,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
    end 
end

%% 
question =2;
% run again and store policy function
%% some graphs
figure

subplot(1,2,1)
hold on
plot(a(:,:,1),'m-')
plot(a(:,:,2),'c-.')
legend('unemployed','employed')
title('decision rules (policy function) of assets, in bad time transition')
hold off

subplot(1,2,2)
hold on
plot(a1(:,:,1),'m-')
plot(a1(:,:,2),'c-.')
legend('unemployed','employed')
title('decision rules (policy function) of assets, in good time transition')
hold off

saveas(gcf,'Q1.png')