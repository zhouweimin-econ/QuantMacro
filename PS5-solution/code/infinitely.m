% =======================================================================
% Quant Macro PS-5 (Final Project)  
% Weimin Zhou
% Due: 1, Nov, 2018
% This file revised to generally cover Question 2,3,4
% Infinitely case 
% =======================================================================
clear;clf;close all;
cd '~/Desktop/PS5/code'  % in order to save png
%% Set up guide 
%  continous method: continous = 1 
%  discrete  method: continous = 0
%
%  transition income matrix parameter: gamma 
%  CRRA risk-averse parameter: sigma
%  Quadratic utility satisfaction level: cbar
%  Productivity shocks deviation: sigma_y
%% parameters set up  
sigma 	  = 2; 
cbar      = 100;

%continous = 1;                  % continous method
continous = 0;                   % discrete method       
r         = 0.04;                % interest rate of guessing
rho       = 0.06;
beta      = inv(1+rho);          % discount factor 

%sigma_y = 0; 			         % certainty case: Y = { 1 }
sigma_y  = 0.1; 				 % uncertainty case: Y = {1-sigma_y, 1+sigma_y}

gamma    = 0;

Pi		 = [(1+gamma)/2, (1-gamma)/2;
           (1-gamma)/2, (1+gamma)/2];	   
       
% method = 1; % Quadratic
% method = 2; % log 
method   = 3;   % CRRA 

if method == 1; % Quadratic utility 
    utility = @(c)  -.5 * (c-cbar)^2;
else 
	if method == 2;  % log utility: \sigma = 1 in this case
		utility = @(c) log(c);
    else   
		utility = @(c) (c^(1-sigma)-1)/(1-sigma);
	end
end
%% Discretize state space A*Y 
% for shocks: 
Y   = [1-sigma_y,1+sigma_y];
ny = size(Y,2);

% form the grid:
nk	= 1000;
dev	= 0.9;
asmin	= (10-dev); % only change the location but not the shape, hence choose 1
asmax	= (10+dev);
devk	= (asmax-asmin)/(ny-1);
agrid	= linspace(asmin,asmax,nk)'; % case 2: a_{t+1} >= 0

% since never reach to natural borrowing constraint, so I ignore this part:
agrid2  = agrid;
agrid2(1,1)= -r/(1+r)*Y(1,1);  % case 1: a_{t+1} >= -A: a1 = -A

% Initialization of policy vector
c		= zeros(nk,ny);
u		= zeros(nk,ny);
v		= zeros(nk,ny);
v0	    = zeros(nk,ny);

crit	= 1;
iter	= 1;
maxit   = 10000;

% piecewise linear splines interpolation nodes
aspline  = agrid; %as previous one

% u(c)
for j=1:ny;
   c		= repmat(Y(j)+(1+r)*agrid',nk,1)-repmat(agrid,1,nk);
   unfes		= find(c<=0);     
   c(unfes)	= NaN;                            % impose constraint
   eval(['util' int2str(j) '= utility(c);'])
   eval(['util' int2str(j) '(unfes)= -inf;']) % store utility as util
   clear c unfes;
end

while crit>1e-6;
   for j=1:ny;
		if continous== 1; %	nodes from continuous methods 
			vi(:,j) = interp1(agrid,v(:,j),aspline,'spline'); % piecewise cubic spline interpolation
			v(:,j)  = vi(:,j);
		else 
			 % discretize
        end   
      % VFI: save time to use repmat for reducing dimension inside loop 
      eval(['[tmpv,tmpdecis]	= max(util' int2str(j) '+beta*repmat(v*Pi(j,:)'',1,nk));']); % [size(v*Pi(j,:),1)*1, size(v*Pi(j,:),2)*nk].
      v0(:,j)=tmpv';
      dr(:,j)=tmpdecis';  % optimal decision 
   end
   crit	= max(max(abs(v0-v)));
   v    = v0;   % Value function
   iter	= iter+1;
end
% dr is the a'(a,y);  
% compute g(c):

% a'(an,y2) < an:
anext = agrid(dr);

alast = agrid(1000);
if dr(1000,2) < alast
    disp('True')
else
    disp('Wrong')
end
% True 

% compute c(a,y)
for j=1:ny;
   c(:,j)	= Y(j)+(1+r)*agrid-anext(:,j);
end
%% Simulate the path 
Amark = 1;
A0    = agrid(Amark,1);    % Initial level of assets
n	  = 1000;              % Number of periods to simulate
s0 	  = 1;                 % Initial level of productivity for markov 
states   = zeros(n-1,2);
controls = zeros(n-1,2);
% markov state computation from https://sites.google.com/site/katrinrabitsch/
[chain,state] = markov(Pi,n,s0); 
for i = 1:n-1;
    if chain(i) == 1; % in state y1 
       Aprime = dr(Amark,1); % a' = g(a,1);
       consum   = ((1+r)*A0+Y(2)-Aprime); 
       Amark  = Aprime; 
    elseif chain(i) == 2; % in state y2
       Aprime = dr(Amark,2);
       consum   = ((1+r)*A0+Y(1)-Aprime); 
       Amark  = Aprime;
    else
      % Nothing 
    end;
    states(i,:) = [ A0 chain(i) ]; 
    controls(i,:) = [ consum Aprime ]; % c, a'
    A0 = Aprime; % update a by a = a'
end

%%  II.4.1.1 and II.4.2.1 gamma = 0, sigma_y = 0 / 0.1, plot consumption functions

% c_41_quad = c;
c_41_CRRA = c;
%%
figure 
subplot(1,2,1)
plot(agrid,c_41_CRRA)
title('C(a,y) for CRRA, T=\infty, certainty with \sigma_y =0')
subplot(1,2,2)
plot(agrid,c_41_quad)
title('C(a,y) for Quadratic, T=\infty, certainty with \sigma_y =0')

saveas(gcf,'41a.png')
%% II.4.2 Q3 sigma = 2, 5, 20
%c2=c;
%c5=c;
c20=c;
%%
figure 
subplot(1,3,1)
plot(agrid,c2(:,1),agrid,c5(:,1),agrid,c20(:,1))
legend('1-\sigma_y','1+\sigma_y')
title('C(a,y), CRRA, T=\infty, uncertainty, \sigma =2')
subplot(1,3,2)
plot(agrid,c5)
legend('1-\sigma_y','1+\sigma_y')
title('C(a,y), CRRA, T=\infty, uncertainty, \sigma =5')
subplot(1,3,3)
plot(agrid,c20)
legend('1-\sigma_y','1+\sigma_y')
title('C(a,y), CRRA, T=\infty, uncertainty, \sigma =20')
saveas(gcf,'43.png')
%%
times=1:1:size(controls,1);
CRRA_controls = controls(:,1); % (:,1) = (:,2) for certainty
%QUAD_controls = controls(:,1); 
figure 
plot(times,CRRA_controls,times,QUAD_controls)
ylim=([-50 50]);
title('time path of consumption for two utility functions')
saveas(gcf,'41q.png')
%%
% based on above codes, we could plot the II4.2 Quesiton 3, 4 and 5 in this
% infinite-lived case.