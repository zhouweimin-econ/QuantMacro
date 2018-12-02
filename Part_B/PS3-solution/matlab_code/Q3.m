%% Calibration Q3
clear all;clc
%% Loading parameters
% Discount factor $\beta=0.95$

betta=0.95;
% 
% Possible values for GDP $y\in\{0.9,1,1.05\}$

y_grid=[0.95, 1.1, 1.5];
% 
% Transition matrix for GDP
% 
% $$\pi_{yy'}=\left[\begin{array}{ccc}0.5 & 0.3 & 0.2\\0.1 & 0.6 & 0.3\\0.2 
% & 0.4 & 0.4\end{array}\right]$$

piy=[0.5, 0.3, 0.2;...
     0.05, 0.65, 0.3;...
     0.02, 0.55, 0.43];
% 
% Risk aversion parameter $\sigma=2$

sig=1.5;
% 
% Utility function $u(c)=\frac{c^{1-\sigma}}{1-\sigma}$

u = @(c) c.^(1-sig)./(1-sig);
%
% Risk-free interest rate $R=1$

R=1;
%
% Probability of regaining access to capital markets next period $\lambda=0.3$

lamda=0.5;
% 
% GDP loss during default. Varies with the level of GDP.
%tau = 0.2;
tau=[0.1,0.4,0.5];
%% simulation 
clear tau betta 
M1=linspace(0.01,0.2,30);
M2=linspace(0.3,0.5,30);
M3=linspace(0.4,0.6,30);
%%
% betta=0.75;
beta = linspace(0.9,0.95,30);
%for iterrr = 1:20
%    betta = beta(iterrr)
%for iterbeta = 1:10
%    disp('beta iteration')
%    betta = beta(iterbeta)
     disp('tau value selection')
     tau=[0.0362, 0.2276, 0.4276];
     
for iterrr=1:30
     iterrr
     betta = beta(iterrr)
%% Initial values and the discretized state space
% Possible levesl of debt issuance $b\in B= \{0,0.05,0.1,...,0.5\}$

B=0:0.05:2.8;

max_iter=1000;
max_iterq=1000;
%% 
% Initial values for the value functions:
% 
% We set the inital guess by assuming that GDP and the level of debt do not 
% change over time. 
% 
% Value function of no default 
% 
% $$V_{0}(b,y)=\frac{u(y)}{1-\beta}$$

V= ones(size(B,2),1)*u(y_grid)/(1-betta);
%% 
% Value function of default
% 
% $${V}_{0}(y)=\frac{u(y)}{1-\beta}$$

Vd = u((1-tau).*y_grid)/(1-betta);
%% 
% 
%% Equilibrium computation
% Our Initial guess for the bond-price schedule is constructed assuming the 
% government never defaults $q(b',y)=\frac{1}{R}$

q=ones(size(B,2),size(y_grid,2));
%% Iterations on the price schedule

for iter_q=1:max_iterq
 
%
%% Value function iterations (taking   as given)


for iter=1:max_iter

    for iy=1:size(y_grid,2)
    
%% 
% $$V_{i+1}(y)=u\left((1-\tau)y\right)+(1-\lambda)\beta E\left\{ V_i(y')\mid 
% y\right\} +\lambda\beta E\left\{ V_{i}(0,y')\mid y\right\} $$

    Vd(iy)=u((1-tau(iy))*y_grid(iy))+(1-lamda)*betta*piy(iy,:)*Vd'+lamda*betta*V(1,:)*piy(iy,:)';
    
        for ib=1:size(B,2)
    
    V_old=V;
    Vd_old=Vd; 
        V(ib,iy)=max(max(u(max(y_grid(iy)+q(:,iy).*B'-B(ib),0))+betta*V*piy(iy,:)'),Vd(iy)) ;
        end
    end
    
 dev= max(max(abs([V_old-V;Vd_old-Vd])));
if dev<=0.00001
    break
end
end

%% 
% Updating the bond price menu

 for iy=1:size(y_grid,2)
    for ib=1:size(B,2)
        
        q(ib,iy)=1-piy(iy,:)*(V(ib,:)<=Vd)';
    end
 end
 
end
%% Bond price menu

%figure
%plot_q(B,q)
%% Simulation of the model
% recovering the policy function

 for iy=1:size(y_grid,2)
       for ib=1:size(B,2)
       % No default value and debt issuance (conditional on no deault)
        [V_ND(ib,iy),bp(ib,iy)]=max(u(max(y_grid(iy)+q(:,iy).*B'-B(ib),0))+betta*V*piy(iy,:)')  ; 
        % Default decision
        Dp(ib,iy)=(V_ND(ib,iy)<=Vd(iy)) ;   
        
    end
 end
%% 
% Simuated sequence of GDP
% 
%  starting value (index)

yt=1;

for t=2:500
    draw_t=rand;
    yt(t)=1+(draw_t>=piy(yt(t-1),1))+(draw_t>=sum(piy(yt(t-1),1:2)));
end

%% 
% initial level of debt
% 
% index in B

bt=ones(500,1);
%% 
% index for the default decision =1 and the default state

Def_b=nan(1,500);
Def_state(500)=0;

for t=2:500
   
   if Def_state(t-1)==0
   
%% 
%     default decision (decided at t)

   Def_b(t)=Dp(bt(t-1),yt(t));
   
   if Def_b(t)==0
 
%% 
%    Debt issuance decision (decided at t)  

 
   bt(t)=bp(bt(t-1),yt(t));
   
       Def_state(t)=0;
   else
       bt(t)=1;
       Def_state(t)=1;
   end
   
   elseif rand<=lamda
   
       Def_b(t)=Dp(bt(t-1),yt(t));
   
   if Def_b(t)==0
  
%% 
%      Debt issuance decision (decided at t)   

   bt(t)=bp(bt(t-1),yt(t));
   Def_state(t)=0;
   else
       bt(t)=1;
       Def_state(t)=1;
   end
   
   else
       
       Def_state(t)=1;
   bt(t)=1;
   end
   
 
%% 
%    0bserved risk spread (1/q-1)

   r_spread(t)=1/q(bt(t),yt(t))-1;
%% 
%    Default probability

   p_model(t)=1-q(bt(t),yt(t));
   
end
       
%% 
%  graph of default and the risk premia

% risk_premia_graph(Def_state, r_spread)
%% Estimating a logit
%  Arranging the data. 
% 
%  We have to be sure that if a default spell lasts for more than one period, 
% we only include the first time default was declared (this is why we
% 
%  distinguish between the default decision and the default state)
% 
% 
% 
% % Number of observations 

N=sum(1-(isnan(Def_b)));
% Regressors. Constant and debt/GDP
X=[B(bt)'./y_grid(yt)'];

% Default decision
Y=Def_b';

% Logit regression. Binomial outcome (0, 1)
[par_est,dev,stats]=glmfit(X(1:end-1),Y(2:end),'binomial');


% Estimated Default probability
X_grid=0:0.1:2;
p_est=1./(1+exp(-par_est(1)-par_est(2)*X_grid));




% default probability in the simulation
% by checking, and intuition:
% this is also the percentage of time in default:
p_est_sim=1./(1+exp(-par_est(1)-par_est(2)*X));

disp('---------------')
disp('Percentage of time in default')
mean(p_est_sim)
disp('---------------')

disp('---------------')
disp('Coefficients for Logit Regression')
disp('True data coeff: [-2.199848***  1.72411***]')
par_est
disp('---------------')
if par_est(1) + 2.199848 < 0.1 && par_est(2) - 1.72411 < 0.1
    break
    disp('----Simulation Results----')
    disp('tau:')
    tau
    disp('betta:')
    betta
    disp('__________________________')
%   else
%    clear all;
%    y_grid=[0.95, 1.1, 1.5];
%    piy=[0.5, 0.3, 0.2;...
%     0.05, 0.65, 0.3;...
%     0.02, 0.55, 0.43];
%    sig=1.5;
%    u = @(c) c.^(1-sig)./(1-sig);
%    R=1;
%    lamda=0.5;
%M1=linspace(0.01,0.2,30);
%M2=linspace(0.3,0.5,30);
%M3=linspace(0.4,0.6,30);
%    beta = linspace(0.5,0.9,10);
end
end

%end
