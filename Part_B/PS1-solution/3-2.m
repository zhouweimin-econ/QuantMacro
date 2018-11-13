% =======================================================================
% Quant Macro Part-2 PS-1    
% Question 3-2
% Only update parameters for the other 50% of agents.
% Weimin Zhou
% Date: 10, Nov, 2018
% =======================================================================
clear;clf;close all;
cd '~/Desktop/PS1'  % in order to save png
%% Code for solving a simplified version of K&S (1998)
% Quantitative Macroeconomics - IDEA programme
%% Initial values and parameters
%%%%%%%%%%%% Finding the transition matrix for the state %%%%%%%%%%%%%%

% Weimin: by the answer in my typed-up pdf, we have the following matrix to
% determine 16 unkowns of pi_{zz'ee'}  by computing A*pi=b  which denotes
% pi = inv(A)*b

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
% Weimin: pi_{zz'ee'} : in other words, the markov chain that describes the
% joint evolution of the exogenous shocks

% Interpretation: 
% finding a job is easier if the economy is exiting from a recession, 
% and losing a job ismore likely when the economy is entering a recession

% Weimin: transtion matrix aggregate state

piZ = [ 7/8  1/8;...
        1/8  7/8];
    
%%%%%%%%%%%%  Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%
betta=0.95;
delta=0.0025;
z=[1.01 0.99];
alfa=0.36;
L=[0.96, 0.9]; % ug = 0.04 and ub = 0.1 

%%%%%%%%%%%%% Starting values for V %%%%%%%%%%%%%%%%%%%     

v1g = @(k,K) log( alfa*z(1)*(K/L(1))^(alfa-1).*k+ (1-alfa)*z(1)*(K/L(1))^(alfa) -delta.*k )/(1-betta);
v1b = @(k,K) log( alfa*z(2)*(K/L(2))^(alfa-1)*k+ (1-alfa)*z(2)*(K/L(2))^(alfa) -delta*k )/(1-betta);
v0g = @(k,K) log( alfa*z(1)*(K/L(1))^(alfa-1)*k -delta*k )/(1-betta);
v0b = @(k,K) log( alfa*z(2)*(K/L(2))^(alfa-1)*k -delta*k )/(1-betta);

%%%%%%%%%%%%% Grid for k and K %%%%%%%%%%%%%%%%%%%%%%%%

k_grid=[0:0.1:5,5.3:0.3:50]; % grid for individual capital: {0,0.1,0.2,...,4.9,5}U{5.3,5.6,5.9,...,50}
K_grid=[16:0.04:18.5];       % aggreagte capital grid: {16,16.04,16.08,16.12,...,18.5}

% Evaluation of the VF
% weimin: for each K, we have value function of all small k, 
% for 1g 1b 0g 0b (zz'ee' state combination)
for j=1:size(K_grid,2)
V1g(:,j)= v1g(k_grid,K_grid(j))';
V1b(:,j)= v1b(k_grid,K_grid(j))';
V0g(:,j)= v0g(k_grid,K_grid(j))';
V0b(:,j)= v0b(k_grid,K_grid(j))';
end

for j=1:size(K_grid,2)
V1g_c(:,j)= v1g(k_grid,K_grid(j))';
V1b_c(:,j)= v1b(k_grid,K_grid(j))';
V0g_c(:,j)= v0g(k_grid,K_grid(j))';
V0b_c(:,j)= v0b(k_grid,K_grid(j))';
end

%%%%%% Perceived law of motion  %%%%%%%%%%%
% initial values
b0g=0;
b1g=1;
b0b=0;
b1b=1;

b0g_c=0;
b1g_c=1;
b0b_c=0;
b1b_c=1;
%% 50 % never update, but notice they influence the Aggregate law of motion
for iter_b=1:1000 
iter_b

% zi is the index for good shock 
H=@(K,zi) exp( (b0g+b1g*log(K))*zi+ (b0b+b1b*log(K))*(1-zi) );
% for 50% population not update their guess  
H_c = @(K,zi) exp( (b0g_c+b1g_c*log(K))*zi+ (b0b_c+b1b_c*log(K))*(1-zi));

% approximation
Ha= @(K,zi) min(abs(K_grid-H(K,zi))); % log(K') - b0 - b1g * log(K) = 0 
Ha_c= @(K,zi) min(abs(K_grid-H_c(K,zi)));
%% Solution of the consumer problem
% Consumption for each possible decision
% e=1 employed
% g=1 good times  =2 bad times

% weimin: c = r * k + w*l*epsilon + (1-delta)*k - k'
c= @(i,I,e,g) max(alfa*z(g)*(K_grid(I)/L(g))^(alfa-1).*k_grid(i)+ ...
             (1-alfa)*z(g)*(K_grid(I)/L(g))^(alfa)*e +(1-delta)*k_grid(i) ...
             - k_grid,0) ;
       
for iter=1:1000  % VFI         
for i=1:size(k_grid,2)
    for I=1:size(K_grid,2)        
        % approximation next period capital
        % log-utility not CRRA
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),1))); % Ha
        V0gt(i,I)= max(log(c(i,I,0,1))' + betta * ([pize(1,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
        V1gt(i,I)= max(log(c(i,I,1,1))' + betta * ([pize(2,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');  
       
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),0)));
        V0bt(i,I)= max(log(c(i,I,0,2))' + betta * ([pize(3,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
        V1bt(i,I)= max(log(c(i,I,1,2))' + betta * ([pize(4,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');  
       
        
        [dif_c,Ip_c]=min(abs(K_grid-H_c(K_grid(I),1))); % Ha_unchanged
        V0gt_c(i,I)= max(log(c(i,I,0,1))' + betta * ([pize(1,:)]*([V0g_c(:,Ip_c),V1g_c(:,Ip_c),V0b_c(:,Ip_c),V1b_c(:,Ip_c)]'))');
        V1gt_c(i,I)= max(log(c(i,I,1,1))' + betta * ([pize(2,:)]*([V0g_c(:,Ip_c),V1g_c(:,Ip_c),V0b_c(:,Ip_c),V1b_c(:,Ip_c)]'))');  
        
        [dif_d,Ip_d]=min(abs(K_grid-H_c(K_grid(I),0)));
        V0bt_c(i,I)= max(log(c(i,I,0,2))' + betta * ([pize(3,:)]*([V0g_c(:,Ip_d),V1g_c(:,Ip_d),V0b_c(:,Ip_d),V1b_c(:,Ip_d)]'))');
        V1bt_c(i,I)= max(log(c(i,I,1,2))' + betta * ([pize(4,:)]*([V0g_c(:,Ip_d),V1g_c(:,Ip_d),V0b_c(:,Ip_d),V1b_c(:,Ip_d)]'))');   
    end
end

dev= max(max(abs( [V0gt-V0g,V1gt-V1g,V0bt-V0b,V1bt-V1b])));
dev_c = max(max(abs( [V0gt_c-V0g_c,V1gt_c-V1g_c,V0bt_c-V0b_c,V1bt_c-V1b_c])));

    if dev<0.001 & dev_c < 0.001
        break
    else
        V0g=V0gt;
        V1g=V1gt;
        V0b=V0bt;
        V1b=V1bt;
    
        V0g_c=V0gt_c;
        V1g_c=V1gt_c;
        V0b_c=V0bt_c;
        V1b_c=V1bt_c;
    end 
end % VFI         
        
% Recover the policy function 
for i=1:size(k_grid,2)
    for I=1:size(K_grid,2)
        
        % approximation next period capital  
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),1))); % H(:,1): good shock
        [V0gt(i,I),a(i,I,2,1)]= max(log(c(i,I,0,1))' + betta * ([pize(1,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
        [V1gt(i,I),a(i,I,1,1)]= max(log(c(i,I,1,1))' + betta * ([pize(2,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');  
       
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),0))); % H(:,1): bad shock
        [V0bt(i,I),a(i,I,2,2)]= max(log(c(i,I,0,2))' + betta * ([pize(3,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
        [V1bt(i,I),a(i,I,1,2)]= max(log(c(i,I,1,2))' + betta * ([pize(4,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
   
        % approximation next period capital  
        [dif_c,Ip_c]=min(abs(K_grid-H_c(K_grid(I),1))); % H(:,1): good shock
        [V0gt_c(i,I),a_c(i,I,2,1)]= max(log(c(i,I,0,1))' + betta * ([pize(1,:)]*([V0g_c(:,Ip_c),V1g_c(:,Ip_c),V0b_c(:,Ip_c),V1b_c(:,Ip_c)]'))');
        [V1gt_c(i,I),a_c(i,I,1,1)]= max(log(c(i,I,1,1))' + betta * ([pize(2,:)]*([V0g_c(:,Ip_c),V1g_c(:,Ip_c),V0b_c(:,Ip_c),V1b_c(:,Ip_c)]'))');  
        
        [dif_d,Ip_d]=min(abs(K_grid-H_c(K_grid(I),0))); % H(:,1): bad shock
        [V0bt_c(i,I),a_c(i,I,2,2)]= max(log(c(i,I,0,2))' + betta * ([pize(3,:)]*([V0g_c(:,Ip_d),V1g_c(:,Ip_d),V0b_c(:,Ip_d),V1b_c(:,Ip_d)]'))');
        [V1bt_c(i,I),a_c(i,I,1,2)]= max(log(c(i,I,1,2))' + betta * ([pize(4,:)]*([V0g_c(:,Ip_d),V1g_c(:,Ip_d),V0b_c(:,Ip_d),V1b_c(:,Ip_d)]'))');
    end 
end
%% Simulation

% A sequence of TFP
% using the index =1 good ,  =2 bad
if iter_b==1  
zt(1)=1;  
for t=2:2000  % store the TFP zt
    draw=rand;
    zt(t)= 1+(rand>=piZ(zt(t-1),1)); % a rand draw of TFP 
end
% Splitting the sample for good and bad times

% "burning" the first 200 periods
ztb=zt;
ztb(1:200)=0;  % ztb unused

% index good times 
i_zg=find(zt==1);

% index bad times
i_zb=find(zt==2);

% initial distribution of assets and employment
N_state(1:460,:,1)=ones(460,1)*[40,1];
% =2 unemployed
N_state(461:500,:,1)=ones(40,1)*[40,2];

N_state_c(1:460,:,1)=ones(460,1)*[40,1];
N_state_c(461:500,:,1)=ones(40,1)*[40,2];
%randInd = randperm(size(N_state,1));
%N_state(:,:,1) = N_state(randInd,:,1); % random ordered 

K_ind(1)=3; 
% K_ind(1)=mean(K_grid);

for t=2:2000
for n=1:500   

% smart people:
% Evolution of assets
    N_state(n,1,t)=a(N_state(n,1,t-1),K_ind(t-1),N_state(n,2,t-1),zt(t-1));
% Evolution of the employment status     
    N_state(n,2,t)= 2-(rand>=pize(1 + zt(t-1)*2 - N_state(n,2,t-1),zt(t)*2-1)/piZ(zt(t-1),zt(t)));  
 
% dumb people:    
    N_state_c(n,1,t)=a_c(N_state(n,1,t-1),K_ind(t-1),N_state_c(n,2,t-1),zt(t-1)); 
    
    % employment should be the same
    N_state_c(n,2,t)= 2-(rand>=pize(1 + zt(t-1)*2 - N_state_c(n,2,t-1),zt(t)*2-1)/piZ(zt(t-1),zt(t))); 
end
N_state_d = [ N_state; N_state_c];
% Storage of the sequence of aggregate capital
[dev2, K_ind(t)]=min(abs(k_grid(round(mean(N_state_d(:,1,t))))-K_grid));
end

% iter_b > 1
else % Question 3, for an update guess on the parameters on Approximation
      
for t=2:2000
for n=1:500
% Evolution of assets for smart people
    N_state(n,1,t)=a(N_state(n,1,t-1),K_ind(t-1),N_state(n,2,t-1),zt(t-1));

% Evolution of assets for dumb people
    N_state_c(n,1,t)=a_c(N_state_c(n,1,t-1),K_ind(t-1),N_state_c(n,2,t-1),zt(t-1));
end
N_state_d = [ N_state; N_state_c];
% Storage of the sequence of aggregate capital
[dev2, K_ind(t)]=min(abs(k_grid(round(mean(N_state_d(:,1,t))))-K_grid));
end
end

% Question 3, for a update guess on the parameters on Approximation

% Regression model for the evolution of aggregate capital
% regression for good times (burning the first 20 periods of g times)

Yg=log(K_grid(K_ind(i_zg(20:end)))');
Xg=[ones(size(i_zg(20:end),2),1),log(K_grid(K_ind(i_zg(20:end)-1))')] ;   
Bg=Xg\Yg
b0gp=Bg(1);
b1gp=Bg(2);

% regression for bad times (burning the first 20 periods of bad times.

Yb=log(K_grid(K_ind(i_zb(20:end)))');
Xb=[ones(size(i_zb(20:end),2),1),log(K_grid(K_ind(i_zb(20:end)-1))')]  ;  
Bb=Xb\Yb
b0bp=Bb(1);
b1bp=Bb(2);

dev_b = max(abs([b0g-b0gp b1g-b1gp b0b-b0bp b1b-b1bp]))

pause(1)
if dev_b<=0.01 % (d) Iterate until convergence on beta for 50% population
    break
end

b0g=0.1*b0gp+0.9*b0g; % updating rules for 50% agents
b1g=0.1*b1gp+0.9*b1g;
b0b=0.1*b0bp+0.9*b0b;
b1b=0.1*b1bp+0.9*b1b;
end
save Solution_to_model;
%% e) R^2
SSTg = sum((Yg - mean(Yg).^2));
SSEg = sum((Yg' - Bg'*Xg').^2);
R2_g= 1- SSEg/SSTg

SSTb = sum((Yb - mean(Yb).^2));
SSEb = sum((Yb' - Bb'*Xb').^2);
R2_b= 1- SSEb/SSTb
%% some graphs
%ind_k = find(k_grid==17);
ind_K = find(K_grid==17.16);

figure

subplot(2,2,1)
hold on
plot(k_grid,a(:,ind_K,2,1),'m-')
plot(k_grid,a(:,ind_K,1,1),'c-.')
legend('unemployed in good time','employed in good time')
title('g(a), Aggregate K = 17.16, in good time')
hold off

subplot(2,2,2)
hold on
plot(k_grid,a(:,ind_K,2,2),'m-')
plot(k_grid,a(:,ind_K,1,2),'c-.')
legend('unemployed in bad time','employed in bad time')
title('g(a) Aggregate K = 17.16, in bad time')
hold off

subplot(2,2,3)
hold on
plot(k_grid,a_c(:,ind_K,2,1),'m-')
plot(k_grid,a_c(:,ind_K,1,1),'c-.') 
legend('unemployed in good time','employed in good time')
title('non-update population, g(a), Aggregate K = 17.16, in good time')
hold off

subplot(2,2,4)
hold on
plot(k_grid,a_c(:,ind_K,2,2),'m-')
plot(k_grid,a_c(:,ind_K,1,2),'c-.')
legend('unemployed in bad time','employed in bad time')
title('non-update population g(a), Aggregate K = 17.16, in bad time')
hold off
saveas(gcf,'Q3_2_a.png')


%% Welfare

for i = 1:size(k_grid,2) 
  W_1g(i) = log(alfa*z(1)*(K_grid(ind_K)/L(1))^(alfa-1).*k_grid(i) + (1-alfa)*z(1)*(K_grid(ind_K)/L(1))^(alfa) +(1-delta).*k_grid(i) - a(i,ind_K,2,1))... 
  + (betta)*(([pize(1,:)]*([V0g(i,ind_K),V1g(i,ind_K),V0b(i,ind_K),V1b(i,ind_K)]'))');   
  W_1b(i) =  log(alfa*z(2)*(K_grid(ind_K)/L(2))^(alfa-1).*k_grid(i)+ (1-alfa)*z(2)*(K_grid(ind_K)/L(2))^(alfa) +(1-delta).*k_grid(i) - a(i,ind_K,2,1) )...
      + (betta)*(([pize(1,:)]*([V0g(i,ind_K),V1g(i,ind_K),V0b(i,ind_K),V1b(i,ind_K)]'))');   
  W_0g(i) = log(alfa*z(1)*(K_grid(ind_K)/L(1))^(alfa-1).*k_grid(i)+ (1-delta).*k_grid(i)- a(i,ind_K,2,1) )...
      + (betta)*(([pize(1,:)]*([V0g(i,ind_K),V1g(i,ind_K),V0b(i,ind_K),V1b(i,ind_K)]'))');   ;
  W_0b(i) = log(alfa*z(2)*(K_grid(ind_K)/L(2))^(alfa-1).*k_grid(i) +(1-delta).*k_grid(i)- a(i,ind_K,2,1))...
      + (betta)*(([pize(1,:)]*([V0g(i,ind_K),V1g(i,ind_K),V0b(i,ind_K),V1b(i,ind_K)]'))');   ;
end
%%
for i = 1:size(k_grid,2) 
  W_1g_c(i) = log(alfa*z(1)*(K_grid(ind_K)/L(1))^(alfa-1).*k_grid(i) + (1-alfa)*z(1)*(K_grid(ind_K)/L(1))^(alfa) +(1-delta).*k_grid(i) - a(i,ind_K,2,1))... 
  + (betta)*(([pize(1,:)]*([V0gt_c(i,ind_K),V1gt_c(i,ind_K),V0bt_c(i,ind_K),V1bt_c(i,ind_K)]'))');   
  W_1b_c(i) =  log(alfa*z(2)*(K_grid(ind_K)/L(2))^(alfa-1).*k_grid(i)+ (1-alfa)*z(2)*(K_grid(ind_K)/L(2))^(alfa) +(1-delta).*k_grid(i) - a(i,ind_K,2,1) )...
      + (betta)*(([pize(1,:)]*([V0gt_c(i,ind_K),V1gt_c(i,ind_K),V0bt_c(i,ind_K),V1bt_c(i,ind_K)]'))');   
  W_0g_c(i) = log(alfa*z(1)*(K_grid(ind_K)/L(1))^(alfa-1).*k_grid(i)+ (1-delta).*k_grid(i)- a(i,ind_K,2,1) )...
      + (betta)*(([pize(1,:)]*([V0gt_c(i,ind_K),V1gt_c(i,ind_K),V0bt_c(i,ind_K),V1bt_c(i,ind_K)]'))');   ;
  W_0b_c(i) = log(alfa*z(2)*(K_grid(ind_K)/L(2))^(alfa-1).*k_grid(i) +(1-delta).*k_grid(i)- a(i,ind_K,2,1))...
      + (betta)*(([pize(1,:)]*([V0gt_c(i,ind_K),V1gt_c(i,ind_K),V0bt_c(i,ind_K),V1bt_c(i,ind_K)]'))');   ;
end
%% compare the welfare of two kinds of agents: 
figure
subplot(2,2,1)
hold on
plot(k_grid,W_0g,'m-')
plot(k_grid,W_0g_c,'c-.')
legend('best-fit agents','never-update agents')
title('Welfare Comparison, unemployed, in good time')
hold off

subplot(2,2,2)
hold on
plot(k_grid,W_1g,'m-')
plot(k_grid,W_1g_c,'c-.')
legend('best-fit agents','never-update agents')
title('Welfare Comparison, employed, in good time')
hold off

subplot(2,2,3)
hold on
plot(k_grid,W_0b,'m-')
plot(k_grid,W_0b_c,'c-.')
legend('best-fit agents','never-update agents')
title('Welfare Comparison, unemployed, in bad time')
hold off

subplot(2,2,4)
hold on
plot(k_grid,W_1b,'m-')
plot(k_grid,W_1b_c,'c-.')
legend('best-fit agents','never-update agents')
title('Welfare Comparison, employed, in bad time')
hold off
saveas(gcf,'Q3_2_b.png')