% =======================================================================
% Quant Macro Part-2 PS-2
% Weimin Zhou
% Date: Nov, 2018
% =======================================================================
clear;clf;close all;
cd '~/Desktop/PS2'  % in order to save png
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
z=[1.6 0.8];
alfa=0.36;
L=[0.96, 0.9]; % ug = 0.04 and ub = 0.1 


v1g = @(k,K,h,H) log( alfa*z(1)*(K/H)^(alfa-1).*k+ (1-alfa)*z(1)*(K/H)^(alfa).*h -delta.*k )/(1-betta);
v1b = @(k,K,h,H) log( alfa*z(2)*(K/H)^(alfa-1)*k+ (1-alfa)*z(2)*(K/H)^(alfa).*h -delta*k )/(1-betta);
v0g = @(k,K,h,H) log( alfa*z(1)*(K/H)^(alfa-1)*k -delta*k )/(1-betta);
v0b = @(k,K,h,H) log( alfa*z(2)*(K/H)^(alfa-1)*k -delta*k )/(1-betta);

%%%%%%%%%%%%% Grid for k and K %%%%%%%%%%%%%%%%%%%%%%%%

% just for a faster computation, technically, we should explore a larger
% grid
k_grid=[0:0.1:5,5.3:0.3:50]; % grid for individual capital: {0,0.1,0.2,...,4.9,5}U{5.3,5.6,5.9,...,50}
K_grid=[16:0.04:18.5];       % aggreagte capital grid: {16,16.04,16.08,16.12,...,18.5}
h_grid=[0.1:0.05:0.9];
H_grid=[0.5:0.05:1];
% K = \integral from 0 to 1 of k_i * di  = E(k) hence, the average of
% individual k


% Evaluation of the VF
% weimin: for each K, we have value function of all small k, 
% for 1g 1b 0g 0b (zz'ee' state combination)
for j=1:size(K_grid,2)
    for m =1:size(H_grid,2)
        for la=1:size(h_grid,2)
V1g(:,j,la,m)= v1g(k_grid,K_grid(j),h_grid(la),H_grid(m))';
V1b(:,j,la,m)= v1b(k_grid,K_grid(j),h_grid(la),H_grid(m))';
V0g(:,j,la,m)= v0g(k_grid,K_grid(j),h_grid(la),H_grid(m))';
V0b(:,j,la,m)= v0b(k_grid,K_grid(j),h_grid(la),H_grid(m))';
        end
    end
end

%%%%%% Perceived law of motion  %%%%%%%%%%%
% initial values
b0g=0;
b1g=1;
b0b=0;
b1b=1;

d0g=0;
d1g=1;
d0b=0;
d1b=1;
%% 2-1 with straight guess version 
H=@(K,zi) exp( (b0g+b1g*log(K))*zi + (b0b+b1b*log(K))*(1-zi) );

G=@(K,zi) exp( (d0g+d1g*log(K))*zi + (d0b+d1b*log(K))*(1-zi) );

c= @(i,I,J,e,g) max(alfa*z(g)*(K_grid(I)/G(K_grid(I),z(g)))^(alfa-1).*k_grid(i)+ ...
              (1-alfa)*z(g)*(K_grid(I)/G(K_grid(I),z(g)))^(alfa)*e*h_grid(J) +(1-delta)*k_grid(i) ...
             - k_grid,0) ;
gamma = 2.0;

for iter=1:500  % VFI         
for i=1:size(k_grid,2)
    for I=1:size(K_grid,2)        
        for J=1:size(h_grid,2)
            for Q=1:size(H_grid,2)
            [dif,Ip]  =min(abs(K_grid - H(K_grid(I),1))); 
            [diff,Ipl]=min(abs(H_grid - G(K_grid(I),1)));
            V0gt(i,I,J,Q)= max(c(i,I,J,0,1).^(1-gamma)./(1-gamma) +...
                betta * ([pize(1,:)]*([V0g(:,Ip,:,Ipl),V1g(:,Ip,:,Ipl),V0b(:,Ip,:,Ipl),V1b(:,Ip,:,Ipl)]'))');
            V1gt(i,I,J,Q)= max(c(i,I,J,0,1).^(1-gamma)./(1-gamma)+...
                betta * ([pize(2,:)]*([V0g(:,Ip,:,Ipl),V1g(:,Ip,:,Ipl),V0b(:,Ip,:,Ipl),V1b(:,Ip,:,Ipl)]'))');  
       
            [dif,Ip]  =min(abs(K_grid - H(K_grid(I),0)));
            [diff,Ipl]=min(abs(H_grid - G(K_grid(I),1)));
            V0bt(i,I,J,Q)= max(c(i,I,J,0,1).^(1-gamma)./(1-gamma) + betta * ([pize(3,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
            V1bt(i,I,J,Q)= max(c(i,I,J,0,1).^(1-gamma)./(1-gamma) + betta * ([pize(4,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');      
            end
        end
    end     
end

dev= max(max(max(abs( [V0gt-V0g,V1gt-V1g,V0bt-V0b,V1bt-V1b]))));

if dev<0.001
    break
else
    V0g=V0gt;
    V1g=V1gt;
    V0b=V0bt;
    V1b=V1bt;
end 
end % VFI         
        
% Recover the policy function only after converging for saving time and memory 
for i=1:size(k_grid,2)
    for I=1:size(K_grid,2)
        
        % approximation next period capital  
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),1))); % H(:,1): good shock
        [V0gt(i,I),a(i,I,2,1)]= max(log(c(i,I,0,1))' + betta * ([pize(1,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
        [V1gt(i,I),a(i,I,1,1)]= max(log(c(i,I,1,1))' + betta * ([pize(2,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');  
       
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),0))); % H(:,1): bad shock
        [V0bt(i,I),a(i,I,2,2)]= max(log(c(i,I,0,2))' + betta * ([pize(3,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
        [V1bt(i,I),a(i,I,1,2)]= max(log(c(i,I,1,2))' + betta * ([pize(4,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
    end 
end
figure

subplot(1,2,1)
hold on
plot(a(:,:,2,1),'m-')
plot(a(:,:,1,1),'c-.')
legend('unemployed in good time','employed in good time')
title('decision rules (policy function) of assets, in good time')
hold off

subplot(1,2,2)
hold on
plot(a(:,:,2,2),'m-')
plot(a(:,:,1,2),'c-.')
legend('unemployed in bad time','employed in bad time')
title('decision rules (policy function) of assets, in bad time')
hold off

saveas(gcf,'Q2_1_ab.png')
%% 
clear all;clc; 
%%   PE version 
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
% v0_{1g} = log(c) + beta* v0_{1g} not changed 
v1g = @(k,K) log( alfa*z(1)*(K/L(1))^(alfa-1).*k+ (1-alfa)*z(1)*(K/L(1))^(alfa) -delta.*k )/(1-betta);
v1b = @(k,K) log( alfa*z(2)*(K/L(2))^(alfa-1)*k+ (1-alfa)*z(2)*(K/L(2))^(alfa) -delta*k )/(1-betta);
v0g = @(k,K) log( alfa*z(1)*(K/L(1))^(alfa-1)*k -delta*k )/(1-betta);
v0b = @(k,K) log( alfa*z(2)*(K/L(2))^(alfa-1)*k -delta*k )/(1-betta);

%%%%%%%%%%%%% Grid for k and K %%%%%%%%%%%%%%%%%%%%%%%%

% just for a faster computation, technically, we should explore a larger
% grid
k_grid=[0:0.1:5,5.3:0.3:50]; % grid for individual capital: {0,0.1,0.2,...,4.9,5}U{5.3,5.6,5.9,...,50}
K_grid=[16:0.04:18.5];       % aggreagte capital grid: {16,16.04,16.08,16.12,...,18.5}
% K = \integral from 0 to 1 of k_i * di  = E(k) hence, the average of
% individual k


% Evaluation of the VF
% weimin: for each K, we have value function of all small k, 
% for 1g 1b 0g 0b (zz'ee' state combination)
for j=1:size(K_grid,2)
V1g(:,j)= v1g(k_grid,K_grid(j))';
V1b(:,j)= v1b(k_grid,K_grid(j))';
V0g(:,j)= v0g(k_grid,K_grid(j))';
V0b(:,j)= v0b(k_grid,K_grid(j))';
end

%%%%%% Perceived law of motion  %%%%%%%%%%%
% initial values
b0g=0;
b1g=1;
b0b=0;
b1b=1;
%B = [0 1 0 1];

for iter_b=1:1000
iter_b

% zi is the index for good shock 
H=@(K,zi) exp( (b0g+b1g*log(K))*zi+ (b0b+b1b*log(K))*(1-zi) );
%log(K') = b0 + b1 * log(K)

% approximation
Ha= @(K,zi) min(abs(K_grid-H(K,zi))); % log(K') - b0 - b1g * log(K) = 0 
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
        
        % take Ip which is the closest grid to the forcast of Aggregate
        % capital from H(which returns the K')
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),1))); % Ha
        
        % pay attention to the proper raw in pize: 
        % 1 - 0g; 2 - 1g; 3 - ob; 4 - 1b
        V0gt(i,I)= max(log(c(i,I,0,1))' + betta * ([pize(1,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
        V1gt(i,I)= max(log(c(i,I,1,1))' + betta * ([pize(2,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');  
       
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),0)));
        V0bt(i,I)= max(log(c(i,I,0,2))' + betta * ([pize(3,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
        V1bt(i,I)= max(log(c(i,I,1,2))' + betta * ([pize(4,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');      
    end     
end

dev= max(max(abs( [V0gt-V0g,V1gt-V1g,V0bt-V0b,V1bt-V1b])));

if dev<0.0001
    break
else
    V0g=V0gt;
    V1g=V1gt;
    V0b=V0bt;
    V1b=V1bt;
end 
end % VFI         
        
% Recover the policy function only after converging for saving time and memory 
for i=1:size(k_grid,2)
    for I=1:size(K_grid,2)
        
        % approximation next period capital  
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),1))); % H(:,1): good shock
        [V0gt(i,I),a(i,I,2,1)]= max(log(c(i,I,0,1))' + betta * ([pize(1,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
        [V1gt(i,I),a(i,I,1,1)]= max(log(c(i,I,1,1))' + betta * ([pize(2,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');  
       
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),0))); % H(:,1): bad shock
        [V0bt(i,I),a(i,I,2,2)]= max(log(c(i,I,0,2))' + betta * ([pize(3,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
        [V1bt(i,I),a(i,I,1,2)]= max(log(c(i,I,1,2))' + betta * ([pize(4,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))');
    end 
end
%% Simulation

% A sequence of TFP
% using the index =1 good ,  =2 bad
if iter_b==1  % Question 2 only one time recovering
zt(1)=1;  
for t=2:2000  % store the TFP zt
    draw=rand;
    zt(t)= 1+(rand>=piZ(zt(t-1),1));
end
% Splitting the sample for good and bad times

% "burning" the first 200 periods
ztb=zt; % no need
ztb(1:200)=0;

% index good times 
i_zg=find(zt==1);

% index bad times
i_zb=find(zt==2);

% initial distribution of assets and employment
indkk = find(k_grid==17);  % 91
N_state(1:960,:,1)=ones(960,1)*[40,1]; % k = 3.9,
%N_state(1:960,:,1)=ones(960,1)*[indkk,1]; % index for k = 17
% =2 unemployed
N_state(961:1000,:,1)=ones(40,1)*[40,2];

K_ind(1)=3;

for t=2:2000
for n=1:1000 
% Evolution of assets
    % optimal saving = a(k,K, employment shock, productivity shock)
    N_state(n,1,t)=a(N_state(n,1,t-1),K_ind(t-1),N_state(n,2,t-1),zt(t-1));
% Evolution of the employment status     
    N_state(n,2,t)= 2-(rand>=pize(1 + zt(t-1)*2 - N_state(n,2,t-1),zt(t)*2-1)/piZ(zt(t-1),zt(t))); 
    % given the probability of moving from zt(t-1) to zt(t)
end
% Storage of the sequence of aggregate capital
[dev2, K_ind(t)]=min(abs(k_grid(round(mean(N_state(:,1,t))))-K_grid));
end

% iter_b larger than 1
else % Question 3, for a update guess on the parameters on Approximation
      
for t=2:2000
for n=1:1000
    
% Evolution of assets: for H changes, a(k,K,e,epsilon) also changes:
    N_state(n,1,t)=a(N_state(n,1,t-1),K_ind(t-1),N_state(n,2,t-1),zt(t-1));
    
   % keep the same evolution of employment
end

% Storage of the sequence of aggregate capital
[dev2, K_ind(t)]=min(abs(k_grid(round(mean(N_state(:,1,t))))-K_grid));


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

% [B1(1:2),s2,s3,s4,s5]=regress(Yb,[ones(size(i_zb(20:end),2),1),log(K_grid(K_ind(i_zb(20:end)-1))')]);R2bad=s5(1);
% [B1(3:4),s2,s3,s4,s5]=regress(Yg,[ones(size(i_zg(20:end),2),1),log(K_grid(K_ind(i_zg(20:end)-1))')]);R2good=s5(1);
% dif_B = max(abs(B-B1));

dev_b = max(abs([b0g-b0gp b1g-b1gp b0b-b0bp b1b-b1bp]))

pause(1)
if dev_b<=0.01 % (d) Iterate until convergence on beta
    break
end
% B = B1*0.1 + B*(0.9);

b0g=0.1*b0gp+0.9*b0g; % updating rules
b1g=0.1*b1gp+0.9*b1g;
b0b=0.1*b0bp+0.9*b0b;
b1b=0.1*b1bp+0.9*b1b;
end
%% e) R^2
SSTg = sum((Yg - mean(Yg).^2));
SSEg = sum((Yg' - Bg'*Xg').^2);
R2_g= 1- SSEg/SSTg

SSTb = sum((Yb - mean(Yb).^2));
SSEb = sum((Yb' - Bb'*Xb').^2);
R2_b= 1- SSEb/SSTb

% another notation:
%disp('R^2 bad aggregate shock:'); R2bad(1)
%disp('R^2 good aggregare shock:'); R2good(1)
%% some graphs

figure

subplot(1,2,1)
hold on
plot(a(:,:,2,1),'m-')
plot(a(:,:,1,1),'c-.')
legend('unemployed in good time','employed in good time')
title('decision rules (policy function) of assets, in good time')
hold off

subplot(1,2,2)
hold on
plot(a(:,:,2,2),'m-')
plot(a(:,:,1,2),'c-.')
legend('unemployed in bad time','employed in bad time')
title('decision rules (policy function) of assets, in bad time')
hold off

saveas(gcf,'Q3_1.png')

% 3-(f) Evolution of the assets distribution
figure
for t_ind=1:100
    
   hist(k_grid(reshape(N_state(:,1,t_ind),1,1000)),40)
   legend(num2str(t_ind))
   pause(1)
end
% the eventual asset distribution store manually, discussion in the pdf.

