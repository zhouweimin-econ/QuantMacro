% =======================================================================
% Quant Macro Part-2 PS-1
% Weimin Zhou
% Date: 10, Nov, 2018
% =======================================================================
clear;clf;close all;
cd '~/Desktop/PS2'  % in order to save png
%% Initial values and parameters
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

piZ = [ 7/8  1/8;...
        1/8  7/8];
    
%  Parameters 
betta=0.95;
delta=0.0025;
z=[1.01 0.99];
alfa=0.36;
L=[0.96, 0.9]; % ug = 0.04 and ub = 0.1 


% just for a faster computation, technically, we should explore a larger
% grid
k_grid=[0:0.1:5,5.3:0.3:50]; % grid for individual capital: {0,0.1,0.2,...,4.9,5}U{5.3,5.6,5.9,...,50}
K_grid=[16:0.04:18.5];       % aggreagte capital grid: {16,16.04,16.08,16.12,...,18.5}
h_grid=[0.2:0.05:0.9];
% K = \integral from 0 to 1 of k_i * di  = E(k) hence, the average of
% individual k
d0g=0;d1g=1;d0b=0;d1b=1;
% initial values
b0g=0;b1g=1;b0b=0;b1b=1;


H=@(K,zi) exp( (b0g+b1g*log(K))*zi + (b0b+b1b*log(K))*(1-zi) );

G=@(K,zi) exp( (d0g+d1g*log(K))*zi + (d0b+d1b*log(K))*(1-zi) );

c= @(i,I,m,e,g) max(alfa*z(g)*(K_grid(I)/G(K_grid(I),z(g)))^(alfa-1).*k_grid(i)+ ...
              (1-alfa)*z(g)*(K_grid(I)/G(K_grid(I),z(g)))^(alfa)*e*h_grid(m) +(1-delta)*k_grid(i) ...
             - k_grid,0) ;

% Starting values for V     
v1g = @(k,K,h) log( alfa*z(1)*(K/G(K,1))^(alfa-1).*k + (1-alfa)*z(1)*(K/G(K,1))^(alfa).*h -delta.*k )./(1-betta);
v1b = @(k,K,h) log( alfa*z(2)*(K/G(K,0))^(alfa-1).*k  + (1-alfa)*z(2)*(K/G(K,0))^(alfa).*h -delta.*k  )./(1-betta);
v0g = @(k,K,h) log( alfa*z(1)*(K/G(K,1))^(alfa-1).*k -delta.*k )./(1-betta);
v0b = @(k,K,h) log( alfa*z(2)*(K/G(K,0))^(alfa-1).*k -delta.*k )./(1-betta);
%____________________________________________________________
%
% Grid for k and K 
%____________________________________________________________

for I=1:size(K_grid,2)
    for m =1:size(h_grid,2)
        V1g(:,I,m)= v1g(k_grid,K_grid(I),h_grid(m))';
        V1b(:,I,m)= v1b(k_grid,K_grid(I),h_grid(m))';
        V0g(:,I,m)= v0g(k_grid,K_grid(I),h_grid(m))';
        V0b(:,I,m)= v0b(k_grid,K_grid(I),h_grid(m))';
    end
end
%% 2-1 with straight guess version 
tic
for iter=1:500  % VFI  
    disp('staring generating')
for i=1:size(k_grid,2)
    for I=1:size(K_grid,2)        
       for m=1:size(h_grid,2)
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),1))); 
        V0gt(i,I,m)= max(log(c(i,I,m,0,1))' + betta * ([pize(1,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');
        V1gt(i,I,m)= max(log(c(i,I,m,1,1))' + betta * ([pize(2,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');  
       
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),0)));
        V0bt(i,I,m)= max(log(c(i,I,m,0,2))' + betta * ([pize(3,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');
        V1bt(i,I,m)= max(log(c(i,I,m,1,2))' + betta * ([pize(4,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');      
       end 
    end
end
disp('deviation: ')
dev = max(max(max(abs( [V0gt-V0g,V1gt-V1g,V0bt-V0b,V1bt-V1b]))))

if dev<0.001
    break
else
    V0g=V0gt;
    V1g=V1gt;
    V0b=V0bt;
    V1b=V1bt;
end 
end % VFI  
toc
% running: 33 mins
% dev: 9.7382e-04
%%   
% Recover the policy function only after converging for saving time and memory 
for i=1:size(k_grid,2)
    for I=1:size(K_grid,2)
        for m =1:size(h_grid,2)
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),1))); 
        [V0gt(i,I,m),a(i,I,m,2,1)]= max(log(c(i,I,m,0,1))' + betta * ([pize(1,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');
        [V1gt(i,I,m),a(i,I,m,1,1)]= max(log(c(i,I,m,1,1))' + betta * ([pize(2,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');  
       
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),0)));
        [V0bt(i,I,m),a(i,I,m,2,2)]= max(log(c(i,I,m,0,2))' + betta * ([pize(3,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');
        [V1bt(i,I,m),a(i,I,m,1,2)]= max(log(c(i,I,m,1,2))' + betta * ([pize(4,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');      
        end 
    end
end

figure

subplot(1,2,1)
hold on
plot(a(:,:,5,2,1),'m-')
plot(a(:,:,5,1,1),'c-.')
legend('unemployed in good time','employed in good time')
title('decision rules (policy function) of assets, in good time (labor = 0.4)')
hold off

subplot(1,2,2)
hold on
plot(a(:,:,5,2,2),'m-')
plot(a(:,:,5,1,2),'c-.')
legend('unemployed in bad time','employed in bad time')
title('decision rules (policy function) of assets, in bad time (labor = 0.4)')
hold off

saveas(gcf,'Q2_1_plot.png')
%% labor 
emp = [1 0];
myfunc = @(c) consumption(c,k_grid(1),a(1,K_grid(1),5,2,1),K_grid(1),emp(1),z(1)); 
c0 = 3; % see Q1alter2.m for a soundful guess.
options = optimset('display','iter','TolX',1e-8,'MaxIter',20);
% c(1,1,1,2,1)=fsolve(myfunc,c0,options);

for i=1:size(k_grid,2)
    for I=1:size(K_grid,2)
        for m =1:size(h_grid,2)          
        myfunc =  @(c) consumption(c,k_grid(i),a(1,K_grid(I),h_grid(m),2,1),K_grid(h),emp(2),z(1)); 
        c(i,I,m,2,1)=fsolve(myfunc, c0, options);
        
        myfunc1 = @(c) consumption(c,k_grid(i),a(1,K_grid(I),h_grid(m),1,1),K_grid(h),emp(1),z(1));
        c(i,I,m,1,1)=fsolve(myfunc, c0, options);
        
        myfunc2 =  @(c) consumption(c,k_grid(i),a(1,K_grid(I),h_grid(m),2,1),K_grid(h),emp(2),z(2)); 
        c(i,I,m,2,2)=fsolve(myfunc, c0, options);
        
        myfunc3 = @(c) consumption(c,k_grid(i),a(1,K_grid(I),h_grid(m),1,1),K_grid(h),emp(1),z(2));
        c(i,I,m,1,2)=fsolve(myfunc, c0, options);
        end 
    end
end
%%   PEA version for updating 
clear v1g v1b v0g v0b
v1g = @(k,K) log( alfa*z(1)*(K/L(1))^(alfa-1).*k+ (1-alfa)*z(1)*(K/L(1))^(alfa) -delta.*k )/(1-betta);
v1b = @(k,K) log( alfa*z(2)*(K/L(2))^(alfa-1)*k+ (1-alfa)*z(2)*(K/L(2))^(alfa) -delta*k )/(1-betta);
v0g = @(k,K) log( alfa*z(1)*(K/L(1))^(alfa-1)*k -delta*k )/(1-betta);
v0b = @(k,K) log( alfa*z(2)*(K/L(2))^(alfa-1)*k -delta*k )/(1-betta);

clear V1g V1b V0g V0b
for j=1:size(K_grid,2)
V1g(:,j)= v1g(k_grid,K_grid(j))';
V1b(:,j)= v1b(k_grid,K_grid(j))';
V0g(:,j)= v0g(k_grid,K_grid(j))';
V0b(:,j)= v0b(k_grid,K_grid(j))';
end
%% this take extremely long time: 33mins * 1000
tic

H=@(K,zi) exp( (b0g+b1g*log(K))*zi + (b0b+b1b*log(K))*(1-zi) );

G=@(K,zi) exp( (d0g+d1g*log(K))*zi + (d0b+d1b*log(K))*(1-zi) );

c= @(i,I,m,e,g) max(alfa*z(g)*(K_grid(I)/G(K_grid(I),z(g)))^(alfa-1).*k_grid(i)+ ...
              (1-alfa)*z(g)*(K_grid(I)/G(K_grid(I),z(g)))^(alfa)*e*h_grid(m) +(1-delta)*k_grid(i) ...
             - k_grid,0) ;

% Starting values for V     
v1g = @(k,K,h) log( alfa*z(1)*(K/G(K,1))^(alfa-1).*k + (1-alfa)*z(1)*(K/G(K,1))^(alfa).*h -delta.*k )./(1-betta);
v1b = @(k,K,h) log( alfa*z(2)*(K/G(K,0))^(alfa-1).*k  + (1-alfa)*z(2)*(K/G(K,0))^(alfa).*h -delta.*k  )./(1-betta);
v0g = @(k,K,h) log( alfa*z(1)*(K/G(K,1))^(alfa-1).*k -delta.*k )./(1-betta);
v0b = @(k,K,h) log( alfa*z(2)*(K/G(K,0))^(alfa-1).*k -delta.*k )./(1-betta);
%____________________________________________________________
%
% Grid for k and K 
%____________________________________________________________

for I=1:size(K_grid,2)
    for m =1:size(h_grid,2)
        V1g(:,I,m)= v1g(k_grid,K_grid(I),h_grid(m))';
        V1b(:,I,m)= v1b(k_grid,K_grid(I),h_grid(m))';
        V0g(:,I,m)= v0g(k_grid,K_grid(I),h_grid(m))';
        V0b(:,I,m)= v0b(k_grid,K_grid(I),h_grid(m))';
    end
end
%% 2-1 with straight guess version 
tic

toc
% running: 33 mins
% dev: 9.7382e-04
%%   


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
       
for iter=1:500  % VFI  
    disp('staring generating')
for i=1:size(k_grid,2)
    for I=1:size(K_grid,2)        
       for m=1:size(h_grid,2)
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),1))); 
        V0gt(i,I,m)= max(log(c(i,I,m,0,1))' + betta * ([pize(1,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');
        V1gt(i,I,m)= max(log(c(i,I,m,1,1))' + betta * ([pize(2,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');  
       
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),0)));
        V0bt(i,I,m)= max(log(c(i,I,m,0,2))' + betta * ([pize(3,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');
        V1bt(i,I,m)= max(log(c(i,I,m,1,2))' + betta * ([pize(4,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');      
       end 
    end
end
disp('deviation: ')
dev = max(max(max(abs( [V0gt-V0g,V1gt-V1g,V0bt-V0b,V1bt-V1b]))))

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
        for m =1:size(h_grid,2)
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),1))); 
        [V0gt(i,I,m),a(i,I,m,2,1)]= max(log(c(i,I,m,0,1))' + betta * ([pize(1,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');
        [V1gt(i,I,m),a(i,I,m,1,1)]= max(log(c(i,I,m,1,1))' + betta * ([pize(2,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');  
       
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),0)));
        [V0bt(i,I,m),a(i,I,m,2,2)]= max(log(c(i,I,m,0,2))' + betta * ([pize(3,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');
        [V1bt(i,I,m),a(i,I,m,1,2)]= max(log(c(i,I,m,1,2))' + betta * ([pize(4,:)]*([V0g(:,Ip,m),V1g(:,Ip,m),V0b(:,Ip,m),V1b(:,Ip,m)]'))');      
        end 
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
H_ind(1)=0.2;
for t=2:2000
for n=1:1000 
% Evolution of assets
    % optimal saving = a(k,K, employment shock, productivity shock)
    N_state(n,1,t)=a(N_state(n,1,t-1),K_ind(t-1),H_ind(t-1),N_state(n,2,t-1),zt(t-1));
% Evolution of the employment status     
    N_state(n,2,t)= 2-(rand>=pize(1 + zt(t-1)*2 - N_state(n,2,t-1),zt(t)*2-1)/piZ(zt(t-1),zt(t))); 
    % given the probability of moving from zt(t-1) to zt(t)
end
% Storage of the sequence of aggregate capital
[dev2, K_ind(t)]=min(abs(k_grid(round(mean(N_state(:,1,t))))-K_grid));
% by Firm's FOC to store sequences of prices
r0(t) =  alfa*zt(t-1)*(K_ind(t-1)/H_ind(t-1))^(alfa-1);
w0(t) = (1-alfa)*zt(t-1)*(K_ind(t-1)/H_ind(t-1))^(alfa);
end

% iter_b larger than 1
else % Question 3, for a update guess on the parameters on Approximation
      
for t=2:2000
for n=1:1000
    
% Evolution of assets: for H changes, a(k,K,e,epsilon) also changes:
    N_state(n,1,t)=a(N_state(n,1,t-1),K_ind(t-1),H_ind(t-1), N_state(n,2,t-1),zt(t-1));
    
   % keep the same evolution of employment
end

% Storage of the sequence of aggregate capital
[dev2, K_ind(t)]=min(abs(k_grid(round(mean(N_state(:,1,t))))-K_grid));
% by Firm's FOC to store sequences of prices
r0(t) =  alfa*zt(t-1)*(K_ind(t-1)/H_ind(t-1))^(alfa-1);
w0(t) = (1-alfa)*zt(t-1)*(K_ind(t-1)/H_ind(t-1))^(alfa);
end

end

Yg=log(K_grid(K_ind(i_zg(20:end)))');
Xg=[ones(size(i_zg(20:end),2),1),log(K_grid(K_ind(i_zg(20:end)-1))')] ;   
Bg=Xg\Yg
b0gp=Bg(1);
b1gp=Bg(2);

YYg=log(H_grid(K_ind(i_zg(20:end)))');
BBg=Xg\YYg
d0gp=BBg(1);
d1gp=BBg(2);

% regression for bad times (burning the first 20 periods of bad times.

Yb=log(K_grid(K_ind(i_zb(20:end)))');
Xb=[ones(size(i_zb(20:end),2),1),log(K_grid(K_ind(i_zb(20:end)-1))')]  ;  
Bb=Xb\Yb
b0bp=Bb(1);
b1bp=Bb(2);

YYb=log(H_grid(K_ind(i_zb(20:end)))');
BBb=Xb\YYb
d0bp=BBb(1);
d1bp=BBb(2);


dev_b = max(abs([b0g-b0gp b1g-b1gp b0b-b0bp b1b-b1bp]))
dev_d = max(abs([d0g-d0gp d1g-d1gp d0b-d0bp d1b-d1bp]))

pause(.5)
if (dev_b<=0.01)&&(dev_d<=0.01) % (d) Iterate until convergence on beta
    break
end

% updating rules for K'
b0g=0.1*b0gp+0.9*b0g; 
b1g=0.1*b1gp+0.9*b1g;
b0b=0.1*b0bp+0.9*b0b;
b1b=0.1*b1bp+0.9*b1b;
% updating rules for H'
d0g=0.1*d0gp+0.9*d0g; 
d1g=0.1*d1gp+0.9*d1g;
d0b=0.1*d0bp+0.9*d0b;
d1b=0.1*d1bp+0.9*d1b;
end