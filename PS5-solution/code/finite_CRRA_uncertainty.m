clear;clf;close all;
cd '~/Desktop/PS5/code'  % in order to save png
% This file solves uncertainty case with CRRA
% Weimin Zhou
%% CRRA, finite, uncertainty
T=45;
gamma=0.95;
prob=[(1+gamma)/2, (1-gamma)/2;
    (1-gamma)/2, (1+gamma)/2]; %transition matrix for income

beta=1/1.06; %discounting factor
r=0.04;
R=1.04;      %interest rate
sigma_y=0.5;
sigma=2;
Y=[1-sigma_y 1+sigma_y];
util=@(x) (x.^(1-sigma)-1)/(1-sigma); %CRRA utility function
v=util([R*[-(1-sigma_y)/R:0.01:3]'+Y(1) R*[-(1-sigma_y)/R:0.01:3]'+Y(2)]);%this is v45(a,y_low,y_high)
for i=1:T  %corresponding to period 45-i
    j=0;
    tempv=[];
    tempind=[];
    for a=[(1-sigma_y)*(R^(-1-i)-1)/(R-1):0.01:3] %state variable: initial wealth "a"
        j=j+1;
        b=[(1-sigma_y)*(R^(-i)-1)/(R-1):0.01:3]';  % possible choice(policy function a'(a)) for T-i+1
        K=util([R*a-b+Y(1) R*a-b+Y(2)]) + beta*(prob*v')'; % return the maximum of all the choice
        
        [M1,I1]=max( K( 1:max(1,min(int16(((min(3,(R*a+Y(1)))-(1-sigma_y)*(R^(-i)-1)/(R-1)))/0.01),max(size(K)))), 1 ));
        tempv(j,1)=M1;
        tempind(j,1)=I1;
        
        [M2,I2]=max( K( 1:max(1,min(int16(((min(3,(R*a+Y(2)))-(1-sigma_y)*(R^(-i)-1)/(R-1)))/0.01),max(size(K)))), 2 ));
        tempv(j,2) =M2;
        tempind(j,2)=I2;
    end
    v=tempv;
    V{2*(T-i+1)-1}=tempv(:,1);
    V{2*(T-i+1)}=tempv(:,2);
    Index{2*(T-i+1)-1}=tempind(:,1);
    Index{2*(T-i+1)}=tempind(:,2);
end
%% policy function for consumption

a_5 = [(1-sigma_y)*(R^(-T-1)-1)/(R-1):0.01:3]';
cl_5 = R*[(1-sigma_y)*(R^(-T-1)-1)/(R-1):0.01:3]'+Y(1)-(.01*tempind(:,1)+(1-sigma_y)*(R^(-T)-1)/(R-1));
ch_5 = R*[(1-sigma_y)*(R^(-T-1)-1)/(R-1):0.01:3]'+Y(2)-(.01*tempind(:,2)+(1-sigma_y)*(R^(-T)-1)/(R-1));
%% Q4. compare sigma_y from .1 to .5
%{
figure
subplot(1,2,1)
plot(a_1, ch_1, a_5, ch_5);
title('cunsumption policy - CRRA (good shock)')
legend('\sigma_y=0.1','\sigma_y=0.5')
xlabel('a') 
ylabel('c(a,y)') 
subplot(1,2,2)
plot(a_1, cl_1, a_5, cl_5);
title('cunsumption policy - CRRA (bad shock)')
legend('\sigma_y=0.1','\sigma_y=0.5')
xlabel('a') 
ylabel('c(a,y)') 
saveas(gcf,'423-c.png')
%}
%% Q5. compare sigma from 0 to .95
%
figure
subplot(1,2,1)
plot(a_1, ch_1, a_5, ch_5);
title('cunsumption policy - CRRA (good shock)')
legend('\gamma=0','\gamma=0.95')
xlabel('a') 
ylabel('c(a,y)') 
subplot(1,2,2)
plot(a_1, cl_1, a_5, cl_5);
title('cunsumption policy - CRRA (bad shock)')
legend('\gamma=0','\gamma=0.95')
xlabel('a') 
ylabel('c(a,y)') 
saveas(gcf,'425.png')
%}
%%
figure
plot(a_5, ch_5, a_5, cl_5);
title('cunsumption policy - CRRA')
legend('1-\sigma_y','1+\sigma_y')
xlabel('a') 
ylabel('c_5 (a,y)') 
%saveas(gcf,'41fa.png')
%% simulated consumption profile
Y_state=[];
Y_realize=[];
Y_state(1)=1;
Y_realize(1)=Y(Y_state(1));
for k=1:T
    Y_state(k+1)=logical(rand(1)<prob(Y_state(k),1))+1; % a Markov income process for y(t)
    Y_realize(k+1)=Y(Y_state(k+1));
end

c=[];
a=2;
for t=0:T-1
    c(1,t+1)=t;
    c(2,t+1)=R*a+Y_realize(t+1)-((Index{2*t+Y_state(t+1)}(int16((a-(1-sigma_y)*(R^(t-T-1)-1)/(R-1))/0.01+1))-1)*0.01+(1-sigma_y)*(R^(t-T)-1)/(R-1));
    a=(Index{2*t+Y_state(t+1)}(int16((a-(1-sigma_y)*(R^(t-T-1)-1)/(R-1))/0.01+1))-1)*0.01+(1-sigma_y)*(R^(t-T)-1)/(R-1);
end
c(1,T+1)=T;
c(2,T+1)=R*a+Y_realize(T+1);

figure
plot(c(1,:),c(2,:))
title('Cunsumption profile for 45 periods-CRRA utility')
xlabel('periods t') 
ylabel('c(t)') 

