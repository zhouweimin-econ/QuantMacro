clear;clf;close all;
cd '~/Desktop/PS5/code'  % in order to save png
% This file solves II.4.2 - Q3  (copied from other files, no comments)
% Weimin Zhou
%% sigma=2, CRRA
T=45;
sigma=2;

gamma=0;
prob=[(1+gamma)/2 (1-gamma)/2;
    (1-gamma)/2 (1+gamma)/2]; 

beta=1/1.06;

r=0.04;
R=1.04;   

sigma_y=0.1;
Y=[1-sigma_y 1+sigma_y];

util=@(x) (x.^(1-sigma)-1)/(1-sigma); 
v=util([R*[-(1-sigma_y)/R:0.01:3]'+Y(1) R*[-(1-sigma_y)/R:0.01:3]'+Y(2)]);
for i=1:T
    j=0;
    tempv=[];
    tempind=[];
    for a=[(1-sigma_y)*(R^(-1-i)-1)/r:0.01:3] 
        j=j+1;
        b=[(1-sigma_y)*(R^(-i)-1)/r:0.01:3]';  
        K=util([R*a-b+Y(1) R*a-b+Y(2)]) + beta*(prob*v')'; 
        [M1,I1]=max( K( 1:max(1,min(int16(((min(3,(R*a+Y(1)))-(1-sigma_y)*(R^(-i)-1)/r))/0.01),max(size(K)))), 1 ));
        tempv(j,1)=M1;
        tempind(j,1)=I1;
        [M2,I2]=max( K( 1:max(1,min(int16(((min(3,(R*a+Y(2)))-(1-sigma_y)*(R^(-i)-1)/r))/0.01),max(size(K)))), 2 ));
        tempv(j,2) =M2;
        tempind(j,2)=I2;
    end
    v=tempv;
    V{2*(T-i+1)-1}=tempv(:,1);
    V{2*(T-i+1)}=tempv(:,2);
    Index{2*(T-i+1)-1}=tempind(:,1);
    Index{2*(T-i+1)}=tempind(:,2);
end
%% simulated consumption profile
Y_state=[];
Y_realize=[];
Y_state(1)=1;
Y_realize(1)=Y(Y_state(1));
for k=1:T
    Y_state(k+1)=logical(rand(1)<prob(Y_state(k),1))+1; 
    Y_realize(k+1)=Y(Y_state(k+1));
end

c2=[];
a=2;
for t=0:T-1
    c2(1,t+1)=t;
    c2(2,t+1)=R*a+Y_realize(t+1)-((Index{2*t+Y_state(t+1)}(int16((a-(1-sigma_y)*(R^(t-T-1)-1)/r)/0.01+1))-1)*0.01+(1-sigma_y)*(R^(t-T)-1)/r);
    a=(Index{2*t+Y_state(t+1)}(int16((a-(1-sigma_y)*(R^(t-T-1)-1)/r)/0.01+1))-1)*0.01+(1-sigma_y)*(R^(t-T)-1)/r;
end
c2(1,T+1)=T;
c2(2,T+1)=R*a+Y_realize(T+1);
%% sigma=5, CRRA
sigma=5;
util=@(x) (x.^(1-sigma)-1)/(1-sigma); 
v=util([R*[-(1-sigma_y)/R:0.01:3]'+Y(1) R*[-(1-sigma_y)/R:0.01:3]'+Y(2)]);
for i=1:T  %corresponding to period 45-i
    j=0;
    tempv=[];
    tempind=[];
    for a=[(1-sigma_y)*(R^(-1-i)-1)/r:0.01:3] 
        j=j+1;
        b=[(1-sigma_y)*(R^(-i)-1)/r:0.01:3]'; 
        K=util([R*a-b+Y(1) R*a-b+Y(2)]) + beta*(prob*v')'; 
        
        [M1,I1]=max( K( 1:max(1,min(int16(((min(3,(R*a+Y(1)))-(1-sigma_y)*(R^(-i)-1)/r))/0.01),max(size(K)))), 1 ));
        tempv(j,1)=M1;
        tempind(j,1)=I1;
        
        [M2,I2]=max( K( 1:max(1,min(int16(((min(3,(R*a+Y(2)))-(1-sigma_y)*(R^(-i)-1)/r))/0.01),max(size(K)))), 2 ));
        tempv(j,2) =M2;
        tempind(j,2)=I2;
    end
    
    v=tempv;
    V{2*(T-i+1)-1}=tempv(:,1);
    V{2*(T-i+1)}=tempv(:,2);
    Index{2*(T-i+1)-1}=tempind(:,1);
    Index{2*(T-i+1)}=tempind(:,2);
end
%% simulated consumption profile
Y_state=[];
Y_realize=[];
Y_state(1)=1;
Y_realize(1)=Y(Y_state(1));
for k=1:T
    Y_state(k+1)=logical(rand(1)<prob(Y_state(k),1))+1;
    Y_realize(k+1)=Y(Y_state(k+1));
end

c5=[];
a=2;
for t=0:T-1
    c5(1,t+1)=t;
    c5(2,t+1)=R*a+Y_realize(t+1)-((Index{2*t+Y_state(t+1)}(int16((a-(1-sigma_y)*(R^(t-T-1)-1)/r)/0.01+1))-1)*0.01+(1-sigma_y)*(R^(t-T)-1)/r);
    a=(Index{2*t+Y_state(t+1)}(int16((a-(1-sigma_y)*(R^(t-T-1)-1)/r)/0.01+1))-1)*0.01+(1-sigma_y)*(R^(t-T)-1)/r;
end
c5(1,T+1)=T;
c5(2,T+1)=R*a+Y_realize(T+1);
%% sigma=20, CRRA
sigma=20;
util=@(x) (x.^(1-sigma)-1)/(1-sigma); 
v=util([R*[-(1-sigma_y)/R:0.01:3]'+Y(1) R*[-(1-sigma_y)/R:0.01:3]'+Y(2)]);
for i=1:T  
    j=0;
    tempv=[];
    tempind=[];
    for a=[(1-sigma_y)*(R^(-1-i)-1)/r:0.01:3] 
        j=j+1;
        b=[(1-sigma_y)*(R^(-i)-1)/r:0.01:3]';  
        K=util([R*a-b+Y(1) R*a-b+Y(2)]) + beta*(prob*v')'; 
        [M1,I1]=max( K( 1:max(1,min(int16(((min(3,(R*a+Y(1)))-(1-sigma_y)*(R^(-i)-1)/r))/0.01),max(size(K)))), 1 ));
        tempv(j,1)=M1;
        tempind(j,1)=I1;
        [M2,I2]=max( K( 1:max(1,min(int16(((min(3,(R*a+Y(2)))-(1-sigma_y)*(R^(-i)-1)/r))/0.01),max(size(K)))), 2 ));
        tempv(j,2) =M2;
        tempind(j,2)=I2;
    end
    v=tempv;
    V{2*(T-i+1)-1}=tempv(:,1);
    V{2*(T-i+1)}=tempv(:,2);
    Index{2*(T-i+1)-1}=tempind(:,1);
    Index{2*(T-i+1)}=tempind(:,2);
end
%% simulated consumption profile
Y_state=[];
Y_realize=[];
Y_state(1)=1;
Y_realize(1)=Y(Y_state(1));
for k=1:T
    Y_state(k+1)=logical(rand(1)<prob(Y_state(k),1))+1; 
    Y_realize(k+1)=Y(Y_state(k+1));
end

c20=[];
a=2;
for t=0:T-1
    c20(1,t+1)=t;
    c20(2,t+1)=R*a+Y_realize(t+1)-((Index{2*t+Y_state(t+1)}(int16((a-(1-sigma_y)*(R^(t-T-1)-1)/r)/0.01+1))-1)*0.01+(1-sigma_y)*(R^(t-T)-1)/r);
    a=(Index{2*t+Y_state(t+1)}(int16((a-(1-sigma_y)*(R^(t-T-1)-1)/r)/0.01+1))-1)*0.01+(1-sigma_y)*(R^(t-T)-1)/r;
end
c20(1,T+1)=T;
c20(2,T+1)=R*a+Y_realize(T+1);

%% quadratic utility case (only one is enough)
util=@(x) -.5*(x-100).^2;
v=util([R*[-(1-sigma_y)/R:0.01:3]'+Y(1) R*[-(1-sigma_y)/R:0.01:3]'+Y(2)]);
for i=1:T  
    j=0;
    tempv=[];
    tempind=[];
    for a=[(1-sigma_y)*(R^(-1-i)-1)/r:0.01:3] 
        j=j+1;
        b=[(1-sigma_y)*(R^(-i)-1)/r:0.01:3]';  
        K=util([R*a-b+Y(1) R*a-b+Y(2)]) + beta*(prob*v')'; 
        [M1,I1]=max( K( 1:max(1,min(int16(((min(3,(R*a+Y(1)))-(1-sigma_y)*(R^(-i)-1)/r))/0.01),max(size(K)))), 1 ));
        tempv(j,1)=M1;
        tempind(j,1)=I1;
        [M2,I2]=max( K( 1:max(1,min(int16(((min(3,(R*a+Y(2)))-(1-sigma_y)*(R^(-i)-1)/r))/0.01),max(size(K)))), 2 ));
        tempv(j,2) =M2;
        tempind(j,2)=I2;
    end
    v=tempv;
    V{2*(T-i+1)-1}=tempv(:,1);
    V{2*(T-i+1)}=tempv(:,2);
    Index{2*(T-i+1)-1}=tempind(:,1);
    Index{2*(T-i+1)}=tempind(:,2);
end
%% simulated consumption profile
Y_state=[];
Y_realize=[];
Y_state(1)=1;
Y_realize(1)=Y(Y_state(1));
for k=1:T
    Y_state(k+1)=logical(rand(1)<prob(Y_state(k),1))+1; 
    Y_realize(k+1)=Y(Y_state(k+1));
end

c_quad=[];
a=2;
for t=0:T-1
    c_quad(1,t+1)=t;
    c_quad(2,t+1)=R*a+Y_realize(t+1)-((Index{2*t+Y_state(t+1)}(int16((a-(1-sigma_y)*(R^(t-T-1)-1)/r)/0.01+1))-1)*0.01+(1-sigma_y)*(R^(t-T)-1)/r);
    a=(Index{2*t+Y_state(t+1)}(int16((a-(1-sigma_y)*(R^(t-T-1)-1)/r)/0.01+1))-1)*0.01+(1-sigma_y)*(R^(t-T)-1)/r;
end
c_quad(1,T+1)=T;
c_quad(2,T+1)=R*a+Y_realize(T+1);
%% 

figure
plot(c2(1,:),c2(2,:));
hold on
plot(c5(1,:),c5(2,:));
hold on
plot(c20(1,:),c20(2,:));
hold on

legend('\sigma=2','\sigma=5','\sigma=20')
title('consumption profile, T = 45, CRRA')
xlabel('time t') 
ylabel('c2(t)') 
saveas(gcf,'sigmacompare1.png')
figure
plot(c_quad(1,:),c_quad(2,:));
title('consumption profile, T = 45, Quadratic')
xlabel('time t') 
ylabel('c2(t)') 
saveas(gcf,'sigmacompare2.png')