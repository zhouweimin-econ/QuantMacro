clear;clf;close all;
cd '~/Desktop/PS5/code'  % in order to save png
% This file solves certainty case with Quadratic
% Weimin Zhou
%% Quadratic, finite, certainty
global a  beta T R;
beta=1/1.06;
R=1.04;
T=5;

y=1;
cbar = 100;

util=@(x) -0.5*(x-cbar).^2; % quad
v=util(R*[-1/R:0.01:3]'+y);

for i=1:T  
    j=0;
    tempv=[];
    
    for a=[(R^(-1-i)-1)/(R-1):0.01:3] %state a
        j=j+1;
        aprime=[(R^(-i)-1)/(R-1):0.01:3]';  % possible choice(policy function a'(a)) for T-i+1
        v0=util(R*a+y-aprime)+beta*v;
        decis = min(int16((min(3,R*a+y)-((1/R)^(i)-1)/(R-1))/0.01), max(size(v0)));
        [M,N]=max(v0(1:max(1,decis))) ;
        tempv(j,1)=M;
        tempv(j,2)=N;
    end
    
    v=tempv(:,1);
    V{T-i+1}=tempv(:,2); % Since backwards 
end

%% consumption profile
c=[];
a=2;
for t=0:T-1
    c(1,t+1)=t;
    c(2,t+1)=R*a+y-((V{t+1}(int16((a-(R^(t-T-1)-1)/(R-1))/0.01+1))-1)*0.01+(R^(t-T)-1)/(R-1));
    a=(V{t+1}(int16((a-(R^(t-T-1)-1)/(R-1))/0.01+1))-1)*0.01+(R^(t-T)-1)/(R-1);
end
c(1,T+1)=T;
c(2,T+1)=R*a+y;

figure
plot(c(1,:),c(2,:))
xlabel('time t') 
ylabel('consumption time profile c(t)') 
saveas(gcf,'42a.png')
%% similar approach 
a_5 = [(R^(-T-1)-1)/(R-1):0.01:3]';
c_5 = R*[(R^(-T-1)-1)/(R-1):0.01:3]'+y-(.01*tempv(:,2)+(R^(-T)-1)/(R-1));
%%
figure
subplot(1,2,1)
plot(a_40,c_40);
xlabel('a') 
ylabel('c_40 (a,y)') 
subplot(1,2,2)
plot(a_5,c_5);
xlabel('a') 
ylabel('c_5 (a,y)') 
saveas(gcf,'42b.png')
