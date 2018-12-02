clear;clf;close all;
cd '~/Desktop/PS5/code'  % in order to save png
% This file solves certainty case with CRRA
% Euler method is used since backwards VFI is not as clean as it.
% Weimin Zhou
%% CRRA, finite, certainty
global a_0  beta T R;

beta=1/1.06;
R=1.04;
T=45;

% initial
a_0=1; 
c_0=0;

% solve c_1
c_1 =[];
c_1(1)= fsolve(@tran,c_0);

for j=2:T+1
    c_1(j)=c_1(j-1)*(beta*R)^.5; % impose euler method: c(t)
end

figure
plot(0:45,c_1)  % T = 45
title('Consumption profile for CRRA utility(initial wealth a_0=1)')
saveas(gcf,'41b.png')

c_0=[];
iter=1;
T = 40;
for a_0=linspace(0,2)  % for a given period, compute policy function w.r.t state a
c_00 = 0;
c_0(iter) = fsolve(@tran,c_00);
iter=iter+1;
end
a_0=linspace(0,2);

figure
plot(a_0,c_0) % T = 40
xlabel('$a$');
ylabel('$c_{40} (a,y)$')
saveas(gcf,'41c.png')

c_0=[];
iter=1;
T = 5;
for a_0=linspace(0,2)  % for a given period, compute policy function w.r.t state a
c_00 = 0;
c_0(iter) = fsolve(@tran,c_00);
iter=iter+1;
end

a_0=linspace(0,2);
figure
plot(a_0,c_0) 
xlabel('$a$');
ylabel('$c_{5} (a,y)$')
saveas(gcf,'41d.png')

