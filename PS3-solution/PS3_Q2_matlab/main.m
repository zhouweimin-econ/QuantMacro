y00 = 0.001 + (0.009-0.001)*rand(100,1);
for i=1:length(y00)
    if y00(i)<=0.0087 
        if y00(i)>=0.0055
            y00(i) = 0.001;
        end
    end
end
y00 = sort(y00);

eta1 = [1, 1.5, 2.5, 3];
options = optimset('display','Iter','TolX',1e-6,'MaxIter',50);
a = [];
c = [];
hph = [];
hpl =[];
h = [];
cph = [];
cpl = [];

R= 1.0057;  %Guess R*beta < 1 means R < 1.0101
x0 = [0.5;.5;.5;.5;0.5; 0.5;0.5]; % Initial Guess
%% 
while  R>1.00 && R<1.30

    for i=1:length(eta1)
        for j=1:length(y00)
        myfun1   = @(x) root1(x,y00(j),eta1(i),R);
        results1 = fsolve(myfun1,x0,options);
        C(i,j)   = results1(1,1);
        Hph(i,j) = results1(2,1);
        Hpl(i,j) = results1(3,1);
        H(i,j)   = results1(4,1);
        A(i,j)   = results1(5,1);
        Cph(i,j) = results1(6,1);
        Cpl(i,j) = results1(7,1);
        end  
    end
    
    if sum(A(1,:))+sum(A(2,:))+sum(A(3,:))+sum(A(4,:))>1   % all are lending, R is too high
        R = R - 0.001; 
        continue;
    
    elseif sum(A(1,:))+sum(A(2,:))+sum(A(3,:))+sum(A(4,:))<-1
        R = R + 0.001;
        continue;
    
    else
        break;
    
    end
    
end

% Above codes takes 1 mins in my laptop, for a small dimension of
% interest rate, takes me half hours, but does not clears the assets
% markets.
%%  Below are the failed part of try, which takes me 2 days, and I don't want to delete them
% use focs.m, which remain the origin focs, unlike root.m changes the expression. 
%R = 1.0057;
%x0 = [0.5;.5;.5;.5;0.5];
%eta1 = [1, 1.5, 2.5, 3];
%options = optimset('display','Iter', 'TolX',1e-6,'MaxIter',30);
%A=zeros(length(eta1), length(y00));
%
%
%for i=1:length(eta1)
%        
%    for j=1:length(y00)
%        myfun1   = @(x) focs(x,y00(j),eta1(i),R);
%        results1 = fsolve(myfun1,x0,options);
%        A(i, j)  = results1(2,1);
%        C(i, j)  = results1(1,1);
%        Hph(i, j) = results1(4,1);
%        Hpl(i,j) = results1(5,1);
%        H(i, j)  = results1(3,1);
%    end  
%end

%for i=1:length(eta1)
%   CH = (eta1(i)+.05)*Hph + R*A;
%    CL = (eta1(i)-.05)*Hpl + R*A;
%end
%% Figure .1 
subplot(2,2,1)
plot(y00,A)
title('Optimal Savings a ');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');

subplot(2,2,2)
plot(y00,C)
title('Optimal Consumption c ');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');

subplot(2,2,3)
plot(y00,Cph)
title('Optimal Consumption c high ');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');

subplot(2,2,4)
plot(y00,Cpl)
title('Optimal Consumption c low ');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');
%%  Figure .2
Savings = A ./(C+A);
plot(y00, Savings)
title('Savings rate');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');
%% Figure .3
hold on 
subplot(2,2,1)
plot(y00,H)
title('hours worked today h');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');

subplot(2,2,2)
plot(y00,Hph,'-.','LineWidth',2)
title('hours worked tomorrow h with good-luck shock');
legend('tmr high \eta_y =1','tmr-high-1.5','tmr-high-2.5','tmr-high-3');

subplot(2,2,3)
plot(y00,Hpl,'-.','LineWidth',2)
title('hours worked tomorrow h with back-luck shock');
legend('tmr-low-1','tmr-low-1.5','tmr-low-2.5','tmr-low-3');
hold off
%% Figure .4
y0 = [y00,y00,y00,y00];
subplot(3,2,1)
hold on
for i=1:1:4
    plot(y00,eta1(i)*H(i,:))
end
title('before-tax labor income today wh');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');
hold on

subplot(3,2,2)
hold on
for i=1:1:4
    plot(y00,eta1(i)*H(i,:)./(eta1(i)*H(i,:)+y0(i)'))
end
title('after-tax labor share (1-tau)wh/((1-tau)wh + y0+T1)');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');
hold off

subplot(3,2,3)
hold on
for i=1:1:4
    plot(y00,(eta1(i)+.05)*Hph(i,:))
end
title('optimal labor income tomorrow (wh)^{prime} for good luck');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');
hold off

subplot(3,2,4)
hold on
for i=1:1:4
    plot(y00,(eta1(i)-.05)*Hpl(i,:))
end
title('optimal labor income tomorrow (wh)^{prime} for back luck');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');
hold off


subplot(3,2,5)
hold on 
for i=1:1:4
    plot(y00,(eta1(i)+.05)*Hph(i,:)./(((eta1(i)+.05)*Hph(i,:))+R*A(i)))
end
title('good-luck after-tax labor share of tomorrow (1-tau)wh/((1-tau)wh+(1+r)a+T2)');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');
hold off

subplot(3,2,6)
hold on 
for i=1:1:4
    plot(y00,(eta1(i)-.05)*Hpl(i,:)./(((eta1(i)-.05)*Hpl(i,:))+R*A(i)))
end
title('bad-luck after-tax labor share of tomorrow (1-tau)wh/((1-tau)wh+(1+r)a+T2)');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');
hold off
%% Figure 5.
gc1 = (Cph-C)./C;
gc2 = (Cpl-C)./C;
Egc = .5*gc1 + .5*gc2;
subplot(3,2,1)
plot(y00,gc1, y00, gc2,'-.')
title('consumption growth ');
legend('gc1 \eta_y =1','gc1,1.5','gc1, 2.5','gc1, 3', 'gc2 \eta_y =1','gc2,1.5','gc2, 2.5','gc2, 3');

subplot(3,2,2)
plot(y00,Egc)
title('Expected consumption growth');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');

for i = 1:1:4
    gwh1(i,:) = ((eta1(i)+.05).*Hph(i,:) - eta1(i).*H(i,:))./ eta1(i).*H(i,:);
    gwh2(i,:) = ((eta1(i)-.05).*Hpl(i,:) - eta1(i).*H(i,:))./ eta1(i).*H(i,:);
end

Egwh = .5*gwh1 + .5*gwh2;

subplot(3,2,3)
plot(y00,gwh1, y00, gwh2,'-.')
title('labor income growth');
legend('gwh1 \eta_y =1','gwh1,1.5','gwh1, 2.5','gwh1, 3', 'gwh2 \eta_y =1','gwh2,1.5','gwh2, 2.5','gwh2, 3');


subplot(3,2,4)
plot(y00,Egwh)
title('Expected labor income growth');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');


subplot(3,2,5)
plot(y00,Egc./Egwh)
title('Expected labor income growth');
legend('\eta_y =1','\eta_y =1.5','\eta_y =2.5','\eta_y =3');

elasticity1 = (gc1./gwh1)./(Egc./Egwh);
elasticity2 = (gc2./gwh2)./(Egc./Egwh);
subplot(3,2,6)
plot(y00,elasticity1,y00, elasticity2)
title('Expected labor income growth');
legend('e1 \eta_y =1','e1,1.5','e1, 2.5','e1, 3', 'e2 \eta_y =1','e2,1.5','e2, 2.5','e2, 3');


%% Figure 6 Capital market clearing


%% Figure 7 Wealfare
nu = 4;
sigma = 3;  % hence, V is negative
kappa = 4; 
beta = 0.99;
U1 = C.^(1-sigma)/(1-sigma) - kappa* H.^(1+1/nu)/(1+1/nu);
U2 = beta* (.5*(Cph.^(1-sigma)/(1-sigma) - kappa* Hph.^(1+1/nu)/(1+1/nu)) + (Cpl.^(1-sigma)/(1-sigma) - kappa* Hpl.^(1+1/nu)/(1+1/nu)));
V =U1 + U2;  
plot(y00,V) 
title('Welfare V as a function of initial wealth');