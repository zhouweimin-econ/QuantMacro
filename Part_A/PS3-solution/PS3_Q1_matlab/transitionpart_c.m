clear

%  structural parameters
z=1.62968;
h=0.31;
A=1.0;             	%   Technology level
beta= 0.980392157;      	   %   Discout factor
sigma=0.999;      	   %   Relative risk aversion
n=0.0;            	%   Population rate
delta=0.25/4;        	%   Depreciation rate
alpha=0.33;     	   %   Capital share? 1- theta
tauc=0.6;          	%   Consumption tax rate
tauy=0.0;      	%   Income tax rate


% steady state

kss0=z*h*(alpha*A*(1-tauy)/(((1+n)/beta)-(1-delta)))^(1/(1-alpha)); % capital stock
yss0=A*kss0^alpha*(z*h)^(1-alpha); % output
css0=(kss0/(1+tauc))*((1-tauy)*(yss0/kss0)-(n+delta)); % consumption
iss0=(n+delta)*kss0; % investment
uss0=(css0^(1-sigma)-1)/(1-sigma); % single-period utility
revss0=tauy*yss0+tauc*css0; % public revenues
revyss0=revss0/yss0; % public revenues as % output
ydss0=(1-tauy)*yss0; % disposable income

    % other way to calculate the steady state
    param=[A beta sigma n delta alpha tauc tauy z h];
    x0=[kss0 kss0 css0 css0]';
    x=fsolve('ss_ck_d',x0,optimset,param);
    kss=x(1);
    css=x(3);
    yss=A*kss0^alpha*(z*h)^(1-alpha); 
    iss=(n+delta)*kss; 
    uss=(css^(1-sigma)-1)/(1-sigma);      
    revss=tauy*yss+tauc*css; 
    revyss=revss/yss; 
    ydss=(1-tauy)*yss;

    %what agents believe after change :
z = z*2;

yss=A*kss^alpha*(z*h)^(1-alpha); % I didn't normalize output to 1 anymore. 
kss=z*h*(alpha*A*(1-tauy)/(((1+n)/beta)-(1-delta)))^(1/(1-alpha)); % capital stock

css=(kss/(1+tauc))*((1-tauy)*(yss/kss)-(n+delta)); % consumption
iss=(n+delta)*kss; % investment
uss=(css^(1-sigma)-1)/(1-sigma); % single-period utility
revss=tauy*yss+tauc*css; % public revenues
revyss=revss0/yss; % public revenues as % output
ydss=(1-tauy)*yss; % disposable income    
    
    J=jacobian('ss_ck_d',x,param); % transition matrix
    MA=[J(1,1) J(1,3);
        J(2,1) J(2,3)];
    MB=[J(1,2) J(1,4);
        J(2,2) J(2,4)];
    MG=-inv(MA)*MB;
    [Mvec,Meig]=eig(MG);
    Meig2=diag(Meig);
    [Meig3,Mord]=sort(abs(Meig2));
    Meig4=Meig2(Mord);
    Mlambda=diag(Meig4);
    MM=Mvec(:,Mord);
    Mm=inv(MM);
    
    stab=-Mm(2,1)/Mm(2,2); % stabilizing constant

% We choose an initial condition for the stock of capital
    
    k00=4 *0.99; 
    
% solution using the model (imposing stability just at the initial period)
    nobs=50;
   
% solution using the exponential convergence (imposing stability just at the initial period)
    k0=k00;
    c0 = css0;
    %c0=css+stab*(k0-kss);
    vk3=[];vk3=[vk3;k0];
    vc3=[];vc3=[vc3;c0];
    for i=1:10
        k1=kss+(Mlambda(1,1)^i)*(k0-kss);
        c1=css+(Mlambda(1,1)^i)*(c0-css);
        vk3=[vk3;k1];
        vc3=[vc3;c1];
    end
 
    for i=10:nobs
        k2=kss0+(Mlambda(1,1)^i)*(k1-kss0);
        c2=css0+(Mlambda(1,1)^i)*(c1-css0);
        vk3=[vk3;k2];
        vc3=[vc3;c2];
    end  
    


% solution using the linear aproximation to the model (imposing stability at all points in time)
    k0=k00;
    c0=css+stab*(k0-kss);
    vk4=[];vk4=[vk4;k0];
    vc4=[];vc4=[vc4;c0];
    for i=1:nobs
        k1=kss+MG(1,1)*(k0-kss)+MG(1,2)*(c0-css);
        c1=css+stab*(k0-kss);
        c0=c1;
        k0=k1;
        vk4=[vk4;k0];
        vc4=[vc4;c0];
    end


figure;
hold on

subplot(2,2,1)
plot(vk3);
title('solution using the exponential convergence (imposing stability just at the initial period)');
ylabel('capital stock 2');
xlabel('periods');
%Axis('tight');

subplot(2,2,2)
plot(vk4);
title('solution using the linear approximation to the model (imposing stability at all points in time)');
ylabel('capital stock 3');
xlabel('periods');  
%Axis('tight');

subplot(2,2,3)
plot(vc3);
title('consumption');
ylabel('Consumption');
xlabel('periods');
% Axis('tight');

subplot(2,2,4)
plot(vc4);
title('consumption2');
ylabel('Consumption');
xlabel('periods');
% Axis('tight');



hold off



