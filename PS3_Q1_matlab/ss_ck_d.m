function f=ss_ck_d(x0,param)

A=param(1);            	%   Technology level
beta=param(2);      	%   Discout factor
sigma=param(3);      	%   Relative risk aversion
n=param(4);            	%   population rate
delta=param(5);        	%   depreciation rate
alpha=param(6);     	%   capital share
tauc=param(7);      	%   consumption tax rate
tauy=param(8);         	%   income tax rate
z = param(9);
h = param(10);

k1=x0(1);
k=x0(2);
c1=x0(3);
c=x0(4);

f(1)=(1-tauy)*A*(k^alpha)*(z*h)^(1-alpha)+(1-delta)*k-(1+tauc)*c-(1+n)*k1;
f(2)=-c1+c*((beta/(1+n))*((1-tauy)*A*alpha*(k1^(alpha-1)*(z*h)^(1-alpha))...
    +(1-delta)))^(1/sigma);
f(3)=k1-k;
f(4)=c1-c;
f=f';