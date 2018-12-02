function [SRagg, Lshare, H]=aggstat(y0,optch,eta,tau,r,T1,T2)
n=length(y0);
wh0=nan(n,2,length(eta));
wh1=nan(n,length(eta));

if nargin==5
    T1=0;
    T2=0;
end
for i=1:length(eta)
    
    if length(tau)==2
    delta=tau(1);
    theta=tau(2);
   
tau1=ones(n,1)-delta*( eta(i)*optch(:,2,i) ).^(-theta);
tau2=ones(n,1)-delta*( eta(i)*optch(:,3,i)).^(-theta);
    else
    tau1=tau;
    tau2=tau;
end
% pretax labor income t=0
wh0(:,1,i)=optch(:,2,i)*eta(i);
% total pre-tax income t=0
wh0(:,2,i)=wh0(:,1,i)+y0'+T1;
% after-tax labor income t=0 
wh0(:,3,i)=optch(:,2,i)*eta(i).*(ones(length(tau1),1)-tau1);
% total after-tax income t=0
wh0(:,4,i)=wh0(:,3,i)+y0'+T1;

% pretax labor income t=1 
wh1(:,1,i)=optch(:,3,i)*(eta(i));
% pretax income t=1 
wh1(:,2,i)=optch(:,3,i)*(eta(i))+optch(:,1,i)*(1+r)+T2;
% after-tax labor income t=1 
wh1(:,3,i)=optch(:,3,i)*(eta(i)).*(ones(length(tau2),1)-tau2);
% total after-tax income t=1 
wh1(:,4,i)=wh1(:,3,i)+optch(:,1,i)*(1+r)+T2;
end
%pretax labor income share t=0
Lshare(1)=sum(sum(wh0(:,1,:)))/sum(sum(wh0(:,2,:)));
%aftertax labor income share t=0
Lshare(2)=sum(sum(wh0(:,3,:)))/sum(sum(wh0(:,4,:)));
%pretax labor income share t=1
Lshare(3)=sum(sum(wh1(:,1,:)))/sum(sum(wh1(:,2,:)));
%aftertax labor income share t=1
Lshare(4)=sum(sum(wh1(:,3,:)))/sum(sum(wh1(:,4,:)));

SRagg=sum(optch(optch(:,1,:)>=0))/(sum(y0)*length(eta));
H(:,1)=sum(sum(optch(:,2,:)))/(length(eta)*n);
H(:,2)=sum(sum(optch(:,3,:)))/(length(eta)*n);
end