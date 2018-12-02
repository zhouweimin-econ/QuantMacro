function [output]=GE(y0,r,eta,sigma,beta,tau,T1,T2,eps,k,nu,etarealized)
global rev
global optch
n=length(y0);
 bln=0;
options=optimoptions('fsolve','MaxFunEvals',300000,'MaxIter',100000);
for j=1:length(eta)
      
   for i=1:length(y0)       
 etaj=eta(j);
 y0i=y0(i);     
 if i==1
     xstart=[0.01 0.5 0.5];
 else
     xstart=[xstar(1)+0.005 xstar(2)-0.03 xstar(3)-0.03];
 end
xstar=fsolve( @(x) SYS(x,y0i,r,etaj,sigma,beta,tau,T1,T2,eps,k,nu),xstart,options);
optch(i,:,j)=xstar; 
 
 % check if the borrowing limit is hit
 bl= -((1-tau)*(eta(j)-eps)*xstar(3)+T2)/(1+r);
         if optch(i,1,j)<= bl
             optch(i,1,j)=bl;
             % get h and hp as a solution of a system of 2 equations in 2
             % unknowns
             xstar=fsolve( @(h) SYSbl(h,y0i,r,etaj,sigma,beta,tau,T1,T2,eps,k,nu,bl),xstart(2:3));
             xstar=[bl xstar];
             bln=bln+1
         else
         end

rev(i,1,j)=xstar(2)*etaj*tau;
rev(i,2,j)=xstar(3)*etarealized(i,j)*tau;


   end   
end

% check aggregate supply and demand for assets
atot=sum(sum(optch(:,1,:)));
% government budget balance (T is distributed for 2 periods)
govbal1=sum(sum(rev(:,1,:)))-T1*n*length(eta);
govbal2=sum(sum(rev(:,2,:)))-T2*n*length(eta);
%govbal=sum(sum(rev(:,1,:)))-T*n*length(eta)+(sum(sum(rev(:,2,:)))-T*n*length(eta))/(1+r);
output=[atot;govbal1;govbal2];
%govbal];

end

  