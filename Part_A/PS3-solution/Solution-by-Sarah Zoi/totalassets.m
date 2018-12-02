function [atot]=totalassets(y0,r,eta,sigma,beta,tau,T1,T2,eps,k,nu)
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
%  syms a h hp
% xstar=solve( beta*0.5*(1+r)*( ((1-tau)*(etaj+eps)*hp+(1+r)*a+T  )^(-sigma)+...
%             ( (1-tau)*(etaj-eps)*hp+(1+r)*a+T  )^(-sigma))-((1-tau)*etaj*h+y0i-a)^(-sigma)==0,...          
%             (1-tau)*etaj*[(1-tau)*etaj*h+y0i-a]^(-sigma)-k*h^(1/nu)==0,...   
%              beta* 0.5* ([(1-tau)*(etaj+eps)*[(1-tau)*(etaj+eps)*hp+(1+r)*a+T]^(-sigma)-k*hp^(1/nu)]+...
%              [(1-tau)*(etaj-eps)*[(1-tau)*(etaj-eps)*hp+(1+r)*a+T]^(-sigma)-k*hp^(1/nu)])==0,...
%               h<1,hp<1,a,h,hp )



% xstar=solve(SYS(xstart,y0i,r,etaj,sigma,beta,tau,T,eps,k,nu),x);
 optch(i,:,j)=xstar; 
 
 % check if the borrowing limit is hit
 bl= -((1-tau)*(eta(j)-eps)*xstar(3)+T2)/(1+r);
         if optch(i,1,j)<= bl
             optch(i,1,j)=bl;
             % get h and hp as a solution of a system of 2 equations in 2
             % unknowns
             xstar=fsolve( @(h) SYSbl(h,y0i,r,etaj,sigma,beta,tau,T1,T2,eps,k,nu,bl),xstart(2:3));
             xstar=[bl xstar];
             bln=bln+1;
         else
         end
 
   end   
end
% check aggregate supply and demand for assets
atot=sum(sum(optch(:,1,:)));

end