function [output]=SYSbl(h,y0i,r,etaj,sigma,beta,tau,T1,T2,eps,k,nu,bl)

output=[                  
         (1-tau)*etaj*[(1-tau)*etaj*h(1)+y0i-bl+T1]^(-sigma)-k*h(1)^(1/nu);                   
         beta*( 0.5*[(1-tau)*(etaj+eps)*[(1-tau)*(etaj+eps)*h(2)+(1+r)*bl+T2]^(-sigma)-k*h(2)^(1/nu)]+0.5*[(1-tau)*(etaj-eps)*[(1-tau)*(etaj-eps)*h(2)+(1+r)*bl+T2]^(-sigma)-k*h(2)^(1/nu)])   ];

end