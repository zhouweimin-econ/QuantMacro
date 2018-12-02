function [output]=SYS(x,y0i,r,etaj,sigma,beta,tau,T1,T2,eps,k,nu)

output=[ beta*0.5*(1+r)*( ((1-tau)*(etaj+eps)*x(3)+(1+r)*x(1)+T2 )^(-sigma) +...
                       +( (1-tau)*(etaj-eps)*x(3)+(1+r)*x(1)+T2  )^(-sigma))-((1-tau)*etaj*x(2)+y0i-x(1)+T1 )^(-sigma);                    
         (1-tau)*etaj*[(1-tau)*etaj*x(2)+y0i-x(1)+T1 ]^(-sigma)-k*x(2)^(1/nu);                   
         beta* 0.5* ([(1-tau)*(etaj+eps)*[(1-tau)*(etaj+eps)*x(3)+(1+r)*x(1)+T2]^(-sigma)-k*x(3)^(1/nu)]+...
              + [(1-tau)*(etaj-eps)*[(1-tau)*(etaj-eps)*x(3)+(1+r)*x(1)+T2 ]^(-sigma)-k*x(3)^(1/nu)])   ];

end