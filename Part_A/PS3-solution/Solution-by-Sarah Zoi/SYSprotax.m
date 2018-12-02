function [output]=SYSprotax(x,y0i,r,etaj,sigma,beta,theta,delta,T1,T2,eps,k,nu)
tau1=1-delta*(etaj*x(2))^(-theta);
tau2h=1-delta*( (etaj+eps)*x(3))^(-theta);
tau2l=1-delta*( (etaj-eps)*x(3))^(-theta);
output=[ beta*0.5*(1+r)*( ((1-tau2h)*(etaj+eps)*x(3)+(1+r)*x(1)+T2  )^(-sigma) +...
                       +( (1-tau2l)*(etaj-eps)*x(3)+(1+r)*x(1)+T2  )^(-sigma))-((1-tau1)*etaj*x(2)+y0i-x(1)+T1)^(-sigma); 
         %FOC wrt h
                     [(1-tau1)*etaj + etaj*x(2)*( -theta*delta*etaj*(etaj*x(2))^(-theta-1) ) ]*[(1-tau1)*etaj*x(2)+y0i-x(1)+T1]^(-sigma)-k*x(2)^(1/nu);                   
         %FOC wrt h'
         beta*0.5* ([ [(1-tau2h)*(etaj+eps)+ (etaj+eps)*x(3)*( -theta*delta*(etaj+eps)*((etaj+eps)*x(3))^(-theta-1) )]*[(1-tau2h)*(etaj+eps)*x(3)+(1+r)*x(1)+T2]^(-sigma)-k*x(3)^(1/nu)]+...
                  + [ [(1-tau2l)*(etaj-eps)+ (etaj-eps)*x(3)*( -theta*delta*(etaj-eps)*((etaj-eps)*x(3))^(-theta-1) )]*[(1-tau2l)*(etaj-eps)*x(3)+(1+r)*x(1)+T2]^(-sigma)-k*x(3)^(1/nu)] )];

end