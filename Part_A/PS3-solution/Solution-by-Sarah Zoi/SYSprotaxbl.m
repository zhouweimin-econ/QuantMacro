function [output]=SYSprotaxbl(x,y0i,r,etaj,sigma,beta,theta,delta,T1,T2,eps,k,nu,bl)
tau1=1-delta*(etaj*x(1))^(-theta);
tau2h=1-delta*( (etaj+eps)*x(2))^(-theta);
tau2l=1-delta*( (etaj-eps)*x(2))^(-theta);
output=[ %FOC wrt h
                     [(1-tau1)*etaj + etaj*x(1)*( -theta*delta*etaj*(etaj*x(1))^(-theta-1) ) ]*[(1-tau1)*etaj*x(1)+y0i-bl+T1]^(-sigma)-k*x(1)^(1/nu);                   
         %FOC wrt h'
         beta*0.5* ([ [(1-tau2h)*(etaj+eps)+ (etaj+eps)*x(2)*( -theta*delta*(etaj+eps)*((etaj+eps)*x(2))^(-theta-1) )]*[(1-tau2h)*(etaj+eps)*x(2)+(1+r)*bl+T2]^(-sigma)-k*x(2)^(1/nu)]+...
                  + [ [(1-tau2l)*(etaj-eps)+ (etaj-eps)*x(2)*( -theta*delta*(etaj-eps)*((etaj-eps)*x(2))^(-theta-1) )]*[(1-tau2l)*(etaj-eps)*x(2)+(1+r)*bl+T2]^(-sigma)-k*x(2)^(1/nu)] )];

end