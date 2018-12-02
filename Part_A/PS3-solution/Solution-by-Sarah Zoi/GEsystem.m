function [output]=GEsystem(sigma,beta,r,eta,y0,tau,T,eps)

n=size(y0)*size(eta);
astar=asavlab(sigma,beta,r,eta,y0,tau,T,eps);
mktclear=[sum(astar);
    T*n]

end