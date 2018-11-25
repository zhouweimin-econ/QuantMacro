function [P_H,grid,par]=ExTransitions(S,grid,mpar,par)

aux   = sqrt(S)*sqrt(1-par.rhoH^2); % Short run std

P_H     = transition(mpar.nh,par.rhoH,aux,grid.boundsH);
P_H=P_H./repmat(sum(P_H,2),[1 mpar.nh]);

end