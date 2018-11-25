function [P_H,grid,par]=ExTransitions(S,grid,mpar,par)

% Generate transition probabilities for given boundsH and short run stdev

aux   = sqrt(S)*sqrt(1-par.rhoH^2); % Short run std

P     = transition(mpar.nh-1,par.rhoH,aux,grid.boundsH);

P_H=[P repmat(mpar.in,[mpar.nh-1 1])];
lastrow=[repmat(0,[1, mpar.nh-1]) 1-mpar.out];
lastrow(ceil(mpar.nh/2))=mpar.out;
P_H=[P_H; lastrow];
P_H=P_H./repmat(sum(P_H,2),[1 mpar.nh]);


end