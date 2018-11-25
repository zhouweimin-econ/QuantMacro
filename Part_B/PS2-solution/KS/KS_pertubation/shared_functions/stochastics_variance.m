function [P_H,grid,par]=stochastics_variance(par, mpar, grid)
% stochastics_variance generates transition probabilities for:
% h: P_H,
% s: P_S,
% sxh: P
% and hgrid, sgrid on which the probability matrices are defined.
% Also an index function INDEX_H is created to read the individual
% productivity state h from the sxh state.

% Copyright (c) 2014-02-28
% Christian Bayer, Ralph Lutticke, Lien Pham-Dao, and Volker Tjaden

% First for human capital
% Generate transition probabilities and grid
% [hgrid,P_H,grid.boundsH] = Tauchen(par.rhoH,mpar.nh,1, 0, 'importance'); % LR variance = 1
% % Correct long run variance for *human capital*
% hgrid               = hgrid*par.sigmaH/sqrt(1-par.rhoH^2);
% grid.h               = exp(hgrid); % Levels instead of Logs
% 
% P_H     = transition(mpar.nh,par.rhoH,sqrt(1-par.rhoH^2),grid.boundsH);
% P_H=P_H./repmat(sum(P_H,2),[1 mpar.nh]);

[P_H, grid.h] = rouwen(par.rhoH, 0, par.sigmaH/sqrt(1-par.rhoH^2), mpar.nh);
grid.h=exp(grid.h');

Paux=P_H^1000^1000;
hh=Paux(1,1:mpar.nh)*grid.h(1:mpar.nh)';
par.H=hh(1); %Total Employment

end