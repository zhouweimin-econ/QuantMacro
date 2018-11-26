function [P_H,grid,par]=stochastics_variance(par, mpar, grid)
% stochastics_variance generates transition probabilities for:
% h: P_H

% Copyright (c) 2017-12-19
% Christian Bayer, Ralph Lutticke

% First for human capital
% Generate transition probabilities and grid
[hgrid,P_H,boundsH] = Tauchen(par.rhoH,mpar.nh-1,1, 0, mpar.tauchen); % LR variance = 1
% Correct long run variance for *human capital*
hgrid               = hgrid*par.sigmaH/sqrt(1-par.rhoH^2);
hgrid               = exp(hgrid); % Levels instead of Logs

grid.h=[hgrid 1];

P_H     = transition(mpar.nh-1,par.rhoH,sqrt(1-par.rhoH^2),boundsH);

% Transitions to entrepreneur state
P_H=[P_H repmat(mpar.in,[mpar.nh-1 1])];
lastrow=[repmat(0,[1, mpar.nh-1]) 1-mpar.out];
lastrow(ceil(mpar.nh/2))=mpar.out;
P_H=[P_H; lastrow];
P_H=P_H./repmat(sum(P_H,2),[1 mpar.nh]);

Paux=P_H^1000^1000;
hh=Paux(1,1:mpar.nh-1)*grid.h(1:mpar.nh-1)';
par.H=hh(1); % Total Employment
par.profitshare=Paux(end,end)^-1; % Profit per household
grid.boundsH=boundsH;
end