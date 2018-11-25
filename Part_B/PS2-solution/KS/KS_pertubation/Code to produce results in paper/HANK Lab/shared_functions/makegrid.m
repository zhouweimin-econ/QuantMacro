function [grid]=makegrid(mpar,grid)

%% Construct basic grids
% This section defines the basic grids for the problem.
k_min = 0;
k_max = 20*grid.K;
grid.k = (exp(linspace(0,log(k_max - k_min+1),mpar.nk))-1+k_min);

m_min = -1.85;
m_max = 10*grid.K;
grid.m = (exp(exp(linspace(0,log(log(m_max - m_min+1)+1),mpar.nm-1))-1)-1+m_min);

grid.m=sort([grid.m 0]);



