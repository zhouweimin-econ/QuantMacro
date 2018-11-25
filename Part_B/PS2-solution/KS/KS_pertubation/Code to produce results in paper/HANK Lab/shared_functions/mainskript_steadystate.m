
%% Set grids
grid.K=40; % Starting guess K

[grid]  = makegrid(mpar,grid);

%% Use Tauchen method to approximate state space

[P_H,grid,par]=stochastics_variance(par, mpar,grid);

[meshes.m,meshes.k,meshes.h] = ndgrid(grid.m,grid.k,grid.h);

%% Solve for steady state capital by bi-section

[c_n_guess,m_n_star,c_a_guess,m_a_star,cap_a_star,psi_guess,joint_distr,R_fc,W_fc,Profits_fc,Output,grid]...
    =steadystate(P_H,grid,meshes,mpar,par);

%% Calculate steady state capital and further statistics

SS_stats

%% Calculate Value Functions
% Calculate Marginal Values of Capital (k) and Liquid Assets(m)
RBRB = par.RB/par.PI + (meshes.m<0).*(par.borrwedge/par.PI);

% Liquid Asset
mutil_c_n = 1./(c_n_guess.^par.xi); % marginal utility at consumption policy no adjustment
mutil_c_a = 1./(c_a_guess.^par.xi); % marginal utility at consumption policy adjustment
mutil_c = par.nu.*mutil_c_a + (1-par.nu).*mutil_c_n; % Expected marginal utility at consumption policy (w &w/o adjustment)
Vm = RBRB.*mutil_c; %take return on money into account
Vm = reshape(reshape(Vm,[mpar.nm*mpar.nk mpar.nh]),[mpar.nm,mpar.nk mpar.nh]);

% Capital
Vk = par.nu.*(par.R+par.Q).*mutil_c_a + (1-par.nu).*par.R.*mutil_c_n + (1-par.nu).*psi_guess; % Expected marginal utility at consumption policy (w &w/o adjustment)
Vk = reshape(reshape(Vk,[mpar.nm*mpar.nk mpar.nh]),[mpar.nm,mpar.nk mpar.nh]);

%% Produce non-parametric Copula
cum_dist = cumsum(cumsum(cumsum(joint_distr,1),2),3);
marginal_m = cumsum(squeeze(sum(sum(joint_distr,2),3)));
marginal_k = cumsum(squeeze(sum(sum(joint_distr,1),3)));
marginal_h = cumsum(squeeze(sum(sum(joint_distr,2),1)));

Copula = griddedInterpolant({marginal_m,marginal_k',marginal_h},cum_dist,'spline');

%% Save
filename=[casenamemaster];
save(filename)