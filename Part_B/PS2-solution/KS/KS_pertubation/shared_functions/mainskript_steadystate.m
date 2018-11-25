
%% Set grids
grid.K  = 90;

[grid]  = makegrid(mpar,grid);

%% Use Tauchen method to approximate state space
switch casenamemaster
    case 'SS_BASELINE_HANC_UNC'
        [P_H,grid,par]=stochastics_variance_UNC(par, mpar,grid);

    otherwise
        [P_H,grid,par]=stochastics_variance(par, mpar,grid);
end
[meshes.m,meshes.h] = ndgrid(grid.m,grid.h);

%% Solve for steady state

[c_guess,m_star,joint_distr,R_fc,W_fc,Output,N,grid]=steadystate(P_H,grid,mpar,par,meshes);

%% Calculate steady state capital and further statistics

SS_stats

%% Calculate Value Functions

% Liquid Asset
mutil_c = 1./(c_guess.^par.xi); % marginal utility at consumption policy no adjustment
Vm = (1+R_fc).*mutil_c; %take return on money into account
Vm = reshape(Vm,[mpar.nm mpar.nh]);

%% Produce non-parametric Copula
cum_dist = cumsum(cumsum(joint_distr,1),2);
marginal_m = cumsum(squeeze((sum(joint_distr,2))));
marginal_h = cumsum(squeeze((sum(joint_distr,1))));
%
Copula = griddedInterpolant({marginal_m,marginal_h},cum_dist,'spline');

%% Save
filename=[casenamemaster];
save(filename)
