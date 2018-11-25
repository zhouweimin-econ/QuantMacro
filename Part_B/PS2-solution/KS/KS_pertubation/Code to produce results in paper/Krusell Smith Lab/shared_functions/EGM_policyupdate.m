function [c_star,m_star]=EGM_policyupdate(EVm,Rminus,Qminus,inc,mesh,grid,par,mpar)

%% EGM Step 1:
EMU = par.beta*reshape(EVm,[mpar.nm,mpar.nh]);
c_new = 1./(EMU.^(1/par.xi));

% Calculate assets consistent with choices being (m')
% Calculate initial asset position from the budget constraint,
% that leads to the optimal consumption choice
m_star_n = (c_new + (Qminus)*mesh.m - inc.labor - inc.profits);
m_star_n = m_star_n./(Rminus+Qminus);

% Identify binding constraints
binding_constraints = mesh.m < repmat(m_star_n(1,:),[mpar.nm 1]);

% Consumption when drawing assets m' to zero: Eat all Resources
Resource = inc.labor + inc.money + inc.profits;

m_star_n = reshape(m_star_n,[mpar.nm mpar.nh]);
c_n_aux  = reshape(c_new,[mpar.nm mpar.nh]);

% Interpolate grid.m and c_n_aux defined on m_star_n over grid.m
% Check monotonicity of m_star_n
if max(sum(abs(diff(sign(diff(m_star_n))))))~=0
    warning('non monotone future liquid asset choice encountered')
end

c_star=zeros(mpar.nm,mpar.nh);
m_star=zeros(mpar.nm,mpar.nh);

for hh=1:mpar.nh
    Savings=griddedInterpolant(m_star_n(:,hh),grid.m); % generate savings function a(s,a*)=a'
    m_star(:,hh)=Savings(grid.m); % Obtain m'(m,h) by Interpolation
    Consumption=griddedInterpolant(m_star_n(:,hh),c_n_aux(:,hh)); % generate consumption function c(s,a*(s,a'))
    c_star(:,hh)=Consumption(grid.m);  % Obtain c(m,h) by interpolation (notice this is out of grid, used linear interpolation)
end

c_star(binding_constraints) = Resource(binding_constraints)-Qminus*grid.m(1);
m_star(binding_constraints) = min(grid.m);

m_star(m_star>grid.m(end)) = grid.m(end);

end
