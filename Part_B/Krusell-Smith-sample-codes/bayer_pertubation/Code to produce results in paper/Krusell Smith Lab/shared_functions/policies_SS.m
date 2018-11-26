function [c_new,m_star,distPOL] = policies_SS(c_guess, grid, inc,P,mpar,par)
% POLICIES solves for the household policies for consumption and
% capital holdings by EGM.

% =========================================================================
% Part of the Matlab code to accompany the paper
% 'Solving heterogeneous agent models in discrete time with
% many idiosyncratic states by perturbation methods', CEPR Discussion Paper
% DP13071
% by Christian Bayer and Ralph Luetticke
% http://www.erc-hump-uni-bonn.net/
% =========================================================================


%% Apply EGM to solve for optimal policies and marginal utilities
money_expense  = repmat(grid.m',[1 mpar.nh]);
distC = 99999;

count=0;
while max([distC])>mpar.crit
    count=count+1;

    %% Step 1: Update policies for only money adjustment
    mutil_c = 1./(c_guess.^par.xi); % marginal utility at consumption policy no adjustment

    mutil_c=par.RB.*mutil_c; %take return on money into account
    aux=reshape(permute(mutil_c,[2 1]),[mpar.nh mpar.nm]);
    % form expectations
    EMU_aux = par.beta*permute(reshape(P*aux,[mpar.nh mpar.nm]),[2 1]);

    c_aux = 1./(EMU_aux.^(1/par.xi));

    % Take borrowing constraint into account
    [c_new,m_star]=EGM_Step1_b(grid,inc,money_expense,c_aux,mpar,par);

    m_star(m_star>grid.m(end)) = grid.m(end); % No extrapolation
   
    %% Step 6: Check convergence of policies
    distC = max((abs(c_guess(:)-c_new(:))));
    
    % Update c policy guesses
    c_guess=c_new;

end
distPOL=[distC];
end

%% SUBFUNCTIONS

function [c_update,m_update]=EGM_Step1_b(grid,inc,money_expense,c_aux,mpar,par)
%%EGM_Step1_b computes the optimal consumption and corresponding optimal money
% holdings in case the capital stock cannot be adjusted by taking the budget
% constraint into account.
% c_update(m,k,s*h,M,K):    Update for consumption policy under no-adj.
% m_update(m,k,s*h,M,K):    Update for money policy under no-adj.

%% EGM: Calculate assets consistent with choices being (m')
% Calculate initial money position from the budget constraint,
% that leads to the optimal consumption choice
m_star = (c_aux + money_expense - inc.labor - inc.profits);
m_star = m_star./par.RB;

% Identify binding constraints
binding_constraints = money_expense < repmat(m_star(1,:),[mpar.nm 1 ]);

% Consumption when drawing assets m' to zero: Eat all Resources
Resource = inc.labor  + inc.money + inc.profits;

%% Next step: Interpolate w_guess and c_guess from new k-grids
% using c(s,h,k',K), k(s,h,k'K

m_star = reshape(m_star,[mpar.nm mpar.nh]);
c_aux= reshape(c_aux,[mpar.nm mpar.nh]);

%Interpolate grid.m and c_n_aux defined on m_star_n over grid.m
c_update=zeros(mpar.nm,mpar.nh);
m_update=zeros(mpar.nm,mpar.nh);

for hh=1:mpar.nh
    Savings=griddedInterpolant(m_star(:,hh),grid.m); % generate savings function a(s,a*)=a'
    m_update(:,hh)=Savings(grid.m); % Obtain m'(m,h) by Interpolation
    Consumption=griddedInterpolant(m_star(:,hh),c_aux(:,hh)); % generate consumption function c(s,a*(s,a'))
    c_update(:,hh)=Consumption(grid.m);  % Obtain c(m,h) by interpolation (notice this is out of grid, used linear interpolation)
end


c_update = reshape(c_update,[mpar.nm, mpar.nh]);
m_update = reshape(m_update,[mpar.nm, mpar.nh]);

c_update(binding_constraints) = Resource(binding_constraints)-grid.m(1);
m_update(binding_constraints) = min(grid.m);

end
