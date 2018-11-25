function [c_n_new,m_n_star,c_a_new,m_a_star,cap_a_star,psi_new,distPOL] = policies_SS(c_a_guess,c_n_guess,psi_guess, grid, inc, RR,RBRB,P,mpar,par,meshes)
% POLICIES solves for the household policies for consumption, money and
% capital holdings in the adjustment and non-adjustment case. Given returns
% {q,R_fc,PI} and interpolation matrix {P_H} a variant of the
% endogeneous gridpoint method is used to iterate on the first order
% conditions until convergence.
% Policies in adjustment case:
% {c_a_new,m_a_star,cap_a_star}
% Policies in non-adjustment case:
% {c_n_new,m_n_star}

% =========================================================================
% Part of the Matlab code to accompany the paper
% 'Solving heterogeneous agent models in discrete time with
% many idiosyncratic states by perturbation methods', CEPR Discussion Paper
% DP13071
% by Christian Bayer and Ralph Luetticke
% http://www.erc-hump-uni-bonn.net/
% =========================================================================
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The most recent version or successor of the above-mentioned paper is
% properly cited in all work and publications that benefit from insights
% derived from the Software.
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% =========================================================================

%% Apply EGM to solve for optimal policies and marginal utilities
distC_n = 99999;
distPSI = distC_n;
distC_a = distPSI;
mutil_c_n = 1./(c_n_guess.^par.xi); % marginal utility at consumption policy no adjustment
mutil_c_a = 1./(c_a_guess.^par.xi); % marginal utility at consumption policy adjustment
mutil_c = par.nu.*mutil_c_a + (1-par.nu).*mutil_c_n; % Expected marginal utility at consumption policy (w &w/o adjustment)

count=0;
while max([distC_n distC_a distPSI])>mpar.crit && count<10000
    count=count+1;
    
    %% Step 1: Update policies for only money adjustment
    mutil_c=RBRB.*mutil_c; %take return on money into account
    aux=reshape(permute(mutil_c,[3 1 2]),[mpar.nh mpar.nm*mpar.nk]);
    % form expectations
    EMU_aux = par.beta*permute(reshape(P*aux,[mpar.nh mpar.nm mpar.nk]),[2 3 1]);
    
    c_n_aux = 1./(EMU_aux.^(1/par.xi));
    
    % Take borrowing constraint into account
    [c_n_new,m_n_star]=EGM_Step1_b(grid,inc,c_n_aux,mpar,par,meshes);    
    m_n_star(m_n_star>grid.m(end)) = grid.m(end); % not extrapolation
    
    %% Step 2: Find for every k on grid some off-grid m*(k')
    m_a_star_aux=EGM_Step2_SS(mutil_c_n,mutil_c_a, psi_guess, grid,P,RBRB,RR,par,mpar);
    
    %% Step 3: Solve for initial resources / consumption in  adjustment case
    [ ~,~,cons_list, res_list, mon_list,cap_list ] = EGM_Step3(EMU_aux,grid,inc,m_a_star_aux,c_n_aux,mpar,par);
    
    %% Step 4: Interpolate Consumption Policy
    [ c_a_new, m_a_star,cap_a_star ] = EGM_Step4( cons_list,res_list, mon_list,cap_list,inc,mpar,grid );
    
    cap_a_star(cap_a_star>grid.k(end)) = grid.k(end); % not extrapolation
    m_a_star(m_a_star>grid.m(end)) = grid.m(end); % not extrapolation
    
    %% Step 5: Update ~psi
    mutil_c_n = 1./(c_n_new.^par.xi); % marginal utility at consumption policy no adjustment
    mutil_c_a = 1./(c_a_new.^par.xi); % marginal utility at consumption policy adjustment
    mutil_c = par.nu.*mutil_c_a + (1-par.nu).*mutil_c_n; % Expected marginal utility at consumption policy (w &w/o adjustment)
    
    %% VFI analogue in updating psi
    term1=((par.nu.* mutil_c_a .*(par.Q + RR))...
        +((1-par.nu).* mutil_c_n .* RR)...
        +(1-par.nu).* psi_guess);
    
    aux = reshape(permute(term1,[3 1 2]),[mpar.nh mpar.nm*mpar.nk]);
    E_rhs_psi = par.beta*permute(reshape(P*aux,[mpar.nh mpar.nm mpar.nk]),[2 3 1]);
    
    E_rhs_psi=reshape(E_rhs_psi, [mpar.nm mpar.nk*mpar.nh]);
    m_n_star=reshape(m_n_star, [mpar.nm mpar.nk*mpar.nh]);
    
    % Interpolation of psi-function at m*_n(m,k)
    [~,idx]=histc(m_n_star,grid.m); % find indexes on grid next smallest to optimal policy
    idx(m_n_star<=grid.m(1))=1; %if below minimum
    idx(m_n_star>=grid.m(end))=mpar.nm-1; %if above maximum
    step=diff(grid.m); %Stepsize on grid
    s = (m_n_star - grid.m(idx))./step(idx);%Distance of optimal policy to next grid point
    
    aux_index=ones(mpar.nm,1)*(0:(mpar.nk*mpar.nh)-1)*mpar.nm; %aux for linear indexes
    aux3=E_rhs_psi(idx(:)+aux_index(:)); %calculate linear indexes
    
    psi_new = aux3 + s(:).*(E_rhs_psi(idx(:)+aux_index(:)+1)-aux3); %linear interpolation
    psi_new = reshape(psi_new,[mpar.nm, mpar.nk, mpar.nh]);
    m_n_star=reshape(m_n_star, [mpar.nm mpar.nk mpar.nh]);
    distPSI = max((abs(psi_guess(:)-psi_new(:))));
    
    %% Step 6: Check convergence of policies
    distC_n = max((abs(c_n_guess(:)-c_n_new(:))));
    distC_a = max((abs(c_a_guess(:)-c_a_new(:))));
    
    % Update c policy guesses
    c_n_guess=c_n_new;
    c_a_guess=c_a_new;
    psi_guess=psi_new;
    
end
distPOL=[distC_n distC_a distPSI];
end

%% SUBFUNCTIONS

function [c_update,m_update]=EGM_Step1_b(grid,inc,c_n_aux,mpar,par,meshes)
%%EGM_Step1_b computes the optimal consumption and corresponding optimal money
% holdings in case the capital stock cannot be adjusted by taking the budget
% constraint into account.
% c_update(m,k,h):    Update for consumption policy under no-adj.
% m_update(m,k,h):    Update for money policy under no-adj.

%% EGM: Calculate assets consistent with choices being (m')
% Calculate initial money position from the budget constraint,
% that leads to the optimal consumption choice
m_star_n = (c_n_aux + meshes.m - inc.labor - inc.rent - inc.profits);
m_star_n = (m_star_n<0).*m_star_n./((par.RB+par.borrwedge)./par.PI)...
    + (m_star_n>=0).*m_star_n./(par.RB./par.PI);

% Identify binding constraints
binding_constraints = meshes.m < repmat(m_star_n(1,:,:),[mpar.nm 1 1]);

% Consumption when drawing assets m' to zero: Eat all Resources
Resource = inc.labor + inc.rent + inc.money + inc.profits;

%% Next step: Interpolate w_guess and c_guess from new k-grids
% using c(s,h,k',K), k(s,h,k'K

m_star_n = reshape(m_star_n,[mpar.nm mpar.nk*mpar.nh]);
c_n_aux= reshape(c_n_aux,[mpar.nm mpar.nk*mpar.nh]);

%Interpolate grid.m and c_n_aux defined on m_star_n over grid.m
c_update=zeros(mpar.nm,mpar.nk*mpar.nh);
m_update=zeros(mpar.nm,mpar.nk*mpar.nh);
for hh=1:mpar.nk*mpar.nh
    Savings=griddedInterpolant(m_star_n(:,hh),grid.m); % generate savings function a(s,a*)=a'
    m_update(:,hh)=Savings(grid.m); % Obtain m'(m,h) by Interpolation
    Consumption=griddedInterpolant(m_star_n(:,hh),c_n_aux(:,hh)); % generate consumption function c(s,a*(s,a'))
    c_update(:,hh)=Consumption(grid.m);  % Obtain c(m,h) by interpolation (notice this is out of grid, used linear interpolation)
end
c_update = reshape(c_update,[mpar.nm, mpar.nk, mpar.nh]);
m_update = reshape(m_update,[mpar.nm, mpar.nk, mpar.nh]);

c_update(binding_constraints) = Resource(binding_constraints)-grid.m(1);
m_update(binding_constraints) = min(grid.m);

end

function  mstar=EGM_Step2_SS(mutil_c_n,mutil_c_a, psi_guess, grid,P,RBRB,RR,par,mpar)
%EGM_Step2 finds the optimal off-grid money choice, m*(k'), given k' (and h)

term1=((par.nu.* mutil_c_a .*(par.Q + RR))+((1-par.nu).* mutil_c_n .* RR)+(1-par.nu).* psi_guess);
aux=reshape(permute(term1,[3 1 2]),[mpar.nh mpar.nm*mpar.nk]);
term1=par.beta*permute(reshape(P*aux,[mpar.nh mpar.nm mpar.nk]),[2 3 1]);

term2=RBRB.*(par.nu.* mutil_c_a +(1-par.nu).* mutil_c_n);
aux=reshape(permute(term2,[3 1 2]),[mpar.nh mpar.nm*mpar.nk]);
term2=par.beta*permute(reshape(P*aux,[mpar.nh mpar.nm mpar.nk]),[2 3 1]);


% Equation (59) in Appedix B.4.
E_return_diff=term1./par.Q-term2;

% Find an m*_n for given k' that solves the difference equation (59)
mstar=Fastroot(grid.m,E_return_diff);
mstar=max(mstar,grid.m(1)); %Use non-negativity constraint and monotonicity
mstar=min(mstar,grid.m(end)); % Do not allow for extrapolation
mstar=reshape(mstar, [mpar.nk,mpar.nh]);
end

function roots = Fastroot(xgrid,fx)
%fast linear interpolation root finding
%(=one Newton step at largest negative function value)
%   stripped down version of interp1 that accepts multiple inputs (max 3)
%   that are interpolated over the same grids x & xi
xgrid=xgrid(:);
fx=reshape(fx,[numel(xgrid),numel(fx)/numel(xgrid)]);

dxgrid=diff(xgrid);
dfx=diff(fx);
idx=ones(1,numel(fx)/numel(xgrid));

% Make use of the fact that the difference equation is monotonically
% increasing in m
idx_min=(fx(1,:)>0); %Corner solutions left (if no solution x* to f(x)=0 exists)
idx_max=(fx(end,:)<0); %Corner solutions right (if no solution x* to f(x)=0 exists)
index=find(and(not(idx_min),not(idx_max))); %interior solutions (if solution x* to f(x)=0 exists)

% Find index of two gridpoints where sign of fx changes from positive to negative,
[~,idx(index)]=max(diff(sign(fx(:,index))));

aux_index=(0:numel(fx)/numel(xgrid)-1)*numel(xgrid); %aux for linear indexes
aux_index2=(0:numel(fx)/numel(xgrid)-1)*(numel(xgrid)-1);
fxx=fx(idx+aux_index);
xl=xgrid(idx)';
dx=dxgrid(idx)';
dfxx=dfx(idx+aux_index2);
% Because function is piecewise linear in gridpoints, one newton step is
% enough to find the solution
roots=xl-fxx.*dx./dfxx;

roots(idx_min)=xgrid(1); %constrained choice
roots(idx_max)=xgrid(end); % no-extrapolation
end

function [ c_star,Resource,cons_list, res_list,mon_list,cap_list ] = EGM_Step3(EMU,grid,inc,m_a_star,c_n_aux,mpar,par)
%EGM_Step3 returns the resources (res_list), consumption (cons_list)
% and money policy (mon_list) for given capital choice (cap_list).
% For k'=0, there doesn't need to be a unique corresponding m*. We get a
% list of consumption choices and resources for money choices m'<m* (mon_list) and cap
% choices k'=0 (cap_list) and merge them with consumption choices and
% resources that we obtain if capital constraint doesn't bind next period.
%
% c_star: optimal consumption policy as function of k',h (both
% constraints do not bind)
%
% Resource: required resource for c_star
%
% cons_list: optimal consumption policy if a) only k>=0 binds and b) both
% constraints do not bind
%
% res_list: Required resorces for cons_list
% c_n_aux: consumption in t as function of t+1 grid (constrained version)

%% Constraints for money and capital are not binding
EMU=reshape(EMU,[mpar.nm, mpar.nk*mpar.nh]);

% Interpolation of psi-function at m*_n(m,k)
[~,idx]=histc(m_a_star,grid.m); % find indexes on grid next smallest to optimal policy
idx(m_a_star<=grid.m(1))=1; %if below minimum
idx(m_a_star>=grid.m(end))=mpar.nm-1; %if above maximum
step=diff(grid.m); %Stepsize on grid
s = (m_a_star - grid.m(idx))./step(idx);%Distance of optimal policy to next grid point

aux_index=(0:(mpar.nk*mpar.nh)-1)*mpar.nm; %aux for linear indexes
aux3=EMU(idx(:)+aux_index(:)); %calculate linear indexes

% Interpolate EMU(m',k',h') over m*_n(k'), m-dim is dropped
EMU_star = aux3 + s(:).*(EMU(idx(:)+aux_index(:)+1)-aux3); %linear interpolation

c_star = 1./(EMU_star.^(1/par.xi));
cap_expenditure = shiftdim(inc.capital(1,:,:));
auxL=shiftdim(inc.labor(1,:,:));
auxP=shiftdim(inc.profits(1,:,:));
% Resources that lead to capital choice k'
% = c + m*(k') + k' - w*h*N = value of todays cap and money holdings
Resource = c_star + m_a_star(:) + cap_expenditure(:) - auxL(:) - auxP(:);

c_star = reshape(c_star, [mpar.nk mpar.nh]);
Resource = reshape(Resource, [mpar.nk mpar.nh]);

%% Money constraint is not binding, but capital constraint is binding
m_star_zero=shiftdim(m_a_star(1,:)); % Money holdings that correspond to k'=0:  m*(k=0)

% Use consumption at k'=0 from constrained problem, when m' is on grid
aux_c=reshape(c_n_aux(:,1,:),[mpar.nm, mpar.nh]);
aux_inc=reshape(inc.labor(1,1,:)+inc.profits(1,1,:),[1 mpar.nh]);
cons_list=cell(mpar.nh,1);
res_list=cell(mpar.nh,1);
mon_list=cell(mpar.nh,1);
cap_list=cell(mpar.nh,1);


for j=1:mpar.nh
    % When choosing zero capital holdings, HHs might still want to choose money
    % holdings smaller than m*(k'=0)
    if m_star_zero(j)>grid.m(1)
        % Calculate consumption policies, when HHs chooses money holdings
        % lower than m*(k'=0) and capital holdings k'=0 and save them in cons_list
        log_index=grid.m<m_star_zero(j);
        % aux_c is the consumption policy under no cap. adj. (fix k�=0), for m�<m_a*(k'=0)
        c_k_cons=aux_c(log_index,j);
        cons_list{j}=c_k_cons; %Consumption at k'=0, m'<m_a*(0)
        % Required Resources: Money choice + Consumption - labor income
        % Resources that lead to k'=0 and m'<m*(k'=0)
        res_list{j}=grid.m(log_index)' + c_k_cons - aux_inc(j);
        mon_list{j}=grid.m(log_index)';
        cap_list{j}=zeros(sum(log_index),1);
    end
end

%% Merge lists
c_star=reshape(c_star,[mpar.nk mpar.nh]);
m_a_star=reshape(m_a_star,[mpar.nk mpar.nh]);
Resource=reshape(Resource,[mpar.nk mpar.nh]);
for j=1:mpar.nh
    cons_list{j}=[cons_list{j} ;c_star(:,j)];
    res_list{j}=[ res_list{j};Resource(:,j)];
    mon_list{j}=[ mon_list{j};m_a_star(:,j)];
    cap_list{j}=[cap_list{j}; grid.k'];
end


end
function [ c_a_new,m_a_new,k_a_new ] = EGM_Step4( cons_list,res_list,mon_list,cap_list,inc,mpar,grid )
%EGM_Step4 obtains consumption, money, and capital policy under adjustment.
%   The function uses the {(cons_list{j},res_list{j})} as measurement
%   points. The consumption function in (m,k) can be obtained from
%   interpolation by using the total resources available at (m,k): R(m,k)=qk+m/pi.
%
%   c_a_new(m,k,h): Update for consumption policy under adjustment
%   m_a_new(m,k,h): Update for money policy under adjustment
%   k_a_new(m,k,h): Update for capital policy under adjustment


c_a_new=NaN([mpar.nm*mpar.nk, mpar.nh]);
m_a_new=NaN([mpar.nm*mpar.nk, mpar.nh]);
k_a_new=NaN([mpar.nm*mpar.nk, mpar.nh]);
Resource_grid=reshape(inc.capital+inc.money+inc.rent,[mpar.nm*mpar.nk, mpar.nh]);
labor_inc_grid=reshape(inc.labor + inc.profits,[mpar.nm*mpar.nk, mpar.nh]);

for j=1:mpar.nh
    log_index=Resource_grid(:,j)<res_list{j}(1);
    
    % when at most one constraint binds:
    cons=griddedInterpolant(res_list{j},cons_list{j});
    c_a_new(:,j)=cons(Resource_grid(:,j));
    mon=griddedInterpolant(res_list{j},mon_list{j});
    m_a_new(:,j)=mon(Resource_grid(:,j));
    cap=griddedInterpolant(res_list{j},cap_list{j});
    k_a_new(:,j)=cap(Resource_grid(:,j));
    
    % Lowest value of res_list corresponds to m_a'=0 and k_a'=0.
    % Any resources on grid smaller then res_list imply that HHs consume all
    % resources plus income.
    % When both constraints are binding:
    c_a_new(log_index,j) = Resource_grid(log_index,j) + labor_inc_grid(log_index,j)-grid.m(1);
    m_a_new(log_index,j) = grid.m(1);
    k_a_new(log_index,j) = 0;
end

c_a_new=reshape(c_a_new,[mpar.nm ,mpar.nk, mpar.nh]);
k_a_new=reshape(k_a_new,[mpar.nm ,mpar.nk, mpar.nh]);
m_a_new=reshape(m_a_new,[mpar.nm ,mpar.nk, mpar.nh]);

end
