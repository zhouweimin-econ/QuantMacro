function [c_a_star,m_a_star,k_a_star,c_n_star,m_n_star]=EGM_policyupdate(EVm,EVk,Qminus,PIminus,RBminus,inc,mesh,grid,par,mpar)

%% EGM Step 1:
EMU = par.beta*reshape(EVm,[mpar.nm,mpar.nk,mpar.nh]);
c_new = 1./(EMU.^(1/par.xi));

% Calculate assets consistent with choices being (m')
% Calculate initial money position from the budget constraint,
% that leads to the optimal consumption choice
m_star_n = (c_new + mesh.m - inc.labor - inc.rent);
m_star_n = (m_star_n<0).*m_star_n./((RBminus+par.borrwedge)/PIminus)...
    + (m_star_n>=0).*m_star_n./(RBminus/PIminus);

% Identify binding constraints
binding_constraints = mesh.m < repmat(m_star_n(1,:,:),[mpar.nm 1 1]);

% Consumption when drawing assets m' to zero: Eat all Resources
Resource = inc.labor + inc.rent + inc.money;

% Next step: Interpolate w_guess and c_guess from new k-grids
% using c(s,h,k'), k(s,h,k'K

m_star_n = reshape(m_star_n,[mpar.nm mpar.nk*mpar.nh]);
c_n_aux  = reshape(c_new,[mpar.nm mpar.nk*mpar.nh]);

% Interpolate grid.m and c_n_aux defined on m_star_n over grid.m
% Check monotonicity of m_star_n
if max(sum(abs(diff(sign(diff(m_star_n))))))~=0
    warning('non monotone future liquid asset choice encountered')
end

c_update=zeros(mpar.nm,mpar.nk*mpar.nh);
m_update=zeros(mpar.nm,mpar.nk*mpar.nh);
for hh=1:mpar.nk*mpar.nh
    Savings=griddedInterpolant(m_star_n(:,hh),grid.m); % generate savings function a(s,a*)=a'
    m_update(:,hh)=Savings(grid.m); % Obtain m'(m,h) by Interpolation
    Consumption=griddedInterpolant(m_star_n(:,hh),c_n_aux(:,hh)); % generate consumption function c(s,a*(s,a'))
    c_update(:,hh)=Consumption(grid.m);  % Obtain c(m,h) by interpolation (notice this is out of grid, used linear interpolation)
end
c_n_star = reshape(c_update,[mpar.nm, mpar.nk, mpar.nh]);
m_n_star = reshape(m_update,[mpar.nm, mpar.nk, mpar.nh]);

c_n_star(binding_constraints) = Resource(binding_constraints)-grid.m(1);
m_n_star(binding_constraints) = min(grid.m);

m_n_star(m_n_star>grid.m(end)) = grid.m(end);
%% EGM Step 2: Find Optimal Portfolio Combinations

term1=par.beta*reshape(EVk,[mpar.nm,mpar.nk,mpar.nh]);

E_return_diff = term1./Qminus-EMU;

% Check quasi-monotonicity of E_return_diff
if max(sum(abs(diff(sign(E_return_diff)))))>2
    warning('multiple roots of portfolio choice encountered')
end

% Find an m_a* for given k' that solves the difference equation (45)
m_a_aux = Fastroot(grid.m,E_return_diff);
m_a_aux = max(m_a_aux,grid.m(1)); %Use non-negativity constraint and monotonicity
m_a_aux = min(m_a_aux,grid.m(end)); % Do not allow for extrapolation
m_a_aux = reshape(m_a_aux, [mpar.nk,mpar.nh]);

%% EGM Step 3:
% Constraints for money and capital are not binding
EMU = reshape(EMU,[mpar.nm, mpar.nk*mpar.nh]);

% Interpolation of psi-function at m*_n(m,k)
[~,idx] = histc(m_a_aux,grid.m); % find indexes on grid next smallest to optimal policy
idx(m_a_aux<=grid.m(1))   = 1; %if below minimum
idx(m_a_aux>=grid.m(end)) = mpar.nm-1; %if above maximum
step = diff(grid.m); %Stepsize on grid
s = (m_a_aux - grid.m(idx))./step(idx);%Distance of optimal policy to next grid point

aux_index = (0:(mpar.nk*mpar.nh)-1)*mpar.nm; %aux for linear indexes
aux3      = EMU(idx(:)+aux_index(:)); %calculate linear indexes

% Interpolate EMU(m',k',s'*h',M',K') over m*_n(k'), m-dim is dropped
EMU_star        = aux3 + s(:).*(EMU(idx(:) + aux_index(:)+1)-aux3); %linear interpolation

c_a_aux         = 1./(EMU_star.^(1/par.xi));
cap_expenditure = shiftdim(inc.capital(1,:,:));
auxL            = shiftdim(inc.labor(1,:,:));
% Resources that lead to capital choice k'
% = c + m*(k') + k' - w*h*N = value of todays cap and money holdings
Resource = c_a_aux + m_a_aux(:) + cap_expenditure(:) - auxL(:);

c_a_aux  = reshape(c_a_aux, [mpar.nk mpar.nh]);
Resource = reshape(Resource, [mpar.nk mpar.nh]);

% Money constraint is not binding, but capital constraint is binding
m_star_zero = shiftdim(m_a_aux(1,:)); % Money holdings that correspond to k'=0:  m*(k=0)

% Use consumption at k'=0 from constrained problem, when m' is on grid
aux_c     = reshape(c_new(:,1,:),[mpar.nm, mpar.nh]);
aux_inc   = reshape(inc.labor(1,1,:),[1 mpar.nh]);
cons_list = cell(mpar.nh,1);
res_list  = cell(mpar.nh,1);
mon_list  = cell(mpar.nh,1);
cap_list  = cell(mpar.nh,1);


for j=1:mpar.nh
    % When choosing zero capital holdings, HHs might still want to choose money
    % holdings smaller than m*(k'=0)
    if m_star_zero(j)>grid.m(1)
        % Calculate consumption policies, when HHs chooses money holdings
        % lower than m*(k'=0) and capital holdings k'=0 and save them in cons_list
        log_index    = grid.m<m_star_zero(j);
        % aux_c is the consumption policy under no cap. adj. (fix k�=0), for m�<m_a*(k'=0)
        c_k_cons     = aux_c(log_index,j);
        cons_list{j} = c_k_cons; %Consumption at k'=0, m'<m_a*(0)
        % Required Resources: Money choice + Consumption - labor income
        % Resources that lead to k'=0 and m'<m*(k'=0)
        res_list{j} = grid.m(log_index)' + c_k_cons - aux_inc(j);
        mon_list{j} = grid.m(log_index)';
        cap_list{j} = zeros(sum(log_index),1);
    end
end

% Merge lists
c_a_aux  = reshape(c_a_aux,[mpar.nk mpar.nh]);
m_a_aux  = reshape(m_a_aux,[mpar.nk mpar.nh]);
Resource = reshape(Resource,[mpar.nk mpar.nh]);
for j=1:mpar.nh
    cons_list{j} = [cons_list{j}; c_a_aux(:,j)];
    res_list{j}  = [res_list{j};  Resource(:,j)];
    mon_list{j}  = [mon_list{j};  m_a_aux(:,j)];
    cap_list{j}  = [cap_list{j};  grid.k'];
end

%% EGM Step 4: Interpolate back to fixed grid
c_a_star = NaN([mpar.nm*mpar.nk, mpar.nh]);
m_a_star = NaN([mpar.nm*mpar.nk, mpar.nh]);
k_a_star = NaN([mpar.nm*mpar.nk, mpar.nh]);
Resource_grid  = reshape(inc.capital+inc.money+inc.rent,[mpar.nm*mpar.nk, mpar.nh]);
labor_inc_grid = reshape(inc.labor,[mpar.nm*mpar.nk, mpar.nh]);

for j=1:mpar.nh
    log_index=Resource_grid(:,j)<res_list{j}(1);
    
    % when at most one constraint binds:
    % Check monotonicity of resources
    if max(sum(abs(diff(sign(diff(res_list{j}))))))~=0
        warning('non monotone resource list encountered')
    end

    cons=griddedInterpolant(res_list{j},cons_list{j});
    c_a_star(:,j)=cons(Resource_grid(:,j));
    mon=griddedInterpolant(res_list{j},mon_list{j});
    m_a_star(:,j)=mon(Resource_grid(:,j));
    cap=griddedInterpolant(res_list{j},cap_list{j});
    k_a_star(:,j)=cap(Resource_grid(:,j));
    % Lowest value of res_list corresponds to m_a'=0 and k_a'=0.
    
    % Any resources on grid smaller then res_list imply that HHs consume all
    % resources plus income.
    % When both constraints are binding:
    c_a_star(log_index,j) = Resource_grid(log_index,j) + labor_inc_grid(log_index,j)-grid.m(1);
    m_a_star(log_index,j) = grid.m(1);
    k_a_star(log_index,j) = 0;
end

c_a_star = reshape(c_a_star,[mpar.nm ,mpar.nk, mpar.nh]);
k_a_star = reshape(k_a_star,[mpar.nm ,mpar.nk, mpar.nh]);
m_a_star = reshape(m_a_star,[mpar.nm ,mpar.nk, mpar.nh]);

k_a_star(k_a_star>grid.k(end)) = grid.k(end);
m_a_star(m_a_star>grid.m(end)) = grid.m(end);
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

aux_index  = (0:numel(fx)/numel(xgrid)-1)*numel(xgrid); %aux for linear indexes
aux_index2 = (0:numel(fx)/numel(xgrid)-1)*(numel(xgrid)-1);
fxx  = fx(idx+aux_index);
xl   = xgrid(idx)';
dx   = dxgrid(idx)';
dfxx = dfx(idx+aux_index2);
% Because function is piecewise linear in gridpoints, one newton step is
% enough to find the solution
roots = xl-fxx.*dx./dfxx;

roots(idx_min)=xgrid(1); %constrained choice
roots(idx_max)=xgrid(end); % no-extrapolation
end
