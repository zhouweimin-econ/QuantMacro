% =======================================================================
%    This program solves the simple neoclassical growth model by value
%                  function iterations. 
% =======================================================================
clear;
clf;
tic
disp('Solving a deterministic growth model by value function iterations')
disp('.')

% set parameter values
alfa   = 0.321;             % elasticity of output to capital 
beta   = 0.988;             % subjective discount factor 
delta  = .031;               % depreciation rate 
Aval  = 1.;               % high value for technology 

% Compute deterministic steady-state of capital stock, consumption and
% investment
kstar  = ((1/beta-(1-delta))/(alfa*Aval))^(1/(alfa-1));
cstar = Aval*(kstar^alfa) - delta*kstar;
istar = delta*kstar;

%   form capital grid 
mink =   0.5*kstar;            %  minimum value of the capital grid 
maxk =   1.5* kstar;            % maximum value of the capital grid 
nk   = 2000;              % choose the number of grid-points
                          
% the grid is now computed 
kgrid = (linspace(mink,maxk,nk))';


%
%   form single period return function
R1=-10000*ones(nk,nk);    % initial values 

%  A. loop through all possible states (i.e. K)
i=1;

while i<=nk;
    k=kgrid(i);
    
    %  B. loop through all possible controls (i.e. K') @
    j=1;    
    while j<=nk; 
        kprime = kgrid(j);
        invest = kprime - (1-delta)*k; 
        consu = Aval*k^(alfa) - invest;
        if consu > 0; 
           R1(j,i)= log(consu);
        end
      j=j+1;
    end
    i=i+1;
end

%%
%  initialize some variables
v       = zeros(nk,1);                 % the value function: row = kgrid 
decis   = zeros(nk,1);                 % the policy function: row = kgrid 
iter    = 0;                           % counter 
metric1 = 10;                          % metric 1 for convergence 
metric2 = 10;                          % metric 2 for convergence

metric = max(metric1,metric2);
maxits = 100;
%  Now iterate on Bellman's equation and get the decision 
%   rules and the value function at the optimum         @
while (metric>1e-6 && iter<=maxits);
%     iter
     z1 = R1 + beta*v(:,1)*ones(1,nk);   % compute hypothetical value
     [tv,tdecis] = max(z1);              % [compute Tv given v
     tv = tv'; tdecis = tdecis';         % ,compute best choice of kprime at each value of k]
    metric1 = max(abs(tdecis-decis));    % check if decision rules have changed @
    metric2 = max(abs(v-tv)./(1+abs(v)));% check if value function has changed @
    metric = max(metric1,metric2);       % compute the max of the two convergence criteria @
    v=tv;                                % value function for next iteration @
    decis=tdecis;                        % decision rule - used only for convergence criteria computation @
    s = sprintf('Iteration: %d\t Metric1:    %d\t Metric2: %g',iter,metric1,metric2); disp(s);
    %pause
    iter=iter+1;
    
end

kdecis = kgrid(decis);                 % optimal policy function for kprime @
cdecis = Aval*kgrid.^alfa - kdecis + (1-delta)*kgrid;                              
                                       % implied consumption function @
idecis = kdecis - (1-delta)*kgrid;     % investment function

%   print out results
disp('Value Function Iterations for the deterministic growth model are completed.')
toc
disp('.')
disp('Number of iterations:'); 
disp(iter-1)
disp('Number of Grid points:'); 
disp(nk)
disp('PARAMETER VALUES')
disp('    alfa      beta      delta  ')
disp([alfa, beta, delta])
disp('Value of Productivity')
disp(Aval)
plot(kgrid,v);
title('Deterministic Growth Model: VALUE FUNCTION');
 
% pause

figure
plot(kgrid,kdecis);
title('Deterministic Growth Model: POLICY FUNCTION FOR KPRIME');
hold on
plot(kgrid,kgrid,'k')
plot(kgrid,kstar,'r.','MarkerSize',0.01)
hold off

% pause

figure
plot(kgrid,cdecis);
title('Deterministic Growth Model: Consumption')
hold on
plot(kgrid,cstar,'r.','MarkerSize',0.01)
hold off

% pause

figure
plot(kgrid,idecis);
title('Deterministic Growth Model: Investment')
hold on
plot(kgrid,istar,'r.','MarkerSize',0.01)
hold off