% =======================================================================
% Quant Macro PS-4  
% Weimin Zhou
% Due: 17, Oct, 2018
% This file is to add a labor choice, and redo item 1.
% =======================================================================
clear;clf;close all;
cd '~/Desktop/PS4'  % in order to save png

disp('Solving deterministic version without labor choice by value function iterations')

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
nk   =   2000;              % choose the number of grid-points
                          
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
tic
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
    %s = sprintf('Iteration: %d\t Metric1:    %d\t Metric2: %g',iter,metric1,metric2); disp(s);
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
saveas(gcf,'a1.png')

% pause

figure
plot(kgrid,cdecis);
title('Deterministic Growth Model: Consumption')
hold on
plot(kgrid,cstar,'r.','MarkerSize',0.01)
hold off
saveas(gcf,'a2.png')
% pause

figure
plot(kgrid,idecis);
title('Deterministic Growth Model: Investment')
hold on
plot(kgrid,istar,'r.','MarkerSize',0.01)
hold off
saveas(gcf,'a3.png')
toc 
% 10.805486 seconds.
%%  (b) taking into account monotonicity of the optimal decision rule.
R1=-10000*ones(nk,nk);    % initial values 
i=1;
tic 
while i<=nk;
    k=kgrid(i);
  
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

ksofar = -100;

while (metric>1e-6 && iter<=maxits);

     z1 = R1 + beta*v(:,1)*ones(1,nk);   
     [tv,tdecis] = max(z1);              
     tv = tv'; tdecis = tdecis';         
    metric1 = max(abs(tdecis-decis));    
    metric2 = max(abs(v-tv)./(1+abs(v)));
    metric = max(metric1,metric2);       
    v=tv;                                
    decis=tdecis;
    % ==== Monotonicity 
    if decis>ksofar % takes an initial -100;
        ksofar = decis;
    else
        break
    end
    % s = sprintf('Iteration: %d\t Metric1:    %d\t Metric2: %g',iter,metric1,metric2); disp(s);
    %pause
    iter=iter+1;
    
end
disp('Value Function Iterations for the deterministic growth model are completed.')
toc
disp('.')
disp('Number of iterations:'); 
disp(iter-1)
plot(kgrid,v);
saveas(gcf,'b.png')  % faster than (a) by checking tic-toc time
% 0.587886 seconds.
%% (c) taking into account concavity of the value function.

R1=-10000*ones(nk,nk);    % initial values 
i=1;
tic 
while i<=nk;
    k=kgrid(i);
  
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

z0 = -100;

while (metric>1e-6 && iter<=maxits);

     z1 = R1 + beta*v(:,1)*ones(1,nk);   
     [tv,tdecis] = max(z1);              
     tv = tv'; tdecis = tdecis';         
    metric1 = max(abs(tdecis-decis));    
    metric2 = max(abs(v-tv)./(1+abs(v)));
    metric = max(metric1,metric2);       
    v=tv;                                
    decis=tdecis;
    % ==== Concavity 
    if z1 > z0 ;
        z0 =z1;  
    else 
        break
    end
    % s = sprintf('Iteration: %d\t Metric1:    %d\t Metric2: %g',iter,metric1,metric2); disp(s);
    %pause
    iter=iter+1;
    
end
disp('Value Function Iterations for the deterministic growth model are completed.')
toc
disp('.')
disp('Number of iterations:'); 
disp(iter-1)

plot(kgrid,v);
saveas(gcf,'c.png')  % faster than (a) by checking tic-toc time
% 0.590535 seconds.

%% (d) taking into account local search on the decision rule
R1=-10000*ones(nk,nk);    % initial values 
i=1;
tic 
while i<=nk;
    k=kgrid(i);
  
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

ksofar = -100;

while (metric>1e-6 && iter<=maxits);

     z1 = R1 + beta*v(:,1)*ones(1,nk);   
     [tv,tdecis] = max(z1);              
     tv = tv'; tdecis = tdecis';         
    metric1 = max(abs(tdecis-decis));    
    metric2 = max(abs(v-tv)./(1+abs(v)));
    metric = max(metric1,metric2);       
    v=tv;                                
    decis=tdecis;
    % ==== Concavity 
    if decis - ksofar < 0.01; 
        break  
    else 
        ksofar = decis;
    end
    % s = sprintf('Iteration: %d\t Metric1:    %d\t Metric2: %g',iter,metric1,metric2); disp(s);
    %pause
    iter=iter+1;
    
end
disp('Value Function Iterations for the deterministic growth model are completed.')
toc
disp('.')
disp('Number of iterations:'); 
disp(iter-1)

plot(kgrid,v);
saveas(gcf,'d.png')  % faster than (a) by checking tic-toc time
% 0.587183 seconds
%% (e) taking into account both concavity of the value function and monotonicity of the decision rule

R1=-10000*ones(nk,nk);    % initial values 
i=1;
tic 
while i<=nk;
    k=kgrid(i);
  
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

z0 = -100;
ksofar = -100;

while (metric>1e-6 && iter<=maxits);

     z1 = R1 + beta*v(:,1)*ones(1,nk);   
     [tv,tdecis] = max(z1);              
     tv = tv'; tdecis = tdecis';         
    metric1 = max(abs(tdecis-decis));    
    metric2 = max(abs(v-tv)./(1+abs(v)));
    metric = max(metric1,metric2);       
    v=tv;                                
    decis=tdecis;
    % ==== Concavity 
    if z1 > z0;
        z0 =z1;  
    else 
        break
    end
    
    if decis>ksofar % takes an initial -100;
        ksofar = decis;
    else
        break
    end
    s = sprintf('Iteration: %d\t Metric1:    %d\t Metric2: %g',iter,metric1,metric2); disp(s);
    %pause
    iter=iter+1;
    
end
disp('Value Function Iterations for the deterministic growth model are completed.')
toc
disp('.')
disp('Number of iterations:'); 
disp(iter-1)

plot(kgrid,v);
saveas(gcf,'e.png')  % faster than (a) by checking tic-toc time
% 0.597594 seconds.
% slightly longer because in Matlab there're two more loops.
%% (f)  Howard?s policy iterations waiting until converged to solve the problem
v       = zeros(nk,1);                 % the value function: row = kgrid 
decis   = zeros(nk,1);                 % the policy function: row = kgrid 
iter    = 0;                           % counter 
metric1 = 10;                          % metric 1 for convergence 
metric2 = 10;                          % metric 2 for convergence

metric = max(metric1,metric2);
maxits = 100;

R1=-10000*ones(nk,nk);    % initial values 
i=1;

tic 
while i<=nk;
    k=kgrid(i);
  
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

while (metric>1e-6 && iter<=maxits);

     z1 = R1 + beta*v(:,1)*ones(1,nk);   
     [tv,tdecis] = max(z1);              
     tv = tv'; tdecis = tdecis';         
    metric1 = max(abs(tdecis-decis));    
    metric2 = max(abs(v-tv)./(1+abs(v)));
    metric = max(metric1,metric2);       
    v=tv;                                
    decis=tdecis;
    
    % Howard Improvement: Update the values:
while i<=nk;
    k=kgrid(i);
  
    j=1;    
    while j<=nk; 
        kprime = kgrid(decis(j));
        invest = kprime - (1-delta)*k; 
        consu = Aval*k^(alfa) - invest;
        if consu > 0; 
           R1(j,i)= log(consu);
        end
      j=j+1;
        
    end
    i=i+1;
end
    
    iter=iter+1;
    
end
disp('Value Function Iterations for the deterministic growth model are completed.')
toc
disp('.')
disp('Number of iterations:'); 
disp(iter-1)

plot(kgrid,v);
saveas(gcf,'f1.png')  % faster than (a) by checking tic-toc time
% 2.130719 seconds
% slightly longer because in Matlab there're two more loops.
%% (g) 
nk = 2000;
v       = zeros(nk,1);                 % the value function: row = kgrid 
decis   = zeros(nk,1);                 % the policy function: row = kgrid 
iter    = 0;                           % counter 
metric1 = 10;                          % metric 1 for convergence 
metric2 = 10;                          % metric 2 for convergence

metric = max(metric1,metric2);
%maxits = 100;

R1=-1000*ones(nk,nk);    % initial values 
i=1;



while i<=nk;
    k=kgrid(i);
  
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

tic 
for maxits = [5 10 20 50];
    
    while (metric>1e-6 && iter<=maxits);

        z1 = R1 + beta*v(:,1)*ones(1,nk);   
        [tv,tdecis] = max(z1);              
        %tv = tv';  %no used anymore 
        tdecis = tdecis';         

        % Only Policy Iterations: 
        while i<=nk;
            k=kgrid(i);
  
            j=1;    
            while j<=nk; 
                kprime = kgrid(tdecis(j));
                invest = kprime - (1-delta)*k; 
                consu = Aval*k^(alfa) - invest;
                if consu > 0; 
                    R1(j,i)= log(consu);
                end
                j=j+1;    
            end
        i=i+1;
        end
    
   Q = sparse(nk,nk);
   for i = 1:nk;
       Q(i,tdecis(i))=1;
   end
   v1 = R1./(speye(nk)-beta*Q);
   %v1 = R1./(1-beta);
   metric = max(abs(tdecis-decis));     
   v=v1;                                
   decis=tdecis;
   iter=iter+1;
          
end
figure
plot(kgrid,v);
end

disp('Value Function Iterations for the deterministic growth model are completed.')
toc
%saveas(gcf,'f1.png')  % faster than (a) by checking tic-toc time
% 2.130719 seconds

%% (g)
nbk     = 1000;
crit    = 1;
epsi    = 1e-6;
ks      = ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));
dev     = 0.7;
kmin    = (1-dev)*ks;
kmax    = (1+dev)*ks;
kgrid   = linspace(kmin,kmax,nbk)'; % builds the grid
v       = zeros(nbk,1);
kp0     = kgrid;
dr      = zeros(nbk,1);
%
% Main loop
%
while crit>epsi;
for i=1:nbk %
% value function
% initial guess on k(t+1)
% decision rule (will contain indices)
% compute indexes for which consumption is positive

  imax    = min(floor((kgrid(i)^alpha+(1-delta)*kgrid(i)-kmin)/dev)+1,nbk);
  %
  % consumption and utility
  %
  c       = kgrid(i)^alpha+(1-delta)*kgrid(i)-kgrid(1:imax);
  util    = log(c);
  %
  % find new policy rule
  %
  [v1,dr(i)]= max(util+beta*v(1:imax));
end 
  %
  % decision rules
  %
  kp  = kgrid(dr);
  c   = kgrid.^alpha+(1-delta)*kgrid-kp;
  %
  % update the value
  %
  util= log(c);
  Q   = sparse(nbk,nbk);
    for i=1:nbk;
        Q(i,dr(i)) = 1;
    end

Tv = (speye(nbk)-beta*Q)\util; 
crit= max(abs(kp-kp0));
v =Tv;
kp0 = kp;
end

plot(kp,v)
saveas(gcf,'g.png')