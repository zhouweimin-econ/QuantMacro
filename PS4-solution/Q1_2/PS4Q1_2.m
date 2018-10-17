% =======================================================================
% Quant Macro PS-4  
% Weimin Zhou
% Due: 17, Oct, 2018
% This file is to add a labor choice, and redo item 1.
% =======================================================================
clear all; clf; cd '~/Desktop/PS4/Q1_2'; 
close all;% in order to save png
% set global parameter values

alfa   = 0.321;             % elasticity of output to capital 
beta   = 0.988;             % subjective discount factor 
delta  = .031;               % depreciation rate 
kappa = 5.24;
nu  = 2.0;   % adding labor choice parameters

% Compute deterministic steady-state variables:  
% Steady states condition stored at 'focss.m'
x0 = [5; 0.3; 0.2];
options = optimset('display','Iter','TolX',1e-8,'MaxIter',100);
myfun   = @(x) focss(x);
results1 = fsolve(myfun,x0,options);
istar = delta*results1(1);
kstar = results1(1);
hstar = results1(2);
cstar = results1(3);
disp('Steady States of the economic system')
fprintf('Kstar = ')
disp(results1(1))
fprintf('Hstar = ')
disp(results1(2))
fprintf('Cstar = ')
disp(results1(3))
fprintf('Istar = ')
disp(istar)
% cstar =(kstar^alfa)*(hstar^(1-alfa)) - delta*kstar;

%   form capital grid 
mink  =   0.5* kstar;            %  minimum value of the capital grid 
maxk  =   1.5* kstar;            % maximum value of the capital grid 
nk    =   500;              % choose the number of grid-points
kgrid =  (linspace(mink,maxk,nk))';  % the grid is now computed 

minh  =   0.5* hstar;           
maxh  =   1.5* hstar;           
nh    =   500;             
hgrid =  (linspace(minh,maxh,nh))';  

mOutput           = zeros(nk,nh);
mValueFunction    = ones(nk,1)*(log(cstar)-kappa*hstar^(1+1/nu)/(1+1/nu));
%mValueFunction =  zeros(nk,1);
mValueFunctionNew = zeros(nk,1);
mPolicyFunction   = zeros(nk,1);
mLaborFunction    = zeros(nk,1);

%==========
%I failed to use fsolve for obtaining consumption and labor instantaneously.
%==========
% Dyfun   = @(x) focdynamic1(x,k,invest);
% results = fsolve(Dyfun,x00,options); Or: 
% results = fsolve(Dyfun,x0,options);
% h      = results(1);    
% consu  = results(2); 

for nCapital=1:nk
    for nLabor=1:nh
            mOutput(nCapital,nLabor) = (kgrid(nCapital))^alfa*(hgrid(nLabor))^(1-alfa)+...
                kgrid(nCapital)*(1-delta);  % yt + (1-delta)*kt
    end
end

%% VFI
maxDifference = 100.0;
tol = 0.0001;
iter = 0;

tic
while (maxDifference>tol)  
        
        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;
        
        for nCapital = 1:nk
                        
            valueHighSoFar = -100.0;
            capitalChoice  = kgrid(1);
            gridLabor=1;
            
            for nCapitalNextPeriod = gridCapitalNextPeriod:nk
                
                valueHighSoFar1=-100.0;
                
                for nLabor = gridLabor:nh
                    
                consumption = mOutput(nCapital,nLabor)-kgrid(nCapitalNextPeriod); %c=yt+(1-delta)*k-k'
                valueProvisional =(log(consumption)-kappa*hgrid(nLabor)^(1+1/nu)/(1+1/nu))+...
                beta*mValueFunction(nCapitalNextPeriod);              
                
                if (valueProvisional>valueHighSoFar1)
                    valueHighSoFar1=valueProvisional;
                    %gridLabor=nLabor;
                    laborChoice=hgrid(nLabor);
                else
                    break;
                end
                
                end
                      
                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    capitalChoice = kgrid(nCapitalNextPeriod);
                    %gridCapitalNextPeriod = nCapitalNextPeriod;
                else
                    break; % We break when we have achieved the max
                end    
                  
            end
            
            mValueFunctionNew(nCapital) = valueHighSoFar;
            mPolicyFunction(nCapital) = capitalChoice;
            mLaborFunction(nCapital)=laborChoice;
        end
    
    maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
    mValueFunction = mValueFunctionNew;
    
    iter = iter+1;
    if (mod(iter,10)==0 || iter ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, maxDifference); 
    end
           
end



disp('Value Function Iterations for the deterministic growth model are completed.')
toc

disp('.')
disp('Number of iterations:'); 
disp(iter-1)
disp('Number of Grid points:'); 
disp(nk)
disp('PARAMETER VALUES')
disp('    alfa      beta      delta   kappa nu')
disp([alfa, beta, delta, kappa, nu])

decis_k = mPolicyFunction;
decis_h = mLaborFunction;
V       = mValueFunction;

subplot(2,1,1)
plot(kgrid,V);
title('Deterministic Growth Model with Labor Choice: VALUE FUNCTION');


subplot(2,1,2)
plot(kgrid,decis_k);
title('Deterministic Growth Model with Labor Choice: g(k)');
hold on
plot(kgrid,kgrid,'k')
plot(kgrid,kstar,'r.','MarkerSize',0.01)
hold off
saveas(gcf,'2a.png')

% ===================================
% Elapsed time is 216.576990 seconds.
%
% Number of iterations:
%   769
% Number of Grid points:
%   500
% PARAMETER VALUES
%    alfa      beta      delta     kappa     nu
%    0.3210    0.9880    0.0310    5.2400    2.0000
% ===================================

% Above I take into account monotonicity of the optimal decision
% notice that, it's two state variables now in computing
% hence, we need two variables to take into account our restrictions.
%% (f) 
% ===============
% below, I take into account Howard's improvement
% the rest are just copy-paste-modify task from above
% ===============
mOutput           = zeros(nk,nh);
mValueFunction    = ones(nk,1)*(log(cstar)-kappa*hstar^(1+1/nu)/(1+1/nu));
%mValueFunction =  zeros(nk,1);
mValueFunctionNew = zeros(nk,1);
mPolicyFunction   = zeros(nk,1);
mLaborFunction    = zeros(nk,1);

%==========
%I failed to use fsolve for obtaining consumption and labor instantaneously.
%==========
% Dyfun   = @(x) focdynamic1(x,k,invest);
% results = fsolve(Dyfun,x00,options); Or: 
% results = fsolve(Dyfun,x0,options);
% h      = results(1);    
% consu  = results(2); 

for nCapital=1:nk
    for nLabor=1:nh
            mOutput(nCapital,nLabor) = (kgrid(nCapital))^alfa*(hgrid(nLabor))^(1-alfa)+...
                kgrid(nCapital)*(1-delta);  % yt + (1-delta)*kt
    end
end

maxDifference = 100.0;
tol = 0.0001;
iter = 0;

tic
while (maxDifference>tol)  
        
        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;
        
        for nCapital = 1:nk
                        
            valueHighSoFar = -100.0;
            capitalChoice  = kgrid(1);
            gridLabor=1;
            
            for nCapitalNextPeriod = gridCapitalNextPeriod:nk
                
                valueHighSoFar1=-100.0;
                
                for nLabor = gridLabor:nh
                    
                consumption = mOutput(nCapital,nLabor)-kgrid(nCapitalNextPeriod); %c=yt+(1-delta)*k-k'
                valueProvisional =(log(consumption)-kappa*hgrid(nLabor)^(1+1/nu)/(1+1/nu))+...
                beta*mValueFunction(nCapitalNextPeriod);              
                
                if (valueProvisional>valueHighSoFar1)
                    valueHighSoFar1=valueProvisional;
                    % Howard's Improvement
                    gridLabor=nLabor; % update the value 
                    laborChoice=hgrid(nLabor);
                else
                    break;
                end
                
                end
                      
                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    capitalChoice = kgrid(nCapitalNextPeriod);
                    % Howard's Improvement
                    gridCapitalNextPeriod = nCapitalNextPeriod; % update the value
                else
                    break; % We break when we have achieved the max
                end
               
            end
            
            mValueFunctionNew(nCapital) = valueHighSoFar;
            mPolicyFunction(nCapital) = capitalChoice;
            mLaborFunction(nCapital)=laborChoice;
        end
    
    maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
    mValueFunction = mValueFunctionNew;
    
    iter = iter+1;
    if (mod(iter,10)==0 || iter ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, maxDifference); 
    end
           
end



disp('Value Function Iterations for the deterministic growth model are completed.')
toc

disp('.')
disp('Number of iterations:'); 
disp(iter-1)
disp('Number of Grid points:'); 
disp(nk)
disp('PARAMETER VALUES')
disp('    alfa      beta      delta   kappa nu')
disp([alfa, beta, delta, kappa, nu])

decis_k = mPolicyFunction;
decis_h = mLaborFunction;
V       = mValueFunction;

subplot(2,1,1)
plot(kgrid,V);
title('Deterministic Growth Model with Labor Choice: VALUE FUNCTION');


subplot(2,1,2)
plot(kgrid,decis_k);
title('Deterministic Growth Model with Labor Choice: g(k)');
hold on
plot(kgrid,kgrid,'k')
plot(kgrid,kstar,'r.','MarkerSize',0.01)
hold off
saveas(gcf,'2f.png')

% Value Function Iterations for the deterministic growth model are completed.
% Elapsed time is 14.298212 seconds.
% 
% Number of iterations:
%   769
% Number of Grid points:
%   500
% PARAMETER VALUES
%    alfa      beta      delta     kappa     nu
%    0.3210    0.9880    0.0310    5.2400    2.0000

% Howard Improvements decrease a lot of time for iterations.