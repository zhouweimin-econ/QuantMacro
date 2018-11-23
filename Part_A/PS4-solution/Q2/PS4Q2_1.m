% Weimin Zhou
% Due: 17, Oct, 2018
% This file is to add a labor choice, and redo item 1.
% =======================================================================
clear;clf;close all; cd '~/Desktop/PS4/Q2';  % in order to save png
disp('Solving deterministic version without labor choice by value function iterations')
%% 1. Add productivity shcoks: 1.01 and 1/1.01 with mean equals to one.
% plot value function which contains labor choice.
% ====================================

aalpha = 0.321;           % Elasticity of output w.r.t. capital
bbeta  = 0.988;          % Discount factor
ddelta = .031;          % Depreciation
kappa = 5.24;
nu = 2.0;

%Grid Points
nGridCapital = 1000;
nGridLabor= 1000;

% Productivity values
vProductivity = [1.01; 1/1.01;]';

% Transition matrix
mTransition   = [1/5, 0;
                 0, 4/5];
             
% Steady State
laborSteadyState = 0.2890; %FROM PSQ1_2.m
capitalSteadyState = ((1/bbeta-1+ddelta)/(aalpha*laborSteadyState^(1-aalpha)))^(1/(aalpha-1));
outputSteadyState = capitalSteadyState^aalpha*laborSteadyState^(1-aalpha);
consumptionSteadyState = outputSteadyState-capitalSteadyState*ddelta;

% We generate the grid of capital
vGridCapital = linspace(0.1*capitalSteadyState,1.9*capitalSteadyState,nGridCapital);
vGridLabor=linspace(0.1*laborSteadyState,1.9*laborSteadyState,nGridLabor);
nGridProductivity = length(vProductivity);

% initializing
mOutput           = zeros(nGridCapital,nGridLabor,nGridProductivity);
mValueFunction  = zeros(nGridCapital,nGridProductivity);
%mValueFunction    = ones(nGridCapital,nGridProductivity)*(log(consumptionSteadyState)-kappa*laborSteadyState^(1+1/nu)/(1+1/nu));
mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
mPolicyFunction   = zeros(nGridCapital,nGridProductivity);
mLaborFunction    = zeros(nGridCapital,nGridProductivity);
expectedValueFunction = zeros(nGridCapital,nGridProductivity);

for nCapital=1:nGridCapital
    for nLabor=1:nGridLabor
        for nProductivity=1:nGridProductivity
            mOutput(nCapital,nLabor,nProductivity) = (vGridCapital(nCapital))^aalpha*...
                vProductivity(nProductivity)*(vGridLabor(nLabor))^(1-aalpha);
        end
    end
end

% iteration
tic 
maxDifference = 100.0;
tolerance = 0.00001;
iteration = 0;

while (maxDifference>tolerance)  
%while (iteration<210)
    
    expectedValueFunction = mValueFunction*mTransition';
    
    if (mod(iteration,10)==0 || iteration ==1)
        
    for nProductivity = 1:nGridProductivity
        
        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;
        
        for nCapital = 1:nGridCapital
                        
            valueHighSoFar = -1000.0;
            capitalChoice  = vGridCapital(1);
            gridLabor=1;
            
            for nCapitalNextPeriod = gridCapitalNextPeriod:nGridCapital
                
                valueHighSoFar1=-1000.0;
                
                for nLabor = gridLabor:nGridLabor
                    
                consumption = mOutput(nCapital,nLabor,nProductivity)-vGridCapital(nCapitalNextPeriod)+vGridCapital(nCapital)*(1-ddelta);
                valueProvisional =log(consumption)-kappa*vGridLabor(nLabor)^(1+1/nu)/(1+1/nu)+...
                bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity);              
                
                if (valueProvisional>valueHighSoFar1)
                    valueHighSoFar1=valueProvisional;
                    gridLabor=nLabor;
                    laborChoice=vGridLabor(nLabor);
                else
                    break;
                end
                
                end
                      
                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    capitalChoice = vGridCapital(nCapitalNextPeriod);
                    gridCapitalNextPeriod = nCapitalNextPeriod;
                else
                    break; % Imposed improvements to speed up
                end    
                  
            end
            
            mValueFunctionNew(nCapital,nProductivity) = valueHighSoFar;
            mPolicyFunction(nCapital,nProductivity) = capitalChoice;
            mLaborFunction(nCapital,nProductivity)=laborChoice;
            mGridPolicyFunction(nCapital,nProductivity)=gridCapitalNextPeriod;
        end
    end
    
    else
             for nProductivity = 1:nGridProductivity
                
                 for nCapital=1:nGridCapital
                
                consumption = (vGridCapital(nCapital))^aalpha*vProductivity(nProductivity)*...
                    (mLaborFunction(nCapital,nProductivity))^(1-aalpha)-mPolicyFunction(nCapital,nProductivity)+...
                    vGridCapital(nCapital)*(1-ddelta); % ct = yt - kp + (1-delta)*kt
                
                gridCapitalNextPeriod1=mGridPolicyFunction(nCapital,nProductivity);
                
                mValueFunctionNew(nCapital,nProductivity)=log(consumption)-...
                    kappa*mLaborFunction(nCapital,nProductivity)^(1+1/nu)/(1+1/nu)+ ...
                    bbeta*expectedValueFunction(gridCapitalNextPeriod1,nProductivity);
                 end
            end
    end  
    
    maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
    mValueFunction = mValueFunctionNew;
    
    iteration = iteration+1;
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
    end
           
end

disp('Value Function Iterations for the growth model are completed.')
toc
disp('.')
disp('Number of iterations:'); 
disp(iteration-1)
%% Plot
figure(1)
subplot(3,1,1)
plot(vGridCapital,mValueFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Value Function')

subplot(3,1,2)
plot(vGridCapital,mPolicyFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Policy Function')

subplot(3,1,3)
plot(vGridCapital,mLaborFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Labor Function')
saveas(gcf,'21.png')
%% 2. Simulated Economy
% tricks I use: make transition function more diverse, and we could see the
% mean is still close to one, but the graph can give difference between two
% different shocks.

% compare the logged and HP-filtered variance of y, c, i, k, w, w*h/y and h in the model, 
% their persistence, and their co-movement with output. 
% Compare them to the data. Is this a good model to understand business cycle fluctuations?


y = mPolicyFunction.^aalpha.*mLaborFunction.^(1-aalpha);
k = mPolicyFunction;
h = mLaborFunction;

for nProductivity = 1:nGridProductivity
                
    for nCapital=1:nGridCapital
                
      consumption(nCapital,nProductivity) = (vGridCapital(nCapital))^aalpha*vProductivity(nProductivity)*...
                    (mLaborFunction(nCapital,nProductivity))^(1-aalpha)-mPolicyFunction(nCapital,nProductivity)+...
                    vGridCapital(nCapital)*(1-ddelta); % ct = yt - kp + (1-delta)*kt


   end
end

c = consumption;
invest = y - c;
logy = log(y);
var(logy); %logged variance
logk = log(k);
var(logk); %logged variance
logh = log(h);
var(logh); %logged variance
logc = log(c);
var(logc); %logged variance
loginvest = log(invest);
var(loginvest); %logged variance
disp('---------------')
disp('logged variance of some variables in the model for good shock')
disp('    output       capital    hours   consumption      investment  ')
disp([var(logy(:,1)), var(logk(:,1)), var(logh(:,1)), var(logc(:,1)), var(loginvest(:,1))])
disp('logged variance of some variables in the model for bad shock')
disp('    output       capital    hours   consumption      investment  ')
disp([var(logy(:,2)), var(logk(:,2)), var(logh(:,2)), var(logc(:,2)), var(loginvest(:,2))])
disp('---------------')


plot(y)
title('i.e. Original simulated output data')
saveas(gcf,'22y.png')

[gnptrend16, Simulatedy16] = hpfilter(y(:,1),1600);
[gnptrend64, Simulatedy64] = hpfilter(y(:,1),6400);
[gnptrendinf, Simulatedyinf] = hpfilter(y(:,1),Inf);
subplot(2,1,1)
plot(Simulatedy16,'b');
hold all
plot(Simulatedyinf - Simulatedy16,'r');
title('\bf HP filtered Output for Good Shocks');
legend('Cyclical output','Difference');
hold off

[gnptrend16, Simulatedy16] = hpfilter(y(:,2),1600);
[gnptrend64, Simulatedy64] = hpfilter(y(:,2),6400);
[gnptrendinf, Simulatedyinf] = hpfilter(y(:,2),Inf);
subplot(2,1,2)
plot(Simulatedy16,'b');
hold all
plot(Simulatedyinf - Simulatedy16,'r');
title('\bf HP filtered Output for Bad Shocks');
legend('Cyclical output','Difference');
hold off

% compute variance of HP filtered data:
[gnptrend16, Simy] = hpfilter(y,1600);
[gnptrend16, Simk] = hpfilter(k,1600);
[gnptrend16, Simc] = hpfilter(c,1600);
[gnptrend16, Simh] = hpfilter(h,1600);
[gnptrend16, Siminvest] = hpfilter(invest,1600);
disp('---------------')
disp('HP filter variance of some variables in the model for good shock')
disp('    output       capital    hours   consumption      investment  ')
disp([var(Simy(:,1)), var(Simk(:,1)), var(Simh(:,1)), var(Simc(:,1)), var(Siminvest(:,1))])
disp('HP filter variance of some variables in the model for bad shock')
disp('    output       capital    hours   consumption      investment  ')
disp([var(Simy(:,2)), var(Simk(:,2)), var(Simh(:,2)), var(Simc(:,2)), var(Siminvest(:,2))])
disp('---------------')


model = arima(1,0,0); %set AR(1) process for obtaining persistence.
estmdly1 = estimate(model, y(:,1));
pery1 = estmdly1.AR;
estmdly1 = estimate(model, y(:,2));
pery2 = estmdly1.AR;

estmdly1 = estimate(model, k(:,1));
perk1 = estmdly1.AR;
estmdly1 = estimate(model, k(:,2));
perk2 = estmdly1.AR;

estmdly1 = estimate(model, h(:,1));
perh1 = estmdly1.AR;
estmdly1 = estimate(model, h(:,2));
perh2 = estmdly1.AR;

estmdly1 = estimate(model, c(:,1));
perc1 = estmdly1.AR;
estmdly1 = estimate(model, c(:,2));
perc2 = estmdly1.AR;

estmdly1 = estimate(model, invest(:,1));
perinvest1 = estmdly1.AR;
estmdly1 = estimate(model, invest(:,2));
perinvest2 = estmdly1.AR;

disp('---------------')
disp('Persistence of some variables in the model for good shock')
disp('    output       capital    hours   consumption      investment  ')
disp([pery1, perk1, perh1, perc1, perinvest1])
disp('Persistence of some variables in the model for bad shock')
disp('    output       capital    hours   consumption      investment  ')
disp([pery2, perk2, perh2, perc2, perinvest2])
disp('---------------')

fid = fopen('test.txt', 'wt');
fprintf(fid, '%f, %f\n', c(:,1), y(:,1));

A = [y(:,1) k(:,1) h(:,1) c(:,1) invest(:,1)];
R = corrcoef(A)

B = [y(:,2) k(:,2) h(:,2) c(:,2) invest(:,2)];
RB = corrcoef(B)

[contribution]=decomp(y(:,1),k(:,1),h(:,1),c(:,1),invest(:,1));

[contributionbad]=decomp(y(:,2),k(:,2),h(:,2),c(:,2),invest(:,2));

%% for the real data, I use US.data from bea.gov
% MODIFY capital and hours from Table 1.12. National Income by Type of Income
% Here, due to the time limit, I only check the invest coonsumption and
% output data


data = xlsread('download.xls','Sheet1','A1:C284');
output = data(1:284,1);
consumption = data(1:284,2);
investment = data(1:284,3);
% do the same procedure as above. 

disp('---------------')
disp('logged variance of some variables of US data')
disp('    output          consumption           investment  ')
disp([var(log(output)), var(log(consumption)), var(log(investment))])
disp('---------------')



[gnptrend16, output16] = hpfilter(output,1600);
[gnptrendinf, outputinf] = hpfilter(output,Inf);
plot(output16,'b');
hold all
plot(outputinf - output16,'r');
title('\bf HP filtered real US Output from 1947 to 2018');
legend('Cyclical output','Difference');
hold off


% compute variance of HP filtered data:
[gnptrend16, Simy] = hpfilter(output,1600);
[gnptrend16, Simc] = hpfilter(consumption,1600);
[gnptrend16, Siminvest] = hpfilter(investment,1600);
disp('---------------')
disp('HP filter variance of some variables in US data')
disp('    output        consumption      investment  ')
disp([var(Simy), var(Simc), var(Siminvest)])
disp('---------------')


model = arima(1,0,0); %set AR(1) process for obtaining persistence.
estmdly1 = estimate(model, output);
pery1 = estmdly1.AR;

estmdly1 = estimate(model, consumption);
perc1 = estmdly1.AR;

estmdly1 = estimate(model, investment);
perinvest1 = estmdly1.AR;

disp('---------------')
disp('Persistence of some variables in US data')
disp('    output     consumption      investment  ')
disp([pery1, perc1, perinvest1])
disp('---------------')


A = [output consumption investment];
R = corrcoef(A)


[contribution]=decomp1(output,consumption,investment);
%% 3.  impulse response function of output, consumption, investment, capital, wages, labor share and hours to a productivity shock

% execute a Dynare file 5
% specify parameters 7
beta = 0.988;
alpha = 0.321;
sigma = 1;  %log-utility of consumption
sigmae = 0.01;  %variance of shocks
rho = 0.9;   %persistence of shocks
delta = 0.031;
kappa = 5.24;
nu = 2.0;
% without labor
save param_nc alpha beta delta rho sigma sigmae
dynare basic_nc_dynare_alt fast noclearall nolog

%% with labor
cd('~/Desktop/PS4');
save param_nc alpha beta delta rho sigma sigmae kappa nu
dynare basic_nc_dynare_alt fast noclearall nolog
