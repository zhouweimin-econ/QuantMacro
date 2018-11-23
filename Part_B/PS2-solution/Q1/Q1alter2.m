%__________________________________________________________________________
%
% Quant Macro Part-2
% PS2-Q1: Solution of the model without aggregate risk
% Author: Weimin Zhou
% Date: Nov, 2018
%__________________________________________________________________________
clear;clc;clf;close all; cd '~/Desktop/PS2';  % in order to save png
disp('Experiment on using PS4-Q2.1 on Raul Part')
disp('Solving Aiyagari model with endogenous labor choice by VFI by generating labor grid and consumption grid separately around their own SS values')
disp('Results: no convergence, due to sensitivity of Steady State, and separate grids')
%% 0 transition compute
A= [ 1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0 ; ...
     0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1 ; ...
     0  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  1  0  0  0  1  0  0 ; ...
     0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0 ; ...
     1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0 5.6 0 -1  0  0  0  0  0 ; ...
    -1 0 28/3 0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
  .02 .48 .05 .45 0 0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0 0 .02 .48 .05 .45 0 0 0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0 ];

b= [7/8; 7/8; 7/8; 7/8; 1/8; 1/8; 1/8; 1/8; 7/24; 21/40; 0; 0; 0.02; 0.005; 0.05; 0.02];

pize = reshape(A^-1*b,4,4); 
% fix z=z_g, normalize the transition of good times:
pie_g0 = pize(1:2,1:2);
pie_b0 = pize(3:4,3:4);   
for i=1:2
    for j=1:2
        pie_g(i,j)=pie_g0(i,j)/sum(pie_g0(i,:));
    end
end
for i=1:2
    for j=1:2
        pie_b(i,j)=pie_b0(i,j)/sum(pie_b0(i,:));
    end
end
%% 1. Add employment shcoks: 1 and 0 with mean equals to one.
% plot value function which contains labor choice.


aalpha = 0.36;           % Elasticity of output w.r.t. capital
bbeta  = 95/100;          % Discount factor
ddelta = 2.5/1000;          % Depreciation

% parameters need to calibrate
kappa = 2.5;
nu = 2.0;

%Grid Points
nGridCapital = 1000;
nGridLabor= 10;

r0  = (1/bbeta-1);
w0  = (1-aalpha)*(aalpha/(r0+ddelta))^(aalpha/(1-aalpha));

% Productivity values
vProductivity = [0; 1]';


question = 1;
% Transition matrix
if question == 1
    mTransition = pie_g;
else
    mTransition = pie_b;
end
 

% Steady State

%x0 = [5; 0.3; 0.2];
%options = optimset('display','Iter','TolX',1e-8,'MaxIter',100);
%myfun   = @(x) focss(x);
%results1 = fsolve(myfun,x0,options);

laborSteadyState = 0.33; 
capitalSteadyState = ((1/bbeta-1+ddelta)/(aalpha*laborSteadyState^(1-aalpha)))^(1/(aalpha-1));
outputSteadyState = capitalSteadyState^aalpha*laborSteadyState^(1-aalpha);
consumptionSteadyState = outputSteadyState-capitalSteadyState*ddelta;

% We generate all state variables' grids based on steady state.
vGridCapital = linspace(0.1*capitalSteadyState,1.9*capitalSteadyState,nGridCapital);
vGridLabor=linspace(0.1*laborSteadyState,1.9*laborSteadyState,nGridLabor);
nGridProductivity = length(vProductivity);

% initializing
mOutput           = zeros(nGridCapital,nGridLabor,nGridProductivity);
mValueFunction    = zeros(nGridCapital,nGridProductivity);
%mValueFunction    = ones(nGridCapital,nGridProductivity)*(log(consumptionSteadyState)-kappa*laborSteadyState^(1+1/nu)/(1+1/nu));
mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
mPolicyFunction   = zeros(nGridCapital,nGridProductivity);
mLaborFunction    = zeros(nGridCapital,nGridProductivity);
expectedValueFunction = zeros(nGridCapital,nGridProductivity);

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
                consumption = (1+r0)*vGridCapital(nCapital) + vGridLabor(nLabor)*vProductivity(nProductivity)*w0 - vGridCapital(nCapitalNextPeriod) + vGridCapital(nCapital)*(1-ddelta);  
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
                consumption =(1+r0)*vGridCapital(nCapital) + mLaborFunction(nCapital,nProductivity)*vProductivity(nProductivity)*w0 - mPolicyFunction(nCapital,nProductivity);

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