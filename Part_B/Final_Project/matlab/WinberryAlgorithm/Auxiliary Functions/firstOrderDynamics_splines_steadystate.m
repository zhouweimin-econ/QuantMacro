function [ys,check] = firstOrderDynamics_polynomials_steadyState(ys,exo)

% Computes stationary equilibrium of the model for Dynare; format is required
% to be called by Dynare (follows example of NK_baseline.mod in Dynare examples)
%
% Thomas Winberry, July 26th, 2016

tStart = tic;
fprintf('\nComputing steady state...\n')

%----------------------------------------------------------------
% Call parameters (the next set of commands will overwrite some)
%----------------------------------------------------------------
setParameters;

%----------------------------------------------------------------
% Read in parameters from Dynare declaration
%----------------------------------------------------------------

% Initialize indicator
check = 0;

% Read parameters from Dynare
global M_ 

% Read out parameters to access them with their name
for iParameter = 1:M_.param_nbr
  paramname = deblank(M_.param_names(iParameter,:));
  eval(['global ' paramname]);
  eval([ paramname ' = M_.params(' int2str(iParameter) ');']);
end

%----------------------------------------------------------------
% Steady State
%----------------------------------------------------------------
displayOpt = 'off';       % 'iter-detailed' or 'off'
coreSteadyState;

% Prices
r = aalpha * (aggregateCapital ^ (aalpha - 1)) * (aggEmployment ^ (1 - aalpha)) - ddelta;
w = (aggregateCapital ^ aalpha) * (1 - aalpha) * (aggEmployment ^ (-aalpha));

%----------------------------------------------------------------
% Save values of steady state variables for Dynare (must be exactly
% as declared in Dynare)
%----------------------------------------------------------------

% Coefficients on conditional expectation function
for iEpsilon = 1 : nEpsilon
	for iAsset = 1 : nAssets
		eval(sprintf('coefficient_%d_%d = mCoefficients(iEpsilon,iAsset);',...
			iEpsilon,iAsset));
	end
end

% Moments and parameters of density away from borrowing constraint
for iEpsilon = 1 : nEpsilon
	for iMoment = 1 : nMeasure
		eval(sprintf('moment_%d_%d = mMoments(iEpsilon,iMoment);',iEpsilon,iMoment));
		eval(sprintf('measureCoefficient_%d_%d = mParameters(iEpsilon,iMoment+1);',iEpsilon,iMoment));
	end
end

% Mass at borrowing constraint
for iEpsilon = 1 : nEpsilon
	eval(sprintf('mHat_%d = mHat(iEpsilon);',iEpsilon));
end

% Other variables
aggregateCapital = (1 - mHat(1,1)) * (1 - aggEmployment) * mMoments(1,1) + (1 - mHat(2,1)) * aggEmployment * mMoments(2,1);
aggregateTFP = 0;
logAggregateOutput = log(exp(aggregateTFP) * (aggregateCapital ^ aalpha) * (aggEmployment ^ (1 - aalpha)));
logAggregateInvestment = log(ddelta * aggregateCapital);
logAggregateConsumption = log(exp(logAggregateOutput) - exp(logAggregateInvestment));
logWage = log(w);

% Save endogenous variables back into ys
for ii = 1 : M_.orig_endo_nbr;
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']); 
end

fprintf('... Done!  Elapsed time: %2.2f seconds \n\n',toc(tStart))
