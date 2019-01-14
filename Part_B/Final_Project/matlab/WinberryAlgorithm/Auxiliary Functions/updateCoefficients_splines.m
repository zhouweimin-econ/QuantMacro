function mAssetsPrimeNew = updateCoefficients_splines(mAssetsPrime)

% Updates linear splines approximating savings decision in steady state
% 
% Inputs
%   (1) mAssetsPrime: nEpsilon x nAssets matrix, storing previous iteration's decision along grid
%
% Outputs
%   (1) mAssetsPrimeNew: nEpsilon x nAssets matrix, storing updated decisions
% 
% Thomas Winberry and Alp Tuncay, July 26th, 2016

% Declare global variables
global bbeta ssigma aalpha ddelta eepsilonBar rrhoEpsilon ssigmaEpsilon aaBar aggEmployment mmu ttau mEpsilonTransition vEpsilonGrid ...
	nEpsilon nAssets nState epsilonMin epsilonMax assetsMin assetsMax nAssetsFine mEpsilonPrimeGridFine mEpsilonGridFine mConsumption ...
	vAssetsGrid vAssetsGridHistogram mEpsilonGrid mAssetsGrid vAssetsGridFine mAssetsGridFine...
	w r mEpsilonPrimeGrid vIndicesAbove vIndicesBelow vWeightAbove vWeightBelow

%---------------------------------------------------------------
% Compute savings policy next period (mAssetsPrimePrime)
%---------------------------------------------------------------

% Compute weights and indices for linear interpolation
[vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = computeLinearWeights(vAssetsGrid,mAssetsPrime(:));

% Compute interpolation for each realization of epsilon
mAssetsPrimePrime = mAssetsPrime(:,vIndicesBelow) .* repmat(vWeightBelow',nEpsilon,1) + ...
	mAssetsPrime(:,vIndicesAbove) .* repmat(vWeightAbove',nEpsilon,1);

%---------------------------------------------------------------
% Compute updated savings policy 
%---------------------------------------------------------------

% Compute consumption in next period using mAssetsPrimePrime and budget constraint
mConsumptionPrime = w * (mmu * (1 - mEpsilonPrimeGrid) + (1 - ttau) * mEpsilonPrimeGrid) + (1 + r) * ...
	repmat(mAssetsPrime(:)',nEpsilon,1) - mAssetsPrimePrime;
	
% Compute savings policy from Euler Equation(i.e. assuming not borrowing constrained)
mConsumption = reshape((sum(repmat(mEpsilonTransition',1,nAssets) .* mConsumptionPrime .^ (-ssigma),1) * ... 
	bbeta * (1 + r)) .^ (-1 / ssigma),nEpsilon,nAssets);
mAssetsPrimeNewStar = (w * (mmu * (1 - mEpsilonGrid) + (1 - ttau) * mEpsilonGrid) + ...
        (1 + r) * mAssetsGrid) - mConsumption;

% Enforce borrowing constraint
mAssetsPrimeNew = max(mAssetsPrimeNewStar,aaBar * ones(nEpsilon,nAssets));