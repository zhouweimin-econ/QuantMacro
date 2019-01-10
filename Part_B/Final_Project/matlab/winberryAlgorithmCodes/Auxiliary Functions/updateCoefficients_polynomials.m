function mCoefficientsNew = updateCoefficients_polynomials(mCoefficients)

% Updates polynomial coefficients approximating the conditional expectation function in steady state
% 
% Inputs
%   (1) mCoefficients: nEpsilon x nAssets matrix, storing previous iteration's coefficients
%
% Outputs
%   (1) mCoefficientsNew: nEpsilon x nAssets matrix, storing updated coefficients
% 
% Thomas Winberry, January 19, 2016

% Declare global variables
global bbeta ssigma aalpha ddelta eepsilonBar rrhoEpsilon ssigmaEpsilon aaBar aggEmployment mmu ttau mEpsilonTransition vEpsilonGrid ...
	nEpsilon nAssets nState epsilonMin epsilonMax assetsMin assetsMax ...
	vAssetsGridZeros vAssetsGrid vAssetsGridHistogram mEpsilonGrid mAssetsGrid ...
	vAssetsGridHistogramZeros vAssetsPoly vAssetsPolySquared mAssetsPolyHistogram w r mEpsilonPrimeGrid

%---------------------------------------------------------------
% Compute current period's savings policy function
%---------------------------------------------------------------

% Compute conditional expectation
mConditionalExpectation = exp(mCoefficients * vAssetsPoly');

% Compute target saving
mAssetsPrimeStar = w * (mmu * (1 - mEpsilonGrid) + (1 - ttau) * mEpsilonGrid) + (1 + r) * mAssetsGrid - ...
	(mConditionalExpectation .^ (-1 / ssigma));

% Compute actual saving
mAssetsPrime = max(mAssetsPrimeStar,aaBar * ones(nEpsilon,nAssets));
mAssetsPrimeGrid = repmat(reshape(mAssetsPrime,1,nState),[nEpsilon 1]);

% Compute next period's polynomials
mAssetsPrimeZeros = scaleDown(mAssetsPrime,assetsMin,assetsMax);
mPolyAssetsPrime = computeChebyshev(nAssets,reshape(mAssetsPrimeZeros,nState,1));

%---------------------------------------------------------------
% Compute next period's savings policy function
%---------------------------------------------------------------

% Compute conditional expectation
mConditionalExpectationPrime = exp(mCoefficients * mPolyAssetsPrime');

% Compute target saving
mAssetsPrimePrimeStar = w * (mmu * (1 - mEpsilonPrimeGrid) + (1 - ttau) * mEpsilonPrimeGrid) + (1 + r) * mAssetsPrimeGrid - ...
	(mConditionalExpectationPrime .^ (-1 / ssigma));

% Compute actual savings
mAssetsPrimePrimeGrid = max(mAssetsPrimePrimeStar,aaBar*ones(nEpsilon,nEpsilon*nAssets));

%---------------------------------------------------------------
% Update conditional expectation function
%---------------------------------------------------------------

% Compute new conditional expectation function
mConsumptionPrime = w * (mmu * (1 - mEpsilonPrimeGrid) + (1 - ttau) * mEpsilonPrimeGrid) + (1 + r) * ...
	mAssetsPrimeGrid - mAssetsPrimePrimeGrid;
aConditionalExpectationTilde = reshape(bbeta * mEpsilonTransition * ((1 + r) * (mConsumptionPrime .^ (-ssigma))),...
	nEpsilon,nEpsilon,nAssets);

% Extract the relevant entries
mConditionalExpectation = zeros(nEpsilon,nAssets);
for iEpsilon = 1:nEpsilon
	mConditionalExpectation(iEpsilon,:) = aConditionalExpectationTilde(iEpsilon,iEpsilon,:);
end

% Update the coefficients
mCoefficientsNew = zeros(nEpsilon,nAssets);
for iEpsilon = 1:nEpsilon
	vCoefficients = sum(vAssetsPoly' .* (ones(nAssets,1) * log(mConditionalExpectation(iEpsilon,:))),2);
	mCoefficientsNew(iEpsilon,:) = (vCoefficients ./ vAssetsPolySquared)';
end

