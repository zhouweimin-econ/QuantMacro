function [residual,mHistogramOptional,mAssetsPrimeOptional,mConsumptionOptional] = computeMCResidualHistogram(capital)

% Computes residual of market-clearing condition, using histogram approximation of distribution
% as in Young (2010); used to compute initial guess for parametric family
% 
% Inputs
%   (1) capital: candidate aggregate capital stock
%
% Outputs
%   (1) residual: residual of market clearing condition (to use in root finder)
%   (2) (optional) mHistogramOptional: mHistogram - histogram distribution
%   (3) (optional) mAssetsPrimeOptional: mAssetsPrime
%   (4) (optional) mConsumptionOptional: mConsumption
% 
% Thomas Winberry, July 26th, 2016

% Declare global variables
global bbeta ssigma aalpha ddelta eepsilonBar rrhoEpsilon ssigmaEpsilon aaBar aggEmployment mmu ttau mEpsilonTransition vEpsilonGrid ...
	nEpsilon nAssets nState epsilonMin epsilonMax assetsMin assetsMax ...
	vAssetsGridZeros vAssetsGrid vAssetsGridHistogram mEpsilonGrid mAssetsGrid ...
	vAssetsPoly vAssetsPolySquared mAssetsPolyHistogram w r mEpsilonPrimeGrid maxIterations tolerance dampening vAssetsPolyFine vAssetsGridFine ...
	mEpsilonGridFine mAssetsGridFine nAssetsFine nStateFine splineOpt
	
% Compute prices
r = aalpha * (capital ^ (aalpha - 1)) * (aggEmployment ^ (1 - aalpha)) - ddelta;
w = (capital ^ aalpha) * (1 - aalpha) * (aggEmployment ^ (-aalpha));

%----------------------------------------------------------------
% Compute individual decisions
%----------------------------------------------------------------

if splineOpt == 0	% approximate conditional expectation function using polynomials

	% Initialize coefficients using rule of thumb savings rule
	mGridInit = log(bbeta * (1 + r) * ((w * (mmu * (1 - mEpsilonGrid) + (1 - ttau) * mEpsilonGrid) + ...
		r * mAssetsGrid) .^ (-ssigma)));
	mCoefficients = zeros(nEpsilon,nAssets);
	for iEpsilon = 1:nEpsilon	% interpolate
		vCoefficients = sum(vAssetsPoly' .* (ones(nAssets,1) * mGridInit(iEpsilon,:)),2);
		mCoefficients(iEpsilon,:) = (vCoefficients ./ vAssetsPolySquared)';
	end

	% Iterate
	err = 100; iteration = 1;
	while err > tolerance && iteration <= maxIterations

		mCoefficientsNew = updateCoefficients_polynomials(mCoefficients);
		err = max(abs(mCoefficientsNew(:) - mCoefficients(:)));
		iteration = iteration + 1;
		mCoefficients = dampening * mCoefficients + (1 - dampening) * mCoefficientsNew;

	end
	
else	% approximate savings decision using linear splines

	% Initialize coefficients
	mAssetsPrime = mAssetsGrid;
	
	% Iterate
	err = 100; iteration = 1;
	while err > tolerance && iteration <= maxIterations
		mAssetsPrimeNew = updateCoefficients_splines(mAssetsPrime);
		err = max(abs(mAssetsPrimeNew(:) - mAssetsPrime(:)));
		iteration = iteration + 1;
		mAssetsPrime = dampening * mAssetsPrime + (1 - dampening) * mAssetsPrimeNew;
	end

end

%----------------------------------------------------------------
% Compute histogram approximation of stationary distribution
%----------------------------------------------------------------

%%%
% Compute policies over histogram grid	
%%%

if splineOpt == 0

	% Compute decision rules along fine grid
	mConditionalExpectation = exp(mCoefficients * vAssetsPolyFine');

	% Compute savings policy
	mAssetsPrimeStar = w * (mmu * (1 - mEpsilonGridFine) + (1 - ttau) * mEpsilonGridFine) + ...
		(1 + r) * mAssetsGridFine - (mConditionalExpectation .^ (-1 / ssigma));
	mAssetsPrimeFine = max(mAssetsPrimeStar,aaBar * ones(nEpsilon,nAssetsFine));
		
else	% linearly interpolate savings rule over fine grid

	% Compute weights
	[vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = computeLinearWeights(vAssetsGrid,vAssetsGridFine);
		
	% Linear interpolation
	mAssetsPrimeFine = mAssetsPrime(:,vIndicesBelow) .* repmat(vWeightBelow',nEpsilon,1) + ...
		mAssetsPrime(:,vIndicesAbove) .* repmat(vWeightAbove',nEpsilon,1);
	
end

% Compute consumption
mConsumptionFine = w * (mmu * (1 - mEpsilonGridFine) + (1 - ttau) * mEpsilonGridFine) + ...
	(1 + r) * mAssetsGridFine - mAssetsPrimeFine;
	
%%%
% Compute transition matrix associated with policy rules
%%%

% Compute weighting matrices
[vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = computeLinearWeights(vAssetsGridFine,mAssetsPrimeFine(:));
	
% Compute transition matrix for assets over full grid
mTransitionAbove = zeros(nStateFine,nAssetsFine);
mTransitionBelow = zeros(nStateFine,nAssetsFine);
for a = 1:nAssetsFine
	mTransitionBelow(vIndicesBelow == a,a) = vWeightBelow(vIndicesBelow == a);
	mTransitionAbove(vIndicesAbove == a,a) = vWeightAbove(vIndicesAbove == a);
end
mAssetsTransition = kron(mTransitionBelow + mTransitionAbove,ones(1,nEpsilon));

% Compute transition matrix for idiosyncratic shocks over full grid
mEpsilonTransitionHistogram = repmat(mEpsilonTransition,nAssetsFine);

% Compute full transition matrix
mTransition = sparse(mAssetsTransition .* mEpsilonTransitionHistogram);

% Compute invariant histogram by iteration
errHistogram = 100;	iterationHistogram = 0;
vHistogram = ones(nStateFine,1) ./ nStateFine;
while errHistogram > 1e-12 && iterationHistogram < 1e4
	
	vHistogramNew = mTransition' * vHistogram;
	errHistogram = max(abs(vHistogramNew - vHistogram)); iterationHistogram = iterationHistogram + 1;
	vHistogram = vHistogramNew;
	
end

% Expand histogram matrix
mHistogram = reshape(full(vHistogramNew),nEpsilon,nAssetsFine);

%----------------------------------------------------------------
% Return market clearing residual
%----------------------------------------------------------------

residual = capital - sum(vAssetsGridFine' .* (mHistogram(1,:) + mHistogram(2,:)));

if nargout > 1 

    mHistogramOptional = mHistogram;
	
    if nargout > 2
        mAssetsPrimeOptional = mAssetsPrimeFine;
        mConsumptionOptional = mConsumptionFine;
    end
	
end