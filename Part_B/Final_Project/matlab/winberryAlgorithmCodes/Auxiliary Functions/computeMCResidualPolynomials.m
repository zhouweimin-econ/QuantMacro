function [residual,mCoefficientsOptional,mParametersOptional,mMomentsOptional,mHatOptional] = ...
	computeMCResidualPolynomials(capital,mMoments,aGridMoments,mHat)

% Computes residual of market-clearing condition, parametric family to approximate distribution
% 
% Inputs
%   (1) capital: candidate aggregate capital stock
%	(2) mMoments: intial guess of moments of distribution (nEpsilon x nMeasure)
%	(3) aGridMoments: grid of centralized moments for computing PDF, corresponding to mMoments
%		(nEpsilon x nAssetsQuadrature x nMoments)
%	(4) mHat: initial guess of mass at borrowing constraint
%
% Outputs
%   (1) residual: residual of market clearing condition
%   (2) (optional) mCoefficientsOptional: coefficients on individual decisions (splines or polynomials)
%   (3) (optional) mParametersOptional: parameters of density away from borrowing constraint
%   (4) (optional) mMomentsOptional: moments of density away from borrowing constraint
%	(5) (optional) mHatOptional: mass at borrowing constraint
% 
% Thomas Winberry, July 26th, 2016

% Declare global variables
global bbeta ssigma aalpha ddelta eepsilonBar rrhoEpsilon ssigmaEpsilon aaBar aggEmployment mmu ttau mEpsilonTransition vEpsilonGrid ...
	nEpsilon nAssets nState epsilonMin epsilonMax assetsMin assetsMax ...
	vAssetsGridZeros vAssetsGrid vAssetsGridHistogram mEpsilonGrid mAssetsGrid ...
	vAssetsPoly vAssetsPolySquared mAssetsPolyHistogram w r mEpsilonPrimeGrid maxIterations tolerance dampening vAssetsPolyFine vAssetsGridFine ...
	mEpsilonGridFine mAssetsGridFine nAssetsFine nStateFine ...
	vAssetsPolyQuadrature vAssetsGridQuadrature mEpsilonGridQuadrature mAssetsGridQuadrature nAssetsQuadrature ...
	nStateQuadrature vQuadratureWeights vEpsilonInvariant nMeasure splineOpt vAssetsPolyBC
	
% Compute prices
r = aalpha * (capital ^ (aalpha - 1)) * (aggEmployment ^ (1 - aalpha)) - ddelta;
w = (capital ^ aalpha) * (1 - aalpha) * (aggEmployment ^ (-aalpha));

%----------------------------------------------------------------
% Set error tolerance & max iteration depending on use
%----------------------------------------------------------------

if nargout == 1 
    err1 = tolerance;
    err2 = 1e-4;
    tol2 = 200;
elseif nargout > 1 
    err1 = 1e-8;
    err2 = 1e-6;
    tol2 = 500;
else
    error('Check # of inputs & outputs for computeMCResidualPolynomials');
end

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
	
	mCoefficientsOptional = mCoefficients;
	
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
	
	mCoefficientsOptional = mAssetsPrime;

end

%----------------------------------------------------------------
% Compute policies over quadrature grid for integration
%----------------------------------------------------------------

if splineOpt == 0

	% Compute conditional expectation
	mConditionalExpectation = exp(mCoefficients * vAssetsPolyQuadrature');

	% Compute savings policy
	mAssetsPrimeStar = w * (mmu * (1 - mEpsilonGridQuadrature) + (1 - ttau) * mEpsilonGridQuadrature) + ...
		(1 + r) * mAssetsGridQuadrature - (mConditionalExpectation .^ (-1 / ssigma));
	mAssetsPrimeQuadrature = max(mAssetsPrimeStar,aaBar * ones(nEpsilon,nAssetsQuadrature));
		
else

	% Compute weights
	[vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = computeLinearWeights(vAssetsGrid,vAssetsGridQuadrature);
		
	% Linear interpolation
	mAssetsPrimeQuadrature = mAssetsPrime(:,vIndicesBelow) .* repmat(vWeightBelow',nEpsilon,1) + ...
		mAssetsPrime(:,vIndicesAbove) .* repmat(vWeightAbove',nEpsilon,1);
	
end

%----------------------------------------------------------------
% Compute policies at borrowing constraint for integration
%----------------------------------------------------------------

if splineOpt == 0

	% Compute conditional expectation
	mConditionalExpectation = exp(mCoefficients * vAssetsPolyBC');

	% Compute savings policy
	mAssetsPrimeStar = w * (mmu * (1 - vEpsilonGrid) + (1 - ttau) * vEpsilonGrid) + ...
		(1 + r) * assetsMin - (mConditionalExpectation .^ (-1 / ssigma));
	mAssetsPrimeBC = max(mAssetsPrimeStar,aaBar * ones(nEpsilon,1));
		
else

	% Compute weights
	[vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = computeLinearWeights(vAssetsGrid,assetsMin);
		
	% Linear interpolation
	mAssetsPrimeBC = mAssetsPrime(:,vIndicesBelow) .* repmat(vWeightBelow',nEpsilon,1) + ...
		mAssetsPrime(:,vIndicesAbove) .* repmat(vWeightAbove',nEpsilon,1);
	
end

%----------------------------------------------------------------
% Compute stationary distribution from these decision rules
%----------------------------------------------------------------

% Initialize iteration
err = 100; iteration = 1; 
options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','notify-detailed',...
	'MaxFunEvals',50000,'TolFun',1e-12,'GradObj','on','MaxIter',1000);
%{ For older versions of MATLAB:
options = optimset('LargeScale','off','Display','notify-detailed',...
	'MaxFunEvals',50000,'TolFun',1e-12,'GradObj','on','MaxIter',1000);
%}

% Iteration
while err > err2 && iteration <= tol2
	
	%%%
	% Update density away from borrowing constraint
	%%%

	% Compute parameters of the distribution by minimization
	mParameters = zeros(nEpsilon,nMeasure+1);
	for iEpsilon = 1 : nEpsilon
		objectiveFunction = @(vParametersTilde) parametersResidual(vParametersTilde,squeeze(aGridMoments(iEpsilon,:,:)));
		[vParameters,normalization] = fminunc(objectiveFunction,zeros(nMeasure,1),options);
		mParameters(iEpsilon,:) = [1 / normalization; vParameters];
	end
	
	% Compute new moments and centered moments grid
	mMomentsNew = zeros(nEpsilon,nMeasure);
	aGridMomentsNew = zeros(nEpsilon,nAssetsQuadrature,nMeasure);
	
	for iEpsilon = 1 : nEpsilon
		
		% Compute first moment (uncentered)
		mMomentsNew(iEpsilon,1) = 0;
		for iEpsilonTilde = 1 : nEpsilon
		
			mMomentsNew(iEpsilon,1) = mMomentsNew(iEpsilon,1) + (1 - mHat(iEpsilonTilde,1)) * vEpsilonInvariant(iEpsilonTilde) * mEpsilonTransition(...
				iEpsilonTilde,iEpsilon) * mParameters(iEpsilonTilde,1) * vQuadratureWeights' * (mAssetsPrimeQuadrature(iEpsilonTilde,:)' .* ...
				exp(squeeze(aGridMoments(iEpsilonTilde,:,:)) * mParameters(iEpsilonTilde,2:nMeasure+1)')) + mHat(iEpsilonTilde,1) * ...
				vEpsilonInvariant(iEpsilonTilde) * mEpsilonTransition(iEpsilonTilde,iEpsilon) * mAssetsPrimeBC(iEpsilonTilde,1);	
				
		end
		
		mMomentsNew(iEpsilon,1) = mMomentsNew(iEpsilon,1) / vEpsilonInvariant(iEpsilon);
		aGridMomentsNew(iEpsilon,:,1) = vAssetsGridQuadrature - mMomentsNew(iEpsilon,1);
		
		% Compute higher order moments (centered)
		for iMoment = 2 : nMeasure
		
			mMomentsNew(iEpsilon,iMoment) = 0;
			
			for iEpsilonTilde = 1 : nEpsilon
			
				mMomentsNew(iEpsilon,iMoment) = mMomentsNew(iEpsilon,iMoment) + (1 - mHat(iEpsilonTilde,1)) * vEpsilonInvariant(iEpsilonTilde) * mEpsilonTransition(...
					iEpsilonTilde,iEpsilon) * mParameters(iEpsilonTilde,1) * vQuadratureWeights' * (((mAssetsPrimeQuadrature(iEpsilonTilde,:)' - ...
					mMomentsNew(iEpsilon,1)) .^ iMoment) .* exp(squeeze(aGridMoments(iEpsilonTilde,:,:)) * ...
					mParameters(iEpsilonTilde,2:nMeasure+1)')) + mHat(iEpsilonTilde,1) * vEpsilonInvariant(iEpsilonTilde) * ...
					mEpsilonTransition(iEpsilonTilde,iEpsilon) * ((mAssetsPrimeBC(iEpsilonTilde,1) - mMomentsNew(iEpsilon,1)) .^ iMoment);
					
			end
			
			mMomentsNew(iEpsilon,iMoment) = mMomentsNew(iEpsilon,iMoment) / vEpsilonInvariant(iEpsilon);
			aGridMomentsNew(iEpsilon,:,iMoment) = (vAssetsGridQuadrature' - mMomentsNew(iEpsilon,1)) .^ iMoment - ...
				mMomentsNew(iEpsilon,iMoment);
				
		end
		
	end

	%%%
	% Update mass at borrowing constraint
	%%%
	
	mHatNew = zeros(nEpsilon,1);
	
	for iEpsilon = 1 : nEpsilon;
	
		for iEpsilonTilde = 1 : nEpsilon
		
			mHatNew(iEpsilon,1) = mHatNew(iEpsilon,1) + (1 - mHat(iEpsilonTilde,1)) * vEpsilonInvariant(iEpsilonTilde) * ...
				mEpsilonTransition(iEpsilonTilde,iEpsilon) * mParameters(iEpsilonTilde,1) * vQuadratureWeights' * ...
				((mAssetsPrimeQuadrature(iEpsilonTilde,:)' <= aaBar + 1e-8) .* exp(squeeze(aGridMoments(iEpsilonTilde,:,:)) * mParameters(iEpsilonTilde,2:nMeasure+1)')) + ...
				mHat(iEpsilonTilde,1) * vEpsilonInvariant(iEpsilonTilde) * mEpsilonTransition(iEpsilonTilde,iEpsilon) * ...
				(mAssetsPrimeBC(iEpsilonTilde,1) <= aaBar + 1e-8);
				
		end
			
		mHatNew(iEpsilon,1) = mHatNew(iEpsilon,1) / vEpsilonInvariant(iEpsilon);
	
   end
	
	%%%
	% Update iteration
	%%%
	
	err = max([max(abs(mMomentsNew(:) - mMoments(:))),max(abs(mHatNew(:) - mHat(:)))]);
	iteration = iteration + 1;
	mMoments = mMomentsNew;
	aGridMoments = aGridMomentsNew;
	mHat = mHatNew;
	
end

%----------------------------------------------------------------
% Return market clearing residual
%----------------------------------------------------------------

capitalNew = (vEpsilonInvariant .* (1 - mHat))' * mMoments(:,1) + aaBar * (vEpsilonInvariant .* mHat)' * ones(nEpsilon,1);
residual = capital - capitalNew;

% Also return optional outputs if requested
if nargout > 2

    mParametersOptional = mParameters;
    mMomentsOptional = mMoments;
	mHatOptional = mHat;
	
end
