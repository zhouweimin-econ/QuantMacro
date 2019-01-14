% Creates grids to use in various approximations
% 
% Thomas Winberry, July 26th, 2016

%---------------------------------------------------------------
% Grids for approximating individual decisions
%---------------------------------------------------------------

if splineOpt == 0	% if using polynomials and approximating conditional expectation

	global vAssetsGridZeros vAssetsGrid mEpsilonGrid mAssetsGrid mEpsilonPrimeGrid

	% Zeros of chebyshev polynomial
	vAssetsGridZeros = -cos(((2 * (1:nAssets)-1)' * pi) / (2 * nAssets));

	% Scale up to state space
	vAssetsGrid = scaleUp(vAssetsGridZeros,assetsMin,assetsMax);
	
else	% if using splines and approximating asset accumulation

	global vAssetsGrid mEpsilonGrid mAssetsGrid mEpsilonPrimeGrid

	% Grid over assets
	vAssetsGrid = exp(linspace(log(assetsMin + .01),log(assetsMax + .01),nAssets)');
	vAssetsGrid = vAssetsGrid - .01;
	
end
	
% Make matrix versions of the grids
mEpsilonGrid = repmat(vEpsilonGrid,[1 nAssets]);
mAssetsGrid = repmat(vAssetsGrid',[nEpsilon 1]);
mEpsilonPrimeGrid = repmat(vEpsilonGrid,[1 nEpsilon*nAssets]);

%---------------------------------------------------------------
% Fine grid, for histogram and plotting functions
%---------------------------------------------------------------

global vAssetsGridFine vAssetsGridFineZeros mEpsilonGridFine mAssetsGridFine mEpsilonPrimeGridFine 

% Assets grid
vAssetsGridFine = linspace(assetsMin,assetsMax,nAssetsFine)';

% Scale down to [-1,1]
vAssetsGridFineZeros = scaleDown(vAssetsGridFine,assetsMin,assetsMax);

% Make matrix versions of grids
mEpsilonGridFine = repmat(vEpsilonGrid,[1 nAssetsFine]);
mAssetsGridFine = repmat(vAssetsGridFine',[nEpsilon 1]);
mEpsilonPrimeGridFine = repmat(vEpsilonGrid,[1 nStateFine]);

%---------------------------------------------------------------
% Quadrature grid, to integrate density (away from borrowing constraint)
%---------------------------------------------------------------

global vQuadratureWeights vAssetsGridQuadratureZeros vAssetsGridQuadrature mEpsilonGridQuadrature ...
	mAssetsGridQuadrature

% Compute grid in the interval [-1, 1]
[vAssetsGridQuadratureZeros,vQuadratureWeights] = computeGaussLegendreQuadrature(nAssetsQuadrature);

% Scale up grid
vAssetsGridQuadrature = scaleUp(vAssetsGridQuadratureZeros,assetsMin+1e-1,assetsMax);

% Make matrix versions of the grids
mEpsilonGridQuadrature = repmat(vEpsilonGrid,[1 nAssetsQuadrature]);
mAssetsGridQuadrature = repmat(vAssetsGridQuadrature', [nEpsilon 1]);
