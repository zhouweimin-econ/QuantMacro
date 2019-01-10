% Compute polynomials for approximations (only relevant if approximating conditional expectation)
% 
% Thomas Winberry, July 26th, 2016

%---------------------------------------------------------------
% Grid for approximating conditional expectation
%---------------------------------------------------------------

global vAssetsPoly vAssetsPolySquared

% Create polynomials
vAssetsPoly = computeChebyshev(nAssets,vAssetsGridZeros);

% Create squared terms for interpolation
vAssetsPolySquared = sum(vAssetsPoly .^ 2)';

%---------------------------------------------------------------
% Fine grid, for histogram and plots
%---------------------------------------------------------------

global vAssetsPolyFine

% Create polynomials
vAssetsPolyFine = computeChebyshev(nAssets,vAssetsGridFineZeros);

%---------------------------------------------------------------
% Quadrature grid, for integrating distribution
%---------------------------------------------------------------

global vAssetsPolyQuadrature

% Create polynomials
vAssetsPolyQuadrature = computeChebyshev(nAssets,vAssetsGridQuadratureZeros);

%---------------------------------------------------------------
% At borrowing constraint, for integrating mass point
%---------------------------------------------------------------

global vAssetsPolyBC

% Create polynomials
vAssetsPolyBC = computeChebyshev(nAssets,-1);	% at edge of state space