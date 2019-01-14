% Shell which declares parameters and calls Dynare to solve model using
% first order approximation of aggregate dynamics
%
% Thomas Winberry, July 26th, 2016

clear all
close all
clc

cd('./Auxiliary Functions');

%----------------------------------------------------------------
% Set parameters
%----------------------------------------------------------------
setParameters;

%----------------------------------------------------------------
% Compute approximation tools
%----------------------------------------------------------------

% Grids
computeGrids;

% Polynomials over grids (only if using polynomials to approximate conditional expectation)
if splineOpt == 0
	computePolynomials;
end

%----------------------------------------------------------------
% Save parameters in .mat files to import into Dynare 
%----------------------------------------------------------------

if splineOpt == 0	% if using polynomials to approximate individual decisions

	% Economic parameters
	save economicParameters.mat bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment ...
		mmu ttau rrhoTFP ssigmaTFP mEpsilonTransition
		
	% Approximation parameters
	save approximationParameters.mat nEpsilon nAssets nState assetsMin assetsMax nAssetsFine nStateFine nAssetsQuadrature nStateQuadrature ...
		nMeasure nMeasureCoefficients kRepSS maxIterations tolerance dampening
		
	% Grids
	save grids.mat vAssetsGridZeros vAssetsGrid mEpsilonGrid mAssetsGrid mEpsilonPrimeGrid vAssetsGridFine ...
		vAssetsGridFineZeros mEpsilonGridFine mAssetsGridFine mEpsilonPrimeGridFine vQuadratureWeights ...
		vAssetsGridQuadratureZeros vAssetsGridQuadrature mEpsilonGridQuadrature mAssetsGridQuadrature
		
	% Polynomials
	save polynomials.mat vAssetsPoly vAssetsPolySquared vAssetsPolyFine vAssetsPolyQuadrature vAssetsPolyBC
	
else	% if using splines to approximate individual decisions

	% Economic parameters
	save economicParameters.mat bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment ...
		mmu ttau rrhoTFP ssigmaTFP mEpsilonTransition
		
	% Approximation parameters
	save approximationParameters.mat nEpsilon nAssets nState assetsMin assetsMax nAssetsFine nStateFine nAssetsQuadrature nStateQuadrature ...
		nMeasure nMeasureCoefficients kRepSS maxIterations tolerance dampening
		
	% Grids
	save grids.mat vAssetsGrid mEpsilonGrid mAssetsGrid mEpsilonPrimeGrid vAssetsGridFine ...
		mEpsilonGridFine mAssetsGridFine mEpsilonPrimeGridFine vQuadratureWeights ...
		vAssetsGridQuadrature mEpsilonGridQuadrature mAssetsGridQuadrature
	
end

%----------------------------------------------------------------
% Run Dynare
%----------------------------------------------------------------

if splineOpt == 0	% if using polynomials to approximate individual decisions

	dynare firstOrderDynamics_polynomials
	
else	% if using splines to approximate individual decisions

	dynare firstOrderDynamics_splines
	
end

cd('../')