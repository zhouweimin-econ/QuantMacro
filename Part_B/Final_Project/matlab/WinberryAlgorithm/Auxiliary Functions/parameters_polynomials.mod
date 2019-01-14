// Declare parameters and load in their values for firstOrderDynamics_polynomials.mod
//
// Thomas Winberry, July 26th, 2016

//----------------------------------------------------------------
// Preliminaries
//----------------------------------------------------------------

// Load in files containing parameters
economicParameters = load('economicParameters');
approximationParameters = load('approximationParameters');
grids = load('grids');
polynomials = load('polynomials');

// Define economic parameters
parameters bbeta ssigma aaBar aalpha ddelta aggEmployment
	mmu ttau rrhoTFP ssigmaTFP;
//Load in their values
@#define nEconomicParameters = 10
for iParam = 1 : @{nEconomicParameters}
	parameterName = deblank(M_.param_names(iParam,:));
	if isfield(economicParameters,parameterName)
		M_.params(iParam) = eval(['economicParameters.' parameterName]);
	end
end

// Epsilon transition matrix
@#for iEpsilon in 1 : 2
	@#for iEpsilonPrime in 1 : 2
		parameters epsilonTransition_@{iEpsilon}_@{iEpsilonPrime};
	@#endfor
@#endfor
for iEpsilon = 1 : 2
	for iEpsilonPrime = 1 : 2
		M_.params(@{nEconomicParameters} + 2 * (iEpsilon - 1) + iEpsilonPrime) = ...
			economicParameters.mEpsilonTransition(iEpsilon,iEpsilonPrime);
	end
end

// Mass of invariant distrbution of idiosyncratic shocks
parameters epsilonMass_1 epsilonMass_2;
epsilonMass_1 = 1 - aggEmployment;
epsilonMass_2 = aggEmployment;

// Define some of the approximation parameters 
parameters nEpsilon nAssets nState assetsMin assetsMax nAssetsFine nStateFine nAssetsQuadrature nStateQuadrature 
	nMeasure nMeasureCoefficients kRepSS maxIterations tolerance dampening;
// Load in their values
@#define nApproximationParameters = 15
for iParam = 1 : @{nApproximationParameters}
	parameterName = deblank(M_.param_names(@{nEconomicParameters} + 6 + iParam,:));
	if isfield(approximationParameters,parameterName)
		M_.params(@{nEconomicParameters} + 6 + iParam) = eval(['approximationParameters.' parameterName]);
	end
end

@#define nCounter = nEconomicParameters + 6 + nApproximationParameters

// Have to impose parameters used as counters directly
@#define nEpsilon = 2
@#define nAssets = 25
@#define nMeasure = 3
@#define nAssetsQuadrature = 8

@#define nState = nEpsilon * nAssets
@#define nStateQuadrature = nEpsilon * nAssetsQuadrature

//
//----------------------------------------------------------------
// Grids for approximating conditional expectation
//----------------------------------------------------------------

// 
// Employment
//

// Define the grids
parameters epsilonGrid_1 epsilonGrid_2;

// Assign values
epsilonGrid_1 = 0;
epsilonGrid_2 = 1;

@#define nCounter = nCounter + 2

//
// Assets
//

// Define the grids
@#for iAssets in 1 : nAssets
	parameters assetsGrid_@{iAssets};
@#endfor

// Assign values (must be in the same order that the parameters were declared)
for iAssets = 1 : @{nAssets}
	M_.params(@{nCounter} + iAssets) = grids.vAssetsGrid(iAssets);
end

// Update counter for future parameter assignments
@#define nCounter = nCounter + nAssets

//----------------------------------------------------------------
// Quadrature grid and weights
//----------------------------------------------------------------

// Define the parameters
@#for iAssets in 1 : nAssetsQuadrature
	parameters quadratureGrid_@{iAssets};
	parameters quadratureWeights_@{iAssets};
@#endfor

// Assign values
for iAssets = 1 : @{nAssetsQuadrature}
	M_.params(@{nCounter} + 2 * (iAssets - 1) + 1) = grids.vAssetsGridQuadrature(iAssets);
	M_.params(@{nCounter} + 2 * (iAssets - 1) + 2) = grids.vQuadratureWeights(iAssets);
end

@#define nCounter = nCounter + 2 * nAssetsQuadrature

//----------------------------------------------------------------
// Conditional expectation polynomials
//----------------------------------------------------------------

//
// Chebyshev polynomials
//

// Define the parameters
@#for iAssets in 1 : nAssets
	@#for iPower in 1 : nAssets
		parameters expectationPoly_@{iAssets}_@{iPower};
	@#endfor
@#endfor

// Assign values
for iAssets = 1 : @{nAssets}
	for iPower = 1 : @{nAssets}
		M_.params(@{nCounter} + @{nAssets} * (iAssets - 1) + iPower) = polynomials.vAssetsPoly(iAssets,iPower);
	end
end

@#define nCounter = nCounter + nAssets * nAssets

// 
// Squared terms in Chebyshev interpolation
//

// Define the parameters
@#for iAssets in 1 : nAssets
	parameters expectationPolySquared_@{iAssets};
@#endfor

// Assign the values
for iAssets = 1 : @{nAssets}
	M_.params(@{nCounter} + iAssets) = polynomials.vAssetsPolySquared(iAssets);	
end

@#define nCounter = nCounter + nAssets

//----------------------------------------------------------------
// Quadrature grid polynomials
//----------------------------------------------------------------

// Define the parameters
@#for iAssets in 1 : nAssetsQuadrature
	@#for iPower in 1 : nAssets
		parameters quadraturePoly_@{iAssets}_@{iPower};		
	@#endfor
@#endfor

// Assign values
for iAssets = 1 : @{nAssetsQuadrature}
	for iPower = 1 : @{nAssets}	
		M_.params(@{nCounter} + @{nAssets} * (iAssets - 1) + iPower) = polynomials.vAssetsPolyQuadrature(iAssets,iPower);		
	end	
end

@#define nCounter = nCounter + nAssetsQuadrature * nAssets

//----------------------------------------------------------------
// Borrowing constraint polynomials
//----------------------------------------------------------------

// Define the parameters
@#for iPower in 1 : nAssets
	parameters bcPoly_@{iPower};
@#endfor

// Assign values
for iPower = 1 : @{nAssets}
	M_.params(@{nCounter} + iPower) = polynomials.vAssetsPolyBC(iPower);
end