// Specifies equations for firstOrderDynamics.mod
//
// Thomas Winberry and Alp Tuncay, July 26th, 2016

//----------------------------------------------------------------
// Individual decisions (#equations = nEpsilon * nAssets)
//----------------------------------------------------------------

@#for iEpsilon in 1 : nEpsilon

	@#for iAssets in 1 : nAssets
	
		// Compute basis representation of spline approximation
		# basis_1_@{iEpsilon}_@{iAssets} = ((assetsGrid_2 - coefficient_@{iEpsilon}_@{iAssets}) / 
			(assetsGrid_2 - assetsGrid_1)) * (coefficient_@{iEpsilon}_@{iAssets} < assetsGrid_2);
			
		@#for iAssetsTilde in 2 : nAssets - 1
			# basis_@{iAssetsTilde}_@{iEpsilon}_@{iAssets} = ((coefficient_@{iEpsilon}_@{iAssets} - 
				assetsGrid_@{iAssetsTilde-1}) / (assetsGrid_@{iAssetsTilde} - assetsGrid_@{iAssetsTilde-1})) * 
				(assetsGrid_@{iAssetsTilde-1} <= coefficient_@{iEpsilon}_@{iAssets}) * (assetsGrid_@{iAssetsTilde} > 
				coefficient_@{iEpsilon}_@{iAssets}) + ((assetsGrid_@{iAssetsTilde+1} - coefficient_@{iEpsilon}_@{iAssets}) / 
				(assetsGrid_@{iAssetsTilde+1} - assetsGrid_@{iAssetsTilde})) * (assetsGrid_@{iAssetsTilde} <=
				coefficient_@{iEpsilon}_@{iAssets}) * (assetsGrid_@{iAssetsTilde+1} > coefficient_@{iEpsilon}_@{iAssets});
		@#endfor
		
		# basis_@{nAssets}_@{iEpsilon}_@{iAssets} = ((coefficient_@{iEpsilon}_@{iAssets} - assetsGrid_@{nAssets-1}) / 
			(assetsGrid_@{nAssets} - assetsGrid_@{nAssets-1})) * (coefficient_@{iEpsilon}_@{iAssets} >= assetsGrid_@{nAssets-1});
			
		// Interpolate asset choices over assetPrime
		@#for iEpsilonPrime in 1 : nEpsilon
		
			# assetsPrimePrime_@{iEpsilonPrime}_@{iEpsilon}_@{iAssets} = (0
			@#for iCoefficient in 1 : nAssets
				+ coefficient_@{iEpsilonPrime}_@{iCoefficient}(+1) * basis_@{iCoefficient}_@{iEpsilon}_@{iAssets}
			@#endfor
			);
			
			# consumptionPrime_@{iEpsilonPrime}_@{iEpsilon}_@{iAssets} = w(+1) * (mmu * (1 - epsilonGrid_@{iEpsilonPrime}) + 
				(1 - ttau) * epsilonGrid_@{iEpsilonPrime}) + (1 + r(+1)) * coefficient_@{iEpsilon}_@{iAssets} - 
				assetsPrimePrime_@{iEpsilonPrime}_@{iEpsilon}_@{iAssets};
				
		@#endfor
			
		# consumption_@{iEpsilon}_@{iAssets} = (bbeta * (1 + r(+1)) * (0
		@#for iEpsilonPrime in 1 : nEpsilon
			+ epsilonTransition_@{iEpsilon}_@{iEpsilonPrime} * (consumptionPrime_@{iEpsilonPrime}_@{iEpsilon}_@{iAssets} ^ 
				(-ssigma))
		@#endfor
		)) ^ (-1 / ssigma);
		
		coefficient_@{iEpsilon}_@{iAssets} = max(aaBar,w * (mmu * (1 - epsilonGrid_@{iEpsilon}) + 
			(1 - ttau) * epsilonGrid_@{iEpsilon}) + (1 + r) * assetsGrid_@{iAssets} - consumption_@{iEpsilon}_@{iAssets});
			
	@#endfor
	
@#endfor

//----------------------------------------------------------------
// Compute various objects over quadrature grid for future use
//----------------------------------------------------------------

// Interpolate savings choice over quadrature grid

@#for iAssets in 1 : nAssetsQuadrature

	// Compute basis representation of spline representation
	# basisQuadrature_1_@{iAssets} = ((assetsGrid_2 - quadratureGrid_@{iAssets}) / 
		(assetsGrid_2 - assetsGrid_1)) * (quadratureGrid_@{iAssets} < assetsGrid_2);
		
	@#for iAssetsTilde in 2 : nAssets - 1
		# basisQuadrature_@{iAssetsTilde}_@{iAssets} = ((quadratureGrid_@{iAssets} - 
				assetsGrid_@{iAssetsTilde-1}) / (assetsGrid_@{iAssetsTilde} - assetsGrid_@{iAssetsTilde-1})) * 
				(assetsGrid_@{iAssetsTilde-1} <= quadratureGrid_@{iAssets}) * (assetsGrid_@{iAssetsTilde} > 
				quadratureGrid_@{iAssets}) + ((assetsGrid_@{iAssetsTilde+1} - quadratureGrid_@{iAssets}) / 
				(assetsGrid_@{iAssetsTilde+1} - assetsGrid_@{iAssetsTilde})) * (assetsGrid_@{iAssetsTilde} <=
				quadratureGrid_@{iAssets}) * (assetsGrid_@{iAssetsTilde+1} > quadratureGrid_@{iAssets});
	@#endfor
		
	# basisQuadrature_@{nAssets}_@{iAssets} = ((quadratureGrid_@{iAssets} - assetsGrid_@{nAssets-1}) / 
		(assetsGrid_@{nAssets} - assetsGrid_@{nAssets-1})) * (quadratureGrid_@{iAssets} >= assetsGrid_@{nAssets-1});	
		
	// Interpolate asset choices over this grid
	@#for iEpsilon in 1 : nEpsilon
		
		# assetsPrimeQuadrature_@{iEpsilon}_@{iAssets} = (0
		@#for iCoefficient in 1 : nAssets
			+ coefficient_@{iEpsilon}_@{iCoefficient} * basisQuadrature_@{iCoefficient}_@{iAssets}
		@#endfor
		);
	
	@#endfor

@#endfor

// Compute density away from borrowing constraint

@#for iEpsilon in 1 : nEpsilon

	@#for iAssets in 1 : nAssetsQuadrature
	
		# measurePDF_@{iEpsilon}_@{iAssets} = exp(0 + measureCoefficient_@{iEpsilon}_1 * (quadratureGrid_@{iAssets} - 
			moment_@{iEpsilon}_1(-1))
			@#for iMoment in 2 : nMeasure
				+ measureCoefficient_@{iEpsilon}_@{iMoment} * ((quadratureGrid_@{iAssets} - moment_@{iEpsilon}_1(-1)) ^ @{iMoment} - 
				moment_@{iEpsilon}_@{iMoment}(-1))
			@#endfor
			);
			
	@#endfor
	
	// Total mass of distribution
	# totalMass_@{iEpsilon} = 0
	@#for iAssets in 1 : nAssetsQuadrature
		+ quadratureWeights_@{iAssets} * measurePDF_@{iEpsilon}_@{iAssets}
	@#endfor
	;
	
@#endfor

//----------------------------------------------------------------
// Compute various objects at borrowing constraint for integrating distribution
//----------------------------------------------------------------

@#for iEpsilon in 1 : nEpsilon

	# assetsPrimeBC_@{iEpsilon} = coefficient_@{iEpsilon}_1;
	
@#endfor

//----------------------------------------------------------------
// Relationship between moments of distribution and parameters
// (#equations = nEpsilon * nMeasure)
//----------------------------------------------------------------

@#for iEpsilon in 1 : nEpsilon
	
	// First moments (uncentered)
	moment_@{iEpsilon}_1(-1) = (0
		@#for iAssets in 1 : nAssetsQuadrature
			+ quadratureWeights_@{iAssets} * quadratureGrid_@{iAssets} * measurePDF_@{iEpsilon}_@{iAssets}
		@#endfor
		) / totalMass_@{iEpsilon};
	
	// Higher order moments (centered)
	@#for iMoment in 2 : nMeasure
	moment_@{iEpsilon}_@{iMoment}(-1) = (0
		@#for iAssets in 1 : nAssetsQuadrature
			+ quadratureWeights_@{iAssets} * measurePDF_@{iEpsilon}_@{iAssets} * ((quadratureGrid_@{iAssets} - 
				moment_@{iEpsilon}_1(-1)) ^ @{iMoment})
		@#endfor
		) / totalMass_@{iEpsilon};
		
	@#endfor
		
@#endfor

//----------------------------------------------------------------
// Law of motion for density away from borrowing constraint 
// (#equations = nEpsilon * nMeasure)
//----------------------------------------------------------------

@#for iEpsilon in 1 : nEpsilon
	
	// First moment (uncentered)
	moment_@{iEpsilon}_1 = (0
	@#for iEpsilonTilde in 1 : nEpsilon
		+ ((1 - mHat_@{iEpsilonTilde}(-1)) * epsilonMass_@{iEpsilonTilde} * epsilonTransition_@{iEpsilonTilde}_@{iEpsilon} * (0
		@#for iAssets in 1 : nAssetsQuadrature
			+ quadratureWeights_@{iAssets} * measurePDF_@{iEpsilonTilde}_@{iAssets} *
				assetsPrimeQuadrature_@{iEpsilonTilde}_@{iAssets}
		@#endfor
		) / totalMass_@{iEpsilonTilde}) + mHat_@{iEpsilonTilde}(-1) * epsilonMass_@{iEpsilonTilde} * epsilonTransition_@{iEpsilonTilde}_@{iEpsilon} 
			* assetsPrimeBC_@{iEpsilonTilde}
	@#endfor
	) / epsilonMass_@{iEpsilon};
	
	// Higher order moments (uncentered)
	@#for iMoment in 2 : nMeasure
		moment_@{iEpsilon}_@{iMoment} = (0
		@#for iEpsilonTilde in 1 : nEpsilon
			+ ((1 - mHat_@{iEpsilonTilde}(-1)) * epsilonMass_@{iEpsilonTilde} * epsilonTransition_@{iEpsilonTilde}_@{iEpsilon} * (0
			@#for iAssets in 1 : nAssetsQuadrature
				+ quadratureWeights_@{iAssets} * measurePDF_@{iEpsilonTilde}_@{iAssets} * 
					(assetsPrimeQuadrature_@{iEpsilonTilde}_@{iAssets} - moment_@{iEpsilon}_1) ^ @{iMoment}
			@#endfor
			) / totalMass_@{iEpsilonTilde}) + mHat_@{iEpsilonTilde}(-1) * epsilonMass_@{iEpsilonTilde} * epsilonTransition_@{iEpsilonTilde}_@{iEpsilon}
				* (assetsPrimeBC_@{iEpsilonTilde} - moment_@{iEpsilon}_1) ^ @{iMoment}
		@#endfor
		) / epsilonMass_@{iEpsilon};
	@#endfor

@#endfor

//----------------------------------------------------------------
// Law of motion for mass at borrowing constraint
// (#equations = nEpsilon)
//----------------------------------------------------------------

@#for iEpsilon in 1 : nEpsilon

	mHat_@{iEpsilon} = (0
	@#for iEpsilonTilde in 1 : nEpsilon
		+ ((1 - mHat_@{iEpsilonTilde}(-1)) * epsilonMass_@{iEpsilonTilde} * epsilonTransition_@{iEpsilonTilde}_@{iEpsilon} * (0
		@#for iAssets in 1 : nAssetsQuadrature
			+ quadratureWeights_@{iAssets} * measurePDF_@{iEpsilonTilde}_@{iAssets} *
				(assetsPrimeQuadrature_@{iEpsilonTilde}_@{iAssets} <= aaBar + 1e-8)
		@#endfor
		) / totalMass_@{iEpsilonTilde}) + mHat_@{iEpsilonTilde}(-1) * epsilonMass_@{iEpsilonTilde} * epsilonTransition_@{iEpsilonTilde}_@{iEpsilon} 
			* (assetsPrimeBC_@{iEpsilonTilde} <= aaBar + 1e-8)
	@#endfor
	) / epsilonMass_@{iEpsilon};

@#endfor

//----------------------------------------------------------------
// Factor prices (# equations = 2)
//----------------------------------------------------------------

# aggregateCapital = (1 - aggEmployment) * moment_1_1(-1) + aggEmployment * moment_2_1(-1);

r = exp(aggregateTFP) * aalpha * (aggregateCapital ^ (aalpha - 1)) * (aggEmployment ^ (1 - aalpha)) - ddelta;
w = exp(aggregateTFP) * (aggregateCapital ^ aalpha) * (1 - aalpha) * (aggEmployment ^ (-aalpha));

//----------------------------------------------------------------
// Law of motion for aggregate TFP (# equations = 1)
//----------------------------------------------------------------

aggregateTFP = rrhoTFP * aggregateTFP(-1) + ssigmaTFP * aggregateTFPShock;

//----------------------------------------------------------------
// Auxiliary variables of interest (# equations = 4)
//----------------------------------------------------------------

// Output
logAggregateOutput = log(exp(aggregateTFP) * (aggregateCapital ^ aalpha) * (aggEmployment ^ (1 - aalpha)));

// Investment
logAggregateInvestment = log((1 - aggEmployment) * moment_1_1 + aggEmployment * moment_2_1 - (1 - ddelta) * aggregateCapital);

// Consumption
logAggregateConsumption = log(exp(logAggregateOutput) - exp(logAggregateInvestment));

// Wage
logWage = log(w);
