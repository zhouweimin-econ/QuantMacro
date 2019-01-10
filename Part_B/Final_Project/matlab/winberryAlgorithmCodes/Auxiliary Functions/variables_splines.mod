// Declare variables for firstOrderDynamics.mod
//
// Thomas Winberry, July 26th, 2016

//----------------------------------------------------------------
// Coefficients on linear spline approximating savings policy
//----------------------------------------------------------------

@#for iCoefficient in 1 : nAssets
	var coefficient_1_@{iCoefficient} coefficient_2_@{iCoefficient};
@#endfor

//----------------------------------------------------------------
// Density of households away from borrowing constraint
//----------------------------------------------------------------

// Moments of the distribution
@#for iMoment in 1 : nMeasure
    var moment_1_@{iMoment} moment_2_@{iMoment};
@#endfor

// Parameters of the distribution
@#for iParameter in 1 : nMeasure
    var measureCoefficient_1_@{iParameter} measureCoefficient_2_@{iParameter};  
@#endfor

//----------------------------------------------------------------
// Mass at borrowing constraint
//----------------------------------------------------------------

var mHat_1 mHat_2;

//----------------------------------------------------------------
// Prices
//----------------------------------------------------------------

var r w;

//----------------------------------------------------------------
// Aggregate TFP
//----------------------------------------------------------------

var aggregateTFP;

//----------------------------------------------------------------
// Auxiliary variables we're interested in
//----------------------------------------------------------------

var logAggregateOutput logAggregateInvestment logAggregateConsumption logWage;

//----------------------------------------------------------------
// Shocks
//----------------------------------------------------------------

varexo aggregateTFPShock;