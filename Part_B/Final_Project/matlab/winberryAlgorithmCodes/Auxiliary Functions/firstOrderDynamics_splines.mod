// Dynare shell which declares model and solves for aggregate dynamics using
// first order approximation (when approximating savings decisions with splines)
//
// Thomas Winberry, July 26th, 2016

//----------------------------------------------------------------
// Load parameters
//----------------------------------------------------------------

@#include "parameters_splines.mod"

//----------------------------------------------------------------
// Define variables
//----------------------------------------------------------------

@#include "variables_splines.mod"

//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;

@#include "equations_splines.mod"

end;

//----------------------------------------------------------------
// 4. Computation
//----------------------------------------------------------------

// Specify shock process

shocks;
    var aggregateTFPShock = 1;
end;

options.steadystate.nocheck = 1;

// Compute steady state (nocheck option ensures that Dynare runs even if steady
// state only computed approximately, i.e., with small numerical error)
//steady(nocheck);

// Check regularity conditions (turn on to check)
//check;
//model_diagnostics;
//model_info;

// Simulate
stoch_simul(order=1,hp_filter=100,irf=40) aggregateTFP logAggregateOutput 
	logAggregateConsumption logAggregateInvestment logWage r;
