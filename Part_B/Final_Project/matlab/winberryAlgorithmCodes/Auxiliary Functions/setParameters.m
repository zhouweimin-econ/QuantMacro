% Sets parameter values 
%
% Thomas Winberry, July 26th, 2016

%----------------------------------------------------------------
% Set economic parameters 
%----------------------------------------------------------------

global bbeta ssigma aaBar aalpha ddelta vEpsilonGrid mEpsilonTransition vEpsilonInvariant aggEmployment ...
	mmu ttau rrhoTFP ssigmaTFP
	
% Preferences
bbeta = .96;										% discount factor (annual calibration)
ssigma = 1;											% coefficient of relative risk aversion
aaBar = 0;											% borrowing constraint

% Technology
aalpha = .36;										% capital share
ddelta = .1;										% depreciation rate (annual calibration)

% Idioynscratic Shocks
vEpsilonGrid = [0;1];
aggEmployment = .93; uDuration = 1;
mEpsilonTransition = [uDuration / (1 + uDuration), 1 - (uDuration / (1 + uDuration));...
					 ((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration))),...
					 1 - ((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration)))];
vEpsilonInvariant = [1 - aggEmployment;aggEmployment];

% Unemployment benefits
mmu = .15;
ttau = mmu * (1 - aggEmployment) / aggEmployment;

% Aggregate Shocks
rrhoTFP = .859;										
ssigmaTFP = .014;

%----------------------------------------------------------------
% Set approximation parameters
%----------------------------------------------------------------

global nEpsilon nAssets nState assetsMin assetsMax nAssetsFine nStateFine nAssetsQuadrature nStateQuadrature ...
	nMeasure nMeasureCoefficients kRepSS maxIterations tolerance dampening splineOpt displayOpt
	
% Whether approximating decision rule with splines or polynomials
splineOpt = 0;	% if splineOpt = 1, use splines to approximate savings policy; if splineOpt = 0, use polynomials
				% to approximate conditional expectation function

% Whether to print out results from steady state computation
displayOpt = 'iter-detailed';       % 'iter-detailed' or 'off'

% Order of approximation
nEpsilon = 2;
nAssets = 25; % number of gridpoints in spline approximation or polynomials in polynomial approximation
nState = nEpsilon * nAssets;

% Bounds on grid space
kRepSS = ((aalpha * (aggEmployment ^ (1 - aalpha))) / ((1 / bbeta) - (1 - ddelta))) ^ (1 / (1 - aalpha));
assetsMin = aaBar;	assetsMax = 3 * kRepSS;

% Finer grid for analyzing policy functions
nAssetsFine = 100;
nStateFine = nEpsilon * nAssetsFine;

% Approximation of distribution
nMeasure = 3;
nAssetsQuadrature = 8;
nStateQuadrature = nEpsilon * nAssetsQuadrature;
nMeasureCoefficients = nEpsilon * nMeasure;

% Iteration on individual decisions
maxIterations = 2e4;
tolerance = 1e-5;
dampening = .95;