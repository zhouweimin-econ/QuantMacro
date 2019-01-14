This archive contains MATLAB and DYNARE software for solving the Krusell and Smith (1998) model, using
the method developed in "A Toolbox for Solving and Estimating Heterogeneous Agent Macro Models" by Thomas
Winberry.  It also contains the PDF file "dynareUserGruide.pdf" explaining how the codes work, and how to 
apply them to other heterogeneous agent models.

Thomas Winberry, July 27th, 2016.  Please email me at Thomas.Winberry@chicagobooth.edu with any feedback,
questions, or bug reports.  If you use these codes, please cite: Winberry (2016), "A Toolbox for Solving 
and Estimating Heterogeneous Agent Macro Models."


The following items are included:

--------------------------------------------------------------------------------------------------------
1.  MATLAB software for computing and analyzing the stationary equilibrium of the model with no aggregate
	shocks 
	
	"steadyState.m"

	This script will solve for the steady state of the model and plot decision rules and the stationary 
	distribution.  It has two options for computing individual decisions: (1) approximate the savings 
	policy rule a'(epsilon,a) using linear splines or (2) approximate the conditional expectation of
	next period's consumption using Chebyshev polynomials (see user guide for details).  The polynomial
	approximation is faster in this particular model, but I include them both as a template for other
	models as well (in some cases, splines may be necessary due to strong kinks).

	The script calls subroutines:

	a. "set Parameters.m": sets parameters.  This includes both parameters of the economic model, and
		parameters controlling the approximation of the model.  Importantly, set splineOpt = 1 to
		approximate the savings decision using linear splines or set splineOpt = 0 to approximate 
		the conditional expectation using polynomials.
	b. "coreSteadyState.m": solves for the market-clearing aggregate capital stock by solving the 
		market clearing condition.
	c. "scaleUp.m": transforms [-1,1] to [a,b]; used for computing Chebyshev polynomials
	d. "scaleDown.m": inverse of scaleUp.m
	e. "computeChebyshev.m": compute Chebyshev polynomials over a pre-defined grid
	f. "computeGaussLegendreQuadrature": compute weights and nodes for Gauss-Legendre 
		quadrature; used to integrate the parametric family
	g. "computeGrids.m": script which computes various grids used in approximation
	h. "computePolynomials.m": script which computes various polynomials used in approximation (only used
		if splineOpt = 0)
	i. "updateCoefficients_polynomials.m": updates polynomial coefficients of conditional expectation function, 
		which is used to approximate individual decisions (only used if splineOpt = 0)
	j. "updateCoefficients_splines.m": updates coefficients of linear splines of savings function, which
		is used to approximate individual decisions (only used if splineOpt = 1)
	k. "computeMCResidualHistogram.m": computes the residual of the market clearing condition that capital 
		supply (from integrating household savings decisions) equals capital demand (from the representative
		firm's first order condition).  Approximates the distribution using a nonparametric histogram and
		serves as initial guess for later approximating the distribution using the parametric family of
		exponential polynomials.
	l. "computeMCResidualPolynomials.m": does the same as computeMCResidualHistogram.m, but approximates the
		distribution using the parametric family of exponential polynomials.
	m. "parametersResidual.m": evaluates the objective function in the minimization problem defining the 
		parameters of the distribution, given moments.  Used in computeMCResidualPolynomials.m.

		
--------------------------------------------------------------------------------------------------------
2.	MATLAB and DYNARE software for computing and analyzing the dynamics of the model with aggregate
	shocks
	
	"dynamics.m"

	This script solves for the first order approximation of the model's dynamics.  You should be able to run only
	 this script and everything else is done automatically.  As with the steady state
	codes, it has two options for computing individual decisions: (1) approximate the savings 
	policy rule a'(epsilon,a) using linear splines or (2) approximate the conditional expectation of
	next period's consumption using Chebyshev polynomials.  

	The files are meant to provide a template for how to solve other models using the method in Dynare as well.
	They are liberally commented, and I also provide some additional information in the descriptions of the sub-
	routines below.  

	In addition to the subroutines described above, the script calls:
	
	a. "firstOrderDynamics_polynomials.mod" or "firstOrderDynamics_splines.mod": Dynare .mod shell which calls other
		files and computes the solution (_polynomials is called if splineOpt = 0 and _splines is called if splineOpt = 1).
		@# signs indicate a macro-processor command (see Dynare documentation).  NB: Dynare will not work unless
		steady(nocheck) is specified; otherwise, the steady state equilibrium conditions are only approximately satisfied, 
		which causes Dynare to shut down.  It is good practice in debugging to first check that these residuals are
		approximately zero.
	b. "parameters_polynomials.mod" or "parameters_splines.mod": declares and loads in parameter values.  NB: this includes
		the grids and if, applicable, polynomials for approximation.  The method treats these as parameters in Dynare because
		they are not endogeneous variables in the model.  "solveDynamicsFirstOrder.m" solves the grids and polynomials in 
		a .mat structure which is read into Dynare in this file.
		!!!NB: some of the approximation parameters (e.g., the size of the grids) must also be set here, in addition to in
		set_parameters.m!!!!
	c. "variables_polynomials.mod" or "variables_splines.mod": declares the variables in the model.  This includes the coefficients
		on individual decisions (either polynomials or splines, depending on splineOpt), the parameters and moments of the
		distribution, the mass of households at the borrowing constraint, prices, aggregate TFP, and aggregate output, 
		investment, and consumption.  Dynare solves for the dynamics of all of these variables.
	d. "equations_polynomials.mod" or "equations_splines.mod:" defines the equilibrium conditions of the model.  NB: compared to 
		a typical Dynare file for a representative agent model, this file libearlly uses the macro-processor (macro-processor
		commands are preceeded by @#; see the Dynare documentation for details).  This allows the user to write loops in Dynare,
		which are useful in characterizing individual decisions and the distribution because one can loop over the gridpoints
		in the individual state space.
	e. "firstOrderDynamics_polynomials_steadystate.m" or "firstOrderDynamics_splines_steadystate.m": an external Matlab .m file
		which computes the model's steady state.  It performs essentially the same commands as steadyState.m, but is
		in the form required to be called by Dynare.
