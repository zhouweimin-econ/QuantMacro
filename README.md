# Quantitative Macroeconomics, 2018-2019

## Economics PhD 2nd Year Course, UAB

### Intro
Our course is separated into two parts.

Part A is from basic numerical methods to
solution of Representative Agent Models (permanent income, life-cycle) with VFI, PFI;
and eventually going to Heterogenous Agents with Incomplete Markets.
In the end, we need to solve Recursive Stationary Equilibrium of the Aiyagari-Bewley-Hugget-Imrohoroglu (ABHI) Model.

Part B accesses to Krusell and Smith (1998) model (Heterogenous agents, incomplete markets and aggregate uncertainty), and defaultable sovereign debt.

### Background Readings:
- Jonathan Heathcote, Kjetil Storesletten, and Gianluca Violante (2009). [Quantitative Macroeconomics with Heterogeneous Households](http://www.annualreviews.org/doi/pdf/10.1146/annurev.economics.050708.142922), Annual Review of Economics
- Guvenen, Fatih (2012). [Macroeconomics with Heterogeneity: A Practical Guide](http://www.richmondfed.org/publications/research/economic_quarterly/2011/q3/pdf/guvenen.pdf), Richmond Fed QR
- Per Krusell, Anthony Smith (2006). [quantitative macroeconomic models with heterogeneous agents](http://aida.wss.yale.edu/smith/paper15.pdf), Advances in Economics and Econometrics: Theory and Applications, Ninth World Congress

### Textbook Reference:
- [Dynamic General Equilibrium Modeling, Computational Methods and Applications](https://www.wiwi.uni-augsburg.de/vwl/maussner/dge_buch/dge_book_2ed/downloads_2nd/): with relative codes for Fortran

- [Applied Computational Economics and Finance](http://www4.ncsu.edu/unity/users/p/pfackler/www/compecon/): with CompEcon Toolbox for Matlab

### Similar Course: 

[Quantitative Economics, NYU by Gianluca Violante](https://sites.google.com/a/nyu.edu/glviolante/teaching/quantmacro)

---

## Main Strategies of solving our models:

### 1. conventional Value Function Iteration (VFI) and Policy Function Iteration (PFI)

VFI and PFI is well-known, the basic algorithms of dynamic programming.
   
eg.1: See [Part_A/PS4-solution](https://github.com/zhouweimin233/QuantMacro/tree/master/Part_A/PS4-solution) where we practiced VFI and PFI using a couple of improvements to speed up (e.g., concavity of value function, previous solution, monotonicity of decision rules).

eg.2: See [Part_B/PS1-solution/main.m](https://github.com/zhouweimin233/QuantMacro/blob/master/Part_B/PS1-solution/main.m), which I use VFI associated with KS algorithm for solving Krusell and Smith model (~300 minutes MATLAB). 

For **Discretization**, we practiced evenly-spaced grids, chebyshev polynomials (finer grids yield better approximations but are costly in terms of computer time, [for algorithm containing numerical optimization and root-finding, Matlab performs worse than Python and Julia](https://web.stanford.edu/~maliars/Files/Files/CEPR-DP13210.pdf))

**Costs and benefits of VFI and PFI**: simple and have a solution for sure, but slow and inaccurate.

### Some extensions: 

- **Euler Equation iteration** 

EEI is faster than VFI (still pretty slow), but much more accurate. For EEI, we update policy function with Euler Equation. 

- **Endogenous Gridpoints Methods**

similar to EEI, but avoid solving non-linear equations. 

eg.1: [Christopher Carroll’s endogenous gridpoint method](http://econ.jhu.edu/people/ccarroll/EndogenousArchive.zip)


### 2. Perturbation method (for more details, go to [Part_B/Final Project](https://github.com/zhouweimin233/QuantMacro/edit/master/Part_B/Final_Project) )

Perturbation method (local) approximate the policy function around the no uncertainty case, and this can deal with large problems but require model equations be differentiable, hence, cannot deal perfectly with zero lower bound and occasionally binding constraints.

To approximate the policy function, we could apply 1st-order taylor expansion, or higher-order, and some more recent resources from Matlab toolbox, the "Automatic Differentiation" method. 

eg.1:  Perturbation methods of solving Krusell&Smith(1998): [Solving the Incomplete Markets Model withAggregate Uncertainty Using a Perturbation Method. Kim, Kim, and Kollmann](http://www.wouterdenhaan.com/suite/finalversion-KKK.pdf). 
 
eg.2 (can be <10 seconds MATLAB): with the borrowing constraint replaced by a penalty function and the finite-state Markov process replaced by a stochastic process with continuous support; or the most recent work done by Christian Bayer, ["Solving heterogeneous agent models in discrete time with many idiosyncratic states by perturbation methods", (with Ralph Luetticke), July 2018, CEPR DP No. 13071](http://www.wiwi.uni-bonn.de/bayer/Working_Papers.html), [MALTAB code](http://www.wiwi.uni-bonn.de/bayer/REPLICATION/Linearize.zip)

### 3. Projection methods  (for more details, go to [Part_B/Final Project](https://github.com/zhouweimin233/QuantMacro/edit/master/Part_B/Final_Project) )

eg.1: **parameterization of the cross-sectional distribution**: do not obtain next period's cross-sectional moments by simulation techniques,but by explicitly integrating the individual choices;

eg.2: **no parameterization of the cross-sectional distribution**: derives the aggregate laws of motion directly from the individual policy rules by simply aggregating them;



### 4. Linearization
We learned this method during our Macro II by using Dynare, billions of resources could be found via Internet. 

---

### Advanced Extensions
**Perturbation + Projection methods**:

I will replicate this method in my [Final Project of Unit 2](https://github.com/zhouweimin233/QuantMacro/tree/master/Part_B/Final_Project). 

A guide to understand [different methods on solving K&S model](http://www.econ.nyu.edu/user/violante/NYUTeaching/QM/Spring15/Lectures/Lecture14_KS_Slides.pdf) and comparisons of different solution methods to solve K&S model: [Comparison of solutions to the incomplete markets model with aggregate uncertainty. Wouter den Haan](http://www.wouterdenhaan.com/papers/comparison.pdf)

eg.1: [Reiter Michael, 2009. "Solving heterogeneous-agent models by projection and perturbation," Journal of Economic Dynamics and Control, Elsevier, vol. 33(3), pages 649-665, March.](https://ideas.repec.org/a/eee/dyncon/v33y2009i3p649-665.html)

eg.2: [Thomas Winberry, 2018 forthcoming. "A Method for Solving and Estimating Heterogeneous Agent Macro Models," Quantitative Economics](http://faculty.chicagobooth.edu/thomas.winberry/research/index.html)

---

**[Parameterized Expectations Algorithm (PEA)](https://web.stanford.edu/~maliars/Files/Codes.html)**

**[Explicit Aggregation (Xpa)](http://www.wouterdenhaan.com/numerical/methodsheteroxpa.pdf)**

**Euler Equation iteration with projection methods**

**...**
