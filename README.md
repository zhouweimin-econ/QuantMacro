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

## Solution Methods for Discrete-Time HA model with Idiosyncratic Shocks and Aggregate Uncertainty:

### 1. conventional VFI and PFI: 

eg.1: See [Part_A/PS4-solution](https://github.com/zhouweimin233/QuantMacro/tree/master/Part_A/PS4-solution) where we practiced VFI and PFI with different improvements.

eg.2: See [Part_B/PS1-solution/main.m](https://github.com/zhouweimin233/QuantMacro/blob/master/Part_B/PS1-solution/main.m), which I use KS algorithm(~300 minutes MATLAB). 

Notice that we use evenly-spaced grids (finer grids yield better approximations but are costly in terms of computer time, [for algorithm containing numerical optimization and root-finding, Matlab performs worse than Python and Julia](https://web.stanford.edu/~maliars/Files/Files/CEPR-DP13210.pdf))

Costs and benefits of VFI:  The benefit of simple VFI is that it is very simple, but it is not very accurate and it is slow.

### 2. Euler Equation iteration 

EEI is faster than VFI, but is still pretty slow, but much more accurate. For EEI, we update policy function with Euler Equation. 

### 3. Endogenous Gridpoints Methods

similar to EEI, but avoid solving non-linear equations. 

eg.1: [Christopher Carroll’s endogenous gridpoint method](http://econ.jhu.edu/people/ccarroll/EndogenousArchive.zip)

---

### 4. Improvements made on above methods: 


**4.1. Linearization**:
We learned this method during our Macro II by using Dynare, billions of resources could be found via Internet. 


**4.2 Perturbation methods**

eg.1:  Perturbation methods of solving Krusell&Smith(1998): [Solving the Incomplete Markets Model withAggregate Uncertainty Using a Perturbation Method. Kim, Kim, and Kollmann](http://www.wouterdenhaan.com/suite/finalversion-KKK.pdf). 

eg.2: see [Part_B/Final Project](https://github.com/zhouweimin233/QuantMacro/edit/master/Part_B/Final_Project).
 
eg.3 (can be <10 seconds MATLAB): with the borrowing constraint replaced by a penalty function and the finite-state Markov process replaced by a stochastic process with continuous support; or the most recent work done by Christian Bayer, ["Solving heterogeneous agent models in discrete time with many idiosyncratic states by perturbation methods", (with Ralph Luetticke), July 2018, CEPR DP No. 13071](http://www.wiwi.uni-bonn.de/bayer/Working_Papers.html), [MALTAB code](http://www.wiwi.uni-bonn.de/bayer/REPLICATION/Linearize.zip)
- See a comparison of different solution methods to solve Incomplete markets with Aggregate Risk (Krusell&Smith(1998) version):
[Comparison of solutions to the incomplete markets model with aggregate uncertainty. Wouter den Haan](http://www.wouterdenhaan.com/papers/comparison.pdf)

**4.3 Projection methods** 

eg.1: **parameterization of the cross-sectional distribution**: do not obtain next period's cross-sectional moments by simulation techniques,but by explicitly integrating the individual choices;

eg.2: **no parameterization of the cross-sectional distribution**: derives the aggregate laws of motion directly from the individual policy rules by simply aggregating them;

**4.4 Mixture of Perturbation and Projection methods**:

I will replicate this method in my Final Project of Unit 2. A guide to understand [different methods on solving K&S model](http://www.econ.nyu.edu/user/violante/NYUTeaching/QM/Spring15/Lectures/Lecture14_KS_Slides.pdf)

eg.1: [Reiter Michael, 2009. "Solving heterogeneous-agent models by projection and perturbation," Journal of Economic Dynamics and Control, Elsevier, vol. 33(3), pages 649-665, March.](https://ideas.repec.org/a/eee/dyncon/v33y2009i3p649-665.html)

eg.2: [Thomas Winberry, 2018 forthcoming. "A Method for Solving and Estimating Heterogeneous Agent Macro Models," Quantitative Economics](http://faculty.chicagobooth.edu/thomas.winberry/research/index.html)

---

**4.5. [Parameterized Expectations Algorithm (PEA)](https://web.stanford.edu/~maliars/Files/Codes.html)**

**4.6. [Explicit Aggregation (Xpa)](http://www.wouterdenhaan.com/numerical/methodsheteroxpa.pdf)**

**4.7. Euler Equation iteration with projection methods**:

1. Just like previous algorithm, we still approximate Law of motion with our Perceived law of motion and update the parameters for simulation. 
2. create state variables nodes: <a href="https://www.codecogs.com/eqnedit.php?latex=s_{i,t}=\{\epsilon_{i,t},&space;k_{i,t},s_{t}\}$&space;with&space;$s_t&space;=&space;\{&space;z_t,&space;K_t,&space;m_t&space;\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s_{i,t}=\{\epsilon_{i,t},&space;k_{i,t},s_{t}\}$&space;with&space;$s_t&space;=&space;\{&space;z_t,&space;K_t,&space;m_t&space;\}" title="s_{i,t}=\{\epsilon_{i,t}, k_{i,t},s_{t}\}$ with $s_t = \{ z_t, K_t, m_t \}" /></a> where <a href="https://www.codecogs.com/eqnedit.php?latex=m_{t&plus;1}&space;=&space;\Gamma&space;(z_{t&plus;1},z_t,&space;s_t;&space;\eta_{\Gamma})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?m_{t&plus;1}&space;=&space;\Gamma&space;(z_{t&plus;1},z_t,&space;s_t;&space;\eta_{\Gamma})" title="m_{t+1} = \Gamma (z_{t+1},z_t, s_t; \eta_{\Gamma})" /></a>


3. approximate <a href="https://www.codecogs.com/eqnedit.php?latex=k_{i,t&plus;1}&space;(\cdot)&space;=&space;P_n&space;(s_{i,t};\eta_{P_n})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k_{i,t&plus;1}&space;(\cdot)&space;=&space;P_n&space;(s_{i,t};\eta_{P_n})" title="k_{i,t+1} (\cdot) = P_n (s_{i,t};\eta_{P_n})" /></a>

4. tranform the FOC: 

<a href="https://www.codecogs.com/eqnedit.php?latex=c_t&space;^{-\gamma}&space;=&space;E&space;[\beta&space;(r(z',K')&plus;(1-\delta))c_{t&plus;1}^{-\gamma}]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?c_t&space;^{-\gamma}&space;=&space;E&space;[\beta&space;(r(z',K')&plus;(1-\delta))c_{t&plus;1}^{-\gamma}]" title="c_t ^{-\gamma} = E [\beta (r(z',K')+(1-\delta))c_{t+1}^{-\gamma}]" /></a>

replace consumption by BC, FOC changes to:

<a href="https://www.codecogs.com/eqnedit.php?latex=\{&space;\begin{align*}&space;r(z,K)&plus;(1-\delta)k&space;&plus;&space;(1-\tau(z))w(z,K)l\epsilon&space;\\&plus;&space;\mu&space;w(z,K)(1-\epsilon)-P_n(s;\eta_{P_n})&space;\end{align*}&space;\}^{-\gamma}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\{&space;\begin{align*}&space;r(z,K)&plus;(1-\delta)k&space;&plus;&space;(1-\tau(z))w(z,K)l\epsilon&space;\\&plus;&space;\mu&space;w(z,K)(1-\epsilon)-P_n(s;\eta_{P_n})&space;\end{align*}&space;\}^{-\gamma}" title="\{ \begin{align*} r(z,K)+(1-\delta)k + (1-\tau(z))w(z,K)l\epsilon \\+ \mu w(z,K)(1-\epsilon)-P_n(s;\eta_{P_n}) \end{align*} \}^{-\gamma}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex==E[\{&space;\begin{align*}&space;\beta(&space;r(z',K')&plus;(1-\delta))\times(r(z',K')&&plus;1-\delta)P_{n}(s;\eta_{P_n})&space;\\&space;&plus;(1-\tau(z'))w(z',K')l\epsilon'&space;&\\&plus;&space;\mu&space;w(z',K')(1-\epsilon')-P_n(s';\eta_{P_n})&space;\end{align*}&space;\}^{-\gamma}]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?=E[\{&space;\begin{align*}&space;\beta(&space;r(z',K')&plus;(1-\delta))\times(r(z',K')&&plus;1-\delta)P_{n}(s;\eta_{P_n})&space;\\&space;&plus;(1-\tau(z'))w(z',K')l\epsilon'&space;&\\&plus;&space;\mu&space;w(z',K')(1-\epsilon')-P_n(s';\eta_{P_n})&space;\end{align*}&space;\}^{-\gamma}]" title="=E[\{ \begin{align*} \beta( r(z',K')+(1-\delta))\times(r(z',K')&+1-\delta)P_{n}(s;\eta_{P_n}) \\ +(1-\tau(z'))w(z',K')l\epsilon' &\\+ \mu w(z',K')(1-\epsilon')-P_n(s';\eta_{P_n}) \end{align*} \}^{-\gamma}]" /></a>


where the expectation is transition matrix times each possible state would happen.

5. determine prices remains the same as before: 

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align*}&space;r(z,K)&space;&=&space;z&space;\alpha&space;(\frac{K}{L(z)})^{\alpha&space;-1&space;}&space;\\&space;w(z,K)&space;&=&space;z&space;(1-\alpha)&space;(\frac{K}{L(z)})^{\alpha}&space;\\&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{align*}&space;r(z,K)&space;&=&space;z&space;\alpha&space;(\frac{K}{L(z)})^{\alpha&space;-1&space;}&space;\\&space;w(z,K)&space;&=&space;z&space;(1-\alpha)&space;(\frac{K}{L(z)})^{\alpha}&space;\\&space;\end{align*}" title="\begin{align*} r(z,K) &= z \alpha (\frac{K}{L(z)})^{\alpha -1 } \\ w(z,K) &= z (1-\alpha) (\frac{K}{L(z)})^{\alpha} \\ \end{align*}" /></a>

6. by minimizing the error term (LHS-RHS) from Step.4: <a href="https://www.codecogs.com/eqnedit.php?latex=\min&space;\sum&space;u^2_k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min&space;\sum&space;u^2_k" title="\min \sum u^2_k" /></a>,we could find the best fitted parameters. 

### ...
