### A Summary

For a complete illustration of the figure interpretation, one could refer to [Dynamic General Equilibrium Modelling, Computational Methods and Applications(2nd)](https://www.wiwi.uni-augsburg.de/vwl/maussner/dge_buch/dge_book_2ed/downloads_2nd/), by Burkhard Heer and Alfred Maussner, in section 8.2 <Transition Dynamics of Krusell and Smith (1998) model>, *page 397-404*. 
 
Besides, by [Comparison of solutions to the incomplete markets model with aggregate uncertainty, Wouter J. DEN HAAN (2009).](http://www.wouterdenhaan.com/datasuite.htm). Obviously, there's many methods to solve **Heterogenous Agents on Incomplete Markets with Idiosyncratic shocks and Aggregate Uncertainty**:
 - **KS algorithm**(~300 minutes MATLAB): also called KS simulation, the method we used for the Problem-Set;
 - **Projection methods I: parameterization of the cross-sectional distribution**: do not obtain next perio's cross-sectional moments by simulation techniques,but by explicitly integrating the individual choices;
 - **Projection methods II: no parameterization of the cross-sectional distribution**: derives the aggregate laws of motion directly from the individual policy rules by simply aggregating them;
 - **Perturbation approach**(can be <10 seconds MATLAB): i.e. with the borrowing constraint replaced by a penalty function and the finite-state Markov process replaced by a stochastic process with continuous support; or the most recent work done by Christian Bayer, ["Solving heterogeneous agent models in discrete time with many idiosyncratic states by perturbation methods", (with Ralph Luetticke), July 2018, CEPR DP No. 13071](http://www.wiwi.uni-bonn.de/bayer/Working_Papers.html), [MALTAB code](http://www.wiwi.uni-bonn.de/bayer/REPLICATION/Linearize.zip)
