## A Summary

### Heterogenous Agent model with endogenous labor supply: under Incomplete Market, with idiosyncratic risk (employment and unemployment), with and without aggregate risk


**For a background study**:

- "Micro and Macro Labor Supply Elasticity" from  [NYU Macro II by Gianluca Violante](https://sites.google.com/a/nyu.edu/glviolante/teaching/macrotheory)

#### *Question 1* is the case without aggregate risk.

Andreas Muller provide a sample matlab code for [Aiyagari model with endogenous labor supply](https://sites.google.com/site/mrandreasmueller/resources?authuser=0).

#### *Question 2* is *Question 1* with aggregate risk (good time and bad time TFP shocks).

- [Incomplete Markets, Labor Supply and Capital Accumulation. Albert, Marcet. Francesc, Obiols-Homs and Philippe, Weil (2003)](http://digital.csic.es/bitstream/10261/57847/1/Incomplete%20Markets.pdf)
gives a complete Algorithm for computing such issues, and also:

- "Solving for the Equilibrium in Models withIdiosyncratic and Aggregate Risk" from [Lecture 12, NYU Quant Macro by Gianluca Violante](http://www.econ.nyu.edu/user/violante/NYUTeaching/QM/Fall15/Lectures/Lecture12_KS_Slides.pdf)

**Notice**: we have a given well-explained transition matrix based on Krusell and Smith (1998), so we don't need to
discretize the stochastic TFP transition process (i.e. Tauchen Method);
We are required to extend from PS1 solution, so I use discretization grid methods for this problem, which takes 33 mins without updating rules. 
Many ways to reduce the time, such as Perturbation methods, i.e.

Christian Bayer proposed a way for solving K&S within seconds in his new working paper (matlab code provided): [Solving heterogeneous agent models in discrete time with many idiosyncratic states by perturbation methods](http://www.cepr.org/active/publications/discussion_papers/dp.php?dpno=13071)
,[MATLAB code](http://www.wiwi.uni-bonn.de/bayer/REPLICATION/Linearize.zip)
