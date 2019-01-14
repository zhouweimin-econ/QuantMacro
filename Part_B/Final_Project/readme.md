### Final Project Summary: 

#### Using Reiter method on solving Krusell&Smith Model with Heterogeneous Agents.

Since the distribution of agents is an infinite dimension object (many state variables), one method is approximating the distribution using only few moments, i.e. Krusell&Smith simulation, etc. And others such as using a perturbation method are more popular recently, i.e. [LeGrand and Ragot](http://www.wouterdenhaan.com/teach/legrandragot.pdf) who truncates the history of idiosyncratic shock; Reiter who solves the individual problem with a projection method, and approximate the law of motion of aggregate variables and the distribution with a perturbation method.

Previously, in PS1, we tried the easy version of KS Algorithm. During the last few lectures, we mentioned Perturbation method, which is more used for Quantitative models, thus, I'm planning to learn this method in my final project. After learning the Literature, I revised my Repository readme.md, to demonstrate a more comprehensive structure. And I notice that Reiter method is more popular on solving KS model, which is a mixture of Projection and Perturbation, I found a few matlab codes from series of Quant Macro courses (I mentioned in my readme.md), and during learning Reiter methods, I download [supplementary matlab files](https://github.com/zhouweimin233/QuantMacro/tree/master/Part_B/Final_Project/matlab/Reiter/supplmentary) from many resources which are mainly for computing:

- "Automatic Differentiation", from MATLAB toolbox: ajac.m; amatinit.m 

- Build Transition Matrix ,from [Mario J. Miranda](miranda.4@osu.edu) and [Paul L. Fackler](paul_fackler@ncsu.edu): BuildTrans.m; lookupfor.m

- Step 2: linear_solution.m(solve for dynamics) require: gensys.m; qzdiv,m; qzswitch.m (find approximation coefficients Using GENSYS algorithm, from Christopher A. Sims (2001)); jacob (compute Jacobian, from CompEcon2016); 

- Convert VAR(1) into Markov-Chain using Tauchen's method, from Github: tauchen.m

- Computes impulse responses for solution of linear RE models, from CompEcon2016: impresp_sims.m