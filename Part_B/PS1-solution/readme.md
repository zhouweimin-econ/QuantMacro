### A Summary

Both PS1 and PS2 are solving **Heterogenous Agents on Incomplete Markets with Idiosyncratic shocks and Aggregate Uncertainty** based on  Krusell&Smith (1998) version

Comment: for a complete illustration of the figure interpretation, one could refer to [Dynamic General Equilibrium Modelling, Computational Methods and Applications(2nd)](https://www.wiwi.uni-augsburg.de/vwl/maussner/dge_buch/dge_book_2ed/downloads_2nd/), by Burkhard Heer and Alfred Maussner, in section 8.2 <Transition Dynamics of Krusell and Smith (1998) model>, *page 397-404*. 

---

**Basic Set-up**

- employment shocks (<a href="https://www.codecogs.com/eqnedit.php?latex=$\epsilon_{i,t}\in&space;\{0,1\}$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\epsilon_{i,t}\in&space;\{0,1\}$" title="$\epsilon_{i,t}\in \{0,1\}$" /></a>)

- Incomplete markets, with borrowing constraint <a href="https://www.codecogs.com/eqnedit.php?latex=k_{i,t&plus;1}\leq&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k_{i,t&plus;1}\leq&space;0" title="k_{i,t+1}\leq 0" /></a>

- transition probabilities calibrated by K&S(1998), where <a href="https://www.codecogs.com/eqnedit.php?latex=z_t$,&space;$\epsilon_{i,t}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?z_t$,&space;$\epsilon_{i,t}" title="z_t$, $\epsilon_{i,t}" /></a> take two values, respectively. And probability of being (un)employed depends on <a href="https://www.codecogs.com/eqnedit.php?latex=z_t" target="_blank"><img src="https://latex.codecogs.com/gif.latex?z_t" title="z_t" /></a>

### Recursive form of problem
<a href="https://www.codecogs.com/eqnedit.php?latex=V(\epsilon,k,Z)&space;=&space;\max_{c,k'}&space;\ln(c_{i,t})&space;&plus;&space;\beta&space;V(\epsilon',k',Z')\\&space;s.t.&space;\\&space;c&space;&plus;&space;k'&space;=&space;rk&space;&plus;&space;(1-\tau)&space;w&space;\bar&space;l&space;\epsilon&plus;\mu&space;w&space;(1-\epsilon)&space;&plus;&space;(1-\delta)k&space;\\&space;k'&space;\leq&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?V(\epsilon,k,Z)&space;=&space;\max_{c,k'}&space;\ln(c_{i,t})&space;&plus;&space;\beta&space;V(\epsilon',k',Z')\\&space;s.t.&space;\\&space;c&space;&plus;&space;k'&space;=&space;rk&space;&plus;&space;(1-\tau)&space;w&space;\bar&space;l&space;\epsilon&plus;\mu&space;w&space;(1-\epsilon)&space;&plus;&space;(1-\delta)k&space;\\&space;k'&space;\leq&space;0" title="V(\epsilon,k,Z) = \max_{c,k'} \ln(c_{i,t}) + \beta V(\epsilon',k',Z')\\ s.t. \\ c + k' = rk + (1-\tau) w \bar l \epsilon+\mu w (1-\epsilon) + (1-\delta)k \\ k' \leq 0" /></a>

where firms maximize profit give:

<a href="https://www.codecogs.com/eqnedit.php?latex=r&space;=&space;z&space;\alpha&space;(\frac{K}{\bar&space;l&space;(1-u(z))})^{\alpha&space;-1&space;}&space;\\&space;w&space;=&space;z&space;(1-\alpha)&space;(\frac{K}{\bar&space;l&space;(1-u(z))})^{\alpha}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?r&space;=&space;z&space;\alpha&space;(\frac{K}{\bar&space;l&space;(1-u(z))})^{\alpha&space;-1&space;}&space;\\&space;w&space;=&space;z&space;(1-\alpha)&space;(\frac{K}{\bar&space;l&space;(1-u(z))})^{\alpha}" title="r = z \alpha (\frac{K}{\bar l (1-u(z))})^{\alpha -1 } \\ w = z (1-\alpha) (\frac{K}{\bar l (1-u(z))})^{\alpha}" /></a>

and government's budget constraint is:

<a href="https://www.codecogs.com/eqnedit.php?latex=\tau&space;w&space;\bar&space;l&space;(1-u(z))&space;=&space;\mu&space;w&space;u(z)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\tau&space;w&space;\bar&space;l&space;(1-u(z))&space;=&space;\mu&space;w&space;u(z)" title="\tau w \bar l (1-u(z)) = \mu w u(z)" /></a>

There's many numerical methods we can choose to solve such problem, in the PS 1 & 2, we're practicing to be familiar with VFI, grid-methods, which is a bit slow. We can also implement many others methods, i.e. in my [final project](https://github.com/zhouweimin233/QuantMacro/tree/master/Part_B/Final_Project), I use Reiter method to solve K&S model, and also talk about some more recent developments.
