### Part B by [Prof. Luis Rojas](https://sites.google.com/site/luiserojasweb)

---
### Useful Links:

- [Wouter den Haan's Website](http://www.wouterdenhaan.com/papers.htm): matlab codes for solving Incomplete Markets model with aggregate uncertainty.

- [Econ 861 by Matteo Iacoviello](https://www2.bc.edu/matteo-iacoviello/teach/0910/EC861.html): a version of the Krusell-Smith model matlab codes.

---
#### PS1 and PS2 is about solving Models with HA - Krusell&Smith (1998)

- employment shocks (<a href="https://www.codecogs.com/eqnedit.php?latex=$\epsilon_{i,t}\in&space;\{0,1\}$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\epsilon_{i,t}\in&space;\{0,1\}$" title="$\epsilon_{i,t}\in \{0,1\}$" /></a>)

- Incomplete markets, with borrowing constraint <a href="https://www.codecogs.com/eqnedit.php?latex=k_{i,t&plus;1}\leq&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k_{i,t&plus;1}\leq&space;0" title="k_{i,t+1}\leq 0" /></a>

- transition probabilities calibrated by K&S(1998), where <a href="https://www.codecogs.com/eqnedit.php?latex=z_t$,&space;$\epsilon_{i,t}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?z_t$,&space;$\epsilon_{i,t}" title="z_t$, $\epsilon_{i,t}" /></a> take two values, respectively. And probability of being (un)employed depends on <a href="https://www.codecogs.com/eqnedit.php?latex=z_t" target="_blank"><img src="https://latex.codecogs.com/gif.latex?z_t" title="z_t" /></a>

### Recursive form of problem
<a href="https://www.codecogs.com/eqnedit.php?latex=V(\epsilon,k,Z)&space;=&space;\max_{c,k'}&space;\ln(c_{i,t})&space;&plus;&space;\beta&space;V(\epsilon',k',Z')\\&space;s.t.&space;\\&space;c&space;&plus;&space;k'&space;=&space;rk&space;&plus;&space;(1-\tau)&space;w&space;\bar&space;l&space;\epsilon&plus;\mu&space;w&space;(1-\epsilon)&space;&plus;&space;(1-\delta)k&space;\\&space;k'&space;\leq&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?V(\epsilon,k,Z)&space;=&space;\max_{c,k'}&space;\ln(c_{i,t})&space;&plus;&space;\beta&space;V(\epsilon',k',Z')\\&space;s.t.&space;\\&space;c&space;&plus;&space;k'&space;=&space;rk&space;&plus;&space;(1-\tau)&space;w&space;\bar&space;l&space;\epsilon&plus;\mu&space;w&space;(1-\epsilon)&space;&plus;&space;(1-\delta)k&space;\\&space;k'&space;\leq&space;0" title="V(\epsilon,k,Z) = \max_{c,k'} \ln(c_{i,t}) + \beta V(\epsilon',k',Z')\\ s.t. \\ c + k' = rk + (1-\tau) w \bar l \epsilon+\mu w (1-\epsilon) + (1-\delta)k \\ k' \leq 0" /></a>

where firms maximize profit give:

<a href="https://www.codecogs.com/eqnedit.php?latex=r&space;=&space;z&space;\alpha&space;(\frac{K}{\bar&space;l&space;(1-u(z))})^{\alpha&space;-1&space;}&space;\\&space;w&space;=&space;z&space;(1-\alpha)&space;(\frac{K}{\bar&space;l&space;(1-u(z))})^{\alpha}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?r&space;=&space;z&space;\alpha&space;(\frac{K}{\bar&space;l&space;(1-u(z))})^{\alpha&space;-1&space;}&space;\\&space;w&space;=&space;z&space;(1-\alpha)&space;(\frac{K}{\bar&space;l&space;(1-u(z))})^{\alpha}" title="r = z \alpha (\frac{K}{\bar l (1-u(z))})^{\alpha -1 } \\ w = z (1-\alpha) (\frac{K}{\bar l (1-u(z))})^{\alpha}" /></a>

and government's budget constraint is:

<a href="https://www.codecogs.com/eqnedit.php?latex=\tau&space;w&space;\bar&space;l&space;(1-u(z))&space;=&space;\mu&space;w&space;u(z)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\tau&space;w&space;\bar&space;l&space;(1-u(z))&space;=&space;\mu&space;w&space;u(z)" title="\tau w \bar l (1-u(z)) = \mu w u(z)" /></a>

There's many numerical methods we can choose to solve such problem, in the PS, we're practicing to be familiar with VFI, grid-methods, which is a bit slow. We can alos implement many others methods: 

### 1. VFI with evenly-spaced grids methods
See [Part_B/PS1-solution/main.m](https://github.com/zhouweimin233/QuantMacro/blob/master/Part_B/PS1-solution/main.m) 

### 2. Euler Equation iteration with projection methods

1. Just like previous algorithm, we still approximate Law of motion with our Perceived law of motion and update the parameters for simulation. 
2. create state variables nodes: <a href="https://www.codecogs.com/eqnedit.php?latex=s_{i,t}=\{\epsilon_{i,t},&space;k_{i,t},s_{t}\}$&space;with&space;$s_t&space;=&space;\{&space;z_t,&space;K_t,&space;m_t&space;\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s_{i,t}=\{\epsilon_{i,t},&space;k_{i,t},s_{t}\}$&space;with&space;$s_t&space;=&space;\{&space;z_t,&space;K_t,&space;m_t&space;\}" title="s_{i,t}=\{\epsilon_{i,t}, k_{i,t},s_{t}\}$ with $s_t = \{ z_t, K_t, m_t \}" /></a> where <a href="https://www.codecogs.com/eqnedit.php?latex=m_{t&plus;1}&space;=&space;\Gamma&space;(z_{t&plus;1},z_t,&space;s_t;&space;\eta_{\Gamma})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?m_{t&plus;1}&space;=&space;\Gamma&space;(z_{t&plus;1},z_t,&space;s_t;&space;\eta_{\Gamma})" title="m_{t+1} = \Gamma (z_{t+1},z_t, s_t; \eta_{\Gamma})" /></a>


3. approximate <a href="https://www.codecogs.com/eqnedit.php?latex=k_{i,t&plus;1}&space;(\cdot)&space;=&space;P_n&space;(s_{i,t};\eta_{P_n})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k_{i,t&plus;1}&space;(\cdot)&space;=&space;P_n&space;(s_{i,t};\eta_{P_n})" title="k_{i,t+1} (\cdot) = P_n (s_{i,t};\eta_{P_n})" /></a>

4. tranform the FOC: 

<a href="https://www.codecogs.com/eqnedit.php?latex=c_t&space;^{-\gamma}&space;=&space;E&space;[\beta&space;(r(z',K')&plus;(1-\delta))c_{t&plus;1}^{-\gamma}]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?c_t&space;^{-\gamma}&space;=&space;E&space;[\beta&space;(r(z',K')&plus;(1-\delta))c_{t&plus;1}^{-\gamma}]" title="c_t ^{-\gamma} = E [\beta (r(z',K')+(1-\delta))c_{t+1}^{-\gamma}]" /></a>

replace consumption by BC, FOC changes to:

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align}&space;\begin{split}&space;\{&space;r(z,K)&plus;(1-\delta)k&space;&plus;&space;(1-\tau(z))w(z,K)l\epsilon&space;\\&plus;&space;\mu&space;w(z,K)(1-\epsilon)-P_n(s;\eta_{P_n})&space;\}^{-\gamma}&space;\end{split}&space;\end{align}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{align}&space;\begin{split}&space;\{&space;r(z,K)&plus;(1-\delta)k&space;&plus;&space;(1-\tau(z))w(z,K)l\epsilon&space;\\&plus;&space;\mu&space;w(z,K)(1-\epsilon)-P_n(s;\eta_{P_n})&space;\}^{-\gamma}&space;\end{split}&space;\end{align}" title="\begin{align} \begin{split} \{ r(z,K)+(1-\delta)k + (1-\tau(z))w(z,K)l\epsilon \\+ \mu w(z,K)(1-\epsilon)-P_n(s;\eta_{P_n}) \}^{-\gamma} \end{split} \end{align}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex==E[\{&space;\begin{align*}&space;\beta(&space;r(z',K')&plus;(1-\delta))\times(r(z',K')&&plus;1-\delta)P_{n}(s;\eta_{P_n})&space;\\&space;&plus;(1-\tau(z'))w(z',K')l\epsilon'&space;&\\&plus;&space;\mu&space;w(z',K')(1-\epsilon')-P_n(s';\eta_{P_n})&space;\end{align*}&space;\}^{-\gamma}]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?=E[\{&space;\begin{align*}&space;\beta(&space;r(z',K')&plus;(1-\delta))\times(r(z',K')&&plus;1-\delta)P_{n}(s;\eta_{P_n})&space;\\&space;&plus;(1-\tau(z'))w(z',K')l\epsilon'&space;&\\&plus;&space;\mu&space;w(z',K')(1-\epsilon')-P_n(s';\eta_{P_n})&space;\end{align*}&space;\}^{-\gamma}]" title="=E[\{ \begin{align*} \beta( r(z',K')+(1-\delta))\times(r(z',K')&+1-\delta)P_{n}(s;\eta_{P_n}) \\ +(1-\tau(z'))w(z',K')l\epsilon' &\\+ \mu w(z',K')(1-\epsilon')-P_n(s';\eta_{P_n}) \end{align*} \}^{-\gamma}]" /></a>


where the expectation is transition matrix times each possible state would happen.

5. determine prices remains the same as before: 
$$ r(z,K) = z \alpha (\frac{K}{L(z)})^{\alpha -1 } \\
w = z (1-\alpha) (\frac{K}{L(z)})^{\alpha}
$$

6. by minimizing the error term (LHS-RHS) from Step.4: $\min \sum u^2_k$,we could find the best fitted parameters. 


### 3. [explicit aggregation](http://www.wouterdenhaan.com/numerical/methodsheteroxpa.pdf)
### ...
