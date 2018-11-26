### Part B by [Prof. Luis Rojas](https://sites.google.com/site/luiserojasweb)

---
### Useful Links:

- [Wouter den Haan's Website](http://www.wouterdenhaan.com/papers.htm): matlab codes for solving Incomplete Markets model with aggregate uncertainty.

- [Econ 861 by Matteo Iacoviello](https://www2.bc.edu/matteo-iacoviello/teach/0910/EC861.html): a version of the Krusell-Smith model matlab codes.

---
#### PS1 and PS2 is about solving Models with HA - Krusell&Smith (1998)

- employment shocks ($\epsilon_{i,t}\in \{0,1\}$)
- Incomplete markets, with borrowing constraint $k_{i,t+1}\leq 0$

- transition probabilities calibrated by K&S(1998), where $z_t$, $\epsilon_{i,t}# take two values, respectively. And probability of being (un)employed depends on $z_t$

### Recursive form of problem
$$V(\epsilon,k,Z) = \max_{c,k'} \ln(c_{i,t}) + \beta V(\epsilon',k',Z')\\
s.t.  \\
c + k' = rk + (1-\tau) w \bar l \epsilon+\mu w (1-\epsilon) + (1-\delta)k \\ 
k' \leq 0 $$
where firms maximize profit give:
$$ 
r = z \alpha (\frac{K}{\bar l (1-u(z))})^{\alpha -1 } \\
w = z (1-\alpha) (\frac{K}{\bar l (1-u(z))})^{\alpha}
$$
and government's budget constraint is:
$$ \tau w \bar l (1-u(z)) = \mu w u(z)$$

There's many numerical methods we can choose to solve such problem, in the PS, we're practicing to be familiar with VFI, grid-methods, which is a bit slow. We can alos implement many others methods: 

### 1. VFI with evenly-spaced grids methods
See [Part_B/PS1-solution/main.m](https://github.com/zhouweimin233/QuantMacro/blob/master/Part_B/PS1-solution/main.m) 

### 2. Euler Equation iteration with projection methods

1. Just like previous algorithm, we still approximate Law of motion with our Perceived law of motion and update the parameters for simulation. 
2. create state variables nodes: $s_{i,t}=\{\epsilon_{i,t}, k_{i,t},s_{t}\}$ with $s_t = \{ z_t, K_t, m_t \}$ where $m_{t+1} = \Gamma (z_{t+1},z_t, s_t; \eta_{\Gamma})$


3. approximate $k_{i,t+1} (\cdot) =  P_n (s_{i,t};\eta_{P_n}) $

4. tranform the FOC: 
$$ c_t ^{-\gamma} = E [\beta (r(z',K')+(1-\delta))c_{t+1}^{-\gamma}]$$
replace consumption by BC, FOC changes to:
$$ \{
\begin{align}
 r(z,K)+(1-\delta)k + (1-\tau(z))w(z,K)l\epsilon \\+ \mu w(z,K)(1-\epsilon)-P_n(s;\eta_{P_n})
\end{align}
\}^{-\gamma} $$
$$
=E[\{
\begin{align}
 \beta( r(z',K')+(1-\delta))\times(r(z',K')&+1-\delta)P_{n}(s;\eta_{P_n}) \\ +(1-\tau(z'))w(z',K')l\epsilon' &\\+ \mu w(z',K')(1-\epsilon')-P_n(s';\eta_{P_n})
\end{align}
 \}^{-\gamma}]
$$

where the expectation is transition matrix times each possible state would happen.

5. determine prices remains the same as before: 
$$ r(z,K) = z \alpha (\frac{K}{L(z)})^{\alpha -1 } \\
w = z (1-\alpha) (\frac{K}{L(z)})^{\alpha}
$$

6. by minimizing the error term (LHS-RHS) from Step.4: $\min \sum u^2_k$,we could find the best fitted parameters. 


### 3. [explicit aggregation](http://www.wouterdenhaan.com/numerical/methodsheteroxpa.pdf)
### ...
