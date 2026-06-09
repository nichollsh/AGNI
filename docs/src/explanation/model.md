# Model description
AGNI models a planetary atmosphere by treating it as a single column (1D) and splitting it up into levels of finite thickness. These levels are defined in pressure-space, and are arranged logarithmically between the surface and the top of the atmosphere. The atmosphere is assumed to be plane-parallel. Quantities such as pressure and temperature are calculated at level-centres and level-edges, while energy fluxes are calculated only at the edges, and thermodynamic properties (e.g. heat capacity) are calculated only at their centres.

## Physics

You can find a detailed description of the physics underpinning AGNI, and the specific implementation within the code, by reading the pages in this section of the documentation. These are linked in the navbar on the left-hand side of the website.

The physics descriptions and implementational approach are largely derived from Harrison Nicholls' [PhD Thesis](https://www.h-nicholls.space/thesis.pdf). Literature references are provided in each page's bibliography section.

In-code references are also provided in the [Bibliography](@ref).

## Obtaining a solution

### Summary
AGNI is designed for modelling planetary atmospheres with high surface pressures and temperatures. This means that the radiative timescale differs by several orders of magnitude across the column, which makes obtaining a solution difficult. To obtain a temperature structure solution that conserves energy more precisely than a time-stepping method, AGNI solves for the temperature structure of the atmosphere as an optimisation problem. This finds the state which conserves energy across all levels and satisfies the required configuration from the user.

### Solution types
It is necessary to tell AGNI what kind of atmospheric solution to solve for. There are currently a few options available set by the `solution_type` variable in the configuration file.
* (1) Aim to conserve energy fluxes throughout the column. The surface temperature is fixed.
* (2) Aim to conserve energy fluxes throughout the column. The surface temperature is set by energy transport through a solid conductive boundary layer of thickness $d$ such that $T_s = T_m - \frac{Fd}{k}$, where $T_m$ is the mantle temperature and $k$ is the thermal conductivity.
* (3) Solve for a state such that the flux carried at each level is equal to $F_\text{int} = \sigma T_\text{int}^4$, representing the rate at which a planet is losing energy into space.
* (4) Solve for a state such that the outgoing longwave radiation (OLR) at the top of atmosphere matches a user-specified target value, while conserving energy throughout the column.

Solution type (1) enforces a fixed surface temperature and instellation, allowing all other temperatures and $F_\text{int}$ to be solved-for as dependent variables. This could be used to model a young planet far from 'radiative equilibrium' (or 'global energy balance') where its surface temperature is very large and the outgoing energy flux is non-zero.

Solution type (2) is appropriate for coupling with a magma ocean model, where a conductive skin of solidified rock forms at the atmosphere-mantle interface. This skin is a kind of boundary layer that must be parameterised as having a particular thickness and thermal conductivity. It's thought that these layers are important for regulating the energy budget of young rocky planets. This is similar to type (1) except that $T_s$ is allowed to change and $T_m$ is fixed.

Solution type (3) is comparable to the implementation inside other atmosphere climate models: the surface temperature is free to change along with the rest of the atmosphere, and AGNI solves for a state of radiative equilibrium.

Solution type (4) is useful when the target OLR is known from observational constraints or when exploring climate states for a specific energy balance condition.


### Construction
The atmosphere is constructed of $N$ levels (cell-centres), corresponding to $N+1$ interfaces (cell-edges). The RT model takes cell-centre temperatures $T_i$, pressures $p_i$, geometric heights, and mixing ratios as input variables at each level $i$, as well as the surface temperature and incoming stellar flux. In return, it provides cell-edges spectral fluxes $F_i$ at all $N+1$ interfaces for LW & SW components and upward & downward streams. Convective fluxes can be estimated using the MLT scheme, condensation fluxes from the condensation scheme, and sensible heat from a simple turbulent kinetic energy (TKE) approximation.


### Definition of residuals
We can use this construction to solve for the temperature profile of the atmosphere as an $N+1$-dimensional optimisation problem. This directly solves for $T(p)$ at radiative-convective equilibrium without having to invoke heating rate calculations, thereby avoiding slow convergence in regions of the atmosphere with long radiative timescales. The residuals vector (length $N+1$)
```math

\bm{r} =

\begin{pmatrix}
r_i     \\
r_{i+1} \\
...     \\
r_N     \\
r_{N+1}
\end{pmatrix}

=

\begin{pmatrix}
F_{i+1} - F_i     \\
F_{i+2} - F_{i+1} \\
...     \\
F_{N+1} - F_N     \\
F_{N+1} - F_t
\end{pmatrix}
```
is what we aim to minimise as our 'objective function', subject to the solution vector of cell-centre temperatures
```math
\bm{x} =

\begin{pmatrix}
T_i     \\
T_{i+1} \\
...     \\
T_N     \\
T_s
\end{pmatrix}
```
where $T_s$ is the surface temperature. Cell-edge temperatures in the bulk atmosphere are interpolated from cell-centres. The bottom- and top-most cell edge temperatures are extrapolated by estimation of $dT/d \log p$. Cell properties (heat capacity, gravity, density, average molecular weight, etc.) are consistently updated at each evaluation of $\bm{r}$. Condensation/rainout are also handled at each evaluation of $\bm{r}$ in order to avoid supersaturation.

The model converges when the cost function $c(\bm{x}) = \sqrt{\sum_i |r_i|}$ satisfies the condition
```math
c(\bm{x}) < c_a + c_r \cdot \underset{i}{\max} \text{ } |F_i|
```
which represents a state where the fluxes are sufficiently conserved. The quantities $c_a$ and $c_r$ are the absolute and relative tolerances provided by the user (parameters `conv_atol` and `conv_rtol` in `solver.solve_energy!()`).

### Iterative steps
The model solves for $\bm{x}$ iteratively, starting from some initial guess. The initial guess should be any reasonable temperature profile which is not significantly cooler than the expected solution. The flowchart below broadly outlines the solution process.

![](fig_model_flowchart.svg)


The Jacobian matrix $\bm{J}$ represents the directional gradient of the residuals with respect to the solution vector. It is a square matrix with elements set according to
```math
J_{uv} = \frac{\partial r_u}{\partial x_v}
```
AGNI estimates $\bm{J}$ using finite-differences, requiring $N+1$ evaluations of $\bm{r}$ in order to fill the matrix. This corresponds to $2(N+1)+1$ objective function calculations under a 2nd order central-difference scheme. Each level $v$ with temperature $x_v$ is perturbed by an amount $\pm \varepsilon x_v$ in order to fill a single column of $\bm{J}$. As such, it can be expensive to construct a full Jacobian, especially when it is discarded at the end of each iteration. To reduce the total number of calculations, AGNI retains some of the columns in $\bm{J}$ between model iterations. This assumes that the second derivative of the residuals is small. A column $v$ is retained only when
```math
\max |r_i| \lt 0.7 \text{ for } i \in \{v-1, v, v+1\}
```
and when $c(\bm{x})/10$ does **not** satisfy the convergence criteria.

With a Jacobian constructed, we can calculate an update $\bm{d}$ to the solution vector $\bm{x} \rightarrow \bm{x} + \bm{d}$. This is primarily done via the Newton-Raphson method
```math
\bm{d} = -\bm{J}^{-1} \bm{r}
```
but can alternatively be performed via the Gauss-Newton and Levenberg–Marquardt methods.

It is possible for the model to become stuck in a local minimum, leading to very small values of $|\bm{d}|$. This is identified when $c(\bm{x})$ has seen little change over the last few iterations. When this occurs, the model is 'nudged' by scaling the update via
```math
\bm{d} \rightarrow 3 \bm{d}
```
In many cases $\bm{d}$ is too large, leading to instabilities. This is due to the non-convexity of the solution space, and the somewhat discontinuous nature of the physics involved (particularly in its temperature derivatives). When $|\bm{d}|>d_{\text{max}}$, the update is crudely scaled via
```math
\bm{d} \rightarrow  d_{\text{max}} \hat{\bm{d}}
```
The update may also be scaled by a linesearch
```math
\bm{d} \rightarrow \alpha \bm{d}
```
on $\alpha$. This is applied if the full step $\bm{d}$ would increase the cost by an unacceptable amount. If the model is close to convergence then a golden-section search method is used to determine the optimal $\alpha$, otherwise a backtracking method is used.  This means that the model is (mostly) able to avoid oscillating around a solution. All three of these scalings to $\bm{d}$ preserve its direction.

## Other features
AGNI can calculate emission spectra, provided with T(p) and the volume mixing ratios of the gases. This is performed using the same RT as the RCE calculations, so is limited in resolution by the choice of correlated-k bands. Similarly, the longwave contribution function can also be calculated.

## Transparent atmospheres
It is useful to run AGNI with a transparent atmosphere in various scenarios. For example, in the calculation of reflectance or emission spectra of 'bare rock' planets. Or alternatively to determine a planet's surface temperature in the absence of an overlying atmosphere. AGNI incorporates this functionality through the configuration variable `composition.transparent=true`. This will set the atmosphere surface pressure to be small and disable the gas opacity, continuum opacity, and other absorption processes in SOCRATES.

When this is enabled, make sure to use the "transparent" solver. The section below does not apply in this case, as there is only one variable (the surface temperature) to solve for - this is done using a Golden-Section search method.

## Julia and Fortran
AGNI is primarily written in Julia, while SOCRATES itself is written in Fortran. Julia was chosen because it allows the SOCRATES binaries to be included in the precompiled code, which significantly improves performance.
