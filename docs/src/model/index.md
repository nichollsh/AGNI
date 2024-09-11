# Model description
AGNI models a planetary atmosphere by treating it as a single column (1D) and splitting it up into levels of finite thickness. These levels are defined in pressure-space, and are arranged logarithmically between the surface and the top of the atmosphere. The atmosphere is assumed to be plane-parallel. Quantities such as pressure and temperature are calculated at level-centres and level-edges, while energy fluxes are calculated only at the edges, and thermodynamic properties (e.g. heat capacity) are calculated only at their centres.

## Radiative transfer
Radiative transfer (RT) refers to the transport of radiation energy through a medium subject to the characteristics of the medium. Radiation passing through an atmosphere is absorbed, emitted, scattered, and reflected. In the context of planetary atmospheres, we also have to handle their surfaces, cloud formation, and radiation from the host star.

AGNI simulates RT using SOCRATES, a numerical code written by the UK Met Office which solves the RT equation using a two-stream solution. SOCRATES is accessed using a Julia interface originally written by Stuart Daines. The atmosphere is assumed to be hydrostatically supported and to behave as an ideal gas. Opacity is handle using the correlated-k approximation, with either random overlap or equivalent extinction used to account for overlapping absorption in mixtures of gases.

The model uses k-terms fitted to spectral absorption cross-section data from [DACE](https://dace.unige.ch/opacityDatabase/?#). The MT_CKD model is used to estimate water continuum absorption cross-sections. Other continuua are derived from the HITRAN tables. Rayleigh scattering and water cloud radiative properties are also included. You can find tools for fitting k-terms and processing line absorption data in my redistribution of [SOCRATES](https://github.com/nichollsh/SOCRATES) on GitHub. The flowchart below outlines how these absorption data are converted into a 'spectral file'.
```@raw html
  <img src="assets/spectral_flowchart.svg" width=100% class="center"/>
```
Surface reflectivity can be modelled as a greybody with an albedo from 0 to 1. Alternatively, it can be modelled from empirical single-scattering data which varies with zenith angle and wavelength.

## Convection
Convection is a turbulent process which occurs across more than one spatial dimension, so it must be parameterised within 1D models like AGNI. In fact, it is typically parameterised inside 3D global circulation models as resolving convection is numerically expensive. AGNI uses mixing length theory (MLT) to parameterise convection. This is in contrast to convective adjustment, which forcibly adjusts a convectively unstable region of the atmosphere to the corresponding adiabat while ensuring that enthalpy is conserved.

MLT directly calculates the energy flux associated with convective heat transport, and thus is the preferred parameterisation within the model. It assumes that parcels of gas are diffused over a characteristic _mixing length_, transporting energy in the process. This requires choosing a scale for this mixing length, but in practice this has very little impact on the results from the model.

Heat capacities are temperature-dependent, using values derived from the JANAF database.

## Phase change
Gases release energy ("latent heat" or "enthalpy") into their surroundings when condensing into a liquid or solid. This is included in the model
through a diffusive condensation scheme, which assumes a fixed condensation timescale. This takes place as follows... firstly,
the mixing ratios of the gases are updated according to the temperature profile, where rainout occurs until all condensibles are
saturated or sub-saturated. The mixing ratios of dry species are increased in order to satisfy the total pressure at condensing
levels. The heat released associated with the change in partial pressure of condensible gases is used to calculate a latent
heating rate. This is then integrated (from the TOA downwards) to provide a latent heat transport flux at cell-edges.

Latent heats are temperature-dependent, using values derived from Coker (2007) and Wagner & Pruß (2001).

## Solar flux
A key input to the radiation model is the shortwave downward-directed flux from the star at the top of the atmosphere. This is quantified by the bolometric instellation flux, a scale factor, an artificial additional albedo factor, and a zenith angle. All of these may be provided to the model through the configuration file. The model also requires a stellar spectrum scaled to the top of the atmosphere.

## Obtaining a solution

### Summary
AGNI is designed for modelling planetary atmospheres with high surface pressures and temperatures. This means that the radiative timescale differs by several orders of magnitude across the column, which makes obtaining a solution difficult. To obtain a temperature structure solution that conserves energy more precisely than a time-stepping method, AGNI solves for the temperature structure of the atmosphere as an optimisation problem. This finds the state which conserves energy across all levels and satisfies the required configuration from the user.

### Solution types
It is necessary to tell AGNI what kind of atmospheric solution to solve for. There are currently a few options available set by the `solution_type` variable in the configuration file.
* (1) Aim to conserve energy fluxes throughout the column. The surface temperature is fixed.
* (2) Aim to conserve energy fluxes throughout the column. The surface temperature is set by energy transport through a solid conductive boundary layer of thickness $d$ such that $T_s = T_m - \frac{Fd}{k}$, where $T_m$ is the mantle temperature and $k$ is the thermal conductivity.
* (3) Solve for a state such that the flux carried at each level is equal to $F_{\text{eff}} = \sigma T_{\text{eff}}^4$, representing the rate at which a planet is losing energy into space.

### Construction
The atmosphere is constructed of $N$ levels (cell-centres), corresponding to $N+1$ interfaces (cell-edges). The RT model takes cell-centre temperatures $T_i$, pressures $p_i$, geometric heights, and mixing ratios as input variables at each level $i$. As well as the surface temperature and incoming stellar flux. In return, it provides cell-edges spectral fluxes $F_i$ at all $N+1$ interfaces for LW & SW components and upward & downward streams. Convective fluxes can be estimated using the MLT scheme, condensation fluxes from the condensation scheme, and sensible heat from a simple TKE approximation.

The total upward-directed energy flux $F_{i}$ describes the total upward-directed energy transport (units of $\text{W m}^{-2}$) from cell $i$ into cell $i-1$ above (or into space for $i=1$). For energy to be conserved throughout the column, it must be true that $F_i = F_t \text{ } \forall \text{ } i$ where $F_t$ is the total amount of energy being transported out of the planet. In global radiative equilibrium, $F_t = 0$.

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
where $T_s$ is the surface temperature. Cell-edge temperatures in the bulk atmosphere are interpolated using the [PCHIP algorithm](https://epubs.siam.org/doi/10.1137/0905021), which avoids Runge's phenomenon and overshooting. The bottom- and top-most cell edge temperatures are extrapolated by estimation of $dT/d \log p$. Cell properties (heat capacity, gravity, density, average molecular weight, etc.) are consistently updated at each evaluation of $\bm{r}$. Condensation/rainout are also handled at each evaluation of $\bm{r}$ in order to avoid supersaturation.

The model converges when the cost function $c(\bm{x}) = \sqrt{\sum_i |r_i|}$ satisfies the condition
```math
c(\bm{x}) < c_a + c_r \cdot \underset{i}{\max} \text{ } |F_i|
```
which represents a state where the fluxes are sufficiently conserved.

### Iterative steps
The model solves for $\bm{x}$ iteratively, starting from some initial guess. The initial guess should be any reasonable temperature profile which is not significantly cooler than the expected solution. The flowchart below broadly outlines the solution process.
```@raw html
  <img src="assets/model_flowchart.svg" width=50% class="center"/>
```
The Jacobian matrix $\bm{J}$ represents the directional gradient of the residuals with respect to the solution vector. It is a square matrix with elements set according to
```math
J_{uv} = \frac{\partial r_u}{\partial x_v}
```
AGNI estimates $\bm{J}$ using finite-differences, requiring $N+1$ evalulations of $\bm{r}$ in order to fill the matrix. This corresponds to $2(N+1)+1$ objective function calculations under a 2nd order central-difference scheme. Each level $v$ with temperature $x_v$ is perturbed by an amount $\pm \varepsilon x_v$ in order to fill a single column of $\bm{J}$. As such, it can be expensive to construct a full Jacobian, especially when it is discarded at the end of each iteration. To reduce the total number of calculations AGNI retains some of the columns in $\bm{J}$ between model iterations. This assumes that the second derivative of the residuals is small. A column $v$ is retained only when
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

## Julia and Fortran
AGNI is primarily written in Julia, while SOCRATES itself is written in Fortran. Julia was chosen because it allows the SOCRATES binaries to be included in the precompiled code, which significantly improves performance.
