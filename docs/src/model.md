# Description 
AGNI models a planetary atmosphere by treating it as a single column (1D) and splitting it up into levels of finite thickness. These levels are defined in pressure-space, and are arranged logarithmically between the surface and the top of the atmosphere. Quantities such as pressure and temperature are calculated at level-centres and level-edges, while energy fluxes are calculated only at the edges, and thermodynamic properties (e.g. heat capacity) are calculated only at their centres.

## Radiative transfer
Radiative transfer (RT) refers to the transport of radiation energy through a medium subject to the characteristics of the medium. Radiation passing through an atmosphere is absorbed, emitted, scattered, and reflected. In the context of planetary atmospheres, we also have to handle their surfaces, cloud formation, and radiation from the host star.

AGNI models RT using SOCRATES, a numerical code written by the UK Met Office which solves the RT equation using a two-stream solution under a plane-parallel approximation. SOCRATES is accessed using a Julia interface originally written by Stuart Daines. The atmosphere is assumed to be hydrostatically supported and to behave as an ideal gas. Opacity is calculated using the correlated-k approximation, with either random overlap or equivalent extinction used to account for overlapping absorption in mixtures of gases. 

For simulating gaseous absorption, the model fits k-terms to spectral absorption cross-section data from DACE. The MT_CKD model is used to estimate continuum absorption cross-sections. Rayleigh scattering and water cloud radiative effects are also included. You can find tools for fitting k-terms and processing line absorption data in my redistribution of [SOCRATES](https://github.com/nichollsh/SOCRATES) on GitHub.

## Convection
Convection is a process that occurs across more than one spatial dimension, so it must be parameterised within 1D models like AGNI. In fact, it's often parameterised in 3D global circulation models, as resolving convection is numerically difficult. AGNI uses 
mixing length theory (MLT) to parameterise convection. This is in contrast to convective adjustment, which forcibly adjusts a convectively unstable region of the atmosphere to the corresponding adiabat, while ensuring that enthalpy is conserved. 
   
MLT directly calculates the energy flux associated with convective heat transport, and thus is the preferred parameterisation within the model. It assumes that parcels of gas are diffused over a characteristic _mixing length_, transporting energy in the process.
   
Heat capacities are temperature-dependent, calculated using the Shomate Equation with coefficients derived from the NIST website.

## Latent heat
Gases release energy (latent heat) into their surroundings when condensing into a liquid or solid. This is included in the model 
through a diffusive condensation scheme, which assumes a fixed condensation timescale. This takes place as follows... firstly, 
the mixing ratios of the gases are updated according to the temperature profile, where rainout occurs until all condensibles are 
saturated or sub-saturated. The mixing ratios of dry species are increased in order to satisfy the total pressure at condensing
levels. The heat released associated with the change in partial pressure of condensible gases is used to calculate a latent
heating rate. This is then integrated to provide a latent heat transport flux.

## Solar flux
The radiation component requires two boundary conditions the energy. The first is the shortwave downward-directed flux from the star at the top of the atmosphere. This is quantified by the instellation, a scale factor, a grey bond albedo, and the solar zenith angle. All of these may be provided to the model through the configuration file.

## Solution types
Depending on the system you wish to model, it is necessary to tell AGNI what kind of solution to solve for. There are currently a few options available set by the `solution_type` variable.   
* (0) Aim to conserve energy fluxes throughout the column. The surface temperature is set at $T_s$ assuming blackbody emission with a fixed surface albedo. The bottom-most temperature value in the column is extrapolated from the rest of the profile.
* (1) Same as 0, but the bottom-most temperature value is fixed equal to $T_s$.
* (2) Aim to conserve energy fluxes throughout the column. The surface temperature is set by energy transport through a solid conductive boundary layer (CBL) such that $T_s = T_m - Fd/k$, where $T_m$ is the interior mantle temperature, while $k$ and $d$ are material properties. 
* (3) Solve for a state such that the flux carried at each level is equal to $\sigma T_{\text{eff}}^4$. In this case, $T_{\text{eff}}$ represents the rate at which a planet is losing energy into space. 


## Obtaining a solution
AGNI is designed for modelling planetary atmospheres with high surface pressures and temperatures. This means that the radiative timescale differs by several orders of magnitude across the column, which makes obtaining a solution difficult. The model contains a suite of methods for obtaining a solution.
    
To obtain a temperature structure solution that conserves energy more precisely than a time-stepping method, it is possible to construct the model as a system of non-linear equations $\vec{r}(\vec{x})$. This algorithm obtains the roots of the nonlinear system formed by the flux divergence $r_i$ at each level $i$, with cell-centre temperatures used as the independent variables $x_i$. Finite-difference methods are used to estimate the Jacobian matrix $\bf{J}$. After the first iteration, columns of $\bf{J}$ are only updated when required. Similar methods have been used in a handful of planetary radiative-convective models (e.g. ATMO), but is more commonly used by the stellar physics community. Currently implemented solvers: Newton-Raphson, Gauss-Newton, and Levenberg-Marquardt. Optionally, the code uses a backtracking linesearch algorithm to determine the optimal step length.

## Other features
AGNI can calculate emission spectra, provided with T(p) and the volume mixing ratios of the gases. This is performed using the same RT as the RCE calculations, so is limited in resolution by the choice of correlated-k bands. Similarly, the longwave contribution function can also be calculated.

## Julia and Fortran
AGNI is primarily written in Julia, while SOCRATES itself is written in Fortran. Julia was chosen because it allows the SOCRATES binaries to be included in the precompiled code, which significantly improves performance.
