# Atmospheric convection

## Convection physics

Convection is a turbulent process that occurs across more than one spatial dimension, so it must be parameterised within 1D models like AGNI. In fact, it is typically parameterised inside 3D global circulation models, as resolving convection is numerically expensive.

Atmospheric convection is driven by locally heating deeper layers of the column, driving a vertical temperature gradient to exceed the critical lapse rate ($d T/d z$) required for triggering convection. In the case of planetary atmospheres, radiation is the primary source of local heating. In the case of Earth: solar radiation passes through the atmosphere and is absorbed at the surface, this energy is transported to the lower layers of the atmosphere, and latent heat from $\rm H_2O$ evaporation., heating it, and driving convection. Energy is transported through the lowermost layers by combination of thermal radiation, sensible heating (usually turbulence and conduction).

The lowest layer of Earth's atmosphere, the troposphere, is shaped by convection occurring from the surface upwards [pierrehumbert_book_2010](@citep). Above Earth's troposphere is the stratosphere, where the atmosphere is thermally inverted (i.e. $d T/d z > 0 $) due to absorption of shortwave radiation by $\rm O_3$ (see [Radiative transfer](@ref)). The stratosphere is defined by convective stability; energy balance is then set mainly by radiative and conductive processes. More generally, stratospheres can also be near-isothermal or slightly uninverted depending on the balance of shortwave and longwave opacities to radiation [pierrehumbert_book_2010, guillot_2010](@citep). The configuration of deep convection below and radiative layers aloft extends to other Solar System bodies with thick atmospheres and may be a general trend for planetary atmospheres more broadly [robinson_common_2014](@citep). It has thus been common in the literature to assume that planetary atmospheres exhibit deep convective zones a-priori; this assumption has come under recent scrutiny [selsis_cool_2023](@citep).

The critical lapse rate required to trigger convection depends on the _type_ of convection involved. Dry convection is taken to be an adiabatic process, meaning that no energy or mass is exchanged between a rising/sinking parcel and the surrounding atmosphere during its motion [pierrehumbert_book_2010](@citep). Radiative heating can then only drive the atmospheric lapse rate up to that of a dry adiabat,
```math
    \nabla_\text{ad} := \frac{d \ln T}{d \ln p}\bigg|_\text{ad} =\frac{R}{c_p}
```
in the case of an ideal gas (see [PhD Thesis](https://www.h-nicholls.space/thesis.pdf) for derivation).  R = 8.314 J/kg/K is the ideal gas constant and $c_p$ is the molar heat capacity of the gas mixture. The Schwarzschild criterion for determining convective instability is then,
```math
    \Bigg|\frac{d \ln T}{d \ln p} \Bigg|_\text{atm} > \Bigg| \frac{d \ln T}{d \ln p} \Bigg|_\text{ad},
```
which also sees wide application in stellar astrophysics [gabriel_schwarz_2014, anders_schwarz_2022](@citep). The temperature profile $T(p)$ in a convective region will then tend towards the integral of the dry adiabatic lapse rate. We thus have an analytic form for $T(p)$.

## Mixing-length implementation

AGNI uses mixing length theory (MLT) to parameterise atmospheric convection. This is in contrast to convective adjustment, which forcibly adjusts a convectively unstable region of the atmosphere to the corresponding adiabat while ensuring that enthalpy is conserved.MLT directly calculates the energy flux associated with convective heat transport, and thus is the preferred parameterisation within the model. Canonical MLT neglects to resolve the full spectrum of turbulence associated with real convection, and instead implicitly considers only a single convective eddy, represented by a rising parcel of air which diffuses energy over some {mixing length} $\lambda$ [joyce_mlt_2023](@citep). This is analogous to only treating the dominant wavenumber of turbulence [canuto_convection_1991](@citep). With this assumption, we can derive an equation which allows a calculation of the energy flux associated with convection that is induced by the local lapse rate $\nabla_T$ exceeding that of an adiabat $\nabla_\text{ad}$.

This requires choosing a scale for this mixing length, but in practice this has very little impact on the results from the model.

When evaluating convective energy fluxes, AGNI first calculates the temperature gradient across each layer of the atmosphere. Convection occurs within each layer that has a lapse rate $dT/dP$ greater than the critical lapse rate for triggering convection. Equations 2 to 6 of [nicholls_convective_2025](@citet) describe the calculation of the convective energy flux under the Schwarzschild criterion. AGNI can also check against the Ledoux criterion for convective stability, which accounts for vertical gradients in the gas MMW. The atmosphere is not explicitly split into convecting and non-convecting regions, thereby allowing disconnected regions of convection.


## Eddy diffusion coefficient $K_{zz}$

AGNI calculates the vertical eddy diffusion coefficient $K_{zz}$ (units of m² s⁻¹) to parameterise vertical mixing processes in the atmosphere. Three parametrisations are available, selectable via the `Kzz_type` configuration parameter:

**Type 1: Constant value**
```math
K_{zz} = K_{\text{break}}
```
where $K_{\text{break}}$ is a user-specified constant diffusion coefficient.

**Type 2: Simple MLT scaling** (default)
```math
K_{zz} = \lambda_{\text{conv}} \, w_{\text{conv}}
```
where $\lambda_{\text{conv}}$ is the mixing length and $w_{\text{conv}}$ is the convective velocity, both calculated from mixing length theory in convective regions.

**Type 3: MLT convective flux scaling**

Following Charnay et al. ([2015](https://iopscience.iop.org/article/10.1088/0004-637X/813/1/15)), Equation 16:
```math
K_{zz} = \frac{H_p}{3} \left(\frac{\lambda_{\text{conv}}}{H_p}\right)^{4/3} \left(\frac{R_{\text{gas}} F_c}{\mu \rho c_p}\right)^{1/3}
```
where $H_p$ is the pressure scale height, $R_{\text{gas}}$ is the gas constant, $F_c$ is the convective energy flux, $\mu$ is the mean molecular weight, $\rho$ is the density, and $c_p$ is the specific heat capacity at constant pressure.

In convective regions, $K_{zz}$ is calculated directly using one of the above parametrisations. In non-convective regions (radiative zones), $K_{zz}$ is extended using a power-law scaling with pressure:
```math
K_{zz}(p) = K_{\text{ref}} \left(\frac{p}{p_{\text{ref}}}\right)^{\alpha}
```
where $K_{\text{ref}}$ and $p_{\text{ref}}$ are the reference diffusion coefficient and pressure at the convective region boundary (or at a user-specified breakpoint pressure), and $\alpha$ is the power-law index (default: $\alpha = -0.4$). This stratospheric scaling follows Tsai et al. ([2020](https://iopscience.iop.org/article/10.3847/1538-4357/ac29bc), Equation 28). Below convective regions, $K_{zz}$ is held constant at its deepest convective value.

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
