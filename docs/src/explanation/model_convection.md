# Atmospheric convection

## Physical background

Convection is a turbulent process that transfers heat by bulk fluid motion. In a planetary atmosphere, a parcel of gas that is warmer (less dense) than its surroundings will rise, expand, and cool adiabatically. If it remains warmer than the environment throughout its ascent, buoyancy drives continued upward motion (the atmosphere is **convectively unstable**) and heat is transported upward. The atmosphere is **convectively stable** in regions where energy transported primarily by radiation.

Convection can be driven by internal heat production or the absorption of shortwave stellar radiation in the deep atmosphere. Atmospheres are not always unstable to convection, leaving a stable radiative layers [selsis_cool_2023, nicholls_convective_2025](@citep). The general pattern of a convective troposphere below a radiative stratosphere is found across Solar System bodies with thick atmospheres and may hold broadly for planetary atmospheres [robinson_common_2014](@citep). AGNI does not enforce this structure *a priori*; the solver determines which layers are convective and which are radiative self-consistently.

## Dry adiabat and Schwarzschild criterion

For an ideal gas parcel undergoing dry (non-condensing) adiabatic ascent, the **dry adiabatic lapse rate** is
```math
\nabla_\text{ad} := \left.\frac{d \ln T}{d \ln p}\right|_\text{ad} = \frac{R}{\mu \, c_{pm}}
```
where $R = 8.314\ \mathrm{J\ mol^{-1}\ K^{-1}}$ is the universal gas constant, $\mu$ is the local mean molecular weight of the gas, and $c_{pm}$ is the molar heat capacity at constant pressure.

The **Schwarzschild criterion** for convective instability is then
```math
\left|\frac{d \ln T}{d \ln p}\right|_\text{atm} > \left|\frac{d \ln T}{d \ln p}\right|_\text{ad}
```

In this case, the atmosphere becomes unstable when its lapse rate exceeds that of an adiabatically rising parcel [gabriel_schwarz_2014, anders_schwarz_2022](@citep). This criterion applies in the absence of composition gradients.

The more general **Ledoux criterion** reduces the tendency to convect due to MMW gradients [gabriel_schwarz_2014](@citep). AGNI supports checking against either criterion (configurable via `physics.convection`).

## Mixing-length theory

Atmospheric convection is fundamentally multidimensional, so it must be parameterised in a 1D column model. AGNI applies **mixing-length theory (MLT)** to directly calculate the convective heat flux [robinson_temperature_2014, lee_dynamically_2024, joyce_mlt_2023](@citep). MLT represents convection by a single characteristic eddy: a rising parcel of gas that diffuses its excess enthalpy over a **mixing length** $\lambda$ before dissolving back into the environment. This is analogous to treating only the dominant wavenumber of turbulence [canuto_convection_1991](@citep). MLT remains a common and well-validated approach in the literature [marley_review_2015](@citep).

The convective heat flux is [nicholls_convective_2025](@citep):
```math
F^\text{cvt} = \frac{1}{2} \rho \, c_{pm} \, w \, T \, \frac{\lambda}{H} \, (\nabla_T - \nabla_\text{ad})
```
where $T$ is the local temperature, $c_{pm}$ is the heat capacity per unit mass, $H = RT/(\mu g)$ is the pressure scale height, $\nabla_T = d\ln T / d\ln p$ is the local atmospheric lapse rate, and
```math
w = \lambda \sqrt{\frac{g}{H} \, (\nabla_T - \nabla_\text{ad})}
```
is the characteristic convective velocity. The convective flux is positive (upward) wherever $\nabla_T > \nabla_\text{ad}$.

The mixing length follows a near-surface formulation:
```math
\lambda = \frac{k_v z}{1 + k_v z / H}
```
where $k_v \approx 1/\sqrt{2\pi}$ is the von Kármán constant and $z$ is the height above the surface [blackadar_mlt_1962, hogstrom_karman_1988](@citep). This ensures $\lambda \to 0$ as $z \to 0$, respecting the Law of the Wall for the near-surface turbulent boundary layer, and $\lambda \to H$ aloft. In practice the resulting temperature structure and convective flux are not sensitive to the specific choice of $\lambda$ parametrisation.

## Eddy diffusion coefficient $K_{zz}$

AGNI calculates the vertical eddy diffusion coefficient $K_{zz}$ (units of m² s⁻¹) to parameterise vertical mixing processes in the atmosphere. This coefficient is used by chemical kinetics codes (e.g. VULCAN). Three parametrisations are available, selectable via the `Kzz_type` configuration parameter:

**Type 1: Constant value**
```math
K_{zz} = K_{\text{break}}
```
where $K_{\text{break}}$ is a user-specified constant diffusion coefficient.

**Type 2: Simple MLT scaling** (default)
```math
K_{zz} = \lambda_{\text{conv}} \, w_{\text{conv}}
```
where $\lambda_{\text{conv}}$ is the mixing length and $w_{\text{conv}}$ is the convective velocity from MLT, calculated only in convective regions.

**Type 3: MLT convective flux scaling**

Following Charnay et al. ([2015](https://iopscience.iop.org/article/10.1088/0004-637X/813/1/15)), Equation 16:
```math
K_{zz} = \frac{H_p}{3} \left(\frac{\lambda_{\text{conv}}}{H_p}\right)^{4/3} \left(\frac{R_{\text{gas}} F_c}{\mu \rho c_p}\right)^{1/3}
```
where $H_p$ is the pressure scale height, $R_{\text{gas}}$ is the gas constant, $F_c$ is the convective energy flux, $\mu$ is the mean molecular weight, $\rho$ is the density, and $c_p$ is the specific heat capacity at constant pressure.

In non-convective (radiative) regions, $K_{zz}$ is extended using a power-law scaling with pressure:
```math
K_{zz}(p) = K_{\text{ref}} \left(\frac{p}{p_{\text{ref}}}\right)^{\alpha}
```
where $K_{\text{ref}}$ and $p_{\text{ref}}$ are the reference diffusion coefficient and pressure at the convective region boundary, and $\alpha$ is the power-law index (default: $\alpha = -0.4$), following [lee_dynamically_2024](@citep). Below convective regions, $K_{zz}$ is held constant at its deepest convective value.

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
