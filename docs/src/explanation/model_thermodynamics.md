# Thermodynamics and EOS

## Gas densities

The density of the gas mixture is calculated using Amagat's additive volume law to combine the densities of the individual components. You can read more about the validity and usage of this formulation here: [magyar_eos_2014, bradley_experimental_2018, magyar_mixing_2013, magyar_ethane_2015](@citet).

Densities of individual gas components are calculated using one of several equations of state:

* The **Van der Waals** equation of state is the default for most gases, accounting for intermolecular attractions and excluded volume [Kontogeorgis_TakingAnot_2019](@citep).
* **AQUA** is implemented as the EOS for water across a wide range of pressures and temperatures, including supercritical states [haldemann_aqua_2020](@citep).
* The [chabrier_eos_2019](@citet) EOS is used for hydrogen.
* For any gas not covered by the above, AGNI falls back to the ideal gas EOS: $\rho = p \mu / (R T)$.

The ideal gas approximation is accurate at low pressures and high temperatures; deviations become important at conditions found deep in the atmospheres of sub-Neptune exoplanets or near phase boundaries.

## Heat capacity

The molar heat capacity at constant pressure $c_p$ determines how much energy is required to raise the temperature of the gas. Its value increases with temperature as rotational and vibrational degrees of freedom in polyatomic molecules become accessible[pierrehumbert_book_2010](@citep). AGNI therefore implements temperature-dependent $c_p(T)$ using tabulated experimental data from the JANAF thermochemical tables [JANAF](@citep), accessed by interpolation.

## Saturation pressure and condensation

The saturation partial pressure $p^\text{sat}(T)$ defines the coexistence curve between vapour and condensate phases. The **Clausius–Clapeyron relation** describes how the saturation pressure varies with temperature along the phase coexistence curve:
```math
\frac{d p^\text{sat}}{d T} = \frac{p^\text{sat} L(T)}{R T^2}
```
where $L(T)$ is the molar latent heat of the relevant phase change. Integrating this from a known reference point gives the saturation curve.

The reference point is anchored using the Antoine equation,
```math
\log_{10}\!\left(\frac{p^\text{sat}}{\text{bar}}\right) = A - \frac{B}{T + C}
```
with coefficients from NIST, evaluated within the Antoine equation's valid temperature range.

## Latent heat

The latent heat $L(T)$ is the enthalpy difference between the vapour and condensate phases at constant pressure. $L$ decreases with increasing temperature and vanishes at the critical point $T_\text{crit}$, where the liquid and vapour phases become indistinguishable.

AGNI implements temperature-dependent $L(T)$ via tabulated empirical reference data. Sources are:
- Water: [IAPWS95](@citep) and [feistel_ice_2006](@citep) (below the triple point)
- Other species: [coker_thermo_2007](@citep)

For gases not covered by these data, $L = 0$ is assumed (i.e. no condensation). Data are interpolated at runtime.

See [Composition and chemistry](@ref) for how these thermodynamic properties are applied to the condensation and phase-change scheme.

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
