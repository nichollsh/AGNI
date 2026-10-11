# Thermodynamics and EOS

## Gas densities

The density of the gas mixture is calculated using Amagat's additive volume law to combine the densities of the individual components. You can read more about the validity and usage of this formulation here: [magyar_eos_2014, bradley_experimental_2018, magyar_mixing_2013, magyar_ethane_2015](@citet).

The mass mixing ratio $q_j$ of each species $j$ is converted from its volume mixing ratio $x_j$ via $q_j = x_j \mu_j / \mu$.

The mixture density is then found from
```math
\frac{1}{\rho_\text{mix}} = \sum_j \frac{q_j}{\rho_j(T,p)}
```
That is, the *specific volumes* of the individual components are added and weighted by their mass fraction, with every component's density $\rho_j$ evaluated at the full mixture temperature and total pressure (**not its own partial pressure**).

Densities of individual gas components are calculated using one of several equations of state:

* The **Van der Waals** equation of state is the default for most gases, accounting for intermolecular attractions and excluded volume [Kontogeorgis_TakingAnot_2019](@citep).
* **AQUA** is implemented as the EOS for water across a wide range of pressures and temperatures, including supercritical states [haldemann_aqua_2020](@citep).
* The [chabrier_eos_2019](@citet) EOS is used for hydrogen.
* For any gas not covered by the above, AGNI falls back to the ideal gas EOS: $\rho = p \mu / (R T)$.

The ideal gas approximation is accurate at low pressures and high temperatures; deviations become important at conditions found deep in the atmospheres of sub-Neptune exoplanets or near phase boundaries.

Real-gas equations of state are tabulated in $(T, \log_{10}p)$ space and interpolated at runtime. Evaluating a tabulated EOS exactly at a phase boundary is numerically discontinuous, since the table will transition abruptly from vapour- and condensate-values for density. AGNI handles this with a configurable scheme (`phs_method`) that either evaluates the table directly, falls back to the ideal-gas law on the condensed side, or uses a 'metastable EOS' approach which evaluates the EOS at a small distance (in $\log_{10}p$ space) from the saturation curve and extrapolates smoothly into the vapour region, so that the resulting density remains numerically well-behaved and differentiable.

The mean molecular weight of the mixture is the mole-fraction-weighted average of the individual gas molecular weights $\mu_j$:
```math
\mu = \sum_j x_j \, \mu_j
```

## Heat capacity

The molar heat capacity at constant pressure $c_p$ determines how much energy is required to raise the temperature of the gas. Its value increases with temperature as rotational and vibrational degrees of freedom in polyatomic molecules become accessible[pierrehumbert_book_2010](@citep). AGNI therefore implements temperature-dependent $c_p(T)$ using tabulated experimental data from the JANAF thermochemical tables [JANAF](@citep), accessed by interpolation.

The heat capacity of the gas mixture is the mass-mixing-ratio-weighted sum of the per-species values:
```math
c_{pm} = \sum_j q_j \, c_{p,j}(T)
```

## Thermal conductivity

Molecular thermal conduction follows Fourier's law. It is negligible in the lower atmosphere but becomes the main way to move heat where the gas is optically thin.

The conductivity of the mixture is combined from mole fractions $x_j$ with the Wassiljewa rule and the Mason-Saxena coefficients [poling_properties_2001](@citep),
```math
\kappa_\text{mix} = \sum_i \frac{x_i \kappa_i}{\sum_j x_j \Phi_{ij}}, \qquad
\Phi_{ij} = \frac{\left[1 + (\eta_i/\eta_j)^{1/2} (M_j/M_i)^{1/4}\right]^2}{\left[8 (1 + M_i/M_j)\right]^{1/2}} .
```

The conductivity of each gas uses the Eucken relation [poling_properties_2001](@citep),
```math
\kappa_j = \eta_j \left( c_{v,j} + \frac{9}{4} \frac{R}{M_j} \right),
```
with $c_{v,j} = c_{p,j} - R/M_j$ per unit mass. For monatomic gases this reduces to $\kappa_j = \tfrac{15}{4} (R/M_j) \, \eta_j$, which AGNI uses directly.

The viscosity of each gas $j$ is taken from first-order Chapman-Enskog theory for the Lennard-Jones (12-6) potential,
```math
\eta_j = \frac{5}{16} \frac{\sqrt{\pi m_j k_B T}}{\pi \sigma_j^2 \, \Omega^{(2,2)*}(T/\varepsilon_j)},
```
where $m_j$ is the molecular mass, $\sigma_j$ the collision diameter, $\varepsilon_j$ the well depth (in K), and $\Omega^{(2,2)*}$ the reduced collision integral, evaluated with the empirical fit of [neufeld_collision_1972](@citet). The Lennard-Jones parameters are taken from the [GRI-Mech 3.0 transport file](http://combustion.berkeley.edu/gri-mech/version30/files30/transport.dat).

Gases missing from GRI-Mech are treated as a hard sphere of diameter 2 Å with $\Omega^{(2,2)*} = 1$.  No laboratory reference exists for atomic hydrogen, so its conductivity scales close to $T^{0.7}$, the exponent commonly used for thermospheres [watson_dynamics_1981](@citep). Dipole moments are ignored, so polar gases are less accurately treated.

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

Antoine coefficients are less readily available for refractory condensate species (e.g. Fe, FeO, MgO, Ti, TiO, TiO$_2$, VO). Their saturation curves are instead compiled from the fits tabulated by [woitke_chemistry_2018](@citet), which collates saturation-vapour-pressure relations for condensates relevant to rock-vapour and dust-forming atmospheres from a range of sources (their Table D2).

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
