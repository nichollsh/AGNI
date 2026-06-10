# Composition and chemistry

## Condensation, evaporation, and oceans

AGNI incorporates a condensation and ocean model that handles the phase changes of volatile gases throughout the atmospheric column. Two reservoirs are tracked for each condensable species:

* An **initial** reservoir of surface condensate provided by the user (e.g. a prescribed ocean depth).
* A **total** reservoir updated at each solver step, based on the calculated chemistry and temperature profile.

At each solver iteration the following steps are performed in order:

1. The surface partial pressures of super-saturated volatiles are reduced to their saturation value $p^\text{sat}(T_s)$. The excess mass is added to the surface ocean reservoir. Sub-saturated volatiles are topped up from the ocean reservoir if condensate is available. The total surface pressure and gas mixing ratios are renormalised, and the pressure grid is regenerated.

2. Equilibrium chemistry is calculated across the column (see [Equilibrium chemistry](@ref) below).

3. At each pressure level in the column, super-saturated volatiles are condensed until their partial pressures equal $p^\text{sat}(T_l)$ at the local temperature $T_l$. Cold-trapping is enforced. Condensate produced at high altitudes is re-evaporated into drier layers below in a downward integration, until all rain either re-evaporates in the atmosphere or reaches the surface, where it contributes to the ocean reservoir.


## Oceans

The ocean layer structure is calculated from the total condensate mass assuming a liquid-phase density. The ocean basin area and continental shelf height determine the filling fraction of the basins. These can result in desert, continental, or aquaplanet scenarios.

## Phase change in the atmosphere

When a volatile condenses or evaporates, it releases or absorbs **latent heat** that contributes to the local heat budget. AGNI handles this through a diffusive condensation scheme with a fixed condensation timescale $t_\text{cond}$ [nicholls_convective_2025](@citep). The scheme is numerically differentiable with respect to temperature, which is necessary for numerical optimisation.

The latent heat transport flux at each cell edge is obtained by integrating the local heating rate downward from the top of the atmosphere:
```math
F_\text{lat}(p) = \int_0^p L_j(T) \frac{d \delta m_j}{dt} \, dp'
```
where $L_j(T)$ is the temperature-dependent latent heat of species $j$, $\delta m_j$ [$\mathrm{kg\ m^{-2}}$] is the mass of condensate produced per layer, and the condensation timescale $t_\text{cond}$ sets the effective rate of phase change [nicholls_convective_2025](@citep).

In condensing regions $F_\text{lat}$ is positive (heat is released); in evaporating regions it is negative. The energy balance is closed by requiring the total condensate mass + gaseous mass budget to be conserved.

The condensation timescale $t_\text{cond}$ represents the microphysics of cloud formation. AGNI uses a representative fixed value. Precipitation in deep atmospheres typically re-evaporates before reaching the surface, analogous to *virga* clouds on Earth.

Latent heats and thermodynamic properties are described in [Thermodynamics and EOS](@ref).

## Equilibrium chemistry

By default, AGNI treats the atmosphere as **well-mixed**: gas mixing ratios are constant with height. Condensation of a super-saturated volatile reduces its mixing ratio locally until it reaches saturation.

AGNI can alternatively couple to [FastChem](https://newstrangeworlds.github.io/FastChem/). FC is a fast numerical code for gas-phase thermochemical equilibrium [kitzmann_fastchem_2024, stock_faschem_2022, kitzmann_fastchem_2026](@citep).

FC solves for the gas-phase speciation of a mixture by minimising the Gibbs free energy subject to conservation of elemental mass and total electric charge. At chemical equilibrium, the number densities $n_j$ of species $j$ satisfy the **law of mass action**:
```math
n_j = K_j(T) \prod_e n_e^{\nu_{j,e}}
```
where $K_j(T)$ is the temperature-dependent equilibrium constant and $\nu_{j,e}$ is the stoichiometric coefficient of element $e$ in species $j$.

When `physics.chemistry = true` is set in the AGNI configuration, FastChem is called at each optimiser iteration. AGNI takes the column's elemental metallicity ratios (e.g. C/H, S/H), layer pressures, and layer temperatures as inputs for FC. The chemical calculation then returns volume mixing ratios at each layer.

Thermochemical equilibrium is a reasonable assumption when reaction timescales are shorter than the dynamical (convective) mixing timescale.

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
