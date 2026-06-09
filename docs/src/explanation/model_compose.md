# Composition and chemistry

## Condensation, evaporation, and oceans

AGNI incorporates a simple ocean model, which is tied to the atmosphere rainout and evaporation schemes. This is divided into two reservoirs for each condensable component:
* An initial reservoir of surface condensate, provided by the user.
* A total reservoir of surface condensate, based on the calculated chemistry and temperature profile.

The condensation and chemistry calculations then operate together. At each solver step, the following actions occur:
1. The *surface* temperature and partial pressures are checked against saturation, for each condensable
    - Super-saturated condensables have their partial pressures decreased, and the mass is added to the ocean reservoir
    - Sub-saturated condensables have their partial pressures increased based on the availability of surface condensate
2. Chemistry is then performed to calculate the gas phase speciation; see [Equilibrium chemistry](@ref) below
3. Rainout is calculated aloft, based on sub/super-saturation at *every* atmosphere level.

In the first step, the model is effectively non-hydrostatic because the surface pressure (and whole pressure grid) is adjusted according to *surface* sub/super-saturation. The sum of condensation at the surface and aloft go towards modulating the surface ocean content. The ocean layer structure is calculated according to liquid density. Two parameters (ocean basin area, continental shelf height) determine the filling fraction of the basins; i.e. whether the planet is a 'desert planet', a 'continental planet', or an 'aqua planet'.

During the third step, where phase change happens aloft, the mixing ratios of dry species are increased in order to satisfy the total pressure at condensing levels. This is treated as a hydrostatic process. The total accumulated amount of condensation (for each volatile) is then re-evaporated in the deeper atmosphere where possible. Rain reaching the surface contributes to the ocean.

## Phase change in the atmosphere

Gases release energy "latent heat" into their surroundings when condensing into a liquid or solid. This is included in the model through a diffusive condensation scheme, which assumes a fixed condensation timescale. Any rain which is not re-evaporated before reaching the surface is considered to contribute towards forming an ocean (secondary reservoir).

The latent heating associated with the change in partial pressure of condensable gases in the atmosphere is used to calculate a latent heating rate at each level of the model (positive where condensing, negative where evaporating). The heating rates in each layer are then integrated (from the TOA downwards) to provide a latent heat transport *flux* at cell-edges, with the assumption being that condensation occurs by updrafts:
```math
F_{\text{latent}}(p) = \int_0^p L_v(T) \frac{dq}{dt} dp
```
where $L_v(T)$ is the temperature-dependent latent heat of vaporisation, $q$ is the mixing ratio of the condensable species, and the integral is performed from the top of atmosphere downwards. The integrated condensable heat flux is balanced by evaporation at deeper layers which closes the energy balance.

This method is conceptually similar to [derras_estimating_2025](@citet). AGNI assumes that no work is done during phase changes aloft, $p dV \approx 0$, so enthalpy and latent heat become equivalent.

Latent heats and other thermodynamic variables are described in [Thermodynamics and EOS](@ref).

## Equilibrium chemistry

By default, AGNI assumes that the atmosphere composition is "well-mixed". This means that the mixing ratios of the species are constant with height. Condensation of a super-saturated volatile will reduce its mixing ratio such that it becomes exactly saturated.

AGNI can couple to [FastChem](https://newstrangeworlds.github.io/FastChem/) — a fast numerical code for gas-phase thermochemical equilibrium [kitzmann_fastchem_2024, stock_faschem_2022](@citep). FastChem takes elemental abundances (metallicities), pressures, and temperatures as input variables and returns the mixing ratios of hundreds of gas-phase species at thermochemical equilibrium. When the configuration variable `physics.chemistry = true` FastChem will be enabled.

At each step of the solver loop, AGNI derives the atmosphere's bulk elemental metallicity from the gas composition (or from user-provided metallicities), then calls FastChem to compute gas-phase equilibrium mixing ratios across the column. AGNI's own rainout and condensation scheme (`_sat_aloft!`) subsequently operates on the post-FastChem composition, handling super-saturation and cold-trapping. The two schemes are therefore **complementary**: FastChem handles gas-phase speciation while AGNI handles condensation independently.

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
