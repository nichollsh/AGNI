# Thermodynamics and EOS

## Gas densities

The density of the gas mixture is calculated using Amagat's additive volume law to combine the densities of the components. You can read more about the validity and usage of this formulation here: [magyar_eos_2014, bradley_experimental_2018, magyar_mixing_2013, magyar_ethane_2015](@citet).

* The densities of each gas component are nominally calculated using the Van der Waals equation of state [Kontogeorgis_TakingAnot_2019](@citep).
* AQUA is implemented as the EOS for water [haldemann_aqua_2020](@citep).
* The [chabrier_eos_2019](@citet) EOS is implemented as the EOS for hydrogen .
* AGNI will fallback to the ideal gas EOS for otherwise unsupported gases.

## Latent heat

Latent heats are temperature-dependent, using values derived from [coker_thermo_2007](@citet) and [IAPWS95](@citet). Heat capacities are also temperature-dependent, using values derived from the JANAF database. See the [ThermoTools repo](https://github.com/nichollsh/ThermoTools) for scripts.

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
