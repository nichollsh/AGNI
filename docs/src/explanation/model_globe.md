# Multicolumn and global mode

## Motivation

Most codes, AGNI included, approach atmospheres as a single 1D column. They thereby represent the full planet with a single zenith angle and stellar flux scale factor (see [Stellar irradiation](@ref)). This naturally cannot capture how climate states vary with longitude and latitude; e.g. the day-night temperature contrast expected on a tidally locked rocky exoplanet. AGNI's multicolumn (globe) mode addresses this by solving several 1D columns, each placed at a different point on the planet's surface, and coupling them through a horizontal heat-redistribution flux. This retains the numerical efficiency of the 1D column solver, without the cost of a fully 3D general circulation model.

## Constructing the globe

A globe is built from a single, already-configured `Atmos_t` instance, which acts as the template ("worker") atmosphere. Multicolumn mode is enabled in the configuration file by adding a `[planet.globe]` section that specifies the longitude and latitude of each column:

```math
\{(\lambda_i, \phi_i)\}_{i=1}^{N_\text{col}}
```

where $\lambda_i \in [0^\circ, 360^\circ]$ and $\phi_i \in [-90^\circ, 90^\circ]$. AGNI constructs one independent copy of the atmosphere for each column, following the pattern described in [`multicol.construct!`](@ref). Each column's zenith angle is derived from its longitude and latitude (`atmosphere.calc_zenith_angle`), and its stellar flux is set accordingly, so that columns nearer the substellar point receive more instellation than those nearer the terminator or nightside. The underlying [SOCRATES](https://proteus-framework.org/SOCRATES/) Fortran library is not thread-safe, so columns are solved sequentially.

## Heat redistribution

Each column receives a heat-redistribution flux $F^\text{redist}_i$, deposited in log-pressure: an [advective heating term](@ref "Advective heating"). This allows a column to act as either a net heat source or a net heat sink [cronin_advective_2016](@citep). Physically-motivated scaling laws can drive the redistribution strength and profile quantities self-consistently from each column's temperature, composition, planetary radius, and rotation rate.

## Solving the globe

`solver.solve_globe!` iterates over the set of columns, calling the standard [Newton-Raphson column solver](@ref "Obtaining a solution") independently on each one, then checks whether all columns have converged to a consistent total flux. Depending on `solution_type` (see [Solution types](@ref)). Each outer iteration re-solves every column to full convergence, multicolumn runs are correspondingly more computationally expensive than single-column runs, scaling roughly with the number of columns.

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
