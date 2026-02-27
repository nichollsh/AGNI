# AGNI — Copilot Instructions

AGNI is a Julia package that simulates radiative-convective equilibrium (RCE) in extreme rocky-exoplanet atmospheres. It wraps the Fortran SOCRATES radiative transfer code via a Julia FFI layer.

## Build, Test, and Run

### Prerequisites
`RAD_DIR` must be set to a compiled SOCRATES installation before any Julia invocation:
```bash
export RAD_DIR=/path/to/SOCRATES
```

### Build (install Julia packages + compile SOCRATES Julia wrappers)
```bash
julia --project=. -e 'using Pkg; Pkg.build()'
```

### Run the model
```bash
julia agni.jl                          # uses res/config/default.toml
julia agni.jl res/config/hotdry.toml  # use a specific config
```

### Run all tests
```bash
cd test/
julia --project=.. runtests.jl
```

### Run fast tests only (skips integration tests)
```bash
cd test/
julia --project=.. runtests.jl fast
```

### Run a single test file
```bash
cd test/
julia --project=.. -e 'include("test_phys.jl")'
```

### Generate coverage
```bash
cd test/
julia --project=.. --code-coverage runtests.jl
julia --project=.. ../test/get_coverage.jl
```

## Architecture

The package entry point is `src/AGNI.jl`, which `include()`s submodule files **in a fixed order** (order matters for dependencies) and then imports them. The executable `agni.jl` activates the project, imports `AGNI`, and calls `AGNI.main()`.

### Central data structure
`atmosphere.Atmos_t` (defined in `src/atmosphere.jl`) is a mutable struct that holds all model state: pressure/temperature grids, flux arrays, gas mixing ratios, solver flags, file paths, and physics options. Nearly every function takes an `atmos::Atmos_t` as its first argument and mutates it.

### Typical execution flow (`AGNI.run_from_config`)
1. Parse and validate the TOML config dict
2. Instantiate `Atmos_t` → `atmosphere.setup!` → `atmosphere.allocate!`
3. Optionally set deep heating profile (`atmosphere.set_deep_heating!`)
4. Initialize T(p) profile via `setpt.request!`
5. Run `chemistry.calc_composition!`
6. Call the requested solver (`solver.solve_energy!` or `solver.solve_transparent!`)
7. Write outputs: `save.write_ncdf`, `save.write_profile`, then plotting functions

### Submodule responsibilities
| Module | Role |
|---|---|
| `atmosphere` | `Atmos_t` struct, setup/allocate/deallocate, hydrostatic integrator |
| `energy` | Flux calculations: SOCRATES RT, MLT convection, conduction, latent/sensible heat |
| `solver` | Nonlinear solvers (Newton, Gauss-Newton, Levenberg-Marquardt, Jacobi-Newton) |
| `setpt` | Setting initial T(p) profiles (dry adiabat, saturation, Guillot, custom) |
| `spectrum` | Spectral file loading and management |
| `chemistry` | Thermochemical equilibrium via FastChem; rainout logic |
| `ocean` | Surface liquid ocean formation and distribution |
| `rfm` | Optional line-by-line RT via the RFM code |
| `phys` | Physical utility functions (EOS, thermodynamics) |
| `consts` | Physical constants |
| `save` | NetCDF and CSV output |
| `plotting` | Makie/GR-based diagnostic plots |
| `load` | Reading profiles and external data |
| `guillot` / `blake` | Analytic T(p) profile implementations |

### External dependencies
- **SOCRATES** (Fortran RT code, loaded at runtime via `include(ENV["RAD_DIR"]/julia/src/SOCRATES.jl)`). The `atmosphere` module hard-includes this path on load.
- **FastChem** (optional, for thermochemical equilibrium). Fetched by `src/get_fastchem.sh`.
- **RFM** (optional, line-by-line RT). Used only when `files.rfm_parfile` is set in config.

## Key Conventions

### Module files must not be run directly
Every `src/*.jl` file begins with:
```julia
if (abspath(PROGRAM_FILE) == @__FILE__)
    error("The file '$thisfile' is not for direct execution")
end
```

### `include()` order in `src/AGNI.jl` is load-order dependent
New source files must be added at the correct position in `AGNI.jl`'s `include` list (dependencies before dependents) and then explicitly `import`ed and `export`ed.

### Bang (`!`) functions mutate `Atmos_t`
Functions like `setup!`, `allocate!`, `deallocate!`, `set_deep_heating!` follow the Julia convention of mutating their primary argument in place and returning `Bool` (true = success).

### Configuration is TOML-driven
All physics options, file paths, solver settings, and output flags come from a single TOML config file. The config dict is validated in `AGNI.open_config` and `AGNI.run_from_config` before any model objects are created. Add new config keys there first.

### Logging
Use Julia's `@info`, `@warn`, `@debug`, `@error` macros (not `println`). Verbosity is controlled by `execution.verbosity` in the config (0 = silent, 1 = normal, 2 = debug).

### Pressure and unit conventions
- Pressures in TOML configs are in **bar**; internally AGNI works in **Pa** (×1e5 conversion).
- Temperatures in K, fluxes in W m⁻², lengths in m throughout the Julia code.

### Test suite structure
Unit tests are individual files (`test_consts.jl`, `test_phys.jl`, etc.) included by `runtests.jl`. Integration tests (`test_integration.jl`) require SOCRATES data and are skipped with the `fast` argument. Add new unit test files by `include()`ing them in `runtests.jl`.
