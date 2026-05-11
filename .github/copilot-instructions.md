# AGNI — Copilot Instructions

AGNI is a Julia package that simulates radiative-convective equilibrium (RCE) in extreme rocky-exoplanet atmospheres. It wraps the Fortran SOCRATES radiative transfer code via a Julia layer.

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
julia --project=. test/runtests.jl
```

### Run fast tests only (skips integration tests, etc.)
```bash
julia --project=. test/runtests.jl fast
```

### Run a single test file
```bash
julia --project=. test/runtests.jl phys
```

### Generate coverage
```bash
julia --project=. --code-coverage test/runtests.jl
julia --project=. test/get_coverage.jl
```

## Architecture

The package entry point is `src/AGNI.jl`, which `include()`s submodule files **in a fixed order** (order matters for dependencies) and then imports them. The executable `agni.jl` activates the project, imports `AGNI`, and calls `AGNI.main()`.

The source tree is split hierarchically by responsibility:
- `src/phys/` for constants, formulae, density/EOS, and analytic profiles
- `src/state/` for model state structs and layer/diagnostic utilities
- `src/compose/` for chemistry/ocean composition logic
- `src/energy/` for RT, fluxes, spectrum, and optional RFM handling
- `src/interface/` for I/O, plotting, config-path helpers, and T(p) setup verbs
- `src/solver/` for solver front-end plus method-specific implementations
- `src/util/` for shared helpers (style, checksums)

### Central data structure
`atmosphere.Atmos_t` (defined in `src/state/atmosphere.jl`) is a mutable struct that holds all model state: pressure/temperature grids, flux arrays, gas mixing ratios, solver flags, file paths, and physics options. Nearly every function takes an `atmos::Atmos_t` as its first argument and mutates it.

### Typical execution flow (`AGNI.run_from_config`)
1. Parse and validate the TOML config dict
2. Instantiate `Atmos_t` → `atmosphere.setup!` → `atmosphere.allocate!`
3. Optionally set deep heating profile (`atmosphere.set_deep_heating!`)
4. Initialize T(p) profile via `setpt.request!`
5. Run `chemistry.calc_composition!`
6. Call the requested solver (`solver.solve_energy!` or `solver.solve_transparent!`)
7. Write outputs: `save.write_ncdf`, `save.write_profile`, then plotting functions

### Submodule responsibilities
| Module path | Role |
|---|---|
| `state/atmosphere.jl` | `Atmos_t` struct, setup/allocate/deallocate, hydrostatic integrator |
| `state/layers.jl`, `state/diagnostics.jl`, `state/multicol.jl` | Layer properties, diagnostics, and globe/multicol state |
| `phys/consts.jl`, `phys/phys.jl`, `phys/formulae.jl`, `phys/density.jl`, `phys/species.jl`, `phys/guillot.jl` | Constants, EOS/thermo utilities, species data, analytic profile pieces |
| `compose/chemistry.jl`, `compose/fastchem.jl`, `compose/ocean.jl` | Thermochemical equilibrium (FastChem), rainout/condensation, surface oceans |
| `energy/energy.jl`, `energy/spectrum.jl`, `energy/rfm.jl` | SOCRATES/RFM RT interfaces, fluxes/heating, spectrum utilities |
| `interface/setpt.jl`, `interface/save.jl`, `interface/load.jl`, `interface/plotting.jl`, `interface/paths.jl` | T(p) setup verbs, I/O, plotting, and path safety helpers |
| `solver/solver.jl` + `solver/*.jl` | Solver front-end and method-specific implementations (energy, transparent, globe, etc.) |
| `util/blake.jl`, `util/style.jl` | Shared utility helpers (hash checks, style/label/color helpers) |

### External dependencies
- **SOCRATES** (Fortran RT code, loaded at runtime via `include(ENV["RAD_DIR"]/julia/src/SOCRATES.jl)`). The `atmosphere` module hard-includes this path on load.
- **FastChem** (optional, for thermochemical equilibrium). Fetched by `src/get_fastchem.sh`.
- **RFM** (optional, line-by-line RT). Used only when `files.rfm_parfile` is set in config.

## Key Conventions

### Module files must not be run directly
The `src/AGNI.jl` file begins with:
```julia
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end
```

### `include()` order in `src/AGNI.jl` is load-order dependent
New source files must be added at the correct position in `AGNI.jl`'s `include` list (dependencies before dependents) and then explicitly `import`ed and `export`ed.

### Bang (`!`) functions mutate `Atmos_t`
Functions like `setup!`, `allocate!`, `deallocate!`, `set_deep_heating!` follow the Julia convention of mutating their primary argument in place.

### Configuration is TOML-driven
All physics options, file paths, solver settings, and output flags come from a single TOML config file. The config dict is validated in `AGNI.open_config` and `AGNI.run_from_config` before any model objects are created. Add new config keys there first.

### Logging
Use Julia's `@info`, `@warn`, `@debug`, `@error` macros (not `println`). Verbosity is controlled by `execution.verbosity` in the config (0 = silent, 1 = normal, 2 = debug).

### Pressure and unit conventions
- Pressures in TOML configs are in **bar**; internally AGNI works in **Pa** (×1e5 conversion).
- Temperatures in K, fluxes in W m⁻², lengths in m throughout the Julia code.

### Test suite structure
Unit tests are individual files (`test_consts.jl`, `test_phys.jl`, etc.) included by `runtests.jl`. Slow tests, such as those requiring SOCRATES data, are skipped with the `fast` argument. Add new unit test files by `include()`ing them in `runtests.jl`.
