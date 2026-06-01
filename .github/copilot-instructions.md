# AGNI — AI Agent Instructions

**Trust these instructions.** Only search if information is incomplete or found to be in error.

**Identity & Mission**: You are an expert Scientific Software Engineer working on the AGNI module of the PROTEUS ecosystem.

## High-Level Instructions

### Rule files you MUST read on every session

AGNI keeps its AI Agent rule files under `.github/copilot-instructions.md`, which is this file.

This file includes a test quality deep-dive: anti-happy-path patterns, discriminating-value guards, physics-invariant tiering, validation certification markers, adversarial-review trigger, buffer-flip propagation, hypothesis seed stability, solver intermediate-state assertions. **Required reading before editing any file under `tests/**` or `src/**`.**. It details the domain-aware physics review (buffer-flip safety, solver intermediate-state types, four PROTEUS-coupling patterns). **Required reading before any code review pass.**

1. **Always** read the two rule files above plus the Testing Standards section below before any code change.

2. **Always** inform the user that you are reading in this file by printing a message at the start of your response: "(Read in copilot-instructions.md...)"

3. When creating a PR, **always** follow the PR template (`.github/pull_request_template.md`) and ensure all sections are filled out with relevant information.


## Model and Ecosystem Context

AGNI is a Julia package that simulates radiative-convective equilibrium (RCE) in extreme rocky-exoplanet atmospheres. It wraps the Fortran SOCRATES radiative transfer code via a Julia layer. AGNI can function as a standalone 1D RCE model, and is accessible via a Command Line Interface, the Julia REPL, Jupyter Notebooks, and via Python through PyCall/JuliaCall.

AGNI is the main atmosphere climate model of the PROTEUS simulation framework software ecosystem. It is called by the main [PROTEUS](https://github.com/FormingWorlds/PROTEUS) coupled atmosphere-interior framework during the climate calculation step.

Sister modules in the ecosystem: CALLIOPE (atmospheric outgassing), Atmodeller (outgassing), SOCRATES (spectral radiative transfer), JANUS (1D convective atmosphere), MORS (stellar evolution), Aragog / SPIDER (interior thermal evolution), VULCAN (atmospheric disequilibrium chemistry), ZEPHYRUS REAS (atmospheric escape), Obliqua (tidal evolution), FastChem (atmospheric thermochemistry).

**Project Type**: Scientific simulation code.

**Languages**: Julia 1.12+ (main), Modern FORTRAN (F90, gfortran compiled).

**Target Platform**: Linux AMD64 (x86-64) or MacOS (Apple Silicon M1/M2). Windows is not currently supported.

## Build, Test, and Run

### Prerequisites
`RAD_DIR` must be set to a compiled SOCRATES installation before any Julia invocation:
```bash
export RAD_DIR=/path/to/SOCRATES
```

### Build (install Julia packages + compile SOCRATES Julia wrappers)
```bash
julia --project -e 'using Pkg; Pkg.build()'
```

### Run the model
```bash
julia agni.jl                          # uses res/config/default.toml
julia agni.jl res/config/hotdry.toml  # use a specific config
```

### Run all tests
```bash
julia --project test/runtests.jl
```

### Run fast tests only (skips integration tests, etc.)
```bash
julia --project test/runtests.jl fast
```

### Run a single test file
```bash
julia --project test/runtests.jl phys
```

### Generate coverage
```bash
julia --project --code-coverage test/runtests.jl
julia --project test/runcoverage.jl
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

Module files must not be run directly. The executable entry point is `agni.jl`, which activates the project and calls `AGNI.main()`. This ensures a single load order and consistent environment for all code paths.

Changes to function signatures should avoid breaking high-level API calls, since AGNI is part of a wider software framework.

### `include()` order in `src/AGNI.jl` is load-order dependent
New source files must be added at the correct position in `AGNI.jl`'s `include` list (dependencies before dependents) and then explicitly `import`ed and `export`ed. All includes are at the top of `AGNI.jl` before any code, and other modules will import from `AGNI` rather than each other. Thie ensures a single load order and avoids circular dependencies.

### Bang (`!`) functions mutate `Atmos_t`
Functions like `setup!`, `allocate!`, `deallocate!`, `set_deep_heating!` follow the Julia convention of mutating their primary argument in place.

### Configuration is TOML-driven
All physics options, file paths, solver settings, and output flags come from a single TOML config file. The config dict is validated in `AGNI.open_config` and `AGNI.run_from_config` before any model objects are created. Add new config keys there first.

### Logging
Use Julia's `@info`, `@warn`, `@debug`, `@error` macros (not `println`). Verbosity is controlled by `execution.verbosity` in the config (0 = silent, 1 = normal, 2 = debug).

### Pressure and unit conventions
- Pressures in TOML configs are in **bar**; internally AGNI works in **Pa** (×1e5 conversion).
- Temperatures in K, fluxes in W m⁻², lengths in m throughout the Julia code.

## Documentation

The documentation lives in the `docs/` folder and is built with Documenter.jl. The main entry point is `docs/src/index.md`, which links to other pages. The documentation includes a user guide, developer guide, and API reference.

When adding new functions or modules, include docstrings in the source code. These will be automatically included in the API reference when the docs are built.

The documentation should remain 1:1 consistent with the codebase. If you change a function signature, update the docstring and the docs page if it is mentioned there. If you add a new config option, document it in the config reference page.

The documentation should also describe the physical assumptions and limitations of the model, as well as any known issues or caveats. This is crucial for users to understand the applicability of AGNI to their research questions. It should outline the physics parameterizations, solver methods, and typical use cases. It should describe explicitly the equations solved, which can be typset in LaTeX and rendered in the docs. Literature references should be included for the physical models and numerical methods used -- this includes URLs which point to scientific papers, textbooks, or online resources and websites that explain the underlying physics and algorithms.

Documentation can include information about SOCRATES, FastChem, and RFM dependencies, but should not require users to understand the internal workings of those codes. The focus should be on how to use AGNI and interpret its outputs, rather than the implementation details of the underlying physics.

Documentation should be clear, concise, and accessible to the target audience of planetary scientists and climate modellers. It can use jargon but should explain any technical terms that are necessary.

## Test suite

Unit tests are individual files (`test_consts.jl`, `test_phys.jl`, etc.) included by `runtests.jl`. Slow tests, such as those requiring SOCRATES data, are skipped with the `fast` argument. Add new unit test files by `include()`ing them in `runtests.jl`.

### Testing Standards

AGNI is scientific simulation code, so the test suite is held to physics-grade rigor. The rules below are the contract.


### Anti-happy-path rules (every new test)

Every new test function MUST include:

1. **At least one edge case** (boundary value, empty input, extreme physical parameter).
2. **At least one path that exercises the error contract** (documented exception, guard return, graceful clamp). If the function under test has no validation, exercise the limit-input behavior and assert the mathematical invariant.
3. **Assertion values that are NOT trivially derivable from the implementation**: discriminating numeric pins or property-based assertions (monotonicity, conservation) preferred over point checks.

**Forbidden patterns itemized below** :

- Single-assert test functions.
- Standalone weak assertions as the only meaningful check (e.g. assert result is not None).
- Tests with no function-level comment (comments should be a clear statement of the physical scenario or contract clause being verified, and formatted with a `#` symbol).
- Tests using `==` adjacent to float literals.
- Tests asserting on a fixture's implicit default.

### Float and numerical comparison

NEVER use `==` for floats. Use functions which check for equality within a tolerance. For pinned numeric values, include a **discrimination guard**: a follow-up `assert` showing the wrong-formula / wrong-buffer / wrong-law value would differ from the correct one by more than the tolerance.

### Documentation per test

- File-level comment: name the source under test, list the invariants and contract clauses the file exercises (formatted with `#` rather than with quotes).
- Function-level comment: state the physical scenario or contract clause being verified. Required (lint-enforced).
- Inline comments: explain **why** a specific input range was chosen.

**Remember**: Trust these instructions. Only search if information is incomplete or found to be in error.

###  Naming tests

- Test names describe behavior, not the called function: `test_olr_increases_with_temperature`, NOT `test_olr`.
- Test names use snake_case and read as full sentences.
- Group related tests in testsuites when they share setup
- Avoid redundant or additional calls to `atmosphere.allocate!` in fast tests, if possible.
- Test file names should aim to mirror source within the `src/**` tree.


### Physics plausibility

- Temperature must be positive everywhere (Kelvin). Flag any code path where T could reach zero or go negative.
- Total pressure must be positive; partial pressures must be non-negative and sum to total pressure within solver tolerance.
- Mole fractions must sum to 1.0. Flag any composition-returning function that doesn't enforce or verify normalization.
- Radiative fluxes must be positive in magnitude. Flag any path where fluxes could flip sign or become NaN.
- Assert that the model's T(p) profiles are physically plausible (e.g. no negative temperatures, no super-adiabatic gradients without convection).

