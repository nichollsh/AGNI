# API reference

## Atmosphere initialisation and variables
Functions from `atmosphere.jl`.
```@autodocs
Modules = [AGNI.atmosphere]
```

## Energy flux calculations
Functions from `energy.jl`.
```@autodocs
Modules = [AGNI.energy]
```

## Aerosol optical properties
Functions from `mie.jl` and `aerosol_optics.jl`.
```@autodocs
Modules = [AGNI.mie, AGNI.aerosol_optics]
```

## Chemistry and phase-change calculations
Functions from `chemistry.jl`.
```@autodocs
Modules = [AGNI.chemistry]
```

## Numerical solvers
Functions from `solver.jl`, including the energy-conserving solver `solve_energy!` (whose
keyword arguments select the linesearch and Broyden updates) and the linesearch helpers.
```@autodocs
Modules = [AGNI.solver, AGNI.solver.solve_energy, AGNI.solver.linesearch]
```

## Multicolumn module
Functions from `multicol.jl`.
```@autodocs
Modules = [AGNI.multicol]
```

## Plotting functions and utilities
Functions from `plotting.jl`.
```@autodocs
Modules = [AGNI.plotting]
```

## Physical properties and constants
Functions from `phys.jl`.
```@autodocs
Modules = [AGNI.phys]
```

## File I/O modules
Functions from `save.jl` and `load.jl`.
```@autodocs
Modules = [AGNI.save, AGNI.load]
```
