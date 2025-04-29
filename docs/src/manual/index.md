# Development manual

## Contributing
If you are interested in contributing to the model, please contact the developers using the information on the main page.

## Coding style
- Indentation uses 4 spaces, no tabs.
- Function names should be lowercase, with words separated by underscores .
- Lines should aim to have a length of no more than 92 characters.
- All functions should have docstrings, ideally with Arguments and Returns listed.
- More comments are always better, even if they seem redundant.
- Use type hinting where possible.
- Print statements should be made through the logger where possible.
- The core package code should not contain global variables, except in the phys module.

## Code reference

### Atmosphere initialisation and variables
Functions from `atmosphere.jl`.
```@autodocs
Modules = [AGNI.atmosphere]
```

### Energy flux evaluation
Functions from `energy.jl`.
```@autodocs
Modules = [AGNI.energy]
```

### Numerical solver
Functions from `solver.jl`.
```@autodocs
Modules = [AGNI.solver]
```

### Plotting functions and utilities
Functions from `plotting.jl`.
```@autodocs
Modules = [AGNI.plotting]
```

### File I/O modules
Functions from `save.jl` and `load.jl`.
```@autodocs
Modules = [AGNI.save, AGNI.load]
```
