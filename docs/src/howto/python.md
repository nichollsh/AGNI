# Accessing AGNI from Python

It is possible to interact with AGNI from Python using the `juliacall` package from
[PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl).

Coupling with Python is done via `juliacall` within the modular
[PROTEUS framework](https://github.com/FormingWorlds/PROTEUS), which couples AGNI
self-consistently to models of planetary interior evolution and volatile outgassing. You
can see this implemented in
[agni.py](https://github.com/FormingWorlds/PROTEUS/blob/main/src/proteus/atmos_clim/agni.py)
within the PROTEUS source code.

## Minimal example

```python
# Import juliacall
from juliacall import Main as jl

# Import AGNI
jl.seval("using Pkg")
jl.Pkg.activate(AGNI_ROOT_DIR)  # <---- set AGNI_ROOT_DIR to your installation path
jl.seval("import AGNI")
jl.AGNI.setup_logging("out.log", 1)

# Setup atmosphere
atmos = jl.AGNI.atmosphere.Atmos_t()
jl.AGNI.atmosphere.setup_b(atmos, ...)   # <--- complete function arguments as per docstring in `AGNI.atmosphere.setup!()`

# Allocate atmosphere
jl.AGNI.atmosphere.allocate_b(atmos, STAR_SPECTRUM_FILE)   # <-- provide path to spectrum

# Solve T(p)
jl.AGNI.solver.solve_energy_b(atmos)

# Write results to a file
jl.AGNI.save.write_ncdf(atmos, "out.nc")
```

See the [API reference](@ref) for full documentation of all available functions.
