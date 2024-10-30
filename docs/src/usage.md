# Running the model
First, follow the [Getting started](@ref) instructions. Only read-on once you
have confirmed that the code is working.

## Input data files
The minimal input data required to run the model will have been downloaded automatically.
If you require more data, such as additional stellar spectra or opacities, then these can
also be easily obtained using the `get_data` script in the AGNI root directory. To see how
to use this script, run it without arguments like so:
```bash
./src/get_data.sh
```

## Tutorials
There are Jupyter notebooks containing tutorials in the `tutorials/` directory
of the repository.

## General execution
The environment variable `RAD_DIR` must point to the SOCRATES installation
directory. This is required for AGNI to find the SOCRATES libraries. The best way to do
this is to add `RAD_DIR=path/to/socrates/folder/` to your shell rc file (e.g. `~/.bashrc`).

Then to use the model, simply run `./agni.jl [cfg]` where `[cfg]` is the path
to the required configuration file. If `[cfg]` is not passed, then the default
configuration will be used.

To calculate equilibrium chemistry self-consistently with the climate, FastChem
is coupled to AGNI. For AGNI to find FastChem, the environment variable `FC_DIR`
must point to the FastChem installation directory.

## Configuration
AGNI configuration files are formatted using [TOML](https://toml.io/en/). There
are examples in `res/config/`. The default configuration file contains comments
explaining the purpose of each parameter, although some are explained in greater
detail below. Take care to format the variables in the TOML file correctly.
There are no 'default values'. Not all parameters are required in all cases,
but the model will return an error naming any parameters which are both
necessary and absent.

Broadly, the configuration files are broken up into four sections:
* `[planet]` -  general characteristics of the planet
* `[files]` - input/output files and other paths
* `[composition]` - atmospheric composition and chemistry
* `[execution]` - what the model should do
* `[plots]` - which plots should be produced

Some parameters:
* `files.input_sf` is the file path to the "spectral file" containing opacity data. Several spectral files are packged with AGNI, but you can find more online via the [Open Science Framework](https://osf.io/vehxg/).

* `execution.solution_type` tells the model which state to solve for. The allowed values (integers) are...
     - 1 : zero flux divergence at fixed `tmp_surf`
     - 2 : zero flux divergence, with `tmp_surf` set such that the conductive skin (CBL) conserves energy flux
     - 3 : the net flux (up minus down) at each layer is equal to `flux_int`

* `execution.solvers` tells the model which solvers to use. This is a list of strings, so multiple solvers can be applied sequentially. An empty string is always appended to the end of this list. Allowed solvers are...
     - [empty string] : no solving takes place, so the model just calculates fluxes using the initial state
     - `newton` : the Newton-Raphson algorithm is used
     - `gauss`  : the Gauss-Newton algorithm is used
     - `levenberg` : the Levenberg–Marquardt algorithm is used

* `execution.initial_state` describes the initial temperature profile applied to the atmosphere. This is a list of strings which are applied in the given order, which allows the user to describe a specific state as required. The descriptors are listed below, some of which take a single argument that needs to immediately follow the descriptor in the list order.
     - `dry`              : integrate the dry adiabatic lapse rate from the surface upwards
     - `str`,       `arg` : apply an isothermal stratosphere at `arg` kelvin
     - `iso`,       `arg` : set the whole atmosphere to be isothermal at `arg` kelvin
     - `csv`,       `arg` : set the temperature profile using the CSV file at the file path `arg`
     - `sat`,       `arg` : apply Clausius-Clapeyron saturation curve for the gas `arg`
     - `ncdf`,      `arg` : load profile from the NetCDF file located at `arg`
     - `loglin`,    `arg` : log-linear profile between `tmp_surf` at the bottom and `arg` at the top

    For example, setting `initial_state = ["dry", "sat", "H2O", "str", "180"]` will set T(p) to follow the dry adiabat from the surface, the water condensation curve above that, and then to be isothermal at 180 K until the top of the model.

* `composition.chem_type` describes the type of chemistry to implement within the model. This is handled externally by FastChem, so you must set the environment variable `FC_DIR` to point to the FastChem directory. The allowed values (integers) are...
     - 0 : Disabled
     - 1 : Equilibrium, gas phase only
     - 2 : Equilibrium, with condensation (condensates retained)
     - 3 : Equilibrium, with condensation (condensates rained out)

## Outputs
Results are optionally plotted and animated, and data will be saved as NetCDF
or CSV files.

## Python
It is possible to interact with the model using Python. This is best done with the `juliacall` package from [PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl), and is implemented this way in the [PROTEUS framework](https://github.com/FormingWorlds/PROTEUS).

