# Using the model
Follow the [Getting started](@ref) instructions. Only read on once you have confirmed that the model runs on your machine.
Otherwise, see the [Troubleshooting](@ref) guide.

## Input data files
The minimal input data required to run the model will have been downloaded automatically from Zenodo. If you require more data, such as additional stellar spectra or opacities, then these can also be easily obtained using the `get_data` script in the AGNI root directory. To see how to use this script, run it without arguments like so:
```bash
./src/get_data.sh
```

Opacities are contained within "spectral files". Use the table within
`res/spectral_files/reference.pdf` to decide which spectral files are best for you.

For example, if you wanted to get the spectral file "Honeyside48" you would run:
```bash
./src/get_data.sh anyspec Honeyside 48
```

## Tutorials
There are Jupyter notebooks containing tutorials in the `tutorials/` directory
of the repository.

## General execution
The environment variable `RAD_DIR` must point to the SOCRATES installation
directory. This is required for AGNI to find the SOCRATES libraries. The best way to do
this is to add `RAD_DIR=path/to/socrates/folder/` to your shell rc file (e.g. `~/.bashrc`).

Then to use the model, simply run
```bash
./agni.jl [cfg]
```
where `[cfg]` is the path to the desired configuration file.
If `[cfg]` is not passed, then the default configuration file will be used.

You should see the following output:
```log
[ INFO  ] Using configuration 'Default'
[ INFO  ] Setting-up a new atmosphere struct
[ INFO  ] Loading thermodyamic data
[ INFO  ] Inserting stellar spectrum and Rayleigh coefficients
[ INFO  ] Allocating atmosphere with composition:
[ INFO  ]       1 H2O     1.00e+00 (EOS_AQUA)
[ INFO  ] Setting T(p): dry, sat
[ INFO  ] Solving with 'none'
[ INFO  ]     done
[ INFO  ] Total RT evalulations: 2
[ INFO  ] Writing results
[ INFO  ] Plotting results
[ INFO  ] Deallocating memory
[ INFO  ] Model runtime: 16.60 seconds
```

The line following "Allocating atmosphere with composition" is a table of gases, their
volume mixing ratios, and flags. In this case there is only one gas.

Potential flags for each species are:
* `EOS_[XX]` - using the `[XX]` equation of state (e.g. ideal gas, AQUA)
* `NO_OPACITY` - no opacity data available, but can contribute to the thermodynamics
* `NO_THERMO` - no thermodynamic data available, so will be treated as a diatomic ideal gas
* `COND` - this gas is allowed to condense

## Configuration
AGNI configuration files are formatted using [TOML](https://toml.io/en/). There
are examples in `res/config/`. The default configuration file contains comments
explaining the purpose of each parameter, although some are explained in greater
detail below. Take care to format the variables in the TOML file correctly.
There are no 'default values'. Not all parameters are required in all cases,
but the model will return an error naming any parameters which are both
necessary and absent.

The configuration files are broken up into four sections (or "tables"), each containing a
number of parameters. These tables of parameters are described below.

### `[planet]`
General properties of the planet.

| Parameter          | Description   |
| -----------------: | :------------ |
| `tmp_surf        ` | Temperature of the surface [kelvin].  |
| `instellation    ` | Stellar flux at planet's orbital distance [W m-2]. |
| `albedo_b        ` | A pseudo bond-albedo which downscales the stellar flux by 1-albedo_b. |
| `s0_fact         ` | Stellar flux scale factor which accounts for planetary rotation (c.f. Cronin+13). |
| `zenith_angle    ` | Characteristic zenith angle for incoming stellar radiation [degrees]. |
| `surface_material` | Surface material (can be "greybody" or path to data file in `res/`). |
| `albedo_s        ` | Spectrally-grey surface albedo used when `material=greybody`. Optional otherwise. |
| `radius          ` | Radius of the planet's surface [m]. |
| `gravity         ` | Gravitational acceleration at the surface [m s-2]. Incompatible with `mass` option. |
| `mass            ` | Total mass of material below the atmosphere [kg]. Incompatible with `gravity` option. |
| `skin_d          ` | Thickness of the conductive boundary layer [m]. Used when `solution_type=2`. |
| `skin_k          ` | Conductivity of the conductive boundary layer [W m-1 K-1]. Used when `solution_type=2`. |
| `tmp_magma       ` | Temperature of the topmost layer of the planet's mantle [K]. Used when `solution_type=2`. |
| `flux_int        ` | Internal flux [W m-2] to be solved-for when `solution_type=3`. |
| `turb_coeff      ` | Turbulent exchange coefficient for sensible heat transport. |
| `wind_speed      ` | Effective wind speed for sensible heat transport [m s-1]. |
| `star_Teff      `  | Stellar photospheric temperature [K] used if `files.input_star=="blackbody"`. |


### `[files]`
Input/output files and other paths.

| Parameter          | Description   |
| -----------------: | :------------ |
| `input_sf       `  | Path to the desired spectral file ending in `.sf`, in `res/spectral_files/`. |
| `input_star     `  | Path to stellar spectrum. If blank, spectrum assumed to be inside spectral file. If "blackbody" must provide `planet.star_Teff`. |
| `output_dir     `  | Path to the output directory. |
| `rfm_parfile  `    | Path to .par linelist file, for running line-by-line calculations with the RFM. |


### `[composition]`
Atmospheric composition and chemistry.

| Parameter         | Description   |
| ----------------: | :------------ |
| `p_top         `  | Total top-of-atmosphere pressure [bar]. |
| `p_dict        `  | Dictionary of gas partial surface pressures [bar]. Summed to obtain `p_surf`. |
| `p_surf        `  | Total surface pressure [bar]. Incompatible with `p_dict`.|
| `vmr_dict      `  | Gas volume mixing ratios (=mole fractions) at the surface. Must be set alongside `p_surf`. |
| `vmr_file      `  | Path to a file containing mixing ratio profiles. Replaces `vmr_dict`. |
| `chemistry     `  | Type of chemistry to be used (see below). |
| `condensates   `  | List of volatiles which are allowed to condense. Incompatible with `chemistry > 0`. |
| `transparent   `  | Make the atmosphere transparent (see below). Replaces all of the above parameters in this table. |


### `[execution]`
Parameters that tell the model what to do.

| Parameter         | Description   |
| ----------------: | :------------ |
| `clean_output  `  | Clean old files from the output folder at model startup (true/false). |
| `verbosity     `  | Logging verbosity (0: quiet, 1: normal, 2: extra logging) |
| `max_steps     `  | Maximum number of steps the solver should take before giving up (typically <200). |
| `max_runtime   `  | Maximum wall-clock runtime [s] before giving up. |
| `num_levels    `  | Number of model levels. Typically ~50, and ideally less than 150.  |
| `continua      `  | Include collisional/continuum absorption in radiative transfer (true/false) |
| `rayleigh      `  | Include Rayleigh scattering in radiative transfer (true/false) |
| `cloud         `  | Include cloud scattering and opacity in radiative transfer (true/false) |
| `overlap_method`  | Method for treating overlapping gas opacities within a given spectral band (see below) |
| `real_gas      `  | Use real-gas equation(s) of state where possible (true/false) |
| `thermo_funct  `  | Use temperature-dependent thermodynamic properties (true/false) |
| `sensible_heat `  | Include turbulent sensible heat transport at the surface (true/false) |
| `convection    `  | Include vertical heat transport associated with convection (true/false) |
| `latent_heat   `  | Include vertical heat transport from condensation and evaporation (true/false) |
| `rainout       `  | Enable compositional rainout of condensables. If disabled, phase change does not impact composition. |
| `initial_state `  | Ordered list of requests describing the initial state of the atmosphere (see below). |
| `solution_type `  | Solution type (see below). |
| `solver        `  | Solver to use (see below). |
| `dx_max        `  | Maximum step size [kelvin] allowed to be taken by the solver during each step. |
| `linesearch    `  | Linesearch method to be used (0: None, 1: Golden section, 2: Backtracking) |
| `easy_start    `  | Initially down-scale convective/condensation fluxes, if initial guess is poor/unknown. **Enable if the model is struggling.** |
| `converge_atol `  | Convergence criterion, absolute amount of energy flux lost [W m-2]. |
| `converge_rtol `  | Convergence criterion, relative amount of energy flux lost [dimensionless]. |
| `perturb_all`     | Perturb all rows of jacobian matrix at each solver iteration? True=stable, False=fast. |
| `rfm_wn_min`      | Line-by-line RFM radiative transfer, minimum wavenumber [cm-1] |
| `rfm_wn_max`      | Line-by-line RFM radiative transfer, maximum wavenumber [cm-1] |

### `[plots]`
Configure plotting routines all of these should be `true` or `false`.

| Parameter         | Description   |
| ----------------: | :------------ |
| `at_runtime     ` | Make some plots at runtime? |
| `temperature    ` | Plot temperature-pressure profile? |
| `fluxes         ` | Plot energy flux profiles? |
| `contribution   ` | Plot spectral contribution function? |
| `emission       ` | Plot outgoing emission spectrum? |
| `albedo         ` | Plot spectral albedo? This is the ratio of upward:downward SW fluxes |
| `mixing_ratios  ` | Plot mixing ratio profiles? |
| `height         ` | Plot radius-pressure profile? |
| `animate        ` | Make an animation of the solver obtaining its solution? |

### Details on specific parameters
* `composition.transparent` configures the atmosphere to be transparent. This works by setting the pressure to be small, and turning off the gas opacity. With this provided, the rest of the parameters in `[configuration]` are redundant. With this enabled, make sure to use the appropriate solver in the `[execution]` table.

* `execution.solution_type` tells the model which state to solve for. The allowed values (integers) are...
     - 1 : zero flux divergence at fixed `tmp_surf`
     - 2 : zero flux divergence such that the conductive skin (CBL) conserves energy flux with fixed `tmp_magma`
     - 3 : the net flux (up minus down) at each layer is equal to `flux_int`

     See the [Model description](@ref) page for an explanation of these solution types.

* `execution.solver` tells the model which solver to use. Allowed solvers are...
     - [empty string] : no solving takes place, so the model just calculates fluxes using the initial state
     - `newton` : the Newton-Raphson algorithm is used
     - `gauss`  : the Gauss-Newton algorithm is used
     - `levenberg` : the Levenberg–Marquardt algorithm is used
     - `transparent` : solver to be used when `composition.transparent=true`.

* `execution.initial_state` describes the initial temperature profile applied to the atmosphere. This is a list of strings which are applied in the given order, which allows the user to describe a specific state as required. The descriptors are listed below, some of which take a single argument that needs to immediately follow the descriptor in the list order.
     - `dry`              : integrate the dry adiabatic lapse rate from the surface upwards
     - `str`,       `arg` : apply an isothermal stratosphere at `arg` kelvin
     - `iso`,       `arg` : set the whole atmosphere to be isothermal at `arg` kelvin
     - `csv`,       `arg` : set the temperature profile using the CSV file at the file path `arg`
     - `sat`,       `arg` : apply Clausius-Clapeyron saturation curve for the gas `arg`
     - `ncdf`,      `arg` : load profile from the NetCDF file located at `arg`
     - `loglin`,    `arg` : log-linear profile between `tmp_surf` at the bottom and `arg` at the top
     - `ana`              : use the Guillot ([2010](https://arxiv.org/abs/1006.4702)) analytical temperature solution

    For example, setting `initial_state = ["dry", "sat", "H2O", "str", "180"]` will set T(p) to follow the dry adiabat from the surface, the water condensation curve above that, and then to be isothermal at 180 K until the top of the model.

* `composition.chemistry` describes the type of chemistry to implement within the model. This is handled externally by FastChem, so you must set the environment variable `FC_DIR` to point to the FastChem directory. The allowed values (integers) are...
     - 0 : Disabled
     - 1 : Equilibrium, gas phase only
     - 2 : Equilibrium, with condensation (condensates retained)
     - 3 : Equilibrium, with condensation (condensates rained out)

     More information on the chemistry is available on the [Equilibrium chemistry](@ref) page

* `execution.overlap_method` tells SOCRATES which algorithm to use to combine gas opacities. The spectral files contain k-tables for pure gases, and combining these coefficients can be done in several ways. See Amundsen+[2017](https://www.aanda.org/articles/aa/full_html/2017/02/aa29322-16/aa29322-16.html) for a nice comparison of overlap methods. Allowed options are...
     - `"ee"`   : equivalent extinction (fastest)
     - `"rorr"` : random overlap, with resorting and re-binning
     - `"ro"`   : random overlap (most accurate, slowest)

## Accessing AGNI from Python
It is possible to interact with AGNI from Python. This is best done with the `juliacall` package from [PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl).

Coupling with Python is done via `juliacall` within the modular [PROTEUS framework](https://github.com/FormingWorlds/PROTEUS), which couples AGNI self-consistently to models of planetary interior evolution and volatile outgassing. You can see this implemented in [agni.py](https://github.com/FormingWorlds/PROTEUS/blob/main/src/proteus/atmos_clim/agni.py) within the PROTEUS source code.

A skeleton example is given below:
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

## Line-by-line radiative transfer
Performed using the [Reference Forward Model](https://eodg.atm.ox.ac.uk/RFM/).
You must provide a HITRAN-formatted `.par` file, setting the path via `files.rfm_parfile`.
This parfile can contain absorption from multiple species, and can be obtained from [hitran.org](https://hitran.org/lbl/).
Alternatively, get the parfiles stored on Zenodo using: `./src/get_data.sh parfiles`.

Then, you must also set the variables `execution.rfm_wn_min` and `execution.rfm_wn_max`.
These two parameters specify the wavenumber [cm-1] range over which to perform the LbL calculations.
The wavenumber resolution is set to 1 cm-1.

Results are saved to the NetCDF file, alongside all the usual data, as `rfm_wn` and `rfm_fl` arrays.
