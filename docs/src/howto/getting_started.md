# Getting started
This page outlines requirements and installation steps for the code. Currently,
GNU/Linux and MacOS (including ARM) are supported.

## Software requirements
* gfortran
* netcdf
* netcdf-fortran
* make
* wget
* curl
* unzip
* cmake

!!! warning
    Do not install Julia using your system package manager. Install only from julialang.org as below.

## Installation

Follow the ordered steps below
1. Install Julia's package manager
    - `curl -fsSL https://install.julialang.org | sh`

2. Download AGNI
    - `git clone https://github.com/nichollsh/AGNI.git && cd AGNI`

3. Setup SOCRATES by doing either **ONE** of the following...
    - a) Follow the instructions on the [SOCRATES GitHub](https://github.com/FormingWorlds/SOCRATES) page
    - b) **OR**, run `./src/get_socrates.sh`

### SOCRATES floating point precision
SOCRATES can be compiled in single or double precision. The Julia interface exposes this
via `SOCRATES_REAL_BYTES`:
- `4` means single precision.
- `8` means double precision.

AGNI records the value in outputs as `SOCRATES_precision`. If the SOCRATES Julia
interface does not define `SOCRATES_REAL_BYTES`, AGNI assumes double precision (`8`).

4. Setup FastChem
    - `./src/get_fastchem.sh`

5. Install AGNI
    - `./src/get_agni.sh`

AGNI is now installed as a package of the Julia environment in the `AGNI/`
directory. These steps will have downloaded some basic input data.

!!! warning
    The `RAD_DIR` environment variable must be set to your SOCRATES path whenever AGNI is being used. This is so that AGNI can locate the SOCRATES libraries.

!!! tip
    Visit the [Troubleshooting](@ref) page if you encounter any problems. This page can
    usually resolve your problem.

## Testing
If you want to run the tests manually, simply use the script in the `test/` folder...
```bash
julia --project test/runtests.jl
```
This will print information on whether tests passed or failed.

## Updating
It's important that you keep AGNI up to date, especially if you are using as part of
the [PROTEUS framework](https://github.com/FormingWorlds/PROTEUS). Use this script to
automatically pull changes from GitHub and download any required data files.
```bash
./src/get_agni.sh
```

## Your first model run

The environment variable `RAD_DIR` must point to the SOCRATES installation
directory. The best way to do this is to add `RAD_DIR=path/to/socrates/folder/` to your
shell rc file (e.g. `~/.bashrc`).

Then run AGNI with the default configuration file:
```bash
./agni.jl
```

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

Output files are written to the directory specified in the configuration file (default: `out/`).
See [Tutorials](@ref) for illustrative results, and [How-to guides](@ref) for
next steps such as configuring the model for your own science case.
