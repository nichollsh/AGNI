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

## What next?

To run AGNI with the default configuration file:
```bash
./agni.jl res/config/default.toml
```

Output files are written to the directory specified in the configuration file (default: `out/`).

See [Tutorials](@ref) for a guided introduction to AGNI with example outputs results. Refer to the other [How-to guides](@ref) for steps on configuring the model for your own science cases.
