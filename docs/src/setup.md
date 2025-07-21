# Getting started
This page outlines requirements and installation steps for the code. Currently,
GNU/Linux and MacOS (including ARM) are supported.

## Software requirements
* gfortran
* netcdf
* netcdf-fortran
* make
* wget
* unzip

!!! warning
    Do not install Julia using your system package manager. Install only from julialang.org as below.

## Installation
Follow the steps below in order to setup the code.
1. Install Julia: `curl -fsSL https://install.julialang.org | sh`
2. Download AGNI: `git clone https://github.com/nichollsh/AGNI.git`
3. Change directory: `cd AGNI`
4. Setup SOCRATES by doing either **ONE** of the following...
    - Follow the instructions on the [SOCRATES GitHub](https://github.com/nichollsh/SOCRATES) page
    - Run `./src/get_socrates.sh`
5. Finally, install AGNI: `./src/get_agni.sh`
AGNI is now installed as a package into the Julia environment of the AGNI
directory. This will also have downloaded some basic input data and run the tests.

!!! warning
    The `RAD_DIR` environment variable must be set to your SOCRATES path whenever AGNI is being used.
    This is so that AGNI can locate the SOCRATES libraries.

## Testing
If you want to run the tests manually, simply use the script in the `test/` folder...
```bash
julia test/runtests.jl
```
This will print information on whether tests passed or failed.

## Updating
It's important that you keep AGNI up to date, especially if you are using as part of
the [PROTEUS framework](https://github.com/FormingWorlds/PROTEUS). Use this script to
automatically pull changes from GitHub and download any required data files.
```bash
./src/get_agni.sh
```

## Using the code
See [Using the model](@ref) for information on using the code.
See [Troubleshooting](@ref) for troubleshooting advice.


## Coupling with FastChem
Coupling with FastChem can be enabled using the configuration parameter `composition.chemistry`.
Of course, it is first necessary to set up FastChem, which can be done by running:
```bash
./src/get_fastchem.sh
```
You **must** then set the `FC_DIR` environment variable to the location of the FastChem
installation folder. Ideally you should also add this variable to your bashrc file.
