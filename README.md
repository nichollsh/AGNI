# AGNI
Radiative-convective solver designed for integration into a coupled atmosphere-interior code.   

AGNI relies on SOCRATES (2311) for calculating radiances. The radiative transfer includes shortwave irradiation from the star, surface emission, gaseous absorption, Rayleigh scattering, parameterised clouds, and continuum absorption. Mixing length theory is used to parametrise convection. Together, energy transport processes allow for an energy-conserving calculation of the atmosphere's temperature profile.      

The model is distributed under a proprietary license. Only once it has been published will the model be distributed under a FOSS license. If you use the model in a publication (once it is open), please cite my paper describing the model.

Consult the [AGNI Wiki](https://github.com/nichollsh/AGNI/wiki) on GitHub for information about the model. 
    
## Repository structure 
* `README.md`       - This file
* `LICENSE.txt`     - License for use and re-use
* `doc/`            - Further documentation
* `out/`            - Output files
* `res/`            - Resources
* `src/`            - AGNI source code
* `socrates/`       - Directory containing SOCRATES and associated files (subject to the license therein)
* `agni.jl`         - AGNI executable
* `demo_steamrun.jl`- Script to demonstrate the pure-steam runaway greenhouse effect
* `demo_tests.jl`   - Script containing quick tests for verifying that the basics of the model are functioning


## Requirements
* Julia (version 1.9.1 or later)
* gfortran (NB: do not use ifort or aocc)
* netcdf-fortran
* make
* OpenMP

## Supported platforms
* MacOS (ARM and x86-64)
* GNU/Linux (x86-64)

## Installation instructions
- `$ cd socrates`
- `$ ./configure`
- `$ ./build_code`
- `$ source set_rad_env`
- `$ cd julia`
- `$ julia`
- `julia> ]`
-  `(@v1.9) pkg> add OffsetArrays Revise PCHIPInterpolation LaTeXStrings Plots NCDatasets DataStructures Glob ArgParse BinnedStatistics LoggingExtras`
-  `(@v1.9) pkg> activate .`
-  Press backspace
-  `julia> cd("src")`
-  `julia> include("generate_wrappers.jl")`
-  `julia> exit()`
-  `$ cd lib`
-  `$ make`
-  `$ cd ../../..`   
You should end up in the root directory of the repository.    

## Running the code
To run the program, execute `./agni.jl [cfg_path]`. If `[cfg_path]` is not provided, then the default configuration file will be used.       
To demo the steam runaway greenhouse effect, execute `./demo_steamrun.jl`.     
To run the unit tests execute `./demo_steamrun.jl`.     

