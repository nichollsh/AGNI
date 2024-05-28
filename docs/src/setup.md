# Getting started
This page outlines requirements and installation steps for the code.

## Requirements
* Julia (NB: install only from julialang.org - do not use your system package manager)
* gfortran (NB: do not use ifort or aocc)
* netcdf-fortran
* make
* OpenMP

## Supported platforms
* MacOS (ARM and x86-64)
* GNU/Linux (x86-64)

## Installation
- `$ export LD_LIBRARY_PATH=""`
- `$ cd socrates`
- `$ ./configure`
- `$ ./build_code`
- `$ source set_rad_env`
- `$ julia`
- `julia> ]` 
- `(@v1.10) pkg> activate ../`
- `(AGNI) pkg> instantiate`
- Press backspace
- `julia> cd("julia/src")`
- `julia> include("generate_wrappers.jl")`
- `julia> exit()`
- `$ cd julia/lib`
- `$ make`
- `$ cd ../../..`   
You will end up in the root directory of the repository.    
You should run the tests next.

## Testing
- `$ julia`
- `julia> ]`
- `(@v1.10) pkg> activate .`
- `(AGNI) pkg> test`
