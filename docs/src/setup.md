# Getting started
This page outlines requirements and installation steps for the code. Currently, 
GNU/Linux and MacOS (including ARM) are supported. 

## Requirements
* gfortran 
* NetCDF library for FORTRAN
* make
* curl

!!! warning
    Do not install Julia using your system package manager. Install only from julialang.org

## Installation
Follow the steps below in order to setup the code.
1. Install Julia: `curl -fsSL https://install.julialang.org | sh`
2. Download AGNI: `git clone https://github.com/nichollsh/AGNI.git`
3. Change directory: `cd AGNI`
4. Setup SOCRATES by doing either **ONE** of the following...
    - Follow the instructions on the SOCRATES GitHub page   
    - Run `source get_socrates.sh`    
5. `julia -e 'using Pkg; Pkg.activate("."); Pkg.build()'`
AGNI is now installed as a package into a Julia environment in the AGNI
directory. This will also have downloaded some basic input data. 
You should run the tests next.

!!! tip 
    The `get_socrates` script automatically adds the radiation code to your
    environment with the variable `RAD_DIR`. This needs to be set whenever 
    AGNI is being used.

## Testing
Now try running the tests in your terminal. 
```bash 
julia ./test/runtests.jl
```
This will print information on whether tests passed or failed.   

## Using the code
See [Running the model](@ref) for information on using the code.    
See [Troubleshooting](@ref) for troubleshooting advice. 
