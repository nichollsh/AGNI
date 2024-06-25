# Getting started
This page outlines requirements and installation steps for the code. Currently, 
GNU/Linux and MacOS (including ARM) are supported. 

## Requirements
* [Julia](https://julialang.org/downloads/) 
* [SOCRATES](https://github.com/nichollsh/SOCRATES) - see instructions below

!!! warning
    Do not install Julia using your system package manager. Install only from julialang.org

## Installation
1. Download AGNI from GitHub, either as a zip from the website or using Git 
    `git clone https://github.com/nichollsh/AGNI.git`
1. Setup SOCRATES by doing either **ONE** of the following...
    - Follow the instructions on the SOCRATES GitHub page   
    - Run `source get_socrates.sh`    
2. `$ julia`
3. `julia> ]` 
4. `(@v1.10) pkg> activate .` ← note the dot
5. `(AGNI) pkg> build`
6. Press backspace to exit the package manager
7. `julia> exit()` 
AGNI is now installed as a package into a Julia environment in the AGNI
directory. You should run the tests next.

!!! tip 
    The `get_socrates` script automatically adds the radiation code to your
    environment with the variable `RAD_DIR`. This needs to be set whenever 
    AGNI is being used.

## Testing
Run `julia ./test/runtests.jl ` in your terminal. This will print information 
on whether tests passed or failed.   

## Using the code
See [Running the model](@ref) for information on using the code.
