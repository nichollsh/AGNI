# Getting started
This page outlines requirements and installation steps for the code.

## Requirements
* Julia (NB: install only from julialang.org - do not use your system package manager)
* [SOCRATES](https://github.com/nichollsh/SOCRATES)

## Supported platforms
* MacOS (ARM and x86-64)
* GNU/Linux (x86-64)

## Installation
1. Setup SOCRATES by doing either **ONE** of the following...
    - Follow the instructions on the SOCRATES GitHub page   
    - Run `./get_socrates.sh`    
2. `$ julia`
3. `julia> ]` 
4. `(@v1.10) pkg> activate .` ← note the dot
5. `(AGNI) pkg> build`
6. Press backspace
7. `julia> exit()` 
AGNI is now installed as a package into a Julia environment in the AGNI
directory. You should run the tests next.

## Testing
Run `./test/runtests.jl ` in your terminal. This will print information 
on whether tests passed or failed.   


## Using the code
See [Running the model](@ref) for information on using the code.
