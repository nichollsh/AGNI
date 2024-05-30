# Getting started
This page outlines requirements and installation steps for the code.

## Requirements
* Julia (NB: install only from julialang.org - do not use your system package manager)
* SOCRATES

## Supported platforms
* MacOS (ARM and x86-64)
* GNU/Linux (x86-64)

## Installation
- Set SOCRATES environment variables.
- `$ julia`
- `julia> ]` 
- `(@v1.10) pkg> activate .`
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

## SOCRATES 
The SOCRATES radiative transfer school is packaged separately. If you haven't already installed it, this can be done by running `./get_socrates.sh`. 


## Testing
- `$ julia`
- `julia> ]`
- `(@v1.10) pkg> activate .`
- `(AGNI) pkg> test`
