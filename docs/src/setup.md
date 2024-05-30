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
- `(@v1.10) pkg> activate .` ← note the dot
- `(AGNI) pkg> instantiate`
- `(AGNI) pkg> build`
- Press backspace
- `julia> exit()` 
You should run the tests next.

## SOCRATES 
The SOCRATES radiative transfer school is packaged separately. If you haven't already installed it, this can be done by running `./get_socrates.sh`. 

## Testing
- `$ julia`
- `julia> ]`
- `(@v1.10) pkg> activate .`
- `(AGNI) pkg> test`
