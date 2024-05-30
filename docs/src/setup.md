# Getting started
This page outlines requirements and installation steps for the code.

## Requirements
* Julia (NB: install only from julialang.org - do not use your system package manager)
* [SOCRATES](https://github.com/nichollsh/SOCRATES)

## Supported platforms
* MacOS (ARM and x86-64)
* GNU/Linux (x86-64)

## Installation
- Setup SOCRATES by doing either **ONE** of the following...
    - Follow the instructions on the SOCRATES GitHub page
    - Run `./get_socrates.sh`
- `$ julia`
- `julia> ]` 
- `(@v1.10) pkg> activate .` ← note the dot
- `(AGNI) pkg> build`
- Press backspace
- `julia> exit()` 
AGNI is now installed as a package into a Julia environment in this directory.   
You should run the tests next.

## Testing
- `$ julia`
- `julia> ]`
- `(@v1.10) pkg> activate .`
- `(AGNI) pkg> test`
These tests may trigger recompilation of depenencies.
