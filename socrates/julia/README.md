# SOCRATES.jl Julia wrapper for SOCRATES radiative transfer codes

This provides a Julia module wrapping the SOCRATES Fortran `radiance_core` code.

SOCRATES Types (control, radout, etc) can be examined and modified from the Julia REPL
(fields appear as Julia properties, ie use . instead of %).

SOCRATES functions (`radiance_core`) can be called using the same syntax as Fortran.

## Building the wrappers
    
1. Check out SOCRATES (revision um13.0 from 2022-07-28 used for these tests)
    
2. Modify SOCRATES make/Mk_cmd_xxx to add -fPIC flag (for gfortran) to
  build 'position independent code' for a shared library. and
  build SOCRATES as usual (we need the radlib.a library)
    
3. Regenerate Fortran and Julia wrappers:
        a) Check  ./julia/src/generate_wrappers.jl
             SOCRATES_DIR  - path to SOCRATES source tree on local machine
             SVN_REV  - (integer) svn revision
        b) Create an empty ./julia/gen folder (to contain generated wrappers)
        c) Activate Julia environment in ./julia folder:
            julia> ] activate .
            (SOCRATES) instantiate
        d) Build wrappers
            julia> cd("src")
            julia> include("generate_wrappers.jl")
  
    
4. Build shared library in julia/lib, using 'make'
    from linux command prompt (tested on Ubuntu linux, Makefile
    will probably need modifying for mac etc)
    This builds a shared library that includes the
    SOCRATES radlib.a + extra C wrapper functions

    NB: RAD_BIN must be set eg by . ./set_rad_env  (see note in Makefile)
    

## Examples
### SOCRATES clear sky minimal test case (LW flux, no absorption)
    
Minimal example to demonstrate Julia <--> Fortran integration

Example from https://code.metoffice.gov.uk/trac/socrates/wiki/SocratesDoc/ClearskyFluxExample
(clear sky fluxes, based on examples/netcdf/CIRC_case6)

Produces expected result (eg flux_down = 0, flux_up = constant for lw, no absorber case)

To run the Julia example script:
- Edit julia/test/test_clearsky.jl to set path
  to SOCRATES source tree on local machine (for data and example files)

- From Julia, starting in the ./julia folder:
        julia> ] activate .
        julia> cd("test/")
        julia> include("test_clearsky.jl")

## Implementation

The overall strategy is to map Fortran types into Julia handle or proxy types that contain
an opaque C pointer holding the memory reference to the Fortran instance.

SOCRATES types are wrapped by autogenerating code for a Julia handle type (holding a pointer
to the Fortran type), and get / set functions. Two files are generated for each Type: a 
Fortran file providing C-callable functions using ISO_C_BINDING, and a Julia file with
a Julia type definition, and functions including `getproperty`, `setproperty` and `propertynames`.

SOCRATES module parameters (for integer and string constants) are translated into Julia modules
with corresponding const variables.

Top-level SOCRATES functions (`calc_radiance`) etc are hand-written.
