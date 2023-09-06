# AGNI
Radiative-convective solver which uses SOCRATES (2306) for radiative-transfer.    
Pronounced: *ag-nee*. Named after the fire deity of Hinduism.    
The main repository for this code is [github.com/nichollsh/AGNI](https://github.com/nichollsh/AGNI).

### Contributors (abbreviations & email addresses)
* HN - Harrison Nicholls (harrison.nicholls@physics.ox.ac.uk)  

### Repository structure 
* `README.md`       - This file
* `LICENSE.txt`     - License for use and re-use
* `doc/`            - Other documentation
* `out/`            - Output files
* `res/`            - Resources
* `src/`            - AGNI source code
* `socrates/`       - Directory containing SOCRATES and associated files
* `agni.jl`         - Main AGNI run file
* `demo_steamrun.jl`- Demo AGNI run file (steam runaway greenhouse effect)


### Requirements
* Julia (version 1.9.1 or later)
* Python (version 3.10 or later)
* NumPy and SciPy
* gfortran
* NetCDF
* netcdf-fortran
* OpenMP

### Supported platforms
* MacOS (ARM and x86-64)
* Ubuntu (x86-64)

### Installation instructions
- `$ cd socrates`
- `$ cp ../res/Mk_cmd_PLAT ./make/Mk_cmd` where PLAT is your platform
- `$ ./build_code`
- `$ cd julia`
- `$ julia`
- `julia> ]`
- `(@v1.9) pkg> add OffsetArrays`
-  `(@v1.9) pkg> add Revise`
-  `(@v1.9) pkg> add PCHIPInterpolation`
-  `(@v1.9) pkg> add LaTeXStrings`
-  `(@v1.9) pkg> add Plots`
-  `(@v1.9) pkg> activate .`
-  Press backspace
-  `julia> cd("src")`
-  `julia> include("generate_wrappers.jl")`
-  `julia> exit()`
-  `$ cd lib`
-  `$ make`
-  `$ cd ..`   
You should end up in a directory called 'julia'.

### Run instructions
Simply run `$ julia agni.jl` in the root directory of the repository.

### License
See [LICENSE.txt](LICENSE.txt) for AGNI.      
See [socrates/COPYRIGHT.txt](socrates/COPYRIGHT.txt) for SOCRATES.   

