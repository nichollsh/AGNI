# THEMIS
Radiative-convective solver which uses SOCRATES (2306) for radiative-transfer.    
Pronounced: *them-iss*. Named after the goddess of justice and order.  

### Contributors (abbreviations & email addresses)
* HN - Harrison Nicholls (harrison.nicholls@physics.ox.ac.uk)  

### Repository structure 
* `README.md`       - This file
* `LICENSE.txt`     - License for use and re-use
* `doc/`            - Other documentation
* `out/`            - Output files
* `res/`            - Resources
* `src/`            - THEMIS source code
* `socrates/`       - Directory containing SOCRATES and associated files
* `themis.jl`       - THEMIS run file


### Requirements
* Julia (version 1.9.1 or later)
* gfortran
* NetCDF
* netcdf-fortran
* OpenMP

### Supported platforms
* MacOS (ARM and x86-64)
* Ubuntu (x86-64)

### Installation instructions
1. `$ cd socrates`
3. `$ cp ../res/Mk_cmd_PLAT ./make/Mk_cmd` where PLAT is your platform
4. `$ ./build_code`
5. `$ source set_rad_env`
6. `$ cd julia`
7. `$ julia`
8. `julia> ]`
9. `(@v1.9) pkg> add OffsetArrays`
10. `(@v1.9) pkg> add Revise`
11. `(@v1.9) pkg> activate .`
12. Press backspace
13. `julia> cd("src")`
14. `julia> include("generate_wrappers.jl")`
15. `julia> exit()`
16. `$ cd lib`
17. `$ make`
18. `$ cd ..`
19. You should end up in a directory called 'julia'.

### Run instructions
Simply run `$ julia themis.jl` in the root directory of the repository.

### License
See [LICENSE.txt](LICENSE.txt) for THEMIS.      
See [socrates/COPYRIGHT.txt](socrates/COPYRIGHT.txt) for SOCRATES.   

