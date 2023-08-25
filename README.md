# Themis
Radiative-convective solver written in Julia. Uses SOCRATES for radiative-transfer.


### Install

1. `$ cd socrates`
2. `$ unzip ../res/socrates_revision*.zip`
3. `$ cp ../res/Mk_cmd_PLAT ./make/Mk_cmd` where PLAT is your platform
4. `$ ./build_code`
5. `$ source set_rad_env`
6. `$ cd julia`
7. `$ julia`
8. `julia> ]`
9. `(@v1.9) pkg> add OffsetArrays`
10. `(@v1.9) pkg> activate .`
11. Press backspace
12. `julia> cd("src")`
13. `julia> include("generate_wrappers.jl")`
14. `julia> exit()`
15. `$ cd lib`
16. `$ make`
17. `$ cd ..`
18. You should end up in a directory called 'julia'.


### Testing the SOCRATES wrapper

Start in the 'julia' folder.
1. `$ julia`
2. `julia> ]`
3. `(@v1.9) pkg> activate .`
4. Press backspace
5. `julia> cd("test")`
6. `julia> include("test_clearsky.jl")`
The code works if you end up with the Julia REPL containing arrays of fluxes.

### Running the radiative-convective model


