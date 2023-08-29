# Themis
Radiative-convective solver written in Julia. Uses SOCRATES for radiative-transfer.

### Install

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

### Running the radiative-convective model

Simply run `$ julia themis.jl` in the root directory of the repository.

### License

See `LICENSE.txt`.



