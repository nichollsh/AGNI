# Troubleshooting
This page may be useful if you are having problems. However, I would suggest that 
you also double check the [Getting started](@ref) instructions.

## Julia errors on start, potentially referencing the CURL library
It is important that the shell environment variable `LD_LIBRARY_PATH` is 
not set when running AGNI. This will cause Julia to use the wrong libraries, 
which will causes problems. You can unset this variable or reset using either of the
following commands 
```bash
unset LD_LIBRARY_PATH
export LD_LIBRARY_PATH=""
```

## Cannot find SOCRATES
Check the installation instructions. Have you set `RAD_DIR`? Try running
`l_run_cdf` in the terminal; if this fails, then SOCRATES has not compiled
or you haven't added it to your `PATH`.


## Finally...
If you are still stuck, or feel that there is a problem with the code, then 
you can contact the authors using the information on the main page.


