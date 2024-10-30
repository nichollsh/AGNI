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
If this does not help, it's possible that you are using a Julia distribution provided by
your system package manager. It's important that you only use Julia distributed from the
official website.


## Cannot find SOCRATES
Check the installation instructions. Have you set `RAD_DIR`? Try running
`l_run_cdf` in the terminal; if this fails, then SOCRATES has not compiled
or you haven't added it to your `PATH`. It is necessary to set the `RAD_DIR` variable
for the environment in which you are running AGNI, so it is best to add it to your shell's
rc file permanently.


## Spectral file does not exist
Check the path in the configuration file. Download additional spectral files using the
`get_data` script. For example, for additional pure-steam spectral files you would run
```bash
./src/get_data.sh steam
```

## Cannot find FastChem
You need to install FastChem. This can be done by running the command:
```bash
./src/get_fastchem.sh
```
and then adding `FC_DIR` to your shell rc file.

## Finally...
If you are still stuck, or feel that there is a problem with the code, then
you can contact the authors using the information on the main page.

