# Troubleshooting
This page may be useful if you are having problems. However, I would suggest that you also double check that you followed all of the [Getting started](@ref) instructions.

## Curl is not installed
You need to install [curl](https://curl.se/). This is a tool for tranferring files over a network. Curl is used by AGNI to obtain lookup data files from Zenodo, and by Julia to obtain the required libraries.

To install curl on Ubuntu:
```bash
sudo apt install curl
```

For other Linux distributions, see the [curl download page](https://curl.se/download.html). MacOS comes with curl pre-installed.

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

## NetCDF is not installed
You need to install NetCDF on your machine. This is a library designed for reading and writing
data files, commonly used in atmospheric sciences. [Wikipedia page](https://en.wikipedia.org/wiki/NetCDF).

To install NetCDF on Ubuntu:
```bash
sudo apt install libnetcdf-dev netcdf-bin ncview libnetcdff-dev
```

To install NetCDF on MacOS:
```bash
sudo brew install netcdf netcdf-fortran
```


## Cannot find SOCRATES
Check the installation instructions. Have you set `RAD_DIR`? Try running
`l_run_cdf` in the terminal; if this fails, then SOCRATES has not compiled
or you haven't added it to your `PATH`. It is necessary to set the `RAD_DIR` variable
for the environment in which you are running AGNI, so it is best to add it to your shell's
rc file permanently.


## Spectral file does not exist
First, check the path in the configuration file.

Download additional spectral files using the `get_data` script.
For example, for additional pure-steam spectral files you would run:
```bash
./src/get_data.sh steam
```

When you downloaded AGNI it should have obtained a "basic" set of data. This will include
a reference guide located at `res/spectral_files/reference.pdf`. Using the table inside
this PDF file, you can decide which set of opacities are appropriate for you.

For example, to download the spectral file `Honeyside16` you would then run:
```bash
./src/get_data.sh anyspec Honeyside 16
```
Note the space between the codename and number of bands.
Other spectral files can be downloaded from the [PROTEUS community on Zenodo](https://zenodo.org/communities/proteus_framework/records?q&f=subject%3Aspectral_files&l=list&p=1&s=10&sort=newest).

## Cannot find FastChem
You need to install FastChem. This can be done by running the command:
```bash
./src/get_fastchem.sh
```
and then adding `FC_DIR` to your shell rc file.

## Finally...
If you are still stuck, or feel that there is a problem with the code, then
you can contact the authors using the information on the main page.

