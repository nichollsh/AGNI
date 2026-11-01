# Obtaining input data

The minimal input data required to run the model will have been downloaded automatically
from Zenodo during installation. If you require more data, such as additional stellar
spectra or opacities, these can also be obtained using the `get_data` script in the AGNI
root directory. To see how to use this script, run it without arguments:
```bash
./src/get_data.sh
```

## Spectral files

Opacities are contained within "spectral files". Each spectral file is identified by a
codename and a number of bands, and consists of two files: `<Codename>.sf` and its
k-coefficient tables `<Codename>.sf_k`.

The spectral files available for AGNI are described in the
[SOCRATES documentation](https://proteus-framework.org/SOCRATES/). Use it to decide which spectral file is best for you:

- [PROTEUS spectral files](https://proteus-framework.org/SOCRATES/Reference/proteus_spectral_file_reference.html)
  lists every codename with its gases, continua, band counts, data sources, references and  Zenodo records, and gives advice on
  [choosing a spectral file](https://proteus-framework.org/SOCRATES/Reference/proteus_spectral_file_reference.html#choosing-a-spectral-file).
- [Spectral files](https://proteus-framework.org/SOCRATES/Explanations/spectral_files.html) explains what each spectral file contains.
- [Validation](https://proteus-framework.org/SOCRATES/Reference/spectral_file_validation.html) compares the Generation III spectral files with published runaway greenhouse calculations, and shows how their results change with the number of wavenumber bands.

Download a spectral file by passing its codename and number of bands to the `get_data` script. For example, to get the spectral file "Honeyside48" you would run:
```bash
./src/get_data.sh anyspec Honeyside 48
```
To get the Generation III water spectral file "Snowline" with 128 bands:
```bash
./src/get_data.sh anyspec Snowline 128
```
The files are placed in `res/spectral_files/<Codename>/<bands>/`. 

Other spectral files can be downloaded directly from the
[PROTEUS community on Zenodo](https://zenodo.org/communities/proteus_framework/records?q&f=subject%3Aspectral_files&l=list&p=1&s=10&sort=newest).

!!! tip "Get missing spectral files"
    See [Spectral file does not exist](@ref) in the troubleshooting guide if a spectral
    file you downloaded cannot be found by AGNI.

## Aerosol refractive indices

Aerosols with `method = "mie"` require refractive index data for the particle material.
These can be obtained by running:
```bash
./src/get_data.sh refractive
```
which places the files in `res/refractive/`. See [Aerosols and clouds](@ref) for their
provenance, and for the list of supported materials.

## Data in another folder

AGNI reads the `thermodynamics`, `scattering`, `refractive` and `blobs` folders of `res/`
through `paths.get_dir`, and each can be read from elsewhere. For each of them, AGNI uses
the first of: the environment variable `AGNI_DIR_<name>` (for example
`AGNI_DIR_refractive`), the folder `<name>` inside `AGNI_DIR_res`, the folder `<name>`
inside `res_dir` from the `[files]` section of the configuration, and finally `res/<name>`.
`get_dir` also accepts `config`, `stellar_spectra` and `spectral_files` for external
callers. Blank variables are ignored, and relative paths are taken from the working
directory. Each folder placed this way is reported once in the log when the atmosphere is
set up. `get_data.sh` always writes to `res/`, and the file paths in `[files]` (such as
`input_sf`) are used as written.

## Mirror on DataverseNL

The Zenodo records used by `get_data.sh` are mirrored on
[DataverseNL](https://dataverse.nl). When Zenodo cannot be reached, or a download from
Zenodo fails, the script takes the same files from the mirror of that record, and it
reads the remaining records of that run from their mirrors too. A record without a mirror
fails with an error. The environment variables `ZENODO_URL` and `DATAVERSE_URL` set the
two servers.
