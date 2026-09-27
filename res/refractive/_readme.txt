Refractive indices (optical constants n, k) of aerosol and cloud condensate materials.

These are used by AGNI to calculate aerosol optical properties at runtime using Mie theory
(see src/energy/aerosol_optics.jl), for aerosols configured with method="mie".

Each file <id>.txt contains three numeric columns: wavelength [micron], n, k. The <id> is
the material name used in the `nk_file` config key. Header lines begin with '#'.

The data are a copy of the refractive index compilation distributed with POSEIDON
(Mullens, Lewis & MacDonald 2024, ApJ 977, 105), which itself collates data from
Wakeford & Sing (2015), Kitzmann & Heng (2018), Burningham et al. (2021), gCMCRT
(Lee et al.) and others. See _manifest.txt and Aerosol-Database-Readme.txt for provenance.

Download these files with `./src/get_data.sh refractive`.
