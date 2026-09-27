# Refractive index and optical data for aerosols

Refractive indices (optical constants n, k) of aerosol and cloud condensate materials.

Download these files with `./src/get_data.sh refractive`.

These are used by AGNI to calculate aerosol optical properties at runtime using Mie theory
(see src/energy/aerosol_optics.jl), for aerosols configured with method="mie".

Each file <id>.txt contains three numeric columns: wavelength [micron], n, k. The <id> is
the material name used in the `nk_file` config key. Header lines begin with '#'.
