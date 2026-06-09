# Radiative transfer

## Radiation physics
Radiative transfer (RT) refers to the transport of radiation energy through a medium subject to the characteristics of the medium [stamnes_radiative_2017](@citep). Radiation passing through an atmosphere is absorbed, emitted, scattered, and reflected. In the context of planetary atmospheres, we also have to handle their surfaces, cloud formation, and radiation from the host star.

## Correlated-k

AGNI nominally simulates RT using a bespoke version of [SOCRATES](https://proteus-framework.org/SOCRATES/): a suite of numerical codes primarily developed by the UK Met Office [manners_fast_2024,sergeev_socrates_2023](@citep).

SOCRATES solves the RT equation using a two-stream solution [edwards_studies_1996, edwards_efficient_1996](@citep). In SOCRATES, bolometric radiation fluxes are determined by summing over some set of spectral regions ('bands'). These bands are then split into sub-bands, which are treated mono-chromatically. Each monochromatic calculation is performed under a two-stream approximation [zdunkowski_twostream_1980](@citep), in which we integrate over the angular component of the radiation field to use a single column with multiple upward- and downward-directed beams [lee_testing_2024, zdunkowski_pifm_1985](@citep).

Opacity is handled using the correlated-k approximation [lacis_corrk_1991](@citep), with either random overlap or equivalent extinction used to account for overlapping absorption in mixtures of gases. The model has generally used k-terms fitted to gas line-absorption opacities from the [DACE](https://dace.unige.ch/opacityDatabase/?#) database [grimm_database_2021](@citep), which largely draws from ExoMol and HITEMP [tennyson_exomol_2018, rothman_hitemp_2010](@citep). These are cross-sections derived from integrating over the linelists, with a line truncation width of 25 cm$^{-1}$.

The MT\_CKD model is used to estimate water continuum absorption cross-sections [mlawer_mtckd_2012, mlawer_mtckd_2023](@citep). Other continua are derived from the HITRAN tables. Rayleigh scattering, water cloud radiative properties, and aerosol parametrisations are also included. For aerosols, generates band-averaged optical properties files at runtime using `scatter_average_90`, then inserts these into the runtime spectral file using `prep_spec`, alongside the stellar spectrum and Rayleigh scattering terms. You can find tools for fitting k-terms and processing line absorption data in my redistribution of [SOCRATES](https://github.com/FormingWorlds/SOCRATES) on GitHub. The flowchart below outlines how these absorption data are converted into a 'spectral file'.

The monochromatic scattering properties are stored in `SOCRATES/data/aerosol/*.mon` files. New `.mon` files can be generated using the `SOCRATES/sbin/Cscatter` script, which can allow the creation of new aerosol data.

![](fig_spectral_flowchart.svg)

## Surface reflectivity

Surface reflectivity can be modelled as a greybody with an albedo from 0 to 1 [Hapke_2012](@citep). Alternatively, the surface can be modelled using empirical reflectance data that varies (spectrally) with wavelength. In the latter case a filepath must be provided via the config. The file can tabulate any one of: spherical reflectance ('r'), hemispherical emissivity ('e'), or single scattering albedo ('w'). These data are compiled on Zenodo [here](https://zenodo.org/communities/proteus_framework/records?q&f=subject%3Asurface_albedos&l=list&p=1&s=10&sort=newest).

## Line by line

AGNI also includes an interface to the [Reference Forward Model](https://eodg.atm.ox.ac.uk/RFM/), which is packaged as a binary
blob since the RFM closed-source [dudhia_rfm_2017](@citep). This interface provides an extremely easy way to validate and benchmark SOCRATES.

## Stellar spectra

A key input to the radiation model is the shortwave downward-directed flux from the star at the top of the atmosphere. This is quantified by the bolometric instellation flux, a scale factor, an artificial additional albedo factor, and a zenith angle [cronin_choice_2014](@citep). All of these may be provided to the model through the configuration file. The model also requires a stellar spectrum scaled to the top of the atmosphere.


## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
