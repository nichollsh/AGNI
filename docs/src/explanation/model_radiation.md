# Radiative transfer

## Radiation physics

Radiative transfer (RT) refers to the transport of radiation energy through a medium subject to the absorption, emission, and scattering properties of that medium [stamnes_radiative_2017](@citep).

Absorption of radiation along a path of length $dl$ through a gas of opacity $\kappa_\nu$ [$\mathrm{m}^2\ \mathrm{kg}^{-1}$] and mass density $\rho$ defines the **optical depth**:
```math
d\tau_\nu = \kappa_\nu \rho \, dl
```

The effective radiating level of the atmosphere (the *photosphere*) is taken to be where some measure of the optical depth $\tau_\nu \sim 1$. This could be a slant or radial measurement from the top-of-atmosphere (TOA), depending on the viewing geometry.

In addition to molecular line absorption, **collision-induced absorption (CIA)** arises from brief multi-molecular complexes formed during collisions and line broadening. CIA can contribute a significant greenhouse effect at high pressures [stamnes_radiative_2017, pierrehumbert_book_2010](@citep).

**Rayleigh scattering** by gas molecules redirects a fraction of incoming stellar (shortwave) radiation upwards, raising the planet's albedo [stamnes_radiative_2017](@citep). The scattering cross-section scales as $\lambda^{-4}$, so it is key when treating the transport of shorter-wavelength stellar radiation.

## SOCRATES correlated-k

AGNI computes radiation fluxes using a two-stream approximation [zdunkowski_twostream_1980](@citep), which reduces the full angular dependence of the radiation field to separate upwelling and downwelling beams. This approximation is standard in both Earth-system and exoplanet atmosphere models [edwards_studies_1996, lee_testing_2024](@citep). The bolometric net upward flux at a given level follows from integrating the up- and down-welling components over spectral regions (bands).

AGNI uses the correlated-k ($k$-distribution) method for approximating and combing gas absorption [lacis_corrk_1991, edwards_studies_1996](@citep). This approach exploits the statistical distribution of absorption cross-sections within spectral regions rather than computing each individual spectral line. Within a given band, the opacity values $\kappa_\nu$ are sorted in ascending order and grouped into $k$-terms. Overlapping absorption by multiple gas species within a band is treated using the random-overlap or equivalent-extinction method [amundsen_treatment_2017](@citep).

AGNI nominally simulates RT using [SOCRATES](https://proteus-framework.org/SOCRATES/): a suite of numerical codes primarily developed by the UK Met Office [manners_fast_2024, sergeev_socrates_2023](@citep). The model sums over spectral bands to yield bolometric fluxes:
```math
F^{\uparrow,\downarrow} = \sum_{\text{band}} \sum_{\text{g-point}} w_g \, F^{\uparrow,\downarrow}_{\nu_g}
```
where $w_g$ are the $k$-term quadrature weights. SOCRATES implements the two-stream solution of [edwards_studies_1996, edwards_efficient_1996](@citep) with the PIFM/hemispheric-mean flux approximation [zdunkowski_pifm_1985](@citep).

## Opacity sources

The $k$-terms are pre-tabulated from gas line-absorption opacities sourced primarily from the [DACE](https://dace.unige.ch/opacityDatabase/) database [grimm_database_2021](@citep), which draws from the ExoMol and HITEMP molecular line databases [tennyson_exomol_2018, rothman_hitemp_2010](@citep). These are cross-sections integrated over the linelists with a line truncation width of 25 cm$^{-1}$.

Water continuum absorption cross-sections are computed using the MTCKD model [mlawer_mtckd_2012, mlawer_mtckd_2023](@citep). Other continua, including CIA between H₂–H₂ and H₂–He, are taken from the HITRAN collision-induced absorption section [karman_hitran_2019](@citep). Rayleigh scattering and water cloud radiative properties are also included.

The flowchart below outlines how these absorption data are converted into a spectral file used at runtime.

![](fig_spectral_flowchart.svg)

You can find tools for fitting k-terms and processing line absorption data in the [SOCRATES repository](https://github.com/FormingWorlds/SOCRATES) on GitHub.

## Stellar irradiation

A key input to the radiation model is the shortwave downward-directed flux from the star at the top of the atmosphere. The instellation flux is calculated from the stellar bolometric luminosity $L_\star$ and the time-averaged planet–star separation $d$:
```math
F^{\text{ins}} = \frac{L_\star}{4\pi d^2}
```
For a planet on an eccentric orbit with semi-major axis $a$ and eccentricity $e$, the time-averaged separation is $d = a(1 + e^2/2)$ [nicholls_redox_2024](@citep).

The model represents the three-dimensional planet with a single 1D column by choosing a **zenith angle** $\theta_z$ and a stellar flux scale factor $f_s$, following [cronin_choice_2014](@citep). Common choices include the substellar point ($\theta_z = 0$, $f_s = 1$), a global average for a rapidly rotating planet ($\cos\theta_z = 1/\sqrt{3}$, $f_s = 1/4$), and a dayside average for a tidally locked planet ($\cos\theta_z = 1/\sqrt{3}$, $f_s = 1/2$).

## Surface reflectivity

Surface reflectivity can be modelled as a greybody with an albedo from 0 to 1 [Hapke_2012](@citep). Alternatively, the surface can be modelled using empirical reflectance data that varies (spectrally) with wavelength. In the latter case a filepath must be provided via the config. The file can tabulate any one of: spherical reflectance ('r'), hemispherical emissivity ('e'), or single scattering albedo ('w'). These data are compiled on Zenodo [here](https://zenodo.org/communities/proteus_framework/records?q&f=subject%3Asurface_albedos&l=list&p=1&s=10&sort=newest).

## RFM line-by-line

AGNI also includes an interface to the [Reference Forward Model](https://eodg.atm.ox.ac.uk/RFM/) (RFM), which is packaged as a binary since the RFM is closed-source [dudhia_rfm_2017](@citep). This line-by-line interface provides an accurate means to validate and benchmark the correlated-k SOCRATES calculations, and has been used to verify spectral features in synthetic emission spectra [nicholls_convective_2025, nicholls_redox_2024](@citep).

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
