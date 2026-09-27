# Aerosols and clouds

AGNI represents the radiative effects of condensed particles (aerosols, hazes, and clouds) within the two-stream radiative transfer calculations performed by SOCRATES [edwards_studies_1996](@citep). Particles extinguish radiation by absorption and scattering, and so modify both the stellar (shortwave) and thermal (longwave) fluxes. This page describes how the optical properties of these particles are obtained.

## Representation in the radiative transfer

Each aerosol species is described by a mass mixing ratio profile $q(p)$ [kg kg$^{-1}$], and by three band-averaged optical properties which are stored in the SOCRATES spectral file: the mass absorption coefficient $\bar{k}_\mathrm{abs}$ and mass scattering coefficient $\bar{k}_\mathrm{sca}$ [m$^2$ kg$^{-1}$ of particle material], and the asymmetry parameter $\bar{g}$. Within a layer of air mass per unit area $\Delta m$, an aerosol contributes an extinction optical depth
```math
\Delta\tau_\mathrm{ext} = q \, (\bar{k}_\mathrm{abs} + \bar{k}_\mathrm{sca}) \, \Delta m ,
```
with single scattering albedo $\omega = \bar{k}_\mathrm{sca}/(\bar{k}_\mathrm{abs}+\bar{k}_\mathrm{sca})$. The phase function is represented by its first moment, $\bar{g}$.

The mixing ratio of each aerosol is either fixed (config key `mmr`) or tied to the condensate yield of a condensable gas (config key `species`), in which case it is updated as the atmospheric structure evolves.

Optical properties are obtained using one of two methods, chosen per species:

- `mon` — pre-computed monochromatic scattering data supplied with SOCRATES (e.g. soot, sulphate, dust, sea salt). These are averaged over the bands of the spectral file using the SOCRATES tool `scatter_average` at runtime.
- `mie` — properties calculated at runtime from the complex refractive index of the particle material, using Mie theory. This allows arbitrary materials, such as the silicate and iron oxide clouds expected in the atmospheres of lava planets.

## Mie theory

For the `mie` method, particles are assumed to be homogeneous spheres. For a sphere of radius $r$ at wavelength $\lambda$, with size parameter $x = 2\pi r/\lambda$ and complex refractive index $m = n + ik$, the extinction and scattering efficiencies are [bohren_absorption_1983](@citep)
```math
Q_\mathrm{ext} = \frac{2}{x^2} \sum_{j=1}^{N} (2j+1) \, \mathrm{Re}(a_j + b_j), \qquad
Q_\mathrm{sca} = \frac{2}{x^2} \sum_{j=1}^{N} (2j+1) \left(|a_j|^2 + |b_j|^2\right),
```
and the asymmetry parameter is
```math
g = \frac{4}{x^2 Q_\mathrm{sca}} \sum_{j=1}^{N} \left[ \frac{j(j+2)}{j+1} \mathrm{Re}(a_j a_{j+1}^* + b_j b_{j+1}^*) + \frac{2j+1}{j(j+1)} \mathrm{Re}(a_j b_j^*) \right],
```
where $a_j$ and $b_j$ are the Mie coefficients. These are evaluated using the BHMIE algorithm of [bohren_absorption_1983](@citet), in which the logarithmic derivative of the Riccati–Bessel function is computed by downward recurrence. This is numerically stable for strongly absorbing particles. The series is truncated after $N = x + 4.05x^{1/3} + 2$ terms [wiscombe_improved_1980](@citep). For very small particles ($x < 10^{-3}$) the Rayleigh limit is used, $Q_\mathrm{sca} = \tfrac{8}{3} x^4 |L|^2$ and $Q_\mathrm{abs} = 4x\,\mathrm{Im}(L)$ with $L = (m^2-1)/(m^2+2)$.

The implementation is verified against an independent Mie code and against the Mie code within SOCRATES (see the [Testing suite](@ref)).

### Size distribution

Particle radii follow a log-normal number distribution,
```math
\frac{dN}{d\ln r} = \frac{N}{\sqrt{2\pi}\,\ln\sigma_g} \exp\left[ -\frac{(\ln r - \ln r_g)^2}{2 \ln^2 \sigma_g} \right],
```
which is configured by its effective (area-weighted mean) radius $r_\mathrm{eff} = \langle r^3 \rangle / \langle r^2 \rangle$ [hansen_light_1974](@citep) and geometric standard deviation $\sigma_g$. The geometric mean radius is then $r_g = r_\mathrm{eff} \exp(-\tfrac{5}{2}\ln^2\sigma_g)$. Setting $\sigma_g = 1$ gives a monodisperse population. Averages over the distribution are computed using Gauss–Hermite quadrature in $\ln r$.

The mass coefficients follow from the distribution-averaged cross-sections and particle volume,
```math
k_\mathrm{abs,sca}(\lambda) = \frac{\langle \pi r^2 Q_\mathrm{abs,sca} \rangle}{\rho_p \, \langle \tfrac{4}{3}\pi r^3 \rangle},
```
where $\rho_p$ is the bulk density of the particle material, and the distribution-averaged asymmetry parameter is weighted by the scattering cross-section.

### Band averaging

The spectral properties are averaged over each band of the spectral file, weighted by the stellar spectrum $F_\star(\lambda)$ ('thin' averaging, as in SOCRATES):
```math
\bar{k} = \frac{\int_\mathrm{band} k(\lambda) F_\star(\lambda) \, d\lambda}{\int_\mathrm{band} F_\star(\lambda) \, d\lambda}, \qquad
\bar{g} = \frac{\int_\mathrm{band} g(\lambda) k_\mathrm{sca}(\lambda) F_\star(\lambda) \, d\lambda}{\int_\mathrm{band} k_\mathrm{sca}(\lambda) F_\star(\lambda) \, d\lambda} .
```
Within each band, the optical properties are evaluated at a minimum of 16 logarithmically spaced wavelengths plus the tabulated wavelengths of the refractive index data, so that narrow features (such as the Si–O stretching resonance of silicates near 9 μm) are resolved. The stellar spectrum is integrated at its native resolution, using power-law interpolation between tabulated points. Where the stellar spectrum contains no flux within a band, uniform weighting is used.

The resultant properties are written into the runtime copy of the spectral file as dry aerosols (block type 11), using aerosol type numbers above those reserved by SOCRATES. They are also written to the output NetCDF file (variables `aer_kabs`, `aer_ksca`, `aer_asym`), and can be plotted with the `plots.aerosol_optics` option.

## Refractive index data

Refractive indices are taken from the compilation distributed with the POSEIDON retrieval code [mullens_implementation_2024](@citep). It collates laboratory measurements from several sources, including the databases of [wakeford_transmission_2015](@citet) and [kitzmann_optical_2018](@citet). Each file is identified by a material name (e.g. `SiO2_amorph`, `FeO`, `MgSiO3`), which is passed to AGNI through the `nk_file` config key. The files can be downloaded using `get_data.sh` (see [Obtaining input data](@ref)). The provenance of each dataset is recorded in the header of its file and in the `Aerosol-Database-Readme.txt` file distributed alongside the data.

The bulk density of each material is tabulated in the source code (`density.condensate_rho`). Only materials which have both a refractive index file and a sourced density can be used. Room-temperature values are used, mostly taken from the unit-cell (X-ray) densities in the Handbook of Mineralogy for crystalline phases. Materials for which a density could not be sourced for the phase of the optical data (for example porous amorphous Al₂O₃ and sol-gel amorphous Mg₂SiO₄) are not currently available.

Outside of the tabulated wavelength range of a material, the refractive index is held at its value at the nearest tabulated wavelength, and a warning is issued. Data rows with unphysical values ($n \le 0$ or $k < 0$) are skipped with a warning.

## Water clouds

Water clouds are treated separately from aerosols, using the droplet parametrisation stored in the spectral file (block type 10). This parametrises the optical properties of water droplets as Padé approximants in the droplet effective radius [edwards_studies_1996](@citep). The cloud water content is set by the condensate yield of H₂O. The spectral file must contain droplet data for clouds to be enabled.

## Assumptions and limitations

- **Spherical, homogeneous particles.** Mie theory neglects the effects of non-spherical shapes, porosity, aggregate structure, and mixed composition.
- **Refractive indices at room temperature.** Most laboratory data were measured on solids at room temperature. Condensates in the atmospheres of lava planets may be hot or molten, and their optical constants may differ substantially. There are few measurements at high temperature.
- **Uniform size distribution.** Each species has a single size distribution which applies at all altitudes. The particle size does not respond to the local microphysics.
- **Dry particles.** Hygroscopic growth is not included for `mie` aerosols.
- **Two-stream scattering.** The phase function is represented only by its asymmetry parameter. Large particles scatter strongly forwards, and this is not captured beyond what the asymmetry parameter describes.
- **Stellar weighting of bands.** The same stellar weighting is applied to all bands, including those dominated by thermal emission. In narrow bands this is a small effect, but it may matter in wide bands where the optical properties vary strongly.
- **Extrapolation.** Bands which lie outside of the tabulated refractive index data use extrapolated values. The fraction of each band affected is reported at runtime.

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
