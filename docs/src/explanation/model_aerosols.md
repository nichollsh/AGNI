# Aerosols and clouds

AGNI represents the radiative effects of condensed particles (aerosols, hazes, and clouds) within the two-stream radiative transfer calculations performed by SOCRATES [edwards_studies_1996](@citep). Particles extinguish radiation by absorption and scattering, and so modify both the stellar (shortwave) and thermal (longwave) fluxes. This page describes how the optical properties of these particles are obtained.

## Representation in the radiative transfer

Each aerosol species is described by a mass mixing ratio profile $q(p)$ (kg kg$^{-1}$), and by three band-averaged optical properties which are stored in the SOCRATES spectral file:

- the mass absorption coefficient $\bar{k}_\mathrm{abs}$
- mass scattering coefficient $\bar{k}_\mathrm{sca}$ (m$^2$ kg$^{-1}$ of particle material)
- and the asymmetry parameter $\bar{g}$.

Within a layer of air mass per unit area $\Delta m$, an aerosol contributes an extinction optical depth
```math
\Delta\tau_\mathrm{ext} = q \, (\bar{k}_\mathrm{abs} + \bar{k}_\mathrm{sca}) \, \Delta m ,
```
with single scattering albedo $\omega = \bar{k}_\mathrm{sca}/(\bar{k}_\mathrm{abs}+\bar{k}_\mathrm{sca})$.

The asymmetry parameter measures how much light is scattered forward versus backward by particles, ranging from -1 (all backward) to 1 (all forward), with 0 indicating equal scattering in all directions. The phase function is represented by its first moment, $\bar{g}$.

The mixing ratio of each aerosol is either fixed (config key `mmr`) or tied to the condensate yield of a condensable gas (config key `species`), in which case it is updated as the atmospheric structure evolves.

Optical properties are obtained using one of two methods, chosen per species:

- `mon`, pre-computed monochromatic scattering data supplied with SOCRATES (e.g. soot, sulphate, dust, sea salt). These are averaged over the bands of the spectral file using the SOCRATES tool `scatter_average` at runtime.
- `mie`, properties calculated at runtime from the complex refractive index of the particle material, using Mie theory. This allows arbitrary materials, such as the silicate and iron oxide clouds expected in the atmospheres of lava planets.

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
where $a_j$ and $b_j$ are the Mie coefficients. These are evaluated using the BHMIE algorithm of [bohren_absorption_1983](@citet), in which the logarithmic derivative of the Riccati–Bessel function is computed by downward recurrence. The series is truncated after $N = x + 4.05x^{1/3} + 2$ terms [wiscombe_improved_1980](@citep). For very small particles the Rayleigh limit is used, $Q_\mathrm{sca} = \tfrac{8}{3} x^4 |L|^2$ and $Q_\mathrm{abs} = 4x\,\mathrm{Im}(L)$ with $L = (m^2-1)/(m^2+2)$.

Bessel functions are a class of special functions that commonly appear in problems involving wave motion, heat conduction, and other physical phenomena with circular or cylindrical symmetry. Riccati-Bessel functions are similar to spherical Bessel functions [du_mie_2004](@citep).

The algorithm implemented in AGNI is almost identical to the classic FORTRAN scheme for Mie
calculations, which you can view here: https://www.astro.princeton.edu/~draine/code/bhmie.f

Mie theory neglects the effects of non-spherical shapes, porosity, aggregate structure, and mixed composition. Most laboratory data were measured on solids at room temperature. Each species has a single size distribution which applies at all altitudes. The particle size does not respond to the local microphysics. The phase function is represented only by its asymmetry parameter.

The implementation is verified against an independent Mie code and against the Mie code within SOCRATES (see the [Testing suite](@ref)).

### Size distribution

Particle radii follow a log-normal number distribution,
```math
\frac{dN}{d\ln r} = \frac{N}{\sqrt{2\pi}\,\ln\sigma_g} \exp\left[ -\frac{(\ln r - \ln r_g)^2}{2 \ln^2 \sigma_g} \right],
```
which is configured by its effective (area-weighted mean) radius $r_\mathrm{eff} = \langle r^3 \rangle / \langle r^2 \rangle$ [hansen_light_1974](@citep) and geometric standard deviation $\sigma_g$. The geometric mean radius is then $r_g = r_\mathrm{eff} \exp(-\tfrac{5}{2}\ln^2\sigma_g)$. Setting $\sigma_g = 1$ gives a monodisperse population.

Averages over the distribution are computed using Gauss–Hermite quadrature in $\ln r$.
This is a numerical integration method (a kind of Gaussian quadrature) used to approximate integrals where the integrand is weighted by a Gaussian function.

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
Within each SOCRATES correlated-k spectral band, the optical properties are evaluated with logarithmically spaced wavelengths within the band. The properties are weighted using the stellar spectrum, which is cumulatively integrated over log-wavelength space.

The resultant properties are written into the runtime copy of the spectral file as dry aerosols (block type 11), using aerosol type numbers above those reserved by SOCRATES. This happens at the same time as inserting Rayleigh scattering data and the stellar spectrum itself. Mie data are also written to the output NetCDF file (variables `aer_kabs`, `aer_ksca`, `aer_asym`), and can be plotted with the `plots.aerosol_optics` option.

## Refractive index data

Refractive indices are largely adopted from the compilation distributed with the POSEIDON retrieval code [mullens_implementation_2024](@citep). It collates laboratory measurements from several sources, including the databases of [wakeford_transmission_2015](@citet) and [kitzmann_optical_2018](@citet). See [Obtaining input data](@ref). The bulk density of each material is required for these calculations. These are set in the source code (`density.condensate_rho`) directly. Only materials which have both a refractive index file and a sourced density can be used for Mie calculations. The refractive index is extrapolated with constant values outside the available WL ranges.

## Water clouds

Water clouds are treated separately from aerosols, using the droplet parametrisation stored in the spectral file (block type 10). This parametrises the optical properties of water droplets with a Padé fit [edwards_studies_1996](@citep). The cloud water content is set by the condensate yield of water.

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
