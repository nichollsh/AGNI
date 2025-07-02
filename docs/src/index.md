```@raw html
    <img class="display-light-only" src="assets/logo_title_light.svg" width=32% alt="AGNI logo, light mode"/>
    <img class="display-dark-only"  src="assets/logo_title_dark.svg"  width=32% alt="AGNI logo, dark mode"/>
    <p align="center">
        <b>An open-source model for extreme atmospheres on rocky exoplanets</b>
    </p>
```

AGNI's primary purpose is to simulate the atmospheric temperature-, height-, and compositional-structures of atmospheres overlying magma oceans. It does this while ensuring that radiative-convective equilibrium is maintained throughout the atmosphere. SOCRATES is used to perform correlated-k radiative transfer including: shortwave irradiation from the star, surface emission, line absorption, Rayleigh scattering, parameterised clouds, and collisional absorption. Mixing length theory is used to parametrise convection. AGNI also supports real gas equations of state, self-gravitation, and various spectral surface compositions. Accounting for these energy transport processes permits an energy-conserving calculation of atmospheric structure, obtained using numerical optimisation, which also yields realistic cooling rates for young rocky planets with magma oceans.

Pronounced as _ag-nee_. Named after the fire deity of Hinduism.

Follow [Getting started](@ref) for information on installing the code and
obtaining results.

Contact: `harrison[dot]nicholls[at]physics.ox.ac.uk`

GitHub: [https://github.com/nichollsh/AGNI](https://github.com/nichollsh/AGNI)

If you use AGNI, please cite the following papers:
* Nicholls et al., (2025a) - [DOI 10.1093/mnras/stae2772](https://doi.org/10.1093/mnras/stae2772)
* Nicholls et al., (2025b) - [DOI 10.21105/joss.07726](https://doi.org/10.21105/joss.07726)
* Nicholls et al., (2025d) - in review at Nature Astronomy

This software is available under the GPLv3. Copyright (C) 2025 Harrison Nicholls.
