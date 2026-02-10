<h1 align="center">
    <div>
        <img src="docs/src/assets/logo_title_light.svg#gh-light-mode-only" style="vertical-align: middle;" width="22%"/>
        <img src="docs/src/assets/logo_title_dark.svg#gh-dark-mode-only" style="vertical-align: middle;" width="22%"/>
    </div>
</h1>

<p align="center">
  <b>An open-source model for extreme atmospheres on rocky exoplanets</b>
</p>

<p align="center" style="margin: 0px">
  <a href="https://doi.org/10.1093/mnras/stae2772"><img src="https://img.shields.io/badge/DOI-10.1093%2Fmnras%2Fstae2772-blue"></a>
  <a href="https://joss.theoj.org/papers/380d8e608e9f863b639af76ceebc7131"><img src="https://joss.theoj.org/papers/380d8e608e9f863b639af76ceebc7131/status.svg"></a>
  <a href="https://ascl.net/2508.020"><img src="https://img.shields.io/badge/ascl-2508.020-navy.svg" alt="ascl:2508.020" /></a>
  <a href="https://doi.org/10.5281/zenodo.15386789"><img src="https://img.shields.io/badge/Zenodo-17431569-blue.svg"></a>
</p>

<p align="center" style="margin: -5px">
  <a href="https://github.com/nichollsh/AGNI/actions/workflows/install_and_test.yml"><img src="https://gist.githubusercontent.com/nichollsh/e20f4fa3c7811c75d34005311fef3696/raw/covbadge.svg"></a>
  <a href="https://www.h-nicholls.space/AGNI/"><img src="https://github.com/nichollsh/AGNI/actions/workflows/documentation.yml/badge.svg"></a>
  <a href="LICENSE.txt"><img src="https://img.shields.io/github/license/nichollsh/AGNI?label=License"></a>
</p>


## Overview
AGNI's primary purpose is to simulate the atmospheric temperature-, height-, and compositional-structures of atmospheres overlying magma oceans. It does this while ensuring that radiative-convective equilibrium is maintained throughout the atmosphere. SOCRATES is used to perform correlated-k radiative transfer including: shortwave irradiation from the star, surface emission, line absorption, Rayleigh scattering, parameterised clouds, and collisional absorption. Mixing length theory is used to parametrise convection. AGNI also supports real gas equations of state, self-gravitation, and various spectral surface compositions. Accounting for these energy transport processes permits an energy-conserving calculation of atmospheric structure, obtained using numerical optimisation, which also yields realistic cooling rates for young rocky planets with magma oceans.

Consult the [AGNI documentation](https://www.h-nicholls.space/AGNI/) for information about the model, including usage, functionality, and how to contribute.

Contact: see information on [my website homepage](https://www.h-nicholls.space/).

## Installation and usage
See the [Getting Started](https://www.h-nicholls.space/AGNI/dev/setup/) page in the documentation for information on installing and using the model.

## Citation
If you use AGNI, please cite the following papers:
* Nicholls et al. (2025a)  - [10.1093/mnras/stae2772](https://doi.org/10.1093/mnras/stae2772)
* Nicholls et al. (2025b)  - [10.21105/joss.07726](https://doi.org/10.21105/joss.07726)
* Nicholls et al. (2026) - [2507.02656](https://arxiv.org/abs/2507.02656)

## Example
Below is an animated example of AGNI solving for a temperature-pressure profile, starting from an isothermal state.
<video autoplay loop muted width="100%" src="https://github.com/user-attachments/assets/759d635e-5de4-410c-8d0e-a0a70ae2ea30"></video>


## Repository structure
* `agni.jl`         - The main AGNI executable
* `LICENSE.txt`     - License for use and re-use
* `deps/`           - Package build scripts
* `docs/`           - Documentation source files
* `misc/`           - Miscellaneous files
* `out/`            - Model output files
* `res/`            - Resources (configs, thermodynamic data, etc.)
* `src/`            - Source code
* `test/`           - Tests for the code
* `tutorials/`      - Notebooks and tutorials

This software is available under GPLv3. Copyright (C) 2023-2026 Harrison Nicholls.
