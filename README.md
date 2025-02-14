<h1 align="center">
    <div>
        <img src="docs/src/assets/logo_title.svg" style="vertical-align: middle;" width="22%"/>
    </div>
</h1>

<p align="center">
  <a href="https://github.com/nichollsh/AGNI/actions/workflows/install_and_test.yml">
    <img src="https://github.com/nichollsh/AGNI/actions/workflows/install_and_test.yml/badge.svg">
  </a>
  <a href="https://nichollsh.github.io/AGNI/dev/">
    <img src="https://github.com/nichollsh/AGNI/actions/workflows/documentation.yml/badge.svg">
  </a>
  <a href="LICENSE.txt">
    <img src="https://img.shields.io/github/license/nichollsh/AGNI?label=License">
  </a>
  <a href="https://doi.org/10.1093/mnras/stae2772">
    <img src="https://img.shields.io/badge/DOI-10.1093%2Fmnras%2Fstae2772-blue">
  </a>
</p>

<p align="center">
  <b>A radiative-convective atmosphere model for lava planets</b>
</p>


## Overview
A numerical model for the atmospheres of hot rocky (exo)planets. The model's primary purpose is to simulate the evolving atmospheres of magma ocean (lava) planets, while ensuring that radiative-convective equilibrium is maintained throughout the atmosphere.

AGNI models correlated-k radiative transfer including shortwave irradiation from the star, surface emission, gaseous absorption, Rayleigh scattering, parameterised clouds, and CIA. Mixing length theory is used to model convection. Together, energy transport processes allow for an energy-conserving calculation of the atmosphere's temperature profile. It also supports real gas equations of state.

Consult the [AGNI documentation](https://nichollsh.github.io/AGNI/) for information about the model.

Contact: `harrison[dot]nicholls[at]physics.ox.ac.uk`

## Installation and usage
See the [Getting Started](https://nichollsh.github.io/AGNI/dev/setup/) page in the documentation for information on installing and using the model.

## Citation
If you use AGNI, please cite the following papers:
* Nicholls et al., (2024b) - https://doi.org/10.1093/mnras/stae2772
* Nicholls et al., (2025) - in review at JOSS

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

This software is available under GPLv3. Copyright (C) 2025 Harrison Nicholls.
