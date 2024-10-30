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
</p>

<p align="center">
  <b>A radiative-convective model for lava planet atmospheres</b>
</p>


## Overview

A numerical model for the atmospheres of hot rocky (exo)planets. AGNI's primary purpose is to simulate the evolving atmospheres of magma ocean planets, while ensuring that radiative-convective equilibrium is maintained throughout the atmosphere.

AGNI models correlated-k radiative transfer including shortwave irradiation from the star, surface emission, gaseous absorption, Rayleigh scattering, parameterised clouds, and CIA. Mixing length theory is used to model convection. Together, energy transport processes allow for an energy-conserving calculation of the atmosphere's temperature profile.

Consult the [AGNI documentation](https://nichollsh.github.io/AGNI/) for information about the model.

Contact: `harrison[dot]nicholls[at]physics.ox.ac.uk`

GitHub: https://github.com/nichollsh/AGNI


## Installation and usage
See the [Getting Started](https://nichollsh.github.io/AGNI/dev/setup/) page in the documentation for information on installing and using the model.

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

This software is available under the GPLv3. Copyright (C) 2024 Harrison Nicholls.
