<h1 align="center">
    <div>
        <img src="docs/src/assets/logo.svg" style="vertical-align: middle;" width="100px"/>
        <span style="vertical-align: middle;">AGNI</span>
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
  <b>A radiative-convective model for lava planet atmospheres.</b>
</p>


## Overview

AGNI relies on [SOCRATES](https://github.com/nichollsh/SOCRATES) for calculating radiances. The radiative transfer includes shortwave irradiation from the star, surface emission, gaseous absorption, Rayleigh scattering, parameterised clouds, and continuum absorption. Mixing length theory is used to parametrise convection. Together, energy transport processes allow for an energy-conserving calculation of the atmosphere's temperature profile.

Consult the [AGNI documentation](https://nichollsh.github.io/AGNI/) for information about the model.

Contact: `harrison[dot]nicholls[at]physics.ox.ac.uk`

GitHub: https://github.com/nichollsh/AGNI


## Installation and usage
See the Getting Started page in the [documentation](https://nichollsh.github.io/AGNI/) for information on installing and using the model.

## Repository structure
* `agni.jl`         - The main AGNI executable
* `LICENSE.txt`     - License for use and re-use
* `get_fastchem.sh` - Download and setup FastChem
* `get_socrates.sh` - Download and setup SOCRATES
* `get_data.sh`     - Download input data files
* `deps/`           - Package build scripts
* `docs/`           - Documentation source files
* `misc/`           - Miscellaneous files
* `out/`            - Model output files
* `res/`            - Resources (configs, thermodynamic data, etc.)
* `src/`            - Package source code
* `test/`           - Package tests
* `tutorials/`      - Notebooks and tutorials

Copyright (C) 2024 Harrison Nicholls
