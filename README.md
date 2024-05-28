# AGNI
Radiative-convective solver designed for integration into a coupled atmosphere-interior code.   

AGNI relies on SOCRATES (2311) for calculating radiances. The radiative transfer includes shortwave irradiation from the star, surface emission, gaseous absorption, Rayleigh scattering, parameterised clouds, and continuum absorption. Mixing length theory is used to parametrise convection. Together, energy transport processes allow for an energy-conserving calculation of the atmosphere's temperature profile.      

The model is distributed under a proprietary license. Only once it has been published will the model be distributed under a FOSS license. If you use the model in a publication (once it is open), please cite my paper describing the model.

Consult the [AGNI Wiki](https://nichollsh.github.io/AGNI/) on GitHub for information about the model. 
    
## Repository structure 
* `README.md`       - This file
* `LICENSE.txt`     - License for use and re-use
* `doc/`            - Further documentation
* `out/`            - Output files
* `res/`            - Resources
* `src/`            - AGNI package source code
* `test/`           - Package tests
* `socrates/`       - Directory containing SOCRATES and associated files (subject to the license therein)
* `.github/`        - GitHub workflows
* `.vscode/`        - Visual Studio Code configuration 
* `agni.jl`         - AGNI executable
* `demo_steamrun.jl`- Script to demonstrate the pure-steam runaway greenhouse effect

## Installation and usage
See the Getting Started page on the AGNI wiki for information on installing and using the model.

## Contributors
* Harrison Nicholls
* Hamish Innes


