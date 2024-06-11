# AGNI
Radiative-convective solver designed for integration into a coupled atmosphere-interior code.   

AGNI relies on SOCRATES (2311) for calculating radiances. The radiative transfer includes shortwave irradiation from the star, surface emission, gaseous absorption, Rayleigh scattering, parameterised clouds, and continuum absorption. Mixing length theory is used to parametrise convection. Together, energy transport processes allow for an energy-conserving calculation of the atmosphere's temperature profile.      

The model is distributed under a proprietary license. Only once it has been published will the model be distributed under a FOSS license. If you use the model in a publication (once it is open), please cite my paper describing the model.

Consult the [AGNI documentation](https://nichollsh.github.io/AGNI/) for information about the model. 

Contact: `harrison.nicholls [a] physics.ox.ac.uk`   
    
GitHub: https://github.com/nichollsh/AGNI    


## Installation and usage
See the Getting Started page in the [documentation](https://nichollsh.github.io/AGNI/) for information on installing and using the model.

    
## Repository structure 
* `LICENSE.txt`     - License for use and re-use
* `doc/`            - Further documentation
* `out/`            - Output files
* `res/`            - Resources
* `src/`            - Package source code
* `test/`           - Package tests
* `build/`          - Package build scripts
* `.github/`        - GitHub workflows
* `.vscode/`        - Visual Studio Code configuration 
* `tutorials/`      - Notebooks and tutorials
* `agni.jl`         - AGNI executable

