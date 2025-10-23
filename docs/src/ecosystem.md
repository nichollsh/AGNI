# Related codes

AGNI is a standalone piece of software designed to model planetary atmospheres. However, its development has always been in the context of the '[PROTEUS](https://github.com/FormingWorlds/PROTEUS)' planetary evolution framework, developed by the [FormingWorlds Lab](https://www.formingworlds.space/). PROTEUS simulates the coupled evolution of the atmospheres and interiors of rocky planets and exoplanets by connecting various codes together as 'modules'. In this manner, AGNI may be used as the atmosphere-climate module within PROTEUS, allowing and evolutionary solution to planetary evolution with realistic atmosphere modelling. PROTEUS accesses AGNI using the `juliacall` package from [PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl). There is more information on this in the [Using the model](@ref) page of this documentation.

Other components of the wider PROTEUS ecosystem may be found in the table below:

|  Code    | Description | Website |
|----------|-------------| ------- |
| PROTEUS  | Model coupling    | https://fwl-proteus.readthedocs.io/ |
| JANUS    | Convective atmosphere | https://fwl-proteus.readthedocs.io/projects/janus/en/latest/ |
| MORS     | Stellar evolution | https://fwl-proteus.readthedocs.io/projects/mors/en/latest/ |
| ZEPHYRUS | Hydrodynamic escape | https://github.com/FormingWorlds/ZEPHYRUS/ |
| CALLIOPE | Volatile outgassing | https://fwl-proteus.readthedocs.io/projects/calliope/en/latest/ |
| LovePy   | Mantle tidal heating | https://github.com/nichollsh/lovepy |
| VULCAN   | Chemical kinetics | https://github.com/FormingWorlds/VULCAN |
| Zalmoxis | Interior structure | https://zalmoxis.readthedocs.io/en/latest/ |
| Aragog   | Interior dynamics | https://github.com/FormingWorlds/aragog |
| SPIDER   | Interior dynamics | https://github.com/djbower/spider |
