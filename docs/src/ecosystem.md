# Related codes

AGNI is a standalone piece of software designed to model planetary atmospheres. However, its development has always been in the context of the '[PROTEUS](https://github.com/FormingWorlds/PROTEUS)' planetary evolution framework, developed by the [FormingWorlds Lab](https://www.formingworlds.space/). PROTEUS simulates the coupled evolution of the atmospheres and interiors of rocky planets and exoplanets by connecting various codes together as 'modules'. In this manner, AGNI may be used as the atmosphere-climate module within PROTEUS, allowing and evolutionary solution to planetary evolution with realistic atmosphere modelling. PROTEUS accesses AGNI using the `juliacall` package from [PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl). There is more information on this in the '[Accessing-AGNI-from-Python](@ref)' page of this documentation.

Other components of the wider PROTEUS ecosystem may be found in the table below:

|  Code    | Description | Website |
|:---------|:------------|:--------|
| PROTEUS  | Model coupling    | [fwl-proteus.readthedocs.io](https://fwl-proteus.readthedocs.io/) |
| JANUS    | Convective atmosphere | [fwl-janus.readthedocs.io](https://fwl-janus.readthedocs.io) |
| MORS     | Stellar evolution | [fwl-mors.readthedocs.io](https://fwl-mors.readthedocs.io/) |
| ZEPHYRUS | Hydrodynamic escape | [fwl-zephyrus.readthedocs.io](https://fwl-zephyrus.readthedocs.io/) |
| CALLIOPE | Volatile outgassing | [fwl-calliope.readthedocs.io](https://fwl-calliope.readthedocs.io/) |
| LovePy   | Mantle tidal heating | [github.com/nichollsh/lovepy](https://github.com/nichollsh/lovepy) |
| VULCAN   | Chemical kinetics | [github.com/FormingWorlds/VULCAN](https://github.com/FormingWorlds/VULCAN) |
| Zalmoxis | Interior structure | [zalmoxis.readthedocs.io](https://zalmoxis.readthedocs.io/) |
| Aragog   | Interior dynamics | [github.com/FormingWorlds/aragog](https://github.com/FormingWorlds/aragog) |
| SPIDER   | Interior dynamics | [github.com/djbower/spider](https://github.com/djbower/spider) |


