# Related codes

AGNI is a standalone piece of software designed to model planetary atmospheres. However, its development has always been in the context of the '[PROTEUS](https://github.com/FormingWorlds/PROTEUS)' planetary evolution framework, developed by the [FormingWorlds Lab](https://www.formingworlds.space/). PROTEUS simulates the coupled evolution of the atmospheres and interiors of rocky planets and exoplanets by connecting various codes together as 'modules'. In this manner, AGNI may be used as the atmosphere-climate module within PROTEUS, allowing and evolutionary solution to planetary evolution with realistic atmosphere modelling. PROTEUS accesses AGNI using the `juliacall` package from [PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl). There is more information on this in the '[Accessing-AGNI-from-Python](@ref)' page of this documentation.

Other components of the wider PROTEUS ecosystem may be found in the table below:

|  Code    | Description | Website |
|:---------|:------------|:--------|
| PROTEUS  | Coupled modelling framework  | [proteus-framework.org/proteus](https://proteus-framework.org/proteus) |
| SOCRATES | Convective atmosphere | [github.com/FormingWorlds/SOCRATES](https://github.com/FormingWorlds/SOCRATES) |
| JANUS    | Convective atmosphere | [proteus-framework.org/janus](https://proteus-framework.org/janus) |
| MORS     | Stellar evolution     | [proteus-framework.org/mors](https://proteus-framework.org/mors) |
| ZEPHYRUS | Hydrodynamic escape   | [proteus-framework.org/zephyrus](https://proteus-framework.org/zephyrus) |
| CALLIOPE | Volatile outgassing   | [proteus-framework.org/calliope](https://proteus-framework.org/calliope) |
| LovePy   | Mantle tidal heating  | [github.com/nichollsh/lovepy](https://github.com/nichollsh/lovepy) |
| Obliqua  | Mantle tidal heating  | [proteus-framework.org/obliqua](https://proteus-framework.org/obliqua) |
| VULCAN   | Chemical kinetics     | [proteus-framework.org/vulcan](https://proteus-framework.org/vulcan) |
| Zalmoxis | Interior structure    | [proteus-framework.org/zalmoxis](https://proteus-framework.org/zalmoxis) |
| Aragog   | Interior dynamics     | [proteus-framework.org/aragog](https://proteus-framework.org/spider) |
| SPIDER   | Interior dynamics     | [proteus-framework.org/spider](https://proteus-framework.org/spider) |

