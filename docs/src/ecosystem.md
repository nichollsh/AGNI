# Related codes

AGNI is a standalone piece of software designed to model planetary atmospheres. However, its development has always been in the context of the [**PROTEUS**](https://github.com/FormingWorlds/PROTEUS) planetary evolution framework.

```@raw html
<center>
<object type="image/png" width="45%" data="https://cdn.jsdelivr.net/gh/FormingWorlds/PROTEUS@main/docs/assets/PROTEUS_white.png"></object>
</center>
```

PROTEUS (/ˈproʊtiəs, PROH-tee-əs) is a modular Python framework that simulates the coupled evolution of the atmospheres and interiors of rocky planets and exoplanets. PROTEUS achieves this by connecting various codes together, which are referred to as 'modules'.  In this manner, AGNI is applied as the primary atmosphere-climate module within PROTEUS, allowing a solution for planetary evolution with realistic atmosphere modelling. It is primarily developed by the [FormingWorlds Lab](https://www.formingworlds.space/) (University of Groningen), and the [Planetary Chemistry](https://www.shorttle.com/people/) group (University of Cambridge).

```@raw html
<center>
<object type="image/svg+xml" width="90%" data="https://cdn.jsdelivr.net/gh/FormingWorlds/PROTEUS@main/docs/assets/proteus_modules_schematic.svg"></object>
</center>
```

PROTEUS' coupling to AGNI is enabled by  the `juliacall` package from [PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl). There is more information on this in the [Accessing AGNI from Python](@ref) page.

Other components of the wider AGNI-adjacent software ecosystem may be found in the table below:

|  Code    | Description | Website |
|:---------|:------------|:--------|
| PyExoCross   | Gas cross-section integrator   | [github.com/nichollsh/PyExoCross](https://github.com/nichollsh/PyExoCross) |
| ThermoTools   | Thermodynamic properties   | [github.com/nichollsh/ThermoTools](https://github.com/nichollsh/ThermoTools) |
| LovePy   | Mantle tidal heating         | [github.com/nichollsh/lovepy](https://github.com/nichollsh/lovepy) |
