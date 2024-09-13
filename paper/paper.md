---
title: 'AGNI: A radiative-convective model for the atmospheres of rocky planets'
tags:
  - Julia
  - Fortran
  - astronomy
  - physics,
  - radiative transfer
  - exoplanets
  - convection
  - radiation
  - planets
  - atmospheres
authors:
  - name: Harrison Nicholls
    orcid: 0000-0002-8368-4641
    affiliation: 1
affiliations:
 - name: Department of Physics, University of Oxford, Parks Road, Oxford OX1 3PU, UK
   index: 1
date: 13 September 2024
bibliography: paper.bib

---

# Summary

It is thought that all rocky planets go through a "magma ocean" stage, where their mantles are fully molten. This is because they are initially very hot. This stage allows for rapid exchange of energy and material between their atmospheres and interiors. Since this phase is likely common to many planets -- including Earth -- it is important that we understand the physical processes involved, and how these processes interact with each other. It is thought that atmospheric blanketing is key to setting the behaviour in this stage: where the opacity of the atmosphere prevents the planet from cooling down, by absorbing energy emitted from the surface and re-radiating it back downwards. Accurate atmosphere models also allow us to connect theoretical models of these young planets to telescope observations of exoplanets, since it is primarily through their atmospheric properties that were are able to characterise these planets.

AGNI[^1] is a Julia program designed to solve for the temperature, composition, and radiation environment within the atmospheres of these planets, with the view of being coupled to an interior model through the PROTEUS framework[^2].

[^1]: AGNI can be found on GitHub [here](https://github.com/nichollsh/AGNI).
[^2]: The PROTEUS framework can be found on GitHub [here](https://github.com/FormingWorlds/PROTEUS).

# Statement of need

Hi

# Acknowledgements

Harrison Nicholls is supported by the Clarendon Fund and the MT Scholarship Trust.

# References
