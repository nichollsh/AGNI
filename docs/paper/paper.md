---
title: 'AGNI: A radiative-convective model for lava planet atmospheres.'
tags:
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
  - name: Raymond Pierrehumbert
    orcid: 0000-0002-5887-1197
    affiliation: 1
  - name: Tim Lichtenberg
    orcid: 0000-0002-3286-7683
    affiliation: 2
affiliations:
 - name: Department of Physics, University of Oxford, Parks Road, Oxford OX1 3PU, UK
   index: 1
 - name: Kapteyn Astronomical Institute, University of Groningen, P.O. Box 800, 9700 AV Groningen, The Netherlands
   index: 2
date: 20 September 2024
bibliography: paper.bib

---

# Summary

It is important that we are able to accurately model the atmospheres of (exo)planets. This is because their atmospheres play a crucial role in setting the environment and conditions on the planet, including how the planet evolves over astronomical timescales. Additionally, it is primarily by observation of their atmospheres that we are able to characterise exoplanets. There is demand for accurate atmosphere models in the context of lava worlds: planets with permanent or fleeting magma oceans.

AGNI[^1] is a Julia program designed to solve for the temperature and radiation environment within the atmospheres of rocky (exo)planets. It leverages a well established FORTRAN code to calculate radiative fluxes from a given atmospheric temperature structure and composition, which -- alongside representations of convection and other processes -- enables an energy-conserving numerical solution for the atmospheric conditions. In contrast to most other numerical atmosphere codes, AGNI uses an Newton-Raphson optimisation method to obtain this solution, which enables improved performance and scalability. AGNI was specifically developed for use alongside planetary interior models within a coupled framework, although it can also be easily applied to scientific problems as a standalone code.

AGNI can be interacted with as a library; it is used in this way within the [Jupyter notebook tutorials](https://github.com/nichollsh/AGNI/tree/main/tutorials) in the repository. It can also be used as an executable program, where it reads TOML configuration files from the disk and outputs figures and NetCDF data to a specified directory.

[^1]: AGNI can be found on GitHub [here](https://github.com/nichollsh/AGNI).

# Statement of need

It is thought that all rocky planets go through a "magma ocean" stage, where their mantles are fully molten [@lichtenberg_review_2024; @elkins_ranges_2008; @schaefer_magma_2016]. For some planets this may be their permanent state, while for others it is fleeting. This stage allows for rapid exchange of energy and material between their atmospheres and interiors. Since this phase is likely common to many planets -- including Earth and Venus -- it is important that we understand the physical processes involved, and how these processes interact with each other [@schaefer_redox_2017; @maurice_volatile_2024]. Accurate atmosphere models can allow us to connect the theory of these young planets to telescope observations, since it is primarily through their atmospheric properties that were are able to characterise them [@piette_rocky_2023; @perryman_book_2018]. Modelling these young planets involves facing several poorly contrained quantites that govern their atmospheric composition [@guimond_mineralogical_2023; @sossi_redox_2020]. Recently, combined observations and modelling of exoplanet TRAPPIST-1 b enabled inference of its past climatic history [@zieba_obs_2023]. @hu_cancri_2024 were able to characterise the atmosphere of 55 Cancri e by using an advanced proprietary atmosphere model.

Several theoretical studies have modelled the atmospheres and evolution of these young planets, but all have made several simplifying assumptions. @lichtenberg_vertically_2021 coupled a simple atmosphere climate model with radiative transfer with an interior evolution model to simulate magma ocean evolution, but they did not account for the possibility of convective stability in their atmospheres. @krissansen_was_2021 and @krissansen_erosion_2024 made similar assumptions. @selsis_cool_2023 used a pure-steam radiative-convective model to investigate convective stability, finding that it is likely to occur, but did not extend their work to coupled scenarios which explore the secular evolution of these planets or gas mixtures; it is unlikely that these atmospheres are exclusively composed of steam [@lichtenberg_vertically_2021]. @piette_rocky_2023 similarly explored potential atmospheres on observable lava planets, but did not consider the physics of atmosphere-interior coupling, and instead chose semi-arbitrary gas compositions. @zilinskas_observability_2023 used the HELIOS atmosphere model to model potential rock-vapour envelopes of sub-Neptune and super-Earth exoplanets, finding that the opacity of various gaseous species (notably SiO) plays a key role in determining the structure and observable properties of exoplanet atmospheres, although they neglect volatile species and apply the ideal gas equations of state. The demand for realistic atmosphere modelling in the context of secular magma ocean evolution is apparent.

Ensuring sufficient spectral resolution is important in modelling the blanketing effect of these atmospheres, as resolving the opacity (and transparency) of their many gases is known to be key in setting the rate at which these planets can cool by radiation to space [@pierrehumbert_book_2010; @boukrouche_beyond_2021]. It is also important that we are able to run grids of models that explore the range of possible (and as-yet poorly constrained) conditions that these planets could exhibit, which demands efficient modelling given finite computational resources. Magma ocean crystallisation could take up to several Gyr because of continuous tidal forcing and atmospheric blanketing [@walterova_thermal_2020; @driscoll_tidal_2015]. The efficiency afforded by AGNI enables simulations of rocky planets over geological timescales as part of a coupled interior-atmosphere planetary evolution framework. AGNI has so far been used in @hammond_photometric_2024, @nicholls_convection_2025, @nicholls_tidal_2025, and @nicholls_escape_2025

![Caption for example figure.](application.svg){ width=20% }

# Comparison with other codes

AGNI is developed with the view of being coupled into the PROTEUS framework[^2] alongside other physical models. In addressing the aforementioned problems, it is able to:

* be coupled to a planet-interior model (PROTEUS) with an appropriate surface boundary condition,
* account for atmospheres of diverse gaseous composition with realistic opacities and equations of state,
* solve for an atmospheric temperature structure that conserves energy and allows for convectively stable regions,
* operate with sufficient speed that it may participate in a wide parameter space,

This is possible due to the method by which AGNI numerically obtains a solution for atmospheric temperature structure and energy transport (@nicholls_convection_2025). Rather than time-stepping each model level according to radiative heating and applying convective adjustment (e.g. @malik_helios_2017, @selsis_cool_2023, @pierrehumbert_book_2010), AGNI uses the Newton-Raphson method to conserve energy fluxes through the column to a required tolerence, similar to the method applied by @drummond_effects_2016 and @goyal_library_2020.

Radiative transfer is performed under the correlated-k and two-stream approximations using SOCRATES[^3], a well-established FORTRAN code developed by the UK Met Office [@manners_socrates_2024; @sergeev_socrates_2023; @amundsen_treatment_2017; @amundsen_radiation_2014]. Convection, condensation, and sensible heat transport are also modelled.

| Model    | Solution method    | Open source?  | Real gas EoS?  | Interior coupling? | Chemistry   | Supported platforms | Reference|
|:---------|:----------------   |:-------------:|:--------------:|:------------------:|:--------    |:------------------- |:----|
| AGNI     | Nonlinear opt.     | ✅            | ✅             |  ✅                | Eqm.        | Linux & MacOS       | @nicholls_convection_2025 |
| HELIOS   | Time-step          | ✅            | ❌             |  ❌                | Diseqm.     | Nvidia devices only | @malik_helios_2017    |
| GENESIS  | Time-step          | ❌            | ❌             |  ✅                | Diseqm.     | ❓                  | @gandhi_genesis_2017    |
| PCM-LBL  | Time-step          | ✅            | ❌             |  ❌                | None        | Linux & MacOS       | @wordsworth_coupled_2021    |
| Exo_k    | Nondim. time-step  | ✅            | ❌             |  ❌                | None        | Linux & MacOS       | @selsis_cool_2023   |
| PICASO   | Nonlinear opt.     | ✅            | ❌             |  ❌                | Diseqm.     | Linux & MacOS       | @mukherjee_picaso_2023    |

[^2]: The PROTEUS framework can be found [here](https://github.com/FormingWorlds/PROTEUS).
[^3]: SOCRATES, packaged with additional tooling, is available [here](https://github.com/nichollsh/SOCRATES).

# Future developments


# Documentation

The documentation for AGNI can be read online at [nichollsh.github.io/AGNI](https://nichollsh.github.io/AGNI/).

# References
