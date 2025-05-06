---
title: 'AGNI: A radiative-convective model for lava planet atmospheres'
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
date: 06 May 2025
bibliography: paper.bib

---

# Summary

It is important that we are able to accurately model the atmospheres of (exo)planets. This is because their atmospheres play a crucial role in setting the environment and conditions on the planet, including how the planet evolves over astronomical timescales. Additionally, it is primarily by observation of their atmospheres that we are able to characterise exoplanets. There is demand for accurate atmosphere models in the context of lava worlds: planets with permanent or fleeting magma oceans.

AGNI[^1] is a Julia program designed to solve for the temperature and radiation environment within the atmospheres of rocky (exo)planets. It leverages a well established FORTRAN code [@edwards_studies_1996; @sergeev_socrates_2023] to calculate radiative fluxes from a given atmospheric temperature structure and composition, which -- alongside representations of convection and other processes -- enables an energy-conserving numerical solution for the atmospheric conditions. In contrast to most other numerical atmosphere codes, AGNI uses an Newton-Raphson optimisation method to obtain this solution, which enables improved performance and scalability. AGNI was specifically developed for use alongside planetary interior models within a coupled framework, although it can also be easily applied to scientific problems as a standalone code.

AGNI can be interacted with as a library; it is used in this way within the [Jupyter notebook tutorials](https://github.com/nichollsh/AGNI/tree/main/tutorials) in the repository. It can also be used as an executable program, where it reads TOML configuration files from the disk and outputs figures and NetCDF data to a specified directory.

[^1]: AGNI can be found on GitHub [here](https://github.com/nichollsh/AGNI).

# Statement of need

It is thought that all rocky planets go through a "magma ocean" stage, where their mantles are fully molten [@lichtenberg_review_2024; @elkins_ranges_2008; @schaefer_magma_2016]. For some planets this may be their permanent state, while for others it is fleeting. This stage allows for rapid exchange of energy and material between their atmospheres and interiors. Since this phase is likely common to many planets -- including Earth and Venus -- it is important that we understand the physical processes involved, and how these processes interact with each other [@schaefer_redox_2017; @maurice_volatile_2024]. Accurate atmosphere models can allow us to connect the theory of these young planets to telescope observations, since it is primarily through their atmospheric properties that were are able to characterise them [@piette_rocky_2023; @perryman_book_2018]. Modelling these young planets involves facing several poorly constrained quantites that govern their atmospheric composition [@guimond_mineralogical_2023; @sossi_redox_2020]. Recently, combined observations and modelling of exoplanet TRAPPIST-1 b enabled inference of its past climatic history [@zieba_obs_2023]. @hu_cancri_2024 were able to characterise the atmosphere of 55 Cancri e by using an advanced proprietary atmosphere model.

Several theoretical studies have modelled the atmospheres and evolution of these young planets, but all have made several simplifying assumptions. @lichtenberg_vertically_2021 coupled a simple atmosphere climate model with an interior evolution model to simulate magma ocean evolution, but they did not account for the possibility of convective stability in their atmospheres. @krissansen_was_2021 and @krissansen_erosion_2024 made similar assumptions. @selsis_cool_2023 used a pure-steam radiative-convective model to investigate convective stability, finding that it is likely to occur, but did not extend their work to coupled scenarios which explore the secular evolution of these planets or gas mixtures; it is unlikely that these atmospheres are exclusively composed of steam [@lichtenberg_vertically_2021]. @piette_rocky_2023 similarly explored potential atmospheres on observable lava planets, but did not consider the physics of atmosphere-interior coupling, and instead chose semi-arbitrary gas compositions. @zilinskas_observability_2023 used the HELIOS atmosphere model to model potential rock-vapour envelopes of sub-Neptune and super-Earth exoplanets, finding that the opacity of various gaseous species (notably SiO) plays a key role in determining the structure and observable properties of exoplanet atmospheres. The demand for realistic atmosphere modelling in the context of secular magma ocean evolution is apparent.

Ensuring sufficient spectral resolution is important in modelling the blanketing effect of these atmospheres, as resolving the opacity (and transparency) of their many gases is known to be key in setting the rate at which these planets can cool by radiation to space [@pierrehumbert_book_2010; @boukrouche_beyond_2021]. It is also important that we are able to run grids of models that explore the range of possible (and as-yet poorly constrained) conditions that these planets could exhibit, which demands efficient modelling given finite computational resources. Magma ocean crystallisation could take up to several Gyr because of continuous tidal forcing and atmospheric blanketing [@walterova_thermal_2020; @driscoll_tidal_2015]. The efficiency afforded by AGNI enables simulations of rocky planets over geological timescales as part of a coupled interior-atmosphere planetary evolution framework. AGNI has so far been used in @hammond_photometric_2024, @nicholls_convection_2025, @nicholls_tidal_2025, and @nicholls_escape_2025

# Comparison with other codes

AGNI is developed with the view of being coupled into the PROTEUS framework[^2] alongside other modules. In addressing the aforementioned problems, it is able to:

* be coupled to a planetary interior modelwith an appropriate surface boundary condition,
* account for atmospheres of diverse gaseous composition with realistic opacities and equations of state,
* solve for an atmospheric temperature structure that conserves energy and allows for convectively stable regions,
* operate with sufficient speed that it may participate in a wide parameter space,

This is possible due to the method by which AGNI numerically obtains a solution for atmospheric temperature structure and energy transport (@nicholls_convection_2025). Our model uses the Newton-Raphson method to conserve energy fluxes through each level of the column to a required tolerence. A typical runtime when applying the model standalone using its command-line interface (Figure 1b) with a poor initial guess of the true temperature profile is 3 minutes. When providing a 'good' guess, such as when AGNI is coupled within the PROTEUS framework (Figure 1a), an atmosphere solution will be obtaind in less than 1 minute. A single radiative transfer calculation takes approximately 30 ms, performed under the correlated-k and two-stream approximations using SOCRATES[^3]: a well-established FORTRAN code developed by the UK Met Office [@edwards_studies_1996; @sergeev_socrates_2023; @amundsen_treatment_2017]. Convection, condensation, and sensible heat transport are also modelled.

HELIOS [@malik_helios_2017] is popular atmosphere model similar to AGNI, but it depends on an Nvidia GPU in order to perform radiative transfer calculations. Whilst this makes each calculation fast, it also means that HELIOS cannot be used on platforms without an Nvidia GPU or with limited resources. GENESIS [@piette_rocky_2023] has been applied to lava planet atmospheres but is closed-source and not publically available. Exo_k [@selsis_cool_2023] is open source and written in pure Python, but is not designed to be coupled with an interior evolution model. These codes have been used to model the atmospheres of static non-evolving planets, so AGNI stands out as being the only open source model currently integrated into a comprehensive interior-atmosphere evolution framework like PROTEUS. No other models of lava planet atmospheres implement a real-gas equation of state.

Coupling with PROTEUS is one primary use-case for AGNI. Our model can also be used standalone (as in @hammond_photometric_2024) through its command-line interface and configuration files, or through Jupyter notebooks (as in the tutorials). Figure 1 below compares the two primary use-cases driving the development of AGNI.

![Visual comparison of the two primary use-cases for AGNI.](application.svg){ width=90% }

[^2]: The PROTEUS framework can be found [here](https://github.com/FormingWorlds/PROTEUS).
[^3]: SOCRATES, packaged with additional tooling, is available [here](https://github.com/nichollsh/SOCRATES).

# Future developments

The landscape of exoplanet science is rapidly evolving. Future updates to AGNI may include:

* Incorporation of aerosols and hazes; supported by SOCRATES in principle, but currently not configurable through AGNI
* An equatorial multi-column mode which parametrises zonal redistribution by atmospheric dynamics
* Accounting for compositional inhibition of convection; e.g. via the Ledoux stability criterion
* Parametrisation of dry convection with a 'full spectrum' model which better represents turbulence in convective fluids

# Documentation

The documentation for AGNI can be read online at [nichollsh.github.io/AGNI](https://nichollsh.github.io/AGNI/).

# References
