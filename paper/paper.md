---
title: 'AGNI: A radiative-convective model for the atmospheres of rocky planets'
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
affiliations:
 - name: Department of Physics, University of Oxford, Parks Road, Oxford OX1 3PU, UK
   index: 1
date: 13 September 2024
bibliography: paper.bib

---

# Summary

AGNI[^1] is a Julia program designed to solve for the temperature, and radiation environment within the atmospheres of rocky (exo)planets. Given basic physical parameters such as planetary radius, atmospheric composition, and the incoming stellar radiation.

[^1]: AGNI can be found on GitHub [here](https://github.com/nichollsh/AGNI).

# Statement of need

It is thought that all rocky planets go through a "magma ocean" stage, where their mantles are fully molten [@lichtenberg_review_2024; @salvador_review_2023]. For some planets this may be their permanent state, while for others it is fleeting. This stage allows for rapid exchange of energy and material between their atmospheres and interiors. Since this phase is likely common to many planets -- including Earth and Venus -- it is important that we understand the physical processes involved, and how these processes interact with each other [@schaefer_redox_2017; @maurice_volatile_2024]. Accurate atmosphere models can allow us to connect the theory of these young planets to telescope observations, since it is primarily through their atmospheric properties that were are able to characterise them [@piette_rocky_2023; @perryman_book_2018]. Modelling these young planets also involves facing several poorly contrained (and unconstrained) quantites which govern their atmospheric composition and dynamics [@guimond_mineralogical_2023; @sossi_redox_2020].

Several theoretical studies have modelled the atmospheres and evolution of these young planets, but have all made several simplifying assumptions. @lichtenberg_vertically_2021 coupled a simple atmosphere model to real-gas radiative transfer and interior evolution codes to simulate magma ocean evolution, but did not account for the possibility of convective stability in their atmospheres. @krissansen_was_2021 made similar assumptions. @selsis_cool_2023 used a pure-steam radiative-convective model to investigate convective stability, finding that it is likely to occur, but did not extend their work to coupled scenarios which explore the secular evolution of these planets or gas mixtures. @zilinskas_observability_2023 used a multi-gas model to simulate the atmospheres of Sub-Neptune and Super-Earth exoplanets an equilibrium, finding that the opacity of various gaseous species (notably SiO) plays a key role in determining the structure and observable properties of exoplanet atmospheres. The demand for realistic atmosphere modelling in the context of secular magma ocean evolution is apparent.

Ensuring sufficient spectral resolution is important in modelling the blanketing effect of these atmospheres, as resolving the opacity (and transparency) of their many gases is known to be key in setting the rate at which these planets can cool by radiation to space [@pierrehumbert_book_2010; @boukrouche_beyond_2021]. It is also important that we are able to run grids of models which explore the range of possible (and as-yet poorly constrained) conditions that these planets could exhibit, which demands efficient modelling given finite computational resources. Performance is paramount.

HELIOS [@malik_helios_2017] is a hydrostatic atmosphere model written in Python and CUDA. HELIOS assumes that the "interior temperature" $T_{\text{int}}$ of a model planet is a known quantity, but coupled time-evolved evolution with an interior model requires that this quantity be an output variable, instead requiring a fixed surface temperature (or something equivalent). It is therefore not possible to apply HELIOS to this problem. The same also applies to Exo_k [@selsis_cool_2023].

# Comparison with other codes

AGNI is a new radiative-convective atmosphere model developed with the view of being coupled into the PROTEUS framework[^2]. In solving the aforementioned problems, it is able to:
* be coupled to an interior model with an appropriate surface boundary condition,
* account for atmospheres of mixed gaseous composition with realistic opacities,
* solve for a temperature structure which conserves energy and allows for convective stability,
* operate with sufficient speed that it may be participate in a wide parameter space.

This is possible due to the method by which AGNI numerically obtains a solution for atmospheric temperature structure and energy transport. Rather than time-stepping each model level according to radiative heating and applying convective adjustment (cf. HELIOS, Exo_k, and various global circulation models), AGNI uses the Newton-Raphson method to find the state which conserves energy fluxes through the column to a required tolerence. This is similar to the method applied by @drummond_effects_2016 and @goyal_library_2020, although with several optimisations. This allows the model to take tens or hundereds of iterations to obtain a solution, in comparison to thousands or tens of thousands with a time-stepping scheme.

In AGNI, radiative transfer is performed under the correlated-k and two-stream approximations with up to 4096 spectral bands distributed between 1 and 35000 cm$^{-1}$ [@lacis_corrk_1991; @stamnes_radiative_2017]. This is done by using SOCRATES[^3], a well established FORTRAN code developed by the UK Met Office [@manners_socrates_2024; @amundsen_treatment_2017; amundsen_radiation_2014]. Convection is parameterised using mixing length theory [@joyce_mlt_2023]. Condensation and sensible heat transport are also modelled.

Alongside the problem of interior-atmosphere coupling, it is useful to be able to leverage SOCRATES through an interactive interface -- in this case thanks to Julia. This was applied in the recent paper by @hammond_photometric_2024.

[^2]: The PROTEUS framework can be found on GitHub [here](https://github.com/FormingWorlds/PROTEUS).
[^3]: Despite SOCRATES being developed under the BSD 3-Clause license, its [main development repository](https://code.metoffice.gov.uk/trac/socrates) is not publically accessible. The code has been re-hosted on GitHub, with additional tools, under the same open-source license [here](https://github.com/nichollsh/SOCRATES).


# Acknowledgements

Harrison Nicholls is supported by the Clarendon Fund and the MT Scholarship Trust.
The author is grateful for the continuing hard work of the Julia developers and those of its many libraries [@julialang].
Additionally, the author thanks James Manners, and others at the Met Office, who have continued to support and develop SOCRATES over the last three decades.

# References
