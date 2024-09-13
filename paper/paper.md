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

AGNI[^1] is a Julia program designed to solve for the temperature, and radiation environment within the atmospheres of these planets. The user is required to specify basic physical parameters such as planetary radius, atmospheric composition, and the incoming stellar radiation. In return, AGNI calculates the temperature p

[^1]: AGNI can be found on GitHub [here](https://github.com/nichollsh/AGNI).

# Statement of need

It is thought that all rocky planets go through a "magma ocean" stage, where their mantles are fully molten. This is because they are initially very hot. This stage allows for rapid exchange of energy and material between their atmospheres and interiors. Since this phase is likely common to many planets -- including Earth -- it is important that we understand the physical processes involved, and how these processes interact with each other. It is thought that atmospheric blanketing and the greenhouse effect are key to setting their behaviour in this stage: where the opacity of the atmosphere can prevent the planet from efficiently cooling. Accurate atmosphere models can also allow us to connect the theory of these young planets to telescope observations, since it is primarily through their atmospheric properties that were are able to characterise them. Modelling these young planets also involves facing several poorly contrained (and unconstrained) quantites which govern their atmospheric composition and dynamics.

Several theoretical studies have modelled the atmospheres and evolution of these young planets, but have all made several simplifying assumptions. @lichtenberg_vertically_2021 coupled a simple atmosphere model to real-gas radiative transfer and interior evolution codes to simulate magma ocean evolution, but did not account for the possibility of convective stability in their atmospheres. @krissansen_was_2021 made similar assumptions. @selsis_cool_2023 used a pure-steam radiative-convective model to investigate convective stability, finding that it is likely to occur, but did not extend their work to coupled scenarios which explore the secular evolution of these planets or gas mixtures. @zilinskas_observability_2023 used a multi-gas model to simulate the atmospheres of Sub-Neptune and Super-Earth exoplanets an equilibrium, finding that the opacity of various gaseous species (notably SiO) plays a key role in determining the structure and observable properties of exoplanet atmospheres. The demand for realistic atmosphere modelling in the context of secular magma ocean evolution is apparent.

Hydrostatic column atmosphere models, such as AGNI and all of those used in the aforementioned works, assume that three-dimensional planetary atmospheres can be reasonably represented by a single column. Making this assumption means that fewer calculations are required to spatially represent a planetary atmosphere, allowing computational resources to be expended on other physics. For AGNI in particular, radiative transfer is performed efficiently at medium and high resolutions with up to 4096 spectral bands distributed between 1 and 35000 cm$^{-1}$. This is done using the correlated-k and two-stream approximations [@amundsen_radiation_2014; @pierrehumbert_book_2010; @lacis_corrk_1991]. Ensuring sufficient spectral resolution is important in modelling the blanketing effect of these atmospheres, as capturing the opacity (and transparency) of the many molecular bands is known to be key in setting the rate at which these planets can cool by radiation to space. It is also important that we are able to run grids of models which explore the range of possible (and as-yet poorly constrained) conditions that these planets could exhibit, which demands efficient modelling given finite computational resources.

HELIOS is a hydrostatic atmosphere model written in Python and CUDA. It uses the same assumptions as AGNI

AGNI was developed with the view of being coupled into the modular PROTEUS framework for planetary evolution[^2].

[^2]: The PROTEUS framework can be found on GitHub [here](https://github.com/FormingWorlds/PROTEUS).



# Acknowledgements

Harrison Nicholls is supported by the Clarendon Fund and the MT Scholarship Trust.
The author is grateful for the continuing hard work of the Julia developers and those of its many libraries [@julialang].
Additionally, the author thanks James Manners, and others at the Met Office, who have continued to support and develop SOCRATES over the last three decades.

# References
