Surface radiative properties lookup data. These values are passed to the variable `rho_alb` in SOCRATES.

AGNI is flexible in taking any form in two columns
 * first column is WL in nm
 * second column describes the optical properties

Options for second column:
 * r = spherical reflectance (AKA Bond albedo), converted to e in code
 * e = spherical emissivity, converted to r in code
 * w = single scattering albedo, converted to r and e in code

----------------

reflectivity: r_hh = r_s = r0 * (1 - γ/(3+3γ) )
emissivity: e_h = 1- r_s

r0 is the diffusive reflectance = (1 − γ)/(1 + γ)
γ is the albedo factor = sqrt(1 − w)

Ecospec and biospec use r_dh = r_h.
Can be converted to gamma if angle is known.
γ = (1 − r_h) / (1 + 2μ0 * rh)

Files can be obtained from Zenodo: https://zenodo.org/communities/proteus_framework/records?q&f=subject%3Asurface_albedos&l=list&p=1&s=10&sort=newest

