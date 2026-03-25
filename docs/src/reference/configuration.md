# Configuration reference

AGNI configuration files are formatted using [TOML](https://toml.io/en/). There are
examples in `res/config/`. There are no default values; the model will return an error
naming any parameters which are both necessary and absent.

See [Configuring AGNI](@ref) for guidance on choosing values for these parameters.

## `[planet]`
General properties of the planet.

| Parameter          | Description   |
| -----------------: | :------------ |
| `tmp_surf        ` | Temperature of the surface [kelvin].  |
| `instellation    ` | Stellar flux at planet's orbital distance [W m-2]. |
| `albedo_b        ` | A pseudo bond-albedo which downscales the stellar flux by 1-albedo_b. |
| `s0_fact         ` | Stellar flux scale factor which accounts for planetary rotation (c.f. Cronin+13). |
| `zenith_angle    ` | Characteristic zenith angle for incoming stellar radiation [degrees]. |
| `surface_material` | Surface material (can be "greybody" or path to data file in `res/`). |
| `albedo_s        ` | Spectrally-grey surface albedo used when `material=greybody`. Optional otherwise. |
| `radius          ` | Radius of the planet's surface [m]. |
| `gravity         ` | Gravitational acceleration at the surface [m s-2]. Incompatible with `mass` option. |
| `mass            ` | Total mass of material below the atmosphere [kg]. Incompatible with `gravity` option. |
| `skin_d          ` | Thickness of the conductive boundary layer [m]. Used when `solution_type=2`. |
| `skin_k          ` | Conductivity of the conductive boundary layer [W m-1 K-1]. Used when `solution_type=2`. |
| `tmp_magma       ` | Temperature of the topmost layer of the planet's mantle [K]. Used when `solution_type=2`. |
| `flux_int        ` | Internal flux [W m-2] to be solved-for when `solution_type=3`. |
| `roughness       ` | Surface roughness length scale [m]. |
| `wind_speed      ` | Effective wind speed for sensible heat transport [m s-1]. |
| `star_Teff       ` | Stellar photospheric temperature [K] used if `files.input_star=="blackbody"`. |


## `[files]`
Input/output files and other paths.

| Parameter          | Description   |
| -----------------: | :------------ |
| `input_sf       `  | Path to the desired spectral file in `res/spectral_files/`. If "greygas", uses double-grey RT scheme. |
| `input_star     `  | Path to stellar spectrum. If blank, spectrum assumed to be inside spectral file. If "blackbody" must provide `planet.star_Teff`. |
| `output_dir     `  | Path to the output directory. |
| `rfm_parfile  `    | Path to .par linelist file, for running line-by-line calculations with the RFM. |


## `[composition]`
Atmospheric composition and chemistry.
There are three main ways to set the composition: partial pressures (`p_dict`), mixing
ratios (`vmr_dict` or `vmr_file`), or by metallicities (`metallicities`).

| Parameter         | Description   |
| ----------------: | :------------ |
| `p_top         `  | Total top-of-atmosphere pressure [bar]. |
| `p_dict        `  | Dictionary of gas partial surface pressures [bar]. Summed to obtain `p_surf`. |
| `p_surf        `  | Total surface pressure [bar]. Incompatible with `p_dict`.|
| `vmr_dict      `  | Gas volume mixing ratios (=mole fractions) at the surface. Must also set `p_surf`. |
| `vmr_file      `  | Path to a file containing mixing ratio profiles. Must also set `p_surf`. |
| `metallicities`   | Dictionary of elemental **mass** abundance ratios relative to hydrogen. Must also set `p_surf`. |
| `condensates   `  | List of volatiles which are allowed to condense. Can be used together with thermochemical equilibrium (`physics.chemistry = true`). |
| `transparent   `  | Make the atmosphere transparent (see below). Replaces all of the above parameters in this table. |


## `[execution]`
Parameters that tell the model what to do.

| Parameter         | Description   |
| ----------------: | :------------ |
| `clean_output  `  | Clean old files from the output folder at model startup (true/false). |
| `verbosity     `  | Logging verbosity (0: quiet, 1: normal, 2: extra logging). |
| `max_steps     `  | Maximum number of steps the solver should take before giving up (typically <200). |
| `max_runtime   `  | Maximum wall-clock runtime [s] before giving up. |
| `num_levels    `  | Number of model levels. Typically ~50, and ideally less than 100.  |
| `converge_atol `  | Convergence criterion absolute tolerance [W m-2]. |
| `converge_rtol `  | Convergence criterion relative tolerance [dimensionless]. |
| `converge_type `  | Definition of convergence criterion (1: cost function, 2: median resid, 3: mean resid). |
| `initial_state `  | Ordered list of requests describing the initial state of the atmosphere (see [Configuring AGNI](@ref)). |
| `solution_type `  | Solution type (see [Configuring AGNI](@ref) and [Model description](@ref)). |
| `solver        `  | Solver to use (see [Configuring AGNI](@ref)). |
| `dx_max        `  | Maximum step size [kelvin] allowed to be taken by the solver during each step. |
| `linesearch    `  | Linesearch method to be used (0: None, 1: Golden section, 2: Backtracking). |
| `easy_start    `  | Initially scale energy fluxes, to help with stability if the model is struggling. |
| `grey_start    `  | Initially solve with double-grey RT scheme, to help with stability if the model is struggling. |
| `perturb_all`     | Perturb all rows of Jacobian matrix at each solver iteration? True=stable, False=fast. |
| `rfm_wn_min`      | Line-by-line RFM radiative transfer, minimum wavenumber [cm-1]. |
| `rfm_wn_max`      | Line-by-line RFM radiative transfer, maximum wavenumber [cm-1]. |

## `[physics]`
Parameters that describe how the model should treat the physics.

| Parameter         | Description   |
| ----------------: | :------------ |
| `chemistry     `  | Include 1D equilibrium chemistry in the atmosphere (true/false). |
| `continua      `  | Include collisional/continuum absorption in radiative transfer (true/false). |
| `rayleigh      `  | Include Rayleigh scattering in radiative transfer (true/false). |
| `aerosol       `  | Include aerosols in radiative transfer  (true/false). |
| `cloud         `  | Include water clouds in radiative transfer (true/false). |
| `overlap_method`  | Method for treating overlapping gas opacities within a given spectral band (see [Configuring AGNI](@ref)). |
| `grey_lw`         | Grey opacity [m2/kg] of longwave (thermal) radiation. Used when `input_sf="greygas"`. |
| `grey_sw`         | Grey opacity [m2/kg] of shortwave (stellar) radiation. Used when `input_sf="greygas"`. |
| `real_gas      `  | Use real-gas equation(s) of state where possible (true/false). |
| `thermo_funct  `  | Use temperature-dependent thermodynamic properties (true/false). |
| `sensible_heat `  | Include turbulent sensible heat transport at the surface (true/false). |
| `convection    `  | Include vertical heat transport associated with convection (true/false). |
| `convection_crit` | Criterion for convective stability. Options: (s)chwarzschild, (l)edoux. |
| `latent_heat   `  | Include vertical heat transport from condensation and evaporation (true/false). |
| `rainout       `  | Enable condensation and evaporation of condensables aloft. Required for `latent_heat=true`. |
| `oceans        `  | Enable condensation and evaporation of condensables at the surface. |

## `[aerosols]`
Parameters controlling SOCRATES 'classic' aerosol parametrizations.

| Parameter           | Description   |
| ------------------: | :------------ |
| `rel_humidity    `  | Mean relative humidity used for moist aerosol optical lookups [0,1]. |
| `avg_phase_moments` | Number of scattering phase-function moments to include. |
| `species_mmr     `  | Dictionary of aerosol mass-mixing-ratio overrides [kg/kg], keyed by SOCRATES aerosol suffix (e.g. `soot`, `dustdiv1`, `biogenic`). |


## `[plots]`
Configure plotting routines; all of these should be `true` or `false`.

| Parameter         | Description   |
| ----------------: | :------------ |
| `at_runtime     ` | Make some plots at runtime? |
| `temperature    ` | Plot temperature-pressure profile? |
| `fluxes         ` | Plot energy flux profiles? |
| `contribution   ` | Plot spectral contribution function? |
| `emission       ` | Plot outgoing emission spectrum? |
| `albedo         ` | Plot spectral albedo? This is the ratio of upward:downward SW fluxes. |
| `mixing_ratios  ` | Plot mixing ratio profiles? |
| `height         ` | Plot radius-pressure profile? |
| `animate        ` | Make an animation of the solver obtaining its solution? |
| `cloud          ` | Plot water cloud mass fraction and area fraction profiles? |
