var documenterSearchIndex = {"docs":
[{"location":"model/#Model-description","page":"Model description","title":"Model description","text":"","category":"section"},{"location":"model/","page":"Model description","title":"Model description","text":"AGNI models a planetary atmosphere by treating it as a single column (1D) and splitting it up into levels of finite thickness. These levels are defined in pressure-space, and are arranged logarithmically between the surface and the top of the atmosphere. The atmosphere is assumed to be plane-parallel. Quantities such as pressure and temperature are calculated at level-centres and level-edges, while energy fluxes are calculated only at the edges, and thermodynamic properties (e.g. heat capacity) are calculated only at their centres.","category":"page"},{"location":"model/#Radiative-transfer","page":"Model description","title":"Radiative transfer","text":"","category":"section"},{"location":"model/","page":"Model description","title":"Model description","text":"Radiative transfer (RT) refers to the transport of radiation energy through a medium subject to the characteristics of the medium. Radiation passing through an atmosphere is absorbed, emitted, scattered, and reflected. In the context of planetary atmospheres, we also have to handle their surfaces, cloud formation, and radiation from the host star.","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"AGNI simulates RT using SOCRATES, a numerical code written by the UK Met Office which solves the RT equation using a two-stream solution. SOCRATES is accessed using a Julia interface originally written by Stuart Daines. Opacity is handled using the correlated-k approximation, with either random overlap or equivalent extinction used to account for overlapping absorption in mixtures of gases.","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"The model uses k-terms fitted to spectral absorption cross-section data from DACE. The MT_CKD model is used to estimate water continuum absorption cross-sections. Other continuua are derived from the HITRAN tables. Rayleigh scattering and water cloud radiative properties are also included. You can find tools for fitting k-terms and processing line absorption data in my redistribution of SOCRATES on GitHub. The flowchart below outlines how these absorption data are converted into a 'spectral file'.","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"  <img src=\"assets/spectral_flowchart.svg\" width=100% class=\"center\"/>","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"Surface reflectivity can be modelled as a greybody with an albedo from 0 to 1. Alternatively, it can be modelled from empirical single-scattering data which varies with zenith angle and wavelength.","category":"page"},{"location":"model/#Convection","page":"Model description","title":"Convection","text":"","category":"section"},{"location":"model/","page":"Model description","title":"Model description","text":"Convection is a turbulent process which occurs across more than one spatial dimension, so it must be parameterised within 1D models like AGNI. In fact, it is typically parameterised inside 3D global circulation models as resolving convection is numerically expensive. AGNI uses mixing length theory (MLT) to parameterise convection. This is in contrast to convective adjustment, which forcibly adjusts a convectively unstable region of the atmosphere to the corresponding adiabat while ensuring that enthalpy is conserved.","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"MLT directly calculates the energy flux associated with convective heat transport, and thus is the preferred parameterisation within the model. It assumes that parcels of gas are diffused over a characteristic mixing length, transporting energy in the process. This requires choosing a scale for this mixing length, but in practice this has very little impact on the results from the model.","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"The atmosphere is assumed to be hydrostatically supported. Gas densities are combined using Amagat's additive volume law. The densities of each gas are nominally calculated using the Van der Walls equation of state (EOS). AQUA is implemented as the EOS for water. The Chabrier+2019 EOS is implemented as the EOS for hydrogen. AGNI will fallback to the ideal gas EOS for otherwise unsupported gases.","category":"page"},{"location":"model/#Phase-change","page":"Model description","title":"Phase change","text":"","category":"section"},{"location":"model/","page":"Model description","title":"Model description","text":"Gases release energy (\"latent heat\" or \"enthalpy\") into their surroundings when condensing into a liquid or solid. This is included in the model through a diffusive condensation scheme, which assumes a fixed condensation timescale. This takes place as follows... firstly, the mixing ratios of the gases are updated according to the temperature profile, where rainout occurs until all condensibles are saturated or sub-saturated. The mixing ratios of dry species are increased in order to satisfy the total pressure at condensing levels. The heat released associated with the change in partial pressure of condensible gases is used to calculate a latent heating rate. This is then integrated (from the TOA downwards) to provide a latent heat transport flux at cell-edges. The integrate condensible heat flux is balanced by evaporation at deeper layers.","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"Latent heats are temperature-dependent, using values derived from Coker (2007) and Wagner & Pruß (2001). Heat capacities are also temperature-dependent, using values derived from the JANAF database. See the ThermoTools repo for scripts.","category":"page"},{"location":"model/#Solar-flux","page":"Model description","title":"Solar flux","text":"","category":"section"},{"location":"model/","page":"Model description","title":"Model description","text":"A key input to the radiation model is the shortwave downward-directed flux from the star at the top of the atmosphere. This is quantified by the bolometric instellation flux, a scale factor, an artificial additional albedo factor, and a zenith angle. All of these may be provided to the model through the configuration file. The model also requires a stellar spectrum scaled to the top of the atmosphere.","category":"page"},{"location":"model/#Obtaining-a-solution","page":"Model description","title":"Obtaining a solution","text":"","category":"section"},{"location":"model/#Summary","page":"Model description","title":"Summary","text":"","category":"section"},{"location":"model/","page":"Model description","title":"Model description","text":"AGNI is designed for modelling planetary atmospheres with high surface pressures and temperatures. This means that the radiative timescale differs by several orders of magnitude across the column, which makes obtaining a solution difficult. To obtain a temperature structure solution that conserves energy more precisely than a time-stepping method, AGNI solves for the temperature structure of the atmosphere as an optimisation problem. This finds the state which conserves energy across all levels and satisfies the required configuration from the user.","category":"page"},{"location":"model/#Solution-types","page":"Model description","title":"Solution types","text":"","category":"section"},{"location":"model/","page":"Model description","title":"Model description","text":"It is necessary to tell AGNI what kind of atmospheric solution to solve for. There are currently a few options available set by the solution_type variable in the configuration file.","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"(1) Aim to conserve energy fluxes throughout the column. The surface temperature is fixed.\n(2) Aim to conserve energy fluxes throughout the column. The surface temperature is set by energy transport through a solid conductive boundary layer of thickness d such that T_s = T_m - fracFdk, where T_m is the mantle temperature and k is the thermal conductivity.\n(3) Solve for a state such that the flux carried at each level is equal to F_textnet = sigma T_textnet^4, representing the rate at which a planet is losing energy into space.","category":"page"},{"location":"model/#Construction","page":"Model description","title":"Construction","text":"","category":"section"},{"location":"model/","page":"Model description","title":"Model description","text":"The atmosphere is constructed of N levels (cell-centres), corresponding to N+1 interfaces (cell-edges). The RT model takes cell-centre temperatures T_i, pressures p_i, geometric heights, and mixing ratios as input variables at each level i. As well as the surface temperature and incoming stellar flux. In return, it provides cell-edges spectral fluxes F_i at all N+1 interfaces for LW & SW components and upward & downward streams. Convective fluxes can be estimated using the MLT scheme, condensation fluxes from the condensation scheme, and sensible heat from a simple TKE approximation.","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"The total upward-directed energy flux F_i describes the total upward-directed energy transport (units of textW m^-2) from cell i into cell i-1 above (or into space for i=1). For energy to be conserved throughout the column, it must be true that F_i = F_t text  forall text  i where F_t is the total amount of energy being transported out of the planet. In global radiative equilibrium, F_t = 0.","category":"page"},{"location":"model/#Definition-of-residuals","page":"Model description","title":"Definition of residuals","text":"","category":"section"},{"location":"model/","page":"Model description","title":"Model description","text":"We can use this construction to solve for the temperature profile of the atmosphere as an N+1-dimensional optimisation problem. This directly solves for T(p) at radiative-convective equilibrium without having to invoke heating rate calculations, thereby avoiding slow convergence in regions of the atmosphere with long radiative timescales. The residuals vector (length N+1)","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"\nbmr =\n\nbeginpmatrix\nr_i     \nr_i+1 \n     \nr_N     \nr_N+1\nendpmatrix\n\n=\n\nbeginpmatrix\nF_i+1 - F_i     \nF_i+2 - F_i+1 \n     \nF_N+1 - F_N     \nF_N+1 - F_t\nendpmatrix","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"is what we aim to minimise as our 'objective function', subject to the solution vector of cell-centre temperatures","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"bmx =\n\nbeginpmatrix\nT_i     \nT_i+1 \n     \nT_N     \nT_s\nendpmatrix","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"where T_s is the surface temperature. Cell-edge temperatures in the bulk atmosphere are interpolated from cell-centres. The bottom- and top-most cell edge temperatures are extrapolated by estimation of dTd log p. Cell properties (heat capacity, gravity, density, average molecular weight, etc.) are consistently updated at each evaluation of bmr. Condensation/rainout are also handled at each evaluation of bmr in order to avoid supersaturation.","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"The model converges when the cost function c(bmx) = sqrtsum_i r_i satisfies the condition","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"c(bmx)  c_a + c_r cdot undersetimax text  F_i","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"which represents a state where the fluxes are sufficiently conserved.","category":"page"},{"location":"model/#Iterative-steps","page":"Model description","title":"Iterative steps","text":"","category":"section"},{"location":"model/","page":"Model description","title":"Model description","text":"The model solves for bmx iteratively, starting from some initial guess. The initial guess should be any reasonable temperature profile which is not significantly cooler than the expected solution. The flowchart below broadly outlines the solution process.","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"  <img src=\"assets/model_flowchart.svg\" width=50% class=\"center\"/>","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"The Jacobian matrix bmJ represents the directional gradient of the residuals with respect to the solution vector. It is a square matrix with elements set according to","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"J_uv = fracpartial r_upartial x_v","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"AGNI estimates bmJ using finite-differences, requiring N+1 evalulations of bmr in order to fill the matrix. This corresponds to 2(N+1)+1 objective function calculations under a 2nd order central-difference scheme. Each level v with temperature x_v is perturbed by an amount pm varepsilon x_v in order to fill a single column of bmJ. As such, it can be expensive to construct a full Jacobian, especially when it is discarded at the end of each iteration. To reduce the total number of calculations AGNI retains some of the columns in bmJ between model iterations. This assumes that the second derivative of the residuals is small. A column v is retained only when","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"max r_i lt 07 text for  i in v-1 v v+1","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"and when c(bmx)10 does not satisfy the convergence criteria.","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"With a Jacobian constructed, we can calculate an update bmd to the solution vector bmx rightarrow bmx + bmd. This is primarily done via the Newton-Raphson method","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"bmd = -bmJ^-1 bmr","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"but can alternatively be performed via the Gauss-Newton and Levenberg–Marquardt methods.","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"It is possible for the model to become stuck in a local minimum, leading to very small values of bmd. This is identified when c(bmx) has seen little change over the last few iterations. When this occurs, the model is 'nudged' by scaling the update via","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"bmd rightarrow 3 bmd","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"In many cases bmd is too large, leading to instabilities. This is due to the non-convexity of the solution space, and the somewhat discontinuous nature of the physics involved (particularly in its temperature derivatives). When bmdd_textmax, the update is crudely scaled via","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"bmd rightarrow  d_textmax hatbmd","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"The update may also be scaled by a linesearch","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"bmd rightarrow alpha bmd","category":"page"},{"location":"model/","page":"Model description","title":"Model description","text":"on alpha. This is applied if the full step bmd would increase the cost by an unacceptable amount. If the model is close to convergence then a golden-section search method is used to determine the optimal alpha, otherwise a backtracking method is used.  This means that the model is (mostly) able to avoid oscillating around a solution. All three of these scalings to bmd preserve its direction.","category":"page"},{"location":"model/#Other-features","page":"Model description","title":"Other features","text":"","category":"section"},{"location":"model/","page":"Model description","title":"Model description","text":"AGNI can calculate emission spectra, provided with T(p) and the volume mixing ratios of the gases. This is performed using the same RT as the RCE calculations, so is limited in resolution by the choice of correlated-k bands. Similarly, the longwave contribution function can also be calculated.","category":"page"},{"location":"model/#Julia-and-Fortran","page":"Model description","title":"Julia and Fortran","text":"","category":"section"},{"location":"model/","page":"Model description","title":"Model description","text":"AGNI is primarily written in Julia, while SOCRATES itself is written in Fortran. Julia was chosen because it allows the SOCRATES binaries to be included in the precompiled code, which significantly improves performance.","category":"page"},{"location":"usage/#Running-the-model","page":"Running the model","title":"Running the model","text":"","category":"section"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"First, follow the Getting started instructions. Only read-on once you have confirmed that the code is working.","category":"page"},{"location":"usage/#Input-data-files","page":"Running the model","title":"Input data files","text":"","category":"section"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"The minimal input data required to run the model will have been downloaded automatically. If you require more data, such as additional stellar spectra or opacities, then these can also be easily obtained using the get_data script in the AGNI root directory. To see how to use this script, run it without arguments like so:","category":"page"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"./src/get_data.sh","category":"page"},{"location":"usage/#Tutorials","page":"Running the model","title":"Tutorials","text":"","category":"section"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"There are Jupyter notebooks containing tutorials in the tutorials/ directory of the repository.","category":"page"},{"location":"usage/#General-execution","page":"Running the model","title":"General execution","text":"","category":"section"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"The environment variable RAD_DIR must point to the SOCRATES installation directory. This is required for AGNI to find the SOCRATES libraries. The best way to do this is to add RAD_DIR=path/to/socrates/folder/ to your shell rc file (e.g. ~/.bashrc).","category":"page"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"Then to use the model, simply run ./agni.jl [cfg] where [cfg] is the path to the required configuration file. If [cfg] is not passed, then the default configuration will be used.","category":"page"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"To calculate equilibrium chemistry self-consistently with the climate, FastChem is coupled to AGNI. For AGNI to find FastChem, the environment variable FC_DIR must point to the FastChem installation directory.","category":"page"},{"location":"usage/#Configuration","page":"Running the model","title":"Configuration","text":"","category":"section"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"AGNI configuration files are formatted using TOML. There are examples in res/config/. The default configuration file contains comments explaining the purpose of each parameter, although some are explained in greater detail below. Take care to format the variables in the TOML file correctly. There are no 'default values'. Not all parameters are required in all cases, but the model will return an error naming any parameters which are both necessary and absent.","category":"page"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"Broadly, the configuration files are broken up into four sections:","category":"page"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"[planet] -  general characteristics of the planet\n[files] - input/output files and other paths\n[composition] - atmospheric composition and chemistry\n[execution] - what the model should do\n[plots] - which plots should be produced","category":"page"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"Some parameters:","category":"page"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"files.input_sf is the file path to the \"spectral file\" containing opacity data. Several spectral files are packged with AGNI, but you can find more online via the Open Science Framework.\nexecution.solution_type tells the model which state to solve for. The allowed values (integers) are...\n1 : zero flux divergence at fixed tmp_surf\n2 : zero flux divergence, with tmp_surf set such that the conductive skin (CBL) conserves energy flux\n3 : the net flux (up minus down) at each layer is equal to flux_int\nexecution.solver tells the model which solver to use. Allowed solvers are...\n[empty string] : no solving takes place, so the model just calculates fluxes using the initial state\nnewton : the Newton-Raphson algorithm is used\ngauss  : the Gauss-Newton algorithm is used\nlevenberg : the Levenberg–Marquardt algorithm is used\nexecution.initial_state describes the initial temperature profile applied to the atmosphere. This is a list of strings which are applied in the given order, which allows the user to describe a specific state as required. The descriptors are listed below, some of which take a single argument that needs to immediately follow the descriptor in the list order.\ndry              : integrate the dry adiabatic lapse rate from the surface upwards\nstr,       arg : apply an isothermal stratosphere at arg kelvin\niso,       arg : set the whole atmosphere to be isothermal at arg kelvin\ncsv,       arg : set the temperature profile using the CSV file at the file path arg\nsat,       arg : apply Clausius-Clapeyron saturation curve for the gas arg\nncdf,      arg : load profile from the NetCDF file located at arg\nloglin,    arg : log-linear profile between tmp_surf at the bottom and arg at the top\nFor example, setting initial_state = [\"dry\", \"sat\", \"H2O\", \"str\", \"180\"] will set T(p) to follow the dry adiabat from the surface, the water condensation curve above that, and then to be isothermal at 180 K until the top of the model.\ncomposition.chem_type describes the type of chemistry to implement within the model. This is handled externally by FastChem, so you must set the environment variable FC_DIR to point to the FastChem directory. The allowed values (integers) are...\n0 : Disabled\n1 : Equilibrium, gas phase only\n2 : Equilibrium, with condensation (condensates retained)\n3 : Equilibrium, with condensation (condensates rained out)","category":"page"},{"location":"usage/#Outputs","page":"Running the model","title":"Outputs","text":"","category":"section"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"Results are optionally plotted and animated, and data will be saved as NetCDF or CSV files.","category":"page"},{"location":"usage/#Python","page":"Running the model","title":"Python","text":"","category":"section"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"It is possible to interact with the model using Python. This is best done with the juliacall package from PythonCall.jl, and is implemented this way in the PROTEUS framework.","category":"page"},{"location":"troubleshooting/#Troubleshooting","page":"Troubleshooting","title":"Troubleshooting","text":"","category":"section"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"This page may be useful if you are having problems. However, I would suggest that you also double check the Getting started instructions.","category":"page"},{"location":"troubleshooting/#Julia-errors-on-start,-potentially-referencing-the-CURL-library","page":"Troubleshooting","title":"Julia errors on start, potentially referencing the CURL library","text":"","category":"section"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"It is important that the shell environment variable LD_LIBRARY_PATH is not set when running AGNI. This will cause Julia to use the wrong libraries, which will causes problems. You can unset this variable or reset using either of the following commands","category":"page"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"unset LD_LIBRARY_PATH\nexport LD_LIBRARY_PATH=\"\"","category":"page"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"If this does not help, it's possible that you are using a Julia distribution provided by your system package manager. It's important that you only use Julia distributed from the official website.","category":"page"},{"location":"troubleshooting/#Cannot-find-SOCRATES","page":"Troubleshooting","title":"Cannot find SOCRATES","text":"","category":"section"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"Check the installation instructions. Have you set RAD_DIR? Try running l_run_cdf in the terminal; if this fails, then SOCRATES has not compiled or you haven't added it to your PATH. It is necessary to set the RAD_DIR variable for the environment in which you are running AGNI, so it is best to add it to your shell's rc file permanently.","category":"page"},{"location":"troubleshooting/#Spectral-file-does-not-exist","page":"Troubleshooting","title":"Spectral file does not exist","text":"","category":"section"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"Check the path in the configuration file. Download additional spectral files using the get_data script. For example, for additional pure-steam spectral files you would run","category":"page"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"./src/get_data.sh steam","category":"page"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"Other spectral files can be downloaded from OSF: https://osf.io/vehxg/.","category":"page"},{"location":"troubleshooting/#Cannot-find-FastChem","page":"Troubleshooting","title":"Cannot find FastChem","text":"","category":"section"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"You need to install FastChem. This can be done by running the command:","category":"page"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"./src/get_fastchem.sh","category":"page"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"and then adding FC_DIR to your shell rc file.","category":"page"},{"location":"troubleshooting/#Finally...","page":"Troubleshooting","title":"Finally...","text":"","category":"section"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"If you are still stuck, or feel that there is a problem with the code, then you can contact the authors using the information on the main page.","category":"page"},{"location":"manual/#Development-manual","page":"Development","title":"Development manual","text":"","category":"section"},{"location":"manual/#Contributing","page":"Development","title":"Contributing","text":"","category":"section"},{"location":"manual/","page":"Development","title":"Development","text":"If you are interested in contributing to the model, please contact the developers using the information on the main page.","category":"page"},{"location":"manual/#Coding-style","page":"Development","title":"Coding style","text":"","category":"section"},{"location":"manual/","page":"Development","title":"Development","text":"Indentation uses 4 spaces, no tabs.\nFunction names should be lowercase, with words separated by underscores .\nLines should aim to have a length of no more than 92 characters.\nAll functions should have docstrings, ideally with Arguments and Returns listed.\nMore comments are always better, even if they seem redundant.\nUse type hinting where possible.\nPrint statements should be made through the logger where possible.\nThe core package code should not contain global variables, except in the phys module.","category":"page"},{"location":"manual/#Code-reference","page":"Development","title":"Code reference","text":"","category":"section"},{"location":"manual/","page":"Development","title":"Development","text":"To be completed.","category":"page"},{"location":"examples/#Example-outputs","page":"Example outputs","title":"Example outputs","text":"","category":"section"},{"location":"examples/#Pure-steam-runaway-greenhouse-effect","page":"Example outputs","title":"Pure steam runaway greenhouse effect","text":"","category":"section"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"By assuming the atmosphere temperature profile follows a dry adiabat and the water vapour-condensate coexistance curve defined by the Clausius-Claperyron relation, we see a characteristic relationship between the outgoing longwave radiation (OLR) and the surface temperature (T_s). Initially OLR increases with T_s, but as the condensing layer (which is independent of T_s) overlaps with the photosphere, OLR and T_s decouple. Eventually the atmosphere reaches a dry post-runaway state, and OLR increases rapidly with T_s.","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"  <img src=\"assets/runaway/curve.png\" width=75% class=\"center\"/>","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"You can find a Jupyter notebook which reproduces this result in the tutorials directory of the repository.","category":"page"},{"location":"examples/#Prescribed-convective-case","page":"Example outputs","title":"Prescribed convective case","text":"","category":"section"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"In this case, a temperature profile is prescribed to follow a dry adiabat from the surface to a moist region, and then a pseudoadiabat to the top of the atmosphere. This is in line with previous works and the OLR curve above.","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"  <img src=\"assets/nosolve/plot_ptprofile.png\" width=62% class=\"center\"/>","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"Radiative fluxes are then calculated according to this  temperature profile. Because the profile is prescribed, the fluxes are not balanced locally or globally across the column.","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"  <img src=\"assets/nosolve/plot_fluxes.png\" width=62% class=\"center\"/>","category":"page"},{"location":"examples/#Radiative-convective-solution","page":"Example outputs","title":"Radiative-convective solution","text":"","category":"section"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"Instead, we can model an atmosphere such that that energy is globally and locally conserved. Convection is parameterised using mixing length theory in this case, allowing the system to be solved using a Newton-Raphson method. In the convective region at ~0.1 bar, we can see that the radiative fluxes and convective fluxes entirely cancel, because AGNI was asked to solve for a case with zero total flux transport.","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"  <img src=\"assets/withsolve/plot_ptprofile.png\" width=62% class=\"center\"/>\n  <br />\n  <img src=\"assets/withsolve/plot_fluxes.png\" width=62% class=\"center\"/>","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"We can also plot the outgoing emission spectrum and normalised longwave contribution function (CF). The spectrum clearly demonstrates complex water absorption features, and exceeds blackbody emission at shorter wavelengths due to Rayleigh scattering. The CF quantifies how much each pressure level contributes to the outgoing emission spectrum at a given wavelength – this is then plotted versus wavelength and pressure.","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"  <img src=\"assets/withsolve/plot_emission.png\" width=90% class=\"center\"/>","category":"page"},{"location":"","page":"Home","title":"Home","text":"    <img src=\"assets/logo_title.svg\" width=32% class=\"center\"/>\n    <p align=\"center\">\n        <b>A radiative-convective model for lava planet atmospheres</b>\n    </p>","category":"page"},{"location":"","page":"Home","title":"Home","text":"A numerical model for the atmospheres of hot rocky (exo)planets. AGNI's primary purpose is to simulate the evolving atmospheres of magma ocean planets, while ensuring that radiative-convective equilibrium is sufficiently maintained. Pronounced as ag-nee. Named after the fire deity of Hinduism.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Follow Getting started for information on installing the code and obtaining results.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Contact: harrison[dot]nicholls[at]physics.ox.ac.uk","category":"page"},{"location":"","page":"Home","title":"Home","text":"GitHub: https://github.com/nichollsh/AGNI","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you use AGNI, please cite the following papers:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Nicholls et al., (2024) - DOI 10.1093/mnras/stae2772\nNicholls et al., (2025) - submitted to JOSS","category":"page"},{"location":"","page":"Home","title":"Home","text":"This software is available under the GPLv3. Copyright (C) 2025 Harrison Nicholls.","category":"page"},{"location":"setup/#Getting-started","page":"Getting started","title":"Getting started","text":"","category":"section"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"This page outlines requirements and installation steps for the code. Currently, GNU/Linux and MacOS (including ARM) are supported.","category":"page"},{"location":"setup/#Requirements","page":"Getting started","title":"Requirements","text":"","category":"section"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"gfortran\nNetCDF library for FORTRAN\nmake\ncurl","category":"page"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"warning: Warning\nDo not install Julia using your system package manager. Install only from julialang.org","category":"page"},{"location":"setup/#Installation","page":"Getting started","title":"Installation","text":"","category":"section"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"Follow the steps below in order to setup the code.","category":"page"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"Install Julia: curl -fsSL https://install.julialang.org | sh\nDownload AGNI: git clone https://github.com/nichollsh/AGNI.git\nChange directory: cd AGNI\nSetup SOCRATES by doing either ONE of the following...\nFollow the instructions on the SOCRATES GitHub page\nRun source src/get_socrates.sh\nbash src/get_agni.sh","category":"page"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"AGNI is now installed as a package into a Julia environment in the AGNI directory. This will also have downloaded some basic input data, and have run the tests.","category":"page"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"tip: Tip\nThe get_socrates.sh script automatically adds the radiation code to your environment with the variable RAD_DIR, which points to the SOCRATES installation. This variable must be set whenever AGNI is being used.","category":"page"},{"location":"setup/#Testing","page":"Getting started","title":"Testing","text":"","category":"section"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"If you want to run the tests manually, simply use the script in the test/ folder...","category":"page"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"julia test/runtests.jl","category":"page"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"This will print information on whether tests passed or failed.","category":"page"},{"location":"setup/#Updating","page":"Getting started","title":"Updating","text":"","category":"section"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"It's important that you keep AGNI up to date, especially if you are using as part of the PROTEUS framework. Use this script to automatically pull changes from GitHub and download any required data files.","category":"page"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"bash src/get_agni.sh","category":"page"},{"location":"setup/#Using-the-code","page":"Getting started","title":"Using the code","text":"","category":"section"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"See Running the model for information on using the code. See Troubleshooting for troubleshooting advice.","category":"page"},{"location":"setup/#Coupling-with-FastChem","page":"Getting started","title":"Coupling with FastChem","text":"","category":"section"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"This can be enabled using the configuration file parameter composition.chem_type. Of course, it is first necessary to setup FastChem, which can be done by running source src/get_fastchem.sh and then setting the FC_DIR environment variable.","category":"page"}]
}
