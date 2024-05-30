var documenterSearchIndex = {"docs":
[{"location":"usage/#Running-the-model","page":"Running the model","title":"Running the model","text":"","category":"section"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"First, follow the Getting started instructions. Only read-on once you have confirmed that the code is working.  ","category":"page"},{"location":"usage/#Execution","page":"Running the model","title":"Execution","text":"","category":"section"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"To run the model, simply execute ./agni.jl [cfg] where [cfg] is the path to the required configuration file. If [cfg] is not passed, then the default configuration will be used.","category":"page"},{"location":"usage/#Configuration","page":"Running the model","title":"Configuration","text":"","category":"section"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"AGNI configuration files are formatted using TOML. There are examples in res/config/.  The default configuration file contains comments explaining the purpose of each parameter, although some are explained in greater detail below. Take care to format the variables in the TOML file correctly. There are no 'default values'. Not all parameters are required in all cases, but the model will return an error naming any parameters which are both necessary and absent.","category":"page"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"Broadly, the configuration files are broken up into four sections:","category":"page"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"[planet] describes the physical characteristics of the planet\n[files] lists input/output files\n[execution] describes what the model should do\n[plots] describes which kind of plot to produce","category":"page"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"Specific parameters:","category":"page"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"solution_type tells the model which state to solve for. The allowed values (integers) are...\n0 : zero flux divergence at fixed tmp_surf, extrapolated tmpl[end]\n1 : zero flux divergence at fixed tmp_surf, constant tmpl[end]\n2 : zero flux divergence, with tmp_surf set such that the conductive skin (CBL) conserves energy flux\n3 : the net upward flux at each layer is equal to flux_eff = sigma * tmp_eff^4\n4 : zero flux divergence and OLR = target_olr\nsolvers tells the model which solvers to use. This is a list of strings, so multiple solvers can be applied sequentially. An empty string is always appended to the end of this list. Allowed solvers are...\n[empty string] : no solving takes place, so the model just calculates fluxes using the initial state\nnewton : the Newton-Raphson algorithm is used\ngauss  : the Gauss-Newton algorithm is used \nlevenberg : the Levenberg–Marquardt algorithm is used \ntimestep : a timestepping method is applied using the Adams-Bashforth method\ninitial_state describes the initial temperature profile applied to the atmosphere. This is a list of strings which are applied in the given order, which allows the user to describe a specific state as required. The descriptors are listed below, some of which take a single argument that needs to immediately follow the descriptor in the list order.\ndry : integrate the dry adiabatic lapse rate from the surface upwards\nstr, arg : apply an isothermal stratosphere at arg kelvin\niso, arg : set the whole atmosphere to be isothermal at arg kelvin\ncsv, arg : set the temperature profile using the CSV file at the file path arg\ncon, arg : apply Clausius-Clapeyron saturation curve for the gas arg\nsat : ensure that no supersaturation occurs at the surface by removing gases as required    \nFor example, setting initial_state = [\"dry\", \"sat\", \"H2O\", \"str\", \"180\"] will set T(p) to follow the dry adiabat from the surface, the water condensation curve above that, and then to be isothermal at 180 K until the top of the model.\nchem_type describes the type of chemistry to implement within the model. This is handled externally by FastChem. You must also provide the path to the FastChem installation directory fastchem_path in the [files] section. The allowed values (integers) are...\n0 : Disabled \n1 : Equilibrium (gas only)\n2 : Equilibrium (condensates retained)\n3 : Equilibrium (condensates rained out)","category":"page"},{"location":"usage/#Outputs","page":"Running the model","title":"Outputs","text":"","category":"section"},{"location":"usage/","page":"Running the model","title":"Running the model","text":"Results are optionally plotted and animated, and data will be saved as NetCDF or CSV files. ","category":"page"},{"location":"model/#Description","page":"Description","title":"Description","text":"","category":"section"},{"location":"model/","page":"Description","title":"Description","text":"AGNI models a planetary atmosphere by treating it as a single column (1D) and splitting it up into levels of finite thickness. These levels are defined in pressure-space, and are arranged logarithmically between the surface and the top of the atmosphere. Quantities such as pressure and temperature are calculated at level-centres and level-edges, while energy fluxes are calculated only at the edges, and thermodynamic properties (e.g. heat capacity) are calculated only at their centres.","category":"page"},{"location":"model/#Radiative-transfer","page":"Description","title":"Radiative transfer","text":"","category":"section"},{"location":"model/","page":"Description","title":"Description","text":"Radiative transfer (RT) refers to the transport of radiation energy through a medium subject to the characteristics of the medium. Radiation passing through an atmosphere is absorbed, emitted, scattered, and reflected. In the context of planetary atmospheres, we also have to handle their surfaces, cloud formation, and radiation from the host star.","category":"page"},{"location":"model/","page":"Description","title":"Description","text":"AGNI models RT using SOCRATES, a numerical code written by the UK Met Office which solves the RT equation using a two-stream solution under a plane-parallel approximation. SOCRATES is accessed using a Julia interface originally written by Stuart Daines. The atmosphere is assumed to be hydrostatically supported and to behave as an ideal gas. Opacity is calculated using the correlated-k approximation, with either random overlap or equivalent extinction used to account for overlapping absorption in mixtures of gases. ","category":"page"},{"location":"model/","page":"Description","title":"Description","text":"For simulating gaseous absorption, the model fits k-terms to spectral absorption cross-section data from DACE. The MT_CKD model is used to estimate continuum absorption cross-sections. Rayleigh scattering and water cloud radiative effects are also included. You can find tools for fitting k-terms and processing line absorption data in my redistribution of SOCRATES on GitHub.","category":"page"},{"location":"model/#Convection","page":"Description","title":"Convection","text":"","category":"section"},{"location":"model/","page":"Description","title":"Description","text":"Convection is a process that occurs across more than one spatial dimension, so it must be parameterised within 1D models like AGNI. In fact, it's often parameterised in 3D global circulation models, as resolving convection is numerically difficult. Two convection parameterisations are included within the model: convective adjustment (CA), and mixing length theory (MLT). ","category":"page"},{"location":"model/","page":"Description","title":"Description","text":"CA forcibly adjusts a convectively unstable region of the atmosphere to the corresponding adiabat, while ensuring that enthalpy is conserved. Only dry adjustment is included in the model. This does not allow the convective energy fluxes to be calculated directly.","category":"page"},{"location":"model/","page":"Description","title":"Description","text":"MLT directly calculates the energy flux associated with convective heat transport, and thus is the preferred parameterisation within the model. It assumes that parcels of gas are diffused over a characteristic mixing length, transporting energy in the process.","category":"page"},{"location":"model/","page":"Description","title":"Description","text":"Heat capacities are temperature-dependent, calculated using the Shomate Equation with coefficients derived from the NIST website.","category":"page"},{"location":"model/#Solar-flux","page":"Description","title":"Solar flux","text":"","category":"section"},{"location":"model/","page":"Description","title":"Description","text":"The radiation component requires two boundary conditions the energy. The first is the shortwave downward-directed flux from the star at the top of the atmosphere. This is quantified by the instellation, a scale factor, a grey bond albedo, and the solar zenith angle. All of these may be provided to the model through the configuration file.","category":"page"},{"location":"model/#Solution-types","page":"Description","title":"Solution types","text":"","category":"section"},{"location":"model/","page":"Description","title":"Description","text":"Depending on the system you wish to model, it is necessary to tell AGNI what kind of solution to solve for. There are currently a few options available set by the solution_type variable.   ","category":"page"},{"location":"model/","page":"Description","title":"Description","text":"(0) Aim to conserve energy fluxes throughout the column. The surface temperature is set at T_s assuming blackbody emission with a fixed surface albedo. The bottom-most temperature value in the column is extrapolated from the rest of the profile.\n(1) Same as 0, but the bottom-most temperature value is fixed equal to T_s.\n(2) Aim to conserve energy fluxes throughout the column. The surface temperature is set by energy transport through a solid conductive boundary layer (CBL) such that T_s = T_m - Fdk, where T_m is the interior mantle temperature, while k and d are material properties. \n(3) Solve for a state such that the flux carried at each level is equal to sigma T_texteff^4. In this case, T_texteff represents the rate at which a planet is losing energy into space. ","category":"page"},{"location":"model/#Obtaining-a-solution","page":"Description","title":"Obtaining a solution","text":"","category":"section"},{"location":"model/","page":"Description","title":"Description","text":"AGNI is designed for modelling planetary atmospheres with high surface pressures and temperatures. This means that the radiative timescale differs by several orders of magnitude across the column, which makes obtaining a solution difficult. The model contains a suite of methods for obtaining a solution.","category":"page"},{"location":"model/#A)-Solving-a-non-linear-system","page":"Description","title":"A) Solving a non-linear system","text":"","category":"section"},{"location":"model/","page":"Description","title":"Description","text":"To obtain a temperature structure solution that conserves energy more precisely than a time-stepping method, it is possible to construct the model as a system of non-linear equations vecr(vecx). This algorithm obtains the roots of the nonlinear system formed by the flux divergence r_i at each level i, with cell-centre temperatures used as the independent variables x_i. Finite-difference methods are used to estimate the jacobian matrix in this case. Similar methods have been used in a handful of planetary radiative-convective models (e.g. ATMO), but is more commonly used by the stellar physics community. Currently implemented solvers: Newton-Raphson, Gauss-Newton, and Levenberg-Marquardt. Optionally, the code uses a golden-section linesearch algorithm to determine the optimal step length. ","category":"page"},{"location":"model/#B)-Time-stepping","page":"Description","title":"B) Time-stepping","text":"","category":"section"},{"location":"model/","page":"Description","title":"Description","text":"AGNI also implements a multistep Adams-Bashforth integrator to integrate the heating rates at each level over time. This is very robust to the initial conditions provided, but is not able to obtain an energy-conserving solution very quickly. Each level to evolves on its own timescale, which provides unphysical intermediate solutions but a physical final state. Time-stepping methods such as these are used in other radiative-convective models (with their own enhancements) such as HELIOS and Exo_k. This method is not the primary solver within AGNI, so it in \"maintainence mode\" rather than active development.","category":"page"},{"location":"model/#Other-features","page":"Description","title":"Other features","text":"","category":"section"},{"location":"model/","page":"Description","title":"Description","text":"AGNI can calculate emission spectra, provided with T(p) and the volume mixing ratios of the gases. This is performed using the same RT as the RCE calculations, so is limited in resolution by the choice of correlated-k bands. Similarly, the the contribution function can also be calculated.","category":"page"},{"location":"model/#Julia-and-Fortran","page":"Description","title":"Julia and Fortran","text":"","category":"section"},{"location":"model/","page":"Description","title":"Description","text":"AGNI is primarily written in Julia, while SOCRATES itself is written in Fortran. Julia was chosen because it allows the SOCRATES binaries to be included in the precompiled code, which significantly improves performance.","category":"page"},{"location":"manual/#Manual","page":"Manual","title":"Manual","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"This page documents functions and data structures within the code.","category":"page"},{"location":"manual/#Index","page":"Manual","title":"Index","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"","category":"page"},{"location":"examples/#Example-outputs","page":"Example outputs","title":"Example outputs","text":"","category":"section"},{"location":"examples/#Pure-steam-runaway-greenhouse-effect","page":"Example outputs","title":"Pure steam runaway greenhouse effect","text":"","category":"section"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"By assuming the atmosphere temperature profile follows a dry adiabat and the water vapour-condensate coexistance curve defined by the Clausius-Claperyron relation, we see a characteristic relationship between the outgoing longwave radiation (OLR) and the surface temperature (T_s). Initially OLR increases with T_s, but as the condensing layer (which is independent of T_s) overlaps with the photosphere, OLR and T_s decouple. Eventually the atmosphere reaches a dry post-runaway state, and OLR increases rapidly with T_s.","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"  <img src=\"assets/runaway/curve.png\" width=70% class=\"center\"/>","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"You can find a Jupyter notebook which reproduces this result in the tutorials directory of the repository.","category":"page"},{"location":"examples/#Prescribed-convective-case","page":"Example outputs","title":"Prescribed convective case","text":"","category":"section"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"In this case, a temperature profile is prescribed to follow a dry adiabat from the surface to a moist region, and then a pseudoadiabat to the top of the atmosphere. This is in line with previous works and the OLR curve above. ","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"  <img src=\"assets/nosolve/plot_ptprofile.png\" width=60% class=\"center\"/> ","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"Radiative fluxes are then calculated according to this  temperature profile. Because the profile is prescribed, the fluxes are not balanced locally or globally across the column.","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"  <img src=\"assets/nosolve/plot_fluxes.png\" width=60% class=\"center\"/> ","category":"page"},{"location":"examples/#Radiative-convective-solution","page":"Example outputs","title":"Radiative-convective solution","text":"","category":"section"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"Instead, we can model an atmosphere such that that energy is globally and locally conserved. Convection is parameterised using mixing length theory in this case, allowing the system to be solved using a Newton-Raphson method. In the convective region at ~0.1 bar, we can see that the radiative fluxes and convective fluxes entirely cancel, because AGNI was asked to solve for a case with zero total flux transport. ","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"  <img src=\"assets/withsolve/plot_ptprofile.png\" width=60% class=\"center\"/> \n  <br />\n  <img src=\"assets/withsolve/plot_fluxes.png\" width=60% class=\"center\"/> ","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"We can also plot the outgoing emission spectrum and normalised contribution function. The spectrum clearly demonstrates complex water absorption features, and exceeds blackbody emission at shorter wavelengths due to Rayleigh scattering. The normalised contribution function quantifies how much each pressure level contributes to the outgoing emission spectrum at a given wavelength – this is then plotted versus wavelength and pressure. ","category":"page"},{"location":"examples/","page":"Example outputs","title":"Example outputs","text":"  <img src=\"assets/withsolve/plot_emission.png\" width=90% class=\"center\"/> \n  <br />\n  <img src=\"assets/withsolve/plot_contfunc.png\" width=90% class=\"center\"/> ","category":"page"},{"location":"#AGNI","page":"Home","title":"AGNI","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A computer code for modelling the atmospheres of hot terrestrial (exo)planets. AGNI's primary purpose is to solve for the atmospheric temperature structure at radiative-convective equilibrium while coupling to an interior model.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The model itself is currently proprietary. It will be made available under the BSD 3-Clause Clear License once published in a journal. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Follow Getting started for information on installing the code and obtaining results.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pronounced: ag-nee. Named after the fire deity of Hinduism.","category":"page"},{"location":"setup/#Getting-started","page":"Getting started","title":"Getting started","text":"","category":"section"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"This page outlines requirements and installation steps for the code.","category":"page"},{"location":"setup/#Requirements","page":"Getting started","title":"Requirements","text":"","category":"section"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"Julia (NB: install only from julialang.org - do not use your system package manager)\nSOCRATES","category":"page"},{"location":"setup/#Supported-platforms","page":"Getting started","title":"Supported platforms","text":"","category":"section"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"MacOS (ARM and x86-64)\nGNU/Linux (x86-64)","category":"page"},{"location":"setup/#Installation","page":"Getting started","title":"Installation","text":"","category":"section"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"Setup SOCRATES by doing either ONE of the following...\nFollow the instructions on the SOCRATES GitHub page\nRun ./get_socrates.sh\n$ julia\njulia> ] \n(@v1.10) pkg> activate . ← note the dot\n(AGNI) pkg> build\nPress backspace\njulia> exit() ","category":"page"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"AGNI is now installed as a package into a Julia environment in this directory.    You should run the tests next.","category":"page"},{"location":"setup/#Testing","page":"Getting started","title":"Testing","text":"","category":"section"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"$ julia\njulia> ]\n(@v1.10) pkg> activate .\n(AGNI) pkg> test","category":"page"},{"location":"setup/","page":"Getting started","title":"Getting started","text":"These tests may trigger recompilation of depenencies.","category":"page"}]
}
