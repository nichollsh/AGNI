# Your first calculation, CLI

So you want to run some standalone climate calculations? Let's start with a simple steam atmosphere with a temperature profile 'prescribed' to follow an adiabat.

Run AGNI with the default configuration file:
```bash
./agni.jl res/config/default.toml
```

You should see the following output:
```log
[ INFO  ] Using configuration 'Default'
[ INFO  ] Setting-up a new atmosphere struct
[ INFO  ] Loading thermodynamic data
[ INFO  ] Inserting stellar spectrum, Rayleigh coefficients
[ INFO  ] Allocating atmosphere with initial composition:
[ INFO  ]       1 H2O      1.00e+00 (EOS_AQUA)
[ INFO  ] Setting T(p): dry, sat
[ INFO  ] Solving with 'none'
[ INFO  ]     done
[ INFO  ] Radiative transfer benchmarking statistics...
[ INFO  ]     total evals:   2, over 0.449157808 secs
[ INFO  ]     time per eval: 224.578904 ms/eval
[ INFO  ] Photosphere defined at τ=0.02 and λ=1.00 μm
[ INFO  ]     pressure: 82.00 mbar, temperature: 370.84 K
[ INFO  ]     planet radius: 1.07 R⊕, bulk density: 4.46e+03 kg/m^3
[ INFO  ] Writing results
[ INFO  ] Plotting results
[ INFO  ] Model runtime: 18.63 seconds

```

We can see a few things here. Firstly, AGNI tells us which configuration file it is using and starts the setup.

The line following "Allocating atmosphere with composition" is a table of gases, their
volume mixing ratios, and flags. In this case there is only one gas. Potential flags for each gas species are:
* `EOS_[XX]` - using the `[XX]` equation of state (e.g. ideal gas, AQUA)
* `NO_OPACITY` - no opacity data available, but can contribute to the thermodynamics
* `NO_THERMO` - no thermodynamic data available, so will be treated as a diatomic ideal gas
* `COND` - this gas is allowed to condense

Note that 'gas' here refers to a chemical species that generally, but not necessarily exists in as gas phase, so it could be supercritical or condensed.

The temperature profile is set the the **dry** adiabat with a **sat**urated upper atmosphere.

No radiative-convective solution is requested here, so there is no solver to be called. Hence, the model outputs "Solving with none".

Output files are written to the directory specified in the configuration file (default: `out/`). This includes a NetCDF data file and various plots.

![](fig_nosolve_ptprofile.png)

Radiative fluxes are then calculated according to this temperature profile. Because the
profile is prescribed, the fluxes are not balanced locally or globally across the column.

![](fig_nosolve_fluxes.png)

This calculation does not solve for radiative-convective equilibrium. In the next tutorial, we will calculate climate structures
in a manner which self-consistently conserves energy.

