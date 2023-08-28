#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# THEMIS main file, for standalone execution
# -------------

# Include local jl files
include("socrates/julia/src/SOCRATES.jl")
include("src/radtrans.jl")
include("src/solver.jl")
include("src/plotting.jl")


# Configuration options
tstar           = 279.39    # LW uflux bottom boundary condition [kelvin]
zenith_degrees  = 45.53     # Zenith angle [degrees from zenith]
toa_heating     = 1381.53   # SW dflux top boundary condition [W m-2]
gravity         = 9.81
nlev_centre     = 69

all_channels = true
# spectral_file = "socrates/data/spectra/ga7/sp_lw_ga7" 
spectral_file = "res/runtime_spectral_file_rscat"
lw = true
output_dir = "out/"

# Create output direct
rm(output_dir,force=true,recursive=true)
mkdir(output_dir)

# Allocate atmos object
atmos = radtrans.Atmos_t()

radtrans.setup_atmos!(atmos, spectral_file, all_channels,
                        false, false, false, false,
                        zenith_degrees, toa_heating, tstar,
                        gravity, nlev_centre
                        )

radtrans.alloc_atmos!(atmos)

radtrans.calc_fluxes!(atmos, true)
radtrans.calc_fluxes!(atmos, false)


plotting.plot_pt(atmos, output_dir)
plotting.plot_fluxes(atmos, output_dir)


println("Done!")
