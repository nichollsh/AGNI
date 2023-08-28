#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# THEMIS main file, for standalone execution
# -------------

# Include local jl files
include("socrates/julia/src/SOCRATES.jl")
include("src/atmosphere.jl")
include("src/setup_pt.jl")
include("src/plotting.jl")


# Configuration options
tstar           = 279.39    # LW uflux bottom boundary condition [kelvin]
zenith_degrees  = 45.53     # Zenith angle [degrees from zenith]
toa_heating     = 1381.53   # SW dflux top boundary condition [W m-2]
gravity         = 9.81
nlev_centre     = 100
p_surf          = 1e+2      # bar
p_top           = 1e-5      # bar 
mixing_ratios   =  Dict([("H2O", 0.7), ("CO2", 0.3)])

# spectral_file = "socrates/data/spectra/ga7/sp_lw_ga7" 
spectral_file = "res/runtime_spectral_file_rscat"
lw = true
output_dir = "out/"

# Create output direct
rm(output_dir,force=true,recursive=true)
mkdir(output_dir)

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup_atmos!(atmos, spectral_file, true,
                        false, false, false, false,
                        zenith_degrees, toa_heating, tstar,
                        gravity, nlev_centre, p_surf, p_top,
                        mixing_ratios
                        )
atmosphere.alloc_atmos!(atmos)

# Set to dry adiabat 
setup_pt.dry_adiabat!(atmos)

# Calculate LW and SW fluxes
atmosphere.calc_fluxes!(atmos, true)
atmosphere.calc_fluxes!(atmos, false)

# Save result
atmosphere.write_pt(atmos, "out/pt.csv")

# Plot result
plotting.plot_pt(atmos, output_dir)
plotting.plot_fluxes(atmos, output_dir)


println("Done!")
