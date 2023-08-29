#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# THEMIS main file, for standalone execution
# -------------

println("Begin THEMIS")

using Revise
using Plots

# Include local jl files
include("socrates/julia/src/SOCRATES.jl")

push!(LOAD_PATH, joinpath(pwd(),"src"))
import atmosphere
import setup_pt
import plotting 
import solver


# Configuration options
tstar           = 3000.0    # LW uflux bottom boundary condition [kelvin]
zenith_degrees  = 45.53     # Zenith angle [degrees from zenith]
toa_heating     = 4.451e+04 # SW dflux top boundary condition [W m-2]
gravity         = 9.81
nlev_centre     = 100
p_surf          = 100.0     # bar
p_top           = 1e-7      # bar 
mixing_ratios   = Dict([("CO", 0.8), ("H2O", 0.1), ("N2", 0.1)])

# spectral_file = "socrates/data/spectra/ga7/sp_lw_ga7" 
spectral_file = "res/runtime_spectral_file_rscat"
output_dir = "out/"

# Create output direct
rm(output_dir,force=true,recursive=true)
mkdir(output_dir)

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, spectral_file, true,
                        false, true, false, false,
                        zenith_degrees, toa_heating, tstar,
                        gravity, nlev_centre, p_surf, p_top,
                        mixing_ratios
                        )
atmosphere.allocate!(atmos)

# Set to dry adiabat 
# setup_pt.dry_adiabat!(atmos)
# setup_pt.condensing!(atmos, "H2O")

# Calculate LW and SW fluxes (once)
# atmosphere.radtrans!(atmos, true)
# atmosphere.radtrans!(atmos, false)

# println("OLR = $(atmos.flux_u_lw[1]) W m-2")

# Call solver 
solver.solve_energy!(atmos, surf_state=0, plot=true)

# Save result
atmosphere.write_pt(atmos, joinpath(output_dir,"pt.csv"))
plotting.plot_pt(atmos, joinpath(output_dir,"pt.pdf"))
plotting.plot_fluxes(atmos, joinpath(output_dir,"fluxes.pdf"))

# Deallocate
atmosphere.deallocate!(atmos)

println("Goodbye")
