#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# THEMIS main file, for standalone execution
# -------------

println("Begin THEMIS")

using Revise

# Include local jl files
include("socrates/julia/src/SOCRATES.jl")

push!(LOAD_PATH, joinpath(pwd(),"src"))
import atmosphere
import setpt
import plotting 
import solver
import phys


# Configuration options
tstar           = 1700.0    # LW uflux bottom boundary condition [kelvin]
zenith_degrees  = 45.53     # Zenith angle [degrees from zenith]
toa_heating     = 4.451e+04 # SW dflux top boundary condition [W m-2]
gravity         = 9.81
nlev_centre     = 100
p_surf          = 300.0     # bar
p_top           = 1e-7      # bar 
mixing_ratios   = Dict([
                        ("CO" , 0.90),
                        ("CO2", 0.05),
                        ("N2" , 0.05)
                        ])

spectral_file = "res/spectral_files/Mallard/Mallard"
star_file     = "res/stellar_spectra/spec_sun.txt"
output_dir = "out/"

# Create output direct
rm(output_dir,force=true,recursive=true)
mkdir(output_dir)

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, spectral_file,
                         star_file, zenith_degrees, 
                         toa_heating, tstar,
                         gravity, nlev_centre, p_surf, p_top,
                         mixing_ratios,
                         flag_gcontinuum=true,
                         flag_rayleigh=true
                         )
atmosphere.allocate!(atmos)

# Set PT profile 
setpt.dry_adiabat!(atmos)
# setpt.condensing!(atmos, "H2O")

# Calculate LW and SW fluxes (once)
atmosphere.radtrans!(atmos, true)
atmosphere.radtrans!(atmos, false)

# Call solver 
solver.solve_energy!(atmos, surf_state=0, plot=true)

# Save result
atmosphere.write_pt(atmos, joinpath(output_dir,"pt.csv"))
plotting.plot_pt(atmos, joinpath(output_dir,"pt.pdf"))
plotting.plot_fluxes(atmos, joinpath(output_dir,"fluxes.pdf"))


# Deallocate atmosphere
atmosphere.deallocate!(atmos)

println("Goodbye")
