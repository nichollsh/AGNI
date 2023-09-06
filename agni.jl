#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# AGNI main file, for standalone execution
# -------------

println("Begin AGNI")

# Get AGNI root directory
ROOT_DIR = dirname(abspath(@__FILE__))

# Include libraries
using Revise

# Include local jl files
include("socrates/julia/src/SOCRATES.jl")
push!(LOAD_PATH, joinpath(ROOT_DIR,"src"))
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

spfile_name   = "Mallard"
star_file     = "res/stellar_spectra/spec_sun.txt"
output_dir    = "out/"

# Create output direct
rm(output_dir,force=true,recursive=true)
mkdir(output_dir)

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name, zenith_degrees, 
                         toa_heating, tstar,
                         gravity, nlev_centre, p_surf, p_top,
                         mixing_ratios,
                         flag_gcontinuum=true,
                         flag_rayleigh=true
                         )
atmosphere.allocate!(atmos;stellar_spectrum=star_file,spfile_noremove=true)

# Set PT profile 
setpt.dry_adiabat!(atmos)

# Calculate LW and SW fluxes (once)
atmosphere.radtrans!(atmos, true)
atmosphere.radtrans!(atmos, false)

# Call solver 
solver.solve_energy!(atmos, surf_state=1, modplot=1)

# Save result
atmosphere.write_pt(atmos,  joinpath(atmos.OUT_DIR,"pt.csv"))
plotting.plot_pt(atmos,     joinpath(atmos.OUT_DIR,"pt.pdf"))
plotting.plot_solver(atmos, joinpath(atmos.OUT_DIR,"pt_hr.pdf"))
plotting.plot_fluxes(atmos, joinpath(atmos.OUT_DIR,"fluxes.pdf"))

# Deallocate atmosphere
atmosphere.deallocate!(atmos)

println("Goodbye")
