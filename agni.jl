#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# AGNI main file for standalone execution
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
tstar           = 1740.0    # LW uflux bottom boundary condition [kelvin]
toa_heating     = 4.391e+04 # SW dflux top boundary condition [W m-2]
gravity         = 10.0
nlev_centre     = 100
p_surf          = 280.4     # bar
p_top           = 1e-6      # bar 
mixing_ratios   = Dict([
                        ("H2O" , 10.36251/280.4),
                        ("CO2" , 94.02385/280.4),
                        ("H2"  , 0.54123/280.4),
                        ("CO"  , 172.91992/280.4),
                        ("N2"  , 2.28491/280.4)
                        ])

spfile_name   = "Mallard"
star_file     = "res/stellar_spectra/trappist-1.txt"
output_dir    = "out/"

# Create output direct
rm(output_dir,force=true,recursive=true)
mkdir(output_dir)

# Setup atmosphere
println("Atmosphere: setting up")
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, tstar,
                         gravity, nlev_centre, p_surf, p_top,
                         mixing_ratios,
                         flag_gcontinuum=true,
                         flag_rayleigh=true
                         )
atmosphere.allocate!(atmos;stellar_spectrum=star_file,spfile_noremove=true)

# Set PT profile 
println("Atmosphere: setting initial T(p)")
setpt.prevent_surfsupersat!(atmos)
setpt.dry_adiabat!(atmos)
setpt.stratosphere!(atmos, 300.0)

# Calculate LW and SW fluxes (once)
println("RadTrans: calculating fluxes")
atmosphere.radtrans!(atmos, true)
atmosphere.radtrans!(atmos, false)

# Call solver 
solver.solve_energy!(atmos, surf_state=1, modplot=1, verbose=true)

# Save result
atmosphere.write_pt(atmos,  joinpath(atmos.OUT_DIR,"pt.csv"))
plotting.plot_pt(atmos,     joinpath(atmos.OUT_DIR,"pt.pdf"))
plotting.anim_solver(atmos)
plotting.plot_fluxes(atmos, joinpath(atmos.OUT_DIR,"fluxes.pdf"))

# Deallocate atmosphere
println("Atmosphere: deallocating arrays")
atmosphere.deallocate!(atmos)

println("Goodbye")
