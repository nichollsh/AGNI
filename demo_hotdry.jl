#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# Example AGNI executable, demonstrating a case similar to Figure1 of Selsis+23
# -------------

tbegin = time()
println("Begin post-runaway demo")

# Get AGNI root directory
ROOT_DIR = dirname(abspath(@__FILE__))
ENV["GKSwstype"] = "100"

# Include libraries
using Revise

# Include local jl files
include("socrates/julia/src/SOCRATES.jl")
push!(LOAD_PATH, joinpath(ROOT_DIR,"src"))
import atmosphere
import setpt
import plotting 
import phys


# Configuration options
tstar           = 2490.0    # Surface temperature [kelvin]
instellation    = 1361.0 * 35.6 
s0_factor       = 3.0/8.0
albedo          = 0.18
zenith          = 48.19
radius          = 6.37e6    # metres
gravity         = 9.81      # m s-2
nlev_centre     = 50  
p_surf          = 270.0    # bar
p_top           = 1e-5      # bar 
mf_dict         = Dict([
                        ("H2O" , 1.0),
                        ])

spfile_name   = "res/spectral_files/Oak/Oak"
star_file     = "res/stellar_spectra/sun.txt"
output_dir    = "out/"

# Setup atmosphere
println("Setting up")
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         instellation, s0_factor, albedo, zenith,
                         tstar,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict,
                         flag_gcontinuum=true,
                         flag_rayleigh=true,
                         overlap_method=4,
                         tmp_floor=2.0,
                         thermo_functions=true,
                 )
atmosphere.allocate!(atmos;stellar_spectrum=star_file,spfile_noremove=true)

# Set PT profile 
println("Setting initial T(p)")
setpt.isothermal!(atmos, tstar-300.0)

# Create output directory
rm(output_dir,force=true,recursive=true)
if !isdir(output_dir) && !isfile(output_dir)
    mkdir(output_dir)
end

println("Running model...")

# Call solver(s)
dry_convect = true
surf_state  = 0

import solver_tstep
solver_tstep.solve_energy!(atmos, surf_state=surf_state, modplot=100, modprop=5, verbose=true, 
                            dry_convect=dry_convect, 
                            accel=true, rtol=1.0e-4, atol=1.0e-2,
                            max_steps=500, min_steps=200, use_mlt=true,
                            dt_max=150.0, F_losspct_conv=1.0)


import solver_nlsol
solver_nlsol.solve_energy!(atmos, surf_state=surf_state,
                            dry_convect=dry_convect,
                            max_steps=100)

# Write arrays
atmosphere.write_pt(atmos,      joinpath(atmos.OUT_DIR,"pt.csv"))
atmosphere.write_ncdf(atmos,    joinpath(atmos.OUT_DIR,"atm.nc"))
atmosphere.write_fluxes(atmos,  joinpath(atmos.OUT_DIR,"fl.csv"))

# Save plots
plotting.plot_x(atmos,      joinpath(atmos.OUT_DIR,"mf.pdf"))
plotting.plot_pt(atmos,     joinpath(atmos.OUT_DIR,"pt.pdf"))
plotting.plot_fluxes(atmos, joinpath(atmos.OUT_DIR,"fl.pdf"))

# Deallocate atmosphere
println("Deallocating arrays")
atmosphere.deallocate!(atmos)

runtime = round(time() - tbegin, digits=2)
println("Runtime: $runtime seconds")
println("Goodbye")
