#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# AGNI executable file for standalone execution
# -------------

tbegin = time()
println("Begin AGNI")

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
tstar           = 2450.0    # Surface temperature [kelvin]
toa_heating     = 3.0/8.0 * 1361.0 * (1-0.18) * 35.6
radius          = 6.37e6    # metres
gravity         = 9.81      # m s-2
nlev_centre     = 50  
p_surf          = 270.0    # bar
p_top           = 1e-5      # bar 
mf_dict         = Dict([
                        ("H2O" , 1.0),
                        # ("CO2" , 0.1),
                        # ("H2" , 1.0),
                        # ("CO" , 90.58514),
                        # ("N2" , 1.41003)
                        ])

spfile_name   = "res/spectral_files/Oak/Oak"
star_file     = "res/stellar_spectra/sun.txt"
output_dir    = "out/"

# Setup atmosphere
println("Setting up")
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, tstar,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict,
                         flag_gcontinuum=true,
                         flag_rayleigh=true,
                         flag_cloud=false,
                         overlap_method=4,
                         zenith_degrees=48.19,
                         skin_d=0.01,
                         skin_k=2.0,
                         tmp_magma=2500.0,
                         tmp_floor=2.0,
                         thermo_functions=true,
                 )
atmosphere.allocate!(atmos;stellar_spectrum=star_file,spfile_noremove=true)

# Set PT profile 
println("Setting initial T(p)")
# setpt.fromcsv!(atmos,"pt.csv")
setpt.isothermal!(atmos, tstar-300.0)
# setpt.prevent_surfsupersat!(atmos)
# setpt.dry_adiabat!(atmos)
# setpt.condensing!(atmos, "H2O")
# setpt.stratosphere!(atmos, 500.0)

# Create output directory
rm(output_dir,force=true,recursive=true)
if !isdir(output_dir) && !isfile(output_dir)
    mkdir(output_dir)
end

atmosphere.write_pt(atmos, joinpath(atmos.OUT_DIR,"pt_ini.csv"))

println("Running model...")

# Calculate LW and SW fluxes (once)
# atmosphere.radtrans!(atmos, true)
# atmosphere.radtrans!(atmos, false)

# Calculate convective fluxes (once)
# println("MLT: calculating fluxes")
# atmosphere.mlt!(atmos)


# Call solver(s)
dry_convect = true
surf_state  = 0

import solver_tstep
solver_tstep.solve_energy!(atmos, surf_state=surf_state, modplot=10, modprop=5, verbose=true, 
                            dry_convect=dry_convect,
                            accel=true, rtol=1.0e-4, atol=1.0e-2,
                            max_steps=600, min_steps=200, use_mlt=true,
                            dt_max=100.0, F_losspct_conv=1.0)


import solver_nlsol
solver_nlsol.solve_energy!(atmos, surf_state=surf_state, dry_convect=dry_convect, max_steps=100)

# Write arrays
atmosphere.write_pt(atmos,      joinpath(atmos.OUT_DIR,"pt.csv"))
atmosphere.write_ncdf(atmos,    joinpath(atmos.OUT_DIR,"atm.nc"))
atmosphere.write_fluxes(atmos,  joinpath(atmos.OUT_DIR,"fl.csv"))

# Save plots
plotting.anim_solver(atmos)
plotting.plot_x(atmos,      joinpath(atmos.OUT_DIR,"mf.pdf"))
plotting.plot_pt(atmos,     joinpath(atmos.OUT_DIR,"pt.pdf"))
plotting.plot_fluxes(atmos, joinpath(atmos.OUT_DIR,"fl.pdf"))

# Deallocate atmosphere
println("Deallocating arrays")
atmosphere.deallocate!(atmos)

runtime = round(time() - tbegin, digits=2)
println("Runtime: $runtime seconds")
println("Goodbye")
