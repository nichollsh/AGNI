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
import solver_euler
import solver_cvode
import phys


# Configuration options
tstar           = 2000.0    # Surface temperature [kelvin]
toa_heating     = 3.745e+04 # Instellation flux [W m-2]
radius          = 7.12e6    # metres
gravity         = 10.8      # m s-2
nlev_centre     = 100  
p_surf          = 120.8     # bar
p_top           = 1e-6      # bar 
mf_dict         = Dict([
                        ("H2O" , 145.64182),
                        ("CO2" , 4.72241),
                        ("H2" , 7.53004),
                        ("CO" , 90.58514),
                        ("N2" , 1.41003)
                        ])

spfile_name   = "res/spectral_files/Mallard/Mallard"
star_file     = "res/stellar_spectra/trappist-1.txt"
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
                         overlap_method=4,
                         zenith_degrees=54.4,
                         skin_d=0.02,
                         skin_k=2.0,
                         tmp_magma=2700.0,
                         tmp_floor=2.0
                 )
atmosphere.allocate!(atmos;stellar_spectrum=star_file,spfile_noremove=true)

# Set PT profile 
println("Setting initial T(p)")
# setpt.fromcsv!(atmos,"out/pt_load.csv")
# setpt.prevent_surfsupersat!(atmos)
# setpt.dry_adiabat!(atmos)
# setpt.stratosphere!(atmos, 300.0)

# Create output directory
# rm(output_dir,force=true,recursive=true)
if !isdir(output_dir) && !isfile(output_dir)
    mkdir(output_dir)
end

atmosphere.write_pt(atmos, joinpath(atmos.OUT_DIR,"pt_ini.csv"))

# Calculate LW and SW fluxes (once)
# println("RadTrans: calculating fluxes")
# atmosphere.radtrans!(atmos, true)
# atmosphere.radtrans!(atmos, false)

# Calculate convective fluxes (once)
# println("MLT: calculating fluxes")
# atmosphere.mlt!(atmos)

# Call solver 
println("Starting solver")
solver_euler.solve_energy!(atmos, surf_state=2, modplot=1, verbose=true, 
                            dry_convect=true, use_mlt=false,
                            max_steps=1000, accel=true, extrap=false, rtol=1.0e-4, atol=1.0e-2)
# solver_cvode.solve_energy!(atmos, surf_state=2,            verbose=true, dry_convect=true,  max_steps=500)

# Write arrays
atmosphere.write_ncdf(atmos,    joinpath(atmos.OUT_DIR,"atm.nc"))
atmosphere.write_pt(atmos,      joinpath(atmos.OUT_DIR,"pt.csv"))
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
