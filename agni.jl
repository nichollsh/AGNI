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
import solver_accel
# import solver_cvode
import phys


# Configuration options
tstar           = 1960.0    # Surface temperature [kelvin]
toa_heating     = 510.0     # Instellation flux [W m-2]
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
                         flag_rayleigh=false,
                         overlap_method=4,
                         zenith_degrees=48.19,
                         skin_d=0.01,
                         skin_k=2.0,
                         tmp_magma=2500.0,
                         tmp_floor=2.0,
                         res_switching=false
                 )
atmosphere.allocate!(atmos;stellar_spectrum=star_file,spfile_noremove=true)

# Set PT profile 
println("Setting initial T(p)")
# setpt.isothermal!(atmos, 1300.0)
# setpt.fromcsv!(atmos,"pt.csv")
# setpt.prevent_surfsupersat!(atmos)
setpt.dry_adiabat!(atmos)
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

# Call solver 
solver_accel.solve_energy!(atmos, surf_state=0, modplot=10, verbose=true, 
                            dry_convect=true, accel=true, extrap=false,
                            max_steps=3000, min_steps=50, use_mlt=false,
                            dt_max=200.0, F_losspct_conv=0.1)

# solver_cvode.solve_energy!(atmos, surf_state=2, verbose=true, dry_convect=true,  max_steps=500)

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
