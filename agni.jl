#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# AGNI executable file for standalone execution
# -------------

tbegin = time()
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
import solver_euler
import solver_cvode
import phys


# Configuration options
tstar           = 2100.0    # Surface temperature [kelvin]
toa_heating     = 3.772e+04 # Instellation flux [W m-2]
radius          = 7.12e6    # metres
gravity         = 10.8      # m s-2
nlev_centre     = 128  
p_surf          = 127.0     # bar
p_top           = 1e-6      # bar 
mf_dict         = Dict([
                        ("H2O" , 0.91805),
                        ("CO2" , 5.98710),
                        ("H2"  , 2.37994),
                        ("CO"  , 115.89786),
                        ("N2"  , 1.77739)
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
                         flag_rayleigh=false,
                         overlap_method=4,
                         zenith_degrees=54.4,
                         skin_d=0.02,
                         skin_k=2.0,
                         tmp_magma=2654.0
                 )
atmosphere.allocate!(atmos;stellar_spectrum=star_file,spfile_noremove=true)

# Set PT profile 
println("Setting initial T(p)")
# setpt.fromcsv!(atmos,"out/pt.csv")
# setpt.prevent_surfsupersat!(atmos)
# setpt.dry_adiabat!(atmos)
# setpt.stratosphere!(atmos, 300.0)

# Create output directory
rm(output_dir,force=true,recursive=true)
if !isdir(output_dir) && !isfile(output_dir)
    mkdir(output_dir)
end

# Calculate LW and SW fluxes (once)
# println("RadTrans: calculating fluxes")
# atmosphere.radtrans!(atmos, true)
# atmosphere.radtrans!(atmos, false)

# Calculate convective fluxes (once)
# println("MLT: calculating fluxes")
# atmosphere.mlt!(atmos)

# Call solver 
println("Starting solver")
solver_euler.solve_energy!(atmos, surf_state=2, modplot=2, verbose=true, dry_convect=true, max_steps=300, min_steps=20, use_mlt=true)
solver_cvode.solve_energy!(atmos, surf_state=2,            verbose=true, dry_convect=true, max_steps=300)

# Write arrays
atmosphere.write_ncdf(atmos,    joinpath(atmos.OUT_DIR,"atm.nc"))
atmosphere.write_pt(atmos,      joinpath(atmos.OUT_DIR,"pt.csv"))
atmosphere.write_fluxes(atmos,  joinpath(atmos.OUT_DIR,"fl.csv"))

# Save plots
plotting.plot_x(atmos,      joinpath(atmos.OUT_DIR,"mf.pdf"))
plotting.plot_pt(atmos,     joinpath(atmos.OUT_DIR,"pt.pdf"))
plotting.plot_fluxes(atmos, joinpath(atmos.OUT_DIR,"fl.pdf"))
plotting.anim_solver(atmos)

# Deallocate atmosphere
println("Deallocating arrays")
atmosphere.deallocate!(atmos)

runtime = round(time() - tbegin, digits=2)
println("Runtime: $runtime seconds")
println("Goodbye")
