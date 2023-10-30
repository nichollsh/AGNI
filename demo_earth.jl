#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# AGNI executable file for standalone execution
# -------------

tbegin = time()
println("Begin Earth demo")

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
import phys


# Configuration options
tstar           = 300.0     # LW uflux bottom boundary condition [kelvin]
toa_heating     = 1367.0    # Daytime instellation flux [W m-2]
radius          = 6.37e6    # metres
gravity         = 9.81      # m s-2
nlev_centre     = 150  
p_surf          = 1.0       # bar
p_top           = 1e-7      # bar 
mf_path         = "doc/example_earth/equ.csv"

spfile_name   = "res/spectral_files/Reach/Reach"
star_file     = "res/stellar_spectra/sun.txt"
output_dir    = "out/"

# Create output directory
rm(output_dir,force=true,recursive=true)
if !isdir(output_dir) && !isfile(output_dir)
    mkdir(output_dir)
end

# Setup atmosphere
println("Setting up")
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, tstar,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_path=mf_path,
                         flag_gcontinuum=true,
                         flag_rayleigh=false,  # the RFM profile does not include scattering
                         zenith_degrees=30.0,
                         overlap_method=4
                 )
atmosphere.allocate!(atmos;stellar_spectrum=star_file)

# Call solver 
println("Starting solver")
solver_euler.solve_energy!(atmos, surf_state=1, modplot=2, verbose=true, extrap=false,
                        dry_convect=true, h2o_convect=true, sens_heat=true, 
                        dt_max=2.0, max_steps=300)

# Write arrays
atmosphere.write_pt(atmos,      joinpath(atmos.OUT_DIR,"pt.csv"))
atmosphere.write_fluxes(atmos,  joinpath(atmos.OUT_DIR,"fl.csv"))

# Save plots
plotting.plot_x(atmos,      joinpath(atmos.OUT_DIR,"mf.pdf"))
plotting.plot_pt(atmos,     joinpath(atmos.OUT_DIR,"pt.pdf"))
plotting.plot_fluxes(atmos, joinpath(atmos.OUT_DIR,"fl.pdf"))

# Deallocate atmosphere
println("Deallocating arrays")
atmosphere.deallocate!(atmos)

println("Goodbye")
