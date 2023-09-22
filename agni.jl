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
import solver
import phys


# Configuration options
tstar           = 2500.0    # LW uflux bottom boundary condition [kelvin]
toa_heating     = 1000.0 # SW dflux top boundary condition [W m-2]
radius          = 6.0e6     # metres
gravity         = 9.81      # m s-2
nlev_centre     = 100  
p_surf          = 270.0     # bar
p_top           = 1e-6      # bar 
mixing_ratios   = Dict([
                        ("H2O" , 1.0),
                        ("CO2" , 0.0),
                        ("CO"  , 0.0),
                        ])

spfile_name   = "res/spectral_files/Oak/Oak"
star_file     = "res/stellar_spectra/sun.txt"
output_dir    = "out/"

# Create output directory
if !isdir(output_dir) && !isfile(output_dir)
    mkdir(output_dir)
end

# Setup atmosphere
println("Atmosphere: setting up")
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, tstar,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mixing_ratios,
                         flag_gcontinuum=true,
                         flag_rayleigh=false
                         )
atmosphere.allocate!(atmos;stellar_spectrum=star_file,spfile_noremove=true)

# Set PT profile 
println("Atmosphere: setting initial T(p)")
# setpt.fromcsv!(atmos,"out/pt.csv")
setpt.prevent_surfsupersat!(atmos)
setpt.dry_adiabat!(atmos)
setpt.stratosphere!(atmos, 300.0)

# Calculate LW and SW fluxes (once)
#println("RadTrans: calculating fluxes")
#atmosphere.radtrans!(atmos, true)
#atmosphere.radtrans!(atmos, false)

# Call solver 
println("RadTrans: starting solver")
solver.solve_energy!(atmos, surf_state=2, modplot=1, verbose=true, dry_adjust=true, h2o_adjust=true, max_steps=500, min_steps=40)

# Write arrays
atmosphere.write_ncdf(atmos,   joinpath(atmos.OUT_DIR,"atm.nc"))
atmosphere.write_pt(atmos,      joinpath(atmos.OUT_DIR,"pt.csv"))
atmosphere.write_fluxes(atmos,  joinpath(atmos.OUT_DIR,"fl.csv"))

# Save plots
plotting.plot_pt(atmos,     joinpath(atmos.OUT_DIR,"pt.pdf"))
plotting.plot_fluxes(atmos, joinpath(atmos.OUT_DIR,"fl.pdf"))
plotting.anim_solver(atmos)

# Deallocate atmosphere
println("Atmosphere: deallocating arrays")
atmosphere.deallocate!(atmos)

runtime = round(time() - tbegin, digits=2)
println("Runtime: $runtime seconds")
println("Goodbye")
