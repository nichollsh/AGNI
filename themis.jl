#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# THEMIS main file, for standalone execution
# -------------

println("Begin THEMIS")

using Revise
using Printf
using Plots

# Include local jl files
include("socrates/julia/src/SOCRATES.jl")

push!(LOAD_PATH, joinpath(pwd(),"src"))
import atmosphere
import setup_pt
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
mixing_ratios   = Dict([("H2O", 1.0)])

# spectral_file = "socrates/data/spectra/ga7/sp_lw_ga7" 
spectral_file = "res/runtime_spectral_file_rscat"
output_dir = "out/"

# Create output direct
rm(output_dir,force=true,recursive=true)
mkdir(output_dir)

# Read T(p) from CSV 
# read_atmos = setup_pt.readcsv("res/tptest.csv")
# p_surf = read_atmos["p_surf"]
# tstar  = read_atmos["T_surf"]

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, spectral_file,
                         zenith_degrees, toa_heating, tstar,
                         gravity, nlev_centre, p_surf, p_top,
                         mixing_ratios,
                         flag_gcontinuum=true
                         )
atmosphere.allocate!(atmos)

run_len = 50
tsurf_arr = range(100,stop=2400,length=run_len)
olr_arr = zeros(Float64, run_len)
for i in 1:run_len

    # Set PT profile 
    atmos.tstar = tsurf_arr[i]
    setup_pt.dry_adiabat!(atmos)
    setup_pt.condensing!(atmos, "H2O")

    # Calculate LW and SW fluxes (once)
    atmosphere.radtrans!(atmos, true)
    atmosphere.radtrans!(atmos, false)

    olr = atmos.flux_u_lw[1]
    @printf("Tsurf = %4.1f K  ,  OLR = %5.1f W m-2 \n",tsurf_arr[i],olr)
    olr_arr[i] = olr 

end

plt = plot(tsurf_arr,olr_arr)
savefig(plt, "runaway.pdf")

# Call solver 
# solver.solve_energy!(atmos, surf_state=0, plot=true)

# Save result
atmosphere.write_pt(atmos, joinpath(output_dir,"pt.csv"))
plotting.plot_pt(atmos, joinpath(output_dir,"pt.pdf"))
plotting.plot_fluxes(atmos, joinpath(output_dir,"fluxes.pdf"))

# Deallocate
atmosphere.deallocate!(atmos)

println("Goodbye")
