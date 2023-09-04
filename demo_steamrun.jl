#!/usr/bin/env -S julia --color=yes --startup-file=no

println("Begin runaway demo")

using Revise
using Printf
using Plots
using DelimitedFiles

# Include local jl files
include("socrates/julia/src/SOCRATES.jl")

push!(LOAD_PATH, joinpath(pwd(),"src"))
import atmosphere
import setpt
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

spectral_file = "res/spectral_files/Oak/Oak"
star_file     = "res/stellar_spectra/spec_sun.txt"
output_dir = "out/"

# Create output direct
rm(output_dir,force=true,recursive=true)
mkdir(output_dir)


# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, spectral_file,
                         star_file, zenith_degrees, 
                         toa_heating, tstar,
                         gravity, nlev_centre, p_surf, p_top,
                         mixing_ratios,
                         flag_gcontinuum=true
                         )
atmosphere.allocate!(atmos)

run_len = 30
tsurf_arr = range(100,stop=2200,length=run_len)
olr_arr = zeros(Float64, run_len)
for i in 1:run_len

    # Set PT profile 
    atmos.tstar = tsurf_arr[i]
    setpt.dry_adiabat!(atmos)
    setpt.condensing!(atmos, "H2O")

    # Calculate LW and SW fluxes (once)
    atmosphere.radtrans!(atmos, true)
    atmosphere.radtrans!(atmos, false)

    olr = atmos.flux_u_lw[1]
    @printf("Tsurf = %4.1f K  ,  OLR = %5.1f W m-2 \n",tsurf_arr[i],olr)
    olr_arr[i] = olr 

end

# Deallocate atmosphere
atmosphere.deallocate!(atmos)

# Make plot
plt = plot()

litdata = readdlm("res/runaway_litdata/Goldblatt13_data.txt", ',', Float64; header=false, skipstart=2)
plot!(plt, litdata[:,1] , litdata[:,2], label="Goldblatt+13")

litdata = readdlm("res/runaway_litdata/Hamano15_data.txt", ',', Float64; header=false, skipstart=2)
plot!(plt, litdata[:,1] , litdata[:,2], label="Hamano+15")

litdata = readdlm("res/runaway_litdata/Kopparapu13_data.txt", ',', Float64; header=false, skipstart=2)
plot!(plt, litdata[:,1] , litdata[:,2], label="Kopparapu+13")

plot!(tsurf_arr,olr_arr,label="THEMIS")

xlabel!(plt, "Surface temperature [K]")
ylabel!(plt, "OLR [W m-2]")
title!(plt, "Steam runaway greenhouse")

savefig(plt, "out/runaway_olr.pdf")

println("Goodbye")
