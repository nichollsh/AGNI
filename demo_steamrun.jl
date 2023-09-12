#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# AGNI executable file for demonstrating pure-steam runaway greenhouse effect
# -------------

# Get AGNI root directory
ROOT_DIR = dirname(abspath(@__FILE__))

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
import plotting

# Configuration options
toa_heating     = 1000.0  # SW dflux top boundary condition [W m-2]
gravity         = 9.81
radius          = 6.0e6
nlev_centre     = 100
p_surf          = 300.0     # bar
p_top           = 1e-8      # bar 
mixing_ratios   = Dict([("H2O", 1.0)])

spectral_file = "Oak"
star_file     = "res/stellar_spectra/sun.txt"
output_dir = "out/"

# Create output direct
rm(output_dir,force=true,recursive=true)
mkdir(output_dir)

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                        spectral_file,
                        toa_heating, 1700.0,
                        gravity, radius,
                        nlev_centre, p_surf, p_top,
                        mixing_ratios,
                        flag_gcontinuum=true
                        )
atmosphere.allocate!(atmos, stellar_spectrum=star_file)

p_boa = atmos.p_boa

plot_frames = false

run_len = 100
tsurf_arr = range(200,stop=2250,length=run_len)
olr_arr = zeros(Float64, run_len)
for i in 1:run_len
    
    # Set PT profile 
    atmos.p_boa = p_boa
    atmosphere.generate_pgrid!(atmos)
    atmos.tstar =       tsurf_arr[i]
    atmos.tmpl[end] =   tsurf_arr[i]

    setpt.prevent_surfsupersat!(atmos)
    setpt.dry_adiabat!(atmos)
    setpt.condensing!(atmos, "H2O")

    # Calculate LW fluxes 
    atmosphere.radtrans!(atmos, true)

    olr = atmos.flux_u_lw[1]
    @printf("Tsurf = %4.1f K  ,  OLR = %5.1f W m-2 \n",atmos.tmpl[end],olr)
    olr_arr[i] = olr 

    if plot_frames 
        tsurf = round(Int,atmos.tmpl[end])
        plotting.plot_pt(atmos, joinpath(output_dir,"pt_$tsurf.png"))
        atmosphere.write_pt(atmos, joinpath(output_dir, "pt_$tsurf.csv") )
    end

end

# Make plot of OLR vs T_surf
plt = plot(framestyle=:box)

lw=2.5

litdata = readdlm("res/runaway_litdata/Goldblatt13_data.txt", ',', Float64; header=false, skipstart=2)
plot!(plt, litdata[:,1] , litdata[:,2], label="Goldblatt+13", lw=lw)

litdata = readdlm("res/runaway_litdata/Hamano15_data.txt", ',', Float64; header=false, skipstart=2)
plot!(plt, litdata[:,1] , litdata[:,2], label="Hamano+15", lw=lw)

litdata = readdlm("res/runaway_litdata/Kopparapu13_data.txt", ',', Float64; header=false, skipstart=2)
plot!(plt, litdata[:,1] , litdata[:,2], label="Kopparapu+13", lw=lw)

plot!(tsurf_arr,olr_arr,label="AGNI", lw=lw)

xlabel!(plt, "Surface temperature [K]")
ylabel!(plt, "OLR [W m-2]")

savefig(plt, "out/runaway_olr.pdf")

# Deallocate atmosphere
atmosphere.deallocate!(atmos)

println("Goodbye")
