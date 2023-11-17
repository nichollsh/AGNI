#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# AGNI executable file for demonstrating the greenhouse effect for multicomponent atmospheres
# -------------

# Get AGNI root directory
ROOT_DIR = dirname(abspath(@__FILE__))
ENV["GKSwstype"] = "100"

println("Begin complex runaway demo")

using Revise
using Printf
using Plots
using DelimitedFiles

# Include local jl files
include("socrates/julia/src/SOCRATES.jl")

push!(LOAD_PATH, joinpath(pwd(),"src"))
import atmosphere
import setpt
import solver_euler
import phys
import plotting


# Configuration options
tstar           = 2235.0    # Surface temperature [kelvin]
toa_heating     = 3.772e+04 # Instellation flux [W m-2]
radius          = 7.12e6    # metres
gravity         = 10.8      # m s-2
nlev_centre     = 100  
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
                         mf_dict=mf_dict,
                         flag_gcontinuum=true,
                         flag_rayleigh=true,
                         overlap_method=4,
                         zenith_degrees=54.4,
                         skin_d=0.02,
                         skin_k=2.0,
                         tmp_magma=2700.
                 )
atmosphere.allocate!(atmos;stellar_spectrum=star_file,spfile_noremove=true)
p_boa = atmos.p_boa

plot_frames = true
run_len = 20
tsurf_arr = range(500,stop=3000,length=run_len)

toa_arr = zeros(Float64, run_len)
boa_arr = zeros(Float64, run_len)
skn_arr = zeros(Float64, run_len)

for i in 1:run_len

    # Set PT profile 
    atmos.p_boa = p_boa
    atmosphere.generate_pgrid!(atmos)
    atmos.tstar   =     tsurf_arr[i]

    atmos.tmpl[:] .=     atmos.tstar # make isothermal
    atmos.tmp[:]  .=     atmos.tstar

    @printf("Running Tsurf = %3.1f K \n",atmos.tmpl[end])

    # Calculate temperature profile
    solver_euler.solve_energy!(atmos, surf_state=1, modplot=0, verbose=false, dry_convect=true, max_steps=400, min_steps=20, use_mlt=true)

    F_toa = atmos.flux_tot[1]
    F_boa = atmos.flux_tot[end]
    F_skn = atmos.skin_k / atmos.skin_d * (atmos.tmp_magma - atmos.tstar)
    
    toa_arr[i] = F_toa
    boa_arr[i] = F_boa 
    skn_arr[i] = F_skn 

    @printf("Completed Tsurf = %3.1f K  ,  F_TOA = %.4e W m-2 ,  F_SKN = %.4e W m-2\n",atmos.tmpl[end],F_toa,F_skn)

    if plot_frames 
        tsurf = round(Int,atmos.tmpl[end])
        plotting.plot_pt(atmos,     joinpath(output_dir, "pt_$tsurf.png"), dpi=220)
        plotting.plot_fluxes(atmos, joinpath(output_dir, "fl_$tsurf.png"), dpi=220)
        atmosphere.write_pt(atmos,  joinpath(output_dir, "pt_$tsurf.csv"))
    end

end

# Make plot of Fluxes vs T_surf
plt = plot(framestyle=:box, size=(500,400))

lw=2.5

plot!(tsurf_arr,toa_arr,label="TOA",  lw=lw)
plot!(tsurf_arr,boa_arr,label="BOA",  lw=lw)
plot!(tsurf_arr,skn_arr,label="Skin", lw=lw)

xlabel!(plt, "Surface temperature [K]")
ylabel!(plt, "Net flux [W m-2]")

savefig(plt, "out/runaway_complex.pdf")

# Deallocate atmosphere
atmosphere.deallocate!(atmos)

println("Goodbye")
