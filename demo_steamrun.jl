#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# AGNI executable file for demonstrating pure-steam runaway greenhouse effect
# -------------

# Get AGNI root directory
ROOT_DIR = dirname(abspath(@__FILE__))
ENV["GKSwstype"] = "100"

using Printf
using Plots
using DelimitedFiles
using LoggingExtras

@info "Begin runaway demo"

# Include local jl files
include("src/phys.jl")
include("src/moving_average.jl")
include("src/spectrum.jl")
include("src/atmosphere.jl")
include("src/setpt.jl")
include("src/plotting.jl")
include("src/energy.jl")
import .phys
import .atmosphere
import .setpt
import .plotting
import .energy

# Configuration options
instellation    = 1000.0  # Solar flux [W m-2]
gravity         = 9.81
radius          = 6.0e6
nlev_centre     = 100
p_surf          = 300.0     # bar
p_top           = 1e-8      # bar 
mole_fractions  = Dict([("H2O", 1.0)])

spectral_file = "res/spectral_files/Oak/Oak.sf"
star_file     = "res/stellar_spectra/sun.txt"
output_dir = "out/"

# Create output direct
rm(output_dir,force=true,recursive=true)
mkdir(output_dir)

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                        spectral_file,
                        instellation, 1.0, 0.0, 48.19,
                        1700.0,
                        gravity, radius,
                        nlev_centre, p_surf, p_top,
                        mf_dict=mole_fractions,
                        flag_gcontinuum=true,
                        thermo_functions=false
                        )
atmosphere.allocate!(atmos, star_file)

p_boa = atmos.p_boa

plot_frames = false

run_len = 200
tsurf_arr = range(200,stop=2400,length=run_len)
olr_arr = zeros(Float64, run_len)
for i in 1:run_len
    
    # Set PT profile 
    atmos.p_boa = p_boa
    atmosphere.generate_pgrid!(atmos)
    atmos.tmp_surf =       tsurf_arr[i]
    atmos.tmpl[end] =   tsurf_arr[i]

    setpt.prevent_surfsupersat!(atmos)
    setpt.dry_adiabat!(atmos)
    setpt.condensing!(atmos, "H2O")

    # Calculate LW fluxes 
    energy.radtrans!(atmos, true)

    olr = atmos.flux_u_lw[1]
    @info @sprintf("Tsurf = %4.1f K  ,  OLR = %5.1f W m-2 ",atmos.tmpl[end],olr)
    olr_arr[i] = olr 

    if plot_frames 
        tsurf = round(Int,atmos.tmpl[end])
        plotting.plot_pt(atmos, joinpath(output_dir,"pt_$tsurf.png"))
        atmosphere.write_pt(atmos, joinpath(output_dir, "pt_$tsurf.csv") )
    end

end

# Make plot of OLR vs T_surf
@info "Making plot"
plt = plot(framestyle=:box, size=(500,400), dpi=300)

lw=2.5

litdata = readdlm("res/literature_data/runaway/Goldblatt13_data.txt", ',', Float64; header=false, skipstart=2)
plot!(plt, litdata[:,1] , litdata[:,2], label="Goldblatt+13", lw=lw)

litdata = readdlm("res/literature_data/runaway/Hamano15_data.txt", ',', Float64; header=false, skipstart=2)
plot!(plt, litdata[:,1] , litdata[:,2], label="Hamano+15", lw=lw)

litdata = readdlm("res/literature_data/runaway/Kopparapu13_data.txt", ',', Float64; header=false, skipstart=2)
plot!(plt, litdata[:,1] , litdata[:,2], label="Kopparapu+13", lw=lw)

plot!(tsurf_arr,olr_arr,label="AGNI", lw=lw)

xlabel!(plt, "Surface temperature [K]")
ylabel!(plt, "OLR [W m-2]")

savefig(plt, "out/runaway_olr.pdf")
savefig(plt, "out/runaway_olr.png")

# Write data 
@info "Writing data"
open("out/OLR_vs_Tsurf.csv","w") do hdl 
    writedlm(hdl, [tsurf_arr olr_arr], ", ")
end 

# Deallocate atmosphere
atmosphere.deallocate!(atmos)

@info "Goodbye"
exit(0)
