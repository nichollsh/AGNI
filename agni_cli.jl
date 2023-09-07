#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# AGNI main file with command line arguments
# -------------

println("AGNI CLI")

# Get AGNI root directory
ROOT_DIR = dirname(abspath(@__FILE__))

# Include libraries
using Revise
using ArgParse

# Command line arguments
s = ArgParseSettings()
@add_arg_table s begin
    "tstar"
        help = "LW uflux bottom boundary condition [K]."
        arg_type = Float64
        required = true
    "toa_heating"
        help = "SW dflux top boundary condition [W m-2]."
        arg_type = Float64
        required = true
    "gravity"
        help = "Surface gravitational acceleration [m s-2]."
        arg_type = Float64
        required = true
    "psurf"
        help = "Total surface pressure [bar]."
        arg_type = Float64
        required = true
    "ptop"
        help = "Total top of atmosphere pressure [bar]."
        arg_type = Float64
        required = true
    "mr_dict"
        help = "Mixing ratios of volatiles formatted as a dictionary e.g. \"H2O=0.8,H2=0.2\" "
        arg_type = String 
        required = true
    "--trppt"
        help = "Initialise with an isothermal stratosphere at this temperature."
        arg_type = Float64
        default = 150.0
    "--cstsurf"
        help = "Fix the surface temperature to be constant."
        action = :store_true
    "--noadjust"
        help = "Disable convective adjustment."
        action = :store_true
    "--once"
        help = "Just calculate fluxes once - don't iterate."
        action = :store_true
    "--output"
        help = "Output directory relative to AGNI directory. This directory will be emptied before being used."
        arg_type = String
        default = "out"
    "--spf"
        help = "Spectral file (name). Default is Mallard."
        arg_type = String
        default = "Mallard"
    "--spf_keep"
        help = "Keep runtime spectral file."
        action = :store_true
    "--stf"
        help = "Path to stellar spectrum txt file. If not provided, spectral file is assumed to already include it."
        arg_type = String
        default = ""
    "--nlevels"
        help = "Number of model levels."
        arg_type = Int64
        default = 100
    "--rscatter"
        help = "Include rayleigh scattering."
        action = :store_true
    "--plot"
        help = "Make plots."
        action = :store_true
    "--animate"
        help = "Make animation."
        action = :store_true
    "--verbose"
        help = "Enable verbose output."
        action = :store_true
end
args = parse_args(s)

# Include local jl files
include("socrates/julia/src/SOCRATES.jl")
push!(LOAD_PATH, joinpath(ROOT_DIR,"src"))
import atmosphere
import setpt
import plotting 

# Set the configuration options
tstar           = args["tstar"]
toa_heating     = args["toa_heating"]
gravity         = args["gravity"]
p_surf          = args["psurf"]
p_top           = args["ptop"]
nlev_centre     = args["nlevels"]
spfile_name     = args["spf"]
output_dir      = args["output"]
oneshot         = args["once"]
plot            = args["plot"]
trppt           = args["trppt"]
star_file       = args["stf"]
rscatter        = args["rscatter"]
verbose         = args["verbose"]
animate         = args["animate"]
fixed_surf      = args["cstsurf"]
no_adjust       = args["noadjust"]

if verbose 
    println("Command line arguments:")
    for (arg,val) in args
        println("    '$arg': $val")
    end
    println("")
end

# Remove old files 
if abspath(output_dir) == abspath(ROOT_DIR)
    error("Output directory cannot be the root directory")
end
rm(output_dir,force=true,recursive=true)
mkdir(output_dir)

spfile_has_star = (star_file == "")
spfile_noremove = args["spf_keep"]

# Parse the provided mr dict 
mixing_ratios = Dict()

mr_strip = strip(args["mr_dict"], [' ', '"'])
mr_split = split(mr_strip, ",")

if length(mr_split) < 1
    error("No input mixing ratios provided. Check the required formatting.")
end 

for pair in mr_split # for each gas:val pair 
    gas_split = split(pair, "=")
    if length(gas_split) != 2
        error("Cannot parse input mixing ratios. Check the required formatting.")
    end 
    gas = String(gas_split[1])
    val = parse(Float64,gas_split[2])
    if (val > 1.0) || (val < 0.0)
        error("Mixing ratios must be between 0 and 1")
    end 
    if gas in keys(mixing_ratios)
        error("Mixing ratio for '$gas' has been provided twice")
    end 
    mixing_ratios[gas] = val
end 

# Setup atmosphere
println("Atmosphere: setting up")
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, tstar,
                         gravity, nlev_centre, p_surf, p_top,
                         mixing_ratios,
                         flag_gcontinuum=true,
                         flag_rayleigh=rscatter
                         )
atmosphere.allocate!(atmos; 
                        spfile_has_star=spfile_has_star, 
                        stellar_spectrum=star_file, 
                        spfile_noremove=spfile_noremove
                        )

# Set PT profile 
setpt.prevent_surfsupersat!(atmos)
setpt.dry_adiabat!(atmos)
setpt.stratosphere!(atmos, trppt)

# Just once or iterate?
if oneshot
    println("RadTrans: Calculating fluxes")
    atmosphere.radtrans!(atmos, true)
    atmosphere.radtrans!(atmos, false)
else
    if animate 
        modplot = 1 
    else
        modplot = 0
    end
    if fixed_surf
        surf_state = 1 
    else 
        surf_state = 0
    end
    import solver
    solver.solve_energy!(atmos, modplot=modplot, verbose=verbose, surf_state=surf_state, dry_adjust=!no_adjust)
end

# Write final PT profile
atmosphere.write_pt(atmos,  joinpath(atmos.OUT_DIR,"pt.csv"))

# Final plots 
if animate && !oneshot
    println("Atmosphere: Making animation")
    plotting.anim_solver(atmos)
end 
if plot 
    println("Atmosphere: Making plots")
    plotting.plot_pt(atmos,     joinpath(atmos.OUT_DIR,"pt.pdf"))
    plotting.plot_fluxes(atmos, joinpath(atmos.OUT_DIR,"fluxes.pdf"))
end 

# Deallocate atmosphere
atmosphere.deallocate!(atmos)

# Done
println("Goodbye")
println(" ")
