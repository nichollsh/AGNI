#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# AGNI executable file with command line arguments
# -------------

# Get AGNI root directory
ROOT_DIR = dirname(abspath(@__FILE__))

# Include libraries
using Revise
using ArgParse

trppt_default = 0.1

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
    "radius"
        help = "Planet radius at the surface [m]."
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
    "--load_csv"
        help = "Path to a CSV file containing a T(p) profile to load. Columns of file should be [Pa, K]"
        arg_type = String 
        default = ""
    "--ini_dry"
        help = "Initialise on a dry adiabat, with an isothermal stratosphere at this temperature."
        action = :store_true
    "--ini_sat"
        help = "Check each level for saturation and set to the saturation coexistance curve when T < T_dew."
        action = :store_true
    "--trppt"
        help = "Apply an isothermal stratosphere at this temperature."
        arg_type = Float64
        default = trppt_default
    "--surface"
        help = "Surface state (0: free, 1: fixed at T(p_surf), 2: fixed at T_star)"
        arg_type = Int
        default = 0
    "--albedo_s"
        help = "Grey surface albedo."
        arg_type = Float64
        default = 0.0
    "--zenith_degrees"
        help = "Direction of solar radiation beam measured from zenith [deg]."
        arg_type = Float64
        default = 54.74
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
        help = "Spectral file path. Default is Mallard."
        arg_type = String
        default = "res/spectral_files/Mallard/Mallard"
    "--stf"
        help = "Path to stellar spectrum txt file. If not provided, spectral file is assumed to already include it."
        arg_type = String
        default = ""
    "--nlevels"
        help = "Number of model levels."
        arg_type = Int
        default = 100
    "--nsteps"
        help = "Number of solver steps (max)."
        arg_type = Int
        default = 250
    "--noaccel"
        help = "Disable model acceleration."
        action = :store_true
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
radius          = args["radius"]
p_surf          = args["psurf"]
p_top           = args["ptop"]
nlev_centre     = args["nlevels"]
spfile_name     = args["spf"]
output_dir      = args["output"]
oneshot         = args["once"]
plot            = args["plot"]
csv_path        = args["load_csv"]
albedo_s        = args["albedo_s"]
zenith_degrees  = args["zenith_degrees"]
ini_dry         = args["ini_dry"]
ini_sat         = args["ini_sat"]
trppt           = args["trppt"]
star_file       = args["stf"]
rscatter        = args["rscatter"]
verbose         = args["verbose"]
animate         = args["animate"]
surf_state      = args["surface"]
max_steps       = args["nsteps"]
no_adjust       = args["noadjust"]
no_accel        = args["noaccel"]

if verbose 
    println("Command line arguments:")
    for (arg,val) in args
        println("    '$arg': $val")
    end
    println("")
end

spfile_has_star = (star_file == "")

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
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mixing_ratios,
                         zenith_degrees=zenith_degrees,
                         albedo_s=albedo_s,
                         flag_gcontinuum=true,
                         flag_rayleigh=rscatter
                         )
atmosphere.allocate!(atmos; 
                        spfile_has_star=spfile_has_star, 
                        stellar_spectrum=star_file
                        )

# Set PT profile 
#    Load CSV if required
if !(csv_path == "")
    setpt.fromcsv!(atmos, abspath(csv_path))
end 
#    Prevent surface supersaturation
setpt.prevent_surfsupersat!(atmos)
#    Apply dry adiabat if required
if ini_dry
    setpt.dry_adiabat!(atmos)
end 
#    Apply stratosphere if required
if trppt > trppt_default+1e-9
    setpt.stratosphere!(atmos, trppt)
end
#    Apply saturation curves if required 
if ini_sat
    for gas in atmos.gases 
        setpt.condensing!(atmos, gas)
    end 
end 

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
    if (surf_state > 2) || (surf_state < 0)
        error("Invalid surface state '$surf_state'")
    end
    import solver
    solver.solve_energy!(atmos, 
                         modplot=modplot, verbose=verbose, 
                         surf_state=surf_state, dry_adjust=!no_adjust, 
                         max_steps=max_steps, gofast=!no_accel, extrap=!no_accel
                         )
end

# Write NetCDF file 
atmosphere.write_ncdf(atmos,    joinpath(atmos.OUT_DIR,"atm.nc"))

# Write final PT profile and final fluxes 
atmosphere.write_pt(atmos,      joinpath(atmos.OUT_DIR,"pt.csv"))
atmosphere.write_fluxes(atmos,  joinpath(atmos.OUT_DIR,"fl.csv"))

# Final plots 
if animate && !oneshot
    println("Atmosphere: Making animation")
    plotting.anim_solver(atmos)
end 
if plot 
    println("Atmosphere: Making plots")
    plotting.plot_pt(atmos,     joinpath(atmos.OUT_DIR,"pt.pdf"))
    plotting.plot_fluxes(atmos, joinpath(atmos.OUT_DIR,"fl.pdf"))
end 

# Deallocate atmosphere
atmosphere.deallocate!(atmos)
