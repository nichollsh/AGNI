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
        help = "Surface temperature [K]."
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
    "--x_dict"
        help = "Mole fractions of volatiles formatted as a dictionary e.g. \"H2O=0.8,H2=0.2\" "
        arg_type = String 
        default = ""
    "--x_path"
        help = "Mole fractions of volatiles formatted as a CSV file (path to file)"
        arg_type = String 
        default = ""
    "--pt_path"
        help = "Path to a CSV file containing a T(p) profile to load. Columns of file should be [Pa, K]"
        arg_type = String 
        default = ""
    "--tstar_enforce"
        help = "Do not allow the PT profile provided by `pt_path`` to overwrite the value of `t_star`"
        action = :store_true
    "--tmp_floor"
        help = "Minimum temperature allowed in the model - prevents numerical issues [K]."
        arg_type = Float64
        default = 10.0
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
        help = "Surface state (0: free, 1: fixed at T(p_surf), 2: conductive skin)."
        arg_type = Int
        default = 0
    "--skin_d"
        help = "Conductive skin thickness [m]."
        arg_type = Float64
        default = 0.1
    "--skin_k"
        help = "Conductive skin thermal conductivity [W m-1 K-1]."
        arg_type = Float64
        default = 2.0
    "--tmp_magma"
        help = "Magma temperature just below the surface [K]."
        arg_type = Float64
        default = 2000.0
    "--albedo_s"
        help = "Grey surface albedo."
        arg_type = Float64
        default = 0.0
    "--zenith_degrees"
        help = "Direction of solar radiation beam measured from zenith [deg]."
        arg_type = Float64
        default = 54.74
    "--output"
        help = "Output directory relative to AGNI directory. This directory will be emptied before being used."
        arg_type = String
        default = "out"
    "--sp_file"
        help = "Spectral file path. Default is Mallard."
        arg_type = String
        default = "res/spectral_files/Mallard/Mallard"
    "--star"
        help = "Path to stellar spectrum txt file. If not provided, spectral file is assumed to already include it."
        arg_type = String
        default = ""
    "--nlevels"
        help = "Number of model levels."
        arg_type = Int
        default = 100
    "--once"
        help = "Just calculate fluxes once - don't iterate."
        action = :store_true
    "--nsteps"
        help = "Number of solver steps (max)."
        arg_type = Int
        default = 250
    "--noadjust"
        help = "Disable convective adjustment."
        action = :store_true
    "--noaccel"
        help = "Disable model acceleration."
        action = :store_true
    "--rscatter"
        help = "Include rayleigh scattering."
        action = :store_true
    "--convcrit_tmpabs"
        help = "Convergence criterion on dtmp [K]."
        arg_type = Float64
        default = 5.0
    "--convcrit_tmprel"
        help = "Convergence criterion on dtmp/tmp/dt [day-1]."
        arg_type = Float64
        default = 5.0
    "--convcrit_fradrel"
        help = "Convergence criterion on dfrad/frad [%]."
        arg_type = Float64
        default = 0.8
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
spfile_name     = args["sp_file"]
output_dir      = args["output"]
x_dict          = args["x_dict"]
x_path          = args["x_path"]
oneshot         = args["once"]
plot            = args["plot"]
pt_path         = args["pt_path"]
tstar_enforce   = args["tstar_enforce"]
tmp_floor       = args["tmp_floor"]
albedo_s        = args["albedo_s"]
zenith_degrees  = args["zenith_degrees"]
ini_dry         = args["ini_dry"]
ini_sat         = args["ini_sat"]
trppt           = args["trppt"]
star_file       = args["star"]
rscatter        = args["rscatter"]
verbose         = args["verbose"]
animate         = args["animate"]
surf_state      = args["surface"]
skin_d          = args["skin_d"]
skin_k          = args["skin_k"]
tmp_magma       = args["tmp_magma"]
max_steps       = args["nsteps"]
no_adjust       = args["noadjust"]
no_accel        = args["noaccel"]
cc_tmpabs       = args["convcrit_tmpabs"]
cc_tmprel       = args["convcrit_tmprel"]
cc_fradrel      = args["convcrit_fradrel"]

if verbose 
    println("Command line arguments:")
    for (arg,val) in args
        println("    '$arg': $val")
    end
    println("")
end

# Create output directory
if !isdir(output_dir) && !isfile(output_dir)
    mkdir(output_dir)
end

spfile_has_star = (star_file == "")

# Handle the provided mole fraction input 
#    Default case
mf_dict = nothing
mf_path = nothing
#    Dictionary case
if x_dict != ""
    mf_dict = Dict()

    x_strip = strip(args["x_dict"], [' ', '"'])
    x_split = split(x_strip, ",")

    if length(x_split) < 1
        error("No input mole fractions provided. Check the required dictionary formatting.")
    end 

    for pair in x_split # for each gas:val pair 
        gas_split = split(pair, "=")
        if length(gas_split) != 2
            error("Cannot parse input mole fractions. Check the required formatting.")
        end 
        gas = String(gas_split[1])
        val = parse(Float64,gas_split[2])
        if (val > 1.0) || (val < 0.0)
            error("Mole fractions must be between 0 and 1")
        end 
        if gas in keys(mole_fractions)
            error("Mole fraction for '$gas' has been provided twice")
        end 
        mole_fractions[gas] = val
    end 
end 
#    File path case
if x_path != ""
    mf_path = x_path
end

# Setup atmosphere
println("Atmosphere: setting up")
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, tstar,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict,
                         mf_path=mf_path,
                         zenith_degrees=zenith_degrees,
                         albedo_s=albedo_s,
                         skin_d=skin_d, skin_k=skin_k, tmp_magma=tmp_magma,
                         tmp_floor=tmp_floor,
                         flag_gcontinuum=true,
                         flag_rayleigh=rscatter
                         )
atmosphere.allocate!(atmos; 
                        spfile_has_star=spfile_has_star, 
                        stellar_spectrum=star_file
                        )

# Set PT profile 
#    Load CSV if required
if !(pt_path == "")
    setpt.fromcsv!(atmos, abspath(pt_path))
end 
#    Do not allow overwriting of tstar
if tstar_enforce:
    atmos.tstar = tstar 
    atmos.tmpl[end] = star
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
                         max_steps=max_steps, accel=!no_accel, extrap=!no_accel,
                         dtmp_conv=cc_tmpabs,drel_dt_conv=cc_tmprel, drel_F_conv=cc_fradrel
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
