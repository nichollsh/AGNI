#!/usr/bin/env -S julia --color=no --startup-file=no

# -------------
# AGNI executable file with command line arguments
# -------------

# Get AGNI root directory
ROOT_DIR = dirname(abspath(@__FILE__))
ENV["GKSwstype"] = "100"

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
    "inst"
        help = "Solar flux at the planet's orbital separation [W m-2]."
        arg_type = Float64
        required = true
    "zenith_degrees"
        help = "Direction of solar radiation beam measured from zenith [deg]."
        arg_type = Float64
        required = true
    "inst_factor"
        help = "Scale factor applied to instellation alongside the zenith angle"
        arg_type = Float64
        required = true
    "albedo_b"
        help = "Grey bond albedo applied to instellation"
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
        help = "Total pressure at the surface [bar]."
        arg_type = Float64
        required = true
    "ptop"
        help = "Total pressure at top of atmosphere  [bar]."
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
        help = "Do not allow the PT profile provided by `pt_path`` to overwrite the value of `tstar`"
        action = :store_true
    "--tmp_floor"
        help = "Minimum temperature allowed in the model - prevents numerical issues [K]."
        arg_type = Float64
        default = 0.5
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
        help = "Surface state (0: free, 1: constant, 2: conductive skin)."
        arg_type = Int
        default = 0
    "--skin_d"
        help = "Conductive skin thickness [m]."
        arg_type = Float64
        default = 0.01
    "--skin_k"
        help = "Conductive skin thermal conductivity [W m-1 K-1]."
        arg_type = Float64
        default = 2.0
    "--tmp_magma"
        help = "Magma temperature just below the surface skin [K]."
        arg_type = Float64
        default = 2000.0
    "--albedo_s"
        help = "Grey surface albedo."
        arg_type = Float64
        default = 0.0
    "--output"
        help = "Output directory relative to AGNI directory. This directory will be emptied before being used."
        arg_type = String
        default = "out"
    "--sp_file"
        help = "Spectral file path. Default is to Mallard."
        arg_type = String
        default = "res/spectral_files/Mallard/Mallard"
    "--star"
        help = "Path to stellar spectrum txt file. If not provided, spectral file is assumed to already include it."
        arg_type = String
        default = ""
    "--nlevels"
        help = "Number of model levels."
        arg_type = Int
        default = 80
    "--dtsolve"
        help = "Use time-stepped solver to obtain a solution."
        action = :store_true
    "--nlsolve"
        help = "Use non-linear solver to obtain solution once time-stepped method has finished."
        action = :store_true
    "--nsteps"
        help = "Maximum number of time-stepped solver steps."
        arg_type = Int
        default = 2000
    "--convect_adj"
        help = "Enable convection via dry adjustment."
        action = :store_true
    "--convect_mlt"
        help = "Enable convection via mixing length theory."
        action = :store_true
    "--noaccel"
        help = "Disable model acceleration."
        action = :store_true
    "--linesearch"
        help = "Enable linesearch with nonlinear solver."
        action = :store_true
    "--roverlap"
        help = "Use random overlap for computing overlapping absorption. Otherwise, equivalent extinction will be used."
        action = :store_true
    "--rscatter"
        help = "Include rayleigh scattering."
        action = :store_true
    "--cloud"
        help = "Include cloud optical properties in the radiative transfer."
        action = :store_true
    "--convcrit_tmprel"
        help = "Convergence criterion on dtmp/tmp/dt [day-1]."
        arg_type = Float64
        default = 4.0
    "--convcrit_fradrel"
        help = "Convergence criterion on dfrad/frad [%]."
        arg_type = Float64
        default = 0.8
    "--convcrit_flosspct"
        help = "Convergence criterion on global flux loss [%]."
        arg_type = Float64
        default = 2.0
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
instellation    = args["inst"]
s0_fact         = args["inst_factor"]
albedo_b        = args["albedo_b"]
zenith_degrees  = args["zenith_degrees"]
gravity         = args["gravity"]
radius          = args["radius"]
p_surf          = args["psurf"]
p_top           = args["ptop"]
nlev_centre     = args["nlevels"]
spfile_name     = args["sp_file"]
output_dir      = args["output"]
x_dict          = args["x_dict"]
x_path          = args["x_path"]
plot            = args["plot"]
pt_path         = args["pt_path"]
tstar_enforce   = args["tstar_enforce"]
tmp_floor       = args["tmp_floor"]
albedo_s        = args["albedo_s"]
ini_dry         = args["ini_dry"]
ini_sat         = args["ini_sat"]
trppt           = args["trppt"]
star_file       = args["star"]
rscatter        = args["rscatter"]
cloud           = args["cloud"]
verbose         = args["verbose"]
animate         = args["animate"]
surf_state      = args["surface"]
skin_d          = args["skin_d"]
skin_k          = args["skin_k"]
tmp_magma       = args["tmp_magma"]
max_steps       = args["nsteps"]
convect_adj     = args["convect_adj"]
convect_mlt     = args["convect_mlt"]
no_accel        = args["noaccel"]
linesearch      = args["linesearch"]
roverlap        = args["roverlap"]
dtsolve         = args["dtsolve"]
nlsolve         = args["nlsolve"]
cc_tmprel       = args["convcrit_tmprel"]
cc_fradrel      = args["convcrit_fradrel"]
cc_floss        = args["convcrit_flosspct"]

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
        if gas in keys(mf_dict)
            error("Mole fraction for '$gas' has been provided twice")
        end 
        mf_dict[gas] = val
    end 
end 
#    File path case
if x_path != ""
    mf_path = x_path
end

# Overlap method
if roverlap
    overlap = 2
else 
    overlap = 4
end

# Convection scheme
dry_convect = false
if convect_adj || convect_mlt
    dry_convect = true 

    if convect_adj && convect_mlt
        error("Both dry convection schemes are enabled! Pick only one at a time")
    end 

    if convect_adj && nlsol 
        error("Non-linear solver isn't compatible with convective adjustment")
    end
end

# Setup atmosphere
println("Setting up")
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         instellation, s0_fact, albedo_b, zenith_degrees,
                         tstar,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict,
                         mf_path=mf_path,
                         albedo_s=albedo_s,
                         skin_d=skin_d, skin_k=skin_k, tmp_magma=tmp_magma,
                         tmp_floor=tmp_floor,
                         overlap_method=overlap,
                         flag_gcontinuum=true,
                         flag_rayleigh=rscatter,
                         flag_cloud=cloud
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
if tstar_enforce
    atmos.tstar = tstar 
    atmos.tmpl[end] = tstar
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

# Oneshot case
if !(dtsolve || nlsolve)
    atmosphere.radtrans!(atmos, true)
    atmosphere.radtrans!(atmos, false)
    if convect_mlt
        atmosphere.mlt!(atmos)
    end
else 
    # Otherwise we will be solving with some method,
    # so check that the surface state is reasonable
    if (surf_state > 2) || (surf_state < 0)
        error("Invalid surface state '$surf_state'")
    end
end    

# Time-stepped solution
if dtsolve
    if animate 
        modplot = 10
    else
        modplot = 0
    end
    
    import solver_tstep
    solver_tstep.solve_energy!(atmos, 
                         modplot=modplot, verbose=verbose, 
                         surf_state=surf_state, dry_convect=dry_convect, use_mlt=convect_mlt,
                         max_steps=max_steps, accel=!no_accel, dt_max=150.0, rtol=1.0e-4, atol=1.0e-2,
                         drel_dt_conv=cc_tmprel, drel_F_conv=cc_fradrel, F_losspct_conv=cc_floss
                         )
end 
    
# Newton-Rapshon solution
if nlsolve
    import solver_nlsol
    solver_nlsol.solve_energy!(atmos, surf_state=surf_state, dry_convect=dry_convect, 
                                max_steps=200, atol=1e-2, use_linesearch=linesearch,
                                calc_cf_end=plot)
end

# Write NetCDF and PT files
atmosphere.write_ncdf(atmos,    joinpath(atmos.OUT_DIR,"atm.nc"))
atmosphere.write_pt(atmos,      joinpath(atmos.OUT_DIR,"pt.csv"))

# Final plots 
if animate
    println("Making animation")
    plotting.anim_solver(atmos)
end 
if plot 
    println("Making plots")
    plotting.plot_pt(atmos,     joinpath(atmos.OUT_DIR,"pt.pdf"))
    plotting.plot_fluxes(atmos, joinpath(atmos.OUT_DIR,"fl.pdf"))
end 

# Deallocate atmosphere
atmosphere.deallocate!(atmos)
