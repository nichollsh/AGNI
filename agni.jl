#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# AGNI executable file for standalone execution
# -------------

# Get AGNI root directory
ROOT_DIR = dirname(abspath(@__FILE__))
ENV["GKSwstype"]="nul"

# Include libraries
using Revise
using LoggingExtras
using Printf
using TOML

# Include local jl files
include("socrates/julia/src/SOCRATES.jl")
push!(LOAD_PATH, joinpath(ROOT_DIR,"src"))
import atmosphere
import setpt
import plotting 
import phys
import solver_tstep
import solver_nlsol 

# Setup terminal + file logging 
function setup_logging(outpath::String)
    # Remove old file 
    rm(outpath, force=true)

    # Setup global logger
     logger_file = FormatLogger(outpath; append=true) do io, args
        @printf(io, "[%-5s] %s",args.level, args.message)
    end;

    # Setup terminal logger 
    logger_term = FormatLogger() do io, args
        color::Int = 39
        if args.level == LoggingExtras.Info
            color = 32
        elseif args.level == LoggingExtras.Warn
            color = 33
        elseif args.level == LoggingExtras.Debug
            color = 36
        end 

        # Set color, set bold, print level, unset bold, unset color, message
        @printf(io, "[\033[%dm\033[1m %-5s \033[21m\033[0m] %s",color, args.level, args.message)
    end;

    # Combine and set 
    logger_both = TeeLogger(logger_file, logger_term);
    global_logger(logger_both)
    disable_logging(Logging.Debug) # disable debug; info only

    return nothing 
end 

# Open and validate config file 
function open_config(cfg_path::String)::Dict

    # open file 
    cfg_dict = TOML.parsefile(cfg_path)

    # check headers 
    headers = ["plots", "planet", "files", "execution", "title"]
    for h in headers 
        if !haskey(cfg_dict, h)
            error("Key $h is missing from configuration file at '$cfg_path'")
        end 
    end 

    # check that output dir is named  
    if !haskey(cfg_dict["files"], "output_dir") || (cfg_dict["files"]["output_dir"] == "")
        error("Output directory is missing from configuration file at '$cfg_path'")
    end 
    out_path = abspath(cfg_dict["files"]["output_dir"])
    
    # check if this is a dangerous path
    if ispath(joinpath(out_path, ".git")) || (joinpath(out_path) == pwd())
        error("Output directory is unsafe")
    end 

    # looks good
    return cfg_dict 
end 

# Main function!
function main()

    # Record start time 
    tbegin = time()

    # Open and validate config file 
    cfg_path::String = joinpath(ROOT_DIR, "res/config/default.toml")
    if length(ARGS)>0
        cfg_path = ARGS[1]
    end
    if !ispath(cfg_path)
        error("Cannot find configuration file at '$cfg_path'")
    end 
    cfg = open_config(cfg_path)

    # Output folder 
    output_dir = abspath(cfg["files"]["output_dir"])
    rm(output_dir,force=true,recursive=true)
    mkdir(output_dir)

    # Copy configuration file 
    cp(cfg_path, joinpath(output_dir, "agni.cfg"))

    # Logging 
    setup_logging(joinpath(output_dir, "agni.log"))

    # Hello
    @info "Hello\n"
    @info "Using configuration '$(cfg["title"])'\n"

    # Read REQUIRED configuration options from dict 
    #    planet stuff 
    t_star::Float64        = cfg["planet"]["t_star"]
    instellation::Float64  = cfg["planet"]["instellation"]
    albedo_b::Float64      = cfg["planet"]["albedo_b"]
    albedo_s::Float64      = cfg["planet"]["albedo_s"]
    asf_sf::Float64        = cfg["planet"]["s0_fact"]
    radius::Float64        = cfg["planet"]["radius"]
    zenith::Float64        = cfg["planet"]["zenith_angle"]
    gravity::Float64       = cfg["planet"]["gravity"]
    p_surf::Float64        = cfg["planet"]["p_surf"]
    p_top::Float64         = cfg["planet"]["p_top"]
    #    solver stuff 
    spfile_name::String    = cfg["files" ]["input_sf"]
    star_file::String      = cfg["files" ]["input_star"]
    nlev_centre::Int       = cfg["execution"]["num_levels"]
    flag_cnt::Bool         = cfg["execution" ]["continua"]
    flag_ray::Bool         = cfg["execution" ]["rayleigh"]
    flag_cld::Bool         = cfg["execution" ]["cloud"]
    flag_aer::Bool         = cfg["execution" ]["aerosol"]
    overlap::Int           = cfg["execution" ]["overlap_method"]
    thermo_funcs::Bool     = cfg["execution" ]["thermo_funcs"]
    dry_type::String       = cfg["execution" ]["dry_convection"]
    surf_state::Int        = cfg["execution" ]["surf_state"]
    solvers_cmd::Array     = cfg["execution" ]["solvers"]
    initial_cmd::Array     = cfg["execution" ]["initial_state"]
    stabilise::Bool        = cfg["execution" ]["stabilise"]
    conv_atol::Float64     = cfg["execution" ]["converge_atol"]
    max_steps::Int         = cfg["execution" ]["max_steps"]
    #    plotting stuff 
    plt_run::Bool          = cfg["plots"     ]["at_runtime"]
    plt_tmp::Bool          = cfg["plots"     ]["temperature"]
    plt_flx::Bool          = cfg["plots"     ]["fluxes"]
    plt_cff::Bool          = cfg["plots"     ]["contribution"]
    plt_ems::Bool          = cfg["plots"     ]["emission"]
    plt_alb::Bool          = cfg["plots"     ]["albedo"]
    plt_vmr::Bool          = cfg["plots"     ]["mixing_ratios"]
    plt_ani::Bool          = cfg["plots"     ]["animate"]

    # Read OPTIONAL configuration options from dict
    #     mixing ratios can be set either way
    if haskey(cfg["planet"],"vmr")
        mf_dict::Dict = cfg["planet"]["vmr"]
        mf_path = nothing
    else 
        mf_dict = nothing
        mf_path = cfg["files"]["input_vmr"]
    end 
    #     sensible heat at the surface 
    turb_coeff::Float64 = 0.0; wind_speed::Float64 = 0.0
    if cfg["execution"]["sensible_heat"]
        turb_coeff = cfg["planet"]["turb_coeff"]
        wind_speed = cfg["planet"]["wind_speed"]
    end 
    #     conductive skin case 
    skin_k::Float64=0.0; skin_d::Float64=0.0; t_magma::Float64=0.0
    if surf_state == 2
        skin_k  = cfg["planet"]["skin_k"]
        skin_d  = cfg["planet"]["skin_d"]
        t_magma = cfg["planet"]["t_magma"]
    end 
    #     interior temperature case 
    t_int::Float64 = 0.0
    if surf_state == 3
        t_int = cfg["planet"]["t_int"]
    end 


    # Setup atmosphere
    @info "Setting up\n"
    atmos = atmosphere.Atmos_t()
    atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                            spfile_name,
                            instellation, asf_sf, albedo_b, zenith,
                            t_star, 
                            gravity, radius,
                            nlev_centre, p_surf, p_top,
                            mf_dict=mf_dict, mf_path=mf_path,
                            flag_gcontinuum=flag_cnt, flag_rayleigh=flag_ray,
                            flag_cloud=flag_cld, flag_aerosol=flag_aer,
                            overlap_method=overlap,
                            skin_d=skin_d, skin_k=skin_k, tmp_magma=t_magma,
                            tmp_floor=5.0,
                            tint=t_int, albedo_s=albedo_s,
                            thermo_functions=thermo_funcs,
                            C_d=turb_coeff, U=wind_speed
                    )
    atmosphere.allocate!(atmos;stellar_spectrum=star_file,spfile_noremove=true,spfile_has_star=isempty(star_file))

    # Set PT profile 
    @info "Setting initial T(p)\n"
    # setpt.fromcsv!(atmos,"pt.csv")
    # setpt.isothermal!(atmos, t_star*0.7)
    # setpt.prevent_surfsupersat!(atmos)
    setpt.dry_adiabat!(atmos)
    setpt.condensing!(atmos, "H2O")
    # setpt.stratosphere!(atmos, 800.0)

    atmosphere.write_pt(atmos, joinpath(atmos.OUT_DIR,"pt_ini.csv"))

    # Solver variables 
    dry_convect::Bool = !isempty(dry_type)
    use_mlt::Bool     = (dry_type == "mlt")
    modplot::Int      = 0
    condensate  = ""

    if plt_run 
        modplot = 1
    end

    # Loop over requested solvers 
    method_map::Array{String,1} = ["newton", "gauss", "levenberg"]
    method::Int = 0
    for sol in solvers_cmd 

        sol = lowercase(sol)

        # No solve - just calc fluxes at the end
        if sol == ""
            @info "Solver = none\n"
            atmosphere.radtrans!(atmos, true, calc_cf=true)
            atmosphere.radtrans!(atmos, false)
            if use_mlt 
                atmosphere.mlt!(atmos)
            end 
        
        # Timestepping
        elseif sol == "timestep"
            @info "Solver = $sol\n"
            solver_tstep.solve_energy!(atmos, surf_state=surf_state, use_physical_dt=false,
                                modplot=modplot, modprop=5, verbose=true, 
                                dry_convect=dry_convect, condensate=condensate,
                                accel=stabilise, step_rtol=1.0e-4, step_atol=1.0e-2, dt_max=1000.0,
                                max_steps=max_steps, min_steps=100, use_mlt=use_mlt)
        
        # Nonlinear methods
        elseif (sol in method_map) 
            @info "Solver = $sol\n"
            method = findfirst(==(sol), method_map)
            solver_nlsol.solve_energy!(atmos, surf_state=surf_state, 
                                dry_convect=dry_convect, condensate=condensate,
                                max_steps=max_steps, conv_atol=conv_atol, method=1,
                                stabilise_mlt=stabilise,modplot=modplot)
        else 
            error("Invalid solver requested '$sol'")
        end 

    end 

    @info "Total RT evalulations: $(atmos.num_rt_eval)\n"

    # Write arrays
    atmosphere.write_pt(atmos,      joinpath(atmos.OUT_DIR,"pt.csv"))
    atmosphere.write_ncdf(atmos,    joinpath(atmos.OUT_DIR,"atm.nc"))
    atmosphere.write_fluxes(atmos,  joinpath(atmos.OUT_DIR,"fl.csv"))

    # Save plots
    @info "Making plots\n"
    plt_ani && plotting.anim_solver(atmos)
    plt_vmr && plotting.plot_x(atmos,          joinpath(atmos.OUT_DIR,"mf.png"))
    plt_cff && plotting.plot_contfunc(atmos,   joinpath(atmos.OUT_DIR,"cf.png"))
    plt_tmp && plotting.plot_pt(atmos,         joinpath(atmos.OUT_DIR,"pt.png"), incl_magma=(surf_state==2))
    plt_flx && plotting.plot_fluxes(atmos,     joinpath(atmos.OUT_DIR,"fl.png"))
    plt_ems && plotting.plot_emission(atmos,   joinpath(atmos.OUT_DIR,"em.png"))
    plt_alb && plotting.plot_albedo(atmos,     joinpath(atmos.OUT_DIR,"al.png"))

    # Deallocate atmosphere
    @info "Deallocating arrays\n"
    atmosphere.deallocate!(atmos)

    # Finish up
    runtime = round(time() - tbegin, digits=2)
    @info "Runtime: $runtime seconds\n"
    @info "Goodbye\n"

    return nothing 
end 

# Call main function 
main()
