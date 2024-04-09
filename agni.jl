#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# AGNI executable file for standalone execution
# -------------

# Get AGNI root directory
ROOT_DIR = dirname(abspath(@__FILE__))
ENV["GKSwstype"] = "100"

# Include libraries
using Revise
using LoggingExtras
using Printf
using TOML

# Include local jl files
push!(LOAD_PATH, joinpath(ROOT_DIR,"src"))
import atmosphere
import setpt
import plotting 
import phys
import solver_tstep
import solver_nlsol 

# Setup terminal + file logging 
function setup_logging(outpath::String, silent::Bool)
    # Remove old file 
    rm(outpath, force=true)

    # If silent 
    if silent 
        global_logger(MinLevelLogger(current_logger(), Logging.Error))
        return nothing
    end 

    # Formatting
    color::Int = 39
    level::String = "UNSET"

    # Setup file logger
    logger_file = FormatLogger(outpath; append=true) do io, args
        if args.level == LoggingExtras.Info
            level = "INFO"
        elseif args.level == LoggingExtras.Warn
            level = "WARN"
        elseif args.level == LoggingExtras.Debug
            level = "DEBUG"
        elseif args.level == LoggingExtras.Error 
            level = "ERROR"
        end 
        @printf(io, "[ %-5s ] %s \n", level, args.message)
    end;

    # Setup terminal logger 
    logger_term = FormatLogger() do io, args
        if args.level == LoggingExtras.Info
            color = 32
            level = "INFO"
        elseif args.level == LoggingExtras.Warn
            color = 93
            level = "WARN"
        elseif args.level == LoggingExtras.Debug
            color = 96
            level = "DEBUG"
        elseif args.level == LoggingExtras.Error 
            color = 91
            level = "ERROR"
        end 
        # Set color, set bold, print level, unset bold, unset color, message
        @printf(io, "[\033[%dm\033[1m %-5s \033[21m\033[0m] %s \n",color, level, args.message)
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
    if cfg["files"]["clean_output"] || !(ispath(output_dir) && isdir(output_dir))
        rm(output_dir,force=true,recursive=true)
        mkdir(output_dir)
    end 

    # Copy configuration file 
    cp(cfg_path, joinpath(output_dir, "agni.cfg"))

    # Logging 
    silent::Bool = cfg["execution"]["silent"]
    setup_logging(joinpath(output_dir, "agni.log"), silent)

    # Hello
    @info "Hello"
    @info "Using configuration '$(cfg["title"])'"

    # Read REQUIRED configuration options from dict 
    #    planet stuff 
    tmp_surf::Float64      = cfg["planet"]["tmp_surf"]
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
    thermo_funct::Bool     = cfg["execution" ]["thermo_funct"]
    dry_type::String       = cfg["execution" ]["dry_convection"]
    condensates::Array     = cfg["execution" ]["condensates"]
    incl_sens::Bool        = cfg["execution" ]["sensible_heat"]
    sol_type::Int          = cfg["execution" ]["solution_type"]
    solvers_cmd::Array     = cfg["execution" ]["solvers"]
    initial_req::Array     = cfg["execution" ]["initial_state"]
    stabilise::Bool        = cfg["execution" ]["stabilise"]
    conv_atol::Float64     = cfg["execution" ]["converge_atol"]
    conv_rtol::Float64     = cfg["execution" ]["converge_rtol"]
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
        mf_dict = cfg["planet"]["vmr"]
        mf_path = nothing
    else 
        mf_dict = nothing
        mf_path = cfg["files"]["input_vmr"]
    end 
    #     sensible heat at the surface 
    turb_coeff::Float64 = 0.0; wind_speed::Float64 = 0.0
    if incl_sens
        turb_coeff = cfg["planet"]["turb_coeff"]
        wind_speed = cfg["planet"]["wind_speed"]
    end 
    #     conductive skin case 
    skin_k::Float64=0.0; skin_d::Float64=0.0; tmp_magma::Float64=0.0
    if sol_type == 2
        skin_k      = cfg["planet"]["skin_k"]
        skin_d      = cfg["planet"]["skin_d"]
        tmp_magma   = cfg["planet"]["tmp_magma"]
    end 
    #     effective temperature case 
    tmp_eff::Float64 = 0.0
    if sol_type == 3
        tmp_eff = cfg["planet"]["tmp_eff"]
    end 
    #     target OLR case 
    target_olr::Float64 = 0.0
    if sol_type == 4
        target_olr = cfg["planet"]["target_olr"]
    end    


    # Setup atmosphere
    @info "Setting up"
    atmos = atmosphere.Atmos_t()
    atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                            spfile_name,
                            instellation, asf_sf, albedo_b, zenith,
                            tmp_surf, 
                            gravity, radius,
                            nlev_centre, p_surf, p_top,
                            mf_dict=mf_dict, mf_path=mf_path,
                            flag_gcontinuum=flag_cnt, flag_rayleigh=flag_ray,
                            flag_cloud=flag_cld, flag_aerosol=flag_aer,
                            overlap_method=overlap,
                            skin_d=skin_d, skin_k=skin_k, tmp_magma=tmp_magma,
                            tmp_floor=5.0, target_olr=target_olr,
                            tmp_eff=tmp_eff, albedo_s=albedo_s,
                            thermo_functions=thermo_funct,
                            C_d=turb_coeff, U=wind_speed
                    )
    atmosphere.allocate!(atmos;stellar_spectrum=star_file,spfile_noremove=true,spfile_has_star=isempty(star_file))

    # Set PT profile by looping over requests
    # Each request may be a command, or an argument following a command
    @info "Setting initial T(p)"
    num_req::Int = length(initial_req)      # Number of requests
    idx_req::Int = 1                        # Index of current request
    str_req::String = "_unset"              # String of current request
    while idx_req <= num_req
        # get command 
        str_req = strip(lowercase(initial_req[idx_req]))

        # handle requests  
        if str_req == "dry"
            # dry adiabat from surface
            setpt.dry_adiabat!(atmos)

        elseif str_req == "str"
            # isothermal stratosphere 
            idx_req += 1
            setpt.stratosphere!(atmos, parse(Float64, initial_req[idx_req]))
            
        elseif str_req == "iso"
            # isothermal profile 
            idx_req += 1
            setpt.isothermal!(atmos, parse(Float64, initial_req[idx_req]))
        
        elseif str_req == "csv"
            # set from csv file 
            idx_req += 1
            setpt.fromcsv!(atmos,initial_req[idx_req])

        elseif str_req == "sat"
            # check surface supersaturation
            setpt.prevent_surfsupersat!(atmos)
        
        elseif str_req == "con"
            # condensing a volatile 
            setpt.condensing!(atmos, initial_req[idx_req])

        else 
            error("Invalid initial state '$str_req'")
        end 
        
        # iterate
        idx_req += 1
    end 

    # Write initial state
    atmosphere.write_pt(atmos, joinpath(atmos.OUT_DIR,"pt_ini.csv"))

    # Solver variables 
    dry_convect::Bool = !isempty(dry_type)
    use_mlt::Bool     = (dry_type == "mlt")
    modplot::Int      = 0

    # Loop over requested solvers 
    method_map::Array{String,1} = ["newton", "gauss", "levenberg"]
    method::Int = 0
    if length(solvers_cmd) == 0  # is empty 
        solvers_cmd = [""]
    end 
    if !isempty(solvers_cmd[end])  # append "no solve" case to end, for calculating cff
        push!(solvers_cmd, "")
    end 
    for sol in solvers_cmd 

        sol = strip(lowercase(sol))

        # No solve - just calc fluxes at the end
        if isempty(sol)
            @info "Solver = none"
            fill!(atmos.flux_tot, 0.0)
            atmosphere.radtrans!(atmos, true, calc_cf=true)
            atmosphere.radtrans!(atmos, false)
            if use_mlt 
                atmosphere.mlt!(atmos)
            end 
            if incl_sens 
                atmosphere.sensible!(atmos)
            end 
            atmosphere.conduct!(atmos)
            atmos.flux_tot = atmos.flux_cdry + atmos.flux_n + atmos.flux_cdct
            atmos.flux_tot[end] += atmos.flux_sens
        
        # Timestepping
        elseif sol == "timestep"
            @info "Solver = $sol"
            # Plotting at runtime
            if plt_run 
                modplot = 10
            end
            solver_tstep.solve_energy!(atmos, sol_type=sol_type, use_physical_dt=false,
                                modplot=modplot, modprop=5, verbose=true,  sens_heat=incl_sens,
                                dry_convect=dry_convect, condensates=condensates,
                                accel=stabilise, step_rtol=1.0e-4, step_atol=1.0e-2, dt_max=1000.0,
                                conv_atol=conv_atol, conv_rtol=conv_rtol,
                                max_steps=max_steps, min_steps=100, use_mlt=use_mlt)
        
        # Nonlinear methods
        elseif (sol in method_map) 
            @info "Solver = $sol"
            if plt_run 
                modplot = 1
            end
            method = findfirst(==(sol), method_map)
            solver_nlsol.solve_energy!(atmos, sol_type=sol_type, 
                                dry_convect=dry_convect, condensates=condensates, sens_heat=incl_sens,
                                max_steps=max_steps, conv_atol=conv_atol, conv_rtol=conv_rtol, method=1,
                                stabilise_mlt=stabilise,modplot=modplot)
        else 
            error("Invalid solver requested '$sol'")
        end 
        @info " "

    end 

    @info "Total RT evalulations: $(atmos.num_rt_eval)"

    # Write arrays
    atmosphere.write_pt(atmos,      joinpath(atmos.OUT_DIR,"pt.csv"))
    atmosphere.write_ncdf(atmos,    joinpath(atmos.OUT_DIR,"atm.nc"))
    atmosphere.write_fluxes(atmos,  joinpath(atmos.OUT_DIR,"fl.csv"))

    # Save plots
    @info "Making plots"
    plt_ani && plotting.anim_solver(atmos)
    plt_vmr && plotting.plot_x(atmos,          joinpath(atmos.OUT_DIR,"plot_vmrs.png"))
    plt_cff && plotting.plot_contfunc(atmos,   joinpath(atmos.OUT_DIR,"plot_contfunc.png"))
    plt_tmp && plotting.plot_pt(atmos,         joinpath(atmos.OUT_DIR,"plot_ptprofile.png"), incl_magma=(sol_type==2))
    plt_flx && plotting.plot_fluxes(atmos,     joinpath(atmos.OUT_DIR,"plot_fluxes.png"), incl_mlt=use_mlt, incl_eff=(sol_type==3))
    plt_ems && plotting.plot_emission(atmos,   joinpath(atmos.OUT_DIR,"plot_emission.png"))
    plt_alb && plotting.plot_albedo(atmos,     joinpath(atmos.OUT_DIR,"plot_albedo.png"))

    # Deallocate atmosphere
    @info "Deallocating arrays"
    atmosphere.deallocate!(atmos)

    # Finish up
    runtime = round(time() - tbegin, digits=2)
    @info "Model runtime: $runtime seconds"
    @info "Goodbye"

    return nothing 
end 

# Call main function 
main()
