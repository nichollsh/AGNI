# Core file containing functions for running the model

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module AGNI

    # Include system libraries
    using LoggingExtras
    using Printf
    using TOML

    # Include local jl files (order matters)
    include("phys.jl")
    include("spectrum.jl")
    include("atmosphere.jl")
    include("setpt.jl")
    include("dump.jl")
    include("plotting.jl")
    include("energy.jl")
    include("solver.jl")

    # Import src files
    import .phys
    import .atmosphere
    import .setpt
    import .dump
    import .plotting 
    import .energy
    import .solver 

    # Export 
    # export atmosphere
    # export solver
    # export plotting
    # export energy

    """
    **Setup terminal logging and file logging**

    Arguments:
    - `outpath::String`     output file (empty to disable file logging)
    - `verbosity::Int`      verbosity (0: silent, 1: normal, 2: debug)
    """
    function setup_logging(outpath::String, verbosity::Int)

        # File logging?
        to_file::Bool = !isempty(outpath)

        # Remove old file 
        if to_file
            rm(outpath, force=true)
        end

        # If silent 
        if verbosity==0 
            global_logger(MinLevelLogger(current_logger(), Logging.Error))
            return nothing
        end 

        # Formatting
        color::Int = 39
        level::String = "UNSET"
        term_io::IO = stdout

        # Setup file logger
        if to_file
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
        end

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
                term_io = stderr
                level = "ERROR"
            end 
            # Set color, set bold, print level, unset bold, unset color, message
            @printf(term_io, "[\033[%dm\033[1m %-5s \033[21m\033[0m] %s \n",color, level, args.message)
        end;

        # Combine and set 
        if to_file
            logger_both = TeeLogger(logger_file, logger_term);
        else
            logger_both = logger_term 
        end 
        global_logger(logger_both)

        if verbosity == 1
            disable_logging(Logging.Debug) # disable debug; info only
        end 

        return nothing 
    end 

    """
    **Open and validate config file.**

    Arguments:
    - `cfg_path::String`        path to configuration file

    Returns:
    - `cfg_dict::Dict`          dictionary containing the configuration
    """
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

    """
    **Main function to be called by executable file.**

    This function parses the cfg file provided by the call arguments, and runs 
    AGNI according to these requirements. It will also save and plot the 
    output as requested by the cfg file.
    
    Returns:
    - `return_success::Bool`        flag for model success
    """
    function main()::Bool

        # Record start time 
        tbegin = time()

        # Folder 
        ROOT_DIR = dirname(abspath(PROGRAM_FILE))

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
        if cfg["execution"]["clean_output"] || !(ispath(output_dir) && isdir(output_dir))
            rm(output_dir,force=true,recursive=true)
            mkdir(output_dir)
        end 

        # Logging 
        verbosity::Int = cfg["execution"]["verbosity"]
        setup_logging(joinpath(output_dir, "agni.log"), verbosity)

        # Hello
        @info "Hello"

        # Temp folders for OUTPUT
        dir_fastchem = joinpath(output_dir,"fastchem/")
        dir_frames   = joinpath(output_dir,"frames/")
        rm(dir_fastchem,force=true,recursive=true)
        rm(dir_frames,force=true,recursive=true)

        # Copy configuration file 
        cp(cfg_path, joinpath(output_dir, "agni.cfg"), force=true)

        # Read REQUIRED configuration options from dict 
        @info "Using configuration '$(cfg["title"])'"
        #    planet stuff 
        tmp_surf::Float64      = cfg["planet"]["tmp_surf"]
        instellation::Float64  = cfg["planet"]["instellation"]
        albedo_b::Float64      = cfg["planet"]["albedo_b"]
        albedo_s::Float64      = cfg["planet"]["albedo_s"]
        asf_sf::Float64        = cfg["planet"]["s0_fact"]
        radius::Float64        = cfg["planet"]["radius"]
        zenith::Float64        = cfg["planet"]["zenith_angle"]
        gravity::Float64       = cfg["planet"]["gravity"]
        
        #    composition stuff 
        p_top::Float64                  = cfg["composition"]["p_top"]
        condensates::Array{String,1}    = cfg["composition"]["condensates"]
        chem_type::Int                  = cfg["composition"]["chemistry"]
        p_surf::Float64                 = 0.0
        mf_dict::Dict{String, Float64}  = Dict{String, Float64}()
        mf_path::String                 = ""
        use_all_gases::Bool             = cfg["composition"]["include_all"]
        if haskey(cfg["composition"],"p_surf")
            # set composition using VMRs + Psurf
            if haskey(cfg["composition"],"vmr_dict")
                # from dict in cfg file 
                mf_dict = cfg["composition"]["vmr_dict"]
            elseif haskey(cfg["composition"], "vmr_path") 
                # from csv file to be read-in
                mf_path = cfg["files"]["input_vmr"]
            else
                @error "Misconfiguration: if providing p_surf, must also provide VMRs"
                exit(1)
            end 
            p_surf = cfg["composition"]["p_surf"]

        elseif haskey(cfg["composition"], "p_dict")
            # set composition from partial pressures (converted to mixing ratios)
            pp_dict::Dict{String, Float64} = cfg["composition"]["p_dict"]
            for k in keys(pp_dict)
                p_surf += pp_dict[k]
            end 
            for k in keys(pp_dict)
                mf_dict[k] = pp_dict[k]/p_surf
            end 

        else
            @error "Misconfiguration: must provide either p_dict or p_surf+VMRs"
            exit(1)
        end 
        if chem_type in [1,2,3] 
            if length(condensates)>0
                @error "Misconfiguration: FastChem coupling is incompatible with AGNI condensation scheme"
                exit(1)
            else
                mkdir(dir_fastchem)
            end 
        end 

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
        conv_type::String      = cfg["execution" ]["convection_type"]
        incl_sens::Bool        = cfg["execution" ]["sensible_heat"]
        incl_latent::Bool      = cfg["execution" ]["latent_heat"]
        sol_type::Int          = cfg["execution" ]["solution_type"]
        solvers_cmd::Array     = cfg["execution" ]["solvers"]
        initial_req::Array     = cfg["execution" ]["initial_state"]
        dx_max::Float64        = cfg["execution" ]["dx_max"]
        linesearch::Bool       = cfg["execution" ]["linesearch"]
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
        plt_ani = plt_ani && plt_tmp

        # Read OPTIONAL configuration options from dict
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
        tmp_int::Float64 = 0.0
        if sol_type == 3
            tmp_int = cfg["planet"]["tmp_int"]
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
                                mf_dict, mf_path, 
                                
                                condensates=condensates,
                                flag_gcontinuum=flag_cnt, flag_rayleigh=flag_ray,
                                flag_cloud=flag_cld, flag_aerosol=flag_aer,
                                overlap_method=overlap,
                                skin_d=skin_d, skin_k=skin_k, tmp_magma=tmp_magma,
                                target_olr=target_olr,
                                tmp_int=tmp_int, albedo_s=albedo_s,
                                thermo_functions=thermo_funct,
                                C_d=turb_coeff, U=wind_speed,
                                use_all_gases=use_all_gases
                        )
        atmosphere.allocate!(atmos,star_file)

        # Set PT profile by looping over requests
        # Each request may be a command, or an argument following a command
        num_req::Int = length(initial_req)      # Number of requests
        idx_req::Int = 1                        # Index of current request
        str_req::String = "_unset"              # String of current request
        prt_req::String = "Setting initial T(p): "
        while idx_req <= num_req
            # get command 
            str_req = strip(lowercase(initial_req[idx_req]))
            prt_req *= str_req*", "

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

            elseif str_req == "ncdf"
                # set from NetCDF file 
                idx_req += 1
                setpt.fromncdf!(atmos,initial_req[idx_req])
            
            elseif str_req == "add"
                # add X kelvin from the currently stored T(p)
                idx_req += 1
                setpt.add!(atmos,parse(Float64, initial_req[idx_req]))

            elseif str_req == "sat"
                # check surface supersaturation
                setpt.prevent_surfsupersat!(atmos)
            
            elseif str_req == "con"
                # condensing a volatile 
                idx_req += 1
                setpt.saturation!(atmos, initial_req[idx_req])
                if flag_cld
                    atmosphere.water_cloud!(atmos)
                end

            else 
                @error "Invalid initial state '$str_req'"
                return false
            end 

            atmosphere.calc_layer_props!(atmos, ignore_errors=true)
            
            # iterate
            idx_req += 1
        end 
        @info prt_req[1:end-2]

        # Write initial state
        dump.write_pt(atmos, joinpath(atmos.OUT_DIR,"pt_ini.csv"))

        # Do chemistry on initial composition
        if chem_type in [1,2,3]
            @debug "Initial chemistry"
            atmosphere.chemistry_eq!(atmos, chem_type, true)
        end 

        # Frame dir
        if plt_ani
            @debug "Will animate"
            mkdir(dir_frames)
        end 

        # Solver variables 
        incl_convect::Bool= !isempty(conv_type)
        use_mlt::Bool     = (conv_type == "mlt")
        modplot::Int      = 0
        incl_conduct::Bool = false

        # Loop over requested solvers 
        return_success::Bool = true
        solver_success::Bool = true
        method_map::Array{String,1} = ["newton", "gauss", "levenberg"]
        
        if length(solvers_cmd) == 0  # is empty 
            solvers_cmd = [""]
        end 
        if !isempty(solvers_cmd[end])  # append "no solve" case to end, for calculating cff
            push!(solvers_cmd, "")
        end

        for sol in solvers_cmd 

            sol = strip(lowercase(sol))
            if isempty(sol)
                sol = "none"
            end 
            @info "Solving with '$sol'"

            # No solve - just calc fluxes at the end
            if sol == "none"
                energy.calc_fluxes!(atmos, incl_latent,  
                                    incl_convect, incl_sens, incl_conduct, 
                                    calc_cf=plt_cff)
                @info "    done"
            
            # Nonlinear solver
            elseif (sol in method_map) 
                if plt_run 
                    modplot = 1
                end
                method_idx = findfirst(==(sol), method_map)            
                solver_success = solver.solve_energy!(atmos, sol_type=sol_type, 
                                    conduct=incl_conduct,  chem_type=chem_type,
                                    convect=incl_convect, latent=incl_latent,
                                    sens_heat=incl_sens, max_steps=max_steps, conv_atol=conv_atol,
                                    conv_rtol=conv_rtol, method=method_idx, 
                                    dx_max=dx_max, linesearch=linesearch,
                                    modplot=modplot,save_frames=plt_ani)
                return_success = return_success && solver_success
            else 
                @error "Invalid solver"
                return_success = false 
                break
            end 
            @info " "
        end 

        @info "Total RT evalulations: $(atmos.num_rt_eval)"

        # Write arrays
        @info "Writing results"
        dump.write_pt(atmos,      joinpath(atmos.OUT_DIR,"pt.csv"))
        dump.write_fluxes(atmos,  joinpath(atmos.OUT_DIR,"fl.csv"))
        dump.write_ncdf(atmos,    joinpath(atmos.OUT_DIR,"atm.nc"))

        # Save plots
        @info "Plotting results"
        plt_alb = plt_alb && (flag_cld || flag_ray)

        flag_cld && plotting.plot_cloud(atmos,     joinpath(atmos.OUT_DIR,"plot_cloud.png"))

        plt_ani && plotting.animate(atmos)
        plt_vmr && plotting.plot_vmr(atmos,        joinpath(atmos.OUT_DIR,"plot_vmrs.png"), size_x=600)
        plt_cff && plotting.plot_contfunc(atmos,   joinpath(atmos.OUT_DIR,"plot_contfunc.png"))
        plt_tmp && plotting.plot_pt(atmos,         joinpath(atmos.OUT_DIR,"plot_ptprofile.png"), incl_magma=(sol_type==2))
        plt_flx && plotting.plot_fluxes(atmos,     joinpath(atmos.OUT_DIR,"plot_fluxes.png"), incl_mlt=use_mlt, incl_eff=(sol_type==3), incl_cdct=incl_conduct, incl_latent=incl_latent)
        plt_ems && plotting.plot_emission(atmos,   joinpath(atmos.OUT_DIR,"plot_emission.png"))
        plt_alb && plotting.plot_albedo(atmos,     joinpath(atmos.OUT_DIR,"plot_albedo.png"))

        # Deallocate atmosphere
        @info "Deallocating memory"
        atmosphere.deallocate!(atmos)

        # Temp folders
        if cfg["execution"]["clean_output"]
            @debug "Cleaning output folder"
            # save fastchem outputs 
            if chem_type in [1,2,3]
                cp(joinpath(output_dir,"fastchem","chemistry.dat"),joinpath(output_dir,"fc_gas.dat"), force=true)
                if chem_type in [2,3]
                    cp(joinpath(output_dir,"fastchem","condensates.dat"),joinpath(output_dir,"fc_con.dat"), force=true)
                end 
            end 
            # remove folders
            rm(dir_fastchem,force=true,recursive=true)
            rm(dir_frames,force=true,recursive=true)
        end

        # Finish up
        runtime = round(time() - tbegin, digits=2)
        @info "Model runtime: $runtime seconds"
        @info "Goodbye"

        return return_success 
    end 

end 

