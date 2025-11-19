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
    import TOML:parsefile

    # Include local jl files (order matters)
    include("phys.jl")
    include("spectrum.jl")
    include("atmosphere.jl")
    include("ocean.jl")
    include("chemistry.jl")
    include("rfm.jl")
    include("setpt.jl")
    include("save.jl")
    include("plotting.jl")
    include("energy.jl")
    include("solver.jl")
    include("load.jl")


    # Import submodules
    import .phys
    import .spectrum
    import .atmosphere
    import .setpt
    import .save
    import .plotting
    import .energy
    import .solver
    import .ocean
    import .chemistry
    import .rfm
    import .load

    # Export submodules (mostly for autodoc purposes)
    export phys
    export spectrum
    export atmosphere
    export setpt
    export save
    export load
    export plotting
    export energy
    export ocean
    export chemistry
    export rfm
    export solver

    const ROOT_DIR::String = abspath(dirname(abspath(@__FILE__)), "../")

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
            logger_file = FormatLogger(outpath; append=false) do io, args
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
            @printf(term_io, "[\033[%dm\033[1m %-5s \033[21m\033[0m] %s \n",
                                color, level, args.message)
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
        cfg_dict = parsefile(cfg_path)

        # check headers
        headers = ["plots", "planet", "files", "execution", "physics", "title"]
        for h in headers
            if !haskey(cfg_dict, h)
                error("Key $h is missing from configuration file at '$cfg_path'")
            end
        end

        # check that output dir is named
        if !haskey(cfg_dict["files"],"output_dir") || (cfg_dict["files"]["output_dir"]=="")
            error("Output directory is missing from configuration file at '$cfg_path'")
        end
        out_path = abspath(cfg_dict["files"]["output_dir"])

        # check if this is a dangerous path
        if ispath(joinpath(out_path, ".git")) || (joinpath(out_path) == pwd()) || samefile(out_path, ROOT_DIR)
            error("Output directory is unsafe")
        end

        # looks good
        return cfg_dict
    end

    """
    **Run AGNI using a configuration dictionary**

    Runs AGNI according to requirements in the configuration dict.
    Will also save and plot the output as requested.

    Arguments:
    - `cfg_dict::Dict`          dictionary containing the configuration

    Returns:
    - `return_success::Bool`    flag for model success
    """
    function run_from_config(cfg::Dict)::Bool
        # Read configuration options from dict
        @info "Using configuration '$(cfg["title"])'"

        # Check that config has these always-required keys
        req_keys = Dict()
        req_keys["files"] = ["output_dir","input_sf"]
        req_keys["planet"] = ["radius","surface_material","instellation","s0_fact","albedo_b","zenith_angle","tmp_surf"]
        req_keys["execution"] = ["num_levels","solver"]
        req_keys["physics"] = ["continua", "rayleigh","cloud","rainout","oceans","overlap_method","thermo_funct","convection_crit"]
        for (k,v) in req_keys
            for kk in v
                if !haskey(cfg[k],kk)
                    @error "Config: missing required key `$k.$kk`"
                    return false
                end
            end
        end

        #    planet structure
        radius::Float64        = cfg["planet"]["radius"]
        gravity::Float64       = 0.0
        mass::Float64          = 0.0
        if haskey(cfg["planet"],"gravity")
            gravity = cfg["planet"]["gravity"]
        end
        if haskey(cfg["planet"],"mass")
            mass = cfg["planet"]["mass"]
            if gravity > 0.0
                # overspecified
                @error "Config: provide `planet.mass` OR `planet.gravity`, not both"
                return false
            end
            # convert mass to gravity
            gravity = phys.grav_accel(mass, radius)
        end

        albedo_s::Float64      = 0.0
        surface_mat::String    = cfg["planet"]["surface_material"]
        if surface_mat == "greybody"
            if !haskey(cfg["planet"],"albedo_s")
                @error "Config: surface is greybody, `albedo_s` must be provided"
                return false
            end
            albedo_s = cfg["planet"]["albedo_s"]
        end

        #    composition stuff (defaults to be overwritten)
        condensates::Array{String,1}        = String[]
        chem::Bool                          = false
        rainout::Bool                       = false
        oceans::Bool                        = false
        p_surf::Float64                     = 0.0
        p_top::Float64                      = 0.0
        pp_dict::Dict{String, Float64}      = Dict{String, Float64}()
        mf_dict::Dict{String, Float64}      = Dict{String, Float64}()
        mf_path::String                     = ""
        metallicities::Dict{String, Float64}= Dict{String, Float64}()
        transparent::Bool                   = false
        real_gas::Bool                      = false

        # transparent atmosphere?
        if haskey(cfg["composition"], "transparent")
            transparent = Bool(cfg["composition"]["transparent"])
        end
        if transparent
            @info "Transparent atmosphere requested"

            # set dummy values
            mf_dict  = Dict("H2"=>1.0)
            p_surf   = 1e-3
            p_top    = 1e-5

        # not transparent
        else
            p_top = Float64(cfg["composition"]["p_top"])
            real_gas = Bool(cfg["physics"]["real_gas"])
            chem = Bool(cfg["physics"]["chemistry"])
            rainout = Bool(cfg["physics"]["rainout"])
            oceans  = Bool(cfg["physics"]["oceans"])
            condensates = cfg["composition"]["condensates"]

            comp_set_by::Int = 0

            # composition set by VMRs and Psurf
            if haskey(cfg["composition"],"p_surf")

                p_surf = cfg["composition"]["p_surf"]

                # from dict in cfg file
                if haskey(cfg["composition"],"vmr_dict")
                    comp_set_by += 1
                    mf_dict = cfg["composition"]["vmr_dict"]
                end

                # from csv file to be read-in
                if haskey(cfg["composition"], "vmr_file")
                    comp_set_by += 1
                    mf_path = cfg["composition"]["vmr_file"]
                end

                # set using metallicity ratios
                if haskey(cfg["composition"], "metallicities")
                    comp_set_by += 1
                    metallicities = cfg["composition"]["metallicities"]

                    # this requires fastchem to be used, to determine the VMRs
                    if !chem
                        @error "Config: must enable FastChem if providing metallicities"
                        return false
                    end

                    # set dummy VMRs for initialising atmosphere struct object
                    mf_dict = Dict("H2"=>0.6, "H2O"=>0.1, "CO2"=>0.1, "N2"=>0.1, "H2S"=>0.1)

                end

                if comp_set_by == 0
                    @error "Config: if providing p_surf, must also provide VMRs"
                    return false
                elseif comp_set_by != 1
                    @error "Config: provide only one of `vmr_dict`, `vmr_file`, `metallicities`"
                    return false
                end

            # composition set by partial pressures
            elseif haskey(cfg["composition"], "p_dict")
                # set composition from partial pressures (converted to mixing ratios)
                pp_dict = cfg["composition"]["p_dict"]
                p_surf = 0.0
                for k in keys(pp_dict)
                    p_surf += pp_dict[k]
                end
                for k in keys(pp_dict)
                    mf_dict[k] = pp_dict[k]/p_surf
                end

                if haskey(cfg["composition"], "vmr_file") || haskey(cfg["composition"], "vmr_dict")
                    @error "Config: cannot provide partial pressures and mixing ratios"
                end

            # most provide either VMR or partial pressures, if not transparent
            else
                @error "Must provide either `p_dict` OR `p_surf` in config"
                @error "    Neglect these keys only when transparent=true"
                return false
            end
        end

        #    chemistry
        if chem &&  transparent
            @error "Config: chemistry is incompatible with transparent atmosphere mode"
            return false
        end

        #    RFM radtrans
        rfm_parfile::String = atmosphere.UNSET_STR
        rfm_wn_min::Float64 = 4000.0
        rfm_wn_max::Float64 = 4020.0
        if haskey(cfg["files"],"rfm_parfile")
            rfm_parfile = cfg["files"]["rfm_parfile"]
            if haskey(cfg["execution"],"rfm_wn_min") && haskey(cfg["execution"],"rfm_wn_max")
                rfm_wn_min = Float64(cfg["execution"]["rfm_wn_min"])
                rfm_wn_max = Float64(cfg["execution"]["rfm_wn_max"])
            else
                @error "Config: RFM calculation enabled (rfm_parfile=$rfm_parfile)"
                @error "        You must also provide `rfm_wn_min` AND `rfm_wn_max`"
                return false
            end
        end

        #    double grey radtrans opacities
        κ_grey_lw::Float64 = 1e-2  # this will be over-written
        κ_grey_sw::Float64 = 1e-2  # ^
        if (lowercase(cfg["files"]["input_sf"]) == "greygas") || cfg["execution"]["grey_start"]
            if all(k in keys(cfg["physics"]) for k in ["grey_sw","grey_lw"])
                κ_grey_lw = Float64(cfg["physics"]["grey_lw"])
                κ_grey_sw = Float64(cfg["physics"]["grey_sw"])
            else
                @error "Config: Grey-gas calculation enabled but opacities are not set"
                @error "        You must also provide `grey_lw` AND `grey_sw`"
                return false
            end
        end

        # star stuff
        star_file::String   = cfg["files" ]["input_star"]
        star_Teff::Float64  = -1.0
        if "star_Teff" in keys(cfg["planet"])
            star_Teff = Float64(cfg["planet"]["star_Teff"])
            if !(lowercase(cfg["files"]["input_star"]) == "blackbody") && (star_Teff>1)
                @warn "Stellar temperature is set by user, but blackbody star isn't enabled"
            end
        end

        #    solver stuff
        incl_convect::Bool     = cfg["physics"]["convection"]
        incl_conduct::Bool     = cfg["physics"]["conduction"]
        incl_advect::Bool      = cfg["physics"]["advection"]
        incl_sens::Bool        = cfg["physics"]["sensible_heat"]
        incl_latent::Bool      = cfg["physics"]["latent_heat"]
        sol_type::Int          = cfg["execution"]["solution_type"]
        conv_atol::Float64     = cfg["execution"]["converge_atol"]
        conv_rtol::Float64     = cfg["execution"]["converge_rtol"]
        perturb_all::Bool      = cfg["execution"]["perturb_all"]

        #    plotting stuff
        plt_tmp::Bool          = cfg["plots"]["temperature"]
        plt_ani::Bool          = cfg["plots"]["animate"]
        plt_ani = plt_ani && plt_tmp

        #    latent heating
        if incl_latent && !rainout
            @error "Config: must enable rainout if also including latent heating"
            return false
        end

        # Read OPTIONAL configuration options from dict
        #     sensible heat at the surface
        roughness::Float64 = 0.001; windspeed::Float64 = 2.0
        if incl_sens && !transparent
            if ! all(k in keys(cfg["planet"]) for k in ["roughness","wind_speed"])
                @error "Config: sensible heating included"
                @error "        you must provide `planet.roughness` and `planet.wind_speed`"
                return false
            end
            roughness  = cfg["planet"]["roughness"]
            windspeed = cfg["planet"]["wind_speed"]
        end
        #     conductive skin case
        skin_k::Float64=2.0; skin_d::Float64=0.1; tmp_magma::Float64=3000.0  # will be overwritten
        if sol_type == 2
            if ! all(k in keys(cfg["planet"]) for k in ["skin_k","skin_d","tmp_magma"])
                @error "Config: solution type $sol_type selected"
                @error "        you must provide `planet.skin_k`, `skin_d`, `tmp_magma`"
                return false
            end
            skin_k      = cfg["planet"]["skin_k"]
            skin_d      = cfg["planet"]["skin_d"]
            tmp_magma   = cfg["planet"]["tmp_magma"]
        end
        #     effective temperature case
        flux_int::Float64 = 0.0
        if sol_type == 3
            if ! haskey(cfg["planet"],"flux_int")
                @error "Config: solution type $sol_type selected"
                @error "        you must provide `planet.flux_int`"
                return false
            end
            flux_int = cfg["planet"]["flux_int"]
        end
        #     target OLR case
        target_olr::Float64 = 250.0
        if sol_type == 4
            if ! haskey(cfg["planet"],"target_olr")
                @error "Config: solution type $sol_type selected"
                @error "        you must provide `planet.target_olr`"
                return false
            end
            target_olr = cfg["planet"]["target_olr"]
        end

        # Output folder
        output_dir = abspath(cfg["files"]["output_dir"])

        # Optional IO folder
        io_dir::String = atmosphere.UNSET_STR
        if haskey(cfg["files"], "io_dir")
            io_dir = cfg["files"]["io_dir"]
        end

        # Create atmosphere structure
        @debug "Instantiate atmosphere"
        atmos = atmosphere.Atmos_t()

        # Setup atmosphere
        @debug "Setup atmosphere "
        atmosphere.setup!(atmos, ROOT_DIR,output_dir,
                                String(cfg["files" ]["input_sf"]),
                                Float64(cfg["planet"]["instellation"]),
                                Float64(cfg["planet"]["s0_fact"]),
                                Float64(cfg["planet"]["albedo_b"]),
                                Float64(cfg["planet"]["zenith_angle"]),
                                Float64(cfg["planet"]["tmp_surf"]),
                                gravity, radius,
                                Int(cfg["execution"]["num_levels"]),
                                p_surf,
                                p_top,
                                mf_dict, mf_path;

                                IO_DIR=io_dir,
                                condensates=condensates,
                                metallicities=metallicities,
                                flag_gcontinuum   = cfg["physics"]["continua"],
                                flag_rayleigh     = cfg["physics"]["rayleigh"],
                                flag_cloud        = cfg["physics"]["cloud"],
                                overlap_method    = cfg["physics"]["overlap_method"],
                                real_gas          = real_gas,
                                thermo_functions  = cfg["physics"]["thermo_funct"],
                                use_all_gases     = chem,
                                surf_roughness=roughness, surf_windspeed=windspeed,
                                skin_d=skin_d, skin_k=skin_k, tmp_magma=tmp_magma,
                                target_olr=target_olr,
                                flux_int=flux_int,
                                surface_material=surface_mat, albedo_s=albedo_s,
                                mlt_criterion=only(cfg["physics"]["convection_crit"][1]),
                                rfm_parfile=rfm_parfile,
                                κ_grey_lw=κ_grey_lw,
                                κ_grey_sw=κ_grey_sw
                        ) || return false

        # Allocate atmosphere
        atmosphere.allocate!(atmos,star_file; stellar_Teff=star_Teff) || return false

        # Set temperatures as appropriate
        if transparent
            # Make the atmosphere transparent. Also sets pressures to be small and turns
            #    off the gas opacity, scattering, continuum, etc.
            @debug "Making atmosphere transparent"
            atmosphere.make_transparent!(atmos)
        else
            # Set T(p) by looping over requests in the config.
            #     Each request may be a command, or an argument following a command
            setpt.request!(atmos, cfg["execution"]["initial_state"]) || return false
        end

        # Write initial state
        save.write_profile(atmos, joinpath(atmos.OUT_DIR,"prof_initial.csv"))

        # Run chemistry in the first instance
        chemistry.calc_composition!(atmos, oceans, chem, rainout)

        # Frame dir
        if plt_ani
            @debug "Will animate"
            rm(atmos.FRAMES_DIR,force=true,recursive=true)
            mkdir(atmos.FRAMES_DIR)
        end

        # Loop over requested solvers
        solver_success::Bool = true
        return_success::Bool = true
        allowed_solvers::Array{String,1} = ["newton", "gauss", "levenberg", "jacobinewton"]
        sol = strip(lowercase(cfg["execution"]["solver"]))
        if isempty(sol)
            sol = "none"
        end
        @info "Solving with '$sol'"

        # No solve - just calc fluxes at the end
        if sol == "none"
            energy.calc_fluxes!(atmos, radiative=true, latent_heat=incl_latent,
                                convective=incl_convect, sens_heat=incl_sens,
                                conductive=incl_conduct, advective=incl_advect,
                                calc_cf=Bool(cfg["plots"]["contribution"]))

        # Transparent atmosphere solver
        elseif sol == "transparent"
            solver_success = solver.solve_transparent!(atmos, sol_type=sol_type,
                                conv_atol=conv_atol,
                                conv_rtol=conv_rtol,
                                max_steps=Int(cfg["execution"]["max_steps"]))

        # Use the nonlinear requested solver
        elseif sol in allowed_solvers
            modplot::Int = 0
            if cfg["plots"]["at_runtime"]
                modplot = 1
            end
            method_idx = findfirst(==(sol), allowed_solvers)
            solver_success = solver.solve_energy!(atmos, sol_type=sol_type,
                                conduct=incl_conduct, chem=chem,
                                convect=incl_convect, latent=incl_latent,
                                sens_heat=incl_sens, advect=incl_advect,
                                max_steps=Int(cfg["execution"]["max_steps"]),
                                max_runtime=Float64(cfg["execution"]["max_runtime"]),
                                conv_atol=conv_atol,
                                conv_rtol=conv_rtol,
                                method=Int(method_idx),
                                rainout=rainout,
                                oceans=oceans,
                                dx_max=Float64(cfg["execution"]["dx_max"]),
                                ls_method=Int(cfg["execution"]["linesearch"]),
                                easy_start=Bool(cfg["execution"]["easy_start"]),
                                grey_start=Bool(cfg["execution"]["grey_start"]),
                                modplot=modplot,
                                save_frames=plt_ani,
                                perturb_all=perturb_all
                                )

        # Invalid selection
        else
            @error "Invalid solver '$sol'"
            return_success = false
        end

        return_success &= solver_success
        @info "    done"
        @info "Total radiative transfer evaluations: $(atmos.num_rt_eval)"

        # RFM calculation?
        if atmos.flag_rfm
            @info "Running RFM line-by-line radiative transfer..."
            rfm.run_rfm(atmos, rfm_wn_min, rfm_wn_max)
            @info "    done"
        end

        # Print information about ocean formation, if any
        for c in atmos.condensates
            if atmos.ocean_tot[c] > eps(0.0)
                @info @sprintf("Surface liquid %s mass: %.2e kg/m^2", c, atmos.ocean_tot[c])
            end
        end
        atmos.ocean_layers = ocean.dist_surf_liq(atmos.ocean_tot,
                                                    atmos.ocean_ob_frac,
                                                    atmos.ocean_cs_height,
                                                    atmos.rp)
        atmos.ocean_topliq = ocean.get_topliq(atmos.ocean_layers)
        atmos.ocean_maxdepth = ocean.get_maxdepth(atmos.ocean_layers)
        atmos.ocean_areacov = ocean.get_areacov(atmos.ocean_layers, atmos.ocean_ob_frac)
        @info @sprintf("Oceans cover %d%% area, max depth %g km",
                            atmos.ocean_areacov*100, atmos.ocean_maxdepth/1e3)

        # Write arrays
        @info "Writing results"
        save.write_ncdf(atmos,    joinpath(atmos.OUT_DIR,"atm.nc"))

        # Save plots
        @info "Plotting results"
        if !transparent
            plt_ani && plotting.animate(atmos.OUT_DIR, atmos.FRAMES_DIR)
            cfg["plots"]["cloud"] && \
                plotting.plot_cloud(atmos,     joinpath(atmos.OUT_DIR,"plot_cloud.png"))
            cfg["plots"]["mixing_ratios"] && \
                plotting.plot_vmr(atmos, joinpath(atmos.OUT_DIR,"plot_vmrs.png"), size_x=600)
            cfg["plots"]["contribution"]  && \
                plotting.plot_contfunc1(atmos, joinpath(atmos.OUT_DIR,"plot_contfunc1.png"))
            cfg["plots"]["height"] && \
                plotting.plot_radius(atmos, joinpath(atmos.OUT_DIR,"plot_radius.png"))
        end
        plt_tmp && \
            plotting.plot_pt(atmos, joinpath(atmos.OUT_DIR,"plot_ptprofile.png"), incl_magma=(sol_type==2))
        cfg["plots"]["fluxes"] && \
            plotting.plot_fluxes(atmos, joinpath(atmos.OUT_DIR,"plot_fluxes.png"),
                                    incl_mlt=incl_convect, incl_eff=(sol_type==3),
                                    incl_cdct=incl_conduct, incl_latent=incl_latent,
                                    incl_advect=incl_advect)
        cfg["plots"]["emission"] && \
            plotting.plot_emission(atmos, joinpath(atmos.OUT_DIR,"plot_emission.png"))
        cfg["plots"]["albedo"] && \
            plotting.plot_albedo(atmos, joinpath(atmos.OUT_DIR,"plot_albedo.png"))

        # Deallocate atmosphere
        @info "Deallocating memory"
        atmosphere.deallocate!(atmos)

        return return_success
    end


    """
    **Main function to be called by executable file.**

    This function parses the cfg file provided by the call arguments, sets up the output
    folders, and configures the logger. It then calls the `run_from_config` function.

    Returns:
    - `return_success::Bool`        flag for model success
    """
    function main()::Bool

        # Record start time
        tbegin = time()

        # Variables
        output_dir::String = atmosphere.UNSET_STR
        clean_output::Bool = false

        # Open and validate config file
        cfg_path::String = joinpath(ROOT_DIR, "res", "config", "default.toml")
        if length(ARGS)>0
            cfg_path = ARGS[1]
        end
        if !ispath(cfg_path)
            error("Cannot find configuration file at '$cfg_path'")
        end
        cfg = open_config(cfg_path)

        # Output folder (make if doesn't exist)
        output_dir = abspath(cfg["files"]["output_dir"])
        clean_output = Bool(cfg["execution"]["clean_output"])
        if clean_output || !isdir(output_dir)
            rm(output_dir,force=true,recursive=true)
            mkdir(output_dir)
        end

        # Logging
        verbosity::Int = cfg["execution"]["verbosity"]
        setup_logging(joinpath(output_dir, "agni.log"), verbosity)

        # Hello
        @debug "Hello"

        # Copy configuration file
        cp(cfg_path, joinpath(output_dir, "agni.toml"), force=true)

        # Run the model
        return_success = run_from_config(cfg)

        # Temp folders
        if clean_output
            @debug "Cleaning output folder"

            # chemistry
            dir_fastchem = joinpath(output_dir,"fastchem")
            fsrc = joinpath(dir_fastchem,"chemistry.dat")
            if isfile(fsrc)
                cp(fsrc, joinpath(output_dir,  "fc_gas.dat"), force=true)
            end
            fsrc = joinpath(dir_fastchem,"condensates.dat")
            if isfile(fsrc)
                cp(fsrc, joinpath(output_dir,  "fc_con.dat"), force=true)
            end
            rm(dir_fastchem,force=true,recursive=true)

            # other
            rm(joinpath(output_dir,"frames"),force=true,recursive=true)
            rm(joinpath(output_dir,"rfm"),   force=true,recursive=true)
        end

        # Finish up
        runtime = round(time() - tbegin, digits=2)
        @info @sprintf("Model runtime: %.2f seconds", runtime)
        @debug "Goodbye"

        return return_success
    end

end

