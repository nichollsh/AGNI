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

# Include local jl files
include("socrates/julia/src/SOCRATES.jl")
push!(LOAD_PATH, joinpath(ROOT_DIR,"src"))
import atmosphere
import setpt
import plotting 
import phys
import solver_tstep
import solver_nlsol 

function main()

    tbegin = time()

    # Configuration options
    tstar::Float64         = 2800.0    # Surface temperature [kelvin]
    instellation::Float64  = 2000.0
    albedo_b::Float64      = 0.1
    asf_sf::Float64        = 3.0/8.0
    radius::Float64        = 6.37e6    # metres
    zenith::Float64        = 48.19
    gravity::Float64       = 9.81      # m s-2
    nlev_centre::Int       = 100  
    p_surf::Float64        = 10.0    # bar
    p_top::Float64         = 1e-5      # bar 
    mf_dict                = Dict([
                                ("H2O" , 1.0),
                                # ("CO2" , 0.8),
                                # ("CO" ,  0.1),
                                # ("N2" ,  0.05)
                            ])

    spfile_name   = "res/spectral_files/Frostflow/48/Frostflow.sf"
    star_file     = "res/stellar_spectra/sun.txt"
    output_dir    = "out/"

    # Create output directory
    rm(output_dir,force=true,recursive=true)
    mkdir(output_dir)

    # Setup global logger
    logger_file = FormatLogger(joinpath(output_dir,"agni.log"); append=true) do io, args
        @printf(io, "[%-5s] %s",args.level, args.message)
    end;
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
    logger_both = TeeLogger(logger_file, logger_term);
    global_logger(logger_both)
    disable_logging(Logging.Debug) # disable debug, info only

    # Hello
    @info "Hello\n"

    # Setup atmosphere
    @info "Setting up\n"
    atmos = atmosphere.Atmos_t()
    atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                            spfile_name,
                            instellation, asf_sf, albedo_b, zenith,
                            tstar,
                            gravity, radius,
                            nlev_centre, p_surf, p_top,
                            mf_dict=mf_dict,
                            flag_gcontinuum=true,
                            flag_rayleigh=true,
                            flag_cloud=false,
                            overlap_method=4,
                            skin_d=0.01,
                            skin_k=2.0,
                            tmp_magma=3000.0,
                            tmp_floor=5.0,
                            tint=0.0,
                            thermo_functions=false,
                    )
    atmosphere.allocate!(atmos;stellar_spectrum=star_file,spfile_noremove=true)

    # Set PT profile 
    @info "Setting initial T(p)\n"
    # setpt.fromcsv!(atmos,"pt.csv")
    # setpt.isothermal!(atmos, tstar*0.7)
    # setpt.prevent_surfsupersat!(atmos)
    setpt.dry_adiabat!(atmos)
    setpt.condensing!(atmos, "H2O")
    # setpt.stratosphere!(atmos, 800.0)

    atmosphere.write_pt(atmos, joinpath(atmos.OUT_DIR,"pt_ini.csv"))

    @info "Running model...\n"

    # Call solver(s)
    dry_convect = true
    condensate  = ""
    surf_state  = 3

    # solver_tstep.solve_energy!(atmos, surf_state=surf_state, use_physical_dt=false,
    #                             modplot=100, modprop=5, verbose=true, 
    #                             dry_convect=dry_convect, condensate=condensate,
    #                             accel=true, step_rtol=1.0e-4, step_atol=1.0e-2, dt_max=1000.0,
    #                             max_steps=400, min_steps=100, use_mlt=true)


    solver_nlsol.solve_energy!(atmos, surf_state=surf_state, 
                                dry_convect=dry_convect, condensate=condensate,
                                max_steps=2000, conv_atol=1.0e-1, method=1,
                                stabilise_mlt=true)

    # Calculate LW and SW fluxes (once)
    atmosphere.radtrans!(atmos, true, calc_cf=true)
    atmosphere.radtrans!(atmos, false)

    # Calculate convective fluxes (once)
    # println("MLT: calculating fluxes")
    # atmosphere.mlt!(atmos)

    @info "Total RT evalulations: $(atmos.num_rt_eval)\n"

    # Write arrays
    atmosphere.write_pt(atmos,      joinpath(atmos.OUT_DIR,"pt.csv"))
    atmosphere.write_ncdf(atmos,    joinpath(atmos.OUT_DIR,"atm.nc"))
    atmosphere.write_fluxes(atmos,  joinpath(atmos.OUT_DIR,"fl.csv"))

    # Save plots
    @info "Making plots\n"
    # plotting.anim_solver(atmos)
    plotting.plot_x(atmos,          joinpath(atmos.OUT_DIR,"mf.png"))
    plotting.plot_contfunc(atmos,   joinpath(atmos.OUT_DIR,"cf.png"))
    plotting.plot_pt(atmos,         joinpath(atmos.OUT_DIR,"pt.png"), incl_magma=(surf_state==2))
    plotting.plot_fluxes(atmos,     joinpath(atmos.OUT_DIR,"fl.png"))
    plotting.plot_emission(atmos,   joinpath(atmos.OUT_DIR,"em.png"))
    plotting.plot_albedo(atmos,     joinpath(atmos.OUT_DIR,"al.png"))

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
