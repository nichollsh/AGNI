# Contains code for iteratively solving for energy balance

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module solver 

    using Printf
    using Statistics
    include("../socrates/julia/src/SOCRATES.jl")
    include("atmosphere.jl")
    include("phys.jl")

    # Calculate heating rates at cell-centres
    function calc_heat(atmos, sens::Bool)

        dF = atmos.flux_n[1:end-1] - atmos.flux_n[2:end]
        dp = atmos.pl[1:end-1]     - atmos.pl[2:end]

        if sens
            C_d = 0.001  # Turbulent exchange coefficient [dimensionless]
            U = 10.0     # Wind speed [m s-1]

            sens = atmos.layer_cp[end]*atmos.p[end]/(phys.R_gas*atmos.tmp[end]) * C_d * U * (atmos.tstar - atmos.tmpl[end])
            dF[end] += sens
        end

        heat = zeros(Float64, atmos.nlev_c)
        heat += (atmos.layer_grav / atmos.layer_cp * atmos.layer_mmw) * dF/dp # K/s
        heat *= 86400.0 # K/day

        return heat
    end

    # Calculates the time-step at each level for a given accuracy
    function calc_stepsize(atmos, heat, dt_min, dt_max, dtmp_step_frac)

        @. expect_heat = abs(heat)
        clamp!(expect_heat, 1e-30, Inf) # Prevent division by zero
        dt = dtmp_step_frac * atmos.tmp / expect_heat

        dt[end] *= 0.5

        return dt
    end

    # Iterates the atmosphere at each level
    function temperature_step!(atmos, 
                                heat::Array, dt::Array, dtmp_clip::Float64, 
                                fixed_bottom::Bool, smooth_width::Int64, 
                                dryadj_steps::Int64, h2oadj_steps::Int64)

        bot_old_e = atmos.tmpl[end]

        # Apply radiative heating temperature change
        dtmp = heat * dt 
        clamp!(dtmp, -1.0 * dtmp_clip, dtmp_clip)
        atm.tmp += dtmp

        adj_changed = 0

        # Dry convective adjustment
        if dryadj_steps > 0
            tmp_before_adj = ones(Float64, atmos.nlev_c) * atmos.tmp
            for _ in 1:dryadj_steps
                DryAdj!(atmos)
            end
            adj_changed += count(x->x>0.0, tmp_before_adj - atmos.tmp) 
        end

        # H2O moist convective adjustment
        # if h2oadj_steps > 0:
        #     tmp_before_adj = copy.deepcopy(atm.tmp)
        #     atm.tmp += moist_adj(atm, 1.0, nb_convsteps=h2oadj_steps)
        #     adj_changed += np.count_nonzero(tmp_before_adj - atm.tmp)
        #     del tmp_before_adj

        # Smooth temperature profile
        # if smooth_width > 1:
        #     atm.tmp = savgol_filter(atm.tmp, smooth_width, 1)

        # Temperature floor (centres)
        clamp!(atmos.tmp, atmos.minT, Inf)
        
        # Interpolate to cell-edge values 
        atmos.tmpl[1]   = atmos.tmp[1]
        atmos.tmpl[end] = atmos.tmp[end]
        for idx in 2:atmos.nlev_l-1
            atmos.tmpl[idx] = 0.5 * (atmos.tmp[idx-1] + atmos.tmp[idx])
        end 

        # Handle boundaries of temperature grid
        # atm.tmpl[0]  = atm.tmp[0]
        if fixed_bottom
            atm.tmpl[end] = bot_old_e
        end

        # Temperature floor (edges)
        clamp!(atmos.tmpl, atmos.minT, Inf)
            
        # Second interpolation back to cell-centres, to prevent grid impriting
        # at medium/low resolutions
        atmos.tmp[:] .= 0.5 .* (atmos.tmpl[2:end] + atmos.tmpl[1:end-1])

        return adj_changed
    end

    # Radiative-convective solver
    function solve_energy!(atmos,
                            rscatter::Bool, 
                            surf_state::Int64=0, surf_value::Float64=350.0, ini_state::Int64=0, 
                            gofast::Bool=true, dry_adjust::Bool=true, h2o_adjust::Bool=false, 
                            sens_heat::Bool=true,
                            verbose::Bool=true, plot::Bool=false )


        println("RCSolver: begin")

        # Run parameters
        steps_max    = 140   # Maximum number of steps
        dtmp_gofast  = 40.0  # Change in temperature below which to stop model acceleration
        wait_adj     = 3     # Wait this many steps before introducing convective adjustment
        modprint     = 10    # Frequency to print when verbose==False

        # Convergence criteria
        dtmp_conv    = 5.0    # Maximum rolling change in temperature for convergence (dtmp) [K]
        drel_dt_conv = 30.0   # Maximum rate of relative change in temperature for convergence (dtmp/tmp/dt) [day-1]
        F_rchng_conv = 0.01   # Maximum relative value of F_loss for convergence [%]
        
        if verbose
            @printf("    convergence criteria       \n")
            @printf("    dtmp_comp   < %.3f K       \n", dtmp_conv)
            @printf("    dtmp/tmp/dt < %.3f K day-1 \n", drel_dt_conv)
            @printf("    F_chng^TOA  < %.3f %%      \n", F_rchng_conv)
            @printf(" \n")
        end

        # Variables
        success = false         # Convergence criteria met
        step = 1                # Current step number
        # atm_hist = []           # Store previous atmosphere states
        F_loss = 1e99           # Flux loss (TOA vs BOA)
        F_TOA_rad = 1e99        # Net upward TOA radiative flux
        dtmp_comp = Inf         # Temperature change comparison
        drel_dt = Inf           # Rate of relative temperature change
        drel_dt_prev  = Inf     # Previous ^
        flag_prev = false       # Previous iteration is meeting convergence
        step_frac = 1e-3        # Step size fraction relative to absolute temperature
        stopfast = false        # Stopping the 'fast' phase?
        # atm_orig = copy.deepcopy(atm)   # Initial atmosphere for plotting

        # Handle surface boundary condition
        if surf_state == 0
            fixed_bottom = false
        elseif surf_state == 1
            fixed_bottom = true
            atmos.tmpl[end] = atmos.T_surf
        elseif surf_state == 2
            fixed_bottom = true
            atmos.tmpl[end] = max(surf_value,atmos.minT)
        else
            error("Invalid surface state for radiative-convective solver")
        end

        # Handle initial state
        if ini_state == 1
            atm.tmp[:]  = atm.tmpl[-1]
            atm.tmpl[:] = atm.tmpl[-1]
        end

        # Main loop
        while (!success) && (step <= steps_max)

            # Validate arrays
            if !(all(isfinite, atmos.tmp) && all(isfinite, atmos.tmpl))
                error("Temperature array contains NaNs")
            end 
            if any(<=(0), atmos.tmp) || any(<=(0), atmos.tmpl)
                error("Temperature array contains negative values")
            end


            # Get heating rates and step size
            atmosphere.radtrans!(atmos, true)
            atmosphere.radtrans!(atmos, false)
            heat = calc_heat(atm, sens_heat && !fixed_bottom)
            dt = calc_stepsize(atm, heat, dt_min, dt_max, step_frac)
            if verbose
                @printf("    dt_max,med  = %.3f, %.3f days", maximum(dt), median(dt))
            end

            # Cancel convective adjustment if disabled
            if ( dry_adjust && (step < wait_adj)) || !dry_adjust 
                dryadj_steps = 0
            end
            if ( h2o_adjust && (step < wait_adj)) || !h2o_adjust
                h2oadj_steps = 0
            end

            # Apply radiative heating rate for full step
            # optionally smooth temperature profile
            # optionally do convective adjustment
            adj_changed = temperature_step!(atmos, heat, dt, 
                                            dtmp_clip, fixed_bottom,
                                            smooth_window,
                                            dryadj_steps, h2oadj_steps)

            

        end # end main loop

        # Print information about the final state
        if !success
            @printf("WARNING: Stopping atmosphere iterations without success")
            @printf("    count_adj   = %d layers   ", adj_changed)
            @printf("    max_change  = %.1f'th lvl ", big_dT_lvl)
            @printf("    dtmp_comp   = %.3f K      ", dtmp_comp)
            @printf("    dtmp/tmp/dt = %.3f day-1  ", drel_dt)
            @printf("    F_chng^TOA  = %.4f %%     ", F_rchng)

        else
            @printf("Convergence criteria met (%d iterations)", step)
        end

    end # end solve_energy

end 
