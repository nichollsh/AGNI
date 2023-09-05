# Contains code for iteratively solving for energy balance

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module solver 

    # Include libraries
    include("../socrates/julia/src/SOCRATES.jl")

    using Printf
    using Statistics
    using Revise
    
    import atmosphere 
    import phys
    import plotting
    import moving_average

    # Dry convective adjustment (single step)
    function adjust_dry!(atmos)

        # Downward pass
        for i in 1:atmos.nlev_c-1
            T1 = atmos.tmp[i]
            p1 = atmos.p[i]

            T2 = atmos.tmp[i+1]
            p2 = atmos.p[i+1]
            
            pfact = (p1/p2)^(phys.R_gas / atmos.layer_cp[i])
            
            # If slope is shallower than adiabat (unstable), adjust to adiabat
            if T1 < T2*pfact
                Tbar = 0.5 * ( T1 + T2 )
                T2 = 2.0 * Tbar / (1.0 + pfact)
                T1 = T2 * pfact
                atmos.tmp[i]   = T1
                atmos.tmp[i+1] = T2
            end
        end

        # Upward pass
        for i in atmos.nlev_c-1:2

            T1 = atmos.tmp[i]
            p1 = atmos.p[i]

            T2 = atmos.tmp[i+1]
            p2 = atmos.p[i+1]
            
            pfact = (p1/p2)^(phys.R_gas / atmos.layer_cp[i])

            if T1 < T2*pfact
                Tbar = 0.5 * ( T1 + T2 )
                T2 = 2.0 * Tbar / ( 1.0 + pfact)
                T1 = T2 * pfact
                atmos.tmp[i]   = T1
                atmos.tmp[i+1] = T2 
            end 
        end
    end

    # Calculate heating rates at cell-centres
    function calc_heat!(atmos, sens::Bool)

        dF = zeros(Float64, atmos.nlev_c)
        dp = zeros(Float64, atmos.nlev_c)

        for i in 1:atmos.nlev_c
            dF[i] = atmos.flux_n[i+1] - atmos.flux_n[i]
            dp[i] = atmos.pl[i+1]     - atmos.pl[i]
        end 

        if sens
            C_d = 0.001  # Turbulent exchange coefficient [dimensionless]
            U = 10.0     # Wind speed [m s-1]
            sens = atmos.layer_cp[end]*atmos.p[end]/(phys.R_gas*atmos.tmp[end]) * C_d * U * (atmos.tstar - atmos.tmpl[end])
            dF[end] += sens
        end

        atmos.heating_rate[:] .= 0.0
        for i in 1:atmos.nlev_c
            atmos.heating_rate[i] = (atmos.layer_grav[i] / atmos.layer_cp[i] * atmos.layer_mmw[i]) * dF[i]/dp[i] # K/s
        end
        atmos.heating_rate *= 86400.0 # K/day
    
    end

    # Calculates the time-step at each level for a given accuracy
    function calc_stepsize(atmos::atmosphere.Atmos_t, 
                            dt_min::Float64, dt_max::Float64, dtmp_step_frac::Float64)
        
        dt = zeros(Float64, length(atmos.heating_rate))

        for i in 1:atmos.nlev_c
            expect_heat = max( abs(atmos.heating_rate[i]) , 1e-30 )
            dt[i] = dtmp_step_frac * atmos.tmp[i] / expect_heat
        end 

        dt[end] *= 0.5

        clamp!(dt, dt_min, dt_max)

        return dt
    end

    # Iterates the atmosphere at each level
    function temperature_step!(atmos::atmosphere.Atmos_t, 
                                dt::Array, dtmp_clip::Float64, 
                                fixed_bottom::Bool, smooth_width::Int64, 
                                dryadj_steps::Int64, h2oadj_steps::Int64)

        bot_old_e = atmos.tmpl[end]

        # Apply radiative heating temperature change
        dtmp = atmos.heating_rate .* dt 
        clamp!(dtmp, -1.0 * dtmp_clip, dtmp_clip)
        atmos.tmp += dtmp

        adj_changed = 0

        # Dry convective adjustment
        if dryadj_steps > 0
            tmp_before_adj = ones(Float64, atmos.nlev_c) .* atmos.tmp
            for _ in 1:dryadj_steps
                adjust_dry!(atmos)
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
        if smooth_width > 2
            if mod(smooth_width,2) == 0
                smooth_width += 1
            end
            atmos.tmp = moving_average.hma(atmos.tmp, smooth_width)
        end 

        # Temperature floor (centres)
        clamp!(atmos.tmp, atmos.minT, Inf)
        
        # Interpolate to cell-edge values 
        for idx in 2:atmos.nlev_l-1
            atmos.tmpl[idx] = 0.5 * (atmos.tmp[idx-1] + atmos.tmp[idx])
        end

        # Extrapolate top boundary
        dt = atmos.tmp[1]-atmos.tmpl[2]
        dp = atmos.p[1]-atmos.pl[2]
        atmos.tmpl[1] = atmos.tmp[1] + dt/dp * (atmos.pl[1] - atmos.p[1])

        # Calculate bottom boundary
        if fixed_bottom
            # Fixed
            atmos.tmpl[end] = bot_old_e
        else
            # Extrapolate
            dt = atmos.tmp[end]-atmos.tmpl[end-1]
            dp = atmos.p[end]-atmos.pl[end-1]
            atmos.tmpl[end] = atmos.tmp[end] + dt/dp * (atmos.pl[end] - atmos.p[end])
        end

        # Second interpolation back to cell-centres 
        atmos.tmp[:] .= 0.5 .* (atmos.tmpl[2:end] + atmos.tmpl[1:end-1])

        # Temperature floor (edges)
        clamp!(atmos.tmpl, atmos.minT, Inf)

        return adj_changed
    end

    # Radiative-convective solver
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            surf_state::Int64=0, surf_value::Float64=350.0, ini_state::Int64=0, 
                            dry_adjust::Bool=true, h2o_adjust::Bool=false, 
                            sens_heat::Bool=true,
                            verbose::Bool=true, plot::Bool=false, gofast::Bool=true )


        println("RCSolver: begin")

        # Run parameters
        steps_max    = 250   # Maximum number of steps
        dtmp_gofast  = 40.0  # Change in temperature below which to stop model acceleration
        wait_adj     = 3   # Wait this many steps before introducing convective adjustment
        modprint     = 10    # Frequency to print when verbose==false
        len_hist     = 5     # Number of previous states to store

        # Convergence criteria
        dtmp_conv    = 5.0    # Maximum rolling change in temperature for convergence (dtmp) [K]
        drel_dt_conv = 5.0   # Maximum rate of relative change in temperature for convergence (dtmp/tmp/dt) [day-1]
        F_rchng_conv = 0.01   # Maximum relative value of F_loss for convergence [%]
        
        if verbose
            @printf("    convergence criteria       \n")
            @printf("    dtmp_comp   < %.3f K       \n", dtmp_conv)
            @printf("    dtmp/tmp/dt < %.3f K day-1 \n", drel_dt_conv)
            @printf("    F_chng^TOA  < %.3f %%      \n", F_rchng_conv)
            @printf(" \n")
        end

        # Tracking
        adj_changed = 0
        F_rchng = Inf
        F_loss = 1.0e99         # Flux loss (TOA vs BOA)
        F_TOA_rad = 1.0e99      # Net upward TOA radiative flux
        F_BOA_rad = 0.0
        F_OLR_rad = 0.0

        # Variables
        success = false         # Convergence criteria met
        step = 1                # Current step number
        dtmp_comp = Inf         # Temperature change comparison
        drel_dt = Inf           # Rate of relative temperature change
        drel_dt_prev  = Inf     # Previous ^
        flag_prev = false       # Previous iteration is meeting convergence
        step_frac = 1e-3        # Step size fraction relative to absolute temperature
        stopfast = false
        dt_min = 1e-4
        dt_max = 1e3

        # Store previous n atmosphere states
        hist_tmp  = zeros(Float64, (len_hist, atmos.nlev_c)) 
        hist_tmpl = zeros(Float64, (len_hist, atmos.nlev_l)) 

        # Handle surface boundary condition
        if surf_state == 0
            fixed_bottom = false
        elseif surf_state == 1
            fixed_bottom = true
            atmos.tmpl[end] = atmos.tstar
        elseif surf_state == 2
            fixed_bottom = true
            atmos.tmpl[end] = max(surf_value,atmos.minT)
        else
            error("Invalid surface state for radiative-convective solver")
        end

        # Handle initial state
        if ini_state == 1
            atmos.tmp[:]  = atmos.tmpl[end]
            atmos.tmpl[:] = atmos.tmpl[end]
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

            # End of the initial fast period
            if gofast && ( (dtmp_comp < dtmp_gofast) || ( step/steps_max > 0.4) )
                gofast = false 
                stopfast = true
            end

             # Set parameters for step size calculation
            if gofast
                # Fast phase
                if verbose || ( mod(step,modprint) == 0)
                    @printf("    step %d (fast) \n", step)
                end
                dtmp_clip = 100.0
                dryadj_steps = 40
                h2oadj_steps = 40
                dt_min = 1e-3
                dt_max = 1e6
                step_frac = 0.1
                smooth_window = Int( max(0.1*atmos.nlev_c,2 )) # 10% of levels
                
            else
                # Slow phase
                if verbose || ( mod(step,modprint) == 0)
                    @printf("    step %d \n", step)
                end
                dtmp_clip = 10.0
                dryadj_steps = 20
                h2oadj_steps = 20
                dt_min = 1e-5
                dt_max = 20.0
                step_frac_max = 5e-3
                smooth_window = 0

                # End of 'fast' period (take average of last two iters)
                if stopfast
                    stopfast = false
                    atmos.tmpl[:] .= 0.5*(hist_tmpl[end,:] .+ hist_tmpl[end-1,:])
                    atmos.tmp[:]  .= 0.5*(hist_tmp[end,:]  .+ hist_tmp[end-1,:])
                    dryadj_steps = 0
                    h2oadj_steps = 0
                end

                # Solver is struggling
                if (step > steps_max * 0.8)
                    step_frac_max *= 0.7
                    dt_max *= 0.5
                    dtmp_clip *= 0.8
                end

                # Adapt the time-stepping accuracy
                if drel_dt < Inf
                    step_frac *= min(max( drel_dt_prev/drel_dt , 0.6 ) , 1.2)
                end 
                step_frac = min(step_frac, step_frac_max)
                if verbose
                    @printf("    step_frac   = %.2e \n", step_frac)
                end
            end

            # Get heating rates and step size
            atmosphere.radtrans!(atmos, true)
            atmosphere.radtrans!(atmos, false)
            calc_heat!(atmos, sens_heat && !fixed_bottom)
            dt = calc_stepsize(atmos,dt_min, dt_max, step_frac)
            if verbose
                @printf("    dt_max,med  = %.3f, %.3f days \n", maximum(dt), median(dt))
            end

            # Apply convective adjustment if ready
            if ( dry_adjust && (step < wait_adj)) || !dry_adjust
                dryadj_steps = 0
            end 
            if ( h2o_adjust && (step < wait_adj)) || !h2o_adjust
                h2oadj_steps = 0
            end

            # Apply radiative heating rate for full step
            # optionally smooth temperature profile
            # optionally do convective adjustment
            adj_changed = temperature_step!(atmos, dt, 
                                            dtmp_clip, fixed_bottom,
                                            smooth_window,
                                            dryadj_steps, h2oadj_steps)

            # Calculate relative rate of change in temperature
            if step > 1
                drel_dt_prev = drel_dt
                drel_dt = maximum(abs.(  ((atmos.tmp[:] .- hist_tmp[end,1:end])./hist_tmp[end,1:end])./dt  ))
            end

            # Calculate maximum average change in temperature (insensitive to oscillations)
            if step > 3
                dtmp_comp = -1
                for i in 1:atmos.nlev_c 
                    tmp_comp_1 = 0.5 * ( atmos.tmp[i]      + hist_tmp[end  ,i] )
                    tmp_comp_2 = 0.5 * ( hist_tmp[end-1,i] + hist_tmp[end-2,i] )
                    dtmp_comp = max(dtmp_comp, abs(tmp_comp_1 - tmp_comp_2))
                end 
            end 

            # Calculate (the change in) flux balance
            F_TOA_rad_prev = F_TOA_rad
            F_TOA_rad = atmos.flux_n[1]
            F_rchng   = abs( (F_TOA_rad - F_TOA_rad_prev) / F_TOA_rad_prev * 100.0)

            F_BOA_rad = atmos.flux_n[end]
            F_OLR_rad = atmos.flux_u_lw[1]
            F_loss = abs(F_TOA_rad-F_BOA_rad)

            # Print debug info to stdout
            if verbose
                @printf("    count_adj   = %d layers   \n", adj_changed)
                @printf("    dtmp_comp   = %.3f K      \n", dtmp_comp)
                @printf("    dtmp/tmp/dt = %.3f day-1  \n", drel_dt)
                @printf("    F_rad^OLR   = %.2e W m-2  \n", F_OLR_rad)
                @printf("    F_rad^TOA   = %.2e W m-2  \n", F_TOA_rad)
                @printf("    F_rad^BOA   = %.2e W m-2  \n", F_BOA_rad)
                @printf("    F_rad^loss  = %.2f W m-2  \n", F_loss)
                @printf("    F_chng^TOA  = %.4f %%     \n", F_rchng)
            end

            # Update history array
            for i in 1:len_hist-1
                hist_tmp[i,:]  .= hist_tmp[i+1,:]
                hist_tmpl[i,:] .= hist_tmpl[i+1,:]
            end
            hist_tmp[end,:]  .= atmos.tmp[:]
            hist_tmpl[end,:] .= atmos.tmpl[:]

            # Plot current state
            # Animate frames with `ffmpeg -framerate 5 -i out/radeqm_monitor_%04d.png -y out/anim.mp4`
            if plot 
                plotting.plot_pt(atmos, @sprintf("out/radeqm_monitor_%04d.png", step))
            end 

            # Convergence check requires that:
            # - minimal temperature change for two iters
            # - minimal rate of temperature change for two iters
            # - minimal change to net radiative flux at TOA for two iters
            # - solver is not in 'fast' mode, as it is unphysical
            flag_this = (dtmp_comp < dtmp_conv) && (drel_dt < drel_dt_conv) && (F_rchng < F_rchng_conv)
            success   = flag_this && flag_prev
            flag_prev = flag_this
            
            # Prepare for next iter
            atexit(exit)
            @printf(" \n")
            step += 1

        end # end main loop

        # Print information about the final state
        if !success
            @printf("RCSolver: Stopping atmosphere iterations without success \n")
            @printf("    count_adj   = %d layers   \n", adj_changed)
            @printf("    dtmp_comp   = %.3f K      \n", dtmp_comp)
            @printf("    dtmp/tmp/dt = %.3f day-1  \n", drel_dt)
            @printf("    F_chng^TOA  = %.4f %%     \n", F_rchng)

        else
            @printf("RCSolver: Convergence criteria met (%d iterations) \n", step)
        end

        @printf("Final radiative fluxes [W m-2] \n")
        @printf("    OLR   = %.2e W m-2         \n", F_OLR_rad)
        @printf("    TOA   = %.2e W m-2         \n", F_TOA_rad)
        @printf("    BOA   = %.2e W m-2         \n", F_BOA_rad)
        @printf("    loss  = %.2f W m-2         \n", F_loss)
    
        # Warn user if there's a sign difference in TOA vs BOA fluxes
        # because this shouldn't be the case
        if F_TOA_rad*F_BOA_rad < 0
            @printf("WARNING: TOA and BOA radiative fluxes have different signs\n")
        end

    end # end solve_energy

end 
