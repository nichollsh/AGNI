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
    using PCHIPInterpolation
    using LinearAlgebra
    
    import atmosphere 
    import phys
    import plotting
    import moving_average


    # Dry convective adjustment, single step
    function adjust_dry!(atmos::atmosphere.Atmos_t)

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

    # Naive steam adjustment, single step (as per AEOLUS)
    function adjust_steam!(atmos::atmosphere.Atmos_t)
        
        # Get mixing ratio of water 
        gas = "H2O"
        mr_h2o = atmosphere.get_mr(atmos, gas)

        # Skip if no water present
        if mr_h2o == 0.0
            return 
        end 

        #Downward pass
        for i in range(1,stop=atmos.nlev_c-1, step=1)
            pp_h2o = atmos.p[i] * mr_h2o
            if (pp_h2o < 1e-10)
                continue
            end 
            Tdew = phys.calc_Tdew(gas, pp_h2o)
            if (atmos.tmp[i] < Tdew)
                atmos.tmp[i] = Tdew
            end
        end

        #Upward pass
        for i in range(atmos.nlev_c-1,stop=2, step=-1)
            pp_h2o = atmos.p[i] * mr_h2o
            if (pp_h2o < 1e-10)
                continue
            end 
            Tdew = phys.calc_Tdew(gas, pp_h2o)
            if (atmos.tmp[i] < Tdew)
                atmos.tmp[i] = Tdew
            end
        end
    
        # Change in temperature is Tmid_cc - Tmid
        dT_conv[:] = (Tmid_cc[:] - atm.tmp[:])/conv_timescale
        return dT_conv
    
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

    # Iterates the atmosphere at each level
    function temperature_step!(atmos::atmosphere.Atmos_t, 
                                dt::Array, dtmp_clip::Float64, 
                                fixed_bottom::Bool, smooth_width::Int64, 
                                dryadj_steps::Int64, h2oadj_steps::Int64)

        bot_old_e = atmos.tmpl[end]
        top_old_e = atmos.tmpl[1]

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
        if h2oadj_steps > 0
            tmp_before_adj = ones(Float64, atmos.nlev_c) .* atmos.tmp
            for _ in 1:h2oadj_steps
                adjust_steam!(atmos)
            end
            adj_changed += count(x->x>0.0, tmp_before_adj - atmos.tmp) 
        end

        # Smooth temperature profile
        if smooth_width > 2
            if mod(smooth_width,2) == 0
                smooth_width += 1
            end
            atmos.tmp = moving_average.hma(atmos.tmp, smooth_width)
        end 

        # Temperature floor (centres)
        clamp!(atmos.tmp, atmos.T_floor, Inf)
        
        # Interpolate to bulk cell-edge values 
        itp = Interpolator(atmos.p, atmos.tmp)
        atmos.tmpl[2:end-1] .= itp.(atmos.pl[2:end-1])

        # Extrapolate top boundary
        dt = atmos.tmp[1]-atmos.tmpl[2]
        dp = atmos.p[1]-atmos.pl[2]
        atmos.tmpl[1] = atmos.tmp[1] + dt/dp * (atmos.pl[1] - atmos.p[1])

        # Limit change at top edge 
        atmos.tmpl[1] = dot( [atmos.tmpl[1]; top_old_e] , [0.6; 0.4] )

        # Calculate bottom boundary
        if fixed_bottom
            # Fixed
            atmos.tmpl[end] = bot_old_e
        else
            # Extrapolate
            dt = atmos.tmp[end]-atmos.tmpl[end-1]
            dp = atmos.p[end]-atmos.pl[end-1]
            atmos.tmpl[end] = atmos.tmp[end] + dt/dp * (atmos.pl[end] - atmos.p[end])
            
            # Limit change at bottom edge 
            atmos.tmpl[1] = dot( [atmos.tmpl[end]; bot_old_e] , [0.6; 0.4] )
        end

        # Second interpolation back to cell-centres 
        itp = Interpolator(atmos.pl, atmos.tmpl)  
        atmos.tmp[:] .= itp.(atmos.p[:])

        # Temperature floor (edges)
        clamp!(atmos.tmpl, atmos.T_floor, Inf)

        return adj_changed
    end

    """
    **Run the radiative-convective solver.**

    Time-steps the temperature profile with the heating rates to obtain global and 
    local energy balance.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `surf_state::Int64=0`             bottom layer temperature, 0: free | 1: T_surf | 2: surf_value
    - `surf_value::Float64=350.0`       bottom layer temperature when `surf_state==2`
    - `dry_adjust::Bool=true`           enable dry convective adjustment
    - `h2o_adjust::Bool=false`          enable naive steam convective adjustment
    - `sens_heat::Bool=true`            include sensible heating 
    - `verbose::Bool=false`             verbose output
    - `modplot::Int=0`                  plot frequency (0 => no plots)
    - `gofast::Bool=true`               enable accelerated fast period at the start 
    - `max_steps::Int64=250`            maximum number of solver steps
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            surf_state::Int64=0, surf_value::Float64=350.0,
                            dry_adjust::Bool=true, h2o_adjust::Bool=false, 
                            sens_heat::Bool=true,
                            verbose::Bool=true, modplot::Int=0, gofast::Bool=true,
                            max_steps::Int64=250)


        println("RCSolver: begin")

        # Run parameters
        dtmp_gofast  = 40.0  # Change in temperature below which to stop model acceleration
        wait_adj     = 10    # Wait this many steps before introducing convective adjustment
        modprint     = 10    # Frequency to print when verbose==false
        len_hist     = 5     # Number of previous states to store

        # Convergence criteria
        dtmp_conv    = 5.0    # Maximum rolling change in temperature for convergence (dtmp) [K]
        drel_dt_conv = 0.1    # Maximum rate of relative change in temperature for convergence (dtmp/tmp/dt) [day-1]
        F_rchng_conv = 0.5   # Maximum relative value of F_loss for convergence [%]
        
        if verbose
            @printf("    convergence criteria       \n")
            @printf("    dtmp_comp   < %.3f K       \n", dtmp_conv)
            @printf("    dtmp/tmp/dt < %.3f K day-1 \n", drel_dt_conv)
            @printf("    F_chng^loss < %.3f %%      \n", F_rchng_conv)
            @printf(" \n")
        end

        # Validate inputs
        wait_adj = max(min(wait_adj, max_steps-2), 1)
        len_hist = max(min(len_hist, max_steps-1), 4)
        modprint = max(modprint, 0)

        # Tracking
        adj_changed = 0
        F_rchng = Inf
        F_loss = 1.0e99         # Flux loss (TOA vs BOA)
        F_loss_prev = 1.0e99    # Previous ^
        F_TOA_rad = 1.0e99      # Net upward TOA radiative flux
        F_BOA_rad = 0.0         # Net upward BOA radiative flux
        F_OLR_rad = 0.0         # Longwave upward TOA radiative flux
        dtmp_comp = Inf         # Temperature change comparison
        drel_dt = Inf           # Rate of relative temperature change
        drel_dt_prev  = Inf     # Previous ^
        ldrel_dt      = ones(Float64, atmos.nlev_c)
        ldrel_dt_prev = ones(Float64, atmos.nlev_c)
        heat_prev = zeros(Float64, atmos.nlev_c) # Previous iteration heating rates
        dt = ones(Float64,  atmos.nlev_c) # time-step [days]

        # Variables
        success = false         # Convergence criteria met
        step = 1                # Current step number

        flag_prev = false       # Previous iteration is meeting convergence
        stopfast = false

        # Store previous n atmosphere states
        hist_tmp  = zeros(Float64, (len_hist, atmos.nlev_c)) 
        hist_tmpl = zeros(Float64, (len_hist, atmos.nlev_l)) 

        # Handle surface boundary condition
        if surf_state == 0
            fixed_bottom = false
        elseif surf_state == 1
            fixed_bottom = true
        elseif surf_state == 2
            fixed_bottom = true
            atmos.tmpl[end] = max(surf_value,atmos.T_floor)
        else
            error("Invalid surface state for radiative-convective solver")
        end

        # Plot initial state 
        if modplot > 0 
            plotting.plot_solver(atmos, @sprintf("%s/solver_monitor_0000.png", atmos.OUT_DIR))
        end 

        # Main loop
        while (!success) && (step <= max_steps)

            # Validate arrays
            if !(all(isfinite, atmos.tmp) && all(isfinite, atmos.tmpl))
                error("Temperature array contains NaNs")
            end 
            if any(<=(0), atmos.tmp) || any(<=(0), atmos.tmpl)
                error("Temperature array contains negative values")
            end

            # End of the initial fast period
            if gofast && ( (dtmp_comp < dtmp_gofast) || ( step/max_steps > 0.4) )
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
                dt_min = 1e1
                dt_max = 1e6
                smooth_window = Int( max(0.1*atmos.nlev_c,2 )) # 10% of levels
                
            else
                # Slow phase
                if verbose || ( mod(step,modprint) == 0)
                    @printf("    step %d \n", step)
                end
                dtmp_clip = 10.0
                dryadj_steps = 20
                h2oadj_steps = 20
                dt_min = 1e-3
                dt_max = 10.0
                smooth_window = 0

                # End of 'fast' period - take average of last two iters
                if stopfast
                    stopfast = false
                    atmos.tmpl[:] .= 0.5*(hist_tmpl[end,:] .+ hist_tmpl[end-1,:])
                    atmos.tmp[:]  .= 0.5*(hist_tmp[end,:]  .+ hist_tmp[end-1,:])
                    dryadj_steps = 0
                    h2oadj_steps = 0
                end

                # Solver is struggling
                if (step > max_steps * 0.8)
                    dt_max *= 0.5
                    dtmp_clip *= 0.8
                end
            end

            # Get fluxes 
            atmosphere.radtrans!(atmos, true)
            atmosphere.radtrans!(atmos, false)

            # Calc heating rates 
            heat_prev[:] .= atmos.heating_rate[:]
            calc_heat!(atmos, sens_heat && !fixed_bottom)

            # Calc step size
            for i in 1:atmos.nlev_c

                # Increase dt
                dt_i = dt[i] * 1.05

                # Set new time-step 
                if gofast 
                    dt[i] = dt_i
                else 
                    # Weight by previous time-step 
                    dt[i] = dot( [dt[i]; dt_i] , [0.2; 0.8] )

                    # oscillating hr => reduce step size
                    if (atmos.heating_rate[i]*heat_prev[i] < 0) 
                        dt[i] *= 0.1
                    end 
                end
            end 

            dt[1] = dt_min
            dt[2] = dt_min
            clamp!(dt, dt_min, dt_max)

            if verbose
                @printf("    dt_max,med  = %.5f, %.5f days \n", maximum(dt), median(dt))
            end

            # Apply convective adjustment if ready
            if ( dry_adjust && (step < wait_adj)) || !dry_adjust
                dryadj_steps = 0
            end 
            if ( h2o_adjust && (step < wait_adj)) || !h2o_adjust
                h2oadj_steps = 0
            end

            # Apply temperature change
            # optionally smooth temperature profile
            # optionally do convective adjustment
            adj_changed = temperature_step!(atmos, dt, 
                                            dtmp_clip, fixed_bottom,
                                            smooth_window,
                                            dryadj_steps, h2oadj_steps)
            
            # Calculate relative rate of change in temperature
            if step > 1
                ldrel_dt_prev[:] .= ldrel_dt[:]
                ldrel_dt[:] .= abs.(  ((atmos.tmp[:] .- hist_tmp[end,1:end])./hist_tmp[end,1:end])./dt  )

                drel_dt_prev = drel_dt
                drel_dt = maximum(ldrel_dt)
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
            F_TOA_rad = atmos.flux_n[1]

            F_BOA_rad = atmos.flux_n[end-1]
            F_OLR_rad = atmos.flux_u_lw[1]

            F_loss_prev = F_loss
            F_loss = abs(F_TOA_rad-F_BOA_rad)
            F_rchng   = abs( (F_loss - F_loss_prev) / F_loss_prev * 100.0)

            # Print debug info to stdout
            if verbose
                @printf("    count_adj   = %d layers   \n", adj_changed)
                @printf("    dtmp_comp   = %.3f K      \n", dtmp_comp)
                @printf("    dtmp/tmp/dt = %.3f day-1  \n", drel_dt)
                @printf("    F_rad^OLR   = %.2e W m-2  \n", F_OLR_rad)
                @printf("    F_rad^TOA   = %.2e W m-2  \n", F_TOA_rad)
                @printf("    F_rad^BOA   = %.2e W m-2  \n", F_BOA_rad)
                @printf("    F_rad^loss  = %.2f W m-2  \n", F_loss)
                @printf("    F_chng^loss = %.4f %%     \n", F_rchng)
                @printf("\n")
            end

            # Update history array
            for i in 1:len_hist-1
                hist_tmp[i,:]  .= hist_tmp[i+1,:]
                hist_tmpl[i,:] .= hist_tmpl[i+1,:]
            end
            hist_tmp[end,:]  .= atmos.tmp[:]
            hist_tmpl[end,:] .= atmos.tmpl[:]

            # Plot current state
            if (modplot > 0) && (mod(step,modplot) == 0)
                plotting.plot_solver(atmos, @sprintf("%s/solver_monitor_%04d.png", atmos.OUT_DIR, step), hist_tmpl=hist_tmpl)
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
            sleep(1e-4)
            atexit(exit)
            step += 1

        end # end main loop

        # Print information about the final state
        if !success
            @printf("RCSolver: Stopping atmosphere iterations without success \n")
            @printf("    count_adj   = %d layers   \n", adj_changed)
            @printf("    dtmp_comp   = %.3f K      \n", dtmp_comp)
            @printf("    dtmp/tmp/dt = %.3f day-1  \n", drel_dt)
            @printf("    F_chng^TOA  = %.4f %%     \n", F_rchng)
            @printf("\n")

        else
            @printf("RCSolver: Convergence criteria met (%d iterations) \n", step)
            @printf("\n")
        end

        @printf("RCSolver: Final radiative fluxes [W m-2] \n")
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
