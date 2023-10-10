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
        for i in 2:atmos.nlev_c

            T1 = atmos.tmp[i-1]    # upper layer
            p1 = atmos.p[i-1]

            T2 = atmos.tmp[i]  # lower layer
            p2 = atmos.p[i]
            
            pfact = (p1/p2)^(phys.R_gas / atmos.layer_cp[i])
            
            # If slope dT/dp is steeper than adiabat (unstable), adjust to adiabat
            if T1 < T2*pfact
                Tbar = 0.5 * ( T1 + T2 )
                T2 = 2.0 * Tbar / (1.0 + pfact)
                T1 = T2 * pfact
                atmos.tmp[i-1]   = T1
                atmos.tmp[i] = T2
            end
        end

        # Upward pass
        for i in atmos.nlev_c:2

            T1 = atmos.tmp[i-1]
            p1 = atmos.p[i-1]

            T2 = atmos.tmp[i]
            p2 = atmos.p[i]
            
            pfact = (p1/p2)^(phys.R_gas / atmos.layer_cp[i])

            if T1 < T2*pfact
                Tbar = 0.5 * ( T1 + T2 )
                T2 = 2.0 * Tbar / ( 1.0 + pfact)
                T1 = T2 * pfact
                atmos.tmp[i-1]   = T1
                atmos.tmp[i] = T2 
            end 
        end
        return nothing
    end

    # Naive steam adjustment, single step (as per AEOLUS)
    function adjust_steam!(atmos::atmosphere.Atmos_t)
        
        # Get mole fraction of water 
        gas = "H2O"

        # Skip if no water present
        if !(gas in atmos.gases)
            return 
        end 

        #Downward pass
        for i in range(1,stop=atmos.nlev_c-1, step=1)
            x = atmosphere.get_x(atmos, gas, i)
            pp = atmos.p[i] * x
            if (pp < 1e-10)
                continue
            end 
            Tdew = phys.calc_Tdew(gas, pp)
            if (atmos.tmp[i] < Tdew)
                atmos.tmp[i] = Tdew
            end
        end

        #Upward pass
        for i in range(atmos.nlev_c-1,stop=2, step=-1)
            x = atmosphere.get_x(atmos, gas, i)
            pp = atmos.p[i] * x
            if (pp < 1e-10)
                continue
            end 
            Tdew = phys.calc_Tdew(gas, pp)
            if (atmos.tmp[i] < Tdew)
                atmos.tmp[i] = Tdew
            end
        end
        return nothing
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
            # sensible heat flux (TKE scheme for this 1D case)
            # transports energy from the bottom level (at tmpl[end]) to the bottom node (at tmp[end])
            atmos.flux_sens = atmos.layer_cp[end]*atmos.p[end]/(phys.R_gas*atmos.tmp[end]) * atmos.C_d * atmos.U * (atmos.tmpl[end] - atmos.tmp[end])
            dF[end]   += atmos.flux_sens
        end

        atmos.heating_rate[:] .= 0.0
        for i in 1:atmos.nlev_c
            atmos.heating_rate[i] = (atmos.layer_grav[i] / atmos.layer_cp[i] * atmos.layer_mmw[i]) * dF[i]/dp[i] # K/s
        end
        atmos.heating_rate *= 86400.0 # K/day
        return nothing
    end

    """
    **Run the radiative-convective solver.**

    Time-steps the temperature profile with the heating rates to obtain global and 
    local energy balance.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `surf_state::Int=0`               bottom layer temperature, 0: free | 1: fixed | 2: skin
    - `dry_adjust::Bool=true`           enable dry convective adjustment
    - `h2o_adjust::Bool=false`          enable naive steam convective adjustment
    - `sens_heat::Bool=true`            include sensible heating 
    - `verbose::Bool=false`             verbose output
    - `modplot::Int=0`                  plot frequency (0 => no plots)
    - `modhydro::Int=5`                 solve hydrostatic equation frequency (0 => never)
    - `accel::Bool=true`                enable accelerated fast period at the start 
    - `extrap::Bool=true`               enable extrapolation forward in time 
    - `dt_max::Float64=10.0`            maximum time-step outside of the accelerated phase
    - `max_steps::Int=200`              maximum number of solver steps
    - `min_steps::Int=15`               minimum number of solver steps
    - `dtmp_conv::Float64=5.0`          convergence: maximum rolling change in temperature  (dtmp) [K]
    - `drel_dt_conv::Float64=1.0`       convergence: maximum rate of relative change in temperature (dtmp/tmp/dt) [day-1]
    - `drel_F_conv::Float64=0.2`        convergence: maximum relative change in F_TOA_rad for convergence [%]
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            surf_state::Int=0,
                            dry_adjust::Bool=true, h2o_adjust::Bool=false, 
                            sens_heat::Bool=true,
                            verbose::Bool=true, modplot::Int=0, modhydro::Int=5,
                            accel::Bool=true, extrap::Bool=true,
                            dt_max::Float64=10.0, max_steps::Int=200, min_steps::Int=15,
                            dtmp_conv::Float64=5.0, drel_dt_conv::Float64=1.0, drel_F_conv::Float64=0.2
                            )


        println("RCSolver: begin")

        # Run parameters
        dtmp_accel   = 40.0  # Change in temperature below which to stop model acceleration (needs to be turned off at small dtmp)
        dtmp_extrap  = 20.0  # Change in temperature above which to use model extrapolation (not worth it at small dtmp)
        wait_adj     = 10    # Wait this many steps before introducing convective adjustment
        modprint     = 10    # Frequency to print when verbose==false
        len_hist     = 10    # Number of previous states to store
        H_small      = 1.0   # A small heating rate [K/day]
        p_large      = 5e5   # A large pressure [Pa] 

        if verbose
            @printf("    convergence criteria       \n")
            @printf("    dtmp_comp   < %.3f K       \n", dtmp_conv)
            @printf("    dtmp/tmp/dt < %.3f K day-1 \n", drel_dt_conv)
            @printf("    drel_F_conv < %.3f %%      \n", drel_F_conv)
            @printf(" \n")
        end

        # Validate inputs
        min_steps = max(min_steps,wait_adj+1)
        max_steps = max(min_steps+10, max_steps)
        wait_adj  = max(min(wait_adj, max_steps-2), 1)
        len_hist  = max(min(len_hist, max_steps-1), 5)
        modprint  = max(modprint, 0)
        
        # Tracking
        adj_changed = 0
        F_loss = 1.0e99         # Radiative flux loss (TOA vs BOA)
        H_stat = 1.0e99         # Heating rate statistic
        F_TOA_rad = 1.0e99      # Net upward TOA radiative flux
        F_TOA_pre = 1.0e99      # Previous ^
        F_TOA_rel = 1.0e99      # Relative change in ^
        F_BOA_rad = 0.0         # Net upward BOA radiative flux
        F_OLR_rad = 0.0         # Longwave upward TOA radiative flux
        dtmp_comp = Inf         # Temperature change comparison
        drel_dt = Inf           # Rate of relative temperature change
        drel_dt_prev  = Inf     # Previous ^

        oscil = falses(atmos.nlev_c)  # layers which are oscillating
        heat_prev = zeros(Float64, atmos.nlev_c) # Previous iteration heating rates
        dt = ones(Float64,  atmos.nlev_c) # time-step [days]
        dtmp = zeros(Float64, atmos.nlev_c) # temperature step [K]

        # Variables
        success = false         # Convergence criteria met
        step = 0                # Current step number
        flag_prev = false       # Previous iteration is meeting convergence
        flag_this = false       # Current iteration is meeting convergence
        plt_magma = false       # Include tmp_magma in plots
        stopaccel = false       # stop accelerated period this iter
        step_stopaccel = 99999  # step at which accelerated period was stopped

        if !accel 
            step_stopaccel = 1
        end

        plt_magma = (surf_state == 2)

        # Store previous n atmosphere states
        hist_tmp  = zeros(Float64, (len_hist, atmos.nlev_c)) 
        hist_tmpl = zeros(Float64, (len_hist, atmos.nlev_l)) 

        # Plot initial state 
        if modplot > 0 
            plotting.plot_solver(atmos, @sprintf("%s/solver_monitor_0000.png", atmos.OUT_DIR), incl_magma=plt_magma)
        end 

        # Main loop
        while (!success) && (step <= max_steps)
            step += 1

            # ----------------------------------------------------------
            # Check that solver is happy 
            # ----------------------------------------------------------
            if !(all(isfinite, atmos.tmp) && all(isfinite, atmos.tmpl))
                error("Temperature array contains NaNs")
            end 
            if any(<=(0), atmos.tmp) || any(<=(0), atmos.tmpl)
                error("Temperature array contains negative values")
            end


            # ----------------------------------------------------------
            # Prepare for new iteration
            # ----------------------------------------------------------

            if verbose || ( mod(step,modprint) == 0)
                @printf("    step %d ", step)
            end

            # End of the initial fast period
            if accel && (step > 2) && ( (dtmp_comp < dtmp_accel) || ( step/max_steps > 0.3)  || flag_prev ) 
                accel = false 
                stopaccel = true
            end

            # Stop extrapolation towards end
            if ((step/max_steps > 0.9) && (step > 2)) || flag_prev
                extrap = false
            end 

            # Handle phases
            if accel
                # Fast phase
                if verbose || ( mod(step,modprint) == 0)
                    @printf("(accel) ")
                end
                dtmp_clip = 80.0
                dryadj_steps = 20
                h2oadj_steps = 20
                dt_min_step = 1.0e1
                dt_max_step = 1.0e4
                smooth_width = Int( max(0.1*atmos.nlev_c,2 )) # 10% of levels
                
            else
                # Slow phase
                dtmp_clip = 50.0
                dryadj_steps = 30
                h2oadj_steps = 30
                dt_min_step = 1.0e-6
                dt_max_step = max(dt_max, 1.0e-3)
                smooth_width = 0

                # End of 'fast' period - take average of last two iters
                if stopaccel
                    stopaccel = false
                    atmos.tmpl[:] .= 0.5*(hist_tmpl[end,:] .+ hist_tmpl[end-1,:])
                    atmos.tmp[:]  .= 0.5*(hist_tmp[end,:]  .+ hist_tmp[end-1,:])
                    dt[:]   .= 0.5 * (dt_min_step + dt_max_step)
                    dryadj_steps = 0
                    h2oadj_steps = 0
                    step_stopaccel = step
                end

                # Solver is struggling
                if (step > max_steps * 0.9)
                    dt_max_step *= 0.5
                    dtmp_clip *= 0.8
                end
            end

            if verbose || ( mod(step,modprint) == 0)
                if extrap
                    @printf("(extrap)")
                end 
                @printf("\n")
            end

            # ----------------------------------------------------------
            # Solve hydrostatic equation + get layer properties
            # ---------------------------------------------------------- 
            if (modhydro > 0) && ( (mod(step,modhydro) == 0) || accel)
                atmosphere.solve_hydro!(atmos)
            end
            
            # ----------------------------------------------------------
            # Get radiative fluxes (LW and SW)
            # ---------------------------------------------------------- 
            atmosphere.radtrans!(atmos, true)
            atmosphere.radtrans!(atmos, false)


            # ----------------------------------------------------------
            # Calculate heating rates
            # ---------------------------------------------------------- 
            heat_prev[:] .= atmos.heating_rate[:]
            calc_heat!(atmos, sens_heat)


            # ----------------------------------------------------------
            # Calculate step size for this step
            # ---------------------------------------------------------- 
            for i in 1:atmos.nlev_c
                oscil[i] = false

                # Increment dt
                if accel 
                    dt[i] *= 1.2
                else
                    # large oscillations => reduce step size
                    if (atmos.heating_rate[i]*heat_prev[i] < 0) && (abs(atmos.heating_rate[i]) > H_small)
                        oscil[i] = true
                        dt[i] *= 0.1
                    
                    # low heating rate at high pressure => increase step size
                    elseif ( abs(atmos.heating_rate[i]) < H_small ) && (atmos.pl[i] > p_large) && !stopaccel
                        dt[i] *= 1.2

                    else
                        dt[i] *= 1.1
                    end
                end
            end 

            dt[1:2] .= dt_min_step
            clamp!(dt, dt_min_step, dt_max_step)

            if verbose
                @printf("    dt_max,med  = %.5f, %.5f days \n", maximum(dt), median(dt))
            end

            # ----------------------------------------------------------
            # Apply temperature step
            # ---------------------------------------------------------- 
            # optionally smooth temperature profile
            # optionally do convective adjustment

            bot_old_e = atmos.tmpl[end]
            top_old_e = atmos.tmpl[1]

            # Apply radiative heating temperature change
            dtmp .= atmos.heating_rate .* dt 
            clamp!(dtmp, -1.0 * dtmp_clip, dtmp_clip)
            atmos.tmp += dtmp

            adj_changed = 0
            
            # Dry convective adjustment
            if dryadj_steps > 0 && dry_adjust && (step >= wait_adj)
                tmp_before_adj = ones(Float64, atmos.nlev_c) .* atmos.tmp

                # do adjustment steps
                for _ in 1:dryadj_steps
                    adjust_dry!(atmos)
                end
                
                # check which levels were changed
                for i in 1:atmos.nlev_c
                    if abs(tmp_before_adj[i] - atmos.tmp[i]) > 0.1 
                        adj_changed += 1
                    end
                end 
            end

            # H2O moist convective adjustment
            if h2oadj_steps > 0 && h2o_adjust && (step >= wait_adj)
                tmp_before_adj = ones(Float64, atmos.nlev_c) .* atmos.tmp
                for _ in 1:h2oadj_steps
                    adjust_steam!(atmos)
                end
                for i in 1:atmos.nlev_c
                    if abs(tmp_before_adj[i] - atmos.tmp[i]) > 0.1 
                        adj_changed += 1
                    end
                end
            end

            # Smooth temperature profile (if required)
            if smooth_width > 2
                if mod(smooth_width,2) == 0
                    smooth_width += 1
                end
                atmos.tmp = moving_average.hma(atmos.tmp, smooth_width)
            end 

            # Temperature floor (centres)
            clamp!(atmos.tmp, atmos.tmp_floor, Inf)

            # Interpolate temperature to bulk cell-edge values 
            itp = Interpolator(atmos.p, atmos.tmp)
            atmos.tmpl[2:end-1] .= itp.(atmos.pl[2:end-1])

            # Extrapolate top boundary temperature
            grad_dt = atmos.tmp[1]-atmos.tmp[2]
            grad_dp = atmos.p[1]-atmos.p[2]
            atmos.tmpl[1] = atmos.tmp[1] + grad_dt/grad_dp * (atmos.pl[1] - atmos.p[1])
            atmos.tmpl[1] = dot( [atmos.tmpl[1]; top_old_e] , [0.9; 0.1] )  # limit change

            # Calculate bottom boundary temperature
            if surf_state == 0
                # Extrapolate
                grad_dt = atmos.tmp[end]-atmos.tmpl[end-1]
                grad_dp = atmos.p[end]-atmos.pl[end-1]
                atmos.tmpl[end] = atmos.tmp[end] + grad_dt/grad_dp * (atmos.pl[end] - atmos.p[end])
                weights = [0.7; 0.3]

            elseif surf_state == 1 || step <= step_stopaccel
                # Fixed
                weights = [0.0; 1.0]

            elseif surf_state == 2
                # Conductive skin
                # Set tmpl[end] such that the skin carries the required flux
                atmos.tmpl[end] = atmos.tmp_magma - atmos.flux_n[1] * atmos.skin_d / atmos.skin_k
                atmos.tmpl[end] = max(atmos.tmp_floor, atmos.tmpl[end])

                # don't update T_star if the change is small (prevents oscillations)
                if abs(atmos.tmpl[end]-bot_old_e) < 5.0
                    weights = [0.0; 1.0]  
                else 
                    weights = [0.15; 0.85]  
                end
            else 
                error("Invalid surface state $(surf_state)")
            end 

            atmos.tmpl[end] = dot( [atmos.tmpl[end]; bot_old_e] , weights )  # limit change at surface

            # Second interpolation back to cell-centres 
            itp = Interpolator(atmos.pl, atmos.tmpl)  
            atmos.tmp[:] .= itp.(atmos.p[:])

            # ----------------------------------------------------------
            # Accelerate solver with time-extrapolation
            # ---------------------------------------------------------- 
            # (at higher pressures, when appropriate) 
            ex_lookback     = 9
            ex_lookback     = min(len_hist-1, ex_lookback)
            ex_dtmp_clip    = dtmp_clip * 2.0
            if !accel && extrap && (step_stopaccel+ex_lookback+1 < step < 0.95*max_steps) && ( mod(step,ex_lookback+1) == 0) && (dtmp_comp > dtmp_extrap)
                # don't extrapolate surface when handling skin
                ex_endlvl = atmos.nlev_c
                if surf_state == 2
                    ex_endlvl -= 1
                end
                # extrapolate each level
                for i in 1:atmos.nlev_c 
                    if (atmos.pl[i] > p_large) && !oscil[i]
                        atmos.tmp[i]  += max(-ex_dtmp_clip, min(ex_dtmp_clip, atmos.tmp[i]  - hist_tmp[end-ex_lookback,i]  ))
                        atmos.tmpl[i] += max(-ex_dtmp_clip, min(ex_dtmp_clip, atmos.tmpl[i] - hist_tmpl[end-ex_lookback,i] ))
                    end 
                end 
            end 

            # Temperature floor
            clamp!(atmos.tmp,  atmos.tmp_floor, Inf)
            clamp!(atmos.tmpl, atmos.tmp_floor, Inf)

            # Set tstar 
            atmos.tstar = atmos.tmpl[end]
            
            # ----------------------------------------------------------
            # Calculate solver statistics
            # ---------------------------------------------------------- 
            # Calculate a measure of the relative rate of change in temperature.
            # Measured only at cell-centres because they are not set by interpolation.
            # 90th percentile is used to exclude 'difficult' layers, which can 
            # can be reasonably neglected if the rest of the column is fine.
            if step > 1
                drel_dt_prev = drel_dt
                drel_dt = quantile(abs.(  ((atmos.tmp[2:end] .- hist_tmp[end,2:end])./hist_tmp[end,2:end])./dt[2:end]  ), 0.9)
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

            # Calculate the (change in) flux balance
            F_TOA_pre = F_TOA_rad 
            F_TOA_rad = atmos.flux_n[1]
            F_TOA_rel = abs((F_TOA_rad-F_TOA_pre)/F_TOA_pre)*100.0

            F_BOA_rad = atmos.flux_n[end-1]
            F_OLR_rad = atmos.flux_u_lw[1]
            F_loss    = abs(F_TOA_rad-F_BOA_rad)

            # Calculate the 'typical' heating rate magnitude
            H_stat    = quantile(abs.(atmos.heating_rate[:]), 0.6)

            # --------------------------------------
            # Print debug info
            # -------------------------------------- 
            if verbose
                @printf("    count_adj   = %d layers   \n", adj_changed)
                @printf("    dtmp_comp   = %.3f K      \n", dtmp_comp)
                @printf("    dtmp/tmp/dt = %.3f day-1  \n", drel_dt)
                @printf("    drel_F_TOA  = %.4f %%     \n", F_TOA_rel)
                @printf("    F_rad^TOA   = %.2e W m-2  \n", F_TOA_rad)
                @printf("    F_rad^BOA   = %.2e W m-2  \n", F_BOA_rad)
                @printf("    F_rad^OLR   = %.2e W m-2  \n", F_OLR_rad)
                @printf("    F_rad^loss  = %.2f W m-2  \n", F_loss)
                @printf("    HR_typical  = %.4f K day-1\n", H_stat)
                @printf("\n")
            end

            # --------------------------------------
            # Store current state for plotting, etc.
            # -------------------------------------- 
            for i in 1:len_hist-1
                hist_tmp[i,:]  .= hist_tmp[i+1,:]
                hist_tmpl[i,:] .= hist_tmpl[i+1,:]
            end
            hist_tmp[end,:]  .= atmos.tmp[:]
            hist_tmpl[end,:] .= atmos.tmpl[:]

            # --------------------------------------
            # Make plots 
            # -------------------------------------- 
            if (modplot > 0) && (mod(step,modplot) == 0)
                plotting.plot_solver(atmos, @sprintf("%s/solver_monitor_%04d.png", atmos.OUT_DIR, step), hist_tmpl=hist_tmpl, incl_magma=plt_magma)
            end 

            # --------------------------------------
            # Convergence check
            # -------------------------------------- 
            # Requires that:
            # - minimal temperature change for two iters
            # - minimal rate of temperature change for two iters
            # - minimal change to net radiative flux at TOA for two iters
            # - solver is not in 'fast' mode, as it is unphysical
            flag_prev = flag_this
            flag_this = (dtmp_comp < dtmp_conv) && (drel_dt < drel_dt_conv) && ( F_TOA_rel < drel_F_conv)
            success   = flag_this && flag_prev && !accel && !stopaccel && (step > min_steps)
            
            # --------------------------------------
            # Sleep in order to capture keyboard interrupt
            # -------------------------------------- 
            sleep(1e-5)
            atexit(exit)

        end # end main loop

        # Print information about the final state
        if !success
            @printf("RCSolver: Stopping atmosphere iterations before convergence \n")
        else
            @printf("RCSolver: Convergence criteria met (%d iterations) \n", step)
        end
        @printf("    dtmp_comp   = %.3f K      \n", dtmp_comp)
        @printf("    dtmp/tmp/dt = %.3f day-1  \n", drel_dt)
        @printf("    drel_F_TOA  = %.4f %%     \n", F_TOA_rel)
        @printf("\n")

        @printf("RCSolver: Final radiative fluxes [W m-2] \n")
        @printf("    OLR   = %.2e W m-2         \n", F_OLR_rad)
        @printf("    TOA   = %.2e W m-2         \n", F_TOA_rad)
        @printf("    BOA   = %.2e W m-2         \n", F_BOA_rad)
        @printf("    loss  = %.2f W m-2         \n", F_loss)
        @printf("    loss  = %.2f %%            \n", F_loss/F_BOA_rad*100.0)
        @printf("\n")
    
        # Warn user if there's a sign difference in TOA vs BOA fluxes
        # because this shouldn't be the case
        if F_TOA_rad*F_BOA_rad < 0
            @printf("WARNING: TOA and BOA radiative fluxes have different signs\n")
        end
        return nothing
    end # end solve_energy

end 
