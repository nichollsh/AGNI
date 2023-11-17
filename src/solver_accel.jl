# Contains code for obtaining energy balance (simplistic euler method)

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module solver_accel

    # Include libraries
    include("../socrates/julia/src/SOCRATES.jl")

    using Printf
    using Statistics
    using Revise
    # using PCHIPInterpolation
    # using LinearAlgebra

    import atmosphere 
    import phys
    import plotting
    import moving_average

   
    """
    **Solve for radiative-convective equilibrium using accelerated time-stepping.**

    Time-steps the temperature profile with the heating rates to obtain global 
    and local energy balance. Uses a first-order Euler method which works well 
    to initialise the integration but may lead to instabilities. Also implments 
    a two-step Adams-Bashforth integrator.

    Convective adjustment is applied to represent convective energy transport in 
    unstable layers, or alternatively MLT is used to calculate the convective 
    fluxes. Adjustment isn't compatible with the Adams-Bashforth integrator.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `surf_state::Int=1`               bottom layer temperature, 0: free | 1: fixed | 2: skin
    - `dry_convect::Bool=true`          enable dry convection
    - `h2o_convect::Bool=false`         enable naive steam convection (not supported by MLT scheme)
    - `use_mlt::Bool=true`              using mixing length theory to represent convection (otherwise use adjustment)
    - `sens_heat::Bool=false`           include sensible heating 
    - `verbose::Bool=false`             verbose output
    - `modplot::Int=0`                  plot frequency (0 => no plots)
    - `accel::Bool=true`                enable accelerated fast period at the start 
    - `extrap::Bool=false`              enable extrapolation forward in time 
    - `adams::Bool=true`                use Adams-Bashforth integrator
    - `rtol::Bool=7.0e-2`               relative tolerence
    - `atol::Bool=1.0e-1`               absolute tolerence
    - `dt_max::Float64=100.0            maximum time-step outside of the accelerated phase
    - `max_steps::Int=300`              maximum number of solver steps
    - `min_steps::Int=15`               minimum number of solver steps
    - `dtmp_conv::Float64=5.0`          convergence: maximum rolling change in temperature  (dtmp) [K]
    - `drel_dt_conv::Float64=0.5`       convergence: maximum rate of relative change in temperature (dtmp/tmp/dt) [day-1]
    - `drel_F_conv::Float64=0.1`        convergence: maximum relative change in F_TOA_rad for convergence [%]
    - `F_losspct_conv::Float64=2.0`     convergence: maximum flux loss throughout column [%]
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            surf_state::Int=1,
                            dry_convect::Bool=true, h2o_convect::Bool=false, use_mlt::Bool=true, 
                            sens_heat::Bool=false,
                            verbose::Bool=true, modplot::Int=0, 
                            accel::Bool=true, extrap::Bool=false, adams::Bool=true,
                            dt_max::Float64=200.0, max_steps::Int=500, min_steps::Int=100,
                            rtol::Float64=9.0e-3, atol::Float64=1.0e-1,
                            dtmp_conv::Float64=5.0, drel_dt_conv::Float64=0.5, drel_F_conv::Float64=0.1, F_losspct_conv::Float64=20.0
                            )

        println("RCSolver: Begin accelerated timestepping")

        # Run parameters
        dtmp_accel   = 5.0   # Change in temperature below which to stop model acceleration (needs to be turned off at small dtmp)
        dtmp_extrap  = 2.0   # Change in temperature above which to use model extrapolation (not worth it at small dtmp)
        wait_con     = 10    # Wait this many steps before introducing convection
        stop_adj     = 300   # After this many steps, disable adjustment if MLT is being used
        modprint     = 25    # Frequency to print when verbose==false
        len_hist     = 12    # Number of previous states to store
        H_small      = 1.0e-1  # A characteristic heating rate [K/day]
        H_large      = 1.0e+5  # A characteristic heating rate [K/day]
        p_large      = 1e6   # A characteristic pressure [Pa] 

        # Validate inputs
        min_steps = max(min_steps,wait_con+1)
        max_steps = max(min_steps+10, max_steps)
        wait_con  = max(min(wait_con, max_steps-2), 1)
        len_hist  = max(min(len_hist, max_steps-1), 5)
        modprint  = max(modprint, 0)

        if use_mlt && h2o_convect 
            error("MLT convection scheme is not compatible with moist convection")
        end

        if (len_hist < 4) && adams
            error("The value of len_hist is too small ($len_hist < 4)")
        end

        if verbose
            @printf("    convergence criteria       \n")
            @printf("    dtmp_comp   < %.3f K       \n", dtmp_conv)
            @printf("    dtmp/tmp/dt < %.3f K day-1 \n", drel_dt_conv)
            @printf("    drel_F_conv < %.3f %%      \n", drel_F_conv)
            if use_mlt
                @printf("    F_loss_conv < %.3f %%      \n", F_losspct_conv)
            end
            @printf(" \n")
        end
        
        # Tracking
        adj_changed::Int   = 0           # Number of convectively adjusted levels
        F_loss::Float64    = 1.0e99      # Total flux loss (TOA vs BOA-1)
        F_losspct::Float64 = 1.0e99      # Total flux loss, relative (TOA vs BOA-1)
        H_stat::Float64    = 1.0e99      # Heating rate statistic
        F_TOA_rad::Float64 = 1.0e99      # Net upward TOA radiative flux
        F_TOA_pre::Float64 = 1.0e99      # Previous ^
        F_TOA_rel::Float64 = 1.0e99      # Relative change in ^
        F_BOA_rad::Float64 = 0.0         # Net upward BOA radiative flux
        F_OLR_rad::Float64 = 0.0         # Longwave upward TOA radiative flux
        F_TOA_tot::Float64 = 0.0         # Total TOA flux
        F_BOA_tot::Float64 = 0.0         # Total BOA flux
        dtmp_comp::Float64 = Inf         # Temperature change comparison
        drel_dt::Float64   = Inf         # Rate of relative temperature change
        drel_dt_prev::Float64 = Inf      # Previous ^

        oscil       = falses(atmos.nlev_c)          # layers which are oscillating
        dt          = ones(Float64,  atmos.nlev_c)  # time-step [days]
        dtmp        = zeros(Float64, atmos.nlev_c)  # temperature step [K]

        # Variables
        step = 0                # Current step number
        success = false         # Convergence criteria met
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
        hist_tmp  = zeros(Float64, (len_hist, atmos.nlev_c))   # temperature (cc)
        hist_tmpl = zeros(Float64, (len_hist, atmos.nlev_l))   # temperature (ce)
        hist_hr   = zeros(Float64, (len_hist, atmos.nlev_c))   # heating rates (cc)

        # Plot initial state 
        if modplot > 0 
            # file name prefixed with zz to make it appear last in directory
            plotting.plot_solver(atmos, @sprintf("%s/zzframe_0000.png", atmos.OUT_DIR), incl_magma=plt_magma)
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
                dtmp_clip = 65.0
                dryadj_steps = 25
                h2oadj_steps = 25
                dt_min_step = 1.0e-1
                dt_max_step = 130.0
                smooth_width = floor(Int, max(0.1*atmos.nlev_c,2 )) # 10% of levels

            else
                # Slow phase
                dtmp_clip = 30.0
                dryadj_steps = 35
                h2oadj_steps = 35
                dt_min_step = 1.0e-4
                dt_max_step = max(dt_max, 1.0e-3)
                smooth_width = 0

                # End of 'fast' period - take average of last two iters
                if stopaccel
                    stopaccel = false
                    atmos.tmpl[:] .= (hist_tmpl[end,:] .+ hist_tmpl[end-1,:] .+ hist_tmpl[end-2,:]) ./ 3.0
                    atmos.tmp[:]  .= (hist_tmp[end,:]  .+ hist_tmp[end-1,:]  .+ hist_tmp[end-2,:] ) ./ 3.0
                    dt[:]         .= 1.0e-1
                    dryadj_steps = 0 # no adjustment in this step
                    h2oadj_steps = 0 # ^^
                    step_stopaccel = step
                end

                # Disable adjustment if mlt is enabled and a sufficient number of steps are reached
                if (step > stop_adj) && use_mlt
                    dryadj_steps = 0
                    h2oadj_steps = 0
                end

                # Adjustment not ready or it's too late
                if (step < wait_con) || (flag_prev && use_mlt)
                    dryadj_steps = 0
                    h2oadj_steps = 0
                end

                # Solver is struggling
                if (step > max_steps * 0.9)
                    dt_max_step *= 0.5
                    dtmp_clip *= 0.8
                end
            end

            if verbose || ( mod(step,modprint) == 0)
                if extrap && (dtmp_comp > dtmp_extrap)
                    @printf("(extrap)")
                end 
                @printf("\n")
            end

            # ----------------------------------------------------------
            # Get fluxes 
            # ---------------------------------------------------------- 
            atmos.flux_tot[:] .= 0.0

            # Radiation
            atmosphere.radtrans!(atmos, true)
            atmosphere.radtrans!(atmos, false)
            atmos.flux_tot += atmos.flux_n

            # Convection
            if use_mlt && (step >= wait_con) 
                atmosphere.mlt!(atmos)
                atmos.flux_tot += atmos.flux_c
            end

            # Turbulence
            if sens_heat
                atmosphere.sensible!(atmos)
                atmos.flux_tot[end] += atmos.flux_sens
            end

            # ----------------------------------------------------------
            # Calculate heating rates
            # ---------------------------------------------------------- 
            atmosphere.calc_hrates!(atmos)

            # ----------------------------------------------------------
            # Calculate step size for this step
            # ---------------------------------------------------------- 
            for i in 1:atmos.nlev_c
                oscil[i] = false

                # Increment dt
                if accel 
                    dt[i] *= 1.2

                elseif (step > step_stopaccel+2)
                    # temperature oscillations  => reduce step size
                    if (atmos.heating_rate[i]*hist_hr[end,i] < 0) || (abs(atmos.heating_rate[i]) > H_large) 
                        oscil[i] = true
                        dt[i] *= 0.02
                    else
                        # Set timestep using tolerences (if not oscillating)
                        if (abs(atmos.tmp[i] - hist_tmp[end,i]) < rtol*abs(hist_tmp[end,i]) + atol) 
                            dt[i] *= 1.05
                        else
                            dt[i] *= 0.95
                        end
                    end
                end

            end 

            # for i in 1:atmos.nlev_c
            #     oscil[i] = false

            #     # Increment dt
            #     if accel 
            #         dt[i] *= 1.2
            #     else
            #         # large oscillations => reduce step size
            #         if (atmos.heating_rate[i]*hist_hr[end,i] < 0)
            #             oscil[i] = true
            #             dt[i] *= 0.1
            #         elseif (abs(atmos.heating_rate[i]) < H_small)
            #             dt[i] *= 1.05
            #         end
            #     end
            # end 

            dt[1:2] .= dt_min_step
            clamp!(dt, dt_min_step, dt_max_step)

            display(dt[1:2:end])

            if verbose
                @printf("    dt_max,med  = %.5f, %.5f days \n", maximum(dt), median(dt))
            end

            # ----------------------------------------------------------
            # Apply temperature step
            # ---------------------------------------------------------- 
            # optionally smooth temperature profile
            # optionally do convective adjustment
            
            dtmp[:] .= 0.0

            # Use Adams-Moulton two-step scheme
            if adams && (step > 3)
                # https://john-s-butler-dit.github.io/NumericalAnalysisBook/Chapter%2004%20-%20Multistep%20Methods/401_Adams%20Bashforth%20Example.html
                for i in 1:atmos.nlev_c 
                    dtmp[i] = 3.0 * atmos.heating_rate[i] - hist_hr[end,i]
                    dtmp[i] *= dt[i] / 2.0
                end 

            # Use Euler one-step scheme
            else 
                dtmp .= atmos.heating_rate .* dt 
            end 
            clamp!(dtmp, -1.0 * dtmp_clip, dtmp_clip)
            atmos.tmp += dtmp

            # Dry convective adjustment
            adj_changed = 0
            if dryadj_steps > 0 && dry_convect
                tmp_before_adj = ones(Float64, atmos.nlev_c) .* atmos.tmp

                # do adjustment steps
                for _ in 1:dryadj_steps
                    atmosphere.adjust_dry!(atmos)
                end
                
                # check which levels were changed
                for i in 1:atmos.nlev_c
                    if abs(tmp_before_adj[i] - atmos.tmp[i]) > 0.1 
                        adj_changed += 1
                    end
                end 
            end

            # H2O moist convective adjustment
            if h2oadj_steps > 0 && h2o_convect
                tmp_before_adj = ones(Float64, atmos.nlev_c) .* atmos.tmp
                for _ in 1:h2oadj_steps
                    atmosphere.adjust_steam!(atmos)
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
            clamp!(atmos.tmp, atmos.tmp_floor, atmos.tmp_ceiling)

            # Set cell-edge values (particularly important for handling conductive skin)
            if (surf_state == 2) && (step <= 10)
                atmosphere.set_tmpl_from_tmp!(atmos, 1, limit_change=true)  # don't allow tmpl to drift during accelerated phase
            else 
                atmosphere.set_tmpl_from_tmp!(atmos, surf_state, limit_change=true)
            end
            

            # ----------------------------------------------------------
            # Accelerate solver with time-extrapolation
            # ---------------------------------------------------------- 
            # (at higher pressures, when appropriate) 
            ex_lookback     = 20
            ex_lookback     = min(len_hist-1, ex_lookback)
            ex_dtmp_clip    = dtmp_clip
            if !accel && extrap && ( mod(step,ex_lookback+1) == 0) && (dtmp_comp > dtmp_extrap)
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
            clamp!(atmos.tmp,  atmos.tmp_floor, atmos.tmp_ceiling)
            clamp!(atmos.tmpl, atmos.tmp_floor, atmos.tmp_ceiling)

            # Set tstar 
            if surf_state != 1
                atmos.tstar = atmos.tmpl[end]
            end
            
            # ----------------------------------------------------------
            # Calculate solver statistics
            # ---------------------------------------------------------- 
            # Calculate a measure of the relative rate of change in temperature.
            # Measured only at cell-centres because they are not set by interpolation.
            # 90th percentile is used to exclude 'difficult' layers, which can 
            # can be reasonably neglected if the rest of the column is fine.
            if step > 1
                drel_dt_prev = drel_dt
                drel_dt = quantile(abs.(  ((atmos.tmp[2:end] .- hist_tmp[end,2:end])./hist_tmp[end,2:end])./dt[2:end]  ), 0.95)
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
            F_BOA_rad = atmos.flux_n[end]
            F_OLR_rad = atmos.flux_u_lw[1]

            F_TOA_tot = atmos.flux_tot[1]
            F_BOA_tot = atmos.flux_tot[end]

            F_loss    = F_TOA_tot-F_BOA_tot
            F_losspct = abs(F_loss/F_TOA_tot*100.0)
            
            # Calculate the 'typical' heating rate magnitude
            H_stat    = quantile(abs.(atmos.heating_rate[:]), 0.95)

            # --------------------------------------
            # Print debug info (some are disabled depending on the convection scheme)
            # -------------------------------------- 
            if verbose
                @printf("    dtmp_comp   = %+.3f K      \n", dtmp_comp)
                @printf("    dtmp/tmp/dt = %+.3f day-1  \n", drel_dt)
                @printf("    drel_F_TOA  = %+.4f %%     \n", F_TOA_rel)

                if use_mlt
                    @printf("    F_pctloss   = %+.4f %%     \n", F_losspct)
                end

                @printf("    F_rad^TOA   = %+.2e W m-2  \n", F_TOA_rad)
                @printf("    F_rad^BOA   = %+.2e W m-2  \n", F_BOA_rad)
                @printf("    F_rad^OLR   = %+.2e W m-2  \n", F_OLR_rad)
                @printf("    F_tot^loss  = %+.2f W m-2  \n", F_loss)
                @printf("    HR_typical  = %+.4f K day-1\n", H_stat)

                if !use_mlt
                    @printf("    count_adj   = %d layers   \n", adj_changed)
                end
                @printf("\n")
            end

            # --------------------------------------
            # Store current state for plotting, etc.
            # -------------------------------------- 
            for i in 1:len_hist-1
                hist_tmp[i,:]  .= hist_tmp[i+1,:]
                hist_tmpl[i,:] .= hist_tmpl[i+1,:]
                hist_hr[i,:]   .= hist_hr[i+1,:]
            end
            hist_tmp[end,:]  .= atmos.tmp[:]
            hist_tmpl[end,:] .= atmos.tmpl[:]
            hist_hr[end,:]   .= atmos.heating_rate[:]

            # --------------------------------------
            # Make plots 
            # -------------------------------------- 
            if (modplot > 0) && (mod(step,modplot) == 0)
                plotting.plot_solver(atmos, @sprintf("%s/solver.png", atmos.OUT_DIR), hist_tmpl=hist_tmpl, incl_magma=plt_magma)
                cp(@sprintf("%s/solver.png", atmos.OUT_DIR), @sprintf("%s/zzframe_%04d.png", atmos.OUT_DIR, step))
                plotting.plot_fluxes(atmos, @sprintf("%s/fluxes.png", atmos.OUT_DIR))
            end 

            # --------------------------------------
            # Convergence check
            # -------------------------------------- 
            # Requires that:
            # - minimal temperature change for two iters
            # - minimal rate of temperature change for two iters
            # - minimal change to net radiative flux at TOA for two iters
            # - minimal flux loss (if MLT is being used for convection)
            # - solver is not in 'fast' mode, as it is unphysical
            flag_prev = flag_this
            flag_this = (dtmp_comp < dtmp_conv) && (drel_dt < drel_dt_conv) && ( F_TOA_rel < drel_F_conv)
            success   = flag_this && flag_prev && !accel && !stopaccel && (step > min_steps) && ( (F_losspct < F_losspct_conv) || !use_mlt)
            
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
        if use_mlt 
            @printf("    loss        = %.2f W m-2    \n", F_loss)
            @printf("    loss        = %.2f %%       \n", F_losspct)
        end
        @printf("\n")

        @printf("RCSolver: Final fluxes [W m-2] \n")
        @printf("    rad_OLR   = %.2e W m-2         \n", F_OLR_rad)
        @printf("    rad_TOA   = %.2e W m-2         \n", F_TOA_rad)
        @printf("    rad_BOA   = %.2e W m-2         \n", F_BOA_rad)
        @printf("    tot_BOA   = %.2e W m-2         \n", F_BOA_tot)
        if (surf_state == 2)
            F_skin = atmos.skin_k / atmos.skin_d * (atmos.tmp_magma - atmos.tstar)
            @printf("    cond_skin = %.2e W m-2         \n", F_skin)
        end
        @printf("\n")
    
        # Warn user if there's a sign difference in TOA vs BOA fluxes
        # because this shouldn't be the case
        if F_TOA_tot*F_BOA_tot < 0
            @printf("WARNING: TOA and BOA total fluxes have different signs\n")
        end
        return nothing

    end # end solve_time


end 
