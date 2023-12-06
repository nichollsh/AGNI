# Contains code for obtaining energy balance (simplistic euler method)

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module solver_tstep

    # Include libraries
    include("../socrates/julia/src/SOCRATES.jl")

    using Printf
    using Statistics
    using Revise

    import atmosphere 
    import phys
    import setpt
    import plotting

   
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
    - `surf_state::Int=1`               bottom layer temperature, 0: free | 1: fixed | 2: conductive skin
    - `update_tstar::Bool=false`        update tstar BC to be equal to tmpl[end]
    - `dry_convect::Bool=true`          enable dry convection
    - `condensate::String=""`           condensate to model (if empty, no condensates are modelled)
    - `use_mlt::Bool=true`              using mixing length theory to represent convection (otherwise use adjustment)
    - `sens_heat::Bool=false`           include sensible heating 
    - `modprop::Int=1`                  frequency at which to update thermodynamic properties (0 => never)
    - `verbose::Bool=false`             verbose output
    - `modplot::Int=0`                  plot frequency (0 => no plots)
    - `accel::Bool=true`                enable accelerated fast period at the start 
    - `adams::Bool=true`                use Adams-Bashforth integrator
    - `rtol::Bool=1.0e-4`               relative tolerence
    - `atol::Bool=1.0e-2`               absolute tolerence
    - `dt_max::Float64=500.0            maximum time-step outside of the accelerated phase
    - `max_steps::Int=1000`             maximum number of solver steps
    - `min_steps::Int=300`              minimum number of solver steps
    - `drel_dt_conv::Float64=1.0`       convergence: maximum rate of relative change in temperature (dtmp/tmp/dt) [day-1]
    - `drel_F_conv::Float64=0.1`        convergence: maximum relative change in F_TOA_rad for convergence [%]
    - `F_losspct_conv::Float64=1.0`     convergence: maximum relative local radiative flux loss [%]
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            surf_state::Int=1, update_tstar::Bool=false,
                            dry_convect::Bool=true, condensate::String="", use_mlt::Bool=true,
                            sens_heat::Bool=false, modprop::Int=1,
                            verbose::Bool=true, modplot::Int=0,
                            accel::Bool=true, adams::Bool=true,
                            dt_max::Float64=500.0, max_steps::Int=1000, min_steps::Int=300,
                            rtol::Float64=1.0e-4, atol::Float64=1.0e-2,
                            drel_dt_conv::Float64=1.0, drel_F_conv::Float64=0.1, F_losspct_conv::Float64=1.0
                            )

        println("TSSolver: begin accelerated timestepping")

        # Run parameters
        dtmp_accel   = 15.0   # Change in temperature below which to stop model acceleration (needs to be turned off at small dtmp)
        smooth_stp   = 150    # Number of steps for which to apply smoothing
        smooth       = true   # Currently smoothing?
        wait_con     = 50     # Introduce convection after this many steps, if ^^ is not already true

        modprint     = 50     # Frequency to print when verbose==false
        len_hist     = 10     # Number of previous states to store
        H_large      = 1.0e5  # A characteristic large heating rate [K/day]

        do_condense  = false  # Allow condensation ever? (overwritten according to condensate)
        is_condense  = false  # Is condensation currently enabled?

        if condensate != "" 
            if condensate in atmos.gases
                do_condense = true 
            else 
                error("Invalid condensate ('$condensate')")
            end 
        end 

        if verbose
            modprint=10
        end

        # Validate inputs
        wait_con  = max(wait_con, 1)
        min_steps = max(min_steps,wait_con+100)
        min_steps = max(min_steps,smooth_stp+50)
        max_steps = max(min_steps+50, max_steps)
        len_hist  = max(min(len_hist, max_steps-1), 5)
        modprint  = max(modprint, 1)

        if (len_hist < 4) && adams
            error("The value of len_hist is too small ($len_hist < 4)")
        end

        if verbose
            @printf("    convergence criteria       \n")
            @printf("    step count  > %d           \n", min_steps)
            @printf("    dtmp/tmp/dt < %.3f K day-1 \n", drel_dt_conv)
            @printf("    dF/F (TOA)  < %.3f %%      \n", drel_F_conv)
            @printf("    F_loc loss  < %.3f %%      \n", F_losspct_conv)
            @printf(" \n")
        end
        
        # Tracking
        adj_changed::Int   = 0           # Number of convectively adjusted levels
        F_loss::Float64    = 1.0e99      # Total flux loss (TOA vs BOA-1)
        F_losspct::Float64 = 1.0e99      # Maximum local relative radiative flux loss [%]
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
        start_con          = false       # Convection has been started at any point

        oscil       = falses(atmos.nlev_c)          # layers which are oscillating
        dt          = ones(Float64,  atmos.nlev_c)  # time-step [days]
        dtmp        = zeros(Float64, atmos.nlev_c)  # temperature step [K]

        # Variables
        step      = 0           # Current step number
        success   = false       # Convergence criteria met
        flag_prev = false       # Previous iteration is meeting convergence
        flag_this = false       # Current iteration is meeting convergence
        plt_magma = false       # Include tmp_magma in plots
        stopaccel = false       # stop accelerated period this iter
        rtol_step = rtol        # rtol used this step
        step_stopaccel = 99999  # step at which accelerated period was stopped

        if !accel 
            step_stopaccel = 1
        end

        plt_magma = (surf_state == 2)
        smooth_width = 1

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
        while (!success) && (step < max_steps)
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

            atmos.mask_c[:] .-= 1.0
            atmos.mask_p[:] .-= 1.0

            if mod(step,modprint) == 0
                @printf("    step %d ", step)
            end

            # End of the initial fast period
            if accel && (step > len_hist) && ( (dtmp_comp < dtmp_accel) || ( step/max_steps > 0.3)  || flag_prev ) 
                accel = false 
                stopaccel = true
            end

            # Handle phases
            if accel
                # Fast phase
                if mod(step,modprint) == 0
                    @printf("(accel) ")
                end
                dtmp_clip    = 50.0
                dryadj_steps = 4
                is_condense  = do_condense
                dt_min_step  = 1.0e-1
                dt_max_step  = 120.0
                rtol_step    = rtol
                smooth       = true 
                smooth_width = floor(Int, max(0.20*atmos.nlev_c,2 )) # 20% of levels

            else
                # Slow phase
                dtmp_clip    = 5.0
                dryadj_steps = 16
                is_condense  = do_condense
                dt_min_step  = 5.0e-5
                dt_max_step  = max(dt_max, 1.0e-2)
                rtol_step    = rtol
                
                if smooth && (step >= smooth_stp)
                    smooth = false
                    clamp!(dt, dt_min_step, 1.0)
                end

                # End of 'fast' period - take average of last two iters
                if stopaccel
                    stopaccel = false
                    atmos.tmpl[:] .= (hist_tmpl[end,:] .+ hist_tmpl[end-1,:] .+ hist_tmpl[end-2,:]) ./ 3.0
                    atmos.tmp[:]  .= (hist_tmp[end,:]  .+ hist_tmp[end-1,:]  .+ hist_tmp[end-2,:] ) ./ 3.0
                    dt[:]         .= 1.0e-1
                    dryadj_steps = 0 # no adjustment in this step
                    is_condense  = false 
                    step_stopaccel = step
                end

                # Solver is struggling
                if (step > max_steps * 0.9)
                    dt_max_step *= 0.5
                    dtmp_clip *= 0.8

                # Not struggling
                else
                    if smooth
                        smooth_width = floor(Int, max(0.08*atmos.nlev_c,2 )) # 10% of levels
                        dt_max_step *= 10.0
                        rtol_step *= 12.0
                    else 
                        smooth_width = 0
                    end
                end
            end

            # Introduce convection and condensation schemes
            if !start_con && (step >= wait_con) && (dry_convect || do_condense)
                start_con = true 
                @printf("(intro convect/condense) ")
            end

            if (mod(step,modprint) == 0) && smooth
                @printf("(smooth) ")
            end

            # ----------------------------------------------------------
            # Recalculate thermodynamic properties at each layer (+ height & gravity)
            # ---------------------------------------------------------- 
            if (modprop > 0) && ( (mod(step, modprop) == 0) || (step < 10))
                if mod(step,modprint) == 0
                    @printf("(props) ")
                end
                atmosphere.calc_layer_props!(atmos)
            end 

            if mod(step,modprint) == 0
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

            # Dry convection (MLT)
            if use_mlt && dry_convect && start_con
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
                    if ( (atmos.heating_rate[i]*hist_hr[end,i] < 0) && (hist_hr[end,i]*hist_hr[end-1,i] < 0))  || (abs(atmos.heating_rate[i]) > H_large) 
                        oscil[i] = true
                        dt[i] *= 0.4
                    else
                        # Set timestep
                        if (abs(atmos.tmp[i] - hist_tmp[end,i]) < (rtol_step*abs(hist_tmp[end,i]) + atol) ) 
                            dt[i] *= 1.03
                        else
                            dt[i] /= 1.02
                        end
                    end
                end
            end 

            dt[end] = min(dt[end], dt_max_step*0.95)
            clamp!(dt, dt_min_step, dt_max_step)

            if verbose && (mod(step,modprint) == 0)
                @printf("    dt_max,med  = %.5f, %.5f days \n", maximum(dt), median(dt))
            end

            # ----------------------------------------------------------
            # Apply temperature step
            # ---------------------------------------------------------- 
            dtmp[:] .= 0.0

            # Use two-step scheme
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

            # Limit change as appropriate
            clamp!(dtmp, -1.0 * dtmp_clip, dtmp_clip)

            # Take step
            atmos.time += dt
            atmos.tmp += dtmp

            # ----------------------------------------------------------
            # (Optionally) Do convective adjustment and condensation
            # ---------------------------------------------------------- 
            # Dry convective adjustment
            adj_changed = 0
            if dryadj_steps > 0 && dry_convect && !use_mlt && start_con

                # do adjustment steps
                tmp_tnd = atmosphere.adjust_dry(atmos, dryadj_steps)
                clamp!(tmp_tnd, -1.0 * dtmp_clip, dtmp_clip)
                
                # check which levels were changed
                for i in 1:atmos.nlev_c
                    if abs(tmp_tnd[i]) > 0.01
                        adj_changed += 1
                        atmos.mask_c[i]  = atmos.mask_decay
                    end
                end 

                atmos.tmp +=  tmp_tnd
            end

            # Apply condensation
            if is_condense && start_con

                lvl_condensing = atmosphere.apply_vlcc!(atmos, condensate)
                
                # check which levels were changed
                for i in 1:atmos.nlev_c
                    if lvl_condensing[i]
                        adj_changed += 1
                        atmos.mask_p[i]  = atmos.mask_decay
                    end
                end 

            end

            # ----------------------------------------------------------
            # (Optionally) Smooth temperature profile for stability
            # ---------------------------------------------------------- 

            if smooth && (smooth_width > 2)
                atmosphere.smooth_centres!(atmos, smooth_width)
            end

            # ----------------------------------------------------------
            # Set remaining temperature values 
            # ---------------------------------------------------------- 
            # Temperature domain
            clamp!(atmos.tmp,  atmos.tmp_floor, atmos.tmp_ceiling)

            # Set cell-edge values (particularly important for handling conductive skin)
            if (surf_state == 2) && (step <= 20)
                atmosphere.set_tmpl_from_tmp!(atmos, 1, limit_change=true)  # don't allow tmpl to drift during accelerated phase
            else 
                atmosphere.set_tmpl_from_tmp!(atmos, surf_state, limit_change=true)
            end

            # Update tstar?
            if update_tstar 
                atmos.tstar = tmpl[end]
            end

            # Temperature domain
            clamp!(atmos.tmpl, atmos.tmp_floor, atmos.tmp_ceiling)
            
            # ----------------------------------------------------------
            # Calculate solver statistics
            # ---------------------------------------------------------- 
            # Calculate a measure of the relative rate of change in temperature.
            # Measured only at cell-centres because they are not set by interpolation.
            # 90th percentile is used to exclude 'difficult' layers, which can 
            # can be reasonably neglected if the rest of the column is fine.
            if step > 1
                drel_dt_prev = drel_dt
                drel_dt = quantile(abs.(  ((atmos.tmp[2:end-1] .- hist_tmp[end,2:end-1])./hist_tmp[end,2:end-1])./dt[2:end-1]  ), 0.95)
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

            F_TOA_tot = atmos.flux_tot[1]
            F_BOA_tot = atmos.flux_tot[end-1]
            F_loss    = abs( F_TOA_tot-F_BOA_tot )

            F_losspct = 0.0
            for i in 1:atmos.nlev_c
                # skip convective layers if adjustment is being used
                # skip condensing layers always
                skipcheck = false
                for j in -3:3  # neighbours
                    if (j+i <= atmos.nlev_c) && (j+i >= 1)  # don't go out of range
                        skipcheck = skipcheck || ( (atmos.mask_c[i+j] > 0) && !use_mlt)|| (atmos.mask_p[i+j] > 0)
                    end
                end
                # if this layer (and its neighbours) are valid for this criterion, check flux loss
                if !skipcheck
                    F_losspct = max(F_losspct, abs( (atmos.flux_n[i]-atmos.flux_n[i+1])/atmos.flux_n[i]))
                end 
            end 
            F_losspct *= 100.0
            
            # Calculate the 'typical' heating rate magnitude
            H_stat    = quantile(abs.(atmos.heating_rate[:]), 0.95)

            # --------------------------------------
            # Print debug info (some are disabled depending on the convection scheme)
            # -------------------------------------- 
            if verbose && (mod(step,modprint) == 0)
                @printf("    dtmp_comp   = %+.3f K      \n", dtmp_comp)
                @printf("    dtmp/tmp/dt = %+.3f day-1  \n", drel_dt)
                @printf("    drel_F_TOA  = %+.4f %%     \n", F_TOA_rel)
                @printf("    F_locloss   = %+.4f %%     \n", F_losspct)

                @printf("    F_rad^TOA   = %+.2e W m-2  \n", F_TOA_rad)
                @printf("    F_rad^BOA   = %+.2e W m-2  \n", F_BOA_rad)
                @printf("    F_rad^OLR   = %+.2e W m-2  \n", F_OLR_rad)
                @printf("    F_tot^loss  = %+.2f W m-2  \n", F_loss)
                @printf("    HR_typical  = %+.4f K day-1\n", H_stat)

                if ( (dryadj_steps > 0) && !use_mlt) || is_condense
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
                plotting.plot_solver(atmos, @sprintf("%s/solver.png", atmos.OUT_DIR), hist_tmpl=hist_tmpl, incl_magma=plt_magma, step=step)
                cp(@sprintf("%s/solver.png", atmos.OUT_DIR), @sprintf("%s/zzframe_%04d.png", atmos.OUT_DIR, step))

                plotting.plot_fluxes(atmos, @sprintf("%s/fluxes.png", atmos.OUT_DIR))
                cp(@sprintf("%s/fluxes.png", atmos.OUT_DIR), @sprintf("%s/zyframe_%04d.png", atmos.OUT_DIR, step))
            end 

            # --------------------------------------
            # Convergence check
            # -------------------------------------- 
            # Requires that:
            # - minimal rate (dT/dt) of temperature change for two iters
            # - minimal change to net radiative flux at TOA for two iters
            # - minimal flux loss between neighbouring layers
            # - solver is not being accelerated in any sense
            flag_prev = flag_this
            flag_this = (drel_dt < drel_dt_conv) && ( F_TOA_rel < drel_F_conv) && (F_losspct < F_losspct_conv) 
            success   = flag_this && flag_prev && !smooth && !accel && !stopaccel && (step > min_steps) && (start_con || (!dry_convect && !do_condense))
            
            # --------------------------------------
            # Sleep in order to capture keyboard interrupt
            # -------------------------------------- 
            sleep(1e-5)
            atexit(exit)

        end # end main loop

        # Print information about the final state
        if !success
            @printf("    stopping atmosphere iterations before convergence \n")
        else
            @printf("    convergence criteria met (%d iterations) \n", step)
        end
        @printf("    dtmp_comp   = %.3f K      \n", dtmp_comp)
        @printf("    dtmp/tmp/dt = %.3f day-1  \n", drel_dt)
        @printf("    dF/F (TOA)  = %.4f %%     \n", F_TOA_rel)
        @printf("    local loss  = %.2f %%     \n", F_losspct)
        @printf("\n")

        @printf("    endpoint fluxes \n")
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

        rm(joinpath(atmos.OUT_DIR,"fluxes.png"), force=true)
        rm(joinpath(atmos.OUT_DIR,"solver.png"), force=true)
        return nothing

    end # end solve_time


end 
