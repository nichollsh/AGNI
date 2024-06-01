# Contains code for obtaining energy balance (simplistic euler method)

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module solver_tstep

    using Printf
    using LoggingExtras
    using Statistics

    import ..atmosphere 
    import ..energy
    import ..phys
    import ..plotting

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
    - `sol_type::Int=1`                 bottom layer temperature, 0: free | 1: fixed | 2: conductive skin
    - `use_physical_dt::Bool=false`     use a single time-step across the entire column (for physical time-evolution)
    - `convect::Bool=true`              enable convection
    - `condensates::Array=[]`           condensates to model (if empty, no condensates are modelled)
    - `chem_type::Int=0`                chemistry type (see wiki)
    - `use_mlt::Bool=true`              using mixing length theory to represent convection (otherwise use adjustment)
    - `sens_heat::Bool=false`           include sensible heating 
    - `conduct::Bool=true`              include conduction
    - `modprop::Int=1`                  frequency at which to update thermodynamic properties (0 => never)
    - `verbose::Bool=false`             verbose output
    - `modplot::Int=0`                  plot frequency (0 => no plots)
    - `save_frames::Bool=true`          save plot frames
    - `accel::Bool=true`                enable accelerated fast period at the start 
    - `adams::Bool=true`                use Adams-Bashforth integrator
    - `dt_max::Float64=1000.0            maximum time-step outside of the accelerated phase
    - `max_steps::Int=1000`             maximum number of solver steps
    - `min_steps::Int=100`              minimum number of solver steps
    - `max_runtime::Float64=400.0`      maximum runtime in wall-clock seconds
    - `step_rtol::Float64=1.0e-4`       step size: relative change in per-level temperature [dimensionless]
    - `step_atol::Float64=1.0e-2`       step size: absolute change in per-level temperature [K]
    - `conv_rtol::Float64=1.0e-4`       convergence: relative tolerance on per-level flux loss [dimensionless]
    - `conv_atol::Float64=1.0e-1`       convergence: absolute tolerance on per-level flux loss [W m-2]
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            sol_type::Int=1, use_physical_dt::Bool=false,
                            convect::Bool=true, condensates::Array{String,1}=String[], chem_type::Int=0,

                            use_mlt::Bool=true,
                            sens_heat::Bool=false, conduct::Bool=true, modprop::Int=1, 
                            verbose::Bool=true, modplot::Int=0, save_frames::Bool=true,
                            accel::Bool=true, adams::Bool=true, dt_max::Float64=1000.0, 
                            max_steps::Int=1000, min_steps::Int=100, max_runtime::Float64=400.0,
                            step_rtol::Float64=1.0e-4, step_atol::Float64=1.0e-2,
                            conv_rtol::Float64=1.0e-4, conv_atol::Float64=1.0e-1
                            )::Bool


        # Start timer 
        wct_start::Float64 = time()

        # Run parameters
        dtmp_accel::Float64 = 15.0   # Change in temperature below which to stop model acceleration (needs to be turned off at small dtmp)
        smooth_stp::Int     = 120    # Number of steps for which to apply smoothing
        do_smooth::Bool     = true   # Is smoothing allowed at ever?
        wait_con::Int       = 30     # Introduce convection after this many steps, if ^^ is not already true
        modprint::Int       = 50     # Frequency to print when verbose==false
        len_hist::Int       = 10     # Number of previous states to store
        H_large::Float64    = 1.0e5  # A characteristic large heating rate [K/day]

        do_condense::Bool   = false  # Allow condensation ever? (overwritten according to condensate)

        if length(condensates) > 0 
            for c in condensates
                if condensate in atmos.gas_all_names
                    do_condense = true 
                else 
                    error("Invalid condensate ('$c')")
                end 
            end
        end 

        if verbose
            modprint=10
        end

        # Validate inputs
        wait_con  = max(wait_con, 1)
        min_steps = max(min_steps,wait_con+50)
        min_steps = max(min_steps,smooth_stp+50)
        max_steps = max(min_steps+50, max_steps)
        len_hist  = max(min(len_hist, max_steps-1), 5)
        modprint  = max(modprint, 1)

        if (len_hist < 4) && adams
            error("The value of len_hist is too small ($len_hist < 4)")
        end

        # Tracking
        adj_changed::Int   = 0           # Number of convectively adjusted levels
        F_loss::Float64    = 1.0e99      # Total flux loss (TOA vs BOA-1)
        H_stat::Float64    = 1.0e99      # Heating rate statistic
        F_TOA_rad::Float64 = 1.0e99      # Net upward TOA radiative flux
        F_BOA_rad::Float64 = 0.0         # Net upward BOA radiative flux
        F_OLR_rad::Float64 = 0.0         # Longwave upward TOA radiative flux
        F_TOA_tot::Float64 = 0.0         # Total TOA flux
        F_BOA_tot::Float64 = 0.0         # Total BOA flux
        dtmp_comp::Float64 = Inf         # Temperature change comparison
        drel_dt::Float64   = Inf         # Rate of relative temperature change
        drel_dt_prev::Float64 = Inf      # Previous ^
        start_con::Bool       = false       # Convection has been started at any point
        runtime::Float64  = 0.0

        oscil       = falses(atmos.nlev_c)          # layers which are oscillating
        dt          = ones(Float64,  atmos.nlev_c)  # time-step [days]
        dtmp        = zeros(Float64, atmos.nlev_c)  # temperature step [K]

        # Convergence 
        F_cri_worst::Float64 = Inf
        F_tol_worst::Float64 = 1.0
        F_rto_worst::Float64 = Inf
        F_rto::Float64       = F_rto_worst
        F_cri::Float64       = F_cri_worst
        F_tol::Float64       = F_tol_worst

        # Variables
        is_condense::Bool   = false       # Is condensation currently enabled?
        is_smooth::Bool     = true        # Currently smoothing?
        step::Int           = 0           # Current step number
        success::Bool       = false       # Convergence criteria met
        flag_prev::Bool     = false       # Previous iteration is meeting convergence
        flag_this::Bool     = false       # Current iteration is meeting convergence
        stopaccel::Bool     = false       # stop accelerated period this iter
        this_rtol::Float64  = step_rtol   # rtol used this step
        step_stopaccel::Int = 99999       # step at which accelerated period was stopped

        info_str::String = ""

        if !accel 
            step_stopaccel = 1
        end

        smooth_width::Int = 5

        # Store previous n atmosphere states
        hist_tmp  = zeros(Float64, (len_hist, atmos.nlev_c))   # temperature (cc)
        hist_tmpl = zeros(Float64, (len_hist, atmos.nlev_l))   # temperature (ce)
        hist_hr   = zeros(Float64, (len_hist, atmos.nlev_c))   # heating rates (cc)

        # Plots
        path_prf::String = @sprintf("%s/solver_prf.png", atmos.OUT_DIR)
        path_flx::String = @sprintf("%s/solver_flx.png", atmos.OUT_DIR)
        function plot_step(i::Int, t::Float64)
            if save_frames
                title_prf = @sprintf("i = %d",i)
                title_flx = @sprintf("t = %.1f s",t)
            end
            plotting.plot_pt(atmos,     path_prf, incl_magma=(sol_type==2), condensates=condensates, title=title_prf)
            plotting.plot_fluxes(atmos, path_flx, incl_eff=(sol_type==3), incl_cdct=conduct, incl_phase=do_condense, title=title_flx)
            if save_frames
                cp(path_prf,@sprintf("%s/frames/%04d_prf.png",atmos.OUT_DIR,i))
                cp(path_flx,@sprintf("%s/frames/%04d_flx.png",atmos.OUT_DIR,i))
            end 
        end 
 

        # Plot initial state 
        if modplot > 0 
            plot_step(0,0.0)
        end 

        # Main loop
        while (!success) && (step < max_steps)
            step += 1
            info_str = ""

            # ----------------------------------------------------------
            # Check that solver is happy 
            # ----------------------------------------------------------
            if !(all(isfinite, atmos.tmp) && all(isfinite, atmos.tmpl))
                error("Temperature array contains NaNs")
            end 
            if any(<=(0), atmos.tmp) || any(<=(0), atmos.tmpl)
                error("Temperature array contains negative values")
            end

            # Check time 
            runtime = time()-wct_start
            if runtime > max_runtime 
                break
            end 

            # ----------------------------------------------------------
            # Prepare for new iteration
            # ----------------------------------------------------------

            atmos.mask_c[:] .-= 1.0
            atmos.mask_p[:] .-= 1.0

            if mod(step,modprint) == 0
                info_str *= @sprintf("    step %d ", step)
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
                    info_str *= @sprintf("(accel) ")
                end
                dtmp_clip    = 50.0
                dryadj_steps = 4
                is_condense  = do_condense
                dt_min_step  = 1.0e-1
                dt_max_step  = 120.0
                this_rtol    = step_rtol
                is_smooth    = do_smooth 
                smooth_width = max(5, floor(Int, 0.2*atmos.nlev_c)) # 20% of levels, and >=5

            else
                # Slow phase
                dtmp_clip    = 5.0
                dryadj_steps = 16
                is_condense  = do_condense
                dt_min_step  = 5.0e-5
                dt_max_step  = max(dt_max, 1.0e-2)
                this_rtol    = step_rtol
                
                if is_smooth && (step >= smooth_stp)
                    is_smooth = false
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

                if is_smooth
                    smooth_width = max(5, floor(Int, 0.08*atmos.nlev_c)) # 8% of levels, and >= 5
                    dt_max_step *= 10.0
                    this_rtol *= 12.0
                end
            end

            # Introduce convection and condensation schemes
            if !start_con && (step >= wait_con) && (convect || do_condense)
                start_con = true 
                info_str *= @sprintf("(intro convect/condense) ")
            end

            if (mod(step,modprint) == 0) && is_smooth
                info_str *= @sprintf("(smooth) ")
            end

            # ----------------------------------------------------------
            # Recalculate thermodynamic properties at each layer (+ height & gravity)
            # ---------------------------------------------------------- 
            if (modprop > 0) && ( (mod(step, modprop) == 0) || (step < 10))
                if mod(step,modprint) == 0
                    info_str *= @sprintf("(props) ")
                end
                atmosphere.calc_layer_props!(atmos)
            end 

            if mod(step,modprint) == 0
                info_str *= ""
                @info info_str
            end

            # ----------------------------------------------------------
            # Get fluxes 
            # ---------------------------------------------------------- 
            energy.calc_fluxes!(atmos, do_condense, use_mlt && convect && start_con, sens_heat, conduct, condensates=condensates)

            # ----------------------------------------------------------
            # Calculate heating rates
            # ---------------------------------------------------------- 
            energy.calc_hrates!(atmos)

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
                        if (abs(atmos.tmp[i] - hist_tmp[end,i]) < (this_rtol*abs(hist_tmp[end,i]) + step_atol) ) 
                            dt[i] *= 1.03
                        else
                            dt[i] /= 1.02
                        end
                    end
                end
            end 

            # Enforce time-stepping limits 
            dt[end] = min(dt[end], dt_max_step*0.95)
            clamp!(dt, dt_min_step, dt_max_step)

            # Enforce a constant dt across column 
            if use_physical_dt
                fill!(dt, minimum(dt))
            end 

            if verbose && (mod(step,modprint) == 0)
                @info @sprintf("    dt_max,med  = %.5f, %.5f days ", maximum(dt), median(dt))
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
            if dryadj_steps > 0 && convect && !use_mlt && start_con

                # do adjustment steps
                tmp_tnd = energy.adjust_dry(atmos, dryadj_steps)
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

                for c in condensates
                    lvl_condensing = atmosphere.apply_vlcc!(atmos, c)
                    
                    # check which levels were changed
                    for i in 1:atmos.nlev_c
                        if lvl_condensing[i]
                            adj_changed += 1
                            atmos.mask_p[i] = atmos.mask_decay
                        end
                    end 
                end
            end

            # ----------------------------------------------------------
            # (Optionally) Smooth temperature profile for stability
            # ---------------------------------------------------------- 
            if is_smooth
                atmosphere.smooth_centres!(atmos, smooth_width)
            end

            # ----------------------------------------------------------
            # Set remaining temperature values 
            # ---------------------------------------------------------- 
            # Temperature domain
            clamp!(atmos.tmp,  atmos.tmp_floor, atmos.tmp_ceiling)

            # Set cell-edge values (particularly important for handling conductive skin)
            atmosphere.set_tmpl_from_tmp!(atmos, limit_change=true)

            # Set bottom edge temperature 
            if (sol_type < 0) || (sol_type > 2)
                # Error cases
                error("Invalid surface state ($sol_type) for timestepping solver")

            elseif (sol_type == 0)
                # Extrapolate (log-linear)
                atmos.tmpl[end] = atmos.tmp[end] +   (atmos.tmp[end]-atmos.tmp[end-1])/(log(atmos.p[end]/atmos.p[end-1]))   * log(atmos.pl[end]/atmos.p[end])

            # elseif (sol_type==1) 
            #   do nothing in this case
                
            elseif (sol_type == 2) && (step > 15)
                # Conductive skin
                atmos.tmpl[end] = atmos.tmp_magma - atmos.flux_tot[1] * atmos.skin_d / atmos.skin_k
                atmos.tmpl[end] = max(atmos.tmp_floor, atmos.tmpl[end])
                atmos.tmp_surf = atmos.tmpl[end]

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
                    tmp_comp_1::Float64 = 0.5 * ( atmos.tmp[i]      + hist_tmp[end  ,i] )
                    tmp_comp_2::Float64 = 0.5 * ( hist_tmp[end-1,i] + hist_tmp[end-2,i] )
                    dtmp_comp = max(dtmp_comp, abs(tmp_comp_1 - tmp_comp_2))
                end 
            end 

            # Calculate the (change in) flux balance
            F_TOA_rad = atmos.flux_n[1]
            F_BOA_rad = atmos.flux_n[end-1]
            F_OLR_rad = atmos.flux_u_lw[1]

            F_TOA_tot = atmos.flux_tot[1]
            F_BOA_tot = atmos.flux_tot[end-1]
            F_loss    = abs( F_TOA_tot-F_BOA_tot )  # flux lost across whole column
            F_rto_worst = 0.0
            for i in 1:atmos.nlev_c
                # skip convective layers if adjustment is being used
                # skip condensing layers always
                skipcheck::Bool = false
                for j in -2:3  # neighbours
                    if (j+i <= atmos.nlev_c) && (j+i >= 1)  # don't go out of range
                        skipcheck = skipcheck || ( (atmos.mask_c[i+j] > 0) && !use_mlt)|| (atmos.mask_p[i+j] > 0)
                    end
                end

                # if this layer (and its neighbours) are valid for this criterion, check flux loss
                if !skipcheck
                    F_tol = abs(atmos.flux_tot[i+1]) * conv_rtol + conv_atol  # required loss for convergence
                    F_cri = abs(atmos.flux_tot[i] - atmos.flux_tot[i+1])
                    F_rto = F_cri/F_tol 
                    if F_rto > F_rto_worst
                        F_rto_worst = F_rto
                        F_cri_worst = F_cri
                        F_tol_worst = F_tol
                    end 
                end 
            end 
            
            # Calculate the 'typical' heating rate magnitude
            H_stat    = quantile(abs.(atmos.heating_rate[:]), 0.95)

            # --------------------------------------
            # Print debug info (some are disabled depending on the convection scheme)
            # -------------------------------------- 
            if verbose && (mod(step,modprint) == 0)
                @info @sprintf("    dtmp_comp   = %+.3f K      ", dtmp_comp)
                @info @sprintf("    dtmp/tmp/dt = %+.3f day-1  ", drel_dt)
                @info @sprintf("    F_cri_worst = %+.2e W m-2  ", F_cri_worst)
                @info @sprintf("    F_tol_worst = %+.2e W m-2  ", F_tol_worst)

                @info @sprintf("    F_rad^TOA   = %+.2e W m-2  ", F_TOA_rad)
                @info @sprintf("    F_rad^BOA   = %+.2e W m-2  ", F_BOA_rad)
                @info @sprintf("    F_rad^OLR   = %+.2e W m-2  ", F_OLR_rad)
                @info @sprintf("    F_tot^loss  = %+.2f W m-2  ", F_loss)
                @info @sprintf("    HR_typical  = %+.4f K day-1", H_stat)

                if ( (dryadj_steps > 0) && !use_mlt) || is_condense
                    @info @sprintf("    count_adj   = %d layers   ", adj_changed)
                end
                @info @sprintf(" ")
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
            if (modplot > 0) && (mod(step, modplot) == 0)
                plot_step(step, runtime)
            end 

            # --------------------------------------
            # Convergence check
            # -------------------------------------- 
            # Requires that:
            # - flux loss is within the rtol,atol tolerances for two consecutive steps
            # - if enabled, convection/condensation are active
            # - solver is not being accelerated or smoothed
            flag_prev = flag_this
            flag_this = (F_rto_worst <= 1.0)
            success   = flag_this && flag_prev && !is_smooth && !accel && !stopaccel && (step > min_steps) && (start_con || (!convect && !do_condense))
            
            # --------------------------------------
            # Sleep in order to capture keyboard interrupt
            # -------------------------------------- 
            sleep(1e-5)
            atexit(exit)

        end # end main loop

        rm(path_prf, force=true)
        rm(path_flx, force=true)

        # Print information about the final state
        atmos.is_solved = true
        if !success
            @warn @sprintf("    stopping atmosphere iterations before convergence ")
            atmos.is_converged = false
        else
            @info @sprintf("    convergence criteria met (%d iterations) ", step)
            atmos.is_converged = true
        end
        @info @sprintf("    dtmp_comp   = %.3f K      ", dtmp_comp)
        @info @sprintf("    dtmp/tmp/dt = %.3f day-1  ", drel_dt)
        @info @sprintf("    F_cri_worst = %+.2e W m-2  ", F_cri_worst)
        @info @sprintf("    F_tol_worst = %+.2e W m-2  ", F_tol_worst)
        if use_physical_dt
            @info @sprintf("    total_time  = %+.2e years  ", atmos.time[1]/365.25)
        else 
            @info @sprintf("    max_time    = %+.2e years  ", maximum(atmos.time)/365.25)
        end 

        @info @sprintf("    endpoint fluxes ")
        @info @sprintf("    rad_OLR   = %.2e W m-2         ", F_OLR_rad)
        @info @sprintf("    rad_TOA   = %.2e W m-2         ", F_TOA_rad)
        @info @sprintf("    rad_BOA   = %.2e W m-2         ", F_BOA_rad)
        @info @sprintf("    tot_BOA   = %.2e W m-2         ", F_BOA_tot)
        if (sol_type == 2)
            F_skin = atmos.skin_k / atmos.skin_d * (atmos.tmp_magma - atmos.tmp_surf)
            @info @sprintf("    cond_skin = %.2e W m-2         ", F_skin)
        end
    
        # Warn user if there's a sign difference in TOA vs BOA fluxes
        # because this shouldn't be the case
        if F_TOA_tot*F_BOA_tot < 0
            @warn @sprintf("TOA and BOA total fluxes have different signs")
        end

        return atmos.is_converged

    end # end solve_time


end 
