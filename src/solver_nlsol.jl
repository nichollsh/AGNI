# Contains code for obtaining energy balance (CVODE method)

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module solver_nlsol

    using Printf
    using Statistics
    using Revise
    using LinearAlgebra

    import atmosphere 
    import phys
    import plotting

    
    """
    **Obtain radiative-convective equilibrium using a matrix method.**

    Solves the non-linear system of equations defined by the flux field
    divergence, minimising flux loss across a cell by iterating the temperature
    profile.

    Not compatible with convective adjustment; MLT must be used.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `surf_state::Int=1`               bottom layer temperature, 0: free | 1: fixed | 2: skin | 3: Tint
    - `condensate::String=""`           condensate to model (if empty, no condensates are modelled)
    - `dry_convect::Bool=true`          enable dry convection
    - `sens_heat::Bool=false`           include sensible heating 
    - `max_steps::Int=2000`             maximum number of solver steps
    - `max_runtime::Float64=600.0`      maximum runtime in wall-clock seconds
    - `fdw::Float64=1.0e-4`             relative width of the "difference" in the finite-difference calculations
    - `use_cendiff::Bool=false`         use central difference for calculating jacobian? If false, use forward difference
    - `method::Int=2`                   numerical method (1: Brute-Force, 2: Newton-Raphson, 3: Gauss-Newton, 4: Levenberg-Marquardt)
    - `linesearch::Bool=true`           use a simple linesearch algorithm to determine the best step size
    - `modplot::Int=0`                  iteration frequency at which to make plots
    - `stabilise_mlt::Bool=true`        stabilise convection by introducing it gradually
    - `step_rtol::Float64=1.0e-2`       step size: relative change in residuals [dimensionless]
    - `step_atol::Float64=1.0e-1`       step size: absolute change in residuals [W m-2]
    - `step_fact::Float64=8.0e4`        step size: scale factor for maximum per-level temperature step [K^2]
    - `conv_atol::Float64=1.0e-2`       convergence: absolute tolerance on per-level flux deviation [W m-2]
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            surf_state::Int=1, condensate::String="",
                            dry_convect::Bool=true, sens_heat::Bool=false,
                            max_steps::Int=2000, max_runtime::Float64=600.0,
                            fdw::Float64=1.0e-4, use_cendiff::Bool=false, 
                            method::Int=2, linesearch::Bool=true,
                            modplot::Int=1, stabilise_mlt::Bool=true,
                            step_rtol::Float64=1.0e-2, step_atol::Float64=1.0e-6, step_fact::Float64=8e4,
                            conv_atol::Float64=1.0e-2
                            )

        # Validate condensation case
        do_condense::Bool  = false 
        i_gas::Int         = -1
        if condensate != "" 
            if condensate in atmos.gases
                do_condense = true 
                i_gas = findfirst(==(condensate), atmos.gases)
            else 
                error("Invalid condensate '$condensate'")
            end 
        end 

        # Validate surf_state
        if (surf_state < 0) || (surf_state > 3)
            error("Invalid surface state ($surf_state)")
        end

        # Start timer 
        wct_start::Float64 = time()

        # Plot paths 
        path_prf::String = @sprintf("%s/solver_prf.png", atmos.OUT_DIR)
        path_flx::String = @sprintf("%s/solver_flx.png", atmos.OUT_DIR)

        # Dimensionality
        arr_len::Int = atmos.nlev_c 
        if (surf_state >= 2)  # states 2 and 3
            arr_len += 1
        end

        # Work arrays 
        prate::Array{Float64,1} = zeros(Float64, atmos.nlev_c)  # condensation production rate [kg /m3 /s]
        rf::Array{Float64,1}    = zeros(Float64, arr_len)  # Forward difference
        rb::Array{Float64,1}    = zeros(Float64, arr_len)  # Backward difference
        x_s::Array{Float64,1}   = zeros(Float64, arr_len)  # Perturbed row, for jacobian
        fd_s::Float64           = 0.0                      # Row perturbation amount

        # Convective flux scale factor 
        convect_sf::Float64 = 1.0e-5

        # Calculate the (remaining) temperatures  
        function _set_tmps!(_x::Array)
            # Read new guess
            for i in 1:atmos.nlev_c
                atmos.tmp[i] = _x[i]
            end
            clamp!(atmos.tmp, atmos.tmp_floor+1.0, atmos.tmp_ceiling-1.0)

            # Interpolate temperature to cell-edge values (not incl. bottommost value)
            atmosphere.set_tmpl_from_tmp!(atmos)

            # Set bottom edge temperature 
            if (surf_state != 1) 
                # For state=1, tmpl[end] is held constant

                # Extrapolate (log-linear)
                grad_dt = atmos.tmp[end]-atmos.tmp[end-1]
                grad_dp = log(atmos.p[end]/atmos.p[end-1])
                atmos.tmpl[end] = atmos.tmp[end] + grad_dt/grad_dp * log(atmos.pl[end]/atmos.p[end])

                if (surf_state >= 2)  # states 2 and 3
                    atmos.tstar = _x[end]  # Surface brightness temperature
                end
            end 
            
            return nothing
        end # end set_tmps

        # Objective function to solve for
        function _fev!(x::Array,resid::Array)

            # Reset values
            atmos.mask_c[:] .= 0
            atmos.mask_p[:] .= 0
           
            # Set temperatures 
            _set_tmps!(x)

            # Reset fluxes
            atmos.flux_tot[:] .= 0.0

            # +Radiation
            atmosphere.radtrans!(atmos, true)
            atmosphere.radtrans!(atmos, false)
            atmos.flux_tot += atmos.flux_n

            # +Dry convection
            if dry_convect
                # Calc flux
                atmosphere.mlt!(atmos)

                # Stabilise?
                atmos.flux_c *= convect_sf

                # Add to total flux
                atmos.flux_tot += atmos.flux_c
            end

            # +Turbulence
            if sens_heat
                atmosphere.sensible!(atmos)
                atmos.flux_tot[end] += atmos.flux_sens
            end

            # Calculate residuals subject to the boundary condition
            if (surf_state == 0) || (surf_state == 1)
                # Conserve fluxes with constant tstar
                resid[1:end] .= atmos.flux_tot[2:end] .- atmos.flux_tot[1:end-1] 
            elseif (surf_state == 2)
                # Conductive boundary layer
                resid[1:end-1] = atmos.flux_tot[2:end] - atmos.flux_tot[1:end-1] 
                resid[end] = atmos.flux_tot[1] - (atmos.tmp_magma - atmos.tmpl[end]) * atmos.skin_k / atmos.skin_d
            elseif (surf_state == 3)
                # Fluxes equal to sigma*Tint^4
                resid[1:end] .= atmos.flux_tot[1:end] .- atmos.flux_int
            end

            # +Condensation
            if do_condense

                a = 1.0
                g = 1.0
                f = 0.0

                for i in 1:atmos.nlev_c

                    pp = atmos.layer_x[i,i_gas] * atmos.p[i]
                    if pp < 1.0e-10 
                        continue
                    end

                    # Layer is condensing if T < T_dew
                    Tsat = phys.calc_Tdew(condensate, pp )
                    g = 1.0 - exp(a * (Tsat - atmos.tmp[i]))
                    f = resid[i]
                    resid[i] = f * g
                    if atmos.tmp[i] < Tsat
                        println("Condensing at level $i")
                        atmos.mask_p[i] = atmos.mask_decay 
                        prate[i]        = -1.0 * f / ( phys.lookup_safe("l_vap",condensate) * (atmos.zl[i] - atmos.zl[i+1])) * 86.4 # g cm-3 day-1
                        atmos.re[i]     = 1.0e-5  # 10 micron droplets
                        atmos.lwm[i]    = 0.8     # 80% of the saturated vapor turns into cloud
                        atmos.clfr[i]   = 1.0     # The cloud takes over the entire cell
                    else 
                        prate[i] = 0.0
                        atmos.re[i]    = 0.0
                        atmos.lwm[i]   = 0.0
                        atmos.clfr[i]  = 0.0
                    end
                end
            end 

            # Check that residuals are real numbers
            if !all(isfinite, resid)
                display(resid)
                error("Residual array contains NaNs and/or Infs")
            end

            return nothing
        end  # end fev

        # Calculate the jacobian and residuals at x using a central-difference method
        function _calc_jac_res_cendiff!(x::Array, jacob::Array, resid::Array)

            # Evalulate residuals at x
            _fev!(x, resid)

            for i in 1:arr_len  # for each x

                # Calculate perturbation
                fd_s = x[i] * fdw

                # Forward
                x_s[:] .= x[:]
                x_s[i] += fd_s * 0.5
                _fev!(x_s, rf)

                # Backward
                x_s[i] -= fd_s    # Only need to modify this level; rest were set during the forward phase
                _fev!(x_s, rb)

                # Set jacobian
                for j in 1:arr_len  # for each r
                    jacob[j,i] = (rf[j] - rb[j]) / fd_s
                end 
            end 

            return nothing
        end # end jr_cd

        # Calculate the jacobian and residuals at x using a forward-difference method
        function _calc_jac_res_fordiff!(x::Array, jacob::Array, resid::Array)

            # Evalulate residuals at x
            _fev!(x, resid)

            for i in 1:arr_len  # for each x

                # Calculate perturbation
                fd_s = x[i] * fdw

                # Forward
                x_s[:] .= x[:]
                x_s[i] += fd_s
                _fev!(x_s, rf)

                # Set jacobian
                for j in 1:arr_len  # for each r
                    jacob[j,i] =  (rf[j] - resid[j]) / fd_s
                end 
            end 

            return nothing
        end # end jr_fd

        # Cost function 
        function _cost(_r::Array)
            return norm(_r)
        end 

        
        # ----------------------------------------------------------
        # Setup initial guess
        # ---------------------------------------------------------- 
        if method == 1
            println("Begin Brute-Force iterations")
        elseif method == 2
            println("Begin Newton-Raphson iterations") 
        elseif method == 3
            println("Begin Gauss-Newton iterations") 
        elseif method == 4
            println("Begin Levenberg-Marquardt iterations") 
        else 
            error("Invalid method choice ($method)")
        end

        @printf("    surf   = %d\n", surf_state)
        if (surf_state == 1)
            @printf("    tstar  = %.2f K\n", atmos.tstar)
        elseif (surf_state == 2)
            @printf("    skin_d = %.2f m\n",         atmos.skin_d)
            @printf("    skin_k = %.2f W K-1 m-1\n", atmos.skin_k)
        elseif (surf_state == 3)
            @printf("    tint   = %.2f K\n",     atmos.tint)
            @printf("    Fint   = %.2f W m-2\n", atmos.flux_int)
        end 
        

        # Allocate initial guess for the x array, as well as a,b arrays
        # Array storage structure:
        #   in the surf_state=2 case
        #       1:end-1 => cell centre temperatures 
        #       end     => bottom cell edge temperature
        #   other cases 
        #       1:end => cell centre temperatures
        x_ini    = zeros(Float64, arr_len) 
        for i in 1:atmos.nlev_c
            x_ini[i]    = clamp(atmos.tmp[i], atmos.tmp_floor , atmos.tmp_ceiling)
        end 
        if (surf_state >= 2)
            x_ini[end] = atmos.tmp[atmos.nlev_c] + 1.0
        end

        # ----------------------------------------------------------
        # Solver loop
        # ---------------------------------------------------------- 
        # Execution variables
        modprint::Int =         1       # Print frequency
        x_dif_clip::Float64 =   300.0   # Maximum allowed step size
        
        # Tracking variables
        step::Int =         0       # Step number
        code::Int =         -1      # Status code 
        lml::Float64 =      2.0     # Levenberg-Marquardt lambda parameter

        # Model statistics tracking
        r_med::Float64 =        9.0     # Median residual
        r_max::Float64 =        9.0     # Maximum residual (sign agnostic)
        x_med::Float64 =        0.0     # Median solution
        x_max::Float64 =        0.0     # Maximum solution (sign agnostic)
        iworst::Int =           0       # Level which is furthest from convergence
        dxmax::Float64 =        9.0     # Maximum change in solution array
        r_cur_2nm::Float64 =    0.01    # Two-norm of residuals 
        r_old_2nm::Float64 =    0.02    # Previous ^


        # Solver variables
        b::Array{Float64,2}      = zeros(Float64, (arr_len, arr_len))    # Approximate jacobian (i)
        x_cur::Array{Float64,1}  = zeros(Float64, arr_len)               # Current best solution (i)
        x_old::Array{Float64,1}  = zeros(Float64, arr_len)               # Previous best solution (i-1)
        x_dif::Array{Float64,1}  = zeros(Float64, arr_len)               # Change in x (i-1 to i)
        r_cur::Array{Float64,1}  = zeros(Float64, arr_len)               # Residuals (i)
        r_old::Array{Float64,1}  = zeros(Float64, arr_len)               # Residuals (i-1)
        r_tst::Array{Float64,1}  = zeros(Float64, arr_len)               # Test for rejection residuals
        dtd::Array{Float64,2}    = zeros(Float64, (arr_len,arr_len))     # Damping matrix for LM method
        c_cur::Float64 = Inf        # current cost (i)
        c_old::Float64 = Inf        # old cost (i-1)
        x_dif_clip_step::Float64 = Inf # maximum step size (IN THIS STEP)

        # Linesearch parameters
        ls_scale::Float64       = 1.0     # best scale factor for linesearch 
        ls_test_cost::Float64   = Inf 
        ls_best_cost::Float64   = Inf 
        ls_best_scale::Float64  = 1.0 

        # Final setup
        x_cur[:] .= x_ini[:]
        for di in 1:arr_len 
            dtd[di,di] = 1.0
        end 
        fill!(r_cur, 1.0e99)
        fill!(r_old, 1.0e98)

        # Stabilise convection?
        stabilise_mlt = stabilise_mlt && dry_convect
        if !stabilise_mlt
            convect_sf = 1.0
        end

        @printf("    step  resid_med  resid_2nm  flux_OLR   xvals_med  xvals_max  |dx|_max   flags\n")
        while true 

            # Update properties (cp, rho, etc.)
            if !all(isfinite, x_cur)
                display(x_cur)
                error("Solution array contains NaNs and/or Infs")
            end 
            _set_tmps!(x_cur)
            atmosphere.calc_layer_props!(atmos)

            # Check time 
            if time()-wct_start > max_runtime 
                code = 3
                break
            end 

            # Update step counter
            step += 1
            if step > max_steps
                code = 1 
                break 
            end 
            if mod(step,modprint) == 0 
                @printf("    %4d  ", step)
            end

            # Reset flags 
            #     Sc,R       = stabilise convection, reduced this step
            #     Cd, Fd     = finite differencing type 
            #     Nr, Gn, Lm = stepping algorithm
            #     Ls         = linesearch active
            stepflags::String = ""

            # Check convective stabilisation
            if stabilise_mlt 
                # We are stabilising
                stepflags *= "Sc"
                
                # Check if sf needs to be increased
                if c_cur < max(100.0*conv_atol, 10.0)
                    if convect_sf < 1.0 
                        # increase sf - reduce stabilisation by 10x
                        stepflags *= "r"
                        convect_sf = min(1.0, convect_sf*10.0)
                        println("convect_sf = $convect_sf")
                        
                        # done stabilising 
                        if convect_sf > 0.99
                            stabilise_mlt = false 
                        end
                    end 
                end
            else
                # No stabilisation at this point
                convect_sf = 1.0 
            end 

            # Evaluate jacobian and residuals
            r_old[:] .= r_cur[:]
            if use_cendiff || (step == 1)
                _calc_jac_res_cendiff!(x_cur, b, r_cur) 
                stepflags *= "Cd"
            else
                _calc_jac_res_fordiff!(x_cur, b, r_cur) 
                stepflags *= "Fd"
            end 

            # Check convergence
            c_old = c_cur
            c_cur = _cost(r_cur)
            if (c_cur < conv_atol) && !stabilise_mlt
                code = 0
                @printf("\n")
                break
            end

            # Check if jacobian is singular 
            if abs(det(b)) < 1.0e-90
                code = 2
                @printf("\n")
                break
            end 

            # Model step 
            x_dif_clip_step = x_dif_clip
            x_old[:] = x_cur[:]
            if method == 1
                # Brute-Force step
                stepflags *= "Bf"
                error("Brute-Force method is not yet implemented")

            elseif (method == 2)
                # Newton-Raphson step 
                x_dif = -b\r_cur
                stepflags *= "Nr"

            elseif method == 3
                # Gauss-Newton step 
                x_dif = -(b'*b) \ (b'*r_cur) 
                stepflags *= "Gn"

            elseif method == 4
                # Levenberg-Marquardt step
                #    Calculate damping parameter ("delayed gratification")
                if r_cur_2nm < r_old_2nm
                    lml /= 5.0
                else 
                    lml *= 1.5
                end

                #    Update our estimate of the solution
                x_dif = -(b'*b + lml * dtd) \ (b' * r_cur)
                stepflags *= "Lm"
            end

            # If not using brute force, limit step size
            if method != 1
                # Limit step size according to inverse temperature (impacts high temperatures)
                # for i in 1:atmos.nlev_c 
                #     x_dif[i] = sign(x_dif[i]) * min(abs(x_dif[i]), step_fact/x_old[i]^1.1)  # slower changes at large temperatures
                # end 

                # Linesearch 
                if linesearch && (step > 2)

                    # Reset
                    stepflags *= "Ls"
                    ls_best_cost = Inf
                    ls_best_scale = 1.0

                    for ls_scale in [0.3, 0.5, 0.8, 1.0] 
                        # try this scale factor 
                        x_cur[:] .= x_old[:] .+ (ls_scale .* x_dif[:])
                        _fev!(x_cur, r_tst)

                        # test improvement
                        ls_test_cost = _cost(r_tst)
                        if ls_test_cost < ls_best_cost 
                            ls_best_scale = ls_scale 
                            ls_best_cost = ls_test_cost 
                        end 
                    end 

                    # apply best linesearch scale 
                    x_dif[:] .*= ls_best_scale

                end # end linesearch 

                # Test new step 
                x_cur[:] .= x_old[:] .+ x_dif[:]
                _fev!(x_cur, r_tst)

                # If convection stabilisation is disabled, check if this step would make the residuals worse.
                # If so, limit the maximum step size to a small value, maintaining the same descent direction.
                # for i in 1:arr_len 
                #     if !stabilise_mlt && (abs(r_tst[i]-r_cur[i]) > abs(r_cur[i]) * step_rtol + step_atol) && ((i>atmos.nlev_c) || (atmos.mask_c[i] > 0))
                #         x_dif_clip_step = 1.0e-2
                #         break
                #     end 
                # end 
            end 

            # Limit step size globally
            if maximum(abs.(x_dif[:])) > x_dif_clip_step
                x_dif[:] .*= x_dif_clip_step ./ maximum(abs.(x_dif[:]))
            end

            # Take the step 
            x_cur[:] .= x_old[:] .+ x_dif[:]
            _fev!(x_cur, r_cur)

            # Model statistics 
            r_med   = median(r_cur)
            iworst  = argmax(abs.(r_cur))
            r_max   = r_cur[iworst]
            x_med   = median(x_cur)
            x_max   = x_cur[argmax(abs.(x_cur))]
            dxmax   = maximum(abs.(x_dif))
            r_old_2nm   = r_cur_2nm
            r_cur_2nm   = norm(r_cur)

            # Plot
            if modplot > 0
                if mod(step, modplot) == 0
                    plotting.plot_pt(atmos,     path_prf, incl_magma=(surf_state==2))
                    plotting.plot_fluxes(atmos, path_flx, incl_int=(surf_state==3))
                end 
            end 
                
            # Inform user
            if mod(step,modprint) == 0 
                @printf("%+.2e  %.3e  %.3e  %+.2e  %+.2e  %.3e  %-9s\n", r_med, r_cur_2nm, atmos.flux_u_lw[1], x_med, x_max, dxmax, stepflags)
            end

        end # end solver loop
        
        rm(path_prf, force=true)
        rm(path_flx, force=true)
        
        # ----------------------------------------------------------
        # Extract solution
        # ---------------------------------------------------------- 
        atmos.is_solved = true
        atmos.is_converged = false 
        if code == 0
            println("    success")
            atmos.is_converged = true
        elseif code == 1
            println("    failure (maximum iterations)")
        elseif code == 2
            println("    failure (singular jacobian)")
        elseif code == 3
            println("    failure (maximum time)")
        else 
            println("    failure (unhandled)")
        end
        println(" ")

        _fev!(x_cur, zeros(Float64, arr_len))
        atmosphere.calc_hrates!(atmos)

        # ----------------------------------------------------------
        # Print info
        # ---------------------------------------------------------- 
        loss = maximum(atmos.flux_tot) - minimum(atmos.flux_tot)
        loss_pct = loss/maximum(atmos.flux_tot)*100.0
        @printf("    summary \n")
        @printf("    outgoing LW flux   = %+.2e W m-2     \n", atmos.flux_u_lw[1])
        if (surf_state == 2)
            F_skin = atmos.skin_k / atmos.skin_d * (atmos.tmp_magma - atmos.tstar)
        @printf("    conduct. skin flux = %+.2e W m-2 \n", F_skin)
        end
        @printf("    total flux at TOA  = %+.2e W m-2     \n", atmos.flux_tot[1])
        @printf("    total flux at BOA  = %+.2e W m-2     \n", atmos.flux_tot[end])
        @printf("    column max. loss   = %+.2e W m-2  (%+.2e %%) \n", loss, loss_pct)
        @printf("    final cost value   = %+.2e W m-2     \n", c_cur)
        @printf("\n")

        return atmos.is_converged
    end # end solve_energy 
end 
