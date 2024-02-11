# Contains code for obtaining energy balance (CVODE method)

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module solver_nlsol

    # Include libraries
    include("../socrates/julia/src/SOCRATES.jl")

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
    - `max_steps::Int=200`              maximum number of solver steps
    - `atol::Float64=1.0e-2`            maximum residual at convergence
    - `fdw::Float64=5.0e-3`             relative width of the "difference" in the finite-difference calculations
    - `use_cendiff::Bool=false`         use central difference for calculating jacobian? If false, use forward difference
    - `method::Int=0`                   numerical method (0: Newton-Raphson, 1: Gauss-Newton, 2: Levenberg-Marquardt)
    - `calc_cf_end::Bool=true`          calculate contribution function?
    - `modplot::Int=0`                  iteration frequency at which to make plots
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            surf_state::Int=1, condensate::String="",
                            dry_convect::Bool=true, sens_heat::Bool=false,
                            max_steps::Int=200, atol::Float64=1.0e-2, 
                            fdw::Float64=1.0e-4, use_cendiff::Bool=false, method::Int=0,
                            calc_cf_end::Bool=true, modplot::Int=1
                            )

        # Validate condensation case
        do_condense  = false 
        i_gas        = -1
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

        # Dimensionality
        arr_len = atmos.nlev_c 
        if (surf_state >= 2)
            arr_len += 1
        end

        # Work arrays 
        calc_cf = false
        prate =   zeros(Float64, atmos.nlev_c)  # condensation production rate [kg /m3 /s]
        rf      = zeros(Float64, arr_len)  # Forward difference
        rb      = zeros(Float64, arr_len)  # Backward difference
        x_s     = zeros(Float64, arr_len)  # Perturbed row, for jacobian
        s       = 0.0                      # Row perturbation amount

        # Calculate the (remaining) temperatures  
        function _set_tmps!(_x::Array)
            # Read new guess
            for i in 1:atmos.nlev_c
                atmos.tmp[i] = _x[i]
            end
            clamp!(atmos.tmp, atmos.tmp_floor, atmos.tmp_ceiling)

            # Interpolate temperature to cell-edge values (not incl. bottommost value)
            atmosphere.set_tmpl_from_tmp!(atmos)

            # Set bottom edge temperature 
            if (surf_state != 1)
                # Extrapolate (log-linear)
                grad_dt = atmos.tmp[end]-atmos.tmp[end-1]
                grad_dp = log(atmos.p[end]/atmos.p[end-1])
                atmos.tmpl[end] = atmos.tmp[end] + grad_dt/grad_dp * log(atmos.pl[end]/atmos.p[end])

                if (surf_state >= 2) # Surface brightness temperature
                    atmos.tstar = _x[end]
                end
            end 
            
            return nothing
        end # end set_tmps

        # Objective function to solve for
        function fev!(x::Array,resid::Array)

            # Reset values
            atmos.mask_c[:] .= 0
            atmos.mask_p[:] .= 0
           
            # Set temperatures 
            _set_tmps!(x)

            # Reset fluxes
            atmos.flux_tot[:] .= 0.0

            # +Radiation
            atmosphere.radtrans!(atmos, true, calc_cf=calc_cf)
            atmosphere.radtrans!(atmos, false)
            atmos.flux_tot += atmos.flux_n

            # +Dry convection
            if dry_convect
                atmosphere.mlt!(atmos)
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
        function calc_jac_res_cendiff!(x::Array, jacob::Array, resid::Array)

            # Evalulate residuals at x
            fev!(x, resid)

            for i in 1:arr_len  # for each x

                # Calculate perturbation
                s = x[i] * fdw

                # Forward
                x_s[:] .= x[:]
                x_s[i] += s * 0.5
                fev!(x_s, rf)

                # Backward
                x_s[i] -= s    # Only need to modify this level; rest were set during the forward phase
                fev!(x_s, rb)

                # Set jacobian
                for j in 1:arr_len  # for each r
                    jacob[j,i] = (rf[j] - rb[j]) / s
                end 
            end 

            return nothing
        end # end jr_cd

        # Calculate the jacobian and residuals at x using a forward-difference method
        function calc_jac_res_fordiff!(x::Array, jacob::Array, resid::Array)

            # Evalulate residuals at x
            fev!(x, resid)

            for i in 1:arr_len  # for each x

                # Calculate perturbation
                s = x[i] * fdw

                # Forward
                x_s[:] .= x[:]
                x_s[i] += s
                fev!(x_s, rf)

                # Set jacobian
                for j in 1:arr_len  # for each r
                    jacob[j,i] =  (rf[j] - resid[j]) / s
                end 
            end 

            return nothing
        end # end jr_fd

        # Cost function 
        function cost(_r::Array)
            return maximum(abs.(_r))
            # return norm(_r)
        end 

        
        # ----------------------------------------------------------
        # Setup initial guess
        # ---------------------------------------------------------- 
        if method == 0
            println("Begin Newton-Raphson iterations") 
        elseif method == 1
            println("Begin Gauss-Newton iterations") 
        elseif method == 2
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
        modprint::Int =     1       # Print frequency

        # Tracking variables
        step::Int =         0       # Step number
        code::Int =         -1      # Status code 
        lml::Float64 =      20.0    # Levenberg-Marquardt lambda parameter

        # Model statistics tracking
        r_med::Float64 =        9.0     # Median residual
        r_max::Float64 =        9.0     # Maximum residual (sign agnostic)
        x_med::Float64 =        0.0     # Median solution
        x_max::Float64 =        0.0     # Maximum solution (sign agnostic)
        dx2nm::Float64 =        9.0     # Current value for 2-norm for the relative change in the solution array
        r_cur_2nm::Float64 =    0.01    # Two-norm of residuals 
        r_old_2nm::Float64 =    0.02    # Previous ^
        count_rej::Int =        0       # Number of rejected steps
        reject::Bool =          false   # step was rejected?

        # Work variables
        b::Array      = zeros(Float64, (arr_len, arr_len))    # Approximate jacobian (i)
        x_cur::Array  = zeros(Float64, arr_len)               # Current best solution (i)
        x_old::Array  = zeros(Float64, arr_len)               # Previous best solution (i-1)
        x_dif::Array  = zeros(Float64, arr_len)               # Change in x (i-1 to i)
        x_dla::Array  = zeros(Float64, arr_len)               # Change in x during the last accepted step
        r_cur::Array  = zeros(Float64, arr_len)               # Residuals (i)
        r_old::Array  = zeros(Float64, arr_len)               # Residuals (i-1)
        r_tst::Array  = zeros(Float64, arr_len)               # Test for rejection residuals
        dtd::Array    = zeros(Float64, (arr_len,arr_len))     # Damping matrix for LM method
        c_old::Float64 = 0.0

        # Final setup
        x_cur[:] .= x_ini[:]
        for di in 1:arr_len 
            dtd[di,di] = 1.0
        end 
        r_cur .+= 1.0e99
        r_old .+= 1.0e98

        @printf("    step  resid_med  resid_max  resid_2nm  xvals_med  xvals_max  deltx_2nm  \n")
        while true 

            # Update properties (cp, rho, etc.)
            if !all(isfinite, x_cur)
                display(x_cur)
                error("Solution array contains NaNs and/or Infs")
            end 
            _set_tmps!(x_cur)
            atmosphere.calc_layer_props!(atmos)

            # Update step counter
            step += 1
            if step > max_steps
                code = 1 
                break 
            end 
            if mod(step,modprint) == 0 
                @printf("    %4d  ", step)
            end

            # Evaluate jacobian and residuals
            r_old[:] .= r_cur[:]
            if use_cendiff || (step == 1)
                calc_jac_res_cendiff!(x_cur, b, r_cur) 
            else
                calc_jac_res_fordiff!(x_cur, b, r_cur) 
            end 

            # Check convergence
            c_old = cost(r_cur)
            if c_old < atol 
                code = 0
                @printf("\n")
                break
            end

            # Check if jacobian is singular 
            if abs(det(b)) < 1.0e-80
                code = 2
                @printf("\n")
                break
            end 

            # Model step 
            x_old[:] = x_cur[:]
            if (method == 0) || (step <= 2)
                # Newton-Raphson step 
                x_dif = -b\r_cur
                x_cur = x_old + x_dif

            elseif method == 1
                # Gauss-Newton step 
                x_dif = -(b'*b) \ (b'*r_cur) 
                x_cur = x_old + x_dif

            elseif method == 2
                # Levenberg-Marquardt step
                #    Calculate damping parameter ("delayed gratification")
                if r_cur_2nm < r_old_2nm
                    lml /= 5.0
                else 
                    lml *= 1.5
                end

                #    Update our estimate of the solution
                x_dif = -(b'*b + lml * dtd) \ (b' * r_cur)
                x_cur = x_old + x_dif
                fev!(x_cur, r_tst)

                #    Accept or reject this step (https://arxiv.org/pdf/1201.5885.pdf)
                # reject = (1.0 - dot(x_dif, x_dla)/(norm(x_dif)*norm(x_dla)) )^2.0 * cost(r_tst) > c_old
                reject = false
                if reject
                    # reject step  
                    x_cur[:] .= x_old[:]
                    lml *= 10.0
                    count_rej += 1
                else 
                    # accept step
                    x_dla = x_cur - x_old
                end
            end

            # Model statistics 
            fev!(x_cur, r_cur)
            r_med   = median(r_cur)
            r_max   = r_cur[argmax(abs.(r_cur))]
            x_med   = median(x_cur)
            x_max   = x_cur[argmax(abs.(x_cur))]
            dx2nm   = norm(x_dif)
            r_old_2nm   = r_cur_2nm
            r_cur_2nm   = norm(r_cur)

            # Plot
            if modplot > 0
                if mod(step, modplot) == 0
                    plotting.plot_pt(atmos,  @sprintf("%s/solver.png", atmos.OUT_DIR), incl_magma=(surf_state==2))
                end 
            end 
                
            # Inform user
            if mod(step,modprint) == 0 
                reject_str = "A"
                if reject
                    reject_str = "R"
                end
                @printf("%+.2e  %+.2e  %.3e  %+.2e  %+.2e  %.3e %s\n", r_med, r_max, r_cur_2nm, x_med, x_max, dx2nm, reject_str)
            end

        end # end solver loop
        
        # Check result
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
            println("    failure (local minimum)")
        else 
            println("    failure (unhandled)")
        end
        @printf("    steps rejected: %d \n", count_rej)
        println(" ")

        # ----------------------------------------------------------
        # Extract solution
        # ---------------------------------------------------------- 
        calc_cf = calc_cf_end
        fev!(x_cur, zeros(Float64, arr_len))
        atmosphere.calc_hrates!(atmos)
        # display(prate)

        # ----------------------------------------------------------
        # Print info
        # ---------------------------------------------------------- 
        loss = atmos.flux_tot[1] - atmos.flux_tot[end]
        loss_pct = loss/atmos.flux_tot[1]*100.0
        @printf("    endpoint fluxes \n")
        @printf("    rad_OLR   = %+.2e W m-2     \n", atmos.flux_u_lw[1])
        if (surf_state == 2)
            F_skin = atmos.skin_k / atmos.skin_d * (atmos.tmp_magma - atmos.tstar)
            @printf("    cond_skin = %+.2e W m-2 \n", F_skin)
        end
        @printf("    tot_TOA   = %+.2e W m-2     \n", atmos.flux_tot[1])
        @printf("    tot_BOA   = %+.2e W m-2     \n", atmos.flux_tot[end])
        @printf("    loss      = %+.2e W m-2     \n", loss)
        @printf("    loss      = %+.2e %%        \n", loss_pct)
        @printf("\n")

        return atmos.is_converged
    end # end solve_energy 

end 