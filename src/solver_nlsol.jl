# Contains code for obtaining energy balance (nonlinear method)

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module solver_nlsol

    using Printf
    using LoggingExtras
    using Statistics
    using Revise
    using LinearAlgebra

    import atmosphere 
    import energy
    import phys
    import plotting

    # Golden section search to minimise a function
    function gs_search(f::Function,a::Float64,b::Float64,dxtol::Float64,max_steps::Int)::Float64
        c::Float64 = (-1+sqrt(5))/2
    
        x1::Float64 = c*a + (1-c)*b
        x2::Float64 = (1-c)*a + c*b

        fx1::Float64 = f(x1)
        fx2::Float64 = f(x2)

        best::Float64 = 0.5*(a+b)
    
        for i = 1:max_steps
            if fx1 < fx2
                b = x2
                x2 = x1
                fx2 = fx1
                x1 = c*a + (1-c)*b
                fx1 = f(x1)
            else
                a = x1
                x1 = x2
                fx1 = fx2
                x2 = (1-c)*a + c*b
                fx2 = f(x2)
            end

            best = 0.5*(a+b)
    
            if abs(b-a) < dxtol
                @debug "GS search succeeded after $(i+2) function evaluations (best = $best)"
                break
            end
        end 
    
        return best
    end 

    """
    **Obtain radiative-convective equilibrium using a matrix method.**

    Solves the non-linear system of equations defined by the flux field
    divergence, minimising flux loss across a cell by iterating the temperature
    profile.

    Not compatible with convective adjustment; MLT must be used.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `sol_type::Int=1`                 solution type, 0: free | 1: fixed | 2: skin | 3: tmp_eff | 4: tgt_olr
    - `condensates::Array=[]`           condensates to model (if empty, no condensates are modelled)
    - `chem_type::Int=0`                chemistry type (see wiki)
    - `convect::Bool=true`              include convection
    - `sens_heat::Bool=false`           include sensible heating 
    - `conduct::Bool=false`             include conductive heat transport within the atmosphere
    - `dx_max::Float64=200.0`           maximum step size [K]
    - `max_steps::Int=2000`             maximum number of solver steps
    - `max_runtime::Float64=600.0`      maximum runtime in wall-clock seconds
    - `fdw::Float64=1.0e-4`             finite difference: relative width (dx/x) of the "difference"
    - `fdc::Bool=false`                 finite difference: always use central difference? 
    - `fdo::Int=2`                      finite difference: scheme order (2nd or 4th)
    - `method::Int=1`                   numerical method (1: Newton-Raphson, 2: Gauss-Newton, 3: Levenberg-Marquardt)
    - `linesearch::Bool=true`           use a golden-section linesearch algorithm to determine the best step size
    - `modplot::Int=0`                  iteration frequency at which to make plots
    - `save_frames::Bool=true`          save plotting frames
    - `modulate_mlt::Bool=true`         improve convergence with convection by introducing MLT gradually
    - `conv_atol::Float64=1.0e-5`       convergence: absolute tolerance on per-level flux deviation [W m-2]
    - `conv_rtol::Float64=1.0e-3`       convergence: relative tolerance on per-level flux deviation [dimensionless]
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                           sol_type::Int=1, condensates::Array=[],
                            condensing::Array=[], chem_type::Int=0,
                            convect::Bool=true, sens_heat::Bool=false,
                            conduct::Bool=false, 
                            dx_max::Float64=200.0, 
                            max_steps::Int=2000, max_runtime::Float64=600.0,
                            fdw::Float64=1.0e-4, fdc::Bool=false, fdo::Int=2,
                            method::Int=1, linesearch::Bool=true,
                            modplot::Int=1, save_frames::Bool=true, 
                            modulate_mlt::Bool=true,
                            conv_atol::Float64=1.0e-5, conv_rtol::Float64=1.0e-3
                            )::Bool

        # Validate condensation case
        do_condense::Bool  = false 
        if length(condensates) > 0
            for c in condensates
                if c in atmos.gas_all_names
                    do_condense = true 
                else 
                    @error "Invalid condensate '$c'"
                    return false
                end 
            end
        end 

        # Validate sol_type
        if (sol_type < 0) || (sol_type > 4)
            @error "Invalid solution type ($sol_type)"
            return false
        end

        # Start timer 
        wct_start::Float64 = time()

        # Plot paths 
        path_prf::String = @sprintf("%s/solver_prf.png", atmos.OUT_DIR)
        path_flx::String = @sprintf("%s/solver_flx.png", atmos.OUT_DIR)
        path_vmr::String = @sprintf("%s/solver_vmr.png", atmos.OUT_DIR)

        # Dimensionality
        arr_len::Int = atmos.nlev_c 
        if (sol_type >= 2)  # states 2,3,4 also solve for tmp_surf
            arr_len += 1
        end

        # Work arrays 
        rf1::Array{Float64,1}    = zeros(Float64, arr_len)  # Forward difference (1 step)
        rb1::Array{Float64,1}    = zeros(Float64, arr_len)  # Backward difference (1 step)
        rf2::Array{Float64,1}    = zeros(Float64, arr_len)  # Forward difference (2 step)
        rb2::Array{Float64,1}    = zeros(Float64, arr_len)  # Backward difference (2 step)
        x_s::Array{Float64,1}    = zeros(Float64, arr_len)  # Perturbed row, for jacobian
        fd_s::Float64            = 0.0                      # Row perturbation amount
        drdx::Float64            = 0.0                      # Jacobian element dr/dx
        do_chemistry::Bool       = (chem_type>0)
        
        # Convective flux scale factor 
        convect_sf::Float64 = 5.0e-5

        # Calculate the (remaining) temperatures  
        function _set_tmps!(_x::Array)
            # Read new guess
            clamp!(_x, atmos.tmp_floor+1.0, atmos.tmp_ceiling-1.0)
            for i in 1:atmos.nlev_c
                atmos.tmp[i] = _x[i]
            end

            # Interpolate temperature to cell-edge values (not incl. bottommost value)
            atmosphere.set_tmpl_from_tmp!(atmos)

            # Set bottom edge temperature 
            if (sol_type != 1) 
                # For state=1, tmpl[end] is held constant

                # Extrapolate (log-linear)
                grad_dt = atmos.tmp[end]-atmos.tmp[end-1]
                grad_dp = log(atmos.p[end]/atmos.p[end-1])
                atmos.tmpl[end] = atmos.tmp[end] + grad_dt/grad_dp * log(atmos.pl[end]/atmos.p[end])

                if (sol_type >= 2)  # states 2,3,4
                    atmos.tmp_surf = _x[end]  # Surface brightness temperature
                end
            end 
            
            return nothing
        end # end set_tmps

        # Objective function to solve for
        function _fev!(x::Array,resid::Array)

            # Reset masks
            fill!(atmos.mask_c, 0.0)
            fill!(atmos.mask_p, 0.0)
           
            # Set temperatures 
            _set_tmps!(x)

            # Calculate layer properties
            if atmos.thermo_funct
                atmosphere.calc_layer_props!(atmos)
            end 

            # Calculate fluxes
            energy.calc_fluxes!(atmos, 
                                do_condense,  convect, sens_heat, conduct, 
                                condensates=condensates, condensing=condensing,
                                convect_sf=convect_sf)

            # Additional energy input
            atmos.flux_dif[:] .+= atmos.ediv_add[:] .* atmos.layer_thick

            # Calculate residuals subject to the solution type
            if (sol_type == 0) || (sol_type == 1)
                # Zero loss with constant tmp_surf
                resid[1:end] .= atmos.flux_dif[1:end]

            elseif (sol_type == 2)
                # Conductive boundary layer
                resid[1:end-1] .= atmos.flux_dif[1:end]
                resid[end] = atmos.flux_tot[end] - (atmos.tmp_magma - atmos.tmpl[end]) * atmos.skin_k / atmos.skin_d

            elseif (sol_type == 3)
                # Zero loss
                resid[2:end] .= atmos.flux_dif[1:end]
                # Total flux at TOA is equal to sigma*tmp_eff^4
                resid[1] = atmos.flux_tot[1] - atmos.flux_eff

            elseif (sol_type == 4)
                # Zero loss
                resid[2:end] .= atmos.flux_dif[1:end]
                # OLR is equal to target_olr
                resid[1] = atmos.target_olr - atmos.flux_u_lw[1]

            end

            # Check that residuals are real numbers
            if !all(isfinite, resid)
                display(resid)
                error("Residual array contains NaNs and/or Infs")
            end

            return nothing
        end  # end fev

        # Calculate the jacobian and residuals at x using a 2nd order central-difference method
        function _calc_jac_res!(x::Array, jacob::Array, resid::Array, central::Bool, order::Int)

            # Evalulate residuals at x
            _fev!(x, resid)

            # For each level...
            for i in 1:arr_len 
                
                # Reset all levels
                x_s[:] .= x[:]

                # Reset residuals
                fill!(rf1, 0.0)
                fill!(rb1, 0.0)
                fill!(rf2, 0.0)
                fill!(rb2, 0.0)

                # Calculate perturbation at this level
                fd_s = x[i] * fdw

                # Forward part (1 step)
                x_s[i] = x[i] + fd_s
                _fev!(x_s, rf1)

                # Forward part (2 step)
                if order == 4
                    x_s[i] = x[i] + 2.0*fd_s
                    _fev!(x_s, rf2)
                end 

                # Backward part
                if central
                    # (1 step)
                    x_s[i] = x[i] - fd_s
                    _fev!(x_s, rb1)

                    # (2 step)
                    if order == 4
                        x_s[i] = x[i] - 2.0*fd_s
                        _fev!(x_s, rb2)
                    end 
                end

                # Set jacobian
                # https://www.dam.brown.edu/people/alcyew/handouts/numdiff.pdf
                # https://en.wikipedia.org/wiki/Finite_difference_coefficient
                for j in 1:arr_len  # for each residual
                    jacob[j,i] = drdx
                    if central 
                        if order == 4
                            # 4th order central difference 
                            jacob[j,i] = (-rf2[j] + 8.0*rf1[j] - 8.0*rb1[j] + rb2[j])/(12.0*fd_s)
                        else 
                            # 2nd order central difference 
                            jacob[j,i] = (rf1[j] - rb1[j])/(2.0*fd_s)
                        end 
                    else
                        if order == 4
                            # 4th order forward difference
                            jacob[j,i] = (-rf2[j] + 4.0*rf1[j] - 3.0*resid[j])/(2.0*fd_s)
                        else 
                            # 2nd order forward difference
                            jacob[j,i] = (rf1[j] - resid[j])/fd_s
                        end  
                    end 
                end 
            end 

            return nothing
        end # end jr_cd

        # Cost function 
        function _cost(_r::Array)
            return norm(_r)
        end 

        # Plot
        title_prf::String = ""
        title_flx::String = ""
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

            if do_chemistry || (length(condensates) > 0)
                plotting.plot_vmr(atmos, path_vmr)
            end 
        end 

        
        # ----------------------------------------------------------
        # Setup initial guess
        # ---------------------------------------------------------- 
            @info @sprintf("    sol_type = %d", sol_type)
        if (sol_type == 1)
            @info @sprintf("    tmp_surf = %.2f K", atmos.tmp_surf)
        elseif (sol_type == 2)
            @info @sprintf("    skin_d   = %.2f m",         atmos.skin_d)
            @info @sprintf("    skin_k   = %.2f W K-1 m-1", atmos.skin_k)
        elseif (sol_type == 3)
            @info @sprintf("    tmp_eff  = %.2f K",     atmos.tmp_eff)
            @info @sprintf("    f_eff    = %.2f W m-2", atmos.flux_eff)
        elseif (sol_type == 4)
            @info @sprintf("    tgt_olr  = %.2f W m-2", atmos.target_olr)
        end 

        # Allocate initial guess for the x array, as well as a,b arrays
        # Array storage structure:
        #   in the sol_type=2 case
        #       1:end-1 => cell centre temperatures 
        #       end     => bottom cell edge temperature
        #   other cases 
        #       1:end => cell centre temperatures
        x_ini    = zeros(Float64, arr_len) 
        for i in 1:atmos.nlev_c
            x_ini[i]    = clamp(atmos.tmp[i], atmos.tmp_floor , atmos.tmp_ceiling)
        end 
        if (sol_type >= 2)
            x_ini[end] = atmos.tmp[atmos.nlev_c] + 1.0
        end

        # ----------------------------------------------------------
        # Solver loop
        # ---------------------------------------------------------- 
        # Execution variables
        modprint::Int =         1       # Print frequency
        convect_incr::Float64 = 6.0     # Factor to increase convect_sf when modulating convection
        
        # Tracking variables
        step::Int =         0       # Step number
        code::Int =         -1      # Status code 
        lml::Float64 =      2.0     # Levenberg-Marquardt lambda parameter
        runtime::Float64  = 0.0     # Model runtime [s]
        fc_retcode::Int  =  0       # Fastchem return code
        step_ok::Bool =     true    # Current step was fine
        struggling::Bool =  false   # Model is (or was) struggling

        # Model statistics tracking
        r_med::Float64 =        9.0     # Median residual
        r_max::Float64 =        9.0     # Maximum residual (sign agnostic)
        c_max::Float64 =        0.0     # Maximum cost (sign agnostic)
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
        dx_max_step::Float64 = Inf  # maximum step size (IN THIS STEP)

        # Linesearch parameters
        ls_best_scale::Float64  = 1.0       # Best found scale
        ls_max_steps::Int       = 12        # Maximum golden section steps 
        ls_begin_step::Int      = 2

        # Final setup
        x_cur[:] .= x_ini[:]
        for di in 1:arr_len 
            dtd[di,di] = 1.0
        end 
        fill!(r_cur, 1.0e99)
        fill!(r_old, 1.0e98)

        # Modulate convection?
        modulate_mlt = modulate_mlt && convect
        if !modulate_mlt
            convect_sf = 1.0
        end

        # Initial plot 
        if modplot > 0
            plot_step(0, 0.0)
        end

        @info @sprintf("    step  resid_med  resid_2nm  flux_OLR   xvals_med  xvals_max  |dx|_max   flags")
        info_str::String = ""
        stepflags::String = ""
        while true 

            # Reset flags 
            #     Cs, Cf     = chemistry model (s)uccess, (f)ailure
            #     Mc,Mr      = Modulate convection: (r)educed this step
            #     Td         = finite differencing type (T) and order (d)
            #     Nr,Gn,Lm   = stepping algorithm
            #     Ls         = linesearch active
            info_str  = ""
            stepflags = ""
            step_ok   = true 

            # Check time 
            runtime = time()-wct_start
            if runtime > max_runtime 
                code = 3
                break
            end 
            struggling = struggling || (runtime > max_runtime*0.5)

            # Update step counter
            step += 1
            if step > max_steps
                code = 1 
                break 
            end 
            if mod(step,modprint) == 0 
                info_str *= @sprintf("    %4d  ", step)
            end
            struggling = struggling || (step > max_steps*0.5)

            # Enable LS now if struggling
            if struggling && linesearch
                ls_begin_step = 0
            end 

            # Check status of guess 
            if !all(isfinite, x_cur)
                code = 4 
                break
            end 
            clamp!(x_cur, atmos.tmp_floor+10.0, atmos.tmp_ceiling-10.0)
            _set_tmps!(x_cur)

            # Run chemistry scheme 
            if chem_type in [1,2,3]
                fc_retcode = atmosphere.chemistry_eq!(atmos, chem_type, false)
                if fc_retcode == 0
                    stepflags *= "Cs-"
                else 
                    stepflags *= "Cf-"
                    step_ok = false
                end
            end 

            # Do condensation here (doesn't affect fluxes currently)

            # Reset condensing array
            condensing = [String[] for x in 1:atmos.nlev_c]
            atmosphere.condense_varyx!(atmos, condensing, condensates=condensates, 
                                        gases=atmos.gas_all_names)
            
            # Update properties (mmw, density, height, etc.)
            step_ok = step_ok && atmosphere.calc_layer_props!(atmos)

            # Check convective modulation
            if modulate_mlt 
                # We are modulating
                stepflags *= "M"
                
                # Check if sf needs to be increased
                if c_cur < max(100.0*conv_rtol, 10.0)
                    if convect_sf < 1.0 
                        # increase sf - reduce modulation by 10x
                        stepflags *= "r"
                        convect_sf = min(1.0, convect_sf*convect_incr)
                        @debug "convect_sf = $convect_sf"
                        
                        # done modulating 
                        if convect_sf > 0.99
                            modulate_mlt = false 
                        end
                    end 
                end

                if stepflags[end] == "M"
                    stepflags *= "c"
                end 
                stepflags *= "-"
            else
                # No modulation at this point
                convect_sf = 1.0 
            end 

            # Evaluate residuals and estimate jacobian with finite-difference 
            r_old[:] .= r_cur[:]
            if fdc || (step == 1) || struggling || (c_cur > c_old)
                # use central difference if: requested, at the start, struggling, or cost increased
                _calc_jac_res!(x_cur, b, r_cur, true, fdo)
                stepflags *= "C$fdo-"
            else
                # otherwise, use forward difference
                _calc_jac_res!(x_cur, b, r_cur, false, fdo)
                stepflags *= "F$fdo-"
            end 

            # Check if jacobian is singular 
            if abs(det(b)) < floatmin()*10.0
                code = 2
                step_ok = false 
                break
            end 

            # Model step 
            dx_max_step = dx_max
            x_old[:] = x_cur[:]
            if (method == 1)
                # Newton-Raphson step 
                x_dif = -b\r_cur
                stepflags *= "Nr-"

            elseif method == 2
                # Gauss-Newton step 
                x_dif = -(b'*b) \ (b'*r_cur) 
                stepflags *= "Gn-"

            elseif method == 3
                # Levenberg-Marquardt step
                #    Calculate damping parameter ("delayed gratification")
                if r_cur_2nm < r_old_2nm
                    lml /= 5.0
                else 
                    lml *= 2.5
                end

                #    Update our estimate of the solution
                x_dif = -(b'*b + lml * dtd) \ (b' * r_cur)
                stepflags *= "Lm-"
            end

            # Limit step size, without changing direction of dx vector
            x_dif[:] .*= min(1.0, dx_max_step / maximum(abs.(x_dif[:])))

            # Linesearch 
            if linesearch && (step >= ls_begin_step)

                # Reset
                stepflags *= "Ls-"

                # Linesearch cost function
                function _costls!(scale::Float64)::Float64
                    x_cur[:] .= x_old[:] .+ scale .* x_dif[:]
                    _fev!(x_cur,r_tst)
                    return _cost(r_tst)
                end 

                # Golden section search
                ls_best_scale = gs_search(_costls!,0.05,1.0,1.0e-2,ls_max_steps)

                # apply best linesearch scale 
                x_dif[:] .*= ls_best_scale

            end # end linesearch 

            # Take the step 
            x_cur[:] .= x_old[:] .+ x_dif[:]
            clamp!(x_cur, atmos.tmp_floor+10.0, atmos.tmp_ceiling-10.0)
            _fev!(x_cur, r_cur)

            # Cost value
            c_old = c_cur
            c_cur = _cost(r_cur)

            # Model statistics 
            r_med   = median(r_cur)
            iworst  = argmax(abs.(r_cur))
            r_max   = r_cur[iworst]
            x_med   = median(x_cur)
            x_max   = x_cur[argmax(abs.(x_cur))]
            dxmax   = maximum(abs.(x_dif))
            r_old_2nm   = r_cur_2nm
            r_cur_2nm   = norm(r_cur)
            c_max = maximum(abs.(atmos.flux_tot))

            # Plot
            if (modplot > 0) && (mod(step, modplot) == 0)
                plot_step(step, runtime)
            end 
                
            # Inform user
            if mod(step,modprint) == 0 
                info_str *= @sprintf("%+.2e  %.3e  %.3e  %+.2e  %+.2e  %.3e  %-s", r_med, r_cur_2nm, atmos.flux_u_lw[1], x_med, x_max, dxmax, stepflags[1:end-1])
                if step_ok
                    @info info_str
                else
                    @warn info_str 
                end
            end

            # Converged?
            if (c_cur < conv_atol + conv_rtol * c_max) && !modulate_mlt
                code = 0
                break
            end

        end # end solver loop
        
        rm(path_prf, force=true)
        rm(path_flx, force=true)
        rm(path_vmr, force=true)
        
        # ----------------------------------------------------------
        # Extract solution
        # ---------------------------------------------------------- 
        atmos.is_solved = true
        atmos.is_converged = false 
        if code == 0
            @info "    success"
            atmos.is_converged = true
        elseif code == 1
            @error "    failure (maximum iterations)"
        elseif code == 2
            @error "    failure (singular jacobian)"
        elseif code == 3
            @error "    failure (maximum time)"
        elseif code == 4
            @error "    failure (NaN temperature)"
        else 
            @error "    failure (unhandled)"
        end

        _fev!(x_cur, zeros(Float64, arr_len))
        energy.calc_hrates!(atmos)

        # ----------------------------------------------------------
        # Print info
        # ---------------------------------------------------------- 
        loss = maximum(atmos.flux_tot) - minimum(atmos.flux_tot)
        loss_pct = 100.0*loss/maximum(atmos.flux_tot)
        @info @sprintf("    outgoing LW flux   = %+.2e W m-2     ", atmos.flux_u_lw[1])
        if (sol_type == 2)
            F_skin = atmos.skin_k / atmos.skin_d * (atmos.tmp_magma - atmos.tmp_surf)
            @info @sprintf("    conduct. skin flux = %+.2e W m-2 ", F_skin)
        end
        @info @sprintf("    total flux at TOA  = %+.2e W m-2     ", atmos.flux_tot[1])
        @info @sprintf("    total flux at BOA  = %+.2e W m-2     ", atmos.flux_tot[end])
        @info @sprintf("    column max. loss   = %+.2e W m-2  (%+.2e %%) ", loss, loss_pct)
        @info @sprintf("    final cost value   = %+.2e W m-2     ", c_cur)
       

        return atmos.is_converged
    end # end solve_energy 
end 
