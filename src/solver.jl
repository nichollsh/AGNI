# Contains code for obtaining energy balance (nonlinear method)

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module solver

    using Printf
    using LoggingExtras
    using Statistics
    using LinearAlgebra

    import ..atmosphere
    import ..setpt
    import ..energy
    import ..phys
    import ..plotting
    import ..chemistry
    import ..ocean

    """
    **Golden section search algorithm**

    Minimises a function `f` between the bounds `a` and `b`.
    The function `f` must return a positive scalar.

    Arguments:
    - `f`           Function to be minimised
    - `a`           Initial bracket, lower value
    - `b`           Initial bracket, upper valuie
    - `dxtol`       Convergence: exit when bracket is smaller than this size
    - `atol`        Convergence: absolute tolerance on minimum
    - `max_steps`   Maximum number of iterations
    - `warnings`    Print warning if search does not converge before `max_steps`` are taken.
    """
    function gs_search(f::Function,a::Float64,b::Float64,
                            dxtol::Float64,atol::Float64,max_steps::Int;
                            warnings::Bool=false)::Float64
        c::Float64 = (-1+sqrt(5))/2

        x1::Float64 = c*a + (1-c)*b
        x2::Float64 = (1-c)*a + c*b

        fx1::Float64 = f(x1)
        fx2::Float64 = f(x2)

        midp::Float64 = 0.5*(a+b)

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


            if fx1 < atol
                # @debug "GS search found minimum. $(i+2) function evaluations, best = $x1"
                return x1
            end

            if fx2 < atol
                # @debug "GS search found minimum. $(i+2) function evaluations, best = $x2"
                return x2
            end

            midp = 0.5*(a+b)

            if abs(b-a) < dxtol
                # @debug "GS search reached minimum bracket size. $(i+2) function evaluations, best = $best"
                return midp
            end
        end

        if warnings
            @warn "GS search did not converge"
        end
        return midp
    end

    # Status codes for solver
    @enum STATUSCODE begin
        CODE_SUC = 0  # success
        CODE_ITE = 1  # max iters
        CODE_SIN = 2  # singular jacobian
        CODE_TIM = 3  # max time
        CODE_NAN = 4  # NaN values returned from objective
        CODE_CFG = 5  # configuration
        CODE_OBJ = 6  # objective function did not succeed
        CODE_STP = 7  # last step not ok
        CODE_HYD = 8  # hydrostatic integration

        CODE_99  = 99 # default code (failure)
    end

    # Solver constants and parameters
    #    plateau parameters
    plateau_n::Int     =    4       # Plateau declared when plateau_i > plateau_n
    plateau_s::Float64 =    3.0     # Scale factor applied to x_dif when plateau_i > plateau_n
    plateau_r::Float64 =    0.98    # Cost ratio for determining whether to increment plateau_i
    #    linesearch
    ls_tau::Float64    =    0.7     # backtracking downscale size

    """
    **Obtain radiative-convective equilibrium using a matrix method.**

    Solves the non-linear system of equations defined by the flux field
    divergence, minimising flux loss across a cell by iterating the temperature
    profile.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.

    Optional arguments:
    - `sol_type::Int`                   solution type, 1: tmp_surf | 2: skin | 3: flux_int | 4: tgt_olr
    - `chem::Bool`                      include eqm thermochemistry when solving for RCE?
    - `convect::Bool`                   include convection
    - `sens_heat::Bool`                 include sensible heating at the surface
    - `conduct::Bool`                   include conductive heat transport within the atmosphere
    - `latent::Bool`                    include latent heat exchange (condensation/evaporation)
    - `advect::Bool`                    include advective heat exchange (dynamics)
    - `rainout::Bool`                   allow rainout (phase change impacts mixing ratios, not just energy fluxes)
    - `dx_min::Float64`                 minimum step size [K]
    - `dx_max::Float64`                 maximum step size [K]
    - `tmp_pad::Float64`                padding around hard limits on temperature floor & ceiling values
    - `max_steps::Int`                  maximum number of solver steps
    - `max_runtime::Float64`            maximum runtime in wall-clock seconds
    - `fd_rel::Float64`                 finite difference: relative width (dx/x) of the difference (rtol)
    - `fd_abs::Float64`                 finite difference: absolute width (dx) of the difference (atol)
    - `fdc::Bool`                       finite difference: central difference? otherwise use forward difference
    - `fdo::Int`                        finite difference: scheme order (2nd or 4th)
    - `method::Int`                     numerical method (1: Newton-Raphson, 2: Gauss-Newton, 3: Levenberg-Marquardt)
    - `easy_start::Bool`                improve convergence reliability; introduce convection and latent heat gradually
    - `grey_start::Bool`                improve convergence reliability; obtain double-grey solution first
    - `ls_method::Int`                  linesearch algorithm (0: None, 1: golden, 2: backtracking)
    - `ls_increase::Bool`               factor by which the cost can increase from last step before triggering linesearch
    - `ls_max_steps::Int`               maximum steps undertaken by linesearch routine
    - `ls_min_scale::Float64`           minimum step scale allowed by linesearch
    - `ls_max_scale::Float64`           maximum step scale allowed by linesearch
    - `detect_plateau::Bool`            assist solver when it is stuck in a region of small dF/dT
    - `perturb_all::Bool`               always recalculate entire Jacobian matrix? Otherwise updates columns only as required
    - `modplot::Int`                    iteration frequency at which to make plots
    - `save_frames::Bool`               save plotting frames
    - `modprint::Int`                   iteration frequency at which to print info
    - `plot_jacobian::Bool`             plot jacobian too?
    - `radiative_Kzz::Bool`             estimate Kzz in radiative zones? Otherwise left at Kzz=0 in these regions.
    - `conv_atol::Float64`              convergence: absolute tolerance on per-level flux deviation [W m-2]
    - `conv_rtol::Float64`              convergence: relative tolerance on per-level flux deviation [dimensionless]

    Returns:
        Nothing
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            sol_type::Int=1,
                            chem::Bool=false,
                            convect::Bool=true, sens_heat::Bool=true,
                            conduct::Bool=true, latent::Bool=true, advect::Bool=true,
                            rainout::Bool=true,
                            dx_min::Float64=1e-7, dx_max::Float64=400.0,
                            tmp_pad::Float64 = 5.0,
                            max_steps::Int=400, max_runtime::Float64=900.0,
                            fd_rel::Float64=3.0e-5, fd_abs::Float64=1e-5,
                            fdc::Bool=true, fdo::Int=2,
                            method::Int=1,
                            easy_start::Bool=false, grey_start::Bool=false,
                            ls_method::Int=1,
                            ls_increase::Float64=1.08,
                            ls_max_steps::Int=10,
                            ls_min_scale::Float64 = 1.0e-5, ls_max_scale::Float64 = 0.99,
                            detect_plateau::Bool=true, perturb_all::Bool=true,
                            modplot::Int=1, save_frames::Bool=true,
                            modprint::Int=1, plot_jacobian::Bool=true,
                            radiative_Kzz::Bool=true,
                            conv_atol::Float64=1.0e-1, conv_rtol::Float64=1.0e-3,
                            )::Bool

        # --------------------
        # Prepare to run the solver
        # --------------------
        @debug "Reticulating splines..."

        # Validate sol_type
        if (sol_type < 1) || (sol_type > 4)
            @error "Invalid solution type ($sol_type)."
            return false
        end

        # Validate finite difference fdo option
        if !(fdo in [2,4])
            @error "Invalid finite-difference order ($fdo). Must be 2 or 4."
            return false
        end

        # Validate if transparent
        if atmos.transparent
            @error "Atmosphere is configured to be transparent."
            @error "    Cannot use `solve_energy`. Use `solve_transparent` instead."
            return false
        end

        # Validate options for rainout (required for latent heat)
        if latent && !rainout
            @error "Must enable rainout if also including latent heating"
            return false
        end

        # Start timer
        wct_start::Float64 = time()

        # Plot paths
        path_plt::String = joinpath(atmos.OUT_DIR,"solver.png")
        path_jac::String = joinpath(atmos.OUT_DIR,"jacobian.png")

        # Dimensionality
        arr_len::Int = atmos.nlev_c
        if sol_type in [2,3,4]  # states 2,3,4 also solve for tmp_surf
            arr_len += 1
        end

        # --------------------
        # Execution parameters
        # --------------------

        #    grey_start
        rt_default::atmosphere.RTSCHEME = atmos.rt_scheme  # record default selection for RT scheme
        if grey_start
            @debug "Starting with double-grey RT scheme"
            atmos.rt_scheme = atmosphere.RT_GREYGAS
        end

        #    easy_start
        easy_start         = Bool(easy_start)
        easy_incr::Float64 = 2.0        # Factor by which to increase easy_sf at each step
        easy_trig::Float64 = 0.1        # Increase sf when cost*easy_trig satisfies convergence

        #    finite difference
        perturb_trig::Float64 = 0.1     # Require full Jacobian update when cost*peturb_trig satisfies convergence
        perturb_crit::Float64 = 0.1     # Require Jacobian update at level i when r_i>perturb_crit
        perturb_mod::Int =      5       # Do full jacobian at least this frequently

        # --------------------
        # Execution variables
        # --------------------
        #     finite difference
        rf1::Array{Float64,1}    = zeros(Float64, arr_len)  # Forward difference  (+1 dx)
        rb1::Array{Float64,1}    = zeros(Float64, arr_len)  # Backward difference (-1 dx)
        rf2::Array{Float64,1}    = zeros(Float64, arr_len)  # Forward difference  (+2 dx)
        rb2::Array{Float64,1}    = zeros(Float64, arr_len)  # Backward difference (-2 dx)
        x_s::Array{Float64,1}    = zeros(Float64, arr_len)  # Perturbed row, for jacobian
        fd_s::Float64            = 0.0                      # Row perturbation amount

        #     solver
        b::Array{Float64,2}      = zeros(Float64, (arr_len, arr_len))   # Approximate jacobian (i)
        x_cur::Array{Float64,1}  = zeros(Float64, arr_len)              # Current best solution (i)
        x_old::Array{Float64,1}  = zeros(Float64, arr_len)              # Previous best solution (i-1)
        x_dif::Array{Float64,1}  = zeros(Float64, arr_len)              # Change in x (i-1 to i)
        r_cur::Array{Float64,1}  = zeros(Float64, arr_len)              # Residuals (i)
        r_old::Array{Float64,1}  = zeros(Float64, arr_len)              # Residuals (i-1)
        r_tst::Array{Float64,1}  = zeros(Float64, arr_len)              # Test for rejection residuals
        dtd::Array{Float64,2}    = zeros(Float64, (arr_len,arr_len))    # Damping matrix for LM method
        perturb::Array{Bool,1}   = falses(arr_len)      # Mask for levels which should be perturbed
        lml::Float64             = 2.0                  # Levenberg-Marquardt lambda parameter
        c_cur::Float64           = Inf                  # current cost (i)
        c_old::Float64           = Inf                  # old cost (i-1)
        dx_max_step::Float64     = Float64(dx_max)      # dx_max allowed in current step
        linesearch::Bool         = Bool(ls_method>0)    # ls enabled?
        ls_alpha::Float64        = 1.0                  # linesearch scale factor
        ls_cost::Float64         = 1.0e99               # linesearch cost
        easy_sf::Float64         = 0.0                  # Convective & phase change flux scale factor
        plateau_apply::Bool      = false                # Plateau declared in this iteration?

        #     tracking
        step::Int =             0       # Step number
        code::STATUSCODE =      CODE_99 # Status code
        runtime::Float64  =     0.0     # Model runtime [s]
        fc_retcode::Int  =      0       # Fastchem return code
        step_ok::Bool =         true    # Current step was fine
        grey_step::Bool = grey_start    # double-grey RT enabled in this step
        easy_step::Bool =       false   # easy_start sf increased in this step
        plateau_i::Int =        0       # Number of iterations for which step was small

        #      statistics
        r_med::Float64 =        9.0     # Median residual
        r_max::Float64 =        9.0     # Maximum residual (sign agnostic)
        c_max::Float64 =        0.0     # Maximum cost (sign agnostic)
        x_med::Float64 =        0.0     # Median solution
        x_max::Float64 =        0.0     # Maximum solution (sign agnostic)
        iworst::Int =           0       # Level which is furthest from convergence
        dx_stat::Float64 =      9.0     # Maximum change in solution array
        r_cur_2nm::Float64 =    0.01    # Two-norm of residuals
        r_old_2nm::Float64 =    0.02    # Previous ^

        # Calculate the (remaining) temperatures from known temperatures
        function _set_tmps!(_x::Array{Float64,1})
            # Read new guess
            clamp!(_x, atmos.tmp_floor+1.0, atmos.tmp_ceiling-1.0)
            for i in 1:atmos.nlev_c
                atmos.tmp[i] = _x[i]
            end

            # Interpolate temperature to cell-edge values
            atmosphere.set_tmpl_from_tmp!(atmos)

            if (sol_type >= 2)  # states 2,3,4
                atmos.tmp_surf = _x[end]  # Surface brightness temperature
            end

            return nothing
        end # end set_tmps

        # Objective function
        function _fev!(x::Array{Float64,1},resid::Array{Float64,1})::Bool

            # Reset masks
            fill!(atmos.mask_c, false)
            fill!(atmos.mask_l, false)

            # Set new temperatures
            _set_tmps!(x)

            # Do condensation and chemistry
            chemistry.calc_composition!(atmos, rainout, chem, rainout)

            # Calculate fluxes
            step_ok &= energy.calc_fluxes!(atmos, radiative=true,
                                latent=latent, convect=convect, sens_heat=sens_heat,
                                conduct=conduct, advect=advect,
                                convect_sf=easy_sf, latent_sf=easy_sf,
                                rainout=rainout)

            # Calculate residuals subject to the solution type
            if (sol_type == 1)
                # Zero loss with constant tmp_surf
                resid[1:end] .= atmos.flux_dif[1:end]

            elseif (sol_type == 2)
                # Zero loss
                resid[1:end-1] .= atmos.flux_dif[1:end]
                # Conductive boundary layer
                resid[end] = atmos.flux_tot[end] - energy.skin_flux(atmos)

            elseif (sol_type == 3)
                # Zero loss
                resid[2:end] .= atmos.flux_dif[1:end]
                resid[1] = atmos.flux_tot[1] - atmos.flux_int

            elseif (sol_type == 4)
                # Zero loss
                resid[2:end] .= atmos.flux_dif[1:end]
                # OLR is equal to target_olr
                resid[1] = atmos.target_olr - atmos.flux_u_lw[1]

            end

            # Check that residuals are real numbers
            if !all(isfinite, resid)
                @error "Residual array contains NaNs and/or Infs"
                @error "resid: $resid"
                @error "flux_n: $(atmos.flux_n)"
                code = CODE_NAN
                return false
            end

            return true
        end  # end fev

        # Calculate the jacobian and residuals at x using a 2nd order central-difference
        function _calc_jac_res!(x::Array{Float64, 1}, jacob::Array{Float64, 2},
                                    resid::Array{Float64 ,1}, central::Bool, order::Int,
                                    which::Array{Bool,1})::Bool

            ok::Bool = true

            # Evalulate residuals at x
            ok = ok && _fev!(x, resid)

            # For each level...
            for i in 1:arr_len
                if !which[i]
                    continue
                end

                # Reset all levels
                @. x_s = x

                # Reset residuals
                fill!(rf1, 0.0)
                fill!(rb1, 0.0)
                fill!(rf2, 0.0)
                fill!(rb2, 0.0)

                # Calculate perturbation at this level
                #    - must be less than tmp_pad/2
                fd_s = min(x[i] * fd_rel + fd_abs, tmp_pad/2)

                # Forward part (1 step)
                x_s[i] = x[i] + fd_s
                ok = ok && _fev!(x_s, rf1)

                # Forward part (2 step)
                if order == 4
                    x_s[i] = x[i] + 2.0*fd_s
                    ok = ok && _fev!(x_s, rf2)
                end

                # Backward part
                if central
                    # (1 step)
                    x_s[i] = x[i] - fd_s
                    ok = ok && _fev!(x_s, rb1)

                    # (2 step)
                    if order == 4
                        x_s[i] = x[i] - 2.0*fd_s
                        ok = ok && _fev!(x_s, rb2)
                    end
                end

                # Set jacobian
                # https://www.dam.brown.edu/people/alcyew/handouts/numdiff.pdf
                # https://en.wikipedia.org/wiki/Finite_difference_coefficient
                if central
                    if order == 4
                        # 4th order central difference
                        @. jacob[:,i] = (-rf2 + 8.0*rf1 - 8.0*rb1 + rb2)/(12.0*fd_s)
                    else
                        # 2nd order central difference
                        @. jacob[:,i] = (rf1 - rb1)/(2.0*fd_s)
                    end
                else
                    if order == 4
                        # 4th order forward difference
                        @. jacob[:,i] = (-rf2 + 4.0*rf1 - 3.0*resid)/(2.0*fd_s)
                    else
                        # 2nd order forward difference
                        @. jacob[:,i] = (rf1 - resid)/fd_s
                    end
                end # end central/forward
            end # end levels

            return ok
        end # end jr_cd

        # Cost function to minimise
        function _cost(_r::Array)
            return norm(_r,3)
        end

        # Plot current state
        function plot_step()

            # plotting.plot_cloud(atmos, "out/cloud.png")

            # Info string
            plt_info::String = ""
            plt_info *= @sprintf("Iteration:%d,  ",step)
            plt_info *= @sprintf("Runtime:%.1f s,  ",runtime)
            plt_info *= @sprintf("Cost:%.2e  ",c_cur)

            # Make subplots (don't save to file)
            plt_pt = plotting.plot_pt(atmos,     "", incl_magma=(sol_type==2))
            plt_fl = plotting.plot_fluxes(atmos, "", incl_eff=(sol_type==3), incl_cdct=conduct, incl_latent=latent)
            plt_mr = plotting.plot_vmr(atmos,    "")
            plt_ra = plotting.plot_radius(atmos, "")

            # Combined plot
            plotting.combined(plt_pt, plt_fl, plt_mr, plt_ra, plt_info, path_plt)

            if save_frames
                cp(path_plt,@sprintf("%s/%04d.png",atmos.FRAMES_DIR,step))
            end
        end


        # ----------------------------------------------------------
        # Setup initial guess
        # ----------------------------------------------------------
            @info @sprintf("    sol_type  = %d", sol_type)
        if (sol_type == 1)
            @info @sprintf("    tmp_surf  = %.2f K", atmos.tmp_surf)
        elseif (sol_type == 2)
            @info @sprintf("    skin_d    = %.2f m",         atmos.skin_d)
            @info @sprintf("    skin_k    = %.2f W K-1 m-1", atmos.skin_k)
        elseif (sol_type == 3)
            @info @sprintf("    flux_int  = %.2f W m-2", atmos.flux_int)
        elseif (sol_type == 4)
            @info @sprintf("    tgt_olr   = %.2f W m-2", atmos.target_olr)
        end

        # Allocate initial guess for the x array, as well as a,b arrays
        # Array storage structure:
        #   in the sol_type>=2 cases
        #       1:end-1 => cell centre temperatures
        #       end     => bottom cell edge temperature
        #   sol_type == 1 case
        #       1:end => cell centre temperatures
        x_ini    = zeros(Float64, arr_len)
        for i in 1:atmos.nlev_c
            x_ini[i] = clamp(atmos.tmp[i], atmos.tmp_floor, atmos.tmp_ceiling)
        end
        if (sol_type >= 2)
            x_ini[end] = atmos.tmp[atmos.nlev_c] + 1.0
        end

        # Final setup
        @. x_cur = x_ini
        for di in 1:arr_len
            dtd[di,di] = 1.0
        end
        fill!(r_cur, 1.0e99)            # reset residual arrays
        fill!(r_old, 1.0e98)            # ^
        energy.reset_fluxes!(atmos)     # reset energy fluxes

        # Modulate convection and phase change?
        easy_start = easy_start && (convect || latent)
        if !easy_start
            easy_sf = 1.0
        end

        # Initial plot
        if modplot > 0
            plot_step()
        end

        if modprint > 0
            @info @sprintf("    step  resid_med    cost     flux_OLR    max(x)    max(|dx|)   flags")
        else
            @info "    please wait..."
        end
        info_str::String = ""
        stepflags::String = ""
        while true

            # Reset flags
            #     Cs, Cf     = chemistry model (s)uccess, (f)ailure
            #     Mc,Mr      = Modulate convection: (r)educed this step
            #     Td         = finite differencing type (T) and order (d)
            #     Nr,Gn,Lm   = stepping algorithm
            #     Ls         = linesearch active
            #     X          = e(x)trapolating along plateau region
            info_str  = ""
            stepflags = ""
            step_ok   = true

            # Check time
            runtime = time()-wct_start
            if runtime > max_runtime
                code = CODE_TIM
                break
            end

            # Update step counter
            @debug "        iterate"
            step += 1
            if step > max_steps
                code = CODE_ITE
                break
            end
            info_str *= @sprintf("    %4d  ", step)

            # Check status of guess
            if !all(isfinite, x_cur)
                code = CODE_NAN
                break
            end
            _set_tmps!(x_cur)

            # Run chemistry scheme
            if chem
                @debug "        chemistry"
                fc_retcode = chemistry.fastchem_eqm!(atmos, false)
                if fc_retcode == 0
                    stepflags *= "Cs-"  # chemistry success
                else
                    stepflags *= "Cf-"  # chemistry failure
                    step_ok = false
                end
            end

            # Switch RT scheme?
            if grey_start & !grey_step
                atmos.rt_scheme = rt_default
            end
            if atmos.rt_scheme == atmosphere.RT_GREYGAS
                stepflags *= "Gg-"
            end

            # Check convective modulation
            easy_step = false
            if easy_start
                # We are modulating
                stepflags *= "M"

                # Check if sf needs to be increased
                if c_cur*easy_trig < conv_atol + conv_rtol * c_max
                    if easy_sf < 1.0
                        # increase sf => reduce modulation by easy_incr
                        stepflags *= "r"

                        # starting from sf=0
                        if easy_sf < 1.0e-10
                            easy_sf = 3e-4
                        end

                        easy_sf = min(1.0, easy_sf*easy_incr)
                        easy_step = true
                        @debug "    updated easy_sf = $easy_sf"

                        # done modulating
                        if easy_sf > 0.99
                            easy_start = false
                        end
                    else
                        easy_sf = 1.0
                    end
                end

                if stepflags[end] == "M"
                    stepflags *= "c"
                end
                stepflags *= "-"
            else
                # No modulation at this point
                easy_sf = 1.0
            end

            # Determine which parts of the Jacobian matrix need to be updated
            @debug "        jacobian"
            if (step <= 2) || perturb_all || easy_step ||
                    (c_cur*perturb_trig < conv_atol + conv_rtol * c_max) ||
                    mod(step,perturb_mod)==0
                # Update whole matrix when any of these are true:
                #    - first step
                #    - it was requested by the user
                #    - we are near global convergence
                fill!(perturb, true)
            else
                # Skip updating Jacobian where the residuals are small, so
                #    that this column of J will be left with the last values
                #    calculated by the finite-difference scheme. This is okay
                #    as long as the jacobian is approx diagonally dominant, or
                #    the layer is near convergence.
                fill!(perturb, true)
                for i in 3:arr_len-4
                    perturb[i] = (sum(abs.(r_cur[i-2:i+2])) > perturb_crit) ||
                                    any(atmos.mask_l[i-2:i+2]) ||
                                    any(atmos.mask_c[i-2:i+2])
                end
            end

            # Evaluate residuals and estimate Jacobian matrix where required
            @. r_old = r_cur
            if fdc || (step == 1)
                # use central difference if:
                #    requested, at the start, or insufficient cost decrease
                if !_calc_jac_res!(x_cur, b, r_cur, true, fdo, perturb)
                    code = CODE_OBJ
                    break
                end
                stepflags *= "C$fdo-"
            else
                # otherwise, use forward difference
                if !_calc_jac_res!(x_cur, b, r_cur, false, fdo, perturb)
                    code = CODE_OBJ
                    break
                end
                stepflags *= "F$fdo-"
            end

            # Check if jacobian is singular
            #    Diagonal elements are usually negative
            #    Off-diagonals are usually positive
            if abs(det(b)) < floatmin()*10.0
                code = CODE_SIN
                step_ok = false
                break
            end

            # Model step
            @. x_old = x_cur
            if (method == 1)
                # Newton-Raphson step
                # @debug "        NR step"
                x_dif = -b\r_cur
                # stepflags *= "Nr-"

            elseif method == 2
                # Gauss-Newton step
                # @debug "        GN step"
                x_dif = -(b'*b) \ (b'*r_cur)
                # stepflags *= "Gn-"

            elseif method == 3
                # Levenberg-Marquardt step
                # @debug "        LM step"
                #    Calculate damping parameter ("delayed gratification")
                if r_cur_2nm < r_old_2nm
                    lml /= 5.0
                else
                    lml *= 2.5
                end

                #    Update our estimate of the solution
                x_dif = -(b'*b + lml * dtd) \ (b' * r_cur)
                # stepflags *= "Lm-"
            end

            # Max step size
            #    limit by user requirement
            dx_max_step = Float64(dx_max)
            #    limit by distance to tmp_ceil
            dx_max_step = min(dx_max_step, abs(atmos.tmp_ceiling-tmp_pad - maximum(x_cur)) )
            #    limit by distance to tmp_floor
            dx_max_step = min(dx_max_step, abs(minimum(x_cur) - atmos.tmp_floor+tmp_pad))

            # Limit step size by dx_max_step, without changing direction of dx vector
            x_dif *= min(1.0, dx_max_step / maximum(abs.(x_dif[:])))

            # Extrapolate step if on plateau.
            #    This acts to give the solver a 'nudge' in (hopefully) the right direction.
            #    Otherwise, this perturbation can still help.
            plateau_apply = (plateau_i > plateau_n)
            if plateau_apply
                @. x_dif *= plateau_s
                plateau_i = 0
                stepflags *= "P-"
            end

            # Linesearch
            # https://people.maths.ox.ac.uk/hauser/hauser_lecture2.pdf
            if linesearch && !plateau_apply
                @debug "        linesearch"

                # Reset
                ls_alpha = ls_max_scale   # Greater than 1 => extension of step
                ls_cost  = 1.0e99   # a big number

                # Internal function minimised by linesearch method
                function _ls_func(scale::Float64)::Float64
                    @. x_cur = x_old + scale * x_dif
                    _fev!(x_cur,r_tst)
                    return _cost(r_tst)
                end

                # Calculate the cost using the full step size
                ls_cost = _ls_func(ls_alpha)

                # Do we need to do linesearch? Triggers due to any of:
                #    - Cost increase from full step is too large
                #    - It is the first step
                if (ls_cost > c_cur*ls_increase ) || (step == 1)

                    # Yes, we do need to do linesearch...
                    stepflags *= "Ls-"

                    if (ls_method == 1) || (ls_cost*0.1 < conv_atol + conv_rtol * c_max)
                        # Use golden-section search method
                        ls_alpha = gs_search(_ls_func, ls_min_scale, ls_max_scale,
                                                1.0e-9, ls_min_scale, ls_max_steps)

                    elseif ls_method == 2
                        # Use backtracking method

                        for il in 1:ls_max_steps
                            # try shrinking it further
                            ls_alpha *= ls_tau

                            if ls_alpha < ls_min_scale
                                # scale too small!
                                ls_alpha = ls_min_scale
                                break
                            end

                            ls_cost = _ls_func(ls_alpha)
                            if ls_cost <= c_cur*ls_increase
                                # this scale is good enough
                                break
                            end
                        end

                    else
                        @error "Invalid linesearch algorithm $ls_method"
                        code = CODE_CFG
                        break
                    end

                    # Apply best scale from linesearch
                    ls_alpha = min(max(ls_alpha, ls_min_scale), ls_max_scale)
                    @debug "    linesearch scale ls_alpha = $ls_alpha"
                    x_dif *= ls_alpha
                end

            end # end linesearch

            # Limit step size by dx_min
            x_dif[x_dif.<0.0] .= clamp.(x_dif[x_dif.<0.0], -dx_max_step, -dx_min)
            x_dif[x_dif.>0.0] .= clamp.(x_dif[x_dif.>0.0], dx_min, dx_max_step)

            # Take the step
            @. x_cur = x_old + x_dif

            # Recalculate layer properties
            if ! atmosphere.calc_layer_props!(atmos)
                code = CODE_HYD
                step_ok = false
                stepflags *= "Ub-"
            end

            # Evaluate fluxes
            _fev!(x_cur, r_cur)

            # New cost value from this step
            c_old = c_cur
            c_cur = _cost(r_cur)

            # If cost ratio is near unity, then the model is struggling
            # to move around because it's on a "plateau" in solution space
            if (plateau_r < c_cur/c_old < 1.0/plateau_r ) && detect_plateau
                plateau_i += 1
            else
                plateau_i = 0
            end

            # Model statistics
            r_med   =   median(r_cur)
            iworst  =   argmax(abs.(r_cur))
            r_max   =   r_cur[iworst]
            x_med   =   median(x_cur)
            x_max   =   x_cur[argmax(abs.(x_cur))]
            dx_stat =   maximum(abs.(x_dif))
            r_old_2nm = r_cur_2nm
            r_cur_2nm = norm(r_cur)
            c_max =     maximum(abs.(atmos.flux_tot))

            # Plot
            if (modplot > 0) && (mod(step, modplot) == 0)
                plot_step()
                if plot_jacobian
                    plotting.jacobian(b, path_jac, perturb=perturb)
                end
            end

            # Inform user
            info_str *= @sprintf("%+.2e  %.3e  %.3e  %.3e  %.3e  %-s",
                                 r_med, c_cur, atmos.flux_u_lw[1],
                                 x_max, dx_stat, stepflags[1:end-1])
            if (modprint>0) && (mod(step, modprint)==0)
                if step_ok
                    @info info_str
                else
                    @warn info_str
                end
            else
                @debug info_str
            end

            # Converged?
            @debug "        check convergence"
            if (c_cur < conv_atol + conv_rtol * c_max) && !easy_start
                # still using grey RT?
                if grey_step
                    # switch to preferred RT scheme
                    grey_step = false
                else
                    # done!
                    code = CODE_SUC
                    break
                end
            end

            # Show benchmark
            if atmos.benchmark
                rt_avg = atmos.tim_rt_eval / atmos.num_rt_eval / 1e9 * 1e3
                @debug @sprintf("Average RT time: %.3f ms", rt_avg)
            end

            # Record that this step not ok
            if (code == 0) && !step_ok
                code = CODE_STP
            end

        end # end solver loop

        # ----------------------------------------------------------
        # Extract solution
        # ----------------------------------------------------------
        atmos.is_solved = true
        atmos.is_converged = false
        if code == CODE_SUC
            @info "    success in $step steps"
            atmos.is_converged = true
            rm(path_plt, force=true)
            rm(path_jac, force=true)
        elseif code == CODE_ITE
            @error "    failure (maximum iterations)"
        elseif code == CODE_SIN
            @error "    failure (singular jacobian)"
            plotting.jacobian(b, path_jac)
        elseif code == CODE_TIM
            @error "    failure (maximum time)"
        elseif code == CODE_NAN
            @error "    failure (NaN values)"
        elseif code == CODE_CFG
            @error "    failure (configuration)"
        elseif code == CODE_OBJ
            @error "    failure (objective function)"
        elseif code == CODE_STP
            @error "    failure (last step not ok)"
        elseif code == CODE_HYD
            @error "    failure (hydrostatic integration)"
        else
            @error "    failure (other)"
        end

        # perform one last evaluation to set `atmos` given the final `x_cur`
        _fev!(x_cur, zeros(Float64, arr_len))

        # calc kzz profile
        if radiative_Kzz
            energy.fill_Kzz!(atmos)
        end

        # calc heating rate profile
        energy.calc_hrates!(atmos)

        # calc LW contribution function
        energy.radtrans!(atmos, true, calc_cf=true)

        # calc diagnostic quantities
        atmosphere.estimate_Ra!(atmos)
        atmosphere.estimate_timescale_conv!(atmos)
        atmosphere.estimate_timescale_rad!(atmos)

        # calc radius of photosphere, and correspondingly the bulk density of the planet
        atmosphere.calc_observed_rho!(atmos)

        # calc ocean scalar quantities
        atmos.ocean_layers = ocean.dist_surf_liq(atmos.cond_total,
                                                    atmos.ocean_ob_frac,
                                                    atmos.ocean_cs_height,
                                                    atmos.rp)
        atmos.ocean_topliq = ocean.get_topliq(atmos.ocean_layers)
        atmos.ocean_maxdepth = ocean.get_maxdepth(atmos.ocean_layers)
        atmos.ocean_areacov = ocean.get_areacov(atmos.ocean_layers, atmos.ocean_ob_frac)


        # ----------------------------------------------------------
        # Print info
        # ----------------------------------------------------------
        loss = maximum(abs.(atmos.flux_tot)) - minimum(abs.(atmos.flux_tot))
        loss_pct = 100.0*loss/maximum(abs.(atmos.flux_tot))
        @info @sprintf("    outgoing LW flux   = %+.2e W m-2     ", atmos.flux_u_lw[1])
        if (sol_type == 2)
            @info @sprintf("    conduct. skin flux = %+.2e W m-2 ", energy.skin_flux(atmos))
        end
        @info @sprintf("    total flux at TOA  = %+.2e W m-2     ", atmos.flux_tot[1])
        @info @sprintf("    total flux at BOA  = %+.2e W m-2     ", atmos.flux_tot[end])
        @info @sprintf("    total flux lost    = %+.2e W m-2  (%+.2e %%) ", loss, loss_pct)
        @info @sprintf("    final cost value   = %+.2e W m-2     ", c_cur)
        @info @sprintf("    surf temperature   = %-9.3f K        ", atmos.tmp_surf)


        return atmos.is_converged
    end # end solve_energy

    """
    **Solve for global (not local) balance with a prescribed atmosphere structure.**

    Comparable to solve_transparent, but with an opaque prescribed atmospheric structure.

    Arguments:
    - `atmos::Atmos_t`         the atmosphere struct instance to be used.
    - `sol_type::Int`          solution types, same as solve_energy
    - `atm_type::Int `         atmosphere prescription (1: isothermal, 2: adiabat, 3: adiabat+stratosphere)
    - `conv_atol::Float64`     convergence: absolute tolerance on global flux [W m-2]
    - `conv_rtol::Float64`     convergence: relative tolerance on global flux [dimensionless]
    - `max_steps::Int`         maximum number of solver steps
    - `tmp_upper::Float64`     upper-bound on Tsurf for golden-section search [K]

    Returns:
    - `Bool` indicating success
    """
    function solve_prescribed!(atmos::atmosphere.Atmos_t;
                                    sol_type::Int=3,
                                    atm_type::Int=1,
                                    conv_atol::Float64=1.0e-3,
                                    conv_rtol::Float64=1.0e-5,
                                    max_steps::Int=300,
                                    tmp_upper::Float64=5000.0)::Bool

        # Validate sol_type (does not allow type=1 here)
        if (sol_type < 1) || (sol_type > 4)
            @error "Invalid solution type ($sol_type)"
            return false
        end
        if sol_type == 1
            @error "Solution type of 1 (fixed Tsurf) is not supported by solve_prescribed"
            return false
        end

        # Validate atm_type
        if (atm_type < 1) || (atm_type > 3)
            @error "Invalid atmosphere prescription ($atm_type)"
            return false
        end


        # Function to set atmosphere according to the desired prescription
        function _prescribe!(atmos::atmosphere.Atmos_t, atm_type::Int, _tsurf::Float64)
            # set tsurf
            atmos.tmp_surf  = deepcopy(_tsurf)
            atmos.tmpl[end] = deepcopy(_tsurf)
            # @debug "    try tmp_surf = $_tsurf K"

            # set profile
            if atm_type == 1
                setpt.isothermal!(atmos, atmos.tmp_surf)
            elseif atm_type == 2
                setpt.dry_adiabat!(atmos)
            elseif atm_type == 3
                setpt.dry_adiabat!(atmos)
                setpt.stratosphere!(atmos, phys.calc_Tskin(atmos.instellation, atmos.albedo_b))
            end
        end

        # Handle different solution types
        if sol_type == 2
            # Conductive boundary layer => find Tsurf based on Tmagma

            function _skinfunc!(_tsurf::Float64)::Float64
                # Cost function, to be minimised. Takes _tsurf and returns the
                # difference between total flux and conductive skin flux

                # Set temperature
                _prescribe!(atmos, atm_type, _tsurf)

                # Residual = radiative flux minus skin flux
                energy.calc_fluxes!(atmos, radiative=true)
                return (atmos.flux_tot[1] - energy.skin_flux(atmos))^2
            end

            # Find solution for T_surf
            tol = conv_atol + conv_rtol * maximum(abs.(atmos.flux_tot))
            T_surf = gs_search(_skinfunc!, atmos.tmp_floor, atmos.tmp_ceiling,
                                    0.0, tol, max_steps; warnings=true)

            # Store final result
            _skinfunc!(T_surf)

        elseif sol_type == 3

            function _intfunc!(_tsurf::Float64)::Float64
                # Cost function, to be minimised. Takes _tsurf and returns the
                # difference between total flux and required total flux

                # Set temperature
                _prescribe!(atmos, atm_type, _tsurf)

                # Residual = radiative flux minus desired flux
                energy.calc_fluxes!(atmos, radiative=true)

                @debug "    flux_tot = $(atmos.flux_tot[1])"
                return (atmos.flux_tot[1] - atmos.flux_int)^2
            end

            # Find solution for T_surf
            tol = conv_atol + conv_rtol * maximum(abs.(atmos.flux_tot))
            T_surf = gs_search(_intfunc!, atmos.tmp_floor, tmp_upper,
                                    0.0, tol, max_steps; warnings=true)

            # Store final result
            _intfunc!(T_surf)

        elseif sol_type == 4

            function _olrfunc!(_tsurf::Float64)::Float64
                # Cost function, to be minimised. Takes _tsurf and returns the
                # difference between calculated OLR and required OLR

                # Set temperature
                _prescribe!(atmos, atm_type, _tsurf)

                # Residual = radiative flux minus desired flux
                energy.calc_fluxes!(atmos, radiative=true)
                return (atmos.flux_u_lw[1] - atmos.target_olr)^2
            end

            # Find solution for T_surf
            tol = conv_atol + conv_rtol * maximum(abs.(atmos.flux_u_lw))
            T_surf = gs_search(_olrfunc!, atmos.tmp_floor, tmp_upper,
                                    0.0, tol, max_steps; warnings=true)

            # Store final result
            _olrfunc!(T_surf)
        end

        # Flag as solved
        atmos.is_converged = true
        atmos.is_solved = true

        # Print info
        @info @sprintf("    outgoing LW flux   = %+.2e W m-2     ", atmos.flux_u_lw[1])
        if (sol_type == 2)
            F_skin = energy.skin_flux(atmos)
            @info @sprintf("    conduct. skin flux = %+.2e W m-2 ", F_skin)
        end
        @info @sprintf("    total flux         = %+.2e W m-2     ", atmos.flux_tot[1])
        @info @sprintf("    surf temperature   = %-9.3f K        ", atmos.tmp_surf)

        return atmos.is_converged
    end # end solve_prescribed

    """
    **Solve for energy balance with a transparent atmosphere.**

    This will use an isothermal temperature profile, and only modify T_surf.

    Arguments:
    - `atmos::Atmos_t`         the atmosphere struct instance to be used.
    - `sol_type::Int`          solution types, same as solve_energy
    - `conv_atol::Float64`     convergence: absolute tolerance on global flux [W m-2]
    - `conv_rtol::Float64`     convergence: relative tolerance on global flux [dimensionless]
    - `max_steps::Int`         maximum number of solver steps
    - `tmp_upper::Float64`     upper-bound on Tsurf for golden-section search [K]

    Returns:
    - `Bool` indicating success
    """
    function solve_transparent!(atmos::atmosphere.Atmos_t;
                                    sol_type::Int=1,
                                    conv_atol::Float64=1.0e-3,
                                    conv_rtol::Float64=1.0e-5,
                                    max_steps::Int=300,
                                    tmp_upper::Float64=5000.0)::Bool


        # Validate sol_type
        if (sol_type < 1) || (sol_type > 4)
            @error "Invalid solution type ($sol_type)"
            return false
        end

        # Check if transparent
        if !atmos.transparent
            @error "Atmosphere is NOT configured to be transparent"
            @error "    Cannot use `solve_transparent`"
            return false
        end

        # Handle different solution types
        if sol_type == 1
            # Fixed temperature case => just calculate radiative fluxes
            energy.calc_fluxes!(atmos, radiative=true)

        elseif sol_type == 2
            # Conductive boundary layer => find Tsurf based on Tmagma

            function _skinfunc!(_tsurf::Float64)::Float64
                # Cost function, to be minimised. Takes _tsurf and returns the
                # difference between total flux and conductive skin flux

                # Set temperature
                atmos.tmp_surf = _tsurf

                # Residual = radiative flux minus skin flux
                energy.calc_fluxes!(atmos, radiative=true)
                return (atmos.flux_tot[end-1] - energy.skin_flux(atmos))^2
            end

            # Find solution for T_surf
            tol = conv_atol + conv_rtol * maximum(abs.(atmos.flux_tot))
            T_surf = gs_search(_skinfunc!, atmos.tmp_floor, atmos.tmp_ceiling,
                                    0.0, tol, max_steps; warnings=true)

            # Store final result
            _skinfunc!(T_surf)

        elseif sol_type == 3

            function _intfunc!(_tsurf::Float64)::Float64
                # Cost function, to be minimised. Takes _tsurf and returns the
                # difference between total flux and required total flux

                # Set temperature
                atmos.tmp_surf = _tsurf

                # Residual = radiative flux minus desired flux
                energy.calc_fluxes!(atmos, radiative=true)
                return (atmos.flux_tot[end-1] - atmos.flux_int)^2
            end

            # Find solution for T_surf
            tol = conv_atol + conv_rtol * maximum(abs.(atmos.flux_tot))
            T_surf = gs_search(_intfunc!, atmos.tmp_floor, tmp_upper,
                                    0.0, tol, max_steps; warnings=true)

            # Store final result
            _intfunc!(T_surf)

        elseif sol_type == 4

            function _olrfunc!(_tsurf::Float64)::Float64
                # Cost function, to be minimised. Takes _tsurf and returns the
                # difference between calculated OLR and required OLR

                # Set temperature
                atmos.tmp_surf = _tsurf

                # Residual = radiative flux minus desired flux
                energy.calc_fluxes!(atmos, radiative=true)
                return (atmos.flux_u_lw[end-1] - atmos.target_olr)^2
            end

            # Find solution for T_surf
            tol = conv_atol + conv_rtol * maximum(abs.(atmos.flux_u_lw))
            T_surf = gs_search(_olrfunc!, atmos.tmp_floor, tmp_upper,
                                    0.0, tol, max_steps; warnings=true)

            # Store final result
            _olrfunc!(T_surf)
        end

        # Flag as solved
        atmos.is_converged = true
        atmos.is_solved = true

        # Print info
        @info @sprintf("    outgoing LW flux   = %+.2e W m-2     ", atmos.flux_u_lw[1])
        if (sol_type == 2)
            F_skin = energy.skin_flux(atmos)
            @info @sprintf("    conduct. skin flux = %+.2e W m-2 ", F_skin)
        end
        @info @sprintf("    total flux         = %+.2e W m-2     ", atmos.flux_tot[1])
        @info @sprintf("    surf temperature   = %-9.3f K        ", atmos.tmp_surf)

        return atmos.is_converged
    end  # end solve_transparent

end
