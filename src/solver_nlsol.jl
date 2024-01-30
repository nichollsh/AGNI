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
    - `surf_state::Int=1`               bottom layer temperature, 0: free | 1: fixed | 2: skin
    - `condensate::String=""`           condensate to model (if empty, no condensates are modelled)
    - `dry_convect::Bool=true`          enable dry convection
    - `sens_heat::Bool=false`           include sensible heating 
    - `max_steps::Int=200`              maximum number of solver steps
    - `atol::Float64=1.0e-2`            maximum residual at convergence
    - `cdw::Float64=5.0e-3`             relative width of the "difference" in the central-difference calculations
    - `calc_cf_end::Bool=true`          calculate contribution function at convergence
    - `modplot::Int=0`                  iteration frequency at which to make plots
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            surf_state::Int=1, condensate::String="",
                            dry_convect::Bool=true, sens_heat::Bool=false,
                            max_steps::Int=200, atol::Float64=1.0e-2, cdw::Float64=5.0e-3,
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

        # Work arrays 
        calc_cf = false
        prate = zeros(Float64, atmos.nlev_c)  # condensation production rate [kg /m3 /s]

        # Objective function to solve for
        function fev!(x::Array,resid::Array)

            # Reset values
            atmos.mask_c[:] .= 0
            atmos.mask_p[:] .= 0
           
            # Read new guess
            for i in 1:atmos.nlev_c
                atmos.tmp[i] = x[i]
            end

            # Check that guess is finite (no NaNs, no Infs)
            if !all(isfinite, atmos.tmp)
                display(atmos.tmp)
                error("Temperature array contains NaNs and/or Infs")
            end 

            # Limit temperature domain
            clamp!(atmos.tmp, atmos.tmp_floor, atmos.tmp_ceiling)

            # Interpolate temperature to cell-edge values (not incl. bottommost value)
            atmosphere.set_tmpl_from_tmp!(atmos)

            # Set bottom edge temperature 
            if (surf_state < 0) || (surf_state > 2)
                # Error cases
                error("Invalid surface state ($surf_state)")

            elseif (surf_state != 1)
                # Extrapolate (log-linear)
                grad_dt = atmos.tmp[end]-atmos.tmp[end-1]
                grad_dp = log(atmos.p[end]/atmos.p[end-1])
                atmos.tmpl[end] = atmos.tmp[end] + grad_dt/grad_dp * log(atmos.pl[end]/atmos.p[end])

                if (surf_state == 2)
                    # Conductive skin
                    atmos.tstar = x[end]
                end
            end 
            atmos.tmpl[end] = clamp(atmos.tmpl[end], atmos.tmp_floor, atmos.tmp_ceiling)

            # Calculate layer properties 
            atmosphere.calc_layer_props!(atmos)
            
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

            # Calculate residuals
            if (surf_state == 2)
                resid[1:end-1] = atmos.flux_tot[2:end] - atmos.flux_tot[1:end-1] 
                resid[end] = atmos.flux_tot[1] - (atmos.tmp_magma - atmos.tmpl[end]) * atmos.skin_k / atmos.skin_d
            else 
                resid[1:end] = atmos.flux_tot[2:end] - atmos.flux_tot[1:end-1] 
            end

            # +Condensation
            if do_condense

                a = 1.0
                g = 1.0
                f = 0.0

                for i in 1:atmos.nlev_c

                    x = atmos.layer_x[i,i_gas]
                    if x < 1.0e-10 
                        continue
                    end

                    # Layer is condensing if T < T_dew
                    Tsat = phys.calc_Tdew(condensate,atmos.p[i] * x )
                    g = 1.0 - exp(a * (Tsat - atmos.tmp[i]))
                    f = resid[i]
                    resid[i] = f * g
                    if atmos.tmp[i] <= Tsat+0.1
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

        # Calculating the jacobian of f at x using a central-difference method
        function jac!(x::Array, jacob::Array)

            # Reset jacobian 
            fill!(jacob, 0.0)

            len_x = length(x)

            # Work variables 
            tmp_s   = zeros(Float64, len_x)  # Perturbed temperature array 
            s       = 0.0                           # Row perturbation
            rf      = zeros(Float64, len_x)  # Forward difference
            rb      = zeros(Float64, len_x)  # Backward difference
            drdt    = zeros(Float64, len_x)  # Jacobian row

            for i in 1:len_x  # for each x

                # Calculate perturbation
                s = x[i] * cdw

                # Forward
                tmp_s[:] .= x[:]
                tmp_s[i] += s
                fev!(tmp_s, rf)

                # Backward
                tmp_s[:] .= x[:]
                tmp_s[i] -= s
                fev!(tmp_s, rb)

                # Central difference
                drdt[:] .= ( (rf[:] .- rb[:]) ./ (2.0 * s) )

                # Set jacobian
                for j in 1:len_x  # for each r
                    jacob[j,i] = drdt[j]
                end 
            end 

            return nothing
        end # end jac

        
        # ----------------------------------------------------------
        # Setup initial guess
        # ---------------------------------------------------------- 

        println("NLSolve: begin Newton-Raphson iterations") 

        arr_len = atmos.nlev_c 
        if (surf_state == 2)
            arr_len += 1
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
        if (surf_state==2)
            x_ini[end] = atmos.tmp[atmos.nlev_c] + 1.0
        end

        # ----------------------------------------------------------
        # Solver loop
        # ---------------------------------------------------------- 
        step::Int =     0       # Step number
        code::Int =     -1      # Status code 
        modprint::Int = 1       # Print frequency
        solving::Bool = true    # Currently solving?
        ssf::Float64  = 1.0     # Step scale factor
        

        b::Array      = zeros(Float64, (arr_len, arr_len))    # Approximate jacobian (i)
        x_cur::Array  = zeros(Float64, arr_len)               # Array to solve for (i)
        x_old::Array  = zeros(Float64, arr_len)               # Old ^ (i-1)
        x_dif::Array  = zeros(Float64, arr_len)               # Change in x (i-1 to i)
        r::Array      = zeros(Float64, arr_len)               # Residuals (i)

        x_cur[:] .= x_ini[:]
        r .+= 1.0e99

        @printf("    step  resid_med  resid_max  xvals_med  xvals_max  \n")
        while solving 
            # Check convergence
            if maximum(abs.(r)) < atol 
                code = 0
                solving = false 
                break
            end

            step += 1
            if step > max_steps
                code = 1 
                solving = false
                break 
            end 
            if mod(step,modprint) == 0 
                @printf("    %4d", step)
            end

            # Evaluate F and J
            fev!(x_cur, r)
            jac!(x_cur, b) 

            # Check if jacobian is singular 
            # if det(b) < 1.0e-99
            #     code = 2
            #     solving = false
            #     break
            # end 

            # Newton-Raphson step 
            x_dif = -b\r
            x_old[:] = x_cur[:]
            x_cur[:] .= x_old[:] .+ (x_dif[:] .* ssf)

            # Plot?
            if modplot > 0
                if mod(step, modplot) == 0
                    plotting.plot_pt(atmos,  @sprintf("%s/solver.png", atmos.OUT_DIR), incl_magma=(surf_state==2))
                end 
            end 
                
            # Inform user
            if mod(step,modprint) == 0 
                r_med = median(r)
                r_max = r[argmax(abs.(r))]
                x_med = median(x_cur)
                x_max = x_cur[argmax(abs.(x_cur))]
                @printf("  %+.2e  %+.2e  %+.2e  %+.2e  \n", r_med, r_max, x_med, x_max)
            end

        end # end solver loop
        
        # Check result
        atmos.is_solved = true
        atmos.is_converged = false 
        if code == 0
            println("    success")
            atmos.is_converged = true
        elseif code == 1
            println("    failure (maximum iterations reached)")
        elseif code == 2
            println("    failure (singular jacobian)")
        else 
            println("    failure (unhandled)")
        end

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