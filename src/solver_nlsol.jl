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

    using NLsolve
    using LineSearches

    import atmosphere 
    import phys
    
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
    - `max_steps::Int=500`              maximum number of solver steps
    - `atol::Int=1.0e-3`                maximum residual at convergence
    - `use_linesearch::Bool=false`      use linesearch to ensure global convergence
    - `calc_cf_end::Bool=true`          calculate contribution function at convergence
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            surf_state::Int=1, condensate::String="",
                            dry_convect::Bool=true, sens_heat::Bool=false,
                            max_steps::Int=500, atol::Float64=1.0e-3,
                            use_linesearch::Bool=false,
                            calc_cf_end::Bool=true
                            )

        # Validate condensation case
        do_condense  = false 
        i_gas        = -1
        if condensate != "" 
            if condensate in atmos.gases
                do_condense = true 
                i_gas = findfirst(==(condensate), atmos.gases)
            else 
                error("Invalid condensate ('$condensate')")
            end 
        end 

        # Work arrays 
        arr_len = atmos.nlev_c 
        if (surf_state == 2)
            arr_len += 1
        end
        resid = zeros(Float64, arr_len)  # residuals
        prate = zeros(Float64, atmos.nlev_c)  # condensation production rate [kg /m3 /s]
        calc_cf = false

        # Objective function to solve for
        function fev!(F,x)

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

            # Pass residuals outward
            F[:] .= resid[:]

            return nothing
        end 

        
        # ----------------------------------------------------------
        # Call the solver
        # ---------------------------------------------------------- 
        the_ls = Static()
        if use_linesearch 
            the_ls = BackTracking()
            println("NLSolver: begin Newton-Raphson iterations (with backtracking linesearch)") 
        else
            println("NLSolver: begin Newton-Raphson iterations") 
        end 

        # Allocate x array
        # Array storage structure:
        #   in the surf_state=2 case
        #       1:end-1 => cell centre temperatures 
        #       end     => bottom cell edge temperature
        #   other cases 
        #       1:end => cell centre temperatures
        x0 = zeros(Float64, arr_len) 
        for i in 1:atmos.nlev_c
            x0[i] = atmos.tmp[i]
        end 
        if (surf_state==2)
            x0[end] = x0[atmos.nlev_c]
        end

        # Start nonlinear solver
        sol = nlsolve(fev!, x0, method = :newton, linesearch = the_ls, 
                        iterations=max_steps, ftol=atol, show_trace=true)

        # ----------------------------------------------------------
        # Extract solution
        # ---------------------------------------------------------- 

        atmos.is_solved = true
        calc_cf = calc_cf_end

        if !converged(sol)
            @printf("    stopping atmosphere iterations before convergence (maybe try enabling linesearch) \n\n")
            atmos.is_converged = false 
        else
            @printf("    convergence criteria met (%d iterations) \n\n", sol.iterations)
            atmos.is_converged = true
        end

        final_x = zeros(Float64, arr_len)
        for i in 1:arr_len
            final_x[i] = sol.zero[i]
        end 
        fev!(zeros(Float64, arr_len), final_x)
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

    end # end solve_energy 

end 