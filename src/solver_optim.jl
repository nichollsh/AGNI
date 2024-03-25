# Contains code for obtaining energy balance (CVODE method)

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module solver_optim

    using Printf
    using Statistics
    using Revise
    using LinearAlgebra
    using NOMAD

    import atmosphere 
    import phys
    import plotting

    
    """
    **Obtain radiative-convective equilibrium using optimisation.**

    Solves the non-linear system of equations defined by the flux field
    divergence, minimising flux loss across a cell by iterating the temperature
    profile.

    Not compatible with convective adjustment; MLT must be used.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `surf_state::Int=1`               bottom layer temperature, 0: free | 1: fixed | 2: skin | 3: tmp_eff
    - `condensate::String=""`           condensate to model (if empty, no condensates are modelled)
    - `dry_convect::Bool=true`          enable dry convection
    - `sens_heat::Bool=false`           include sensible heating 
    - `max_steps::Int=200`              maximum number of solver steps
    - `atol::Float64=1.0e-2`            maximum residual at convergence
    - `modplot::Int=0`                  iteration frequency at which to make plots
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            surf_state::Int=1,
                            dry_convect::Bool=true, sens_heat::Bool=false,
                            max_steps::Int=200, atol::Float64=1.0e-2
                            )

        # Validate surf_state
        if (surf_state < 0) || (surf_state > 3)
            error("Invalid surface state ($surf_state)")
        end

        # Dimensionality
        arr_len = atmos.nlev_c 
        if (surf_state >= 2)
            arr_len += 1
        end

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
                    atmos.tmp_surf = _x[end]
                end
            end 
            
            return nothing
        end # end set_tmps

        # Objective function to solve for
        function fev(x::Array)

            resid = zeros(Float64, arr_len)

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
                # Conserve fluxes with constant tmp_surf
                resid[1:end] .= atmos.flux_tot[2:end] .- atmos.flux_tot[1:end-1] 
            elseif (surf_state == 2)
                # Conductive boundary layer
                resid[1:end-1] = atmos.flux_tot[2:end] - atmos.flux_tot[1:end-1] 
                resid[end] = atmos.flux_tot[1] - (atmos.tmp_magma - atmos.tmpl[end]) * atmos.skin_k / atmos.skin_d
            elseif (surf_state == 3)
                # Fluxes equal to sigma*tmp_eff^4
                resid[1:end] .= atmos.flux_tot[1:end] .- atmos.flux_eff
            end

            # Check that residuals are real numbers
            if !all(isfinite, resid)
                display(resid)
                error("Residual array contains NaNs and/or Infs")
            end


            return resid
        end  # end fev

        # Cost function 
        function cost(x::Array)
            cst = norm(fev(x))
            # @printf("Cost: %.3e \n", cst)
            return cst
        end 


        # ----------------------------------------------------------
        # Setup initial guess
        # ---------------------------------------------------------- 
        println("Begin optimisation") 

        @printf("    surf   = %d\n", surf_state)
        if (surf_state == 1)
            @printf("    tmp_surf  = %.2f K\n", atmos.tmp_surf)
        elseif (surf_state == 2)
            @printf("    skin_d = %.2f m\n",         atmos.skin_d)
            @printf("    skin_k = %.2f W K-1 m-1\n", atmos.skin_k)
        elseif (surf_state == 3)
            @printf("    tmp_eff   = %.2f K\n",     atmos.tmp_eff)
            @printf("    Fint   = %.2f W m-2\n", atmos.flux_eff)
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
        # Solver
        # ---------------------------------------------------------- 

        function eval_fct(x)
            return (true, true, [cost(x)])
        end
        pb = NomadProblem(  arr_len, # number of inputs of the blackbox
                            1, # number of outputs of the blackbox
                            ["OBJ"], # type of outputs of the blackbox
                            eval_fct;
                            lower_bound=fill(atmos.tmp_floor,arr_len),
                            upper_bound=fill(atmos.tmp_ceiling,arr_len),
                            options=NOMAD.NomadOptions(
                                display_stats= ["TIME", "BBE", "OBJ"],
                                speculative_search=false
                            )
                        )
        result = solve(pb, x_ini) 
        x_cur = result.x_best_feas

        _set_tmps!(x_cur)

        # Check result
        atmos.is_solved = true
        atmos.is_converged = false 
        if result[0]
            println("    success")
            atmos.is_converged = true
        else 
            println("    failure")
        end
        println(" ")

        fev(x_cur)
        atmosphere.calc_hrates!(atmos)

        # ----------------------------------------------------------
        # Print info
        # ---------------------------------------------------------- 
        loss = atmos.flux_tot[1] - atmos.flux_tot[end]
        loss_pct = loss/atmos.flux_tot[1]*100.0
        @printf("    endpoint fluxes \n")
        @printf("    rad_OLR   = %+.2e W m-2     \n", atmos.flux_u_lw[1])
        if (surf_state == 2)
            F_skin = atmos.skin_k / atmos.skin_d * (atmos.tmp_magma - atmos.tmp_surf)
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