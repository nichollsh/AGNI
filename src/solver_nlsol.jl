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
    - `use_linesearch::Bool=false`      use linesearch to ensure global convergence
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            surf_state::Int=1, condensate::String="",
                            dry_convect::Bool=true, sens_heat::Bool=false,
                            max_steps::Int=500, atol::Float64=1.0e-3,
                            use_linesearch::Bool=false
                            )

        # Validate condensation case
        do_condense  = false 
        i_gas        = -1
        if condensate != "" 
            if condensate in atmos.gases
                do_condense = true 
                i_gas = findfirst(==(gas), atmos.gases)
            else 
                error("Invalid condensate ('$condensate')")
            end 
        end 

        # Work arrays 
        pwr = zeros(Float64, atmos.nlev_c)  # power delivery into each cell

        # Objective function to solve for
        function fev!(F,x)

            # Reset values
            atmos.mask_c[:] .= 0.0
            atmos.mask_p[:] .= 0.0
           
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

            # Interpolate temperature to cell-edge values 
            atmosphere.set_tmpl_from_tmp!(atmos, surf_state)

            # Calculate layer properties 
            atmosphere.calc_layer_props!(atmos)
            
            # Reset fluxes
            atmos.flux_tot[:] .= 0.0

            # +Radiation
            atmosphere.radtrans!(atmos, true)
            atmosphere.radtrans!(atmos, false)
            atmos.flux_tot += atmos.flux_n

            # +Dry convection
            if dry_convect
                atmosphere.mlt!(atmos, pmin=1.0)
                atmos.flux_tot += atmos.flux_c
            end

            # +Turbulence
            if sens_heat
                atmosphere.sensible!(atmos)
                atmos.flux_tot[end] += atmos.flux_sens
            end

            # Calculate power into in each cell
            pwr[1:end] = atmos.flux_tot[2:end] - atmos.flux_tot[1:end-1] 

            # +Condensation
            if do_condense
                
                # This handled by requiring that the residuals (i.e. the power
                # being delivered into a cell) can only be positive in regions 
                # where condensation is occuring, because they cannot be allowed 
                # to cool down. All internal production (i.e. negative power 
                # delivery) yields extra condensate, not a change in temp.

                # Check if each level is condensing. 
                for i in 1:atmos.nlev_c

                    x = atmos.layer_x[i,i_gas]
                    if x < 1.0e-10 
                        continue
                    end

                    Tsat = phys.calc_Tdew(gas,atmos.p[i] * x )
                    if atmos.tmp[i] < Tsat
                        pwr[i] = max(pwr[i], 0.0)

                        atmos.mask_p[i] = atmos.mask_decay 
                        atmos.re[i]   = 1.0e-5  # 10 micron droplets
                        atmos.lwm[i]  = 0.8     # 80% of the saturated vapor turns into cloud
                        atmos.clfr[i] = 1.0     # The cloud takes over the entire cell
                    else 
                        atmos.re[i]   = 0.0
                        atmos.lwm[i]  = 0.0
                        atmos.clfr[i] = 0.0
                    end
                end
            end 

            # Check fluxes are real numbers
            if !all(isfinite, atmos.flux_tot)
                display(atmos.flux_tot)
                error("Flux array contains NaNs and/or Infs")
            end

            # Pass residuals outward
            F[:] .= pwr[:]

            return nothing
        end 

        
        # ----------------------------------------------------------
        # Call solver
        # ---------------------------------------------------------- 
        the_ls = Static()
        if use_linesearch 
            the_ls = BackTracking()
            println("NLSolver: begin Newton-Raphson iterations (with backtracking linesearch)") 
        else
            println("NLSolver: begin Newton-Raphson iterations") 
        end 

        atmosphere.set_tmpl_from_tmp!(atmos, surf_state)

        x0 = zeros(Float64, atmos.nlev_c)
        x0[:] .= atmos.tmp[:]

        sol = nlsolve(fev!, x0, method = :newton, linesearch = the_ls, 
                        iterations=max_steps, ftol=atol, show_trace=true)

        # ----------------------------------------------------------
        # Extract solution
        # ---------------------------------------------------------- 

        if !converged(sol)
            @printf("    stopping atmosphere iterations before convergence (maybe try enabling linesearch) \n\n")
        else
            @printf("    convergence criteria met (%d iterations) \n\n", sol.iterations)
        end

        atmos.tmp[:] .= sol.zero[:]
        fev!(zeros(Float64, atmos.nlev_c), atmos.tmp)
        atmosphere.calc_hrates!(atmos)

        # ----------------------------------------------------------
        # Print info
        # ---------------------------------------------------------- 
        loss = atmos.flux_tot[1] - atmos.flux_tot[end]
        loss_pct = loss/atmos.flux_tot[1]*100.0
        @printf("    endpoint fluxes \n")
        @printf("    rad_OLR   = %+.2e W m-2     \n", atmos.flux_u_lw[1])
        @printf("    tot_TOA   = %+.2e W m-2     \n", atmos.flux_tot[1])
        @printf("    tot_BOA   = %+.2e W m-2     \n", atmos.flux_tot[end])
        @printf("    loss      = %+.2e W m-2     \n", loss)
        @printf("    loss      = %+.2e %%        \n", loss_pct)
        @printf("\n")

    end # end solve_energy 

end 