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
    using PCHIPInterpolation
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
    - `verbose::Bool=false`             verbose output?
    - `dry_convect::Bool=true`          enable dry convection
    - `sens_heat::Bool=false`           include sensible heating 
    - `max_steps::Int=500`              maximum number of solver steps
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            surf_state::Int=1, verbose::Bool=false,
                            dry_convect::Bool=true, sens_heat::Bool=false,
                            max_steps::Int=500, atol::Float64=1.0e-4
                            )

        modprint = 10
        call = 0

        # Objective function
        function fev!(F,x)

            call += 1
           
            # ----------------------------------------------------------
            # Set atmosphere
            # ----------------------------------------------------------
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
            
            atmos.flux_tot[:] .= 0.0

            # Radiation
            atmosphere.radtrans!(atmos, true)
            atmosphere.radtrans!(atmos, false)
            atmos.flux_tot += atmos.flux_n

            # Convection
            if dry_convect
                atmosphere.mlt!(atmos, pmin=100.0)
                atmos.flux_tot += atmos.flux_c
            end

            # Turbulence
            if sens_heat
                atmosphere.sensible!(atmos)
                atmos.flux_tot[end] += atmos.flux_sens
            end

            # Check fluxes are real numbers
            if !all(isfinite, atmos.flux_tot)
                display(atmos.flux_tot)
                error("Flux array contains NaNs")
            end

            # Pass outward
            for i in 1:atmos.nlev_c 
                F[i] = atmos.flux_tot[i] - atmos.flux_tot[i+1]
            end

            if mod(call, modprint) == 0
                F_OLR_rad = atmos.flux_u_lw[1]
            
                F_TOA_tot = atmos.flux_tot[1]
                F_BOA_tot = atmos.flux_tot[end-1]
                F_loss    = abs( F_TOA_tot-F_BOA_tot )

                r_max = maximum(F) 
                r_med = median(F)
                
                @printf("    call %d \n",  call )
                @printf("    F_rad^OLR   = %+.2e W m-2  \n", F_OLR_rad)
                @printf("    F_tot^TOA   = %+.2e W m-2  \n", F_TOA_tot)
                @printf("    F_tot^BOA   = %+.2e W m-2  \n", F_BOA_tot)
                @printf("    F_tot^loss  = %+.2e W m-2  \n", F_loss)
                @printf("    r_max       = %+.2e W m-2  \n", r_max)
                @printf("    r_med       = %+.2e W m-2  \n", r_med)
                println(" ")
            end

            return nothing
        end 

        
        # ----------------------------------------------------------
        # Call solver
        # ---------------------------------------------------------- 

        println("NLSolver: begin Newton-Raphson-Backtracking iterations")


        atmosphere.smooth_centres!(atmos, 5)
        atmosphere.set_tmpl_from_tmp!(atmos, surf_state)

        x0 = zeros(Float64, atmos.nlev_c)
        x0[:] .= atmos.tmp[:]

        sol = nlsolve(fev!, x0, method = :newton, linesearch = BackTracking(), iterations=max_steps, ftol=atol)

        # ----------------------------------------------------------
        # Extract solution
        # ---------------------------------------------------------- 

        if !converged(sol)
            @printf("    stopping atmosphere iterations before convergence \n\n")
        else
            @printf("    convergence criteria met (%d iterations) \n\n", sol.iterations)
        end

        atmos.tmp[:] .= sol.zero[:]
        fev!(zeros(Float64, atmos.nlev_c), atmos.tmp)

        # ----------------------------------------------------------
        # Print info
        # ---------------------------------------------------------- 
        loss = atmos.flux_tot[1] - atmos.flux_tot[end]
        loss_pct = loss/atmos.flux_tot[1]*100.0
        @printf("    endpoint fluxes \n")
        @printf("    rad_OLR   = %.2e W m-2         \n", atmos.flux_u_lw[1])
        @printf("    tot_TOA   = %.2e W m-2         \n", atmos.flux_tot[1])
        @printf("    tot_BOA   = %.2e W m-2         \n", atmos.flux_tot[end])
        @printf("    loss      = %.2e W m-2         \n", loss)
        @printf("    loss      = %.2e %%            \n", loss_pct)
        @printf("\n")

    end # end solve_root
end 