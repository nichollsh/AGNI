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

    using SciMLBase
    using NonlinearSolve

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
                            max_steps::Int=500
                            )


        bot_orig_e = atmos.tmpl[end]
        call::Int = 0
        modprint::Int=25

        if verbose 
            modprint = 10
        end
        if max_steps < modprint
            modprint = 2
        end

        function objective(du, u, p)

            call += 1
            if mod(call,modprint) == 0 
                println("    call $call")
            end
            

            # ----------------------------------------------------------
            # Set atmosphere
            # ---------------------------------------------------------- 
            for i in 1:atmos.nlev_c
                atmos.tmp[i] = u[i]
            end
            clamp!(atmos.tmp, atmos.tmp_floor, atmos.tmp_ceiling)

            # Interpolate temperature to cell-edge values 
            atmosphere.set_tmpl_from_tmp!(atmos, surf_state)
            
            # ----------------------------------------------------------
            # Get fluxes 
            # ---------------------------------------------------------- 
            atmos.flux_tot[:] .= 0.0

            # Radiation
            atmosphere.radtrans!(atmos, true)
            atmosphere.radtrans!(atmos, false)
            atmos.flux_tot += atmos.flux_n

            # Convection
            if dry_convect
                atmosphere.mlt!(atmos)
                atmos.flux_tot += atmos.flux_c
            end

            # Turbulence
            if sens_heat
                atmosphere.sensible!(atmos)
                atmos.flux_tot[end] += atmos.flux_sens
            end

            for i in 1:atmos.nlev_c 
                du[i] = atmos.flux_tot[i] - atmos.flux_tot[i+1]
            end
            
            if mod(call,modprint) == 0 
                println("    Max df = $(maximum(abs.(du)))")
                println("    Avg df = $(mean(abs.(du)))")
                println(" ")
            end

            return nothing
        end # end objective function
        
        # ----------------------------------------------------------
        # Call solver
        # ---------------------------------------------------------- 

        println("NLSolver: Begin chi-squared minimisation")

        u0 = zeros(Float64, atmos.nlev_c)
        u0[:] .= atmos.tmp[:]

        p = 2.0

        # prob_ss = SteadyStateProblem(objective, u0, p)
        # sol = solve(prob_ss, DynamicSS(CVODE_BDF()), dt=1.0e-3,  abstol=1e-3, reltol=1e-5, maxiters=max_steps)

        prob_nl = NonlinearProblem(objective, u0, p)
        sol = solve(prob_nl, NewtonRaphson(autodiff=false), dt=50.0, abstol=1e-1, reltol=1e-3, maxiters=max_steps)

        if sol.retcode == :Success
            println("NLSolver: Iterations completed (converged)")
        else
            println("NLSolver: Iterations completed (maximum iterations or failure)")
        end

        # ----------------------------------------------------------
        # Extract solution
        # ---------------------------------------------------------- 
        atmos.tmp[:] .= sol.u[:]
        objective(zeros(Float64, atmos.nlev_c), atmos.tmp, p)

        # ----------------------------------------------------------
        # Print info
        # ---------------------------------------------------------- 
        loss = atmos.flux_tot[1] - atmos.flux_tot[end]
        loss_pct = loss/atmos.flux_tot[1]*100.0
        @printf("NLSolver: Final total fluxes [W m-2] \n")
        @printf("    rad_OLR   = %.2e W m-2         \n", atmos.flux_u_lw[1])
        @printf("    tot_TOA   = %.2e W m-2         \n", atmos.flux_tot[1])
        @printf("    tot_BOA   = %.2e W m-2         \n", atmos.flux_tot[end])
        @printf("    loss      = %.4f W m-2         \n", loss)
        @printf("    loss      = %.4f %%            \n", loss_pct)
        @printf("\n")

    end # end solve_root
end 