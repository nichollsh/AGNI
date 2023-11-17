# Contains code for obtaining energy balance (CVODE method)

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module solver_cvode

    # Include libraries
    include("../socrates/julia/src/SOCRATES.jl")

    using Printf
    using Statistics
    using Revise
    using PCHIPInterpolation
    using LinearAlgebra

    using SciMLBase
    using SteadyStateDiffEq
    using Sundials

    import atmosphere 
    import phys

    """
    **Obtain radiative-convective equilibrium using a CVODE solver.**

    Time-steps the temperature profile with the heating rates to obtain global and 
    local energy balance. Uses a high-order adaptive stiff ODE integrator from the
    Sundials library (CVode Adams-Moulton) to obtain an accurate solution. 
    
    This is very expensive, so it's best to get close to the solution by using
    the functions in solver_euler.jl, and then using this one afterwards.

    Not compatible with convective adjustment; MLT must be used.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `surf_state::Int=1`               bottom layer temperature, 0: free | 1: fixed | 2: skin
    - `verbose::Bool=false`             verbose output?
    - `dry_convect::Bool=true`          enable dry convection
    - `sens_heat::Bool=false`           include sensible heating 
    - `max_steps::Int=50`               maximum number of solver steps
    """
    function solve_energy!(atmos::atmosphere.Atmos_t;
                            surf_state::Int=1, verbose::Bool=false,
                            dry_convect::Bool=true, sens_heat::Bool=false,
                            max_steps::Int=200
                            )


        bot_orig_e = atmos.tmpl[end]
        call::Int = 0
        modprint::Int=25

        if verbose 
            modprint = 2
        end
        if max_steps < modprint
            modprint = 2
        end

        function objective(du, u, p, t)

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
            clamp!(atmos.tmp, atmos.tmp_floor, Inf)

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

            # ----------------------------------------------------------
            # Calculate heating rate in each layer (should be minimised)
            # ---------------------------------------------------------- 
            atmosphere.calc_hrates!(atmos)
            for i in 1:atmos.nlev_c 
                du[i] = atmos.heating_rate[i]
            end
            
            if mod(call,modprint) == 0 
                println("    Max hr = $(maximum(abs.(du)))")
                println("    Avg hr = $(mean(abs.(du)))")
                println(" ")
            end

            return nothing
        end # end objective function
        
        # ----------------------------------------------------------
        # Call solver
        # ---------------------------------------------------------- 

        println("RCSolver: Begin SUNDIALS CVODE integration")

        u0 = zeros(Float64, atmos.nlev_c)
        u0[:] .= atmos.tmp[:]

        p = 2.0

        prob_ss = SteadyStateProblem(objective, u0, p)
        sol = solve(prob_ss, DynamicSS(CVODE_BDF()), dt=1.0e-3,  abstol=1e-3, reltol=1e-5, maxiters=max_steps)

        # prob_nl = NonlinearProblem(objective, u0, p)
        # sol = solve(prob_nl, LevenbergMarquardt(autodiff=false), dt=200.0, abstol=1e-1, reltol=1e-2, maxiters=max_steps)

        if sol.retcode == :Success
            println("RCSolver: Iterations completed (converged)")
        else
            println("RCSolver: Iterations completed (maximum iterations or failure)")
        end

        # ----------------------------------------------------------
        # Extract solution
        # ---------------------------------------------------------- 
        atmos.tmp[:] .= sol.u[:]
        objective(zeros(Float64, atmos.nlev_c), atmos.tmp, p, 2.0e20)

        # ----------------------------------------------------------
        # Print info
        # ---------------------------------------------------------- 
        loss = atmos.flux_tot[1] - atmos.flux_tot[end]
        loss_pct = loss/atmos.flux_tot[1]*100.0
        @printf("RCSolver: Final total fluxes [W m-2] \n")
        @printf("    rad_OLR   = %.2e W m-2         \n", atmos.flux_u_lw[1])
        @printf("    tot_TOA   = %.2e W m-2         \n", atmos.flux_tot[1])
        @printf("    tot_BOA   = %.2e W m-2         \n", atmos.flux_tot[end])
        @printf("    loss      = %.2f W m-2         \n", loss)
        @printf("    loss      = %.2f %%            \n", loss_pct)
        @printf("\n")

    end # end solve_root
end 