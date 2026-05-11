module solve_globe

    using Printf
    using LoggingExtras

    import ..atmosphere
    import ..multicol
    include("energy.jl"); import .solve_energy: solve_energy!

    """
    **Solve multi-column problem for radiative-convective equilibrium across the globe**

    Repeatedly call `solve_energy!` on each column, solving them as separate systems of
    equations, until all columns are converged to the specified tolerance.

    Arguments:
    - `globe::Globe_t`              globe struct containing the columns to solve
    - `globe_iters::Int`            maximum number of iterations to perform
    - `kwargs...`                   keyword arguments to pass to `solve_energy!`

    Returns:
    - `succ::Bool`                 whether the solve was successful for all columns
    """
    function solve_globe!(globe::multicol.Globe_t, globe_iters::Int; kwargs...)::Bool

        # extract convergence tolerances from kwargs
        globe_atol::Float64 = kwargs[:conv_atol]
        globe_rtol::Float64 = kwargs[:conv_rtol]

        # Track state of globe solver
        succ::Bool = true       # success during this specific iteration
        conv::Bool = false      # globe converged

        # Iterate until convergence, or until giving up
        for iter in 1:globe_iters
            @info "Iteration $iter/$globe_iters of globe solver (ncol=$(globe.ncol))"

            # Assume success during this iteration
            succ = true

            # Update boundary conditions, after first iteration
            if iter>1
                succ &= multicol.set_surface_bc!(globe, kwargs[:sol_type])
                @info "Updated columns' surface boundary conditions"
            end

            # Update heating profiles, after first iteration
            if iter>1
                succ &= multicol.set_redist!(globe)
                @info "Updated columns' heat redistribution profiles"
            end

            # Solve each column independently
            #    This function runs `solve_energy!` on each column independently
            succ &= multicol.call_for_globe!(globe, solve_energy!; kwargs...)

            # Check if any of the solvers failed
            succ || @warn "Global solver partially failed during iteration $iter"

            # Check that globe converged on state where all columns get same flux_tot.
            #   For sol_type=1, system closed by updating heat_redist
            #   For sol_type=2, system closed by updating surface BC and heat_redist
            #   For sol_type=3, system closed by updating heat_redist
            @info "Checking globe for convergence"
            conv = succ

            #   Determine a reasonable target value
            conv_val = median([atmos.flux_tot[1] for atmos in globe.atmos_arr])
            @info @sprintf("    Convergence target = %+.2e W m-2", conv_val)

            #    Check all columns
            for (iatmos,atmos) in enumerate(globe.atmos_arr)
                resid_val = atmos.flux_tot[1]-conv_val
                conv_val  = globe_atol + globe_rtol * abs(conv_val)
                conv &= (abs(resid_val) < conv_val)
                @info @sprintf("    Column %d: resid = %+.2e W m-2 (%+.2f%%)",
                                iatmos, resid_val, resid_val/conv_val*100
                                )
            end

            # Converged?
            conv && break

            @info "-------------------------------"
            @info " "
        end

        if conv
            @info "Globe solve converged"
        else
            @warn "Globe solve did not converge ($globe_iters iterations)"
        end
        return conv
    end

end
