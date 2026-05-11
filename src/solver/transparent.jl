module solve_transparent

    using Printf
    using LoggingExtras

    import ..atmosphere
    import ..energy
    import ..golden: gs_search

    """
    **Solve for energy balance with a transparent atmosphere.**

    This will use an isothermal temperature profile, and only modify T_surf.

    Arguments:
    - `atmos::Atmos_t`         the atmosphere struct instance to be used.
    - `sol_type::Int64`        solution types, same as solve_energy
    - `conv_atol::Float64`     convergence: absolute tolerance on global flux [W m-2]
    - `conv_rtol::Float64`     convergence: relative tolerance on global flux [dimensionless]
    - `max_steps::Int64`       maximum number of solver steps
    - `tmp_upper::Float64`     upper-bound on Tsurf for golden-section search [K]

    Returns:
    - `Bool` indicating success
    """
    function solve_transparent!(atmos::atmosphere.Atmos_t;
                                    sol_type::Int64=1,
                                    conv_atol::Float64=1.0e-3,
                                    conv_rtol::Float64=1.0e-5,
                                    max_steps::Int64=300,
                                    tmp_upper::Float64=5000.0)::Bool


        # Validate sol_type
        if (sol_type < 1) || (sol_type > 4)
            @warn "Invalid solution type ($sol_type)"
            return false
        end

        # Check if transparent
        if !atmos.transparent
            @warn "Atmosphere is NOT configured to be transparent"
            @warn "    Cannot use `solve_transparent`"
            return false
        end

        succ::Bool = false

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
                return (atmos.flux_tot[1] - energy.skin_flux(atmos))^2
            end

            # Find solution for T_surf
            tol = conv_atol + conv_rtol * maximum(abs.(atmos.flux_tot))
            (T_surf,succ) = gs_search(_skinfunc!, atmos.tmp_floor, atmos.tmp_ceiling,
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
                return (atmos.flux_tot[1] - atmos.flux_int)^2
            end

            # Find solution for T_surf
            tol = conv_atol + conv_rtol * maximum(abs.(atmos.flux_tot))
            (T_surf,succ) = gs_search(_intfunc!, atmos.tmp_floor, tmp_upper,
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
                return (atmos.flux_u_lw[1] - atmos.target_olr)^2
            end

            # Find solution for T_surf
            tol = conv_atol + conv_rtol * maximum(abs.(atmos.flux_u_lw))
            (T_surf,succ) = gs_search(_olrfunc!, atmos.tmp_floor, tmp_upper,
                                        0.0, tol, max_steps; warnings=true)

            # Store final result
            _olrfunc!(T_surf)
        end

        # Flag as solved
        atmos.is_converged = succ
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
    export solve_transparent!


end
