# This file is part of AGNI. License is GPL-3.0: https://www.gnu.org/licenses

"""
**Solves for global (not local) balance with a prescribed atmosphere structure.**

For example, this can be used to find the surface temperature of an atmosphere with a
prescribed temperature profile, such that the outgoing flux matches the incoming flux, or
such that the outgoing longwave radiation matches a target value.
"""
module solve_prescribed

    using Printf
    using LoggingExtras
    using Statistics: median, mean
    using LinearAlgebra

    import ..atmosphere
    import ..energy
    import ..setpt
    import ..phys
    import ..golden: gs_search

    """
    **Solve for global (not local) balance with a prescribed atmosphere structure.**

    Comparable to solve_transparent, but with an opaque prescribed atmospheric structure.

    Arguments:
    - `atmos::Atmos_t`         the atmosphere struct instance to be used.
    - `sol_type::Int64`        solution types, same as solve_energy
    - `atm_type::Int64 `         atmosphere prescription (1: isothermal, 2: adiabat, 3: adiabat+stratosphere)
    - `conv_atol::Float64`     convergence: absolute tolerance on global flux [W m-2]
    - `conv_rtol::Float64`     convergence: relative tolerance on global flux [dimensionless]
    - `max_steps::Int64`       maximum number of solver steps
    - `tmp_upper::Float64`     upper-bound on Tsurf for golden-section search [K]

    Returns:
    - `Bool` indicating success
    """
    function solve_prescribed!(atmos::atmosphere.Atmos_t;
                                    sol_type::Int64=3,
                                    atm_type::Int64=1,
                                    conv_atol::Float64=1.0e-3,
                                    conv_rtol::Float64=1.0e-5,
                                    max_steps::Int64=300,
                                    tmp_upper::Float64=5000.0)::Bool

        # Validate sol_type (does not allow type=1 here)
        if (sol_type < 1) || (sol_type > 4)
            @warn "Invalid solution type ($sol_type)"
            return false
        end
        # Validate atm_type
        if (atm_type < 1) || (atm_type > 3)
            @warn "Invalid atmosphere prescription ($atm_type)"
            return false
        end

        succ::Bool = false

        # Function to set atmosphere according to the desired prescription
        function _prescribe!(atmos::atmosphere.Atmos_t, atm_type::Int64, _tsurf::Float64)
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

            # layer properties
            atmosphere.calc_layer_props!(atmos)
        end

        # Handle different solution types
        if sol_type == 1
            # Constant surface temperature, which is trivial
            _prescribe!(atmos, atm_type, atmos.tmp_surf)
            energy.calc_fluxes!(atmos, radiative=true)
            succ = true

        elseif sol_type == 2
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
            (T_surf,succ) = gs_search(_skinfunc!, atmos.tmp_floor, atmos.tmp_ceiling,
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
            (T_surf,succ) = gs_search(_intfunc!, atmos.tmp_floor, tmp_upper,
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
            (T_surf,succ) = gs_search(_olrfunc!, atmos.tmp_floor, tmp_upper,
                                        0.0, tol, max_steps; warnings=true)

            # Store final result
            _olrfunc!(T_surf)
        end

        # calc LW contribution function
        energy.calc_fluxes!(atmos,  radiative=true,
                                latent_heat=true, convective=true, sens_heat=true,
                                conductive=true, deep=true,
                                calc_cf=true, calc_hr=true)

        # calc diagnostic quantities
        diagnostics.estimate_Ra!(atmos)
        diagnostics.estimate_timescale_conv!(atmos)
        diagnostics.estimate_timescale_rad!(atmos)

        # Flag as solved
        atmos.is_converged = succ
        atmos.is_solved = true

        # Print info
        # @info @sprintf("    outgoing LW flux   = %+.2e W m-2     ", atmos.flux_u_lw[1])
        if (sol_type == 2)
            F_skin = energy.skin_flux(atmos)
            @info @sprintf("    conduct. skin flux = %+.2e W m-2 ", F_skin)
        end
        @info @sprintf("    total flux         = %+.2e W m-2     ", atmos.flux_tot[1])
        @info @sprintf("    surf temperature   = %-9.3f K        ", atmos.tmp_surf)
        @info @sprintf("    surf pressure      = %.1e bar      ", atmos.p_boa/1e5)

        return atmos.is_converged
    end # end solve_prescribed
    export solve_prescribed!

end
