module layers

    import ..atmosphere
    import ..phys
    import ..density
    import ..species

    # Hydrostatic+gravity+mass calculation (constants and limits)
    HYDROGRAV_steps::Int64   = 2000      # total number of steps in height integration
    HYDROGRAV_maxdr::Float64 = 1e8       # maximum dz across each layer [m]
    HYDROGRAV_mindr::Float64 = 1e-5      # minimum dz across each layer [m]
    HYDROGRAV_ming::Float64  = 1e-4      # minimum allowed gravity [m/s^2]
    HYDROGRAV_constg::Bool   = false     # constant gravity with height?
    HYDROGRAV_selfg::Bool    = true      # include self-gravity of the atmosphere?

    """
    **Calculate properties within each layer of the atmosphere (e.g. density, mmw).**

    Function will return false if hydrostatic calculcation fails. This is usually when
    the atmosphere becomes unbound.

    Arguments:
        - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
        - `ok::Bool`                function result is ok
    """
    function calc_layer_props!(atmos::atmosphere.Atmos_t)::Bool
        if !atmos.is_param
            @error("Atmosphere struct has not been setup")
            return false
        end

        # Status
        ok::Bool = true

        # MMW
        calc_profile_mmw!(atmos)

        # Heat capacity and thermal conductivity
        calc_profile_cpkc!(atmos)

        # Set density at each level
        calc_profile_density!(atmos)

        # Perform hydrostatic integration
        ok = ok && calc_profile_radius!(atmos)

        # Calculate scale height at each layer
        calc_profile_Hp!(atmos)

        return ok
    end

    """
    **Calculate radii, gravities, and masses for all layers.**

    Performs hydrostatic integration from the ground upwards.
    Requires density, temperature, pressure to have already been set.

    Does not account for surface ocean height.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
    - `bound::Bool`             atmosphere is strongly bound by gravity
    """
    function calc_profile_radius!(atmos::atmosphere.Atmos_t)::Bool

        # Calculate net surface acceleration [m s-2]
        a_surf::Float64 = atmos.grav_surf -
                            phys.cent_accel(atmos.axial_period, atmos.rp, atmos.col_lat)

        # Reset arrays
        fill!(atmos.r         ,    atmos.rp)
        fill!(atmos.rl        ,    atmos.rp)
        fill!(atmos.g         ,    atmos.grav_surf)
        fill!(atmos.gl        ,    atmos.grav_surf)
        fill!(atmos.a         ,    a_surf)
        fill!(atmos.al        ,    a_surf)
        fill!(atmos.m         ,    atmos.interior_mass)
        fill!(atmos.ml        ,    atmos.interior_mass)
        fill!(atmos.layer_thick,   1.0)
        fill!(atmos.layer_σ ,      1.0)
        fill!(atmos.layer_isbound, true)

        # Check config...
        if HYDROGRAV_constg && HYDROGRAV_selfg
            @warn "Incompatible gravity parameters have been set:"
            @warn "    constant with height (HYDROGRAV_constg=$HYDROGRAV_constg)"
            @warn "    atmos self-attraction (HYDROGRAV_selfg=$HYDROGRAV_selfg)"
        end

        # Temporary values
        nsub::Int64 = round(Int64, HYDROGRAV_steps/atmos.nlev_c, RoundUp)

        # Integrate from surface upwards
        for i in range(start=atmos.nlev_c, stop=1, step=-1)

            # ------------
            # Integrate from lower edge to centre
            atmos.r[i], atmos.g[i], atmos.m[i] =
                integ_hydrograv(atmos.rl[i+1],
                                    atmos.gl[i+1], atmos.al[i+1],
                                    atmos.ml[i+1], atmos.pl[i+1],
                                    atmos.p[i], atmos.layer_ρ[i], nsub)

            #   apply radius limiter
            atmos.r[i] = max(atmos.r[i], atmos.rl[i+1] + HYDROGRAV_mindr)
            if atmos.r[i] > atmos.rl[i+1] + HYDROGRAV_maxdr/2
                atmos.r[i] = atmos.rl[i+1] + HYDROGRAV_maxdr/2
                atmos.layer_isbound[i] = false
            end

            #   apply gravity limiter
            if HYDROGRAV_constg
                atmos.g[i]  = atmos.grav_surf
            end
            if atmos.g[i] < HYDROGRAV_ming
                atmos.g[i] = HYDROGRAV_ming
                atmos.layer_isbound[i] = false
            end

            # calculate net acceleration at layer centre
            atmos.a[i] = atmos.g[i] -
                            phys.cent_accel(atmos.axial_period, atmos.r[i], atmos.col_lat)

            # ------------
            # Integrate from centre to upper edge
            atmos.rl[i], atmos.gl[i], atmos.ml[i] =
                integ_hydrograv(atmos.r[i],
                                    atmos.g[i], atmos.a[i],
                                    atmos.m[i], atmos.p[i],
                                    atmos.pl[i], atmos.layer_ρ[i], nsub)

            #   apply radius limiter
            atmos.rl[i] = max(atmos.rl[i], atmos.r[i] + HYDROGRAV_mindr)
            if atmos.rl[i] > atmos.r[i] + HYDROGRAV_maxdr/2
                atmos.rl[i] = atmos.r[i] + HYDROGRAV_maxdr/2
                atmos.layer_isbound[i] = false
            end

            #   apply gravity limiter
            if HYDROGRAV_constg
                atmos.gl[i] = atmos.grav_surf
            end
            if atmos.gl[i] < HYDROGRAV_ming
                atmos.gl[i] = HYDROGRAV_ming
                atmos.layer_isbound[i] = false
            end

            #  calculate net acceleration at layer upper edge
            atmos.al[i] = atmos.gl[i] -
                            phys.cent_accel(atmos.axial_period, atmos.rl[i], atmos.col_lat)

            # Store: Layer geometrical thickness [m]
            atmos.layer_thick[i] = atmos.rl[i] - atmos.rl[i+1]

            # Mass of layer, per unit area at layer-centre [kg m-2]
            atmos.layer_σ[i] = (atmos.ml[i] - atmos.ml[i+1])/(4 * pi * atmos.r[i]^2)
        end

        return all(atmos.layer_isbound)
    end

    """
    **Integrate hydrostatic and gravity equations across a pressure interval.**

    Uses the classic fourth-order Runge-Kutta method.

    Internally, this function integrates from p0 to p1, which means that the pressure is
    decreasing across the interval. The gravity is updated at each integration step,
    and is offset by the centrifugal acceleration.

    Arguments:
    - `r0::Float64`     radius   at start of interval [m]
    - `g0::Float64`     gravity  at start of interval [m s-2]
    - `a0::Float64`     net accel at start of interval (grav - cent) [m s-2]
    - `m0::Float64`     mass enc at start of interval [kg]
    - `p0::Float64`     pressure at start of interval [Pa]
    - `p1::Float64`     pressure at end   of interval [Pa]
    - `rho::Float64`    density throughout interval, constant [kg m-3]
    - `n::Int64`        number of steps for integration (n >= 2)

    Returns:
    - `rj::Float64`     radius   at end of interval [m]
    - `gj::Float64`     gravity  at end of interval [kg]
    - `mj::Float64`     mass enc at end of interval [kg]
    """
    function integ_hydrograv(r0::Float64, g0::Float64, a0::Float64, m0::Float64, p0::Float64,
                                    p1::Float64, rho::Float64, n::Int64)::Tuple{Float64,Float64,Float64}

        # Work variables
        pj::Float64 = p0    # rolling pressure (decreasing)
        rj::Float64 = r0    # rolling radius   (increasing)
        gj::Float64 = g0    # rolling gravity  (incr, decr, or constant)
        mj::Float64 = m0    # rolling mass     (increasing)

        # Get centripetal acceleration's rotation factor, (2 pi cosθ / period)^2
        # Backing-out the factor this way avoids extra trig function calls.
        cent_fact::Float64 = (g0 - a0)/r0  # units: s-2

        # Get gravitational acceleration at r
        function _grav(r)
            if HYDROGRAV_constg
                # gravity is constant
                return g0
            else
                if HYDROGRAV_selfg
                    # gravity changes with mass and radius
                    # approximate mass as constant for this part of the integration
                    return phys.grav_accel(mj, r)
                else
                    # gravity changes with radius only
                    return g0 * (r0/r)^2
                end
            end
        end

        # Get net acceleration at r
        function _accel(r)
            return _grav(r) - r*cent_fact
        end

        # Derivative to integrate
        #   dr/dp = -1 / (rho * a(r))
        function _drdp(p,r)
            return -1 / (rho * _accel(r))
        end

        # Parameters
        dp::Float64  = (p1-p0)/max(2,n) # this will be negative
        dp2::Float64 = dp/2
        k1::Float64  = 0.0; k2::Float64 = 0.0
        k3::Float64  = 0.0; k4::Float64 = 0.0

        # Loop over sub-levels between p0 and p1
        for _ in range(p0, stop=p1, step=dp)

            # Integrate radius ...
            k1 = _drdp(pj,       rj)
            k2 = _drdp(pj + dp2, rj + k1*dp2)
            k3 = _drdp(pj + dp2, rj + k2*dp2)
            k4 = _drdp(pj + dp,  rj + k3*dp)
            rj += dp/6 * (k1 + 2*k2 + 2*k3 + k4)

            # Integrate mass enclosed ...
            mj += 4 * pi * rj^2 * (-1 * dp) / _accel(rj)

            # Integrate pressure (negative change )
            pj += dp

            # Update gravity
            gj = _grav(rj)
        end

        return (rj, gj, mj)
    end

    """
    **Calculate mean molecular weight for all layers.**

    MMW stored as kg/mol.

    Arguments:
        - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    """
    function calc_profile_mmw!(atmos::atmosphere.Atmos_t)
        fill!(atmos.layer_μ, 0.0)
        for gas in atmos.gas_names
            @. atmos.layer_μ += atmos.gas_vmr[gas] * atmos.gas_dat[gas].mmw
        end
        return nothing
    end

    """
    **Calculate specific heat capacity and thermal conductivity for all layers.**

    Specific heat per unit mass: J K-1 kg-1.
    Thermal conductivity: W m-1 K-1.

    Arguments:
        - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    """
    function calc_profile_cpkc!(atmos::atmosphere.Atmos_t)
        # Loop over layers
        @inbounds for i in 1:atmos.nlev_c
            calc_single_cpkc!(atmos, i)
        end
        return nothing
    end

    """
    **Calculate specific heat capacity and thermal conductivity of a single layer.**

    Specific heat per unit mass: J K-1 kg-1.
    Thermal conductivity: W m-1 K-1.

    Arguments:
        - `atmos::Atmos_t`      the atmosphere struct instance to be used.
        - `idx::Int64`          index of the layer
    """
    function calc_single_cpkc!(atmos::atmosphere.Atmos_t, idx::Int64)
        # Reset
        mmr::Float64 = 0.0
        atmos.layer_cp[idx] = 0.0
        atmos.layer_kc[idx] = 0.0
        # Loop over gases
        for gas in atmos.gas_names
            mmr = atmos.gas_vmr[gas][idx] * atmos.gas_dat[gas].mmw/atmos.layer_μ[idx]
            atmos.layer_cp[idx] += mmr * species.get_Cp(atmos.gas_dat[gas], atmos.tmp[idx])
            atmos.layer_kc[idx] += mmr * species.get_Kc(atmos.gas_dat[gas], atmos.tmp[idx])
        end
        return nothing
    end

    """
    **Calculate the mass-density for all layers.**

    Requires temperature, pressure, mmw to be already have been set.

    Arguments:
        - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    """
    function calc_profile_density!(atmos::atmosphere.Atmos_t)
        # Loop over levels
        @inbounds for i in 1:atmos.nlev_c
            calc_single_density!(atmos, i)
        end
        return nothing
    end

    """
    **Calculate the mass-density of a single layer.**

    Requires temperature, pressure, mmw to be already have been set.

    Arguments:
    - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    - `idx::Int64`          index of the layer
    """
    function calc_single_density!(atmos::atmosphere.Atmos_t, idx::Int64)
        atmos.layer_ρ[idx] = density.calc_rho_mix(
                                    [atmos.gas_dat[gas] for gas in atmos.gas_names],
                                    [atmos.gas_vmr[gas][idx] for gas in atmos.gas_names],
                                    atmos.tmp[idx], atmos.p[idx],
                                    atmos.layer_μ[idx]
                                )
        return nothing
    end

        """
    **Calculate the scale height for all layers.**

    Requires temperature, pressure, mmw to be already have been set.

    Arguments:
        - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    """
    function calc_profile_Hp!(atmos::atmosphere.Atmos_t)
        @inbounds for i in 1:atmos.nlev_c
            calc_single_Hp!(atmos, i)
        end
        return nothing
    end

    """
    **Calculate the scale height of a single layer.**

    Requires temperature, pressure, mmw to be already have been set.

    Arguments:
    - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    - `idx::Int64`          index of the layer
    """
    function calc_single_Hp!(atmos::atmosphere.Atmos_t, idx::Int64)
        if atmos.real_gas
            atmos.layer_Hp[idx] = atmos.p[idx] / (atmos.layer_ρ[idx] * atmos.g[idx])
        else
            atmos.layer_Hp[idx] = phys.R_gas * atmos.tmp[idx] / (atmos.layer_μ[idx] * atmos.g[idx])
        end
        return nothing
    end


end
