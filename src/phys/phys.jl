# Contains physical data

module phys

    # Import packages
    using LoggingExtras

    # Include internal modules
    using ..consts

    """
    **Evaluate the Planck function at a given wavelength and temperature.**

    Integrated over a hemisphere.

    Arguments:
    - `wave::Float64`       Wavelength [nm]
    - `tmp::Float64`        Temperature [K]

    Returns:
    - `flx::Float64`        Spectral flux density [W m-2 nm-1]
    """
    function evaluate_planck(wav::Float64, tmp::Float64)::Float64

        # Output value
        flx::Float64 = 0.0

        # Convert nm to m
        wav = wav * 1.0e-9

        # Optimisation variables
        wav5::Float64 = wav*wav*wav*wav*wav
        hc::Float64   = h_pl * c_vac

        # Calculate planck function value [W m-2 sr-1 m-1]
        # http://spiff.rit.edu/classes/phys317/lectures/planck.html
        flx = 2.0 * hc * (c_vac / wav5) / ( exp(hc / (wav * k_B * tmp)) - 1.0)

        # Integrate solid angle (hemisphere), convert units
        flx = flx * pi * 1.0e-9 # [W m-2 nm-1]

        return flx
    end

    """
    **Calculate gravitational acceleration.**

    Using the Newtonian formula for universal gravitation in spherical geometry.

    Arguments:
    - `mass::Float64`       Enclosed mass [kg]
    - `radius::Float64`     Enclosed radius [m]

    Returns:
    - `grav::Float64`       Gravitational acceleration [m s-2]
    """
    function grav_accel(mass::Float64, radius::Float64)::Float64
        return G_grav * mass / (radius * radius)
    end

    """
    **Calculate radial component of centripetal acceleration.**

    Note that this is only the component of the centripetal acceleration which
    acts in the radial direction, and is relevant for calculating the effective
    gravity at the surface of a rotating planet. This means `a_c` is zero at the poles.

    The full centripetal acceleration is given by the equation below, where `d` is the
    distance from the rotation axis and `p` is the rotation period.

    `a_c = 4 * pi^2 * d / p^2`

    The distance `d` is equal to `r * cos(θ)`, so the equation can be rewritten as:

    `a_c = 4 * pi^2 * d * cos(θ) / p^2`

    And then to get the radial component, we multiply by `cos(θ)` again:

    `a_c = r * ( 2 * pi * cos(θ) / p )^2`

    Arguments:
    - `p::Float64`       Axial rotation period [s]
    - `r::Float64`       Radius of this layer [m]
    - `θ::Float64`       Latitude of this column [degrees]

    Returns:
    - `a_c::Float64`     Centripetal acceleration [m s-2]
    """
    function cent_accel(p::Float64, r::Float64, θ::Float64)::Float64
        return r * ( 2 * pi * cosd(θ) / p )^2
    end


    """
    **Calculate thermal diffusivity.**

    https://en.wikipedia.org/wiki/Thermal_diffusivity?useskin=vector

    Arguments:
    - `k::Float64`       Thermal conductivity [W m-1 K-1]
    - `ρ::Float64`       Density [kg m-3]
    - `cp::Float64`      Specific heat capacity [J K-1 kg-1]

    Returns:
    - `α::Float64`       Thermal diffusivity [m2 s-1]
    """
    function calc_therm_diffus(k::Float64, ρ::Float64, cp::Float64)::Float64
        return k / (ρ * cp)
    end
    export calc_therm_diffus

    """
    **Calculate planetary equilibrium temperature.**

    https://en.wikipedia.org/wiki/Planetary_equilibrium_temperature?useskin=vector

    Arguments:
    - `S::Float64`       Bolometric instellation [W m-2]
    - `α::Float64`       Bond albedo

    Returns:
    - `Teq::Float64`     Planetary equilibrium temperature [K]
    """
    function calc_Teq(S::Float64, α::Float64)::Float64
        return (S*(1-α)/(4*consts.σSB))^0.25
    end
    export calc_Teq

    """
    **Calculate planetary *skin* temperature.**

    https://en.wikipedia.org/wiki/Skin_temperature_(atmosphere)?useskin=vector

    Arguments:
    - `S::Float64`       Bolometric instellation [W m-2]
    - `α::Float64`       Bond albedo

    Returns:
    - `Tskin::Float64`   Planetary skin temperature [K]
    """
    function calc_Tskin(S::Float64, α::Float64)::Float64
        return calc_Teq(S, α) * (0.5^0.25)
    end
    export calc_Tskin

end # end module
