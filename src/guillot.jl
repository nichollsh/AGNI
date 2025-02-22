# Contains module for Guillot analytical T(p) profiles

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module guillot

    import SpecialFunctions:expinti

    # Scalars
    κ_vs::Float64 = 4e-3 * 10  # m2/kg
    κ_th::Float64 = 1e-2 * 10  # m2/kg
    γ::Float64    = κ_vs/κ_th  # opacity ratio

    # Evalulate exponential integral of order 1
    function _E1(z::Float64)::Float64
        return -1*expinti(-z)
    end

    # Evaluate exponential integral of order 2
    function _E2(z::Float64)::Float64
        return exp(-z) - z*_E1(z)
    end

    # Evalulate irradiation temperature
    function eval_Tirr(Tstar::Float64, Rstar::Float64, sep::Float64)::Float64
        return Tstar * sqrt(Rstar/sep)
    end

    # Evalulate planetary equilibrium temperature
    function eval_Teqm(Tstar::Float64, Rstar::Float64, sep::Float64)::Float64
        return Tstar * sqrt(Rstar/(2*sep))
    end

    """
    **Evalulate planet-average T^4 at a given optical depth.**

    Equation 49 in Guillot (2010).
    """
    function _eval_T4_avg(τ::Float64, Tint::Float64, Teqm::Float64)::Float64
        pt1 = (3.0 * Tint^4 / 4) * (2.0/3 + τ)
        pt2 = (2.0/(3*γ)) * (1+  (γ*τ/2 - 1)* exp(-γ*τ) )
        pt3 = (2.0*γ/3) * (1- τ*τ/2)*_E2(γ*τ)
        return pt1 + (3.0 * Teqm^4 / 4) * (2.0/3 + pt2 + pt3)
    end

    """
    **Evalulate collimated-beam T^4 at a given optical depth.**

    Equation 27 in Guillot (2010).
    """
    function _eval_T4_cos(τ::Float64, Tint::Float64, Tirr::Float64, θ::Float64)::Float64
        mu  = cos(θ * pi / 180)
        pt1 = (3.0 * Tint^4 / 4) * (2.0/3 + τ)
        pt2 = mu/γ
        pt3 = (γ/(3*mu) - mu/γ) * exp(-γ*τ/mu)
        return pt1 + (3* Tirr^4 / 4)*mu*(2.0/3 + pt2 + pt3)
    end

    """
    **Evalulate LW optical depth as a function of pressure.**

    Assuming constant gravity.
    """
    function eval_tau(p::Float64, grav::Float64)::Float64
        return p * κ_th / grav
    end

    """
    **Calculate planet-average temperature profile.**
    """
    function calc_profile_avg(p_arr::Array{Float64,1}, grav::Float64,
                                Tint::Float64, Teqm::Float64)::Array{Float64, 1}

        t_arr = zeros(Float64, length(p_arr))
        for (i,p) in enumerate(p_arr)
            t_arr[i] = eval_T4_avg(eval_tau(p, grav), Tint, Teqm)^0.25
        end
        return t_arr
    end


    """
    **Calculate temperature profile for a given zenith angle.**
    """
    function calc_profile_cos(p_arr::Array{Float64,1}, grav::Float64,
                                Tint::Float64, Tirr::Float64, θ::Float64)::Array{Float64, 1}

        t_arr = zeros(Float64, length(p_arr))
        for (i,p) in enumerate(p_arr)
            t_arr[i] = eval_T4_cos(eval_tau(p, grav), Tint, Tirr, θ)^0.25
        end
        return t_arr
    end

end
