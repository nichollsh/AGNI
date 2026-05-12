# This file is part of AGNI. License is GPL-3.0: https://www.gnu.org/licenses

"""
**Module for calculating the density of liquids and gases.**

Contains functions for evaluating equation of state for liquids and gases, and calculating
the density of a mixture of gases using Amagat's law. Also contains lookup tables and
ideal gas calculations for liquid and gas densities.
"""
module density

    import ..consts: BIGFLOAT, R_gas
    import ..species: Gas_t, is_vapour, EOS_IDEAL, EOS

    # Table of liquid-phase density for ocean calculation [kg/m^3]
    #     All taken from this website:
    #     https://encyclopedia.airliquide.com/water#properties
    const _lookup_liquid_rho::Dict{String, Float64} = Dict([
        ("H2O", 958.37 ),  # boiling
        ("CO2", 1178.4 ),  # triple
        ("H2" , 70.516 ),  # boiling
        ("CH4", 422.36 ),  # boiling
        ("CO" , 793.2  ),  # boiling
        ("N2" , 806.11 ),  # boiling
        ("NH3", 681.97 ),  # boiling
        ("SO2", 1461.1 ),  # boiling
    ])

    """
    **Evaluate the density of a liquid phase.**

    Returns BIGFLOAT density for unsupported phases, to avoid divide-by-zero error

    Arguments:
    - `name::String`    Name of liquid

    Returns:
    - `rho::Float64`    Density of liquid phase [kg m-3]
    """
    function liquid_rho(name::String)::Float64
        if name in keys(_lookup_liquid_rho)
            return _lookup_liquid_rho[name]
        else
            return BIGFLOAT
        end
    end

    """
    **Evaluate the density of a single gas using the ideal gas EOS.**

    Arguments:
    - `tmp::Float64`        temperature [K]
    - `prs::Float64`        pressure [Pa]
    - `mmw::Float64`        mean molecular weight [kg mol-1]

    Returns:
    - `rho::Float64`        mass density [kg m-3]
    """
    function _rho_ideal(tmp::Float64, prs::Float64, mmw::Float64)::Float64
        return prs * mmw / (tmp * R_gas)
    end

    """
    **Calculate the density of a gas using the most appropriate equation of state.**

    Arguments:
    - `tmp::Float64`        temperature [K]
    - `prs::Float64`        pressure [Pa]
    - `gas::Gas_t`          the gas struct to be used

    Returns:
    - `rho::Float64`        mass density [kg m-3]
    """
    function calc_rho_gas(tmp::Float64, prs::Float64, gas::Gas_t)::Float64
        if isequal(gas.eos, EOS_IDEAL) || !is_vapour(gas, tmp, prs)
            # analytical form of ideal gas equation of state
            return _rho_ideal(tmp, prs, gas.mmw)
        else
            # otherwise, will use tabulated real-gas EOS to evaluate the density
            return gas.eos_I(tmp, log10(prs))
        end
    end

    """
    **Calculate the density of a mixture of gases using Amagat's law.**

    Arguments:
    - `gas::Array{Gas_t,1}`     array of gases
    - `vmr::Array{Float64,1}`   array of volume mixing ratios
    - `tmp::Float64`            temperature [K]
    - `prs::Float64`            pressure [Pa]

    Returns:
    - `rho::Float64`            mass density [kg m-3]
    """
    function calc_rho_mix(gas::Array{Gas_t,1}, vmr::Array{Float64,1},
                            tmp::Float64, prs::Float64, mmw::Float64)::Float64

        ngas::Int64 = length(gas)

        # single gas case
        if ngas == 1
            return calc_rho_gas(tmp, prs, gas[1])
        end

        # calculate the density (and mass-mixing ratio) of each gas
        rho::Array{Float64, 1} = zeros(Float64, ngas)
        mmr::Array{Float64, 1} = zeros(Float64, ngas)
        for i in 1:ngas
            rho[i] = calc_rho_gas(tmp, prs, gas[i])
            mmr[i] = vmr[i] * gas[i].mmw / mmw
        end

        # add them together, assuming ideal additive volumes (inverse density)
        return 1.0 / sum(mmr[:] ./ rho[:])
    end

end
