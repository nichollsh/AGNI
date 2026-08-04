# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0

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

    Evaluates `prs` without knowledge of whether it is the total or partial pressure.
    The density of mixtures is calculated in `calc_rho_mix` using Amagat's law, in which
    case this function should be evaluated at the total pressure (not partial pressure).

    The flag `phs_method` is an integer specifying how to handle EOS evaluations with phase
    boundaries, which is tricky for non-ideal equations of state. A check is made against
    phase boundaries with a small tolerance, set within the call to `is_vapour()`. Note that
    the tolerance `phs_εlogp` in `is_vapour()` differs from the offset `phs_dlogp` on psat.

    Options for `phs_method` flag:
    - `1`, naively evaluate the density at the requested T-P.
    - `2`, switch to ideal gas when `is_vapour() == false`.
    - `3`, [default] evaluate real gas EOS with clamp applied to log10 pressure in the form
            `eval_log10prs = min(log10(prs), log10(prs_sat) - phs_dlogp)`, to ensure that
            density is evaluated in the vapour region.

    Arguments:
    - `tmp::Float64`        temperature [K]
    - `prs::Float64`        pressure [Pa]
    - `gas::Gas_t`          the gas struct to be used
    - `phs_method::Int64`   method for handling P-T conditions near non-vapour regions
    - `phs_dlogp::Float64`  log10 pressure offset relative to saturation pressure

    Returns:
    - `rho::Float64`        mass density [kg m-3]
    """
    function calc_rho_gas(tmp::Float64, prs::Float64, gas::Gas_t;
                            phs_method::Int64=2,
                            phs_dlogp::Float64=0.2)::Float64

        # log10 pressure for evaluating EOS (Pa)
        eval_log10prs::Float64 = log10(prs)

        # determine which EOS to use, based on the `vap_enforce` flag...
        #     firstly, check if the gas is specified to be ideal
        eval_ideal::Bool = isequal(gas.eos, EOS_IDEAL)

        #     next, check if we are within a non-vapour regime (see docstring)
        #     this uses a small negative tolerance to ensure that exactly-saturated cases
        #     are treated as condensates, enabling special treatment if phs_method=2,3.
        if !eval_ideal && !is_vapour(gas, tmp, prs; phs_εlogp=-phs_dlogp)
            if phs_method == 1
                # don't need to do anything
            elseif phs_method == 2
                # switch to ideal gas
                eval_ideal = true
            elseif phs_method == 3
                # use shifted log10pressure (same temperature) for evaluating density
                # this decreases the pressure such that we fall within the vapour region
                # @debug "Applying clamp at phase boundary (old logP=$eval_log10prs)"
                eval_log10prs = min(eval_log10prs, gas.sat_I(tmp) - phs_dlogp)
                # @debug "    new logP=$eval_log10prs => vapour=$(is_vapour(gas, tmp, 10.0^eval_log10prs))"
            end
        else
            @debug "Is vapour: T=$tmp, P=$prs, logP=$eval_log10prs"
        end

        # evaluate EOS
        if eval_ideal
            # analytical form of ideal gas equation of state
            # this doesn't care about phase boundaries
            # @debug "Evaluating ideal gas: T=$tmp, P=$prs"
            return _rho_ideal(tmp, prs, gas.mmw)
        else
            # otherwise, will use tabulated real-gas EOS to evaluate the density
            # this requires careful handling of phase boundaries
            # @debug "Evaluating real gas:  T=$tmp, new logP=$eval_log10prs"
            return gas.eos_I(tmp, eval_log10prs)
        end
    end

    """
    **Calculate the density of a mixture of gases using Amagat's law.**

    This evaluates the density of each component at a given temperature and pressure. It is
    important that the *total* pressure is used for each species, since we then weight the
    density of each species by its mass mixing ratio.

    Arguments:
    - `gas::Array{Gas_t,1}`     array of gases
    - `vmr::Array{Float64,1}`   array of volume mixing ratios
    - `tmp::Float64`            temperature [K]
    - `prs::Float64`            total pressure [Pa]

    Returns:
    - `rho::Float64`            mass density [kg m-3]
    """
    function calc_rho_mix(gas::Array{Gas_t,1}, vmr::Array{Float64,1},
                            tmp::Float64, prs::Float64, mmw::Float64)::Float64

        ngas::Int64 = length(gas)

        # single gas case
        # (the total pressure is identical to the partial pressure)
        if ngas == 1
            return calc_rho_gas(tmp, prs, gas[1])
        end

        # calculate the density (and mass-mixing ratio) of each gas
        rho::Array{Float64, 1} = zeros(Float64, ngas)
        mmr::Array{Float64, 1} = zeros(Float64, ngas)
        for i in 1:ngas
            rho[i] = calc_rho_gas(tmp, prs, gas[i]) # total temperature and pressure
            mmr[i] = vmr[i] * gas[i].mmw / mmw # convert VMR to MMR
        end

        # add them together, assuming ideal additive volumes (inverse density)
        return 1.0 / sum(mmr[:] ./ rho[:])
    end

end
