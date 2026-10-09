# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0

"""
**Module for handling gas/vapour species data and properties.**

Defines the `Gas_t` struct, which contains all relevant data for a single
species, and functions to load this data from files and evaluate properties such as
saturation pressure, heat capacity, and thermal conductivity.

Note that 'gas' here refers to any vapour species, including super-critical phases.
"""
module species

    # Import packages
    using NCDatasets
    using LoggingExtras
    import Interpolations: interpolate, Gridded, Linear, Flat, extrapolate, Extrapolation

    # Import local modules
    using ..consts
    import ..consts: _lookup_lj
    using ..formulae
    import ..blake: valid_file
    import ..style: pretty_colour, pretty_name

    # Minimum data file version [YYYYMMDD, as integer]
    const MIN_DATA_VERSION::Int64 = 20260201

    # Pressure limits for EOS evaluation (should be consistent with NetCDF data files)
    const EOS_LOGPMIN::Float64 = 0.0   # log10 Pa
    const EOS_LOGPMAX::Float64 = 11.0  # log10 Pa

    # Fallback value for particle size [m]
    const FALLBACK_SIZE::Float64 = 2e-10 # generic hard sphere of 2 Å

    # Enumerate potential equations of state
    @enum EOS EOS_IDEAL=1 EOS_VDW=2 EOS_AQUA=3 EOS_CMS19=4
    export EOS, EOS_IDEAL, EOS_VDW, EOS_AQUA, EOS_CMS19

    # Enable/disable flags
    ENABLE_CHECKSUM::Bool = true  # can still be disabled when function is called
    ENABLE_AQUA::Bool     = true
    ENABLE_CMS19::Bool    = true


    # Structure containing data for a single gas
    mutable struct Gas_t

        # Names
        formula::String         # Formula used by SOCRATES
        JANAF_name::String      # JANAF name
        fastchem_name::String   # FastChem name (to be determined from FC output file)

        # Is this a stub?
        stub::Bool

        # Fail if file found but cannot be parsed
        fail::Bool

        # Should evaluations be temperature-dependent or use constant values?
        tmp_dep::Bool

        # Maximum valid range for T,P
        tmp_max::Float64
        log10prs_max::Float64
        log10prs_min::Float64

        # Constituent atoms (dictionary of numbers)
        atoms::Dict{String, Int64}

        # Mean molecular weight [kg mol-1]
        mmw::Float64

        # Triple and critical points [K]
        T_trip::Float64
        T_crit::Float64

        # Saturation curve
        no_sat::Bool                # No saturation data
        sat_T::Array{Float64,1}     # Reference temperatures [K]
        sat_P::Array{Float64,1}     # Corresponding saturation pressures [log10 Pa]
        sat_I::Extrapolation        # log10 Psat(T), 1D linear interpolator-extrapolator

        # Latent heat (enthalpy) of phase change
        lat_T::Array{Float64,1}     # Reference temperatures [K]
        lat_H::Array{Float64,1}     # Corresponding heats [J kg-1]
        lat_I::Extrapolation        # 1D linear interpolator-extrapolator

        # Specific heat capacity
        cap_T::Array{Float64,1}     # Reference temperatures [K]
        cap_C::Array{Float64,1}     # Corresponding Cp values [J K-1 kg-1]
        cap_I::Extrapolation        # 1D linear interpolator-extrapolator

        # Particle mass [kg], collision diameter [m], and Lennard-Jones well depth [K]
        particle_m::Float64
        particle_d::Float64
        lj_eps::Float64             # ε/k_B; NaN when Lennard-Jones data are not available

        # Plotting colour (hex code) and label
        plot_color::String
        plot_label::String

        # Which equation of state should be used for this gas?
        eos::EOS

        # EOS original grid (flattened 2D arrays)
        eos_T::Array{Float64,1}     # temperature [K]
        eos_P::Array{Float64,1}     # log pressure [log10 Pa]
        eos_ρ::Array{Float64,2}     # log density [kg m-3]

        # EOS interpolator with constant-value extrapolation
        eos_I::Extrapolation        # log10 rho(T,P),2D linear interpolator-extrapolator

        Gas_t() = new()
    end # end gas struct
    export Gas_t

    """
    **Load gas data into a new struct.**

    Arguments:
    - `thermo_dir::String`      directory containing thermodynamic data files
    - `formula::String`         molecular formula of the gas (e.g. "H2O")
    - `tmp_dep::Bool`           enable temperature-dependent thermodynamic evaluations?
    - `real_gas::Bool`          use a real-gas EOS if available
    - `check_integrity::Bool`   check the integrity of the data file

    Returns:
    - `gas::Gas_t`              struct containing gas data
    """
    function load_gas(thermo_dir::String, formula::String,
                            tmp_dep::Bool, real_gas::Bool;
                            check_integrity::Bool=true)::Gas_t

        @debug ("Loading data for gas $formula")

        # Clean input and get file path
        formula = String(strip(formula))
        fpath = joinpath(thermo_dir, "$formula.nc" )

        # Initialise struct
        gas = Gas_t()
        gas.formula = formula
        gas.tmp_dep = tmp_dep
        gas.fail = false

        # Count atoms
        gas.atoms = count_atoms(formula)
        for e in keys(gas.atoms)
            if !(e in elems_standard)
                @error "Gas '$formula' contains unsupported element '$e'"
                gas.fail = true
                return gas
            end
        end

        # Set plotting color and label
        gas.plot_color = pretty_colour(formula)
        gas.plot_label = pretty_name(formula)

        # Fastchem name (to be learned later)
        gas.fastchem_name = "_unknown"

        # Default parameters, assuming we have no data...
        gas.mmw = get_mmw(formula)
        gas.JANAF_name = "_unknown"

        # Collision diameter [m] and well depth [K] from the Lennard-Jones table,
        if haskey(_lookup_lj, formula)
            gas.particle_d, gas.lj_eps = _lookup_lj[formula]
        else
            gas.particle_d = FALLBACK_SIZE
            gas.lj_eps = NaN
        end

        # heat capacity set to zero
        gas.cap_T = [0.0, BIGFLOAT]
        gas.cap_C = [Cp_ideal/gas.mmw, Cp_ideal/gas.mmw]

        # latent heat set to zero
        gas.lat_T = [0.0, BIGFLOAT]
        gas.lat_H = [0.0, 0.0]

        # saturation pressure set to large value (ensures always gas phase)
        gas.sat_T = [0.0, BIGFLOAT]
        gas.sat_P = [BIGLOGFLOAT, BIGLOGFLOAT]
        gas.no_sat = true

        # critical set to small value (always supercritical)
        gas.T_crit = 0.0
        gas.T_trip = 0.0

        # set EOS to ideal gas
        gas.eos = EOS_IDEAL
        gas.tmp_max = BIGFLOAT
        gas.log10prs_max = BIGLOGFLOAT
        gas.log10prs_min = SMALLLOGFLOAT
        eos_name = "ideal gas"

        # Check if we have data from file
        gas.stub = !isfile(fpath)
        if gas.stub
            # no data
            @debug("    stub")

        elseif ENABLE_CHECKSUM && check_integrity && !valid_file(fpath)
            # file exists - check its integrity
            @warn("    ncdf file is corrupt: '$fpath'")
            gas.fail = true
            return gas

        else
            # have data => load what we can find inside the file
            @debug("    ncdf")

            # open the file
            with_logger(MinLevelLogger(current_logger(), Logging.Info)) do
            Dataset(fpath,"r") do ds

                # check date created
                created::Int64 = 0
                if !haskey(ds,"created")
                    @warn("Data file ($formula) has no creation date")
                    gas.fail = true
                    return gas
                end
                created = ds["created"][1]
                if created < MIN_DATA_VERSION
                    @warn("Data file ($formula) is outdated ($created < $MIN_DATA_VERSION)")
                    gas.fail = true
                    return gas
                end

                # we always have these
                gas.mmw = ds["mmw"][1]
                gas.JANAF_name = String(ds["JANAF"][1])

                # triple point and critical point
                if haskey(ds, "T_trip")
                    gas.T_trip = ds["T_trip"][1]
                end
                if haskey(ds, "T_crit")
                    gas.T_crit = ds["T_crit"][1]
                end

                # heat capacity
                if haskey(ds, "cap_T")
                    gas.cap_T = ds["cap_T"][:]
                    gas.cap_C = ds["cap_C"][:]
                end

                # latent heat of phase change
                if haskey(ds, "lat_T")
                    gas.lat_T = ds["lat_T"][:]
                    gas.lat_H = ds["lat_H"][:]
                end

                # saturation pressure
                if haskey(ds, "sat_T")
                    gas.sat_T = ds["sat_T"][:] # K
                    gas.sat_P = ds["sat_P"][:] # log10 Pa
                    gas.no_sat = false
                end

                # work out which is the best available equation of state
                if real_gas
                    # try to use van der waals EOS
                    if haskey(ds, "vdw_T")
                        gas.eos = EOS_VDW
                    end

                    #  try to use aqua EOS for water
                    if (formula == "H2O") && ENABLE_AQUA
                        if haskey(ds, "aqua_T")
                            gas.eos = EOS_AQUA
                            # aqua data found -  this is preferred
                        else
                            @warn("Could not find AQUA table for H2O equation of state")
                            @warn("    Using ideal gas EOS for H2O")
                        end
                    end

                    # try to use cms19 EOS for dihydrogen
                    if (formula == "H2") && ENABLE_CMS19
                        if haskey(ds, "cms19_T")
                            gas.eos = EOS_CMS19
                            # cms19 data found -  this is preferred
                        else
                            @warn("Could not find CMS19 table for H2 equation of state")
                            @warn("    Using ideal gas EOS for H2")
                        end
                    end

                end

                # prepare eos data if necessary
                if gas.eos != EOS_IDEAL
                    # load data (T and P are 1d, rho is 2d)
                    if gas.eos == EOS_VDW
                        eos_name = "Van der Waals"
                        gas.eos_P = ds["vdw_P"][:]      # log Pa
                        gas.eos_T = ds["vdw_T"][:]      # K
                        gas.eos_ρ = ds["vdw_rho"][:,:]  # log kg m-3 (converted later)
                    elseif gas.eos == EOS_AQUA
                        eos_name = "AQUA"
                        gas.eos_P = ds["aqua_P"][:]
                        gas.eos_T = ds["aqua_T"][:]
                        gas.eos_ρ = ds["aqua_rho"][:,:]
                    elseif gas.eos == EOS_CMS19
                        eos_name = "CMS19"
                        gas.eos_P = ds["cms19_P"][:]
                        gas.eos_T = ds["cms19_T"][:]
                        gas.eos_ρ = ds["cms19_rho"][:,:]
                    end

                    # check shape
                    if !(length(gas.eos_ρ) == length(gas.eos_P) * length(gas.eos_T))
                        @warn("Could not parse $formula EOS data from file")
                        @warn("    temp. length = $(length(gas.eos_T))")
                        @warn("    pres. length = $(length(gas.eos_P))")
                        @warn("    dens. length = $(length(gas.eos_ρ))")
                        gas.fail = true
                        return gas
                    end

                    # check ascending
                    if !issorted(gas.eos_P)
                        @warn("Could not parse $formula EOS data from file")
                        @warn("    Pressure array must be strictly ascending")
                        gas.fail = true
                        return gas
                    end
                    if !issorted(gas.eos_T)
                        @warn("Could not parse $formula EOS data from file")
                        @warn("    Temperature array must be strictly ascending")
                        gas.fail = true
                        return gas
                    end

                    # record valid T,P range
                    gas.tmp_max      = maximum(gas.eos_T)
                    gas.log10prs_max = min(maximum(gas.eos_P), EOS_LOGPMAX)
                    gas.log10prs_min = max(minimum(gas.eos_P), EOS_LOGPMIN)

                    # ensure min/max are compatible
                    if gas.log10prs_min >= gas.log10prs_max
                        @warn("Could not parse $formula EOS data from file")
                        @warn("    The valid pressure domain is too small")
                        @warn("    log10(Pmin/Pa) = $(gas.log10prs_min)")
                        @warn("    log10(Pmax/Pa) = $(gas.log10prs_max)")
                        gas.fail = true
                        return gas
                    end

                    # interpolate to 2D grid
                    gas.eos_I = extrapolate(interpolate(
                                                        (gas.eos_T,gas.eos_P), gas.eos_ρ,
                                                        Gridded(Linear())), # linear interp.
                                            Flat()) # constant-value extrap.
                end # /EOS

            end # /NetCDF
            end # /MinLevelLogger

            # Some data files have a saturation curve but no known critical point.
            # Take the critical point to be the highest T of the saturation curve.
            if !gas.no_sat && (gas.T_crit <= minimum(gas.sat_T))
                gas.T_crit = maximum(gas.sat_T)
                @debug "    $formula: updated T_crit=$(gas.T_crit) K from saturation curve"
            end

            # Setup 1D interpolators for Cp, Lv, and Psat
            gas.cap_I = extrapolate(interpolate((gas.cap_T,), gas.cap_C, Gridded(Linear())), Flat())
            gas.lat_I = extrapolate(interpolate((gas.lat_T,), gas.lat_H, Gridded(Linear())), Flat())
            gas.sat_I = extrapolate(interpolate((gas.sat_T,), gas.sat_P, Gridded(Linear())), Flat())
        end

        # Particle mass [kg] from the final molar mass (N_A = R_gas / k_B)
        gas.particle_m = gas.mmw * k_B / R_gas

        @debug("    using '$eos_name' equation of state")
        @debug("    done")
        return gas
    end # end load_gas
    export load_gas

    """
    **Get gas saturation pressure for a given temperature.**

    If the temperature is above the critical point, then a large value
    is returned.

    Arguments:
    - `gas::Gas_t`              the gas struct to be used
    - `t::Float64`              temperature [K]

    Returns:
    - `p::Float64`              saturation pressure [Pa]
    """
    function get_Psat(gas::Gas_t, t::Float64)::Float64

        # Handle stub cases
        if gas.stub
            return BIGFLOAT
        end
        if gas.no_sat
            return BIGFLOAT
        end

        # Above critical point. In practice, a check for this should be made
        #    before any attempt to evaluate this function.
        if t > gas.T_crit + 1.0e-5
            return BIGFLOAT
        end

        # Get value from interpolator
        return 10.0 ^ gas.sat_I(t)
    end
    export get_Psat

    """
    **Get gas dew point temperature for a given partial pressure.**

    If the pressure is below the critical point pressure, then T_crit is returned.
    This function is horrendous, and should be avoided at all costs.

    Arguments:
    - `gas::Gas_t`              the gas struct to be used
    - `p::Float64`              pressure [Pa]

    Returns:
    - `t::Float64`              dew point temperature [K]
    """
    function get_Tdew(gas::Gas_t, p::Float64)::Float64

        # Handle stub case
        if gas.stub
            return 0.0
        end

        p = log10(p)

        # Find closest value in array
        i::Int64 = argmin(abs.(gas.sat_P .- p))
        return min(gas.sat_T[i], gas.T_crit)
    end
    export get_Tdew

    """
    **Check if pressure-temperature coordinate is within the vapour regime.**

    Returns true if p < p_sat and t < t_crit.
    Also returns true if no phase-change data are available for this gas.

    This function performs the saturation check in log10 pressure, and applies small
    numerical tolerance on the inequality so that exactly-saturated cases are
    considered to be vapours.

    Arguments:
    - `gas::Gas_t`              the gas struct to be used
    - `t::Float64`              temperature [K]
    - `p::Float64`              partial pressure [Pa]
    - `phs_εlogp::Float64`      tolerance for log10(pressure) inequality comparison

    Returns:
    - `vapour::Bool`           is within vapour regime?
    """
    function is_vapour(gas::Gas_t, t::Float64, p::Float64;
                        phs_εlogp::Float64=1e-5)::Bool

        # Handle stub cases
        if gas.stub || gas.no_sat
            return true
        end

        # Above critical point?
        if t >= gas.T_crit
            return true
        end

        # Is vapour when p < p_sat (see docstring)
        # Comparison is made in log10-space with a small tolerance
        # A `phs_εlogp>0` means that super-saturated pressures are considered to be vapour
        # A `phs_εlogp<0` means that sub-saturated pressures are considered to be vapour
        return log10(p) - phs_εlogp < gas.sat_I(t)
    end
    export is_vapour

    """
    **Get gas enthalpy (latent heat) of phase change.**

    If the temperature is above the critical point, then a zero value
    is returned. Evaluates at 0 Celcius if `gas.tmp_dep=false`.

    Arguments:
    - `gas::Gas_t`              the gas struct to be used
    - `t::Float64`              temperature [K]

    Returns:
    - `h::Float64`              enthalpy of phase change [J kg-1]
    """
    function get_Lv(gas::Gas_t, t::Float64)::Float64

        # Handle stub case
        if gas.stub
            return gas.lat_H[1]
        end

        # Above critical point
        if t > gas.T_crit
            return 0.0
        end

        # Constant value
        if !gas.tmp_dep
            t = zero_celcius
        end

        # Get value from interpolator
        return gas.lat_I(t)
    end
    export get_Lv

    """
    **Get gas heat capacity for a given temperature.**

    Evaluates at 0 Celcius if `gas.tmp_dep=false`.

    Arguments:
    - `gas::Gas_t`              the gas struct to be used
    - `t::Float64`              temperature [K]

    Returns:
    - `cp::Float64`             heat capacity of gas [J K-1 kg-1]
    """
    function get_Cp(gas::Gas_t, t::Float64)::Float64

        # Handle stub case
        if gas.stub
            return gas.cap_C[1]
        end

        # Constant value
        if !gas.tmp_dep
            t = zero_celcius
        end

        # Temperature floor, since we can get weird behaviour as Cp -> 0.
        t = max(t, 0.5)

        # Get value from interpolator
        return gas.cap_I(t)
    end
    export get_Cp

    """
    **Reduced collision integral Ω(2,2)* for the Lennard-Jones 12-6 potential.**

    Empirical fit of Neufeld, Janzen & Aziz (1972), J. Chem. Phys. 57, 1100, Table I,
    without its small sinusoidal term. The fit was made over 0.3 ≤ T* ≤ 100.

    Source: https://archive.org/download/wikipedia-scholarly-sources-corpus/10.1063%252F1.100061.zip/10.1063%252F1.1678363.pdf

    Arguments:
    - `t_red::Float64`          reduced temperature T* = T / (ε/k_B)

    Returns:
    - `omega::Float64`          collision integral, dimensionless
    """
    function omega22(t_red::Float64)::Float64
        return 1.16145 / t_red^0.14874 +
               0.52487 * exp(-0.77320 * t_red) +
               2.16178 * exp(-2.43787 * t_red)
    end
    export omega22

    """
    **Get dilute-gas dynamic viscosity at a given temperature.**

    First-order Chapman-Enskog viscosity,
    `η = (5/16) sqrt(π m k_B T) / (π σ² Ω(2,2)*(T*))`, with the Lennard-Jones
    collision diameter σ and well depth ε of `_lookup_lj`.

    Gases without LJ coeffs are treated as a hard sphere. Dipole moments are ignored.

    - Chapman & Cowling (1970), The Mathematical Theory of Non-Uniform Gases
    - https://doi.org/10.1063/1.1678363 (Neufeld et al. 1972)

    Arguments:
    - `gas::Gas_t`              the gas struct to be used
    - `t::Float64`              temperature [K]

    Returns:
    - `eta::Float64`            dynamic viscosity [Pa s]
    """
    function get_Visc(gas::Gas_t, t::Float64)::Float64

        # Constant value
        if !gas.tmp_dep
            t = zero_celcius
        end

        # Temperature floor, keeping the reduced temperature positive
        t = max(t, 0.5)

        # Collision integral
        omega::Float64 = isnan(gas.lj_eps) ? 1.0 : omega22(t / gas.lj_eps)

        return 5.0/16.0 * sqrt(pi * gas.particle_m * k_B * t) /
                    (pi * gas.particle_d^2 * omega)
    end
    export get_Visc

    """
    **Get single gas thermal conductivity at a given temperature.**

    Eucken relation `k = η (c_v + (9/4) R / M)` with `c_v = c_p - R/M` and M the molar mass.

    Monatomic gases use Chapman-Enskog result `k = (15/4) (R/M) η`, which the Eucken
    relation reduces to  when `c_v = (3/2) R/M`. Dipole effects are ignored

    - Poling, Prausnitz, O'Connell (2001), The Properties of Gases and Liquids (5th ed.)
    - Watson et al. 1981, Icarus 48, 150
    - https://webbook.nist.gov/chemistry/fluid/

    Arguments:
    - `gas::Gas_t`              the gas struct to be used
    - `t::Float64`              temperature [K]

    Returns:
    - `kc::Float64`             thermal conductivity [W m-1 K-1]
    """
    function get_Kc(gas::Gas_t, t::Float64)::Float64

        # Constant value
        if !gas.tmp_dep
            t = zero_celcius
        end

        # Specific gas constant [J K-1 kg-1]
        r_spec::Float64 = R_gas / gas.mmw

        # Monatomic: translation only, independent of stub heat capacity
        if sum(values(gas.atoms)) == 1
            return get_Visc(gas, t) * 3.75 * r_spec
        end

        # Heat capacity at constant volume, bounded below by translational DOF: 3/2=1.5
        cv::Float64 = max(get_Cp(gas, t) - r_spec, 1.5 * r_spec)

        # Eucken relation (9/4 = 2.25)
        return get_Visc(gas, t) * (cv + 2.25 * r_spec)
    end
    export get_Kc


    """
    **Calculate the demixing temperature for a given pressure and H2O molar fraction.**

    Source: https://www.aanda.org/articles/aa/pdf/2025/11/aa56322-25.pdf (Appendix A).

    Arguments:
    - p::Float64            Pressure in Pa (converted to kbar internally)
    - x::Float64            Molar fraction of H2O in the mixture (0-1)

    Returns:
    - Tdemix::Float64      Demixing temperature in K
    """
    function _Tdemix_H2O(p::Float64, x::Float64)::Float64

        Pkbar::Float64 = p * 1e-8  # Pa -> kbar

        # Table A1 Coefficients
        a::Float64 = 1.2035e-4
        b::Float64 = 0.5501
        c::Float64 = 1.9163e-2
        d::Float64 = 0.4498
        e::Float64 = -6.2253e-2
        f::Float64 = 74.5041
        g::Float64 = -3.1495e-4
        h::Float64 = 5.0828e6
        i::Float64 = 4.0719

        # Fitting function
        pt1 = (a / pi) * 0.5 * (b + c * Pkbar) / ((x - d)^2.0 + (0.5 * b)^2.0)
        pt2 = e * Pkbar^3.0 + f * Pkbar^2.0 + g * Pkbar + h
        return pt1 * pt2 + i * Pkbar
    end

    """
    **Calculate the gas demixing temperature, for a given pressure and molar fraction.**

    This is a wrapper function which calls the appropriate demixing fit.

    Arguments:
    - `gas::Gas_t`          the gas struct
    - `p::Float64`          pressure [Pa]
    - `x::Float64`          molar fraction of the gas in the mixture

    Returns:
    - `Tdemix::Float64`     demixing temperature [K]
    """
    function get_Tdemix(gas::Gas_t, p::Float64, x::Float64)::Float64
        if gas.formula == "H2O"
            return _Tdemix_H2O(p, x)
        else
            return -1.0 * BIGFLOAT # always above this temperature
        end
    end
    export get_Tdemix


end
