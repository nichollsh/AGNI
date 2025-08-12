# Contains physical data

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module phys

    # Import external modules
    using NCDatasets
    using LoggingExtras
    import Interpolations: interpolate, Gridded, Linear, Flat, extrapolate, Extrapolation

    # Import internal modules
    include("consts.jl")
    include("blake.jl")
    using .consts
    import .blake

    # A large floating point number
    const BIGFLOAT::Float64     = 1e99
    const BIGLOGFLOAT::Float64  = 99.0

    # Minimum data file version [YYYYMMDD, as integer]
    const MIN_DATA_VERSION::Int64 = 20250220

    # Enable/disable flags
    ENABLE_CHECKSUM::Bool = true
    ENABLE_AQUA::Bool     = true
    ENABLE_CMS19::Bool    = true

    # Enumerate potential equations of state
    @enum EOS EOS_IDEAL=1 EOS_VDW=2 EOS_AQUA=3 EOS_CMS19=4

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
        prs_max::Float64

        # Constituent atoms (dictionary of numbers)
        atoms::Dict{String, Int}

        # Mean molecular weight [kg mol-1]
        mmw::Float64

        # Triple and critical points [K]
        T_trip::Float64
        T_crit::Float64

        # Saturation curve
        no_sat::Bool                # No saturation data
        sat_T::Array{Float64,1}     # Reference temperatures [K]
        sat_P::Array{Float64,1}     # Corresponding saturation pressures [log10 Pa]
        sat_I::Extrapolation        # Psat(T), 1D linear interpolator-extrapolator

        # Latent heat (enthalpy) of phase change
        lat_T::Array{Float64,1}     # Reference temperatures [K]
        lat_H::Array{Float64,1}     # Corresponding heats [J kg-1]
        lat_I::Extrapolation        # 1D linear interpolator-extrapolator

        # Specific heat capacity
        cap_T::Array{Float64,1}     # Reference temperatures [K]
        cap_C::Array{Float64,1}     # Corresponding Cp values [J K-1 kg-1]
        cap_I::Extrapolation        # 1D linear interpolator-extrapolator

        # Thermal conductivity
        kc::Float64                 # Constant conductivity [J K-1 kg-1]

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
        eos_I::Extrapolation        # 2D linear interpolator-extrapolator

        Gas_t() = new()
    end # end gas struct

    """
    Load gas data into a new struct
    """
    function load_gas(thermo_dir::String, formula::String,
                            tmp_dep::Bool, real_gas::Bool)::Gas_t

        @debug ("Loading data for gas $formula")

        # Clean input and get file path
        formula = String(strip(formula))
        fpath = joinpath(thermo_dir, "$formula.nc" )

        # Initialise struct
        gas = Gas_t()
        gas.formula = formula
        gas.tmp_dep = tmp_dep
        fail = false

        # Count atoms
        gas.atoms = count_atoms(formula)
        for e in keys(gas.atoms)
            if !(e in elems_standard)
                error("Gas '$formula' contains unsupported element '$e'")
            end
        end

        # Set plotting color and label
        gas.plot_color = _pretty_color(formula)
        gas.plot_label = _pretty_name(formula)

        # Thermal conductivity
        gas.kc = 0.0

        # Fastchem name (to be learned later)
        gas.fastchem_name = "_unknown"

        # Default parameters, assuming we have no data...
        gas.mmw = _get_mmw(formula)
        gas.JANAF_name = "_unknown"

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
        gas.prs_max = BIGFLOAT
        eos_name = "ideal gas"

        # Check if we have data from file
        gas.stub = !isfile(fpath)
        if gas.stub
            # no data
            @debug("    stub")

        elseif ENABLE_CHECKSUM && !blake.valid_file(fpath)
            # file exists - check its integrity
            @debug("    ncdf file is corrupt")
            fail = true

        else
            # have data => load what we can find inside the file
            @debug("    ncdf")

            # open the file
            with_logger(MinLevelLogger(current_logger(), Logging.Info)) do
            Dataset(fpath,"r") do ds

                # check date created
                created::Int64 = 0
                if !haskey(ds,"created")
                    @error("Data file ($formula) has no creation date")
                    fail = true
                    return gas
                end
                created = ds["created"][1]
                if created < MIN_DATA_VERSION
                    @error("Data file ($formula) is outdated ($created < $MIN_DATA_VERSION)")
                    fail = true
                    return gas
                end

                # we always have these
                gas.mmw = ds["mmw"][1]
                gas.JANAF_name = String(ds["JANAF"][:])

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
                        end
                    end

                    # try to use cms19 EOS for dihydrogen
                    if (formula == "H2") && ENABLE_CMS19
                        if haskey(ds, "cms19_T")
                            gas.eos = EOS_CMS19
                            # cms19 data found -  this is preferred
                        else
                            @warn("Could not find CMS19 table for H2 equation of state")
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
                        @error("Could not parse $formula EOS data from file")
                        @error("    temp. length = $(length(gas.eos_T))")
                        @error("    pres. length = $(length(gas.eos_P))")
                        @error("    dens. length = $(length(gas.eos_ρ))")
                        fail = true
                        return gas
                    end

                    # check ascending
                    if !issorted(gas.eos_P)
                        @error("Could not parse $formula EOS data from file")
                        @error("    Pressure array must be strictly ascending")
                        fail = true
                        return gas
                    end
                    if !issorted(gas.eos_T)
                        @error("Could not parse $formula EOS data from file")
                        @error("    Temperature array must be strictly ascending")
                        fail = true
                        return gas
                    end

                    # record valid T,P range
                    gas.tmp_max = maximum(gas.eos_T)
                    gas.prs_max = 10.0 ^ maximum(gas.eos_P)

                    # convert density to SI units
                    @. gas.eos_ρ = 10.0 ^ gas.eos_ρ

                    # interpolate to 2D grid
                    gas.eos_I = extrapolate(interpolate(
                                                        (gas.eos_T,gas.eos_P), gas.eos_ρ,
                                                        Gridded(Linear())), # linear interp.
                                            Flat()) # constant-value extrap.
                end # /EOS

            end # /NetCDF
            end # /MinLevelLogger

            # Setup 1D interpolators for Cp, Lv, and Psat
            gas.cap_I = extrapolate(interpolate((gas.cap_T,), gas.cap_C, Gridded(Linear())), Flat())
            gas.lat_I = extrapolate(interpolate((gas.lat_T,), gas.lat_H, Gridded(Linear())), Flat())
            gas.sat_I = extrapolate(interpolate((gas.sat_T,), gas.sat_P, Gridded(Linear())), Flat())
        end

        @debug("    using '$eos_name' equation of state")
        @debug("    done")
        gas.fail = fail
        return gas
    end # end load_gas


    """
    Get number of atoms from formula, returning a dictionary
    """
    function count_atoms(m::String)::Dict
        # Setup
        out = Dict()
        nchar::Int = length(m)
        i::Int = 1
        elem::String = ""
        count::Int=-1
        last::Bool=false

        # Loop through string
        while i <= nchar
            last = (i == nchar)

            # new element
            if isuppercase(m[i])
                count = 0
                elem = string(m[i])
                if !last && islowercase(m[i+1])  # two letter element name
                    elem = elem*string(m[i+1])
                    i += 1
                end
            end

            last = (i == nchar)

            # get count
            if count == 0   # expecting number
                # number of atoms
                if last || isletter(m[i+1]) # got letter => count=1
                    count = 1
                else
                    count = parse(Int, m[i+1])
                end
                # repeated element
                if elem in keys(out)
                    out[elem] += count
                else
                    out[elem] = count
                end
                # reset
                elem = ""
                count = -1
            end
            i += 1
        end

        return out
    end

    """
    Check if two gas atom dicts are equivalent
    """
    function same_atoms(d1::Dict, d2::Dict)::Bool

        # check if have same atoms at all
        for k in keys(d1)
            if !(k in keys(d2))
                return false
            end
        end

        # ^^ reverse combination
        for k in keys(d2)
            if !(k in keys(d1))
                return false
            end
        end

        # check counts
        for k in keys(d1)
            if d1[k] != d2[k]
                return false
            end
        end

        # if we haven't returned false so far, then it must be true
        return true
    end

    """
    Calculate species mean molecular weight [kg mol-1] from formula or use known value
    """
    function _get_mmw(m::String)::Float64

        # already defined?
        if m in keys(consts._lookup_mmw)
           return consts._lookup_mmw[m]
        end

        # get atoms
        atoms::Dict{String, Int} = count_atoms(m)

        # add up atoms
        mmw::Float64 = 0.0
        for k in keys(atoms)
            mmw += consts._lookup_mmw[k]*atoms[k]
        end

        return mmw
    end

    """
    Convert formula to pretty unicode string
    """
    function _pretty_name(gas::String)::String
        out::String = ""
        for c in gas
            if isnumeric(c)
                d = parse(Int, c)
                out *= Char(parse(Int,"208$d", base=16))
            else
                out *= c
            end
        end
        return out
    end

    """
    Generate a colour hex code from a molecular formula
    """
    function _pretty_color(gas::String)::String
        # Defined
        if gas in keys(consts._lookup_color)
            return consts._lookup_color[gas]
        end

        # Else, generate colour from atoms
        atoms = count_atoms(gas)
        r::Float64 = 0.0
        g::Float64 = 0.0
        b::Float64 = 0.0
        for e in keys(atoms)
            r += parse(Int,consts._lookup_color[e][2:3],base=16)*atoms[e]
            g += parse(Int,consts._lookup_color[e][4:5],base=16)*atoms[e]
            b += parse(Int,consts._lookup_color[e][6:7],base=16)*atoms[e]
        end
        m::Float64 = max(r,g,b)

        # prevents the colour getting too close to white
        if r+g+b > 705
            m *= 255.0/235.0
        end

        # convert to hex code
        out::String = "#"
        out *= string(floor(Int,255 * r/m),base=16,pad=2)
        out *= string(floor(Int,255 * g/m),base=16,pad=2)
        out *= string(floor(Int,255 * b/m),base=16,pad=2)
        return out
    end

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
        i::Int = argmin(abs.(gas.sat_P .- p))
        return min(gas.sat_T[i], gas.T_crit)
    end

    """
    **Check if pressure-temperature coordinate is within the vapour regime.**

    Returns true if p > p_sat and t < t_crit.

    Arguments:
    - `gas::Gas_t`              the gas struct to be used
    - `t::Float64`              temperature [K]
    - `p::Float64`              temperature [K]

    Returns:
    - `vapour::Bool`           is within vapour regime?
    """
    function is_vapour(gas::Gas_t, t::Float64, p::Float64)::Bool

        # Handle stub cases
        if gas.stub || gas.no_sat
            return true
        end

        # Above critical point?
        if t >= gas.T_crit
            return true
        end

        # Saturated by pressure? (with offset to account for transition)
        return Bool(p > 10.0 ^ gas.sat_I(t+0.2))
    end

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

    """
    **Get gas thermal conductivity at a given temperature.**

    This is always set to zero - not yet implemented.

    Arguments:
    - `gas::Gas_t`              the gas struct to be used
    - `t::Float64`              temperature [K]

    Returns:
    - `kc::Float64`             thermal conductivity [W m-1 K-1]
    """
    function get_Kc(gas::Gas_t, t::Float64=-1.0)::Float64
        return gas.kc
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

        # validate lengths
        ngas::Int = length(gas)
        if ngas != length(vmr)
            @error "The number of gases and the number of mixing ratios are different"
            exit(1)
        end

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
        if (gas.eos == EOS_IDEAL) || !is_vapour(gas, tmp, prs)
            # analytical form of ideal gas equation of state
            return _rho_ideal(tmp, prs, gas.mmw)
        else
            # otherwise, will use tabulated real-gas EOS to evaluate the density
            return gas.eos_I(tmp, log10(prs))
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
    **Evaluate the density of a liquid phase.**

    Returns BIGFLOAT density for unsupported phases, to avoid divide-by-zero error

    Arguments:
    - `name::String`     Name of liquid

    Returns:
    - `rho::Float64`
    """
    function liquid_rho(name::String)::Float64
        if name in keys(_lookup_liquid_rho)
            return _lookup_liquid_rho[name]
        else
            return BIGFLOAT
        end
    end

end # end module
