# Contains physical data

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module phys 

    using NCDatasets
    using PCHIPInterpolation
    using LoggingExtras

    # Sources:
    # - Pierrehumbert (2010)
    # - NIST
    # - Sources outlined in Thermo files (https://github.com/nichollsh/Thermo)

    # Universal gas constant, J K-1 mol-1
    const R_gas::Float64 = 8.314462618 # NIST CODATA

    # Stefan-boltzmann constant, W m−2 K−4
    const sigma::Float64 =  5.670374419e-8 # NIST CODATA

    # Von Karman constant, dimensionless
    const k_vk::Float64 = 0.40  # Hogstrom 1988

    # Planck's constant, J s
    const h_pl::Float64 = 6.62607015e-34  # NIST CODATA

    # Boltzmann constant, J K-1 
    const k_B::Float64 = 1.380649e-23 # NIST CODATA

    # Speed of light, m s-1
    const c_vac::Float64 = 299792458.0 # NIST CODATA

    # A large floating point number 
    const fbig::Float64 = 1e90

    # List of elements included in the model 
    const elements_list::Array{String,1} = ["H","C","N","O","S","P", 
                                            "Fe","Mg","Si","Ca","Al"]

    # 0 degrees celcius, in kelvin 
    const zero_celcius::Float64 = 273.15

    # Newton's gravitational constant [m3 kg-1 s-2]
    const G_grav::Float64 = 6.67430e-11 # NIST CODATA

    # Molecule mean molecular weight, kg mol-1
    const lookup_mmw::Dict{String, Float64} = Dict{String, Float64}([
        # molecules
        ("H2O", 1.801530E-02 ), 
        ("CO2", 4.401000E-02 ), 
        ("O3", 4.799820E-02 ), 
        ("N2O", 4.401280E-02 ), 
        ("CO", 2.801060E-02 ), 
        ("CH4", 1.604300E-02 ), 
        ("O2", 3.199880E-02 ), 
        ("NO", 3.000610E-02 ), 
        ("SO2", 6.406280E-02 ), 
        ("NO2", 4.600550E-02 ), 
        ("NH3", 1.703060E-02 ), 
        ("HNO3", 6.301290E-02 ), 
        ("N2", 2.801340E-02 ), 
        ("TiO", 6.386600E-02 ), 
        ("VO", 6.694090E-02 ), 
        ("H2", 2.015880E-03 ), 
        ("OCS", 6.007500E-02 ), 
        ("FeH", 5.685300E-02 ), 
        ("CrH", 5.300400E-02 ), 
        ("PH3", 3.399758E-02 ), 
        ("C2H2", 2.603730E-02 ), 
        ("HCN", 2.702530E-02 ), 
        ("H2S", 3.408100E-02 ), 
        ("NO3", 6.301280E-02 ), 
        ("N2O5", 1.080104E-01 ), 
        ("HONO", 4.701340E-02 ), 
        ("HO2NO2", 7.901220E-02 ), 
        ("H2O2", 3.401470E-02 ), 
        ("C2H6", 3.006900E-02 ), 
        ("CH3", 1.503450E-02 ), 
        ("H2CO", 3.002600E-02 ), 
        ("HO2", 3.300670E-02 ), 
        ("HDO", 1.902140E-02 ), 
        ("HCl", 3.646100E-02 ), 
        ("HF", 2.000689E-02 ),

        # elements from https://iupac.qmul.ac.uk/AtWt/
        ("H",  1.008000000e-03 ),
        ("He",  4.002000000e-03 ),
        ("Li",  6.940000000e-03 ),
        ("Be",  9.012000000e-03 ),
        ("B",  1.081000000e-02 ),
        ("C",  1.201100000e-02 ),
        ("N",  1.400700000e-02 ),
        ("O",  1.599900000e-02 ),
        ("F",  1.899800000e-02 ),
        ("Ne",  2.017970000e-02 ),
        ("Na",  2.298900000e-02 ),
        ("Mg",  2.430500000e-02 ),
        ("Al",  2.698100000e-02 ),
        ("Si",  2.808500000e-02 ),
        ("P",  3.097300000e-02 ),
        ("S",  3.206000000e-02 ),
        ("Cl",  3.545000000e-02 ),
        ("Ar",  3.995000000e-02 ),
        ("K",  3.909830000e-02 ),
        ("Ca",  4.007800000e-02 ),
        ("Sc",  4.495500000e-02 ),
        ("Ti",  4.786700000e-02 ),
        ("V",  5.094150000e-02 ),
        ("Cr",  5.199610000e-02 ),
        ("Mn",  5.493800000e-02 ),
        ("Fe",  5.584500000e-02 ),
        ("Co",  5.893300000e-02 ),
        ("Ni",  5.869340000e-02 ),
        ("Cu",  6.354600000e-02 ),
        ("Zn",  6.538000000e-02 ),
        ("Ga",  6.972300000e-02 ),
        ("Ge",  7.263000000e-02 ),
        ("As",  7.492100000e-02 ),
        ("Se",  7.897100000e-02 ),
        ("Br",  7.990400000e-02 ),
        ("Kr",  8.379800000e-02 ),
        ("Rb",  8.546780000e-02 ),
        ("Sr",  8.762000000e-02 ),
        ("Y",  8.890500000e-02 ),
        ("Zr",  9.122400000e-02 ),
        ("Nb",  9.290600000e-02 ),
        ("Mo",  9.595000000e-02 ),
        ("Ru",  1.010700000e-01 ),
        ("Rh",  1.029050000e-01 ),
        ("Pd",  1.064200000e-01 ),
        ("Ag",  1.078682000e-01 ),
        ("Cd",  1.124140000e-01 ),
        ("In",  1.148180000e-01 ),
        ("Sn",  1.187100000e-01 ),
        ("Sb",  1.217600000e-01 ),
        ("Te",  1.276000000e-01 ),
        ("I",  1.269040000e-01 ),
        ("Xe",  1.312930000e-01 ),
        ("Cs",  1.329050000e-01 ),
        ("Ba",  1.373270000e-01 ),
        ("La",  1.389050000e-01 ),
        ("Ce",  1.401160000e-01 ),
        ("Pr",  1.409070000e-01 ),
        ("Nd",  1.442420000e-01 ),
        ("Sm",  1.503600000e-01 ),
        ("Eu",  1.519640000e-01 ),
        ("Gd",  1.572500000e-01 ),
        ("Tb",  1.589250000e-01 ),
        ("Dy",  1.625000000e-01 ),
        ("Ho",  1.649300000e-01 ),
        ("Er",  1.672590000e-01 ),
        ("Tm",  1.689340000e-01 ),
        ("Yb",  1.730450000e-01 ),
        ("Lu",  1.749668000e-01 ),
        ("Hf",  1.784860000e-01 ),
        ("Ta",  1.809470000e-01 ),
        ("W",  1.838400000e-01 ),
        ("Re",  1.862070000e-01 ),
        ("Os",  1.902300000e-01 ),
        ("Ir",  1.922170000e-01 ),
        ("Pt",  1.950840000e-01 ),
        ("Au",  1.969660000e-01 ),
        ("Hg",  2.005920000e-01 ),
        ("Tl",  2.043800000e-01 ),
        ("Pb",  2.072000000e-01 ),
        ("Bi",  2.089800000e-01 ),
        ("Th",  2.320377000e-01 ),
        ("Pa",  2.310350000e-01 ),
        ("U",  2.380280000e-01),
    ])

    # Structure containing data for a single gas
    mutable struct Gas_t 

        # Names 
        formula::String         # Formula used by SOCRATES
        JANAF_name::String      # JANAF name 
        fastchem_name::String   # FastChem name (to be determined from FC output file)

        # Is this a stub?
        #   This is the case when we cannot find an appropriate data file
        stub::Bool

        # Should evaluations be temperature-dependent or use constant values?
        tmp_dep::Bool

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
        sat_P::Array{Float64,1}     # Corresponding saturation pressures [Pa]
        sat_I::Interpolator         # Interpolator struct

        # Latent heat (enthalpy) of phase change
        lat_T::Array{Float64,1}     # Reference temperatures [K]
        lat_H::Array{Float64,1}     # Corresponding heats [J kg-1]
        lat_I::Interpolator         # Interpolator struct

        # Specific heat capacity
        cap_T::Array{Float64,1}     # Reference temperatures [K]
        cap_C::Array{Float64,1}     # Corresponding Cp values [J K-1 kg-1]
        cap_I::Interpolator         # Interpolator struct

        # Plotting colour (hex code) and label
        plot_color::String 
        plot_label::String 

        Gas_t() = new()
    end # end gas struct 

    # Pre-defined colors 
    const lookup_color::Dict{String, String} = Dict{String, String}([
        # common volatiles 
        ("H2O", "#027FB1" ),
        ("CO2", "#D24901" ),
        ("H2" , "#008C01" ),
        ("CH4", "#C720DD" ),
        ("CO" , "#D1AC02" ),
        ("N2" , "#870036" ),
        ("NH3", "#675200" ),
        ("S2" , "#FF8FA1" ),
        ("SO2", "#00008B" ),

        # volatile elements 
        ("H"  , "#0000ff"),
        ("C"  , "#ff0000"),
        ("O"  , "#00dd00"),
        ("N"  , "#eeaa00"),
        ("S" ,  "#ff22ff"),
        ("P" ,  "#33ccff"),
        ("He" , "#30FF71" ),

        # refractory elements 
        ("Fe" , "#888888"),
        ("Si" , "#aa2277"),
        ("Mg" , "#996633"),
    ])

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
    function get_mmw(m::String)::Float64

        # already defined?
        if m in keys(lookup_mmw)
           return lookup_mmw[m]
        end 

        # get atoms 
        atoms::Dict{String, Int} = count_atoms(m)

        # add up atoms 
        mmw::Float64 = 0.0
        for k in keys(atoms)
            mmw += lookup_mmw[k]*atoms[k]
        end

        return mmw
    end

    """
    Convert formula to pretty unicode string 
    """
    function pretty_name(gas::String)::String 
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
    function pretty_color(gas::String)::String 
        # Defined 
        if gas in keys(lookup_color)
            return lookup_color[gas]
        end 

        # Else, generate colour from atoms 
        atoms = count_atoms(gas)
        r::Float64 = 0.0
        g::Float64 = 0.0
        b::Float64 = 0.0
        for e in keys(atoms)
            r += parse(Int,lookup_color[e][2:3],base=16)*atoms[e]
            g += parse(Int,lookup_color[e][4:5],base=16)*atoms[e]
            b += parse(Int,lookup_color[e][6:7],base=16)*atoms[e]
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
    Load gas data into a new struct  
    """
    function load_gas(thermo_dir::String, formula::String, tmp_dep::Bool)::Gas_t

        @debug ("Loading data for gas $formula")

        # Clean input and get file path
        formula = String(strip(formula))
        fpath = joinpath(thermo_dir, "$formula.nc" )

        # Initialise struct 
        gas = Gas_t()
        gas.formula = formula 
        gas.tmp_dep = tmp_dep

        # Count atoms 
        gas.atoms = count_atoms(formula)
        for e in keys(gas.atoms)
            if !(e in elements_list)
                error("Gas '$formula' contains unsupported element '$e'")
            end 
        end 

        # Fastchem name (to be learned later)
        gas.fastchem_name = "_unknown"
        
        # Set plotting color and label 
        gas.plot_color = pretty_color(formula)
        gas.plot_label = pretty_name(formula)

        get_mmw(formula)

        # Check if we have data from file 
        gas.stub = !isfile(fpath)
        gas.no_sat = false
        if gas.stub
            # no data => generate stub
            @debug("    stub")

            gas.mmw = get_mmw(formula)
            gas.JANAF_name = "_unknown"

            # heat capacity set to zero
            gas.cap_T = [0.0, fbig]
            gas.cap_C = [0.0, 0.0]

            # latent heat set to zero 
            gas.lat_T = [0.0, fbig]
            gas.lat_H = [0.0, 0.0]

            # saturation pressure set to large value (ensures always gas phase)
            gas.sat_T = [0.0, fbig]
            gas.sat_P = [fbig, fbig]
            gas.no_sat = true

            # critical set to small value (always supercritical)
            gas.T_crit = 0.0
            gas.T_trip = 0.0

        else
            # have data => load from file 
            @debug("    ncdf")
            @debug "ALL DEBUG SUPPRESSED"
            with_logger(MinLevelLogger(current_logger(), Logging.Info-200)) do
                ds = Dataset(fpath,"r")

                # scalar properties  
                gas.mmw = ds["mmw"][1]
                gas.T_trip = ds["T_trip"][1]
                gas.T_crit = ds["T_crit"][1]
                gas.JANAF_name = String(ds["JANAF"][:])

                # variable properties 
                gas.cap_T = ds["cap_T"][:]
                gas.cap_C = ds["cap_C"][:]

                gas.lat_T = ds["lat_T"][:]
                gas.lat_H = ds["lat_H"][:]

                gas.sat_T = ds["sat_T"][:]
                gas.sat_P = ds["sat_P"][:]
                if length(gas.sat_P) < 3
                    gas.no_sat = true
                end 

                # close file 
                close(ds)
            end
            @debug "ALL DEBUG RESTORED"

            # extrapolate to high temperatures 
            push!(gas.cap_T, 9000.0); push!(gas.cap_C, gas.cap_C[end])
            push!(gas.lat_T, 9000.0); push!(gas.lat_H, gas.lat_H[end])
            push!(gas.sat_T, 9000.0); push!(gas.sat_P, gas.sat_P[end])

            # setup interpolators 
            gas.cap_I = Interpolator(gas.cap_T, gas.cap_C)
            gas.lat_I = Interpolator(gas.lat_T, gas.lat_H)
            gas.sat_I = Interpolator(gas.sat_T, gas.sat_P)

        end 

        @debug("    done")
        return gas
    end # end load_gas 

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
            return fbig
        end 
        if gas.no_sat
            return fbig 
        end 

        # Above critical point. In practice, a check for this should be made 
        #    before any attempt to evaluate this function.
        if t > gas.T_crit + 1.0e-5
            return fbig 
        end 

        # Get value from interpolator
        return gas.sat_I(t)
    end 

    """
    **Approximate the dew point temperature without interpolation**

    This should be avoided as much as possible.

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

        # Find closest value in array 
        i::Int = argmin(abs.(gas.sat_P .- p))
        return min(gas.sat_T[i], gas.T_crit)
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

    Arguments:
    - `gas::Gas_t`              the gas struct to be used
    - `t::Float64`              temperature [K]

    Returns:
    - `kc::Float64`             thermal conductivity [W m-1 K-1]
    """
    function get_Kc(gas::Gas_t, t::Float64=-1.0)::Float64 
        return 0.0
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

        # Calculate planck function value [W m-2 sr-1 m-1]
        # http://spiff.rit.edu/classes/phys317/lectures/planck.html
        flx = 2.0 * h_pl * c_vac * (c_vac / wav^5.0) / 
              ( exp(h_pl * c_vac / (wav * k_B * tmp)) - 1.0)

        # Integrate solid angle (hemisphere), convert units
        flx = flx * pi * 1.0e-9 # [W m-2 nm-1]

        return flx
    end 
    
end # end module 
