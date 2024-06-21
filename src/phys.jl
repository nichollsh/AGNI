# Contains physical data

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module phys 

    using NCDatasets
    using PCHIPInterpolation

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
    const elements_list = ["H","C","N","O","S","P", "Fe","Mg","Si","Ca","Al"]

    # Molecule mean molecular weight, kg mol-1
    const lookup_mmw = Dict{String, Float64}([
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
        ("CFC11", 1.373686E-01 ), 
        ("CFC12", 1.209140E-01 ), 
        ("CFC113", 1.873765E-01 ), 
        ("HCFC22", 8.646892E-02 ), 
        ("HFC125", 1.200223E-01 ), 
        ("HFC134A", 1.020318E-01 ), 
        ("CFC114", 1.709210E-01 ), 
        ("TiO", 6.386600E-02 ), 
        ("VO", 6.694090E-02 ), 
        ("H2", 2.015880E-03 ), 
        ("He", 4.002602E-03 ), 
        ("OCS", 6.007500E-02 ), 
        ("Na", 2.298977E-02 ), 
        ("K", 3.909830E-02 ), 
        ("FeH", 5.685300E-02 ), 
        ("CrH", 5.300400E-02 ), 
        ("Li", 6.941000E-03 ), 
        ("Rb", 8.546780E-02 ), 
        ("Cs", 1.329055E-01 ), 
        ("PH3", 3.399758E-02 ), 
        ("C2H2", 2.603730E-02 ), 
        ("HCN", 2.702530E-02 ), 
        ("H2S", 3.408100E-02 ), 
        ("Ar", 3.994800E-02 ), 
        ("O", 1.599940E-02 ), 
        ("N", 1.400674E-02 ), 
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
        ("HF", 2.000689E-02 )
    ])

    # Structure containing data for a single gas
    mutable struct Gas_t 

        # Names 
        formula::String     # Formula used by SOCRATES
        JANAF_name::String  # JANAF name 

        # Constituent atoms (dictionary of numbers)
        atoms::Dict{String, Int}

        # Mean molecular weight [kg mol-1]
        mmw::Float64

        # Triple and critical points [K]
        T_trip::Float64 
        T_crit::Float64 

        # Saturation curve 
        sat_T::Array{Float64,1}     # Reference temperatures [K]
        sat_P::Array{Float64,1}     # Corresponding saturation pressures [Pa]
        sat_i::Interpolator         # Interpolator struct

        # Latent heat (enthalpy) of phase change
        lat_T::Array{Float64,1}     # Reference temperatures [K]
        lat_H::Array{Float64,1}     # Corresponding heats [J kg-1]
        lat_i::Interpolator         # Interpolator struct

        # Specific heat capacity
        cap_T::Array{Float64,1}     # Reference temperatures [K]
        cap_C::Array{Float64,1}     # Corresponding Cp values [J mol-1 K-1]
        cap_i::Interpolator         # Interpolator struct

        # Plotting colour (hex code) and label
        plot_color::String 
        plot_label::String 

        Gas_t() = new()
    end # end gas struct 

    # Pre-defined colors 
    const lookup_color = Dict{String, String}([
        # common volatiles 
        ("H2O", "#027FB1" ),
        ("CO2", "#D24901" ),
        ("H2" , "#008C01" ),
        ("CH4", "#C720DD" ),
        ("CO" , "#D1AC02" ),
        ("N2" , "#870036" ),
        ("O2" , "#00008B" ),
        ("NH3", "#675200" ),

        # volatile elements 
        ("H"  , "#0000ff"),
        ("C"  , "#ff0000"),
        ("O"  , "#00ff00"),
        ("N"  , "#ffff00"),
        ("S" ,  "#00ffee"),

        # refractory elements 
        ("Fe" , "#ff22ff"),
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
    Calculate species mean molecular weight [kg mol-1] from formula or use known value
    """
    function get_mmw(m::String)::Float64

        # already defined?
        if m in keys(lookup_mmw)
           return lookup_mmw[formula]
        end 

        # get atoms 
        atoms = count_atoms(m)

        # add up atoms 
        mmw::Float64 = 0.0
        for k in keys(atoms)
            mmw += lookup_mmw[k]*atoms[k]
        end

        return 
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
    Generate a color for a formula 
    """
    function pretty_color(gas::String)::String 
        # Defined 
        if gas in keys(lookup_color)
            return lookup_color[gas]
        end 

        # Else, calculate color from atoms 
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

        out::String = "#"
        out *= string(floor(Int,255 * r/m),base=16,pad=2)
        out *= string(floor(Int,255 * g/m),base=16,pad=2)
        out *= string(floor(Int,255 * b/m),base=16,pad=2)
        return out
    end 

    """
    Load gas data into a new struct  
    """
    function load_gas(formula::String)::Gas_t

        # Clean input and get file path
        formula = strip(formula)
        fpath = joinpath(abspath(@__FILE__), "..", "res", "thermo", "$formula.ncdf" )

        # Initialise struct 
        gas = Gas_t()
        gas.formula = formula 

        # Count atoms 
        gas.atoms = count_atoms(formula)
        
        # Set plotting color and label 
        gas.plot_color = pretty_color(formula)
        gas.plot_label = pretty_name(formula)

        # Check if we have data from file 
        if !isfile(fpath)
            # no data => generate stub
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

            # critical set to large value (never supercritical)
            gas.T_crit = fbig
            gas.T_trip = 0.0

        else
            # have data => load from file 
            ds = Dataset(fpath,"r")

            # scalar properties  
            gas.mmw = ds["mmw"][1]
            gas.T_trip = ds["T_trip"][1]
            gas.T_crit = ds["T_crit"][1]
            gas.JANAF_name = ds["JANAF"][1]

            # variable properties 
            gas.cap_T = ds["cap_T"][:]
            gas.cap_C = ds["cap_C"][:]

            gas.lat_T = ds["lat_T"][:]
            gas.lat_H = ds["lat_H"][:]

            gas.sat_T = ds["sat_T"][:]
            gas.sat_P = ds["sat_P"][:]

            # close file 
            close(ds)

        end 

        return gas
    end # end load_gas 

    """
    Helper function to find the nearest index in an array to a value 
    """
    function _findnearest(val::Float64, arr:Array{Float64,1})::Int
        return findmin(abs.(arr-val))[2]
    end 

    """
    Helper function to find the two elements of an array which bound a value.
    
    Assumes that the array is sorted.
    """
    function _findbounding(val::Float64, arr:Array{Float64,1})::Tuple{Int,Int}

        # Get length of the search array
        l::Int = length(arr)

        # Invalid case 
        if l < 2
            error("Cannot find bounding elements in array of length <2")
        end 

        # Trivial case 
        if l == 2
            return (1,2)
        end 

        # Closest element
        i::Int = _findnearest(val, arr)

        # Handle cases where we are at either end of the array
        if i == 1
            return (1,2)
        elseif i == l
            return (l-1,l)
        end 

        # Check either side
        if arr[i] > val 
            # left side 
            return (i-1, i)
        else 
            # right side 
            return (i, i+1)
        end 

    end 

    """
    **Get gas saturation pressure for a given temperature.**

    If temperature is ommitted, then the triple point is used. If the
    temperature is above the critical point, then a large value is returned.

    Arguments:
    - `gas::Gas_t`              the gas struct to be used
    - `t::Float64`              temperature [K]

    Returns:
    - `p::Float64`              saturation pressure [Pa]
    """
    function get_Psat(gas::Gas_t, t::Float64=-1.0)::Float64 

        # Above critical point. In practice, a check for this should be made 
        #    before any attempts to evaluate this function.
        if t > gas.T_crit
            return fbig 
        end 

        # Temperature not provided 
        if t < 0.0
            t = gas.T_trip + 1.0e-2
        end 

        # Find nearest lookup points 
        b::Tuple{Int, Int} = _findbounding(gas.sat_T, t)

        # Interpolate between points 
        out::Float64 = ()/(gas.sat_Tb[1]-b[0])
    end 
    
end # end module 
