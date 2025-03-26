# Contains physical constants

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module consts
    # Sources:
    # - Pierrehumbert (2010)
    # - NIST
    # - Sources outlined in Thermo files (https://github.com/nichollsh/Thermo)

    # Universal gas constant, J K-1 mol-1
    const R_gas::Float64 = 8.314462618 # NIST CODATA
    export R_gas

    # Stefan-boltzmann constant, W m−2 K−4
    const σSB::Float64 =  5.670374419e-8 # NIST CODATA
    export σSB

    # Von Karman constant, dimensionless
    const k_vk::Float64 = 0.40  # Hogstrom 1988
    export k_vk

    # Planck's constant, J s
    const h_pl::Float64 = 6.62607015e-34  # NIST CODATA
    export h_pl

    # Boltzmann constant, J K-1
    const k_B::Float64 = 1.380649e-23 # NIST CODATA
    export k_B

    # Speed of light, m s-1
    const c_vac::Float64 = 299792458.0 # NIST CODATA
    export c_vac

    # 0 degrees celcius, in kelvin
    const zero_celcius::Float64 = 273.15
    export zero_celcius

    # Mixing length parameter, dimensionless
    const αMLT::Float64 = 1.0
    export αMLT

    # Newton's gravitational constant [m3 kg-1 s-2]
    const G_grav::Float64 = 6.67430e-11 # NIST CODATA
    export G_grav

    # List of elements included in the model
    const elems_standard::Array{String,1} = ["H","C","N","O","S","P",
                                                "Fe","Mg","Si","Ca","Al"]
    export elems_standard

    # Standard species
    const gases_standard::Array{String, 1} = [
        "C2H4", "CO", "H2O", "H2SO4", "N2", "O2", "OH", "SO", "SiO", "CH4", "CO2", "H2",
        "H2S", "HCN", "NH3", "OCS", "S2", "S8", "SO2", "SiO2", "O3", "N2O", "NO", "NO2",
        "HNO3", "FeH", "PH3", "C2H2", "NO3", "N2O5", "HONO", "HO2NO2", "H2O2", "C2H6",
        "CH3", "H2CO", "HO2", "C", "Fe", "FeO", "H", "Mg2", "Mg", "MgO", "N", "O", "Si",
        "S", "SO", "CS2"
    ]
    export gases_standard

    # Table of gas molecule mean molecular weights, kg mol-1
    const _lookup_mmw::Dict{String, Float64} = Dict([
        # molecules
        ("H2O",     1.801530E-02 ),
        ("CO2",     4.401000E-02 ),
        ("O3",      4.799820E-02 ),
        ("N2O",     4.401280E-02 ),
        ("CO",      2.801060E-02 ),
        ("CH4",     1.604300E-02 ),
        ("O2",      3.199880E-02 ),
        ("NO",      3.000610E-02 ),
        ("SO2",     6.406280E-02 ),
        ("NO2",     4.600550E-02 ),
        ("NH3",     1.703060E-02 ),
        ("HNO3",    6.301290E-02 ),
        ("N2",      2.801340E-02 ),
        ("TiO",     6.386600E-02 ),
        ("VO",      6.694090E-02 ),
        ("H2",      2.015880E-03 ),
        ("OCS",     6.007500E-02 ),
        ("FeH",     5.685300E-02 ),
        ("CrH",     5.300400E-02 ),
        ("PH3",     3.399758E-02 ),
        ("C2H2",    2.603730E-02 ),
        ("HCN",     2.702530E-02 ),
        ("H2S",     3.408100E-02 ),
        ("NO3",     6.301280E-02 ),
        ("N2O5",    1.080104E-01 ),
        ("HONO",    4.701340E-02 ),
        ("HO2NO2",  7.901220E-02 ),
        ("H2O2",    3.401470E-02 ),
        ("C2H6",    3.006900E-02 ),
        ("CH3",     1.503450E-02 ),
        ("H2CO",    3.002600E-02 ),
        ("HO2",     3.300670E-02 ),
        ("HDO",     1.902140E-02 ),
        ("HCl",     3.646100E-02 ),
        ("HF",      2.000689E-02 ),

        # elements from https://iupac.qmul.ac.uk/AtWt/
        ("H",   1.008000000e-03 ),
        ("He",  4.002000000e-03 ),
        ("Li",  6.940000000e-03 ),
        ("Be",  9.012000000e-03 ),
        ("B",   1.081000000e-02 ),
        ("C",   1.201100000e-02 ),
        ("N",   1.400700000e-02 ),
        ("O",   1.599900000e-02 ),
        ("F",   1.899800000e-02 ),
        ("Ne",  2.017970000e-02 ),
        ("Na",  2.298900000e-02 ),
        ("Mg",  2.430500000e-02 ),
        ("Al",  2.698100000e-02 ),
        ("Si",  2.808500000e-02 ),
        ("P",   3.097300000e-02 ),
        ("S",   3.206000000e-02 ),
        ("Cl",  3.545000000e-02 ),
        ("Ar",  3.995000000e-02 ),
        ("K",   3.909830000e-02 ),
        ("Ca",  4.007800000e-02 ),
        ("Sc",  4.495500000e-02 ),
        ("Ti",  4.786700000e-02 ),
        ("V",   5.094150000e-02 ),
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
        ("Y",   8.890500000e-02 ),
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
        ("I",   1.269040000e-01 ),
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
        ("W",   1.838400000e-01 ),
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
        ("U",   2.380280000e-01),
    ])

    # Table of pre-defined colors for plotting
    const _lookup_color::Dict{String, String} = Dict([
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
        ("H"  , "#0000aa"),
        ("C"  , "#ff0000"),
        ("O"  , "#00cc00"),
        ("N"  , "#ffaa00"),
        ("S" ,  "#ff22ff"),
        ("P" ,  "#33ccff"),
        ("He" , "#30FF71" ),

        # refractory elements
        ("Fe" , "#888888"),
        ("Si" , "#aa2277"),
        ("Mg" , "#996633"),
        ("Na" , "#1144ff")
    ])

end
