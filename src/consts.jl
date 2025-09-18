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

    # Convection beta parameter for obtaining Ra, dimensionless
    const βRa::Float64 = 0.3
    export βRa

    # Newton's gravitational constant [m3 kg-1 s-2]
    const G_grav::Float64 = 6.67430e-11 # NIST CODATA
    export G_grav

    # Specific heat capacity for ideal gas [J mol-1 K-1]
    const Cp_ideal::Float64 = R_gas * 7/2 # assuming that it is diatomic
    export Cp_ideal

    # Proton mass [kg]
    const proton_mass::Float64 = 1.67262192595e-27 # NIST CODATA
    export proton_mass

    # List of elements included in the model
    const elems_standard::Array{String,1} = ["H","C","N","O","S","P","He",
                                                "Fe","Mg","Si","Ca","Al"]
    export elems_standard

    # Standard species
    const gases_standard::Array{String, 1} = [
        "C2H4", "CO", "H2O", "H2SO4", "N2", "O2", "OH", "SO", "SiO", "CH4", "CO2", "H2",
        "H2S", "HCN", "NH3", "OCS", "S2", "S8", "SO2", "SiO2", "O3", "N2O", "NO", "NO2",
        "HNO3", "FeH", "PH3", "C2H2", "NO3", "N2O5", "HONO", "HO2NO2", "H2O2", "C2H6",
        "CH3", "H2CO", "HO2", "C", "Fe", "FeO", "H", "Mg2", "Mg", "MgO", "N", "O", "Si",
        "S", "SO", "CS2","He",
    ]
    export gases_standard

    # Atom counts in some standard species
    const _lookup_count_atoms::Dict{String, Dict} = Dict([
        "H2O" => Dict("H"=>2, "O"=>1),
        "H2"  => Dict("H"=>2),
        "O2"  => Dict("O"=>2),
        "O3"  => Dict("O"=>3),
        "CO"  => Dict("C"=>1, "O"=>1),
        "CO2" => Dict("C"=>1, "O"=>2),
        "CH4" => Dict("C"=>1, "H"=>4),
        "N2"  => Dict("N"=>2),
        "NH3" => Dict("N"=>1, "H"=>3),
        "SO2" => Dict("S"=>1, "O"=>2),
        "H2S" => Dict("S"=>1, "H"=>2),
        "S2"  => Dict("S"=>2),
    ]
    )
    export _lookup_count_atoms

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

    # Solar metallicities taken from FastChem source files (Asplund et al. 2020, A&A)
    const _solar_metallicity::Dict{String, Float64} = Dict([
        ("Al" , 6.43),
        ("Ar" , 6.38),
        ("C"  , 8.46),
        ("Ca" , 6.30),
        ("Cl" , 5.31),
        ("Co" , 4.94),
        ("Cr" , 5.62),
        ("Cu" , 4.18),
        ("F"  , 4.40),
        ("Fe" , 7.46),
        ("Ge" , 3.62),
        ("H"  , 12.00),
        ("He" , 10.914),
        ("K"  , 5.07),
        ("Mg" , 7.55),
        ("Mn" , 5.42),
        ("N"  , 7.83),
        ("Na" , 6.22),
        ("Ne" , 8.06),
        ("Ni" , 6.20),
        ("O"  , 8.69),
        ("P"  , 5.41),
        ("S"  , 7.12),
        ("Si" , 7.51),
        ("Ti" , 4.97),
        ("V"  , 3.90),
        ("Zn" , 4.56),
        ("Li" , 0.96),
        ("Be" , 1.38),
        ("B"  , 2.70),
        ("Sc" , 3.14),
        ("Ga" , 3.02),
        ("As" , 2.30),
        ("Se" , 3.34),
        ("Rb" , 2.32),
        ("Sr" , 2.83),
        ("Y"  , 2.21),
        ("Zr" , 2.59),
        ("Nb" , 1.47),
        ("Mo" , 1.88),
        ("Ru" , 1.75),
        ("Rh" , 0.78),
        ("Pd" , 1.57),
        ("Ag" , 0.96),
        ("Cd" , 1.71),
        ("In" , 0.80),
        ("Sn" , 2.02),
        ("Sb" , 1.01),
        ("Te" , 2.18),
        ("Cs" , 1.08),
        ("Ba" , 2.27),
        ("La" , 1.11),
        ("Ce" , 1.58),
        ("Pr" , 0.75),
        ("Nd" , 1.42),
        ("Sm" , 0.95),
        ("Eu" , 0.52),
        ("Gd" , 1.08),
        ("Tb" , 0.31),
        ("Dy" , 1.10),
        ("Ho" , 0.48),
        ("Er" , 0.93),
        ("Tm" , 0.11),
        ("Yb" , 0.85),
        ("Lu" , 0.10),
        ("Hf" , 0.85),
        ("Ta" , -0.15),
        ("W"  , 0.79),
        ("Re" , 0.26),
        ("Os" , 1.35),
        ("Ir" , 1.32),
        ("Pt" , 1.61),
        ("Au" , 0.91),
        ("Hg" , 1.17),
        ("Tl" , 0.92),
        ("Pb" , 1.95),
        ("Bi" , 0.65),
        ("Th" , 0.03),
        ("U" , -0.54),
    ])

end
