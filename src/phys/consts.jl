# This file is part of AGNI. License is GPL-3.0: https://www.gnu.org/licenses

"""
**Module for defining physical and numerical constants.**
"""
module consts

    # Misc constants
    const UNSET_STR::String = "__AGNI_UNSET_STR"
    export UNSET_STR

    # Code versions
    const AGNI_VERSION::String     = "1.10.0"  # current agni version
    export AGNI_VERSION
    const SOCVER_minimum::Float64  = 2407.2    # minimum required socrates version
    export SOCVER_minimum

    # A large floating point number
    const BIGFLOAT::Float64     = 1e99
    export BIGFLOAT
    const BIGLOGFLOAT::Float64  = 99.0
    export BIGLOGFLOAT

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

    # Length of an hour in seconds
    const u_hour_sec::Float64 = 3600.0
    export u_hour_sec

    # Length of a day in seconds
    const u_day_sec::Float64 = 24.0 * u_hour_sec
    export u_day_sec

    # Length of a year in seconds
    const u_year_sec::Float64 = 365.25 * u_day_sec
    export u_year_sec

    # Specific heat capacity for ideal gas [J mol-1 K-1]
    const Cp_ideal::Float64 = R_gas * 7/2 # assuming that it is diatomic
    export Cp_ideal

    # Proton mass [kg]
    const proton_mass::Float64 = 1.67262192595e-27 # NIST CODATA
    export proton_mass

    # List of elements included in the model
    const elems_standard::Array{String,1} = ["H","D","C","N","O","S","P",
                                                "He","Ar","F","Cl","Br",
                                                "Ti", "V", "Fe", "Si","Al", "Cr",
                                                "Mg","Ca","Na","Li","K"]
    export elems_standard

    # Standard species
    const vols_standard::Array{String,1} = [
        # volatile atoms
        "H", "O", "C", "N", "S", "P", "He",
        # basic
        "CH4", "CO2", "CO", "H2", "H2O", "O2", "OH", "O3",
        # carbon
        "H2O2", "C2H6", "C2H4","C2H2",  "CH3", "H2CO", "HO2", "C2", "C3", "C4", "C5",
        "CH3CHO", "CH3OOH", "CH3COCH3", "CH3COCHO", "CHOCHO",
        "C2H5CHO", "HOCH2CHO", "C2H5COCH3", "CH3ONO2", "C2H3", "C3H4", "C4H3",
        "C2N2", "HCO",
        # sulfur
        "SO", "S2", "S3", "S6", "S8", "SO2", "SO3", "H2SO4", "H2S", "CS2", "OCS",
        "CH3SH", "CH3S", "C2H6S", "C2H6S2",
        # nitrogen
        "N2", "HCN", "NH3",  "HNO3", "N2O5", "HONO", "HO2NO2",
        "NO3", "N2O", "NO", "NO2", "N2O4", "N2H4", "N2O3", "CN",
        # halogens and biosignatures
        "HCl", "HF", "CH3Cl", "CH3F", "CH3Br", "SF6",
        # phosphorous
        "PH3", "PS", "PO", "PN",
        # isotopologues
        # "HDO",
    ]
    const vaps_standard::Array{String,1} = [
        # (semi)refractory atoms
        "Na", "Si", "Ti", "V", "Mg", "K", "Fe", "Li", "Ca", "Al", "Cr",
        # ions
        # "H-",
        # rock vapours
        "SiO2", "SiO", "SiH", "SiH2", "SiH4",
        "FeO", "FeH",
        "Mg2", "MgO",
        "TiO", "TiO2", "VO", "CrH",
        "CaO", "AlO", "Na2", "NaO", "NaOH", "KOH",
        "HAlO2"
    ]
    const gases_standard::Array{String, 1} = vcat(vols_standard, vaps_standard)
    export vols_standard
    export vaps_standard
    export gases_standard

    # Solar metallicities taken from FastChem source files
    # These are log10 molar (number) ratios relative to hydrogen, offset by 12 dex
    # https://www.aanda.org/articles/aa/pdf/2021/09/aa40445-21.pdf
    const solar_metallicity::Dict{String, Float64} = Dict([
        ("Al" , 6.43),
        ("Ar" , 6.38),
        ("C"  , 8.46),
        ("Ca" , 6.30),
        ("Cl" , 5.31),
        ("Co" , 4.94),
        ("Cr" , 5.62),
        ("Cu" , 4.18),
        ("F"  , 4.40),
        ("Br" , -9.0),  # no value, placeholder
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
