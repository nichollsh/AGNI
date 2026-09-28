# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0

"""
**Module for defining physical and numerical constants.**
"""
module consts

    # Misc constants
    const UNSET_STR::String = "__AGNI_UNSET_STR"
    export UNSET_STR

    # Code versions
    const AGNI_VERSION::String    = "1.13.0"  # current agni version
    export AGNI_VERSION
    const SOCVER_minimum::String  = "2603.9"    # minimum required socrates version
    export SOCVER_minimum

    # A large floating point number
    const BIGFLOAT::Float64     = floatmax(Float32)/2.0
    export BIGFLOAT
    const BIGLOGFLOAT::Float64  = log(BIGFLOAT)
    export BIGLOGFLOAT

    # A small positive floating point number
    const SMALLFLOAT::Float64   = floatmin(Float32)*2.0
    export SMALLFLOAT
    const SMALLLOGFLOAT::Float64 = log(SMALLFLOAT)
    export SMALLLOGFLOAT

    # Sources:
    # - Pierrehumbert (2010) textbook "Principles of Planetary Climate"
    # - NIST CODATA (https://physics.nist.gov/cuu/Constants/index.html)
    # - Sources outlined in Thermo files (https://github.com/nichollsh/ThermoTools)

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

    # Earth radius [m]
    const R_earth::Float64 = 6378137.0 # m
    export R_earth

    # Proton mass [kg]
    const proton_mass::Float64 = 1.67262192595e-27 # NIST CODATA
    export proton_mass

    # List of elements included in the model
    const elems_standard::Array{String,1} = ["H","D","C","N","O","S","P",
                                                "He","Ne","Ar","Kr","Xe",
                                                "F","Cl","Br",
                                                "Ti", "V", "Fe", "Si","Al", "Cr",
                                                "Mg","Ca","Na","Li","K"]
    export elems_standard

    # Standard species
    const vols_standard::Array{String,1} = [
        # volatile atoms
        "H", "O", "C", "N", "S", "P",
        # Noble gases
        "He", "Ne", "Ar", "Kr", "Xe",
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
        "Na", "Si", "Ti", "V", "Mg", "K", "Fe", "Ca", "Al",
        # ions
        # "H-",
        # rock vapours
        "SiO2", "SiO", "SiH", "SiH2", "SiH4",
        "FeO", "FeH", "FeO2H2",
        "Mg2", "MgO",
        "TiO", "TiO2", "VO", "CrH",
        "CaO", "AlO", "Na2", "NaO", "NaOH", "KOH",
        "HAlO2",
        # aerosol and cloud condensates
        "Al2O3", "CaTiO3", "Fe2O3", "Fe2SiO4", "FeS", "MgSiO3", "Mg2SiO4",
        "Na2S", "NaCl", "KCl",
    ]

    const gases_standard::Array{String, 1} = vcat(vols_standard, vaps_standard)
    export vols_standard
    export vaps_standard
    export gases_standard

    # Solar abundances of each element, relative to hydrogen.
    # These are log10 molar (number) ratios relative to hydrogen, offset by 12 dex
    # Current source is Table 3 from: https://arxiv.org/pdf/2608.23155
    # Alternative source: https://www.aanda.org/articles/aa/pdf/2021/09/aa40445-21.pdf
    const solar_metallicity::Dict{String, Float64} = Dict([
        ("H"  , 12.0),
        ("He" , 10.934),
        ("Li" , 0.96),
        ("Be" , 1.21),
        ("B"  , 2.7),
        ("C"  , 8.47),
        ("N"  , 7.84),
        ("O"  , 8.7),
        ("F"  , 4.4),
        ("Ne" , 8.07),
        ("Na" , 6.22),
        ("Mg" , 7.56),
        ("Al" , 6.43),
        ("Si" , 7.55),
        ("P"  , 5.35),
        ("S"  , 7.06),
        ("Cl" , 5.31),
        ("Ar" , 6.38),
        ("K"  , 5.07),
        ("Ca" , 6.3),
        ("Sc" , 3.13),
        ("Ti" , 4.97),
        ("V"  , 3.9),
        ("Cr" , 5.65),
        ("Mn" , 5.47),
        ("Fe" , 7.46),
        ("Co" , 4.94),
        ("Ni" , 6.2),
        ("Cu" , 4.18),
        ("Zn" , 4.56),
        ("Ga" , 3.02),
        ("Ge" , 3.62),
        ("As" , 2.29),
        ("Se" , 3.32),
        ("Br" , 2.55),
        ("Kr" , 3.19),
        ("Rb" , 2.34),
        ("Sr" , 2.84),
        ("Y"  , 2.29),
        ("Zr" , 2.59),
        ("Nb" , 1.47),
        ("Mo" , 1.88),
        ("Tc" , -12.0), # no value, placeholder
        ("Ru" , 1.75),
        ("Rh" , 1.06),
        ("Pd" , 1.57),
        ("Ag" , 1.15),
        ("Cd" , 1.66),
        ("In" , 0.8),
        ("Sn" , 2.02),
        ("Sb" , 1.01),
        ("Te" , 2.14),
        ("I"  , 1.66),
        ("Xe" , 2.26),
        ("Cs" , 1.03),
        ("Ba" , 2.27),
        ("La" , 1.12),
        ("Ce" , 1.58),
        ("Pr" , 0.75),
        ("Nd" , 1.42),
        ("Pm" , -12.0), # no value, placeholder
        ("Sm" , 0.95),
        ("Eu" , 0.49),
        ("Gd" , 1.08),
        ("Tb" , 0.31),
        ("Dy" , 1.1),
        ("Ho" , 0.48),
        ("Er" , 0.93),
        ("Tm" , 0.11),
        ("Yb" , 0.85),
        ("Lu" , 0.1),
        ("Hf" , 0.86),
        ("Ta" , -0.14),
        ("W"  , 0.79),
        ("Re" , 0.29),
        ("Os" , 1.31),
        ("Ir" , 1.32),
        ("Pt" , 1.6),
        ("Au" , 0.91),
        ("Hg" , 1.04),
        ("Tl" , 0.92),
        ("Pb" , 1.95),
        ("Bi" , 0.61),
        ("Po" , -12.0), # no value, placeholder
        ("At" , -12.0), # no value, placeholder
        ("Rn" , -12.0), # no value, placeholder
        ("Fr" , -12.0), # no value, placeholder
        ("Ra" , -12.0), # no value, placeholder
        ("Ac" , -12.0), # no value, placeholder
        ("Th" , 0.03),
        ("Pa" , -12.0), # no value, placeholder
        ("U"  , -0.52),
    ])

    # Table of condensed-phase density for ocean and aerosol calculations [kg/m^3]
    const _lookup_rho::Dict{String, Float64} = Dict([
        # https://encyclopedia.airliquide.com/water#properties
        ("H2O", 958.37 ),  # boiling
        ("CO2", 1178.4 ),  # triple
        ("H2" , 70.516 ),  # boiling
        ("CH4", 422.36 ),  # boiling
        ("CO" , 793.2  ),  # boiling
        ("N2" , 806.11 ),  # boiling
        ("NH3", 681.97 ),  # boiling
        ("SO2", 1461.1 ),  # boiling

        # Handbook of Mineralogy (Anthony et al., https://www.handbookofmineralogy.org/)
        ("SiO2",                2660.0),  # alpha quartz
        ("SiO2_amorph",         2201.0),  # fused silica
        ("SiO",                 2140.0),  # amorphous
        ("FeO",                 5970.0),  # wustite
        ("Fe",                  7874.0),  # alpha iron
        ("Fe2O3",               5255.0),  # hematite
        ("FeS",                 4850.0),  # troilite
        ("Fe2SiO4_KH",          4400.0),  # fayalite
        ("MgO",                 3580.0),  # periclase
        ("MgSiO3",              3189.0),  # enstatite
        ("MgSiO3_amorph_glass", 2710.0),  # glass
        ("Mg2SiO4_crystalline", 3271.0),  # forsterite
        ("CaTiO3_KH",           4020.0),  # perovskite (synthetic)
        ("TiO2_anatase",        3890.0),  # anatase
        ("TiO2_rutile",         4250.0),  # rutile
        ("Na2S",                1856.0),  # PubChem CID 14804 (Merck Index)
        ("KCl",                 1987.0),  # sylvite
        ("NaCl",                2165.0),  # halite
        ("C",                   2260.0),  # graphite

        # https://webmineral.com/data/Corundum.shtml
        ("Al2O3", 4050.0), # Corundum

        # https://www.webmineral.com/data/Forsterite.shtml
        ("Mg2SiO4", 3270.0), # Forsterite

        # Bond & Bergstrom (2006), Aerosol Sci. Technol. 40, 27
        ("Soot",                1800.0),

        # Imanaka et al. (2012) via Horst & Tolbert (2013)
        ("Tholin",              1350.0),  # 1.3-1.4 g/cm3
    ])

end
