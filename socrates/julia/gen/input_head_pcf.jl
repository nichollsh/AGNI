# Autogenerated from src/modules_gen/input_head_pcf.f90
# svn revision 1226

module input_head_pcf



const NP_max_length_line = 132

const NPD_in_profile = 9

const NPD_data_column = 35

const NPD_phys_type = 62

const NPD_unit = 20

const len_col_header = 10

const len_long_title = 30

const len_file_suffix = 12

const IP_physical_data = 1

const IP_gaseous_data = 2

const IP_aerosol_data = 3

const IP_pressure = 1

const IP_height = 2

const IP_temperature = 3

const IP_ozone_density = 4

const IP_water_vapour_density = 5

const IP_lwc = 6

const IP_re = 7

const IP_null_field = 8

const IP_temperature_ground = 9

const IP_cloud_fraction = 10

const IP_flux_up = 11

const IP_flux_down = 12

const IP_flux_direct = 13

const IP_flux_net = 14

const IP_heating_rate = 15

const IP_flux_total_down = 16

const IP_lwm = 17

const IP_t_dew = 18

const IP_spec_humidity = 19

const IP_hum_mix_ratio = 20

const IP_iwc = 21

const IP_iwm = 22

const IP_ire = 23

const IP_tau = 24

const IP_ssa = 25

const IP_gsc = 26

const IP_fsc = 27

const IP_planck = 28

const IP_cloud_fraction_cv = 29

const IP_lwccv = 30

const IP_lwmcv = 31

const IP_recv = 32

const IP_iwccv = 33

const IP_iwmcv = 34

const IP_irecv = 35

const IP_pressure_ground = 36

const IP_t_level = 37

const IP_rel_hum = 38

const IP_discard = 39

const IP_cloud_fraction_w = 40

const IP_cloud_fraction_i = 41

const IP_cloud_fraction_w_cv = 42

const IP_cloud_fraction_i_cv = 43

const IP_solar_zenith_angle = 44

const IP_solar_azimuth = 45

const IP_solar_toa = 46

const IP_surface_type = 47

const IP_view_polar = 48

const IP_view_azim = 49

const IP_radiance = 50

const IP_surface_char = 51

const IP_optical_water = 52

const IP_optical_ice = 53

const IP_optical_ss = 54

const IP_iso_inc = 55

const IP_view_geom = 56

const IP_photolysis = 57

const IP_p_level = 58

const IP_contrib_funci = 59

const IP_contrib_funcf = 60

const IP_actinic_flux = 61

const IP_photolysis_rate = 62

const header_phys = [
    "PRESS",
    "HGT",
    "TEMP",
    "O3DEN",
    "H2ODEN",
    "LWC",
    "RE",
    "",
    "TSTAR",
    "CLFRAC",
    "UFLX",
    "DFLX",
    "SFLX",
    "NFLX",
    "HRTS",
    "VFLX",
    "LWM",
    "TDEW",
    "SPH",
    "HMR",
    "IWC",
    "IWM",
    "IRE",
    "TAU",
    "SSA",
    "ASYM",
    "FRWSC",
    "PLANCK",
    "CVFRAC",
    "LWCCV",
    "LWMCV",
    "RECV",
    "IWCCV",
    "IWMCV",
    "IRECV",
    "PSTAR",
    "TLEV",
    "RH",
    "DISCRD",
    "WCLFRC",
    "ICLFRC",
    "CWFRAC",
    "CIFRAC",
    "SZEN",
    "SAZIM",
    "STOA",
    "SURF",
    "POLAR",
    "AZIM",
    "RADN",
    "SRFCHR",
    "OPWT",
    "OPICE",
    "OPSS",
    "ISOS",
    "GEOM",
    "PHOTOL",
    "PLEV",
    "CFI",
    "CFF",
    "AFLX",
    "PHRATE",
]

const phys_suffix = [
    "p",
    "hgt",
    "t",
    "o3d",
    "qd",
    "lwc",
    "re",
    "null",
    "tstar",
    "clfr",
    "uflx",
    "dflx",
    "sflx",
    "nflx",
    "hrts",
    "vflx",
    "lwm",
    "tdw",
    "q",
    "hmr",
    "iwc",
    "iwm",
    "ire",
    "tau",
    "ssa",
    "gsc",
    "fsc",
    "plk",
    "ccfr",
    "lwccv",
    "lwmcv",
    "recv",
    "iwccv",
    "iwmcv",
    "irecv",
    "pstar",
    "tl",
    "rh",
    "",
    "wclfr",
    "iclfr",
    "wccfr",
    "iccfr",
    "szen",
    "sazim",
    "stoa",
    "surf",
    "vwpol",
    "vwazim",
    "radn",
    "surf",
    "op_water",
    "op_ice",
    "ss",
    "isos",
    "view",
    "photol",
    "pl",
    "cfi",
    "cff",
    "aflx",
    "ph_rate",
]

const phys_title = [
    " Pressure",
    " Height",
    " Temperature",
    " Ozone Mass Density",
    " Water Vapour Mass Density",
    " Liquid Water Content",
    " Effective Radius",
    " Null Grid",
    " Temperature of surface",
    " Cloud Fraction",
    " Upward Flux",
    " Downward Flux (diffuse)",
    " Direct Flux",
    " Net Downward Flux",
    " Heating Rates (K, day)",
    " Total Downward Flux",
    " Liquid Water Mass Fraction",
    " Dew point Temperature",
    " Specific Humidity",
    " Humidity Mixing Ratio",
    " Ice Content",
    " Ice Mass Fraction",
    " Ice Effective Radius",
    " Optical Depth",
    " Albedo of Single Scattering",
    " Asymmetry Factor",
    " Forward Scattering Factor",
    " Planck Function",
    " Convective Cloud Fraction",
    " Conv. Liquid Water Content",
    " Conv. Liquid Water Mass Frac",
    " Conv. Effective Radius",
    " Convective Ice Content",
    " Convective Ice Mass Fraction",
    " Conv. Ice Effective Radius",
    " Pressure at surface",
    " Temperature on Levels",
    " Relative Humidity",
    "",
    " Water Cloud Fraction",
    " Ice Cloud Fraction",
    " Conv. Water Cloud Fraction",
    " Conv. Ice Cloud Fraction",
    " Solar zenith angle",
    " Solar azimuthal angle",
    " Solar Irradiance at TOA",
    " Type of surface",
    " Polar viewing angle",
    " Azimuthal viewing angle",
    " Radiance",
    " Surface characteristics",
    " Optical data for droplets",
    " Optical data for ice crystal",
    " Single scattering properties",
    " Isotropic source",
    " Viewing Geometry",
    " Rate of photolysis",
    " Pressure on Levels",
    " Contribution function (inty)",
    " Contribution function (flux)",
    " Actinic flux",
    " Photolysis rate",
]

const header_gas = [
    "H2O",
    "CO2",
    "O3",
    "N2O",
    "CO",
    "CH4",
    "O2",
    "NO",
    "SO2",
    "NO2",
    "NH3",
    "HNO3",
    "N2",
    "CFC11",
    "CFC12",
    "CFC113",
    "HCFC22",
    "HFC125",
    "HFC134A",
    "CFC114",
    "TiO",
    "VO",
    "H2",
    "He",
    "OCS",
    "Na",
    "K",
    "FeH",
    "CrH",
    "Li",
    "Rb",
    "Cs",
    "PH3",
    "C2H2",
    "HCN",
    "H2S",
    "Ar",
    "AIR",
    "O",
    "N",
    "NO3",
    "N2O5",
    "HONO",
    "HO2NO2",
    "H2O2",
    "C2H6",
    "CH3",
    "H2CO",
    "HO2",
    "HDO",
    "HCl",
    "HF",
    "cOSSO",
    "tOSSO",
    "yOSOS",
]

const gas_suffix = [
    "q",
    "co2",
    "o3",
    "n2o",
    "co",
    "ch4",
    "o2",
    "no",
    "so2",
    "no2",
    "nh3",
    "hno3",
    "n2",
    "cfc11",
    "cfc12",
    "cfc113",
    "hcfc22",
    "hfc125",
    "hfc134a",
    "cfc114",
    "tio",
    "vo",
    "h2",
    "he",
    "ocs",
    "na",
    "k",
    "feh",
    "crh",
    "li",
    "rb",
    "cs",
    "ph3",
    "c2h2",
    "hcn",
    "h2s",
    "ar",
    "air",
    "o",
    "n",
    "no3",
    "n2o5",
    "hono",
    "ho2no2",
    "h2o2",
    "c2h6",
    "ch3",
    "h2co",
    "ho2",
    "hdo",
    "hcl",
    "hf",
    "cosso",
    "tosso",
    "yosos",
]

const gas_title = [
    " Water Vapour",
    " Carbon Dioxide",
    " Ozone",
    " Dinitrogen oxide",
    " Carbon Monoxide",
    " Methane",
    " Oxygen",
    " Nitrogen Monoxide",
    " Sulphur Dioxide",
    " Nitrogen Dioxide",
    " Ammonia",
    " Nitric acid",
    " Nitrogen",
    " CFC-11",
    " CFC-12",
    " CFC-113",
    " HCFC-22",
    " HFC-125",
    " HFC-134a",
    " CFC-114",
    " Titanium oxide",
    " Vanadium oxide",
    " Hydrogen",
    " Helium",
    " Carbonyl sulphide",
    " Sodium",
    " Potassium",
    " Iron hydride",
    " Chromium hydride",
    " Lithium",
    " Rubidium",
    " Cesium",
    " Phosphine",
    " Acetylene",
    " Hydrogen cyanide",
    " Hydrogen sulphide",
    " Argon",
    " Dry air",
    " Atomic Oxygen",
    " Atomic Nitrogen",
    " Nitrate radical",
    " Dinitrogen Pentoxide",
    " Nitrous acid",
    " Peroxynitric acid",
    " Hydrogen peroxide",
    " Ethane",
    " Methyl radical",
    " Formaldehyde",
    " Hydroperoxy radical",
    " Semiheavy water",
    " Hydrogen chloride",
    " Hydrogen fluoride",
    " cis-OSSO",
    " trans-OSSO",
    " OSO-S",
]

const header_aerosol = [
    "WTSOL",
    "DUST",
    "OCN",
    "SOOT",
    "ASH",
    "SULPH",
    "NH4SO4",
    "AUNCH",
    "SAHARA",
    "ACCUM",
    "AITKEN",
    "FRSOOT",
    "AGSOOT",
    "NACL",
    "NACLFLM",
    "NACLJET",
    "DUSTDIV1",
    "DUSTDIV2",
    "DUSTDIV3",
    "DUSTDIV4",
    "DUSTDIV5",
    "DUSTDIV6",
    "BIOMS1",
    "BIOMS2",
    "BIOGENIC",
    "FROCFF",
    "AGOCFF",
    "DELTA",
    "",
    "NITRATE",
    "DUST2BIN1",
    "DUST2BIN2",
]

const aerosol_suffix = [
    "wtsol",
    "dust",
    "ocn",
    "soot",
    "ash",
    "sulph",
    "nh4so4",
    "aunch",
    "sahara",
    "accum",
    "aitken",
    "frsoot",
    "agsoot",
    "nacl",
    "naclflm",
    "nacljet",
    "dustdiv1",
    "dustdiv2",
    "dustdiv3",
    "dustdiv4",
    "dustdiv5",
    "dustdiv6",
    "bioms1",
    "bioms2",
    "biogenic",
    "frocff",
    "agocff",
    "delta",
    "",
    "nitrate",
    "dust2bin1",
    "dust2bin2",
]

const aerosol_opt_suffix = [
    "op_wtsol",
    "op_dust",
    "op_ocn",
    "op_soot",
    "op_ash",
    "op_sulph",
    "op_nh4so4",
    "op_aunch",
    "op_sahara",
    "op_accum",
    "op_aitken",
    "op_frsoot",
    "op_agsoot",
    "op_nacl",
    "op_naclflm",
    "op_nacljet",
    "op_dustdiv1",
    "op_dustdiv2",
    "op_dustdiv3",
    "op_dustdiv4",
    "op_dustdiv5",
    "op_dustdiv6",
    "op_bioms1",
    "op_bioms2",
    "op_biogenic",
    "op_frocff",
    "op_agocff",
    "op_delta",
    "",
    "op_nitrate",
    "op_dust2bin",
    "op_dust2bin",
]

const aerosol_title = [
    " Water-soluble Aerosol",
    " Dust-like Aerosol",
    " Oceanic Aerosol",
    " Soot Aerosol",
    " Ash Aerosol",
    " 75% Sulphuric acid Aerosol",
    " Ammonium Sulphate Aerosol",
    " Uncharacterized Aerosol",
    " Saharan Dust Aerosol",
    " Accumulation SO4 Aerosol",
    " Aitken-mode SO4 Aerosol",
    " Fresh Soot Aerosol",
    " Aged Soot Aerosol",
    " Generic Sodium Chloride",
    " Sodium Chloride (Film mode)",
    " Sodium Chloride (Jet mode)",
    " Dust aerosol (division 1)",
    " Dust aerosol (division 2)",
    " Dust aerosol (division 3)",
    " Dust aerosol (division 4)",
    " Dust aerosol (division 5)",
    " Dust aerosol (division 6)",
    " Biomass aerosol (division 1)",
    " Biomass aerosol (division 2)",
    " Biogenic aerosol",
    " Fresh fossil-fuel org. carbo",
    " Aged fossil-fuel org. carbon",
    " Unspecified (delta) aerosol",
    "",
    " Ammonium nitrate aerosol",
    " Two-bin Dust aerosol (div 1)",
    " Two-bin Dust aerosol (div 2)",
]

const IP_unit_pa = 1

const IP_unit_mb = 2

const IP_unit_k = 3

const IP_unit_m = 4

const IP_unit_km = 5

const IP_unit_um = 6

const IP_unit_kgm_3 = 7

const IP_unit_gm_3 = 8

const IP_unit_kgm_2 = 9

const IP_unit_gm_2 = 10

const IP_unit_gkg = 11

const IP_unit_g_g = 12

const IP_unit_none = 13

const IP_unit_m_3 = 14

const IP_unit_cm_3 = 15

const IP_unit_wm_2 = 16

const IP_unit_c = 17

const name_unit = [
    "PA",
    "MB",
    "K",
    "M",
    "KM",
    "UM",
    "KGM-3",
    "GM-3",
    "KGM-2",
    "GM-2",
    "GKG-1",
    "GG-1",
    "NONE",
    "M-3",
    "CM-3",
    "WM-2",
    "C",
    "PPMV",
    "VOL-MARS",
    "VOL-VENUS",
]

const factor_unit = [
    1.0,
    1.0E+02,
    1.0,
    1.0,
    1.0E+03,
    1.0E-06,
    1.0,
    1.0E-03,
    1.0,
    1.0E-03,
    1.0E-03,
    1.0,
    1.0,
    1.0,
    1.0E+06,
    1.0,
    1.0,
    -1.0E-06,
    -1.0,
    -1.0,
]

const offset_unit = [
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    2.7315E+02,
    0.0,
    0.0,
    0.0,
]

end
