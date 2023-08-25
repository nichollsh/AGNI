# 
# Example from https://code.metoffice.gov.uk/trac/socrates/wiki/SocratesDoc/ClearskyFluxExample
# (clear sky fluxes, based on examples/netcdf/CIRC_case6)
# 
# Currently hardcoded to run clear sky, LW flux, no absorption
# (ie expect constant flux_up, zero flux_down)

# Code based on src/aux/l_run_cdf.F90, sbin/Cl_run_cdf

# NB: File suffixes for physical quantities etc are defined in
# src/modules_gen/input_head_pcf.f90


import SOCRATES
import NCDatasets
using Plots

# Edit for path on local computer
SOCRATES_DIR = "../.."

# SOCRATES data and example folders relative to SOCRATES_DIR
RAD_DATA = joinpath(SOCRATES_DIR, "data")

#####################################
# create Fortran Types 
#######################################

dimen = SOCRATES.StrDim()
control = SOCRATES.StrCtrl()
spectrum = SOCRATES.StrSpecData()
atm = SOCRATES.StrAtm()
cld = SOCRATES.StrCld()
aer = SOCRATES.StrAer()
bound = SOCRATES.StrBound()
radout = SOCRATES.StrOut()

#####################################
# Configure options
#####################################

# SOCRATES example to use
example_dir = joinpath(SOCRATES_DIR, "examples/netcdf/CIRC_case6")
base_name = "case6"

# channel for every band
all_channels = true

# output plots
do_plot = true

# deallocate arrays (false to allow inspection from REPL)
# do_deallocate = true
do_deallocate = false

spec_opt = "sp_lw_ga7"
example_opts = ["gas"]

# spec_opt = "sp_sw_ga7"
# example_opts = ["rayleigh", "gas"]

if spec_opt == "sp_lw_ga7"
    control.spectral_file = joinpath(RAD_DATA, "spectra/ga7/sp_lw_ga7")
    control.isolir = SOCRATES.rad_pcf.ip_infra_red
    # n_channel = 1

elseif spec_opt == "sp_sw_ga7"
    control.spectral_file = joinpath(RAD_DATA, "spectra/ga7/sp_sw_ga7")
    control.isolir = SOCRATES.rad_pcf.ip_solar
    # n_channel = 1 
else
    error("no recognized spectral file option $spec_opt")
end

control.l_rayleigh = ("rayleigh" in example_opts)
control.l_gas = ("gas" in example_opts)
control.l_continuum = ("continuum" in example_opts)
surface_albedo = 0.0

##########################
# Initialisation
#############################

control.i_cloud_representation = SOCRATES.rad_pcf.ip_cloud_type_homogen
atm.n_profile = 0
cld.n_condensed = 0

#########################################
# spectral data
#########################################

# SOCRATES.read_spectrum(control.spectral_file, spectrum)
SOCRATES.set_spectrum(spectrum=spectrum, spectral_file=control.spectral_file, l_all_gasses=true)

gas_index = zeros(Int, SOCRATES.gas_list_pcf.npd_gases) # pointers to gases in spectral file
for i in 1:spectrum.Gas.n_absorb
    ti = spectrum.Gas.type_absorb[i]
    gas_index[ti] = i
    println("spectrum.Gas $i type_absorb $ti $(SOCRATES.gas_list_pcf.name_absorb[ti])")
end

#########################################
# input files
#########################################

nc_tstar = NCDatasets.NCDataset(joinpath(example_dir, base_name*".tstar"))
# no .surf file ??
nc_surflw = NCDatasets.NCDataset(joinpath(example_dir, base_name*".surflw"))
nc_t = NCDatasets.NCDataset(joinpath(example_dir, base_name*".t"))
nc_tl = NCDatasets.NCDataset(joinpath(example_dir, base_name*".tl"))

nc_open = [nc_tstar, nc_surflw, nc_t, nc_tl]

#########################################
# diagnostics
#########################################
control.l_actinic_flux = Bool(spectrum.Basic.l_present[2])
control.l_photolysis_rate = spectrum.Photol.n_pathway > 0
control.l_flux_div = spectrum.Photol.n_pathway > 0

#########################################
# set array sizes and allocate arrays
########################################

if all_channels
    n_channel = spectrum.Basic.n_band
else
    n_channel = 1
end

# modules_gen/dimensions_field_cdf_ucf.f90
npd_direction = 1    # Maximum number of directions for radiances
npd_layer = nc_t.dim["plev"]
npd_latitude = nc_t.dim["lat"]
npd_longitude = nc_t.dim["lon"]

npd_column = 24 # Maximum number of cloudy subcolumns
# npd_profile = 131 # Number of atmospheric profiles
npd_profile = npd_latitude
npd_profile == 1 ||
    error("only single column supported: requested $npd_profile")
npd_max_order = 101 #       Maximum order of spherical harmonics used
npd_brdf_basis_fnc = 2 #       Number of BRDF basis functions
npd_brdf_trunc = 5 #       Order of BRDF truncation
npd_profile_aerosol_prsc = 9 # Size allocated for profiles of prescribed aerosol optical properties
npd_profile_cloud_prsc = 9 # Size allocated for profiles of prescribed cloudy optical properties
npd_opt_level_aerosol_prsc = 170 # Size allocated for levels of prescribed aerosol optical properties
npd_opt_level_cloud_prsc = 170 # Size allocated for levels of prescribed cloudy optical properties

# modules_gen/dimensioms_fixed_pcf.f90
npd_cloud_component        =  4 #   Number of components of clouds.
npd_cloud_type             =  4 #   Number of permitted types of clouds.
npd_overlap_coeff          = 18 #   Number of overlap coefficients for cloud
npd_source_coeff           =  2 #   Number of coefficients for two-stream sources
npd_region                 =  3 # Number of regions in a layer


dimen.nd_profile                = npd_profile
dimen.nd_flux_profile           = npd_profile
dimen.nd_2sg_profile            = npd_profile
dimen.nd_radiance_profile       = npd_profile
dimen.nd_j_profile              = 1
dimen.nd_layer                  = npd_layer
dimen.nd_layer_clr              = npd_layer
dimen.id_cloud_top              = 1
dimen.nd_channel                = n_channel
dimen.nd_column                 = npd_column
dimen.nd_max_order              = npd_max_order
dimen.nd_direction              = npd_direction
dimen.nd_viewing_level          = npd_layer
dimen.nd_brdf_basis_fnc         = npd_brdf_basis_fnc
dimen.nd_brdf_trunc             = npd_brdf_trunc
dimen.nd_profile_aerosol_prsc   = npd_profile_aerosol_prsc
dimen.nd_profile_cloud_prsc     = npd_profile_cloud_prsc
dimen.nd_opt_level_aerosol_prsc = npd_opt_level_aerosol_prsc
dimen.nd_opt_level_cloud_prsc   = npd_opt_level_cloud_prsc
dimen.nd_cloud_component        = npd_cloud_component
dimen.nd_cloud_type             = npd_cloud_type
dimen.nd_overlap_coeff          = npd_overlap_coeff
dimen.nd_source_coeff           = npd_source_coeff
dimen.nd_region                 = npd_region
dimen.nd_point_tile             = 1
dimen.nd_tile                   = 1
dimen.nd_subcol_gen             = 1
dimen.nd_subcol_req             = 1
dimen.nd_aerosol_mode           = 1

SOCRATES.allocate_atm(atm, dimen, spectrum)
SOCRATES.allocate_cld(cld, dimen, spectrum)
SOCRATES.allocate_aer(aer, dimen, spectrum)
SOCRATES.allocate_bound(bound, dimen, spectrum)

# atm sizes and coordinates 
atm.n_layer = nc_t.dim["plev"]
atm.n_layer <= dimen.nd_layer ||
    error("nd_layer $(dimen.nd_layer) <= n_layer $(atm.n_layer)")

n_latitude = nc_t.dim["lat"]
n_longitude = nc_t.dim["lon"]
atm.n_profile = n_latitude*n_longitude
atm.n_profile <= dimen.nd_profile ||
    error("nd_profile $(dimen.nd_profile) < n_profile = $(n_latitude)*$(n_longitude)")
atm.lat[1:n_latitude] .= nc_t["lat"][:]
atm.lon[1:n_longitude] .= nc_t["lon"][:]

###########################################
# Spectral region
###########################################


if control.isolir == SOCRATES.rad_pcf.ip_solar
    Bool(spectrum.Basic.l_present[2]) ||
        error("The spectral file contains no solar spectral data.")

    #   Assign the solar zenith angles from the input file.
    #   They will be converted to trigonometric functions later.
    NCDatasets.NCDataset(joinpath(example_dir, base_name*".szen")) do nc_szen
        nc_szen.dim["lat"]*nc_szen.dim["lon"] == 1 ||
            error("only single column supported")
        bound.zen_0[1] = nc_szen["szen"][1, 1]
    end

    #   The solar azimuthal angles.
    #   TODO - not needed for two stream

    #   The file of solar irradiances.
    NCDatasets.NCDataset(joinpath(example_dir, base_name*".stoa")) do nc_stoa
        nc_stoa.dim["lat"]*nc_stoa.dim["lon"] == 1 ||
            error("only single column supported")
        bound.solar_irrad[1] = nc_stoa["stoa"][1, 1]
    end

elseif control.isolir == SOCRATES.rad_pcf.ip_infra_red
    Bool(spectrum.Basic.l_present[6]) ||
        error("The spectral file contains no data for the Planckian function." )

    if Bool(spectrum.Basic.l_present[2])
        control.l_solar_tail_flux = true
    end
else
    error("unknown isolir $(control.isolir)")
end

###########################################
# Range of bands
###########################################

control.last_band = spectrum.Basic.n_band
control.last_band  <= spectrum.Basic.n_band ||
    error("last_band out of range")

control.first_band = 1
control.first_band  >=1 ||
    error("first_band out of range")

n_band_active = control.last_band - control.first_band + 1

# Map spectral bands into output channels
if ( (n_channel*floor(n_band_active/n_channel) != n_band_active)  &&
    (spectrum.Var.n_sub_band >= n_channel) )
    # Number of bands not a multiple of channels so use sub-bands
    control.l_map_sub_bands = true
end

SOCRATES.allocate_control(control, spectrum)

if n_channel == 1
    control.map_channel[1:spectrum.Basic.n_band] .= 1
elseif n_channel == spectrum.Basic.n_band
    control.map_channel[1:spectrum.Basic.n_band] .= 1:n_channel
else
    error("n_channel $n_channel != 1 and != $n_band_active not supported ")
end

# Calculate the weighting for the bands.
control.weight_band .= 1.0

##############################################
# Scaling of optical depth for direct flux
#############################################

# 'Entre treatment of optical depth for direct solar flux (0/1/2)'
# '0: no scaling; 1: delta-scaling; 2: circumsolar scaling'
control.i_direct_tau = 1

############################################
# Check Options
############################################

if control.l_rayleigh
    Bool(spectrum.Basic.l_present[3]) ||
        error("The spectral file contains no rayleigh scattering data.")
end

# -a
control.l_aerosol = false
if control.l_aerosol
    Bool(spectrum.Basic.l_present[11]) ||
        error("The spectral file contains no aerosol data.")
end

if control.l_gas
    Bool(spectrum.Basic.l_present[5]) ||
        error("The spectral file contains no gaseous absorption data.")
end

if control.l_continuum
    Bool(spectrum.Basic.l_present[9]) ||
        error("The spectral file contains no continuum absorption data.")
end

################################
# Gaseous absorption
#################################

if control.l_gas
    control.i_gas_overlap = SOCRATES.rad_pcf.ip_overlap_random # = 2
    for j in control.first_band:control.last_band
        control.i_gas_overlap_band[j] = control.i_gas_overlap
    end

    mr_gases = Dict()
    for i_gas in 1:spectrum.Gas.n_absorb
        # Read gas mixing ratios
        ti = spectrum.Gas.type_absorb[i_gas]
        sfx = SOCRATES.input_head_pcf.gas_suffix[ti]
        fn = joinpath(example_dir, base_name*"."*sfx)
        println("Reading mixing ratio for gas $ti $(SOCRATES.gas_list_pcf.name_absorb[ti])    $fn")

        if isfile(fn)
            NCDatasets.NCDataset(fn, "r") do nc
                mr_gases[sfx] = nc[sfx][1, 1, :]
                atm.gas_mix_ratio[1, :, i_gas] .= mr_gases[sfx]
            end
        else
            println("  no file found - setting mixing ratio to 0.0")
            atm.gas_mix_ratio[:, :, i_gas] .= 0.0
        end
    end
end

################################
# Aerosol processes
#################################

if control.l_aerosol 
    error("aerosols not implemented")
else
    dimen.nd_profile_aerosol_prsc   = 1
    dimen.nd_opt_level_aerosol_prsc = 1
    dimen.nd_phf_term_aerosol_prsc  = 1
    aer.mr_source .= SOCRATES.rad_pcf.ip_aersrc_classic_roff
end

SOCRATES.allocate_aer_prsc(aer, dimen, spectrum)
for i = 1:spectrum.Dim.nd_aerosol_species
    aer.mr_type_index[i] = i
end

#######################################
# Clouds
# see src/aux/input_cloud_cdf.f
#######################################

control.i_cloud = SOCRATES.rad_pcf.ip_cloud_off # 5 (clear sky)
control.l_cloud = false

dimen.nd_profile_cloud_prsc   = 1
dimen.nd_opt_level_cloud_prsc = 1
dimen.nd_phf_term_cloud_prsc  = 1

SOCRATES.allocate_cld_prsc(cld, dimen, spectrum)

#####################################
# Angular integration
# see src/aux/angular_control_cdf.f
#####################################

control.i_angular_integration = SOCRATES.rad_pcf.ip_two_stream

if control.i_angular_integration == SOCRATES.rad_pcf.ip_two_stream
    # see src/aux/angular_control_cdf.f

    if control.isolir == SOCRATES.rad_pcf.ip_infra_red
        control.i_2stream = 12 # -t 12: 
        # the two-stream approximation to use. 
        # 12 is the recommended approximation for the LW as used within the Met Office. 
        # (Practical improved flux method with diffusivity factor = 1.66.)
    else
        control.i_2stream = 16 # Cl_run_cdf 
    end

    control.l_rescale = false # Cl_run_cdf default  (+R for true)
    if control.l_rescale
        control.l_henyey_greenstein_pf = true
    end

    control.i_solver = 13 # -v 13: 
    # the solver used for the two-stream calculations. 
    # 13 is recommended for clear-sky, 
    # 16 is recommended for cloudy-sky,
    # 17 is recommended for cloud with separate stratiform and convective regions.

    #      Arrays of fluxes must be of the full size.
    dimen.nd_2sg_profile = dimen.nd_profile
    dimen.nd_flux_profile = dimen.nd_profile
    dimen.nd_radiance_profile = 1
    dimen.nd_j_profile = 1
    dimen.nd_viewing_level = 1
    dimen.nd_sph_coeff = 1

    #   Convert the zenith angles to secants.
    if control.isolir == SOCRATES.rad_pcf.ip_solar
        for l in 1:atm.n_profile
            bound.zen_0[l] = 1.0/cosd(bound.zen_0[l])
        end
    end

    # Reset dimen.nd_max_order to reduce memory requirements
    dimen.nd_max_order = 1

else
    error("unsupported i_angular_integration=$(control.i_angular_integration)")
end

#####################################
# surface properties
# see src/aux/assign_surface_char_cdf.f
# IP_surface_char  = 51, file suffix 'surf'
#####################################

# sets bound.rho_alb  from file
# rho_alb(l, IP_surf_alb_dir, i_band)
# rho_alb(l, IP_surf_alb_diff, i_band)

bound.rho_alb[:, SOCRATES.rad_pcf.ip_surf_alb_diff, :] .= surface_albedo

if control.i_angular_integration == SOCRATES.rad_pcf.ip_two_stream

#       Explicit basis functions are not required, but if there is
#       no separate direct albedo the diffuse value must be copied
#       into the direct field.

    if control.isolir == SOCRATES.rad_pcf.ip_solar
        bound.rho_alb[:, SOCRATES.rad_pcf.ip_surf_alb_dir, :] .=
            bound.rho_alb[:, SOCRATES.rad_pcf.ip_surf_alb_diff, :]
    end
end

###################################################
# Treatment of scattering
###################################################

control.i_scatter_method = SOCRATES.rad_pcf.ip_scatter_full
for i in control.first_band:control.last_band
    control.i_scatter_method_band[i] = control.i_scatter_method
end

###################################################
# Temperatures
###################################################

# vertical (assuming single column)
atm.p[1, :] .= nc_t["plev"][:]
atm.t[1, :] .= nc_t["t"][1, 1, :]

atm.p_level[1, 0:end] .= nc_tl["plev"][:]
atm.t_level[1, 0:end] .= nc_tl["tl"][1, 1, :]

####################################################
# Pressures
###################################################

# assume defined along with temperatures

####################################################
# Surface temperatures
# IP_temperature_ground = 9
# suffix .tstar
####################################################

if control.isolir == SOCRATES.rad_pcf.ip_infra_red
    bound.t_ground[1] = nc_tstar["tstar"][1, 1, 1]
end

#####################################################
# Variation of the temperature within layers
####################################################

if control.isolir == SOCRATES.rad_pcf.ip_infra_red
    # 'Is the ir-source function to be ' //      &
    # 'taken as linear or quadratic in tau? (l/q)'
    control.l_ir_source_quad = false
    # control.l_ir_source_quad = true  # Cl_run_cdf -q
end

######################################################
# Determination of the mass and density in each layer
######################################################

# Physical constants from 
# src/modules_core/rad_ccf.F90

grav_acc           = 9.80665

# Gas constant for dry air
r_gas_dry          = 287.026

# Molecular weight of dry air/ molecular weight of water
ratio_molar_weight = 28.966 / 18.0153

r = r_gas_dry
c_virtual = ratio_molar_weight - 1.0

for i in range(atm.n_layer, 1, step=-1)
    for l in 1:atm.n_profile
        atm.mass[l, i] = (atm.p_level[l, i] - atm.p_level[l, i-1])/grav_acc
    end
end

# dry air density (TODO H2O)
for i in 1:atm.n_layer
    for l in 1:atm.n_profile
        atm.density[l, i] = atm.p[l, i]/(r*atm.t[l, i])
    end
end

##########################################
# Determination of heights, and lit layers
###########################################

if control.l_spherical_solar
    error("l_spherical_solar not supported")
end

#########################################
# Subgrid cloud generator
#########################################

if control.i_cloud == SOCRATES.rad_pcf.ip_cloud_mcica
    error("MCICA subgrid clouds not supported")
end

##############################################
# Calculate the clear-sky mean relative humidity for moist aerosols.
###################################################

!control.l_aerosol ||
    error("(moist) aerosols not supported")

##########################################################
# Calculation of radiances or irradiances.
############################################################

SOCRATES.radiance_calc(control, dimen, spectrum, atm, cld, aer, bound, radout)

# display output
if do_plot
    p = plot(atm.t_level[1, :], atm.p_level[1, :], label="t_level", yflip=true, xlabel="Temperature (K)", ylabel="pressure (Pa)")
    display(p)
end

if do_plot && control.l_gas
    p = plot(yflip=true, xscale=:log10, xlabel="Mass mixing ratio (g/g)", ylabel="pressure (Pa)")
    for (sfx, mr) in mr_gases
        plot!(mr, atm.p[1, :], label=sfx)
    end
    display(p)
end

println()
println("radout.flux_down (total):")
display(radout.flux_down)
if do_plot
    p = plot(radout.flux_down[1, :, 1], atm.p_level[1, :], label="flux_down", yflip=true, xlabel="W m-2", ylabel="pressure (Pa)")
end

println()
println("radout.flux_up:")
display(radout.flux_up)
if do_plot
    plot!(radout.flux_up[1, :, 1], atm.p_level[1, :], label="flux_up")
end

if control.isolir == SOCRATES.rad_pcf.ip_solar
    println()
    println("radout.flux_direct:")
    display(radout.flux_direct)
    if do_plot
        plot!(radout.flux_direct[1, :, 1], atm.p_level[1, :], label="flux_direct")
    end

    flux_diffuse_down = radout.flux_down[:, :, :] .- radout.flux_direct[:, :, :]
    println()
    println("flux_diffuse_down: (flux_down - flux_direct)")
    display(flux_diffuse_down)
    if do_plot
        plot!(flux_diffuse_down[1, :, 1], atm.p_level[1, :], label="flux_diffuse_down")
    end
end
if control.l_flux_div
    println()
    println("radout.flux_div:")
    display(radout.flux_div)
end

if do_plot
    display(p)
end

##############################################
# Deallocate arrays
#############################################

# close netcdf files
close.(nc_open)

if do_deallocate
    SOCRATES.deallocate_atm(atm)
    SOCRATES.deallocate_cld(cld)
    SOCRATES.deallocate_cld_prsc(cld)
    SOCRATES.deallocate_aer(aer)
    SOCRATES.deallocate_aer_prsc(aer)
    SOCRATES.deallocate_bound(bound)
    SOCRATES.deallocate_out(radout)
end
