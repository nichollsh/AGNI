#!/usr/bin/env -S julia --color=yes --startup-file=no

# Get AGNI root directory
ROOT_DIR = abspath(joinpath(dirname(abspath(PROGRAM_FILE)),"../"))

# Activate environment
ENV["GKSwstype"] = "100"
import Pkg
Pkg.activate(ROOT_DIR)

# Include libraries
using LoggingExtras

@info "Begin tests"

# Include local jl files
include("../src/phys.jl")
include("../src/spectrum.jl")
include("../src/atmosphere.jl")
include("../src/setpt.jl")
include("../src/energy.jl")
include("../src/plotting.jl")
import .phys
import .atmosphere
import .setpt
import .energy
import .plotting


# Prepare
output_dir    = "out/"
p_top           = 1e-8
nlev_centre     = 100
radius          = 1.0e7    # metres
gravity         = 10.0      # m s-2

rm(output_dir,force=true,recursive=true)
if !isdir(output_dir) && !isfile(output_dir)
    mkdir(output_dir)
end

passing = true


# -------------
# Test shomate
# -------------
@info " "
@info "Testing heat capacity functions"
data_H2O::phys.Gas_t = phys.load_gas("res/thermodynamics/", "H2O", true)
c_expt::Array{Float64, 1} = [4.975, 35.22, 41.27 , 51.20 , 55.74 ]     # Expected values of cp [J mol-1 K-1]
t_test::Array{Float64, 1} = [10.0,  500.0, 1000.0, 2000.0, 3000.0]     # Tested values of temperature
c_obs::Array{Float64,1} = zeros(Float64, 5)
cp_pass = true
for i in 1:5
    c_obs[i] = phys.get_Cp(data_H2O, t_test[i]) * data_H2O.mmw # get value and convert units
    if abs(c_expt[i]- c_obs[i])/c_expt[i] > 0.01  # error must be <1%
        global cp_pass = false
    end
end
@info "Expected values = $(c_expt) J mol-1 K-1"
@info "Modelled values = $(c_obs) J mol-1 K-1"
if cp_pass
    @info "Pass"
else
    @warn "Fail"
    passing = false
end
@info "--------------------------"




# -------------
# Test mixing ratios
# -------------
@info " "
@info "Testing composition"

tmp_surf        = 200.0     # Surface temperature [kelvin]
toa_heating     = 10000.00  # Instellation flux [W m-2]
p_surf          = 1.0       # bar
theta           = 65.0
mf_dict         = Dict([
                        ("H2O" , 0.5),
                        ("CO2" , 0.2),
                        ("N2"  , 0.1),
                        ("H2"  , 0.2)
                        ])
spfile_name   = joinpath(ROOT_DIR,"res/spectral_files/Dayspring/48/Dayspring.sf")

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir,
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tmp_surf,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict, ""
                 )
atmosphere.allocate!(atmos,"")

dct_e::Dict{String, Float64} = mf_dict
dct_o::Dict{String, Float64} = Dict()
sp_pass = true
for k in keys(dct_e)
    dct_o[k] = atmos.gas_vmr[k][20]
    global sp_pass = sp_pass && ( abs(dct_o[k]-dct_e[k]) < 1.0e-6 )
end
@info "Expected values = $(dct_e)"
@info "Modelled values = $(dct_o)"
if sp_pass
    @info "Pass"
else
    @warn "Fail"
    passing = false
end
atmosphere.deallocate!(atmos)
@info "--------------------------"




# -------------
# Test instellation
# -------------
@info " "
@info "Testing instellation"

tmp_surf        = 200.0    # Surface temperature [kelvin]
toa_heating     = 1000.00     # Instellation flux [W m-2]
p_surf          = 50.0    # bar
theta           = 65.0
mf_dict         = Dict([
                        ("H2O" , 0.8),
                        ("CO2" , 0.2),
                        ])
spfile_name   = joinpath(ROOT_DIR,"res/spectral_files/Dayspring/48/Dayspring.sf")

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir,
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tmp_surf,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict,"",
                         flag_gcontinuum=false,
                         flag_rayleigh=false,
                         overlap_method="ro"
                 )
atmosphere.allocate!(atmos,joinpath(ROOT_DIR,"res/stellar_spectra/sun.txt"))
energy.radtrans!(atmos, false)

val_e = toa_heating * cosd(theta)
val_o = atmos.flux_d_sw[1]
@info "Expected value = $(val_e) W m-2"
@info "Modelled value = $(val_o) W m-2"
if abs(val_e-val_o) < 2
    @info "Pass"
else
    @warn "Fail"
    passing = false
end
@info "--------------------------"


# -------------
# Test no-scattering case
# -------------
@info " "
@info "Testing zero scattering"

val_o = atmos.flux_u_sw[2]
val_e = 0.0
@info "Expected value = $(val_e) W m-2"
@info "Modelled value = $(val_o) W m-2"
if abs(val_e-val_o) < 1.0e-10
    @info "Pass"
else
    @warn "Fail"
    passing = false
end
atmosphere.deallocate!(atmos)
@info "--------------------------"



# -------------
# Test greenhouse effect
# -------------
@info " "
@info "Testing greenhouse effect"

tmp_surf        = 1300.0    # Surface temperature [kelvin]
toa_heating     = 1000.00    # Instellation flux [W m-2]
p_surf          = 300.0    # bar
theta           = 45.0
mf_dict         = Dict([
                        ("H2O" , 1.0),
                        ("N2"  , 1.0e-9)
                        ])
spfile_name   = joinpath(ROOT_DIR,"res/spectral_files/Oak/318/Oak.sf")

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir,
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tmp_surf,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict,"",
                         flag_gcontinuum=true,
                         flag_rayleigh=false,
                         overlap_method="ee",
                         condensates=["H2O"],
                         surface_material="greybody",
                         albedo_s=0.5
                 )
atmosphere.allocate!(atmos,joinpath(ROOT_DIR,"res/stellar_spectra/sun.txt"))
setpt.prevent_surfsupersat!(atmos)
setpt.dry_adiabat!(atmos)
setpt.saturation!(atmos, "H2O")
atmosphere.calc_layer_props!(atmos)
energy.radtrans!(atmos, true)
energy.radtrans!(atmos, false)

val_e = [270.0, 280.0]
val_o = atmos.flux_u_lw[1]
@info "Expected range = $(val_e) W m-2"
@info "Modelled value = $(val_o) W m-2"
if ( val_o > val_e[1]) && (val_o < val_e[2])
    @info "Pass"
else
    @warn "Fail"
    passing = false
end
@info "--------------------------"


# -------------
# Test hydrostatic integrator
# -------------
@info " "
@info "Testing hydrostatic integration"

val_e = 431792.0902977038  # known from previous tests
val_o = atmos.z[1] # top level
@info "Expected value = $(val_e) m"
@info "Modelled value = $(val_o) m"
if abs(val_o - val_e) < 0.1
    @info "Pass"
else
    @warn "Fail"
    passing = false
end
@info "--------------------------"


# -------------
# Test plot T(p)
# -------------
@info " "
@info "Testing plot temperatures"
plt_path::String = "/tmp/agni_plot_tmp.png"
rm(plt_path, force=true)
plotting.plot_pt(atmos, plt_path)
@info "Expecting file at $plt_path"
if isfile(plt_path)
    @info "Found file at $plt_path"
    @info "Pass"
    rm(plt_path, force=true)
else
    @warn "Fail"
    passing = false
end
@info "--------------------------"


# -------------
# Test plot height
# -------------
@info " "
@info "Testing plot height"
plt_path = "/tmp/agni_plot_hei.png"
plotting.plot_height(atmos, plt_path)
@info "Expecting file at $plt_path"
if isfile(plt_path)
    @info "Found file at $plt_path"
    @info "Pass"
    rm(plt_path, force=true)
else
    @warn "Fail"
    passing = false
end
@info "--------------------------"


# -------------
# Test plot albedo
# -------------
@info " "
@info "Testing plot albedo"
plt_path = "/tmp/agni_plot_alb.png"
plotting.plot_albedo(atmos, plt_path)
@info "Expecting file at $plt_path"
if isfile(plt_path)
    @info "Found file at $plt_path"
    @info "Pass"
    rm(plt_path, force=true)
else
    @warn "Fail"
    passing = false
end
@info "--------------------------"


# -------------
# Test surface albedo
# -------------
@info " "
@info "Testing surface albedo "
val_e = 30.31364083727528  # known from previous tests
val_o = atmos.flux_u_sw[end] # bottom level
@info "Expected value = $(val_e) W m-2"
@info "Modelled value = $(val_o) W m-2"
if abs(val_o - val_e) < 1e-6
    @info "Pass"
else
    @warn "Fail"
    passing = false
end
atmosphere.deallocate!(atmos)
@info "--------------------------"


# -------------
# Test Rayleigh scattering
# -------------
@info " "
@info "Testing Rayleigh scattering"

tmp_surf        = 400.0    # Surface temperature [kelvin]
toa_heating     = 1000.00    # Instellation flux [W m-2]
p_surf          = 10.0    # bar
theta           = 75.0
mf_dict         = Dict([
                        ("H2O" , 0.6),
                        ("CO2" , 0.4),
                        ])
spfile_name   = joinpath(ROOT_DIR,"res/spectral_files/Dayspring/48/Dayspring.sf")

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir,
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tmp_surf,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict, "",
                         flag_gcontinuum=true,
                         flag_rayleigh=true,
                         overlap_method="ro"
                 )
atmosphere.allocate!(atmos,joinpath(ROOT_DIR,"res/stellar_spectra/sun.txt"))
energy.radtrans!(atmos, true)
energy.radtrans!(atmos, false)

val_e = [1.0e-4, 1e9]
val_o = atmos.flux_u_sw[20]
@info "Expected range = $(val_e) W m-2"
@info "Modelled value = $(val_o) W m-2"
if ( val_o > val_e[1]) && (val_o < val_e[2])
    @info "Pass"
else
    @warn "Fail"
    passing = false
end
atmosphere.deallocate!(atmos)
@info "--------------------------"


# -------------
# Test heating rate calculation
# -------------
@info " "
@info "Testing heating rates"

tmp_surf        = 2500.0    # Surface temperature [kelvin]
toa_heating     = 1000.0    # Instellation flux [W m-2]
p_surf          = 5.0     # bar
theta           = 45.0
mf_dict         = Dict([
                        ("H2O" , 1.0)
                        ])
spfile_name   = joinpath(ROOT_DIR,"res/spectral_files/Oak/318/Oak.sf")

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir,
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tmp_surf,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict, "",
                         flag_gcontinuum=true,
                         flag_rayleigh=false,
                         overlap_method="ro",
                         thermo_functions=false
                 )
atmosphere.allocate!(atmos,joinpath(ROOT_DIR,"res/stellar_spectra/sun.txt"))
setpt.isothermal!(atmos, 300.0)
atmosphere.calc_layer_props!(atmos)
atmos.flux_tot[:] .= 0.0
energy.radtrans!(atmos, true)
energy.radtrans!(atmos, false)
atmos.flux_tot += atmos.flux_n
energy.calc_hrates!(atmos)

val_e = 6.366831453838685  # from previous tests
val_o = atmos.heating_rate[atmos.nlev_c-10]
@info "Expected value = $(val_e) K/day"
@info "Modelled value = $(val_o) K/day"
if abs(val_o - val_e) < 1e-6
    @info "Pass"
else
    @warn "Fail"
    passing = false
end
atmosphere.deallocate!(atmos)
@info "--------------------------"


# -------------
# Inform at end
# -------------
@info " "
if passing
    @info "All tests passed"
    exit(0)
else
    @warn "Some tests failed"
    exit(1)
end

