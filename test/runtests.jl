#!/usr/bin/env -S julia --color=yes --startup-file=no

# Get AGNI root directory
ROOT_DIR = abspath(joinpath(dirname(abspath(PROGRAM_FILE)),"../"))

# Activate environment
ENV["GKSwstype"] = "100"
import Pkg
Pkg.activate(ROOT_DIR)

# Include libraries
using LoggingExtras
using AGNI

@info "Begin tests"


# Prepare
output_dir    = "out/"
p_top           = 1e-8
nlev_centre     = 100
radius          = 1.0e7    # metres
gravity         = 10.0      # m s-2
total  = 0
failed = 0

rm(output_dir,force=true,recursive=true)
if !isdir(output_dir) && !isfile(output_dir)
    mkdir(output_dir)
end


rtol   = 1e-3


# -------------
# Test heat capacity lookup tables
# -------------
@info " "
@info "Testing heat capacity functions"
ideal_H2O::phys.Gas_t = phys.load_gas("res/thermodynamics/", "H2O", true, false)
t_test  = [10.0,  500.0, 1000.0, 2000.0, 3000.0]     # Tested values of temperature
v_expt  = [4.975, 35.22, 41.27 , 51.20 , 55.74 ]     # Expected values of cp [J mol-1 K-1]
v_obs   = zero(t_test)
test_pass = true
for i in 1:5
    v_obs[i] = phys.get_Cp(ideal_H2O, t_test[i]) * ideal_H2O.mmw # get value and convert units
    if abs(v_expt[i]- v_obs[i])/v_expt[i] > rtol
        global test_pass = false
    end
end
@info "Expected values = $(v_expt) J mol-1 K-1"
@info "Modelled values = $(v_obs) J mol-1 K-1"
if test_pass
    @info "Pass"
else
    @warn "Fail"
    failed += 1
end
total += 1
@info "--------------------------"


# -------------
# Test ideal gas equation of state
# -------------
@info " "
@info "Testing H2O ideal gas equation of state"
t_test = [200.0,  300.0, 500.0,   1273.0,  3200.0] # Tested values of temperature [K]
p_test = [1e0,    1e3,   1e5,     1e7,     1e8]    # Tested values of pressure [Pa]
v_expt = [1.0833532e-5, 7.2223549e-3, 4.33341295e-1, 1.7020475e1, 6.7709577e1]  # Expected rho [kg m-3]
v_obs  = zero(p_test)
test_pass = true
for i in 1:5
    v_obs[i] = phys.calc_rho_gas(t_test[i], p_test[i], ideal_H2O)
    if abs(v_expt[i] - v_obs[i])/v_expt[i] > rtol
        global test_pass = false
    end
end
@info "Expected values = $(v_expt) kg m-3"
@info "Modelled values = $(v_obs) kg m-3"
if test_pass
    @info "Pass"
else
    @warn "Fail"
    failed += 1
end
total += 1
@info "--------------------------"

# -------------
# Test AQUA equation of state
# -------------
@info " "
@info "Testing H2O AQUA equation of state"
aqua_H2O::phys.Gas_t = phys.load_gas("res/thermodynamics/", "H2O", true, true)
v_expt = [9.26121571e2, 7.2269991e-3, 4.35193299e-1, 1.7022212e1, 6.68198662e1]
v_obs  = zero(p_test)
test_pass = true
for i in 1:5
    v_obs[i] = phys.calc_rho_gas(t_test[i], p_test[i], aqua_H2O)
    if abs(v_expt[i] - v_obs[i])/v_expt[i] > rtol
        global test_pass = false
    end
end
@info "Expected values = $(v_expt) kg m-3"
@info "Modelled values = $(v_obs) kg m-3"
if test_pass
    @info "Pass"
else
    @warn "Fail"
    failed += 1
end
total += 1
@info "--------------------------"


# -------------
# Test VdW equation of state
# -------------
@info " "
@info "Testing CO2 VdW equation of state"
vdw_CO2::phys.Gas_t = phys.load_gas("res/thermodynamics/", "CO2", true, true)
v_expt = [2.6465333e-5, 1.7664486595e-2, 1.0612271156, 4.12385376e1, 1.47631823e2]
v_obs  = zero(p_test)
test_pass = true
for i in 1:5
    v_obs[i] = phys.calc_rho_gas(t_test[i], p_test[i], vdw_CO2)
    if abs(v_expt[i] - v_obs[i])/v_expt[i] > rtol
        global test_pass = false
    end
end
@info "Expected values = $(v_expt) kg m-3"
@info "Modelled values = $(v_obs) kg m-3"
if test_pass
    @info "Pass"
else
    @warn "Fail"
    failed += 1
end
total += 1
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
test_pass = true
for k in keys(dct_e)
    dct_o[k] = atmos.gas_vmr[k][20]
    global test_pass = test_pass && ( abs(dct_o[k]-dct_e[k]) < 1.0e-6 )
end
@info "Expected values = $(dct_e)"
@info "Modelled values = $(dct_o)"
if test_pass
    @info "Pass"
else
    @warn "Fail"
    failed += 1
end
total += 1
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
                         overlap_method="ro",
                         real_gas=false
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
    failed += 1
end
total += 1
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
    failed += 1
end
total += 1
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
                         300, p_surf, p_top,
                         mf_dict,"",
                         flag_gcontinuum=true,
                         flag_rayleigh=false,
                         overlap_method="ee",
                         condensates=["H2O"],
                         surface_material="greybody",
                         albedo_s=0.5,
                         real_gas=false
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
    failed += 1
end
total += 1
@info "--------------------------"


# -------------
# Test hydrostatic integrator
# -------------
@info " "
@info "Testing hydrostatic integration"

val_e = 424238.22072713776   # known from previous tests
val_o = atmos.z[1] # top level
@info "Expected value = $(val_e) m"
@info "Modelled value = $(val_o) m"
if abs(val_o - val_e)/val_e < rtol
    @info "Pass"
else
    @warn "Fail"
    failed += 1
end
total += 1
@info "--------------------------"


# -------------
# Test write NetCDF
# -------------
@info " "
@info "Testing write NetCDF"
out_path::String = "/tmp/agni_atm.nc"
rm(out_path, force=true)
dump.write_ncdf(atmos, out_path)
@info "Expecting file at $out_path"
if isfile(out_path)
    @info "Found file at $out_path"
    @info "Pass"
    rm(out_path, force=true)
else
    @warn "Fail"
    failed += 1
end
total += 1
@info "--------------------------"

# -------------
# Test plot T(p)
# -------------
@info " "
@info "Testing plot temperatures"
out_path = "/tmp/agni_plot_tmp.png"
rm(out_path, force=true)
plotting.plot_pt(atmos, out_path)
@info "Expecting file at $out_path"
if isfile(out_path)
    @info "Found file at $out_path"
    @info "Pass"
    rm(out_path, force=true)
else
    @warn "Fail"
    failed += 1
end
total += 1
@info "--------------------------"


# -------------
# Test plot height
# -------------
@info " "
@info "Testing plot height"
out_path = "/tmp/agni_plot_hei.png"
plotting.plot_height(atmos, out_path)
@info "Expecting file at $out_path"
if isfile(out_path)
    @info "Found file at $out_path"
    @info "Pass"
    rm(out_path, force=true)
else
    @warn "Fail"
    failed += 1
end
total += 1
@info "--------------------------"


# -------------
# Test plot albedo
# -------------
@info " "
@info "Testing plot albedo"
out_path = "/tmp/agni_plot_alb.png"
plotting.plot_albedo(atmos, out_path)
@info "Expecting file at $out_path"
if isfile(out_path)
    @info "Found file at $out_path"
    @info "Pass"
    rm(out_path, force=true)
else
    @warn "Fail"
    failed += 1
end
total += 1
@info "--------------------------"


# -------------
# Test surface albedo
# -------------
@info " "
@info "Testing surface albedo "
val_e = 30.24053638241024  # known from previous tests
val_o = atmos.flux_u_sw[end] # bottom level
@info "Expected value = $(val_e) W m-2"
@info "Modelled value = $(val_o) W m-2"
if abs(val_o - val_e)/val_e < rtol
    @info "Pass"
else
    @warn "Fail"
    failed += 1
end
total += 1
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
                         overlap_method="ro",
                         real_gas=false
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
    failed  += 1
end
total += 1
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
                         thermo_functions=false,
                         real_gas=false
                 )
atmosphere.allocate!(atmos,joinpath(ROOT_DIR,"res/stellar_spectra/sun.txt"))
setpt.isothermal!(atmos, 300.0)
atmosphere.calc_layer_props!(atmos)
atmos.flux_tot[:] .= 0.0
energy.radtrans!(atmos, true)
energy.radtrans!(atmos, false)
atmos.flux_tot += atmos.flux_n
energy.calc_hrates!(atmos)

val_e = 6.366596000871719  # from previous tests
val_o = atmos.heating_rate[atmos.nlev_c-10]
@info "Expected value = $(val_e) K/day"
@info "Modelled value = $(val_o) K/day"
if abs(val_o - val_e)/val_e < rtol
    @info "Pass"
else
    @warn "Fail"
    failed += 1
end
total += 1
@info "--------------------------"


# -------------
# Test flux calculation
# -------------
@info " "
@info "Testing fluxes"
energy.calc_fluxes!(atmos, true, true, true, true)
val_e = 8.60308315276e3  # from previous tests
val_o = atmos.flux_tot[atmos.nlev_c-10]
@info "Expected value = $(val_e) W m-2"
@info "Modelled value = $(val_o) W m-2"
if abs(val_o - val_e)/val_e < rtol
    @info "Pass"
else
    @warn "Fail"
    failed += 1
end
total += 1
atmosphere.deallocate!(atmos)
@info "--------------------------"


# -------------
# Inform at end
# -------------
@info " "
if failed == 0
    @info "All $total tests have passed"
    exit(0)
else
    @warn "Some tests failed ($failed/$total)"
    exit(1)
end

