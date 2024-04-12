#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# This file runs a series of tests, to make sure that the model is behaving properly.
# -------------


# Get AGNI root directory
ROOT_DIR = dirname(abspath(@__FILE__))
ENV["GKSwstype"] = "100"

# Include libraries
using Revise
using LoggingExtras

@info "Begin tests"


# Include local jl files
include("socrates/julia/src/SOCRATES.jl")
push!(LOAD_PATH, joinpath(ROOT_DIR,"src"))
import atmosphere
import setpt
import plotting 
import phys

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
c_expt::Array{Float64, 1} = [35.22, 41.27 , 51.20 , 55.74 , 59.40 ]     # Expected values of cp [J mol-1 K-1]
t_test::Array{Float64, 1} = [500.0, 1000.0, 2000.0, 3000.0, 5000.0]     # Tested values of temperature 
cp_pass = true 
for i in 1:5
    c_this = phys.shomate_cp("H2O", t_test[i]) * 18.0153 * 1.0e-3 # get value and convert units 
    if abs(c_expt[i]-c_this)/c_expt[i] > 0.01  # error must be <1%
        @warn "At tmp=$(t_test[i]) K \nModelled value = $c_this J K-1 mol-1 \nExpected value = $(c_expt[i]) J K-1 mol-1"
        global cp_pass = false 
    end 
end 
if cp_pass
    @info "Pass"
else
    @warn "Fail"
    passing = false
end




# -------------
# Test mixing ratios
# -------------
@info " "
@info "Testing composition"

tmp_surf           = 200.0    # Surface temperature [kelvin]
toa_heating     = 1000.00  # Instellation flux [W m-2]
p_surf          = 50.0     # bar
theta           = 65.0
mf_dict         = Dict([
                        ("H2O" , 0.5),
                        ("CO2" , 0.1),
                        ("O2"  , 0.1),
                        ("He" , 0.1),
                        ("H2"  , 0.2)
                        ])
spfile_name   = "res/spectral_files/Mallard/Mallard.sf"

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tmp_surf,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict
                 )
atmosphere.allocate!(atmos,"")

dct_e::Dict{String, Float64} = mf_dict
dct_o::Dict{String, Float64} = Dict()
sp_pass = true
for k in keys(dct_e)
    dct_o[k] = atmosphere.get_x(atmos, k, 25)
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





# -------------
# Test instellation
# -------------
@info " "
@info "Testing instellation"

tmp_surf           = 200.0    # Surface temperature [kelvin]
toa_heating     = 1000.00     # Instellation flux [W m-2]
p_surf          = 50.0    # bar
theta           = 65.0
mf_dict         = Dict([
                        ("H2O" , 0.8),
                        ("CO2" , 0.2),
                        ])
spfile_name   = "res/spectral_files/Mallard/Mallard.sf"

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tmp_surf,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict,
                         flag_gcontinuum=false,
                         flag_rayleigh=false,
                         overlap_method=2
                 )
atmosphere.allocate!(atmos,"res/stellar_spectra/sun.txt")
atmosphere.radtrans!(atmos, false)

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
atmosphere.deallocate!(atmos)




# -------------
# Test greenhouse effect
# -------------
@info " "
@info "Testing greenhouse effect"

tmp_surf           = 1300.0    # Surface temperature [kelvin]
toa_heating     = 1000.00    # Instellation flux [W m-2]
p_surf          = 300.0    # bar
theta           = 45.0
mf_dict         = Dict([
                        ("H2O" , 1.0),
                        ])
spfile_name   = "res/spectral_files/Oak/Oak.sf"

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tmp_surf,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict,
                         flag_gcontinuum=true,
                         flag_rayleigh=false,
                         overlap_method=4
                 )
atmosphere.allocate!(atmos,"res/stellar_spectra/sun.txt")
setpt.prevent_surfsupersat!(atmos)
setpt.dry_adiabat!(atmos)
setpt.condensing!(atmos, "H2O")
atmosphere.radtrans!(atmos, true)

val_e = [270.0, 290.0]
val_o = atmos.flux_u_lw[1]
@info "Expected range = $(val_e) W m-2"
@info "Modelled value = $(val_o) W m-2"
if ( val_o > val_e[1]) && (val_o < val_e[2]) 
    @info "Pass"
else
    @warn "Fail"
    passing = false
end
atmosphere.deallocate!(atmos)



# -------------
# Test Rayleigh scattering
# -------------
@info " "
@info "Testing Rayleigh scattering"

tmp_surf           = 400.0    # Surface temperature [kelvin]
toa_heating     = 1000.00    # Instellation flux [W m-2]
p_surf          = 10.0    # bar
theta           = 75.0
mf_dict         = Dict([
                        ("H2O" , 0.6),
                        ("CO2" , 0.4),
                        ])
spfile_name   = "res/spectral_files/Mallard/Mallard.sf"

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tmp_surf,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict,
                         flag_gcontinuum=true,
                         flag_rayleigh=true,
                         overlap_method=2
                 )
atmosphere.allocate!(atmos,"res/stellar_spectra/sun.txt")
atmosphere.radtrans!(atmos, true)
atmosphere.radtrans!(atmos, false)

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




# -------------
# Test heating rate calculation
# -------------
@info " "
@info "Testing heating rates"

tmp_surf           = 2500.0    # Surface temperature [kelvin]
toa_heating     = 1000.0    # Instellation flux [W m-2]
p_surf          = 5.0     # bar
theta           = 45.0
mf_dict         = Dict([
                        ("H2O" , 1.0)
                        ])
spfile_name   = "res/spectral_files/Oak/Oak.sf"

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tmp_surf,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict,
                         flag_gcontinuum=true,
                         flag_rayleigh=false,
                         overlap_method=2,
                         thermo_functions=false
                 )
atmosphere.allocate!(atmos,"res/stellar_spectra/trappist-1.txt")
setpt.isothermal!(atmos, 300.0)
atmos.flux_tot[:] .= 0.0
atmosphere.radtrans!(atmos, true)
atmosphere.radtrans!(atmos, false)
atmos.flux_tot += atmos.flux_n
atmosphere.calc_hrates!(atmos)
val_e = [20.0, 70.0]  # tests have found ~54 K/day for this setup
val_o = atmos.heating_rate[atmos.nlev_c-2]
@info "Expected range = $(val_e) K/day"
@info "Modelled value = $(val_o) K/day"
if ( val_o > val_e[1]) && (val_o < val_e[2]) 
    @info "Pass"
else
    @warn "Fail"
    passing = false
end
atmosphere.deallocate!(atmos)


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

