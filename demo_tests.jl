#!/usr/bin/env -S julia --color=yes --startup-file=no

# -------------
# This file runs a series of tests, to make sure that the model is behaving properly.
# -------------

println("Begin tests")


# Get AGNI root directory
ROOT_DIR = dirname(abspath(@__FILE__))
ENV["GKSwstype"] = "100"

# Include libraries
using Revise

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
# Test mixing ratios
# -------------
println(" ")
println("Testing composition")

tstar           = 200.0    # Surface temperature [kelvin]
toa_heating     = 1000.00  # Instellation flux [W m-2]
p_surf          = 50.0     # bar
theta           = 65.0
mf_dict         = Dict([
                        ("H2O" , 0.5),
                        ("CO2" , 0.1),
                        ("He"  , 0.1),
                        ("SO2" , 0.1),
                        ("O3"  , 0.2)
                        ])
spfile_name   = "res/spectral_files/Reach/Reach"
star_file     = "res/stellar_spectra/sun.txt"

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tstar,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict
                 )
atmosphere.allocate!(atmos;stellar_spectrum=star_file)

val_e = mf_dict
val_o = Dict()
sp_pass = true
for k in keys(val_e)
    val_o[k] = atmosphere.get_x(atmos, k, 25)
    global sp_pass = sp_pass && ( abs(val_o[k]-val_e[k]) < 1.0e-6 )
end
println("Expected values = $(val_e)")
println("Modelled values = $(val_o)")
if sp_pass
    println("Pass")
else
    println("Fail")
    passing = false
end
atmosphere.deallocate!(atmos)





# -------------
# Test instellation
# -------------
println(" ")
println("Testing instellation")

tstar           = 200.0    # Surface temperature [kelvin]
toa_heating     = 1000.00     # Instellation flux [W m-2]
p_surf          = 50.0    # bar
theta           = 65.0
mf_dict         = Dict([
                        ("H2O" , 0.8),
                        ("CO2" , 0.2),
                        ])
spfile_name   = "res/spectral_files/Mallard/Mallard"
star_file     = "res/stellar_spectra/sun.txt"

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tstar,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict,
                         flag_gcontinuum=false,
                         flag_rayleigh=false,
                         overlap_method=2
                 )
atmosphere.allocate!(atmos;stellar_spectrum=star_file)
atmosphere.radtrans!(atmos, false)

val_e = toa_heating * cosd(theta)
val_o = atmos.flux_d_sw[1]
println("Expected value = $(val_e) W m-2")
println("Modelled value = $(val_o) W m-2")
if abs(val_e-val_o) < 2
    println("Pass")
else
    println("Fail")
    passing = false
end
atmosphere.deallocate!(atmos)




# -------------
# Test greenhouse effect
# -------------
println(" ")
println("Testing greenhouse effect")

tstar           = 1300.0    # Surface temperature [kelvin]
toa_heating     = 1000.00    # Instellation flux [W m-2]
p_surf          = 300.0    # bar
theta           = 45.0
mf_dict         = Dict([
                        ("H2O" , 1.0),
                        ])
spfile_name   = "res/spectral_files/Oak/Oak"
star_file     = "res/stellar_spectra/sun.txt"

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tstar,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict,
                         flag_gcontinuum=true,
                         flag_rayleigh=false,
                         overlap_method=4
                 )
atmosphere.allocate!(atmos;stellar_spectrum=star_file)
setpt.prevent_surfsupersat!(atmos)
setpt.dry_adiabat!(atmos)
setpt.condensing!(atmos, "H2O")
atmosphere.radtrans!(atmos, true)

val_e = [275.0, 290.0]
val_o = atmos.flux_u_lw[1]
println("Expected range = $(val_e) W m-2")
println("Modelled value = $(val_o) W m-2")
if ( val_o > val_e[1]) && (val_o < val_e[2]) 
    println("Pass")
else
    println("Fail")
    passing = false
end
atmosphere.deallocate!(atmos)



# -------------
# Test Rayleigh scattering
# -------------
println(" ")
println("Testing Rayleigh scattering")

tstar           = 400.0    # Surface temperature [kelvin]
toa_heating     = 1000.00    # Instellation flux [W m-2]
p_surf          = 10.0    # bar
theta           = 75.0
mf_dict         = Dict([
                        ("H2O" , 0.6),
                        ("CO2" , 0.4),
                        ])
spfile_name   = "res/spectral_files/Mallard/Mallard"
star_file     = "res/stellar_spectra/sun.txt"

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tstar,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict,
                         flag_gcontinuum=true,
                         flag_rayleigh=true,
                         overlap_method=2
                 )
atmosphere.allocate!(atmos;stellar_spectrum=star_file, spfile_noremove=true)
atmosphere.radtrans!(atmos, true)
atmosphere.radtrans!(atmos, false)

val_e = [1.0e-4, 1e9]
val_o = atmos.flux_u_sw[20]
println("Expected range = $(val_e) W m-2")
println("Modelled value = $(val_o) W m-2")
if ( val_o > val_e[1]) && (val_o < val_e[2]) 
    println("Pass")
else
    println("Fail")
    passing = false
end
atmosphere.deallocate!(atmos)




# -------------
# Test heating rate calculation
# -------------
println(" ")
println("Testing heating rates")

tstar           = 2500.0    # Surface temperature [kelvin]
toa_heating     = 1000.0    # Instellation flux [W m-2]
p_surf          = 5.0     # bar
theta           = 45.0
mf_dict         = Dict([
                        ("H2O" , 1.0)
                        ])
spfile_name   = "res/spectral_files/Oak/Oak"
star_file     = "res/stellar_spectra/trappist-1.txt"

# Setup atmosphere
atmos = atmosphere.Atmos_t()
atmosphere.setup!(atmos, ROOT_DIR, output_dir, 
                         spfile_name,
                         toa_heating, 1.0, 0.0, theta,
                         tstar,
                         gravity, radius,
                         nlev_centre, p_surf, p_top,
                         mf_dict=mf_dict,
                         flag_gcontinuum=true,
                         flag_rayleigh=false,
                         overlap_method=2,
                         thermo_functions=false
                 )
atmosphere.allocate!(atmos;stellar_spectrum=star_file, spfile_noremove=true)
setpt.isothermal!(atmos, 300.0)
atmos.flux_tot[:] .= 0.0
atmosphere.radtrans!(atmos, true)
atmosphere.radtrans!(atmos, false)
atmos.flux_tot += atmos.flux_n
atmosphere.calc_hrates!(atmos)
val_e = [20.0, 30.0]  # tests have found ~26 K/day for this setup
val_o = atmos.heating_rate[atmos.nlev_c-2]
println("Expected range = $(val_e) K/day")
println("Modelled value = $(val_o) K/day")
if ( val_o > val_e[1]) && (val_o < val_e[2]) 
    println("Pass")
else
    println("Fail")
    passing = false
end
atmosphere.deallocate!(atmos)


# -------------
# Inform at end
# -------------
println(" ")
if passing
    println("SUCCESS - all tests passed")
    exit(0)
else 
    println("WARNING - some tests failed")
    exit(1)
end

