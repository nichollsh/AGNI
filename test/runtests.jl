#!/usr/bin/env -S julia --color=yes --startup-file=no
# Run this function from inside the `test/` folder

# Get AGNI root directory
ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))

# Activate environment
ENV["GKSwstype"] = "100"
import Pkg
Pkg.activate(ROOT_DIR)

# Include libraries
using LoggingExtras
using NCDatasets
using AGNI

@info "Begin AGNI tests"

# Prepare
RES_DIR         = joinpath(ROOT_DIR,"res/")
OUT_DIR         = joinpath(ROOT_DIR,"out/")
TEST_DIR        = joinpath(ROOT_DIR,"test/")
p_top           = 1e-8
nlev_centre     = 100
radius          = 1.0e7    # metres
gravity         = 10.0      # m s-2
total  = 0
failed = 0

# which test suite to run?
# 0 - none
# 10 - fast
# 20 - all
suite::Int64 = 20
if length(ARGS)>0
    if ARGS[1] == "all"
        suite = 20
    elseif ARGS[1] == "fast"
        suite = 10
    elseif ARGS[1] == "none"
        suite = 0
    else
        if !isnothing(tryparse(Int64, ARGS[1]))
            suite = parse(Int64, ARGS[1])
        else
            @warn "Invalid test suite option '$(ARGS[1])'. Running all tests."
            suite = 20
        end
    end
end
suite = min(max(suite, 0), 20)
@info "Using suite $suite"

rm(OUT_DIR,force=true,recursive=true)
if !isdir(OUT_DIR) && !isfile(OUT_DIR)
    mkdir(OUT_DIR)
end

rtol   = 1e-3

# -------------
# Test module imported
# -------------
if suite >= 0
    @info " "
    @info "Testing module imported"
    if isdefined(AGNI.atmosphere, :setup!)
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"
end


if suite > 2
    # -------------
    # Test thermodynamic lookup data validity
    # -------------
    @info " "
    @info "Testing thermodynamics data validity"
    ideal_H2O::phys.Gas_t = phys.load_gas("$RES_DIR/thermodynamics/", "H2O", true, false)
    if !ideal_H2O.fail
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"


    # -------------
    # Test heat capacity lookup tables
    # -------------
    @info " "
    @info "Testing heat capacity functions for H2O"
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
end

# -------------
# Test AQUA equation of state
# -------------
if suite > 5
    @info " "
    @info "Testing H2O AQUA equation of state"
    aqua_H2O::phys.Gas_t = phys.load_gas("$RES_DIR/thermodynamics/", "H2O", true, true)
    v_expt = [926.12116198786, 0.007222354920, 0.4333412952269, 17.038999553692, 66.87150907049]
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
end


# -------------
# Test VdW equation of state
# -------------
if suite > 4
    @info " "
    @info "Testing CO2 VdW equation of state"
    vdw_CO2::phys.Gas_t = phys.load_gas("$RES_DIR/thermodynamics/", "CO2", true, true)
    v_expt = [2.646533036586e-5, 0.01764355357724, 1.061227115599, 41.2385376189, 147.631823888]
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
end


# -------------
# Test mixing ratios
# -------------
if suite > 6
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
    spfile_name   ="$RES_DIR/spectral_files/Dayspring/48/Dayspring.sf"

    # Setup atmosphere
    atmos = atmosphere.Atmos_t()
    atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
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
end



if suite > 7
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
    spfile_name   = "$RES_DIR/spectral_files/Dayspring/48/Dayspring.sf"

    # Setup atmosphere
    atmos = atmosphere.Atmos_t()
    atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
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
    atmosphere.allocate!(atmos,"$RES_DIR/stellar_spectra/sun.txt")
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
end


if suite > 7
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
    spfile_name   = "$RES_DIR/spectral_files/Oak/318/Oak.sf"

    # Setup atmosphere
    atmos = atmosphere.Atmos_t()
    atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            spfile_name,
                            toa_heating, 1.0, 0.0, theta,
                            tmp_surf,
                            gravity, radius,
                            300, p_surf, p_top,
                            mf_dict,"",
                            flag_gcontinuum=true,
                            flag_rayleigh=false,
                            overlap_method="ee",
                            surface_material="greybody",
                            albedo_s=0.5,
                            real_gas=false
                    )
    atmosphere.allocate!(atmos,"$RES_DIR/stellar_spectra/sun.txt")
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

    val_e = 1.0438945347722374e7   # known from previous tests
    val_o = atmos.r[1] # height of topmost layer-centre
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
    out_path::String = joinpath(OUT_DIR,"agni_atm.nc")
    rm(out_path, force=true)
    save.write_ncdf(atmos, out_path)
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
    out_path = joinpath(OUT_DIR,"agni_plot_tmp.png")
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
    # Test plot radius
    # -------------
    @info " "
    @info "Testing plot radius"
    out_path = joinpath(OUT_DIR,"agni_plot_rad.png")
    plotting.plot_radius(atmos, out_path)
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
    out_path = joinpath(OUT_DIR,"agni_plot_alb.png")
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
end


# -------------
# Test Rayleigh scattering
# -------------
if suite > 12
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
    spfile_name   = "$RES_DIR/spectral_files/Dayspring/48/Dayspring.sf"

    # Setup atmosphere
    atmos = atmosphere.Atmos_t()
    atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
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
    atmosphere.allocate!(atmos,"$RES_DIR/stellar_spectra/sun.txt")
    energy.radtrans!(atmos, true)
    energy.radtrans!(atmos, false)

    val_e = 37.125717835788286
    val_o = atmos.flux_u_sw[20]
    @info "Expected value = $(val_e) W m-2"
    @info "Modelled value = $(val_o) W m-2"
    if abs(val_o - val_e)/val_e < rtol
        @info "Pass"
    else
        @warn "Fail"
        failed  += 1
    end
    total += 1
    atmosphere.deallocate!(atmos)
    @info "--------------------------"
end


if suite > 9
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
    spfile_name   = "$RES_DIR/spectral_files/Oak/318/Oak.sf"

    # Setup atmosphere
    atmos = atmosphere.Atmos_t()
    atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
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
    atmosphere.allocate!(atmos,"$RES_DIR/stellar_spectra/sun.txt")
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
    energy.calc_fluxes!(atmos, radiative=true, convective=true, conductive=true, sens_heat=true, latent_heat=true, advective=true)
    val_e = 8602.78747109532  # from previous tests
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
    @info "--------------------------"

    # -------------
    # Transparent atmosphere solver (sol_type = 4)
    # -------------
    @info " "
    @info "Testing transparent solver (sol_type=4)"
    atmosphere.make_transparent!(atmos)
    atmos.target_olr = 5000.0
    solver.solve_transparent!(atmos; sol_type=4)
    val_e = atmos.target_olr
    val_o = atmos.flux_u_lw[1]
    @info "Expected value = $(val_e) W m-2"
    @info "Modelled value = $(val_o) W m-2"
    if abs(val_o - val_e)/val_e < rtol
        @info "Pass"
    else
        @warn "Fail"
        show(atmos.flux_u_lw)
        failed += 1
    end
    total += 1
    @info "--------------------------"

    # -------------
    # Transparent atmosphere solver (sol_type = 3)
    # -------------
    if suite > 10
        @info " "
        @info "Testing transparent solver (sol_type=3)"
        atmosphere.make_transparent!(atmos)
        atmos.flux_int = 1200.0
        solver.solve_transparent!(atmos; sol_type=3)
        val_e = atmos.flux_int
        val_o = atmos.flux_tot[1]
        @info "Expected value = $(val_e) W m-2"
        @info "Modelled value = $(val_o) W m-2"
        if abs(val_o - val_e)/val_e < rtol
            @info "Pass"
        else
            @warn "Fail"
            show(atmos.flux_tot)
            failed += 1
        end
        total += 1
        @info "--------------------------"
    end
end

if suite > 15
    # -------------
    # Run example TOML config file (writes to output folder)
    # -------------
    @info " "
    @info "Testing model with energy-conserving TP solver"
    cfg = AGNI.open_config(joinpath(TEST_DIR, "test.toml"))

    # check return code is fine
    succ = AGNI.run_from_config(cfg)
    if succ
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    # -------------
    # Compare result from NetCDF
    # -------------
    @info " "
    @info "Testing values from model solution"

    # Variable to store expected values
    arr_e::Array{Float64,1}     = zeros(Float64, cfg["execution"]["num_levels"]+1)
    arr_o_tmp::Array{Float64,1} = zero(arr_e)
    arr_o_prs::Array{Float64,1} = zero(arr_e)
    arr_o_flN::Array{Float64,1} = zero(arr_e)
    arr_o_flC::Array{Float64,1} = zero(arr_e)
    arr_o_Kzz::Array{Float64,1} = zero(arr_e)

    atol = Float64(cfg["execution"]["converge_atol"])

    # Read result
    @debug "ALL DEBUG SUPPRESSED"
    with_logger(MinLevelLogger(current_logger(), Logging.Info-200)) do

        # open netcdf
        ds = Dataset(joinpath(OUT_DIR,"atm.nc"),"r")

        # read profiles from netCDF file
        arr_o_tmp[:] .= ds["tmpl"][:]
        arr_o_prs[:] .= ds["pl"][:]
        arr_o_flN[:] .= ds["fl_N"][:]
        arr_o_flC[:] .= ds["fl_cnvct"][:]
        arr_o_Kzz[:] .= ds["Kzz"][:]

        # close netcdf
        close(ds)
    end
    @debug "ALL DEBUG RESTORED"

    # pressure profile
    @info "Pressure..."
    arr_e[:] .= [1.0, 1.61372830078497, 2.60411902875434, 4.20234057531355, 6.78143591592047, 10.9433950574805, 17.6596663109266, 28.4979033083612, 45.9878730817361, 74.2119322849048, 119.757895384089, 193.256705023749, 311.863814213277, 503.263462986711, 812.130492972705, 1310.5579604405, 2114.88447058187, 3412.8489230686, 5507.4108934593, 8887.46482282671, 14341.9535068263, 23144.0162625079, 37348.1540366365, 60269.7731509967, 97259.0386156535, 156949.663121218, 253274.113177377, 408715.604290549, 659555.937616089, 1064344.08248185, 1717562.16767397, 2771678.67833306, 4472736.32390834, 7217781.18783982, 11647537.7716905, 18795961.3366388, 30331574.7493941, 48946920.5804721, 78987030.9769821, 127463607.282535, 205691630.391968, 331930405.19812, 535645488.759229, 646902370.004265, 648216250.228125 ]
    if all(abs.(arr_e[:] .- arr_o_prs[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        @info "Pass"
    else
        @warn "Fail \n Observed: $(arr_o_prs)"
        failed += 1
    end
    total += 1

    # temperature profile
    @info "Temperature..."
    arr_e[:] .= [160.663417882603, 172.097128812409, 178.036068253899, 179.208908569656, 179.910164720311, 179.281972493295, 179.483379941445, 181.607889256078, 183.112356682626, 182.74338943682, 183.440372062087, 186.233959168495, 190.338756315127, 195.355039534191, 201.103646189501, 207.489881120782, 214.2815536312, 221.705657433829, 231.687457669314, 246.481687372921, 263.356595002906, 287.5630763156, 320.301515497962, 354.852774535303, 392.118464544988, 429.26259724383, 460.588392363893, 487.062314373208, 514.455301328919, 543.168345970076, 570.387595711342, 593.846621505299, 610.931353333751, 620.771283792757, 625.59828953617, 628.122297051458, 629.66973541692, 630.570749507789, 630.967916271267, 631.091540930136, 631.121574633071, 631.128221749457, 631.12965996539, 631.129771581873, 631.129771953692]
    if all(abs.(arr_e[:] .- arr_o_tmp[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        @info "Pass"
    else
        @warn "Fail \n Observed: $(arr_o_tmp)"
        failed += 1
    end
    total += 1

    # radiative flux
    @info "Radiative flux..."
    arr_e[:] .= [-4.9872131683059706e-5, -0.00021685517390324094, -5.446596077263166e-5, -6.956502420507604e-5, -4.464714481855481e-5, -4.1937900391531e-5, -6.189397629441373e-5, -8.169547623992912e-5, -4.90136637267824e-5, -4.301526104200093e-5, -7.095810627788524e-5, -0.00037362132775342616, -0.0016387527005576885, -0.005163568319403566, -0.015416830214064703, -0.04250362273199926, -0.10932251298316942, -0.30647951633102366, -1.4784710919412305, -9.927186737124345, -38.40847329339533, -39.38202174156268, -40.03873661212435, -24.729291535171853, -8.090326572761398, -0.0002255708594987027, -0.0001643575964322963, -0.00015236288642483942, -0.00017095811467982003, -0.00016738894379386693, -0.0001507433066763042, -0.00012077610864480448, -7.536869504320975e-5, -3.700287694208271e-5, -1.7974538663700912e-5, -1.075028051067406e-5, -6.867279150668537e-6, -3.3999808010110044e-6, -1.1308877208085488e-6, -2.810934676129688e-7, -6.233177590498777e-8, -1.3682293434019401e-8, -4.055073132242046e-8, 2.369779394782707e-8, -2.6217860067846884e-5]
    if all(abs.(arr_e[:] .- arr_o_flN[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        @info "Pass"
    else
        @warn "Fail \n Observed: $(arr_o_flN)"
        failed += 1
    end
    total += 1

    # convective flux
    @info "Convective flux..."
    arr_e[:] .= [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1758088369679132, 21.882272290150052, 19.522283264947543, 8.769662805655303, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if all(abs.(arr_e[:] .- arr_o_flC[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        @info "Pass"
    else
        @warn "Fail \n Observed: $(arr_o_flC)"
        failed += 1
    end


    # eddy diffusion coefficients
    @info "Kzz..."
    arr_e[:] .= [293515.330261667, 242381.028037358, 200155.006214063, 165285.323017843, 136490.40572033, 112711.948729345, 93076.0174630662, 76860.9284503441, 63470.7251477871, 52413.2746248124, 43282.1800995569, 35741.8445533194, 29515.1364634408, 24373.2043307382, 20127.0656526946, 16620.6612101877, 13725.1193905086, 11334.0197421416, 9359.48168174701, 7728.93460077965, 6382.45066279981, 5270.54226322015, 13528.8039526634, 12294.5144847605, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378]
    if all(abs.(arr_e[:] .- arr_o_Kzz[:]) .< abs.(arr_e[:].*rtol) .+ 1e4)
        @info "Pass"
    else
        @warn "Fail \n Observed: $(arr_o_Kzz)"
        failed += 1
    end
    total += 1
    @info "--------------------------"
end


# -------------
# Inform at end
# -------------
@info " "
if failed == 0
    @info "All $total tests have passed (suite $suite)"
    exit(0)
else
    @warn "Some tests failed ($failed/$total)"
    exit(1)
end

