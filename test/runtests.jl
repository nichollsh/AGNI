#!/usr/bin/env -S julia --color=yes --startup-file=no

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
                            condensates=["H2O"],
                            surface_material="greybody",
                            albedo_s=0.5,
                            real_gas=false
                    )
    atmosphere.allocate!(atmos,"$RES_DIR/stellar_spectra/sun.txt")
    chemistry.regrid_saturated_surf!(atmos)
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
    energy.calc_fluxes!(atmos, radiative=true, convect=true, conduct=true, sens_heat=true, latent_heat=true, advect=true)
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
    arr_e::Array{Float64,1} = zeros(Float64, cfg["execution"]["num_levels"]+1)
    arr_o_tmp::Array{Float64,1} = zero(arr_e)
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
        arr_o_flN[:] .= ds["fl_N"][:]
        arr_o_flC[:] .= ds["fl_cnvct"][:]
        arr_o_Kzz[:] .= ds["Kzz"][:]

        # close netcdf
        close(ds)
    end
    @debug "ALL DEBUG RESTORED"

    # temperature profile
    @info "Temperature..."
    arr_e[:] .= [160.5675648884, 171.994064313637, 177.938575310175, 179.113900851404,  179.801478507666, 179.167818376284, 179.40614826252, 181.556769798685,  183.007095076024, 182.624208591575, 183.399201292413, 186.261911326026,  190.412634382101, 195.47048174129, 201.251679339825, 207.654691790094,  214.451046157178, 221.848979602947, 231.737439525989, 246.412351141343,  263.513451113498, 288.073705588398, 320.936573630417, 355.668949498865,  393.166000909574, 430.332142915339, 461.765909449321, 488.839407157943,  516.896086763742, 545.950964338871, 573.226076102865, 596.498543067942,  613.203600564291, 622.65575000879, 627.239385241944, 629.635777119968,  631.102258523775, 631.936413339924, 632.289738849541, 632.396879253831,  632.422477552567, 632.428125540173, 632.429337440917, 632.429439784017,  632.429440335934]
    if all(abs.(arr_e[:] .- arr_o_tmp[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        @info "Pass"
    else
        @warn "Fail \n Observed: $(arr_o_tmp)"
        failed += 1
    end
    total += 1

    # radiative flux
    @info "Radiative flux..."
    arr_e[:] .= [-1.70936355061713e-05, -0.0001848507873774, -2.19587959691125e-05, -3.68052697012899e-05, -1.17535213348674e-05, -9.21098234130113e-06, -2.99803821803835e-05, -4.89456000991595e-05, -1.5358869575266e-05, -1.09188006263139e-05, -3.93173775137257e-05, -0.000354796292072024, -0.00165987941761614, -0.00528027329954739, -0.0158055913856288, -0.0433923529017761, -0.11140604524968, -0.308141836228856, -1.46423434726591, -9.77401066606802, -37.3100359467623, -38.7488777778792, -39.2155408656643, -23.6328050591985, -7.23626887173242, -0.00023511365176887, -0.000176458639629118, -0.000168000405722069, -0.000180597452342113, -0.000171841531511063, -0.000151788906691763, -0.00011949688596502, -7.29475105316624e-05, -3.52666797951429e-05, -1.70943386361699e-05, -1.02762551650315e-05, -6.48362776556333e-06, -3.05875096534569e-06, -9.87984617728532e-07, -2.40686721166353e-07, -5.31945902250514e-08, -1.16863944068938e-08, -1.6438207365123e-08, 1.3533281159821e-09, -1.25751721498093e-05 ]
    if all(abs.(arr_e[:] .- arr_o_flN[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        @info "Pass"
    else
        @warn "Fail \n Observed: $(arr_o_flN)"
        failed += 1
    end
    total += 1

    # convective flux
    @info "Convective flux..."
    arr_e[:] .= [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.70370297292099, 22.4258287268845, 17.7152047690876, 7.20391869209392, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    if all(abs.(arr_e[:] .- arr_o_flC[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        @info "Pass"
    else
        @warn "Fail \n Observed: $(arr_o_flC)"
        failed += 1
    end


    # eddy diffusion coefficients
    @info "Kzz..."
    arr_e[:] .= [343861.803277211, 283750.654425665, 234147.652108032, 193215.846845093, 159439.40985939, 131567.497342447, 108567.927918327, 89588.9578396376, 73927.7383356734, 61004.2869927183, 50340.0092478434, 41539.9745820394, 34278.2910464086, 28286.0364958031, 23341.2995869213, 19260.9617288462, 15893.9156467514, 13115.4694216393, 10822.7287707491, 8930.78656048291, 7369.57843796924, 6081.2881357725, 13638.9561240586, 11895.0041403355, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682, 8254.94954692682]
    if all(abs.(arr_e[:] .- arr_o_Kzz[:]) .< abs.(arr_e[:].*rtol) .+ 1e3)
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

