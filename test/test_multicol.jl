using Test
using AGNI


ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR         = joinpath(ROOT_DIR,"res/")
OUT_DIR         = joinpath(ROOT_DIR,"out/")
TEST_DIR        = joinpath(ROOT_DIR,"test/")

function _make_grey_atmos(; nlev_c::Int64=100,
                                instellation::Float64=900.0,
                                p_surf::Float64=1.0,
                                p_top::Float64=1e-6,
                                transspec_ref_tau::Float64=1.1,
                                transspec_ref_wl::Float64=1e-6,
                                transspec_ref_p::Float64=1e-3)
    atmos = atmosphere.Atmos_t()
    ok = atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            "greygas",
                            instellation, 1.0, 0.0, 0.0,
                            400.0,
                            10.0, 1.0e7,
                            nlev_c, p_surf, p_top,
                            Dict("N2" => 1.0), "";
                            real_gas=false,
                            thermo_functions=false,
                            flag_rayleigh=false,
                            tmp_magma=1200.0,
                            flag_cloud=false,
                            transspec_ref_tau=transspec_ref_tau,
                            transspec_ref_wl=transspec_ref_wl,
                            transspec_ref_p=transspec_ref_p)
    ok || error("Failed to setup transspec atmosphere")
    atmosphere.allocate!(atmos, ""; check_safe_gas=false) || error("Failed to allocate transspec atmosphere")
    setpt.isothermal!(atmos, 500.0)
    atmosphere.calc_layer_props!(atmos)
    return atmos
end

@testset "multicol" begin

    @testset "test_init_globe" begin

        # make atmosphere object first
        atmos = _make_grey_atmos()

        # then make globe from it
        globe = multicol.Globe_t()
        lons = Float64[0.0, 45.0, 160.0]
        lats = Float64[0.0, 0.0, 0.0]
        redist_flux = Float64[-1.0, -2.0, 3.0]
        redist_Pmid = Float64[1.0, 1.0, 1.0] .* 1e5
        redist_Pwid = Float64[4.0, 4.0, 4.0]


        # construct it
        ok = multicol.construct!(
                                globe, atmos,
                                lons,
                                lats,
                                redist_flux,
                                redist_Pmid,
                                redist_Pwid,
                                )
        @test ok
        @test globe.is_constructed

        # test variables set
        @test all(isapprox.(globe.redist_Pmid, redist_Pmid; rtol=rtol))
        @test all(isapprox.(globe.redist_flux, redist_flux; rtol=rtol))
        @test all(isapprox.(globe.redist_Pwid, redist_Pwid; rtol=rtol))

        # test number of columns
        @test length(globe.atmos_arr) == globe.ncol
        @test length(globe.atmos_arr) == length(lons)

        # test inheret boundary conditions
        @test isapprox(globe.tmp_surf, atmos.tmp_surf; rtol=0.0, atol=1e-6)
        @test isapprox(globe.tmp_magma, atmos.tmp_magma; rtol=0.0, atol=1e-6)
        @test isapprox(globe.flux_int, atmos.flux_int; rtol=0.0, atol=1e-6)

        # test atmosphere columns have correct lat/lon
        for i in eachindex(lons)
            @test isapprox(globe.atmos_arr[i].col_lon, lons[i]; rtol=0.0, atol=1e-6)
            @test isapprox(globe.atmos_arr[i].col_lat, lats[i]; rtol=0.0, atol=1e-6)
        end

        # test zenith angle positive
        for i in eachindex(lons)
            @test globe.atmos_arr[i].zenith_degrees >= 0.0
        end

        # test tear-down
        multicol.deconstruct!(globe)
        @test !globe.is_constructed
    end


    # copy from column to worker
    @testset "copy_atmos_fields" begin

        # init
        atmos = _make_grey_atmos()
        globe = multicol.Globe_t()
        lons = Float64[0.0, 45.0, 160.0]
        lats = Float64[0.0, 0.0, 0.0]
        redist_flux = Float64[-1.0, -2.0, 3.0]
        redist_Pmid = Float64[1.0, 1.0, 1.0] .* 1e5
        redist_Pwid = Float64[4.0, 4.0, 4.0]
        ok = multicol.construct!(
                                globe, atmos,
                                lons,
                                lats,
                                redist_flux,
                                redist_Pmid,
                                redist_Pwid,
                                )
        @test ok

        # copy to worker
        ok = multicol.copy_atmos_fields!(atmos, globe.atmos_arr[1])
        @test ok

        # modify and copy
        globe.atmos_arr[1].tmp_surf += 10.0
        ok = multicol.copy_atmos_fields!(atmos, globe.atmos_arr[1])
        @test ok
        @test isapprox(atmos.tmp_surf, globe.atmos_arr[1].tmp_surf; rtol=0.0, atol=1e-6)

        # tear down
        multicol.deconstruct!(globe)
    end

    @testset "call_agnostic" begin

        # init
        atmos = _make_grey_atmos()
        globe = multicol.Globe_t()
        lons = Float64[0.0, 45.0, 160.0]
        lats = Float64[0.0, 0.0, 0.0]
        redist_flux = Float64[-1.0, -2.0, 3.0]
        redist_Pmid = Float64[1.0, 1.0, 1.0] .* 1e5
        redist_Pwid = Float64[4.0, 4.0, 4.0]
        ok = multicol.construct!(
                                globe, atmos,
                                lons,
                                lats,
                                redist_flux,
                                redist_Pmid,
                                redist_Pwid,
                                )
        @test ok

        # call for globe
        function _set_toa!(atm)
            atm.toa_heating = 123.0
        end
        ok = multicol.call_for_globe!(globe, _set_toa!)
        @test ok
        for i in eachindex(lons)
            @test isapprox(globe.atmos_arr[i].toa_heating, 123.0; rtol=0.0, atol=1e-6)
        end

        # call agnostic on globe
        function _set_albedo!(atm; albedo::Float64=0.5)
            atm.albedo_b = albedo
        end
        ok = multicol.call_agnostic!(globe, _set_albedo!; albedo=0.3)
        @test ok
        for i in eachindex(lons)
            @test isapprox(globe.atmos_arr[i].albedo_b, 0.3; rtol=0.0, atol=1e-6)
        end

        # call agnostic on atmos
        ok = multicol.call_agnostic!(atmos, _set_albedo!; albedo=0.2)
        @test ok
        @test isapprox(atmos.albedo_b, 0.2; rtol=0.0, atol=1e-6)

        # globe columns should be unchanged by call on atmos
        for i in eachindex(lons)
            @test isapprox(globe.atmos_arr[i].albedo_b, 0.3; rtol=0.0, atol=1e-6)
        end

        # tear down
        multicol.deconstruct!(globe)
    end

    @testset "set_surface_bc" begin

        # init
        atmos = _make_grey_atmos()
        globe = multicol.Globe_t()
        lons = Float64[0.0, 45.0, 160.0]
        lats = Float64[0.0, 0.0, 0.0]
        redist_flux = Float64[-1.0, -2.0, 3.0]
        redist_Pmid = Float64[1.0, 1.0, 1.0] .* 1e5
        redist_Pwid = Float64[4.0, 4.0, 4.0]
        ok = multicol.construct!(
                                globe, atmos,
                                lons,
                                lats,
                                redist_flux,
                                redist_Pmid,
                                redist_Pwid,
                                )
        @test ok

        # set new bc: tsurf
        globe.tmp_surf = 401.0
        ok = multicol.set_surface_bc!(globe, 1) # bctype = 1
        @test ok
        for i in eachindex(lons)
            @test isapprox(globe.atmos_arr[i].tmp_surf, 401.0; rtol=0.0, atol=1e-6)
        end

        # set new bc: magma temp
        globe.tmp_magma = 1300.0
        ok = multicol.set_surface_bc!(globe, 2) # bctype = 2
        @test ok
        for i in eachindex(lons)
            @test isapprox(globe.atmos_arr[i].tmp_magma, 1300.0; rtol=0.0, atol=1e-6)
        end

        # set new bc: flux_int
        globe.flux_int = 500.0
        ok = multicol.set_surface_bc!(globe, 3) # bctype = 3
        @test ok
        for i in eachindex(lons)
            @test isapprox(globe.atmos_arr[i].flux_int, 500.0; rtol=0.0, atol=1e-6)
        end

        # failure state: set invalid bc type
        ok = multicol.set_surface_bc!(globe, 999) # invalid bctype
        @test !ok

        # tear down
        multicol.deconstruct!(globe)
    end

    @testset "set_redist" begin

        # init
        atmos = _make_grey_atmos()
        globe = multicol.Globe_t()
        lons = Float64[0.0, 45.0, 160.0]
        lats = Float64[0.0, 0.0, 0.0]
        redist_flux = Float64[-1.0, -2.0, 3.0]
        redist_Pmid = Float64[1.0, 1.0, 1.0] .* 1e5
        redist_Pwid = Float64[4.0, 4.0, 4.0]
        ok = multicol.construct!(
                                globe, atmos,
                                lons,
                                lats,
                                redist_flux,
                                redist_Pmid,
                                redist_Pwid,
                                )
        @test ok

        # set new redist
        ok = multicol.set_redist!(globe)
        @test ok

        # tear down
        multicol.deconstruct!(globe)
    end

    # _is_succ is a pure classifier, no Atmos_t/Globe_t required at all.
    @testset "is_succ" begin
        @test multicol._is_succ(true) == true
        @test multicol._is_succ(false) == false
        @test multicol._is_succ(1.0) == true       # finite Float
        @test multicol._is_succ(NaN) == false      # non-finite Float
        # fallback branch: any type other than Bool/AbstractFloat is treated as
        # success (e.g. a function returning `nothing`, or a plain value)
        @test multicol._is_succ("ok") == true
        @test multicol._is_succ(nothing) == true
    end

    # construct!'s four early-return guards all fire before deepcopy/allocation work
    @testset "construct_guards" begin

        # allocate atmosphere cheaply
        atmos = atmosphere.Atmos_t()
        ok0 = atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                                "greygas",
                                900.0, 1.0, 0.0, 0.0,
                                400.0,
                                10.0, 1.0e7,
                                30, 1.0, 1e-6,
                                Dict("N2" => 1.0), "";
                                real_gas=false,
                                thermo_functions=false,
                                flag_rayleigh=false,
                                flag_cloud=false)
        ok0 || error("Failed to setup test atmosphere")
        @test !atmos.is_alloc

        globe = multicol.Globe_t()
        logs, ok = Test.collect_test_logs() do
            multicol.construct!(globe, atmos, [0.0], [0.0], [0.0], [1.0e5], [4.0])
        end
        @test ok == false
        @test !globe.is_constructed
        @test any(occursin("not allocated", l.message) for l in logs if l.level == Logging.Warn)

        atmos.is_alloc = true  # white-box: only to get past the guard above

        logs, ok = Test.collect_test_logs() do
            multicol.construct!(globe, atmos, [0.0, 45.0], [0.0], [0.0, 0.0], [1.0e5, 1.0e5], [4.0, 4.0])
        end
        @test ok == false
        @test any(occursin("must match num of columns", l.message) for l in logs if l.level == Logging.Warn)

        logs, ok = Test.collect_test_logs() do
            multicol.construct!(globe, atmos, [400.0], [0.0], [0.0], [1.0e5], [4.0])
        end
        @test ok == false
        @test any(occursin("Invalid longitude", l.message) for l in logs if l.level == Logging.Warn)

        logs, ok = Test.collect_test_logs() do
            multicol.construct!(globe, atmos, [0.0], [100.0], [0.0], [1.0e5], [4.0])
        end
        @test ok == false
        @test any(occursin("Invalid latitude", l.message) for l in logs if l.level == Logging.Warn)

        # discrimination guard: none of these failure paths must set is_constructed=true
        @test !globe.is_constructed
    end

end
