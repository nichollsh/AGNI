using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR  = joinpath(ROOT_DIR,"res/")
OUT_DIR  = joinpath(ROOT_DIR,"out/")

const DH_GRAVITY = 10.0
const DH_RADIUS = 1.0e7

function _make_greygas_atmos(; instellation::Float64=1000.0,
                               nlev_c::Int=30,
                               p_surf::Float64=10.0,
                               p_top::Float64=1e-6)
    atmos = atmosphere.Atmos_t()
    mf_dict = Dict([("N2", 1.0)])

    ok = atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            "greygas",
                            instellation, 1.0, 0.0, 0.0,
                            300.0,
                            DH_GRAVITY, DH_RADIUS,
                            nlev_c, p_surf, p_top,
                            mf_dict, "";
                            real_gas=false,
                            thermo_functions=false,
                            flag_rayleigh=false,
                            flag_cloud=false)
    ok || error("Failed to setup greygas atmosphere for deep heating tests")

    ok = atmosphere.allocate!(atmos, ""; check_safe_gas=false)
    ok || error("Failed to allocate greygas atmosphere for deep heating tests")

    return atmos
end

@testset "deep_heating" begin
    @testset "mass_norm_rel" begin
        atmos = _make_greygas_atmos(instellation=1000.0)

        ok = atmosphere.set_deep_heating!(atmos,
                                            1.0e5, 0.7,
                                            0.1, 0.0,
                                            "mass", "clamp", "rel")
        @test ok

        energy.deep_heating!(atmos)
        @test isapprox(atmos.flux_deep[end], -100.0; rtol=1e-6, atol=1e-8)
        @test all(diff(atmos.flux_deep) .<= 0.0)
        @test atmos.flux_deep[1] == 0.0

        atmosphere.deallocate!(atmos)
    end

    @testset "boundary_flux_abs" begin
        atmos = _make_greygas_atmos(instellation=500.0)

        ok = atmosphere.set_deep_heating!(atmos,
                                            1.0e9, 0.5,
                                            0.0, 250.0,
                                            "pressure", "boundary_flux", "abs")
        @test ok

        energy.calc_fluxes!(atmos; deep=true)
        @test all(isapprox.(atmos.flux_deep, -250.0; rtol=0.0, atol=1e-10))

        diff = atmos.flux_deep .- atmos.flux_tot
        @test all(isapprox.(diff, 0.0; atol=1e-10))

        atmosphere.deallocate!(atmos)
    end

    @testset "power_off" begin
        atmos = _make_greygas_atmos(instellation=800.0)

        ok = atmosphere.set_deep_heating!(atmos,
                                            1.0e5, 0.5,
                                            1.0, 999.0,
                                            "mass", "clamp", "off")
        @test ok

        energy.deep_heating!(atmos)
        @test all(atmos.flux_deep .== 0.0)

        atmosphere.deallocate!(atmos)
    end
end
