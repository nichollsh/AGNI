using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)), "../"))
OUT_DIR = joinpath(ROOT_DIR, "out/")

# helper to create a simple grey gas atmosphere for testing diagnostics and energy utilities
function _make_greygas_atmos(; instellation::Float64=1200.0,
                               nlev_c::Int64=30,
                               p_surf::Float64=10.0,
                               p_top::Float64=1e-6)
    atmos = atmosphere.Atmos_t()
    ok = atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            "greygas",
                            instellation, 1.0, 0.0, 0.0,
                            350.0,
                            10.0, 1.0e7,
                            nlev_c, p_surf, p_top,
                            Dict("N2" => 1.0), "";
                            real_gas=false,
                            thermo_functions=false,
                            flag_rayleigh=false,
                            flag_cloud=false)
    ok || error("Failed to setup test atmosphere")
    atmosphere.allocate!(atmos, ""; check_safe_gas=false) || error("Failed to allocate test atmosphere")
    return atmos
end

@testset "transparent" begin
    atmos = _make_greygas_atmos()
    # transparent solver with type=0, returns false
    @test !solver.solve_transparent!(atmos; sol_type=0)

    # test with transparent flag off, returns false
    atmos.transparent = false
    @test !solver.solve_transparent!(atmos; sol_type=1)

    # now make transparent and try solving with type=1
    atmosphere.make_transparent!(atmos)
    @test !solver.solve_transparent!(atmos; sol_type=1)
    @test atmos.is_solved
    @test !atmos.is_converged

    # try with type=3 and check converges
    atmosphere.make_transparent!(atmos)
    atmos.flux_int = 1200.0
    @test solver.solve_transparent!(atmos; sol_type=3, max_steps=150)
    @test isapprox(atmos.flux_tot[1], atmos.flux_int; rtol=1e-2, atol=0.5)

    # try with type=4 and check converges to target OLR
    atmosphere.make_transparent!(atmos)
    atmos.target_olr = 5000.0
    @test solver.solve_transparent!(atmos; sol_type=4, max_steps=150)
    @test isapprox(atmos.flux_u_lw[1], atmos.target_olr; rtol=1e-2, atol=0.5)
end
