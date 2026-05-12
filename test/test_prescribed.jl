using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)), "../"))
OUT_DIR = joinpath(ROOT_DIR, "out/")

# helper function for setting up new atmosphere with greygas
function _make_prescribed_atmos(; instellation::Float64=1200.0,
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
    ok || error("Failed to setup prescribed-solver atmosphere")
    atmosphere.allocate!(atmos, ""; check_safe_gas=false) || error("Failed to allocate prescribed-solver atmosphere")
    return atmos
end

@testset "prescribed_solver" begin
    atmos = _make_prescribed_atmos()
    try
        # validation branches
        @test !solver.solve_prescribed!(atmos; sol_type=0, atm_type=1)
        @test !solver.solve_prescribed!(atmos; sol_type=1, atm_type=99)

        # prescribed solver basic solve path
        atmos.tmp_surf = 700.0
        @test solver.solve_prescribed!(atmos; sol_type=1, atm_type=1)
        @test atmos.is_solved
        @test atmos.is_converged
        @test all(atmos.tmp .≈ atmos.tmp_surf)
        @test all(atmos.tmpl .≈ atmos.tmp_surf)
    finally
        atmosphere.deallocate!(atmos)
    end
end
