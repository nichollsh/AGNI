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
                            check_integrity=false,
                            flag_rayleigh=false,
                            flag_cloud=false)
    ok || error("Failed to setup prescribed-solver atmosphere")
    atmosphere.allocate!(atmos, ""; check_safe_gas=false) || error("Failed to allocate prescribed-solver atmosphere")
    return atmos
end

@testset "prescribed" begin
    atmos = _make_prescribed_atmos()
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

    atmosphere.deallocate!(atmos)

    # sol_type=2: conductive skin-flux boundary condition. Solver should find a Tsurf
    # such that the radiative flux balances the conductive flux through the magma skin
    # layer (energy-balance invariant), rather than just returning without erroring.
    atmos2 = _make_prescribed_atmos()
    @test solver.solve_prescribed!(atmos2; sol_type=2, atm_type=1)
    @test atmos2.is_solved
    @test atmos2.is_converged
    F_skin = energy.skin_flux(atmos2)
    @test isapprox(atmos2.flux_tot[1], F_skin; rtol=1e-2, atol=0.5)
    atmosphere.deallocate!(atmos2)

    # sol_type=3: prescribed total internal flux boundary condition. Solver should find
    # a Tsurf such that the total outgoing flux matches the requested internal flux.
    atmos3 = _make_prescribed_atmos()
    atmos3.flux_int = 1200.0
    @test solver.solve_prescribed!(atmos3; sol_type=3, atm_type=1)
    @test atmos3.is_solved
    @test atmos3.is_converged
    @test isapprox(atmos3.flux_tot[1], atmos3.flux_int; rtol=1e-2, atol=0.5)
    atmosphere.deallocate!(atmos3)
end
