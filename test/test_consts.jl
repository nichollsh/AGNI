using Test
using AGNI

@testset "consts" begin
    @test isapprox(AGNI.consts.R_gas, 8.314462618; atol=1e-12)
    @test isapprox(AGNI.consts._lookup_mmw["H2O"], 1.801530E-02; rtol=1e-8)
    @test AGNI.consts._lookup_count_atoms["H2O"]["H"] == 2
    @test "H2O" in keys(AGNI.consts._lookup_color)
end
