using Test
using AGNI

# Ocean module: pure geometric functions for distributing condensed liquids
# across ocean basins and continental shelves.
# Liquid densities used below come from consts._lookup_liquid_rho:
#   H2O = 958.37 kg/m^3
#   CO2 = 1178.4 kg/m^3  (denser → sinks to bottom)

@testset "ocean" begin

    # Shared planetary parameters (easy to update for a different test planet)
    fOB = 0.3      # ocean basin area fraction
    hCS = 1000.0   # continental shelf height [m]
    Rpl = 6.371e6  # planet radius [m]

    # -------------
    # No liquids → empty layers, no coverage
    # -------------
    @testset "empty" begin
        sigs   = Dict{String, Float64}()
        layers = AGNI.ocean.dist_surf_liq(sigs, fOB, hCS, Rpl)
        @test length(layers) == 0
        @test AGNI.ocean.get_topliq(layers)          == "none"
        @test AGNI.ocean.get_maxdepth(layers)        == 0.0
        @test AGNI.ocean.get_areacov(layers, fOB)    == 0.0
    end

    # -------------
    # Small amount of H2O fits entirely within the ocean basins
    # Expected: one layer, area coverage = fOB (continental planet)
    # -------------
    @testset "basin_fill" begin
        sigs   = Dict("H2O" => 100.0)  # 100 kg/m^2 of rainout
        layers = AGNI.ocean.dist_surf_liq(sigs, fOB, hCS, Rpl)
        @test length(layers) == 1
        @test AGNI.ocean.get_topliq(layers)          == "H2O"
        @test AGNI.ocean.get_maxdepth(layers)        > 0.0
        @test isapprox(AGNI.ocean.get_areacov(layers, fOB), fOB; atol=1e-10)
    end

    # -------------
    # Very large amount of H2O overflows basins onto the continental shelf
    # Expected: area coverage = 1.0 (aqua planet)
    # -------------
    @testset "planet_flood" begin
        sigs   = Dict("H2O" => 1.0e6)  # huge amount; greatly exceeds basin capacity
        layers = AGNI.ocean.dist_surf_liq(sigs, fOB, hCS, Rpl)
        @test AGNI.ocean.get_topliq(layers)       == "H2O"
        @test AGNI.ocean.get_areacov(layers, fOB) == 1.0
    end

    # -------------
    # Two liquids: CO2 is denser (1178 kg/m^3) so it sinks below H2O (958 kg/m^3).
    # The top of the ocean is the least-dense liquid exposed to the atmosphere.
    # -------------
    @testset "two_liquids" begin
        sigs   = Dict("H2O" => 200.0, "CO2" => 100.0)
        layers = AGNI.ocean.dist_surf_liq(sigs, fOB, hCS, Rpl)
        @test length(layers) == 2
        @test AGNI.ocean.get_topliq(layers) == "H2O"   # H2O floats on CO2
        @test AGNI.ocean.get_maxdepth(layers) > 0.0
    end

    # -------------
    # Depth is additive over layers (basin + shelf components)
    # -------------
    @testset "depth_additive" begin
        # Place enough liquid that both in-basin and shelf portions are non-zero
        sigs   = Dict("H2O" => 1.0e5)
        layers = AGNI.ocean.dist_surf_liq(sigs, fOB, hCS, Rpl)
        depth  = AGNI.ocean.get_maxdepth(layers)
        manual = sum(la[3] + la[4] for la in layers)
        @test isapprox(depth, manual; atol=1e-10)
    end

end
