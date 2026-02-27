using Test
using AGNI

@testset "consts" begin

    # -------------
    # Fundamental physical constants (NIST CODATA values)
    # Update these if the source values are changed.
    # -------------
    @test isapprox(AGNI.consts.R_gas,        8.314462618;      atol=1e-12)
    @test isapprox(AGNI.consts.σSB,          5.670374419e-8;   rtol=1e-9)
    @test isapprox(AGNI.consts.k_B,          1.380649e-23;     rtol=1e-9)
    @test isapprox(AGNI.consts.h_pl,         6.62607015e-34;   rtol=1e-9)
    @test isapprox(AGNI.consts.c_vac,        299792458.0;      atol=0.1)
    @test isapprox(AGNI.consts.G_grav,       6.67430e-11;      rtol=1e-5)
    @test isapprox(AGNI.consts.zero_celcius, 273.15;           atol=1e-10)

    # -------------
    # Mean molecular weights [kg mol-1] for key gases
    # -------------
    mmw_cases = [
        ("H2O", 1.801530e-02),
        ("CO2", 4.401000e-02),
        ("N2",  2.801340e-02),
        ("H2",  2.015880e-03),
    ]
    for (gas, expected) in mmw_cases
        @test isapprox(AGNI.consts._lookup_mmw[gas], expected; rtol=1e-6)
    end

    # -------------
    # Atom counts for key molecules
    # -------------
    @test AGNI.consts._lookup_count_atoms["H2O"]["H"] == 2
    @test AGNI.consts._lookup_count_atoms["H2O"]["O"] == 1
    @test AGNI.consts._lookup_count_atoms["CO2"]["C"] == 1
    @test AGNI.consts._lookup_count_atoms["CO2"]["O"] == 2

    # -------------
    # Plotting colours: known gases and key elements must be present
    # -------------
    for gas in ["H2O", "CO2", "H2", "N2", "CH4"]
        @test haskey(AGNI.consts._lookup_color, gas)
    end
    for elem in ["H", "C", "O", "N", "S", "Fe", "Si"]
        @test haskey(AGNI.consts._lookup_color, elem)
    end

    # -------------
    # Standard species lists are non-empty and self-consistent
    # -------------
    @test length(AGNI.consts.vols_standard) > 0
    @test length(AGNI.consts.vaps_standard) > 0
    @test length(AGNI.consts.gases_standard) == length(AGNI.consts.vols_standard) +
                                                length(AGNI.consts.vaps_standard)
    @test "H2O" in AGNI.consts.vols_standard
    @test "H"   in AGNI.consts.elems_standard

    # -------------
    # Solar metallicities: hydrogen is anchored at 12.00 by definition
    # -------------
    @test isapprox(AGNI.consts._solar_metallicity["H"], 12.00; atol=1e-10)
    @test haskey(AGNI.consts._solar_metallicity, "Fe")

    # -------------
    # Liquid densities for ocean module
    # -------------
    @test isapprox(AGNI.consts._lookup_liquid_rho["H2O"], 958.37; rtol=1e-5)
    @test AGNI.consts._lookup_liquid_rho["CO2"] > AGNI.consts._lookup_liquid_rho["H2O"]

end
