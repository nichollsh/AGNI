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
    @test isapprox(AGNI.consts.proton_mass,  1.67262192595e-27; rtol=1e-9)

    # -------------
    # Dimensionless parameters for atmospheric modeling
    # -------------
    @test AGNI.consts.k_vk > 0.0  # Von Karman constant
    @test AGNI.consts.αMLT > 0.0  # Mixing length parameter
    @test AGNI.consts.βRa  > 0.0  # Convection beta parameter

    # -------------
    # Derived constants
    # -------------
    @test isapprox(AGNI.consts.Cp_ideal, AGNI.consts.R_gas * 7/2; rtol=1e-12)
    @test AGNI.consts.Cp_ideal > 0.0  # Sanity check: positive heat capacity

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

    # Test elements MMW
    element_mmw_cases = [
        ("H",   1.008000e-03),
        ("C",   1.201100e-02),
        ("O",   1.599900e-02),
        ("Fe",  5.584500e-02),
    ]
    for (elem, expected) in element_mmw_cases
        @test isapprox(AGNI.consts._lookup_mmw[elem], expected; rtol=1e-6)
    end

    # All MMW values should be positive
    for (species, mmw) in AGNI.consts._lookup_mmw
        @test mmw > 0.0  # MMW for all species should be positive
    end

    # -------------
    # Atom counts for key molecules
    # -------------
    @test AGNI.consts._lookup_count_atoms["H2O"]["H"] == 2
    @test AGNI.consts._lookup_count_atoms["H2O"]["O"] == 1
    @test AGNI.consts._lookup_count_atoms["CO2"]["C"] == 1
    @test AGNI.consts._lookup_count_atoms["CO2"]["O"] == 2
    @test AGNI.consts._lookup_count_atoms["S8"]["S"] == 8

    # All atom counts should be positive integers
    for (molecule, atoms) in AGNI.consts._lookup_count_atoms
        for (element, count) in atoms
            @test count > 0  # Atom counts should be positive
            @test count == floor(count)  # Atom counts should be integers
        end
    end

    # -------------
    # Plotting colours: known gases and key elements must be present
    # -------------
    for gas in ["H2O", "CO2", "H2"]
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
    # Molecular composition validation
    # For all molecules with defined atom counts, ensure constituent
    # elements are in elems_standard and have assigned colors
    # -------------
    for (molecule, atoms_dict) in AGNI.consts._lookup_count_atoms
        for element in keys(atoms_dict)
            @test element in AGNI.consts.elems_standard
            @test haskey(AGNI.consts._lookup_color, element)
        end
    end

    # Additional check: All single-element species in vols_standard and vaps_standard
    # should be in elems_standard and have colors
    for species in vcat(AGNI.consts.vols_standard, AGNI.consts.vaps_standard)
        # Check if it's a single element (1-2 character string, starts with uppercase)
        if length(species) <= 2 && occursin(r"^[A-Z][a-z]?$", species)
            @test species in AGNI.consts.elems_standard
            @test haskey(AGNI.consts._lookup_color, species)
        end
    end

    # -------------
    # Solar metallicities: hydrogen is anchored at 12.00 by definition
    # -------------
    @test isapprox(AGNI.consts._solar_metallicity["H"], 12.00; atol=1e-10)
    @test haskey(AGNI.consts._solar_metallicity, "Fe")

    # Test key solar metallicity values (important for atmospheric chemistry)
    solar_metal_cases = [
        ("C",  8.46),
        ("N",  7.83),
        ("O",  8.69),
        ("He", 10.914),
        ("Fe", 7.46),
        ("Si", 7.51),
        ("Mg", 7.55),
    ]
    for (elem, expected) in solar_metal_cases
        @test isapprox(AGNI.consts._solar_metallicity[elem], expected; atol=1e-10)
    end

    # Most elements in elems_standard should have solar metallicity values
    # Note: D is currently missing from _solar_metallicity
    elem_with_metals = filter(e -> haskey(AGNI.consts._solar_metallicity, e), AGNI.consts.elems_standard)
    @test length(elem_with_metals) >= 21  # Most elements should have solar metallicity data

    # -------------
    # Liquid densities for ocean module
    # -------------
    @test isapprox(AGNI.consts._lookup_liquid_rho["H2O"], 958.37; rtol=1e-5)
    @test AGNI.consts._lookup_liquid_rho["CO2"] > AGNI.consts._lookup_liquid_rho["H2O"]

    # All liquid densities should be positive
    for (species, rho) in AGNI.consts._lookup_liquid_rho
        @test rho > 0.0  # All liquid densities should be positive
    end

    # -------------
    # Color assignments: all colors should be valid hex codes
    # -------------
    for (species, color) in AGNI.consts._lookup_color
        @test occursin(r"^#[0-9A-Fa-f]{6}$", color)  # Valid hex color format
    end

    # -------------
    # Data consistency checks
    # -------------
    # Species with MMW should include key atmospheric gases
    for gas in ["H2O", "CO2", "H2", "N2", "CH4", "O2"]
        @test haskey(AGNI.consts._lookup_mmw, gas)
    end

    # gases_standard should contain most species without significant duplication
    @test length(AGNI.consts.gases_standard) >= 50  # Should have many species
    duplicate_count = length(AGNI.consts.gases_standard) - length(unique(AGNI.consts.gases_standard))
    @test duplicate_count == 0

    # vols_standard and vaps_standard should not overlap
    vols_set = Set(AGNI.consts.vols_standard)
    vaps_set = Set(AGNI.consts.vaps_standard)
    @test length(intersect(vols_set, vaps_set)) == 0

end
