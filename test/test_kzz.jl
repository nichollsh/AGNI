using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR         = joinpath(ROOT_DIR,"res/")
OUT_DIR         = joinpath(ROOT_DIR,"out/")
TEST_DIR        = joinpath(ROOT_DIR,"test/")

@testset "Kzz" begin
    # -------------
    # Test Kzz type 1: Constant value
    # -------------
    tmp_surf        = 1500.0
    toa_heating     = 5000.0
    p_surf          = 100.0
    p_top           = 1e-6
    theta           = 60.0
    gravity         = 10.0
    nlev_centre     = 100
    radius          = 1.0e7
    Kzz_kbreak      = 1e6  # 1e6 cm2/s
    mf_dict         = Dict([("H2O", 1.0)])
    spfile_name     = "$RES_DIR/spectral_files/Frostflow/48/Frostflow.sf"

    atmos = atmosphere.Atmos_t()
    atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            spfile_name,
                            toa_heating, 1.0, 0.0, theta,
                            tmp_surf,
                            gravity, radius,
                            nlev_centre, p_surf, p_top,
                            mf_dict, "";
                            Kzz_kbreak=Kzz_kbreak,
                            Kzz_type=1
                    )
    atmosphere.allocate!(atmos,"$RES_DIR/stellar_spectra/sun.txt")

    # Set up a simple temperature profile
    setpt.dry_adiabat!(atmos)
    setpt.stratosphere!(atmos, 200.0)

    # Calculate layer properties (density, cp, etc)
    atmosphere.calc_layer_props!(atmos)

    # Run MLT to calculate convection and Kzz
    energy.convection!(atmos)

    # In convective regions, Kzz should be constant = Kzz_kbreak
    for i in 1:atmos.nlev_l
        if atmos.mask_c[i]
            @test isapprox(atmos.Kzz[i], Kzz_kbreak; rtol=1e-6)
        end
    end

    # -------------
    # Test Kzz type 2: Simple scaling (w * λ)
    # -------------

    atmos.Kzz_type = 2

    # Run MLT
    energy.convection!(atmos)

    # In convective regions, Kzz = w_conv * λ_conv
    for i in 1:atmos.nlev_l
        if atmos.mask_c[i]
            expected_Kzz = atmos.λ_conv[i] * atmos.w_conv[i]
            @test isapprox(atmos.Kzz[i], expected_Kzz; rtol=1e-6)
            @test atmos.Kzz[i] > 0.0
        end
    end

    # -------------
    # Test Kzz type 3: Charnay+15 parametrisation
    # -------------

    atmos.Kzz_type = 3

    # Calculate hydrostatic structure
    atmosphere.calc_layer_props!(atmos)

    # Run MLT
    energy.convection!(atmos)

    # In convective regions, check Kzz is positive and finite
    for i in 1:atmos.nlev_l
        if atmos.mask_c[i]
            @test isfinite(atmos.Kzz[i])
            @test atmos.Kzz[i] > 0.0

            # Verify formula: Kzz = (Hp/3) * (λ/Hp)^(4/3) * (R*Fc/(mu*rho*cp))^(1/3)
            # Get interpolated values
            m1 = (atmos.p[i] - atmos.pl[i]) / (atmos.p[i] - atmos.p[i-1])
            m2 = (atmos.pl[i] - atmos.p[i-1]) / (atmos.p[i] - atmos.p[i-1])

            mu   = atmos.layer_μ[i]  * m2 + atmos.layer_μ[i-1]  * m1
            c_p  = atmos.layer_cp[i] * m2 + atmos.layer_cp[i-1] * m1
            rho  = atmos.layer_ρ[i]  * m2 + atmos.layer_ρ[i-1]  * m1
            Hp   = atmos.layer_Hp[i] * m2 + atmos.layer_Hp[i-1] * m1

            expected_Kzz = (Hp/3.0) * (atmos.λ_conv[i]/Hp)^(4.0/3.0) *
                            (phys.R_gas*atmos.flux_cdry[i]/(mu*rho*c_p))^(1.0/3.0)

            @test isapprox(atmos.Kzz[i], expected_Kzz; rtol=0.05)
        end
    end

    # -------------
    # Test surface values are zero after MLT
    # -------------

    # Surface values should be zero
    @test atmos.w_conv[end] == 0.0
    @test atmos.λ_conv[end] == 0.0
    @test atmos.flux_cdry[end] == 0.0



    # -------------
    # Test fill_Kzz!: Extension above convective zone
    # -------------

    # Calculate hydrostatic structure
    atmosphere.calc_layer_props!(atmos)

    # Run MLT and fill_Kzz
    energy.convection!(atmos)
    ok = energy.fill_Kzz!(atmos)
    @test ok

    # Find top of convective region
    Kzz_eps = 1.0e-10
    if any(atmos.Kzz .> Kzz_eps)
        i_top = findfirst(x -> x > Kzz_eps, atmos.Kzz)
        i_bot = findlast(x -> x > Kzz_eps, atmos.Kzz)

        # All Kzz values should be positive
        @test all(atmos.Kzz .> 0.0)

        # Above convective zone: power-law scaling
        # Kzz[i] = Kzz[i_top] * (p[i]/p[i_top])^Kzz_power
        for i in 1:i_top-1
            expected = atmos.Kzz[i_top] * (atmos.pl[i]/atmos.pl[i_top])^atmos.Kzz_power
            @test isapprox(atmos.Kzz[i], expected; rtol=1e-6)
        end

        # Below convective zone: constant
        for i in i_bot+1:atmos.nlev_l
            @test isapprox(atmos.Kzz[i], atmos.Kzz[i_bot]; rtol=1e-10)
        end
    end

    # -------------
    # Test fill_Kzz!: No convection case (use Kzz_pbreak)
    # -------------

    # Set up isothermal profile (no convection)
    setpt.isothermal!(atmos, 500.0)

    # Calculate hydrostatic structure
    atmosphere.calc_layer_props!(atmos)

    # Run MLT (should find no convection)
    energy.convection!(atmos)

    # Manually zero out Kzz to simulate no convection
    atmos.Kzz .= 0.0

    # Run fill_Kzz - should use Kzz_pbreak as reference
    ok = energy.fill_Kzz!(atmos)
    @test ok

    # All Kzz values should be positive after fill
    @test all(atmos.Kzz .> 0.0)

    # Find index near Kzz_pbreak (1 bar = 1e5 Pa by default)
    i_ref = findmin(abs.(atmos.pl .- atmos.Kzz_pbreak))[2]
    @test isapprox(atmos.Kzz[i_ref], Kzz_kbreak; rtol=1e-6)

    atmosphere.deallocate!(atmos)
end
