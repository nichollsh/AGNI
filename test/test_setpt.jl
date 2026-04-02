using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR = joinpath(ROOT_DIR,"res/")
OUT_DIR = joinpath(ROOT_DIR,"out/")

@testset "setpt" begin
    
    # Test _parse_tmp_str with different inputs
    @testset "parse_tmp_str" begin
        # Setup minimal atmosphere for testing
        p_surf = 1.0       # bar
        p_top = 1e-8
        theta = 65.0
        gravity = 10.0
        nlev_centre = 50
        radius = 1.0e7
        tmp_surf = 300.0
        toa_heating = 1000.0
        mf_dict = Dict("H2O" => 1.0)
        spfile = joinpath(RES_DIR, "spectral_files", "Dayspring", "48", "Dayspring.sf")
        
        atmos = AGNI.atmosphere.Atmos_t()
        AGNI.atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                              spfile,
                              toa_heating, 1.0, 0.0, theta,
                              tmp_surf,
                              gravity, radius,
                              nlev_centre, p_surf, p_top,
                              mf_dict, ""
                      )
        
        # Test with numeric string
        result = AGNI.setpt._parse_tmp_str(atmos, "500.0")
        @test isapprox(result, 500.0; atol=1e-10)
        
        # Test with numeric value
        result = AGNI.setpt._parse_tmp_str(atmos, 600.0)
        @test isapprox(result, 600.0; atol=1e-10)
        
        # Test with "tsurf" keyword
        result = AGNI.setpt._parse_tmp_str(atmos, "tsurf")
        @test isapprox(result, tmp_surf; atol=1e-10)
        
        # Test with "teq" keyword
        result = AGNI.setpt._parse_tmp_str(atmos, "teq")
        expected_teq = AGNI.phys.calc_Teq(toa_heating, 0.0)
        @test isapprox(result, expected_teq; rtol=1e-3)
        
        # Test with invalid input (should error)
        @test_throws ErrorException AGNI.setpt._parse_tmp_str(atmos, Dict())
    end
    
end
