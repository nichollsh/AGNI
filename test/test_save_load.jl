using Test
using AGNI
using NCDatasets

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR = joinpath(ROOT_DIR,"res/")
OUT_DIR = joinpath(ROOT_DIR,"out/")

@testset "save" begin
    
    # Setup a minimal atmosphere for testing
    # This is needed once and reused across tests
    p_surf = 1.0       # bar
    p_top = 1e-8
    theta = 65.0
    gravity = 10.0
    nlev_centre = 50   # Increase from 20 to meet minimum requirement
    radius = 1.0e7
    tmp_surf = 500.0
    toa_heating = 1000.0
    mf_dict = Dict("H2O" => 0.8, "CO2" => 0.2)
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
    AGNI.atmosphere.allocate!(atmos, "")
    
    # Set some T(p) profile for testing
    AGNI.setpt.isothermal!(atmos, tmp_surf)
    
    # Initialize flux arrays to zero (normally done by energy module)
    fill!(atmos.flux_u_lw, 0.0)
    fill!(atmos.flux_d_lw, 0.0)
    fill!(atmos.flux_n_lw, 0.0)
    fill!(atmos.flux_u_sw, 0.0)
    fill!(atmos.flux_d_sw, 0.0)
    fill!(atmos.flux_n_sw, 0.0)
    fill!(atmos.flux_u, 0.0)
    fill!(atmos.flux_d, 0.0)
    fill!(atmos.flux_n, 0.0)
    fill!(atmos.flux_cdry, 0.0)
    fill!(atmos.flux_l, 0.0)
    fill!(atmos.flux_deep, 0.0)
    fill!(atmos.flux_tot, 0.0)
    
    # Test write_profile (CSV output)
    @testset "write_profile" begin
        csv_file = joinpath(OUT_DIR, "test_profile.csv")
        
        # Remove file if it exists
        rm(csv_file, force=true)
        
        # Write profile
        result = AGNI.save.write_profile(atmos, csv_file)
        
        # Check that function returns nothing
        @test result === nothing
        
        # Check that file was created
        @test isfile(csv_file)
        
        # Read and validate file content
        lines = readlines(csv_file)
        @test length(lines) > 2  # At least header + some data
        @test occursin("pressure", lines[1])
        @test occursin("temperature", lines[1])
        @test occursin("radius", lines[1])
        
        # Check that we have the right number of data lines
        # (header comments + data = nlev_l + nlev_c + 2 header lines)
        @test length(lines) == atmos.nlev_l + atmos.nlev_c + 2
        
        # Validate a data line format (should have 3 comma-separated values)
        data_line = lines[3]  # First data line after headers
        @test occursin(",", data_line)
        values = split(data_line, ",")
        @test length(values) == 3
        
        # Clean up
        rm(csv_file, force=true)
    end
    
    # Test write_fluxes (CSV output)
    @testset "write_fluxes" begin
        csv_file = joinpath(OUT_DIR, "test_fluxes.csv")
        
        # Remove file if it exists
        rm(csv_file, force=true)
        
        # Write fluxes
        result = AGNI.save.write_fluxes(atmos, csv_file)
        
        # Check that function returns nothing
        @test result === nothing
        
        # Check that file was created
        @test isfile(csv_file)
        
        # Read and validate file content
        lines = readlines(csv_file)
        @test length(lines) > 2
        @test occursin("pressure", lines[1])
        @test occursin("U_LW", lines[1])
        @test occursin("D_LW", lines[1])
        @test occursin("N_LW", lines[1])
        @test occursin("U_SW", lines[1])
        @test occursin("D_SW", lines[1])
        
        # Check number of data lines (nlev_l + 2 header lines)
        @test length(lines) == atmos.nlev_l + 2
        
        # Validate a data line has the expected number of columns
        data_line = lines[3]
        values = split(data_line, ",")
        @test length(values) == 14  # pressure + 13 flux columns
        
        # Clean up
        rm(csv_file, force=true)
    end
    
    # Test write_ncdf (NetCDF output)
    @testset "write_ncdf" begin
        nc_file = joinpath(OUT_DIR, "test_atmos.nc")
        
        # Remove file if it exists
        rm(nc_file, force=true)
        
        # Write NetCDF
        result = AGNI.save.write_ncdf(atmos, nc_file)
        
        # Check that function returns nothing
        @test result === nothing
        
        # Check that file was created
        @test isfile(nc_file)
        
        # Open and validate NetCDF structure
        Dataset(nc_file, "r") do ds
            # Check global attributes
            @test haskey(ds.attrib, "description")
            @test haskey(ds.attrib, "date")
            @test haskey(ds.attrib, "AGNI_version")
            @test haskey(ds.attrib, "SOCRATES_version")
            
            # Check dimensions
            @test haskey(ds.dim, "nlev_c")
            @test haskey(ds.dim, "nlev_l")
            @test haskey(ds.dim, "ngases")
            @test ds.dim["nlev_c"] == atmos.nlev_c
            @test ds.dim["nlev_l"] == atmos.nlev_l
            
            # Check scalar variables
            @test haskey(ds, "tmp_surf")
            @test haskey(ds, "flux_int")
            @test haskey(ds, "instellation")
            @test haskey(ds, "planet_radius")
            @test haskey(ds, "surf_gravity")
            
            # Check vector variables
            @test haskey(ds, "p")
            @test haskey(ds, "pl")
            @test haskey(ds, "tmp")
            @test haskey(ds, "tmpl")
            @test haskey(ds, "gases")
            @test haskey(ds, "x_gas")
            
            # Validate data values match atmosphere
            @test ds["tmp_surf"][] ≈ atmos.tmp_surf
            @test ds["planet_radius"][] ≈ atmos.rp
            @test ds["surf_gravity"][] ≈ atmos.grav_surf
            
            # Validate array dimensions
            @test length(ds["p"][:]) == atmos.nlev_c
            @test length(ds["pl"][:]) == atmos.nlev_l
            @test length(ds["tmp"][:]) == atmos.nlev_c
        end
        
        # Clean up
        rm(nc_file, force=true)
    end
    
    # Clean up atmosphere
    AGNI.atmosphere.deallocate!(atmos)
end


@testset "load" begin
    
    # Setup atmosphere and save it first
    p_surf = 1.0
    p_top = 1e-8
    theta = 65.0
    gravity = 10.0
    nlev_centre = 50  # Match the save test
    radius = 1.0e7
    tmp_surf = 500.0
    toa_heating = 1000.0
    mf_dict = Dict("H2O" => 0.8, "CO2" => 0.2)
    spfile = joinpath(RES_DIR, "spectral_files", "Dayspring", "48", "Dayspring.sf")
    
    # Create and save an atmosphere
    atmos1 = AGNI.atmosphere.Atmos_t()
    AGNI.atmosphere.setup!(atmos1, ROOT_DIR, OUT_DIR,
                          spfile,
                          toa_heating, 1.0, 0.0, theta,
                          tmp_surf,
                          gravity, radius,
                          nlev_centre, p_surf, p_top,
                          mf_dict, ""
                  )
    AGNI.atmosphere.allocate!(atmos1, "")
    AGNI.setpt.isothermal!(atmos1, tmp_surf)
    
    # Initialize flux arrays
    fill!(atmos1.flux_u_lw, 100.0)  # Use non-zero values to test
    fill!(atmos1.flux_d_lw, 200.0)
    fill!(atmos1.flux_n_lw, 300.0)
    fill!(atmos1.flux_u_sw, 10.0)
    fill!(atmos1.flux_d_sw, 20.0)
    fill!(atmos1.flux_n_sw, 30.0)
    fill!(atmos1.flux_u, 110.0)
    fill!(atmos1.flux_d, 220.0)
    fill!(atmos1.flux_n, 330.0)
    fill!(atmos1.flux_cdry, 50.0)
    fill!(atmos1.flux_l, 0.0)
    fill!(atmos1.flux_deep, 0.0)
    fill!(atmos1.flux_tot, 330.0)
    
    # Modify some values to make them unique
    atmos1.tmp_surf = 600.0
    atmos1.tmp_magma = 2000.0
    atmos1.skin_d = 5.0
    atmos1.skin_k = 3.0
    
    # Save to NetCDF
    nc_file = joinpath(OUT_DIR, "test_load.nc")
    rm(nc_file, force=true)
    AGNI.save.write_ncdf(atmos1, nc_file)
    
    # Test load_ncdf! (round-trip)
    @testset "load_ncdf_roundtrip" begin
        # Create a new atmosphere with same config
        atmos2 = AGNI.atmosphere.Atmos_t()
        AGNI.atmosphere.setup!(atmos2, ROOT_DIR, OUT_DIR,
                              spfile,
                              toa_heating, 1.0, 0.0, theta,
                              400.0,  # Different initial temp
                              gravity, radius,
                              nlev_centre, p_surf, p_top,
                              mf_dict, ""
                      )
        AGNI.atmosphere.allocate!(atmos2, "")
        
        # Load from NetCDF
        success = AGNI.load.load_ncdf!(atmos2, nc_file)
        
        # Check that load succeeded
        @test success == true
        
        # Validate that loaded values match original
        @test atmos2.tmp_surf ≈ atmos1.tmp_surf
        @test atmos2.tmp_magma ≈ atmos1.tmp_magma
        @test atmos2.skin_d ≈ atmos1.skin_d
        @test atmos2.skin_k ≈ atmos1.skin_k
        @test atmos2.rp ≈ atmos1.rp
        @test atmos2.grav_surf ≈ atmos1.grav_surf
        @test atmos2.instellation ≈ atmos1.instellation
        
        # Check array values
        @test all(atmos2.p .≈ atmos1.p)
        @test all(atmos2.pl .≈ atmos1.pl)
        @test all(atmos2.tmp .≈ atmos1.tmp)
        @test all(atmos2.tmpl .≈ atmos1.tmpl)
        
        # Check gas VMRs
        for gas in atmos1.gas_names
            @test all(atmos2.gas_vmr[gas] .≈ atmos1.gas_vmr[gas])
        end
        
        # Check flux values
        @test all(atmos2.flux_u .≈ atmos1.flux_u)
        @test all(atmos2.flux_d .≈ atmos1.flux_d)
        @test all(atmos2.flux_n .≈ atmos1.flux_n)
        @test all(atmos2.flux_cdry .≈ atmos1.flux_cdry)
        
        AGNI.atmosphere.deallocate!(atmos2)
    end
    
    # Test load_ncdf! error: mismatched levels
    @testset "load_ncdf_mismatch_levels" begin
        # Create atmosphere with different number of levels
        atmos3 = AGNI.atmosphere.Atmos_t()
        AGNI.atmosphere.setup!(atmos3, ROOT_DIR, OUT_DIR,
                              spfile,
                              toa_heating, 1.0, 0.0, theta,
                              tmp_surf,
                              gravity, radius,
                              60,  # Different nlev_centre
                              p_surf, p_top,
                              mf_dict, ""
                      )
        AGNI.atmosphere.allocate!(atmos3, "")
        
        # Store original values to ensure they weren't overwritten
        orig_tmp_surf = atmos3.tmp_surf
        
        # Try to load - the function will log an error but may still return
        # The important thing is that the data wasn't corrupted
        AGNI.load.load_ncdf!(atmos3, nc_file)
        
        # Verify that original atmosphere wasn't modified (indicates load failed properly)
        # If load succeeded incorrectly, tmp_surf would change to 600.0
        @test atmos3.tmp_surf == orig_tmp_surf
        
        AGNI.atmosphere.deallocate!(atmos3)
    end
    
    # Test load_ncdf! error: mismatched gases
    @testset "load_ncdf_mismatch_gases" begin
        # Create atmosphere with different number of gases
        atmos4 = AGNI.atmosphere.Atmos_t()
        AGNI.atmosphere.setup!(atmos4, ROOT_DIR, OUT_DIR,
                              spfile,
                              toa_heating, 1.0, 0.0, theta,
                              tmp_surf,
                              gravity, radius,
                              nlev_centre, p_surf, p_top,
                              Dict("H2O" => 0.5, "CO2" => 0.3, "N2" => 0.2),  # Different gases
                              ""
                      )
        AGNI.atmosphere.allocate!(atmos4, "")
        
        # Store original values
        orig_tmp_surf = atmos4.tmp_surf
        
        # Try to load - should not modify atmosphere
        AGNI.load.load_ncdf!(atmos4, nc_file)
        
        # Verify original atmosphere wasn't modified
        @test atmos4.tmp_surf == orig_tmp_surf
        
        AGNI.atmosphere.deallocate!(atmos4)
    end
    
    # Clean up
    AGNI.atmosphere.deallocate!(atmos1)
    rm(nc_file, force=true)
end
