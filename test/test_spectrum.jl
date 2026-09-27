"""
Tests for `src/energy/spectrum.jl`.

Invariants exercised:
- Spectral file parsing and block insertion succeed with valid inputs.
- Aerosol headers and aerosol averaging produce non-empty outputs.
- Guard paths return gracefully on invalid inputs.
"""
const _SPECTRUM_TESTS_DOC = nothing
using Test
using AGNI
using Printf

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR = joinpath(ROOT_DIR,"res/")
RAD_DIR = AGNI.paths.RAD_DIR

temp_sf = tempname() * ".sf"

# Use 4 bytes as test of single precision
module DummySocratesPrecision
    const SOCRATES_REAL_BYTES = 4
end

# Missing precision simulates older SOCRATES Julia builds.
module DummySocratesNoPrecision
end


@testset "spectrum" begin

    # Check SOCRATES metainfo
    @testset "socrates_meta" begin

        # Check socrates was found
        @test isfile(joinpath(RAD_DIR,"version"))

        # Check SOCRATES precision getter
        @testset "socrates_precision" begin

            # test defined as single
            precision = AGNI.spectrum.get_socrates_precision(DummySocratesPrecision)
            @test !isempty(precision)
            @test precision == "single"

            # test fallback to double
            fallback = AGNI.spectrum.get_socrates_precision(DummySocratesNoPrecision)
            @test !isempty(fallback)
            @test fallback == "double"
        end

        # Test get_socrates_version
        @testset "socrates_version" begin
            version = AGNI.spectrum.get_socrates_version(RAD_DIR)
            @test version isa String
            @test startswith(version, "2") # this millenium
        end
    end

    # Test count_gases with a real spectral file
    @testset "count_gases_valid" begin
        # Use one of the test spectral files
        spfile = joinpath(RES_DIR, "spectral_files", "Dayspring", "48", "Dayspring.sf")
        @test isfile(spfile)
        num_gases = AGNI.spectrum.count_gases(spfile)
        @test num_gases > 0
        @test num_gases isa Int
    end

    # Test count_gases with non-existent file
    @testset "count_gases_missing" begin
        rm(temp_sf, force=true)  # Ensure the file does not exist
        num_gases = AGNI.spectrum.count_gases(temp_sf)
        @test num_gases == -1
    end

    # Test count_gases with invalid file (no gas count line)
    @testset "count_gases_invalid" begin
        # Create a temporary file with invalid content
        write(temp_sf, "This is not a valid spectral file\nNo gas information here\n")
        num_gases = AGNI.spectrum.count_gases(temp_sf)
        @test num_gases == -1
        rm(temp_sf, force=true)
    end

    # Test count_gases with malformed gas count line (wrong format)
    @testset "count_gases_malformed_line" begin
        write(temp_sf, "Total number of gaseous absorbers not_a_number\n")
        num_gases = AGNI.spectrum.count_gases(temp_sf)
        @test num_gases == -1
        rm(temp_sf, force=true)
    end

    # Test count_gases with wrong split count
    @testset "count_gases_wrong_split" begin
        write(temp_sf, "Total number of gaseous absorbers = 5 = extra\n")
        num_gases = AGNI.spectrum.count_gases(temp_sf)
        @test num_gases == -1
        rm(temp_sf, force=true)
    end

    # Test count_gases with unparseable number (to hit catch block)
    @testset "count_gases_unparseable" begin
        write(temp_sf, "Total number of gaseous absorbers = not_a_number\n")
        num_gases = AGNI.spectrum.count_gases(temp_sf)
        @test num_gases == -1
        rm(temp_sf, force=true)
    end

    # Test insert_aerosol_header
    @testset "insert_aerosol_header" begin
        # Create a minimal spectral file
        sf_content = """
Line 1
Line 2
Line 3
*END
Some more content
*BLOCK: TYPE =    1
"""
        write(temp_sf, sf_content)

        # Try inserting aerosol header with a valid aerosol name
        success = AGNI.spectrum.insert_aerosol_header(temp_sf, ["dust"])
        @test success

        # Check that file was modified
        content = read(temp_sf, String)
        @test contains(content, "Total number of aerosols")
        @test contains(content, "Dust-like Aerosol")

        rm(temp_sf, force=true)
    end

    # Test insert_aerosol_header with existing aerosol data
    @testset "insert_aerosol_header_existing" begin
        write(temp_sf, "This file already has aerosols in it\n*END\n")
        success = AGNI.spectrum.insert_aerosol_header(temp_sf, ["dust"])
        @test !success
        rm(temp_sf, force=true)
    end

    # Test insert_aerosol_header with missing *END marker
    @testset "insert_aerosol_header_no_end" begin
        write(temp_sf, "Line 1\nLine 2\nNo END marker here\n")
        success = AGNI.spectrum.insert_aerosol_header(temp_sf, ["sulphuric_acid"])
        @test !success
        rm(temp_sf, force=true)
    end

    # Test blackbody_star
    @testset "blackbody_star" begin
        Teff = 5800.0  # Sun-like star
        S0 = 1361.0     # Solar constant

        wl, fl = AGNI.spectrum.blackbody_star(Teff, S0)

        @test length(wl) == length(fl)
        @test length(wl) > 0
        @test all(fl .>= AGNI.spectrum.SMALLFLOAT)
        @test all(fl .<= AGNI.spectrum.BIGFLOAT)
        @test wl[1] == 1.0
        @test wl[end] == 100e3
    end

    # Test load_from_file with valid file
    @testset "load_from_file_valid" begin
        # Create a test stellar spectrum file
        test_data = """
# Header line 1
# Header line 2
100.0 1.5e10
200.0 2.3e10
300.0 1.8e10
"""
        write(temp_sf, test_data)
        wl, fl = AGNI.spectrum.load_from_file(temp_sf)

        @test length(wl) == 3
        @test length(fl) == 3
        @test wl[1] == 100.0
        @test fl[1] == 1.5e10
        rm(temp_sf, force=true)
    end

    # Test write_to_socrates_format
    @testset "write_to_socrates_format" begin
        # Create test wavelength and flux arrays
        wl = collect(range(100.0, 1000.0, length=1000))
        fl = ones(Float64, 1000) .* 1e10

        temp_star = tempname() * ".txt"
        success = AGNI.spectrum.write_to_socrates_format(wl, fl, temp_star, 500)
        @test success
        @test isfile(temp_star)

        # Check file content
        content = read(temp_star, String)
        @test contains(content, "Star spectrum at TOA")
        @test contains(content, "*BEGIN_DATA")
        @test contains(content, "*END")

        rm(temp_star, force=true)
    end

    # Test write_to_socrates_format with mismatched arrays
    @testset "write_to_socrates_format_mismatch" begin
        wl = collect(range(100.0, 1000.0, length=1000))
        fl = ones(Float64, 999)  # Different length

        temp_star = tempname() * ".txt"
        success = AGNI.spectrum.write_to_socrates_format(wl, fl, temp_star)
        @test !success
        rm(temp_star, force=true)
    end

    # Test write_to_socrates_format with short spectrum
    @testset "write_to_socrates_format_short" begin
        wl = collect(range(100.0, 1000.0, length=100))
        fl = ones(Float64, 100) .* 1e10

        temp_star = tempname() * ".txt"
        success = AGNI.spectrum.write_to_socrates_format(wl, fl, temp_star)
        @test success
        rm(temp_star, force=true)
    end

    # Test write_to_socrates_format with invalid wavelength (too small)
    @testset "write_to_socrates_format_wl_too_small" begin
        wl = [1e-50, 100.0, 200.0]
        fl = ones(Float64, 3) .* 1e10

        temp_star = tempname() * ".txt"
        success = AGNI.spectrum.write_to_socrates_format(wl, fl, temp_star)
        @test !success
        rm(temp_star, force=true)
    end

    # Test write_to_socrates_format with descending wavelength array
    @testset "write_to_socrates_format_descending" begin
        wl = [1000.0, 500.0, 100.0]  # Descending order
        fl = ones(Float64, 3) .* 1e10

        temp_star = tempname() * ".txt"
        success = AGNI.spectrum.write_to_socrates_format(wl, fl, temp_star)
        @test !success
        rm(temp_star, force=true)
    end

    # Test write_to_socrates_format with duplicate wavelengths
    @testset "write_to_socrates_format_duplicates" begin
        wl = [100.0, 200.0, 200.0, 300.0]  # Has duplicate
        fl = ones(Float64, 4) .* 1e10

        temp_star = tempname() * ".txt"
        success = AGNI.spectrum.write_to_socrates_format(wl, fl, temp_star)
        @test success  # Should succeed after removing duplicates
        rm(temp_star, force=true)
    end

    # Clean up any remaining temp files
    rm(temp_sf, force=true)


    @testset "aerosol_guard_invalid_input" begin
        tmpdir = mktempdir()
        orig = joinpath(tmpdir, "base.sf")
        star = joinpath(tmpdir, "star.dat")
        outp = joinpath(tmpdir, "out.sf")

        # insert_blocks: missing original file
        @test !AGNI.spectrum.insert_blocks( RAD_DIR,
            orig, star, outp, false, false; aerosol_avg_files=Dict{String,String}()
        )

        # write minimal original + _k to pass cp step, but missing star file should fail
        write(orig, "Total number of gaseous absorbers = 1\n*END\n")
        write(orig * "_k", "dummy k-table\n")
        @test !AGNI.spectrum.insert_blocks( RAD_DIR,
            orig, star, outp, false, false; aerosol_avg_files=Dict{String,String}()
        )

        # write star file, but malformed spectral file => gas count parse failure
        write(star, "header\nheader\n1.0 1.0\n")
        write(orig, "No gas count line here\n*END\n")
        @test !AGNI.spectrum.insert_blocks( RAD_DIR,
            orig, star, outp, false, false; aerosol_avg_files=Dict{String,String}()
        )

        # valid gas-count line, but missing prep binaries / execution path should still be handled
        write(orig, "Total number of gaseous absorbers = 1\n*END\n")
        @test !AGNI.spectrum.insert_blocks( RAD_DIR,
            orig, star, outp, true, false; aerosol_avg_files=Dict{String,String}()
        )

        # generate_aerosol_avg_files guard paths
        missing_orig = joinpath(tmpdir, "missing.sf")
        @test isempty(AGNI.spectrum.generate_aerosol_avg_files( RAD_DIR,
            missing_orig, ["dust"], tmpdir, 2, star, tmpdir
        ))

        @test isempty(AGNI.spectrum.generate_aerosol_avg_files(RAD_DIR,
            orig, ["dust"], tmpdir, 0, star, tmpdir
        ))

        @test isempty(AGNI.spectrum.generate_aerosol_avg_files(RAD_DIR,
            orig, ["dust"], tmpdir, 2, "", tmpdir
        ))

        @test isempty(AGNI.spectrum.generate_aerosol_avg_files(RAD_DIR,
            orig, ["dust"], tmpdir, 2, joinpath(tmpdir, "missing_star.dat"), tmpdir
        ))

        # no species => early empty return without external tool
        @test isempty(AGNI.spectrum.generate_aerosol_avg_files(RAD_DIR,
            orig, String[], tmpdir, 2, star, tmpdir
        ))

        # species provided but missing .mon file
        @test isempty(AGNI.spectrum.generate_aerosol_avg_files(RAD_DIR,
            orig, ["dust"], tmpdir, 2, star, tmpdir
        ))

        rm(tmpdir; force=true, recursive=true)
    end

    # Check successful data insersion
    @testset "aerosol_insertion_success" begin

        # Setup paths
        tmpdir = mktempdir()
        spfile = joinpath(RES_DIR, "spectral_files", "Dayspring", "48", "Dayspring.sf")

        # A real spectral file
        @test isfile(spfile)
        @test isfile(spfile * "_k")

        # Some SOCRATES scattering data is pre-computed
        scattering_dir = joinpath(RES_DIR, "scattering")
        @test isdir(scattering_dir)

        # Choose the 1,2 aerosol species with a matching .mon file in the repo data.
        available_species = [
            s for s in AGNI.spectrum.input_head_pcf.aerosol_suffix
            if isfile(joinpath(scattering_dir, s * ".mon"))
        ]
        @test !isempty(available_species)
        species = [available_species[1], available_species[2]] # first two

        # Get stellar spectrum
        star_file = joinpath(tmpdir, "star.dat")
        wl = collect(range(100.0, 1000.0, length=1000))
        fl = ones(Float64, 1000) .* 1e10
        @test AGNI.spectrum.write_to_socrates_format(wl, fl, star_file, 500)
        @test isfile(star_file)

        # Generate aerosol average-properties file(s)
        avg_files = AGNI.spectrum.generate_aerosol_avg_files(
            RAD_DIR, spfile, species, tmpdir, 1, star_file, scattering_dir
        )

        # should be 2
        @test length(avg_files) == length(species)

        # check output from input
        @test haskey(avg_files, species[1])
        avg_file = avg_files[species[1]]
        @test isfile(avg_file)
        @test filesize(avg_file) > 0
        @test length(readlines(avg_file)) > 1

        # insert data into spectral file
        outp_file = joinpath(tmpdir, "runtime.sf")
        success = AGNI.spectrum.insert_blocks(
            RAD_DIR, spfile, star_file, outp_file, false, true;
            aerosol_avg_files=avg_files
        )

        # check the spectral file
        @test success
        @test isfile(outp_file)
        @test isfile(outp_file * "_k")
        @test filesize(outp_file) > 0
        @test filesize(outp_file * "_k") > 0

        # check the spectral file contains aerosol data
        out_lines = readlines(outp_file)
        @test any(contains.(out_lines, "Total number of aerosols"))
        @test any(contains.(out_lines, "List of indexing numbers of aerosols"))

        # remove this folder
        rm(tmpdir; force=true, recursive=true)
    end


    # -------------
    # Block-0 aerosol rows are read by SOCRATES with the Fortran format (i5, 7x, i5, 7x, a),
    # so the index must occupy columns 1-5 and the type number columns 13-17. Type numbers
    # above 99 (used for runtime-calculated aerosols) must fit this layout.
    # -------------
    @testset "aerosol_row_matches_fortran_layout" begin
        for (idx, typ) in ((1, 4), (3, 101), (12, 32), (99999, 99999))
            row = AGNI.spectrum.format_aerosol_row(idx, typ, "name")
            @test parse(Int, row[1:5]) == idx
            @test parse(Int, row[13:17]) == typ
            @test strip(row[18:end]) == "name"
        end
        # Discrimination guard: the previous layout ("   %2d          %2d") put 3-digit
        # type numbers outside columns 13-17
        old = "    1          101       name "
        @test tryparse(Int, strip(old[13:17])) != 101
        # Error contract: out-of-range numbers
        @test_throws ErrorException AGNI.spectrum.format_aerosol_row(0, 101, "x")
        @test_throws ErrorException AGNI.spectrum.format_aerosol_row(1, 100000, "x")
    end

    # -------------
    # Band edges parsed from block 1 of a real spectral file are positive, contiguous, and
    # increasing; the first band starts at 0.2857 μm in Dayspring/48.
    # -------------
    @testset "band_edges_parsed_from_block_1" begin
        spfile = joinpath(RES_DIR, "spectral_files", "Dayspring", "48", "Dayspring.sf")
        bands = AGNI.spectrum.read_band_edges(spfile)
        @test size(bands) == (48, 2)
        @test all(bands .> 0.0)
        @test all(bands[:,2] .> bands[:,1])
        @test all(isapprox.(bands[2:end,1], bands[1:end-1,2]; rtol=1e-8))
        @test isapprox(bands[1,1], 2.857143673e-7; rtol=1e-8)
        # Error contract: missing file, and a file without block 1
        @test_throws ErrorException AGNI.spectrum.read_band_edges(tempname())
        write(temp_sf, "*BLOCK: TYPE =    0\n*END\n")
        @test_throws ErrorException AGNI.spectrum.read_band_edges(temp_sf)
        rm(temp_sf, force=true)
    end

    # -------------
    # Runtime-calculated aerosols appended to a spectral file are read back by SOCRATES with
    # the written values, in the dry-aerosol parametrisation, alongside any aerosols already
    # present. Invalid properties are rejected without modifying the file.
    # -------------
    @testset "custom_aerosols_round_trip_through_socrates" begin
        SOC = AGNI.atmosphere.SOCRATES
        tmpdir = mktempdir()
        spfile = joinpath(RES_DIR, "spectral_files", "Dayspring", "48", "Dayspring.sf")
        work = joinpath(tmpdir, "work.sf")
        cp(spfile, work; force=true)
        cp(spfile*"_k", work*"_k"; force=true)
        nb = 48

        # Values chosen to vary across bands, so that a band offset would be detected
        ka = [[10.0*b for b in 1:nb], fill(5.0, nb)]
        ks = [[1000.0 + b for b in 1:nb], fill(0.0, nb)]     # purely absorbing second aerosol
        g  = [[0.5 + 0.004*b for b in 1:nb], fill(-0.2, nb)]

        # Guard paths: mismatched lengths, unphysical values, NaN, wrong band count
        orig = read(work, String)
        @test !AGNI.spectrum.append_custom_aerosols!(work, ["a"], [101, 102], ka[1:1], ks[1:1], g[1:1])
        @test !AGNI.spectrum.append_custom_aerosols!(work, ["a"], [101], [-ka[1]], ks[1:1], g[1:1])
        @test !AGNI.spectrum.append_custom_aerosols!(work, ["a"], [101], ka[1:1], ks[1:1], [fill(1.5, nb)])
        @test !AGNI.spectrum.append_custom_aerosols!(work, ["a"], [101], [fill(NaN, nb)], ks[1:1], g[1:1])
        @test !AGNI.spectrum.append_custom_aerosols!(work, ["a"], [101], [ones(nb-1)], [ones(nb-1)], [zeros(nb-1)])
        @test read(work, String) == orig
        @test AGNI.spectrum.append_custom_aerosols!(work, String[], Int[], Vector{Float64}[],
                                                    Vector{Float64}[], Vector{Float64}[])

        # Append two aerosols to a file with no aerosols
        @test AGNI.spectrum.append_custom_aerosols!(work, ["sio2", "feo"], [101, 102], ka, ks, g)

        sp = SOC.StrSpecData()
        kw = Dict{Symbol,Any}(:spectrum=>sp, :spectral_file=>work)
        if startswith(AGNI.spectrum.get_socrates_version(RAD_DIR), "24")
            kw[:l_all_gasses] = true
        else
            kw[:l_all_gases] = true
        end
        SOC.set_spectrum(; kw...)
        A = sp.Aerosol
        @test A.n_aerosol == 2
        @test Int.(A.type_aerosol[1:2]) == [101, 102]
        @test all(Int.(A.i_aerosol_parm[1:2]) .== SOC.rad_pcf.ip_aerosol_param_dry)
        @test Bool(sp.Basic.l_present[11])
        for i in 1:2, b in 1:nb
            @test isapprox(A.abs[1, i, b],        ka[i][b]; rtol=1e-8)
            @test isapprox(A.scat[1, i, b],       ks[i][b]; rtol=1e-8, atol=1e-12)
            @test isapprox(A.phf_fnc[1, 1, i, b], g[i][b];  rtol=1e-8)
        end
        # Discrimination guard: a one-band offset would be detected
        @test abs(A.abs[1, 1, 10] - ka[1][11]) > 1.0

        rm(tmpdir; force=true, recursive=true)
    end

    # -------------
    # Runtime-calculated aerosols are appended after aerosols inserted by prep_spec, taking
    # the next indices in the block-0 list, and both are read back by SOCRATES.
    # -------------
    @testset "custom_aerosols_follow_prep_spec_aerosols" begin
        SOC = AGNI.atmosphere.SOCRATES
        tmpdir = mktempdir()
        spfile = joinpath(RES_DIR, "spectral_files", "Dayspring", "48", "Dayspring.sf")
        star_file = joinpath(tmpdir, "star.dat")
        wl = collect(range(100.0, 1000.0, length=1000))
        @test AGNI.spectrum.write_to_socrates_format(wl, ones(1000) .* 1e10, star_file, 500)
        avg = AGNI.spectrum.generate_aerosol_avg_files(RAD_DIR, spfile, ["soot"], tmpdir, 1,
                                                        star_file, joinpath(RES_DIR, "scattering"))
        outp = joinpath(tmpdir, "runtime.sf")
        @test AGNI.spectrum.insert_blocks(RAD_DIR, spfile, star_file, outp, false, true;
                                            aerosol_avg_files=avg)
        nb = 48
        @test AGNI.spectrum.append_custom_aerosols!(outp, ["sio2"], [101],
                                                    [fill(2.0, nb)], [fill(3.0, nb)], [fill(0.1, nb)])
        lines = readlines(outp)
        @test any(l -> startswith(l, "Total number of aerosols =     2"), lines)

        sp = SOC.StrSpecData()
        kw = Dict{Symbol,Any}(:spectrum=>sp, :spectral_file=>outp)
        if startswith(AGNI.spectrum.get_socrates_version(RAD_DIR), "24")
            kw[:l_all_gasses] = true
        else
            kw[:l_all_gases] = true
        end
        SOC.set_spectrum(; kw...)
        A = sp.Aerosol
        @test A.n_aerosol == 2
        @test Int.(A.type_aerosol[1:2]) == [4, 101]    # soot, then custom
        @test isapprox(A.abs[1, 2, 7], 2.0; rtol=1e-8)
        @test A.abs[1, 1, 7] > 0.0                     # soot data from prep_spec intact
        @test abs(A.abs[1, 1, 7] - 2.0) > 1e-3

        rm(tmpdir; force=true, recursive=true)
    end


    # -------------
    # Independent check of AGNI's Mie code against SOCRATES' own Mie implementation
    # (Cscatter -M), for soot with a log-normal distribution (r_g = 0.5 μm, σ_g = 1.65).
    # Wavelengths are taken from the SOCRATES refractive index table so that interpolation
    # does not enter. Compared quantity is ⟨σ⟩/⟨V⟩ (independent of density). Agreement is
    # limited by size-distribution quadrature and 6-digit output (observed ≈ 1e-3).
    # -------------
    @testset "mie_agrees_with_socrates_cscatter" begin

        setenv_file = joinpath(RAD_DIR, "set_rad_env")
        refract = joinpath(RAD_DIR, "data", "aerosol", "refract_soot")

        tmpdir = mktempdir()
        rows = Vector{Float64}[]
        on = false
        for l in readlines(refract)
            startswith(l, "*BEGIN_DATA") && (on = true; continue)
            startswith(l, "*END") && break
            on && push!(rows, parse.(Float64, split(l)))
        end
        sel = [r for r in rows if 0.25e-6 <= r[1] <= 30e-6]
        open(joinpath(tmpdir, "wl"), "w") do f
            write(f, "Wavelengths\n*BEGIN_DATA\n")
            foreach(r -> write(f, @sprintf("  %.6e\n", r[1])), sel)
            write(f, "*END\n")
        end
        cmd = "cd $tmpdir && source $setenv_file >/dev/null 2>&1 && " *
                "$(joinpath(RAD_DIR, "sbin", "Cscatter")) -w wl -r $refract -l -t 1 -C 4 " *
                "-g 1.0 0.5e-6 1.65 -n 1.0e8 -M -o soot.mon >/dev/null 2>&1"
        run(`bash -c $cmd`; wait=true)
        mon = joinpath(tmpdir, "soot.mon")
        @test isfile(mon)

        L = readlines(mon)
        φ = parse(Float64, split(L[findfirst(l->contains(l, "Volume fraction"), L)], "=")[2])
        reff_soc = parse(Float64, split(split(L[findfirst(l->contains(l, "Effective radius"), L)], "=")[2])[1])
        i0 = findfirst(l->contains(l, "Wavelength (m)"), L)
        dat = [parse.(Float64, split(l)) for l in L[i0+1:end] if length(split(l)) == 4]
        @test length(dat) == length(sel)

        # SOCRATES reports r_eff = r_g exp(2.5 ln²σ_g), matching AGNI's conversion. The
        # tolerance reflects SOCRATES' own size quadrature, which recovers the input
        # number density to only ~1e-4.
        r_eff = 0.5e-6 * exp(2.5 * log(1.65)^2)
        @test isapprox(reff_soc, r_eff; rtol=1e-4)
        @test isapprox(AGNI.mie.geometric_radius(r_eff, 1.65), 0.5e-6; rtol=1e-12)

        λ = [r[1] for r in sel]
        m = [complex(r[2], r[3]) for r in sel]
        ka, ks, g = AGNI.mie.mass_coefficients(λ, m, r_eff, 1.65, 1.0)
        for i in eachindex(λ)
            @test isapprox(ka[i], dat[i][2] / φ; rtol=5e-3)
            @test isapprox(ks[i], dat[i][3] / φ; rtol=5e-3)
            @test isapprox(g[i],  dat[i][4];     atol=2e-3)
        end
        # Discrimination guard: using r_g in place of r_eff changes k_sca by >> 0.5%
        _, ks_wrong, _ = AGNI.mie.mass_coefficients(λ[1:1], m[1:1], 0.5e-6, 1.65, 1.0)
        @test abs(ks_wrong[1] / (dat[1][3] / φ) - 1) > 0.05

        rm(tmpdir; force=true, recursive=true)
    end

    # -------------
    # Independent check of AGNI's stellar-weighted 'thin' band averaging against SOCRATES'
    # scatter_average_90, applied to the same monochromatic soot data (soot.mon). Bands
    # beyond 40 μm are excluded: soot.mon has no data between 40 μm and 10 mm, so results
    # there depend only on how the gap is interpolated. Agreement elsewhere is limited by
    # the sparse .mon sampling (62 points) and the down-binned stellar spectrum. This checks
    # units, normalisation, and band mapping; in these narrow bands the choice of weighting
    # changes results by less than the tolerance, so the weighting itself is tested with
    # synthetic data in test_aerosol_optics.jl.
    # -------------
    @testset "band_average_agrees_with_socrates_scatter_average" begin
        tmpdir = mktempdir()
        spfile = joinpath(RES_DIR, "spectral_files", "Dayspring", "48", "Dayspring.sf")
        monfile = joinpath(RES_DIR, "scattering", "soot.mon")
        wl, fl = AGNI.spectrum.load_from_file(joinpath(RES_DIR, "stellar_spectra", "sun.txt"))
        star_file = joinpath(tmpdir, "star.dat")
        @test AGNI.spectrum.write_to_socrates_format(wl, fl, star_file)
        avg = AGNI.spectrum.generate_aerosol_avg_files(RAD_DIR, spfile, ["soot"], tmpdir, 1,
                                                        star_file, joinpath(RES_DIR, "scattering"))
        @test isfile(avg["soot"])

        function readblock(path)
            L = readlines(path)
            φ = parse(Float64, split(split(L[findfirst(l->contains(l, "Volume fraction"), L)], "=")[2])[1])
            i0 = findfirst(l->contains(l, "Absorption"), L)
            d = [parse.(Float64, split(l)) for l in L[i0+1:end]
                    if length(split(l)) == 4 && !isnothing(tryparse(Float64, split(l)[1]))]
            return φ, permutedims(hcat(d...))
        end
        φm, mon = readblock(monfile)
        φa, soc = readblock(avg["soot"])

        λm = mon[:,1]
        lin(y, x) = x <= λm[1] ? y[1] : x >= λm[end] ? y[end] :
                    (i = searchsortedlast(λm, x); y[i] + (x - λm[i]) / (λm[i+1] - λm[i]) * (y[i+1] - y[i]))
        props(λ) = ([lin(mon[:,2] ./ φm, x) for x in λ], [lin(mon[:,3] ./ φm, x) for x in λ],
                    [lin(mon[:,4], x) for x in λ], fill(false, length(λ)))
        bands = AGNI.spectrum.read_band_edges(spfile)
        star = AGNI.aerosol_optics.star_cumulative(wl, fl)
        ka, ks, g, _ = AGNI.aerosol_optics.band_average(bands, star, λm, props)

        covered = findall(bands[:,2] .<= 40e-6)
        @test length(covered) >= 40
        for b in covered
            @test isapprox(ka[b], soc[b,2] / φa; rtol=0.05)
            @test isapprox(ks[b], soc[b,3] / φa; rtol=0.05)
            @test isapprox(g[b],  soc[b,4];      atol=0.01)
        end
        rm(tmpdir; force=true, recursive=true)
    end

end
