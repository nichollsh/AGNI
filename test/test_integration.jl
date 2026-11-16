
using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR         = joinpath(ROOT_DIR,"res/")
OUT_DIR         = joinpath(ROOT_DIR,"out/")
TEST_DIR        = joinpath(ROOT_DIR,"test/")
CONFIG_PATH     = joinpath(TEST_DIR,"test.toml")


@testset "integration" begin

    # -------------
    # Run example TOML config file (writes to output folder)
    # -------------
    @testset "full_run" begin

        # Test opening test config
        cfg = AGNI.open_config(CONFIG_PATH)
        @test haskey(cfg, "planet")
        @test haskey(cfg["planet"], "tmp_surf")
        @test haskey(cfg, "execution")

        # set to use spectral surface albedo
        cfg["planet"]["surface_material"] = "res/surface_albedos/lunarmarebasalt.dat"

        # check return code is fine
        succ = AGNI.main(cfg_path=CONFIG_PATH)
        @test succ

        # check io directory exists
        @test isdir(cfg["files"]["io_dir"])
        @test isfile(joinpath(OUT_DIR, "plot_contfunc1.png"))
        @test isfile(joinpath(OUT_DIR, "plot_contfunc2.png"))
        @test isfile(joinpath(OUT_DIR, "plot_tau.png"))

        # with [plots] at_runtime=true and animate=true, solve_energy!'s
        # in-loop plot_step() (plot_pt/plot_fluxes/plot_vmr/plot_radius/
        # plot_cloud/jacobian/combined) must have run at least once per
        # solver step, saving one frame per step plus a combined animation -
        # this is real solver-loop plotting coverage, not just the final
        # post-solve plots checked above
        frames_dir = joinpath(cfg["files"]["io_dir"], "frames")
        @test isdir(frames_dir)
        frame_files = filter(f -> endswith(f, ".png"), readdir(frames_dir))
        @test length(frame_files) > 1
        @test isfile(joinpath(frames_dir, "0001.png"))
        @test isfile(joinpath(OUT_DIR, "animation.mp4"))

        # on a successful solve, the transient solver.png (the source frames
        # are copied from) is cleaned up rather than left behind
        @test !isfile(joinpath(OUT_DIR, "solver.png"))

        # -------------
        # Compare result from NetCDF
        # -------------

        # Variable to store expected values
        arr_e::Array{Float64,1}     = zeros(Float64, cfg["execution"]["num_levels"]+1)
        arr_o_tmp::Array{Float64,1} = zero(arr_e)
        arr_o_prs::Array{Float64,1} = zero(arr_e)
        arr_o_flN::Array{Float64,1} = zero(arr_e)
        arr_o_flC::Array{Float64,1} = zero(arr_e)
        arr_o_Kzz::Array{Float64,1} = zero(arr_e)

        atol = Float64(cfg["execution"]["converge_atol"])
        rtol = 1e-3

        # Read result
        @debug "ALL DEBUG SUPPRESSED"
        with_logger(MinLevelLogger(current_logger(), Logging.Info-200)) do

            # open netcdf
            ds = Dataset(joinpath(OUT_DIR,"atm.nc"),"r")

            # read profiles from netCDF file
            arr_o_tmp[:] .= ds["tmpl"][:]
            arr_o_prs[:] .= ds["pl"][:]
            arr_o_flN[:] .= ds["fl_N"][:]
            arr_o_flC[:] .= ds["fl_cnvct"][:]
            arr_o_Kzz[:] .= ds["Kzz"][:]

            # close netcdf
            close(ds)
        end
        @debug "ALL DEBUG RESTORED"

        # pressure profile
        @debug ("Pressure...")
        arr_e[:] .= [0.001, 0.00180007897342944, 0.00324028431058278, 0.00583276765541337,
                    0.0104994424134089, 0.0188998255211106, 0.0340211785220364,
                    0.0612408081088068, 0.11023829099249, 0.198437629682378,
                    0.357203404728425, 0.642994338089044, 1.15744058802827, 2.08348446550349,
                    3.7504365778197, 6.7510820249139, 12.152480800945, 21.875425164786,
                    39.3774928739604, 70.8825969487837, 127.594272349579, 229.679766786507,
                    413.441718814568, 744.227744776629, 1339.66871481522, 2411.50948490012,
                    4340.90751799436, 7813.97634874342, 14065.774524248, 25319.5049660982,
                    45577.1085071156, 82042.3946933709, 147682.789617336, 265840.684327569,
                    478534.226140151, 861399.398541213, 1550586.9450388, 2791178.95623854,
                    5024342.55020371, 9044213.37992856, 16280298.3364186, 29305822.7165454,
                    52752795.2711041, 94959197.5571423, 170934054.856344, 307694797.989939,
                    553874936.095308, 698575789.445476, 700000000.0]
        test_check = all(abs.(arr_e[:] .- arr_o_prs[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        if !test_check
            @error ("Fail 'prs' \n Observed: $(arr_o_prs) \n Expected: $(arr_e)")
        end
        @test test_check

        # temperature profile
        @debug ("Temperature...")
        arr_e[:] .= [155.317074485401, 173.62957588613, 181.405471247682, 169.92952555296,
                    173.257122481789, 183.040926359106, 180.670654850232, 180.622656442572,
                    177.430681118349, 173.284536547052, 170.995843392193, 172.198252234912,
                    173.424491029447, 173.256360960756, 170.755697193101, 170.343736762165,
                    172.154396502088, 172.47456312956, 170.422743860886, 169.203225432218,
                    171.467818949909, 174.647625064909, 177.977206747958, 183.736411066939,
                    194.751202028798, 207.534596302865, 230.035158537662, 260.204111692055,
                    289.587987044993, 326.148744617617, 370.496752102442, 418.781460349971,
                    466.026849534799, 505.802821057431, 538.441570484272, 564.273678704435,
                    582.165642674503, 595.324986967272, 609.343595544882, 626.615564940366,
                    644.506263405301, 659.430444021865, 671.303149701877, 681.788496461628,
                    691.383708201906, 700.941316835694, 714.090321341753, 723.296955294404,
                    723.402090893429 ]
        test_check = all(abs.(arr_e[:] .- arr_o_tmp[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        if !test_check
            @error ("Fail 'tmp' \n Observed: $(arr_o_tmp) \n Expected: $(arr_e)")
        end
        @test test_check

        # radiative flux
        @debug ("Radiative flux...")
        arr_e[:] .= [-0.0148119407702438, -0.0148755271200685, -0.0148007292347643,
                    -0.0147311189086281, -0.0149350325094133, -0.0147822744714858,
                    -0.0148217780746336, -0.0148009499044974, -0.0147964715718558,
                    -0.0147933151854431, -0.014811401774466, -0.0148213040046699,
                    -0.0148115336518231, -0.0148100209825657, -0.0147926515066388,
                    -0.014826790453867, -0.0148109780162713, -0.0148158275396781,
                    -0.0147928679559186, -0.01482300551055, -0.014820552814399,
                    -0.0148285820117735, -0.0148123597745098, -0.0148274566036548,
                    -0.0147569803602892, -0.0146775524008262, -0.538614347956695,
                    -33.9579611882832, -65.1951073063441, -59.018497249933,
                    -44.2590713193362, -20.2680832136679, -0.000144744003193864,
                    -0.000141332438403197, -0.000137619349473539, -0.000108122734189919,
                    -7.17754206647214e-05, -6.32583288151523e-05, -8.13925929321613e-05,
                    -9.74487593623508e-05, -8.84064578059451e-05, -6.71268932208591e-05,
                    -5.69923697817387e-05, -5.29325436939843e-05, -4.79065534442213e-05,
                    -5.27487149497574e-05, -7.91829501594066e-05, -0.000153221427861543,
                    -3.60088035775821 ]
        test_check = all(abs.(arr_e[:] .- arr_o_flN[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        if !test_check
            @error ("Fail 'flN' \n Observed: $(arr_o_flN) \n Expected: $(arr_e)")
        end
        @test test_check

        # convective flux
        arr_e[:] .= [0, 0, 0, 0, 1.62252792680597e-05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.747672235533771, 34.0091829002652,
                    64.9196691838259, 59.0546053518388, 44.59572614017, 20.0908199499065, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
        test_check = all(abs.(arr_e[:] .- arr_o_flC[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        if !test_check
            @error ("Fail 'flC' \n Observed: $(arr_o_flC) \n Expected: $(arr_e)")
        end
        @test test_check

        # eddy diffusion coefficients
        @debug ("Kzz...")
        arr_e[:] .= [26047.1151219131, 20589.3675244423, 16275.2018053592,
                    12864.9990579224, 10169.3486040735, 6581.61713940047, 6581.61713940047,
                    6581.61713940047, 6581.61713940047, 6581.61713940047, 6581.61713940047,
                    6581.61713940047, 6581.61713940047, 6581.61713940047, 6581.61713940047,
                    6581.61713940047, 6581.61713940047, 6581.61713940047, 6581.61713940047,
                    6581.61713940047, 6581.61713940047, 6581.61713940047, 6581.61713940047,
                    6581.61713940047, 6581.61713940047, 6581.61713940047, 6581.61713940047,
                    22561.424214886, 26693.8271438194, 24409.5295421561, 20715.4161920978,
                    14644.7517065071, 14644.7517065071, 14644.7517065071, 14644.7517065071,
                    14644.7517065071, 14644.7517065071, 14644.7517065071, 14644.7517065071,
                    14644.7517065071, 14644.7517065071, 14644.7517065071, 14644.7517065071,
                    14644.7517065071, 14644.7517065071, 14644.7517065071, 14644.7517065071,
                    14644.7517065071, 14644.7517065071 ]
        test_check = all(abs.(arr_e[:] .- arr_o_Kzz[:]) .< abs.(arr_e[:].*rtol) .+ 1e4)
        if !test_check
            @error ("Fail 'Kzz' \n Observed: $(arr_o_Kzz) \n Expected: $(arr_e)")
        end
        @test test_check
    end
end
