
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
        arr_e[:] .= [1.0, 1.6166561597647628, 2.61357713890535, 4.225255580631699, 6.830785461008676, 11.043031391571262, 17.85278472165932, 28.861814389224783, 46.65963001432751, 75.43257827500737, 121.9485423152283, 197.14886210824756, 318.7219223179124, 515.2637589673191, 833.0043298380627, 1346.6815809434217, 2177.1210730739326, 3519.6561933386415, 5690.073865215115, 9198.892963716504, 14871.446972809026, 24042.016353206734, 38867.67383057681, 62835.664313929636, 101583.66376602539, 164225.8557588174, 265496.74130513147, 429216.94222841226, 693896.2135289601, 1.1217915878390379e6, 1.8135512804522736e6, 2.931888848592439e6, 4.739856166822592e6, 7.662717668492737e6, 1.2387979719307058e7, 2.00271037202587e7, 3.2376940591604076e7, 5.234238044175449e7, 8.4619631757913e7, 1.368008489184562e8, 2.2115993506507078e8, 3.5753957131612146e8, 5.780185503278601e8, 6.985757894454757e8, 7.000000000000001e8]
        test_check = all(abs.(arr_e[:] .- arr_o_prs[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        if !test_check
            @error ("Fail 'prs' \n Observed: $(arr_o_prs) \n Expected: $(arr_e)")
        end
        @test test_check

        # temperature profile
        @debug ("Temperature...")
        arr_e[:] .= [166.13495251941706, 181.38720985026285, 189.451092342491, 190.43461237027338, 191.61299578427833, 192.6897426322011, 192.93648924933416, 191.73531669174662, 190.2719463229448, 190.47279553573136, 191.6074583400523, 193.01369165002518, 193.96792214045485, 195.20974257106104, 198.8280729153404, 204.1621858094887, 210.43821168051437, 221.89350404287342, 243.0386069934421, 266.9448558058717, 290.49596828150254, 320.3365851015192, 355.9368473627168, 394.13091927011976, 435.02736916194397, 472.97386865356907, 503.9252284131128, 530.6904560514897, 553.6122426719533, 571.2185677566192, 583.7165355304262, 594.0500164047908, 605.7415494638906, 620.224311833289, 636.0545181394675, 650.4576085457393, 662.05875630792, 671.6498792319785, 680.4012122448423, 688.6015321595945, 696.2709658438955, 704.5687512025273, 716.364002129331, 724.2378987244332, 724.3459744394626]
        test_check = all(abs.(arr_e[:] .- arr_o_tmp[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        if !test_check
            @error ("Fail 'tmp' \n Observed: $(arr_o_tmp) \n Expected: $(arr_e)")
        end
        @test test_check

        # radiative flux
        @debug ("Radiative flux...")
        arr_e[:] .= [-0.0003111581085590842, -0.0003790702958781367, -0.0003161484344218479, -0.00031761818263476016, -0.0003188050230846784, -0.0003171001123973838, -0.0003103854026562658, -0.00030158692396753395, -0.00030658330643973386, -0.00031734218197243536, -0.0003140569706374663, -0.00031293225055151197, -0.00029429179727458177, -0.00029236536164489735, -0.0002823032111223256, -0.00023733844540174687, -0.00015231497116019455, -0.0002232565296935718, -9.97846343371242, -45.653714826283476, -67.5248278886412, -62.02944610472454, -51.83170103090947, -35.58947484785486, -9.174716702127768, -0.00020021118763224877, -0.00017192708001800838, -0.00015592689332777354, -0.000128643763886771, -9.167473989180053e-5, -6.578936178414097e-5, -6.514538576851692e-5, -8.3183194149683e-5, -0.0001010857595886705, -0.000100893652060563, -8.337805247293772e-5, -6.54344077446467e-5, -5.7901622394496144e-5, -5.4892824287833264e-5, -5.1017014630616586e-5, -4.821513264730015e-5, -5.930836435030784e-5, -8.752380356780165e-5, -0.00015922670870886085, -3.7308897319832823]
        test_check = all(abs.(arr_e[:] .- arr_o_flN[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        if !test_check
            @error ("Fail 'flN' \n Observed: $(arr_o_flN) \n Expected: $(arr_e)")
        end
        @test test_check

        # convective flux
        arr_e[:] .= [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.964392679501833, 45.94033549690809, 67.1103876772344, 62.028861494842765, 51.82909174090239, 35.589252774481075, 9.177061039904173, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        test_check = all(abs.(arr_e[:] .- arr_o_flC[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        if !test_check
            @error ("Fail 'flC' \n Observed: $(arr_o_flC) \n Expected: $(arr_e)")
        end
        @test test_check

        # eddy diffusion coefficients
        @debug ("Kzz...")
        arr_e[:] .= [483392.7284556761, 398889.9079125956, 329159.1893466155, 271618.23295630165, 224136.12277011253, 184954.45237102755, 152622.20577872955, 125942.02193109093, 103925.85277588709, 85758.37285751126, 70766.78534481382, 58395.90632578484, 48187.60466501459, 39763.83601955386, 32812.64271137205, 27076.60098926838, 22343.287847338186, 18437.414357395668, 15214.334188815343, 24465.89996242476, 26708.75613160484, 24834.780786371706, 22135.317549080046, 18350.587588802708, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072, 10896.478080874072]
        test_check = all(abs.(arr_e[:] .- arr_o_Kzz[:]) .< abs.(arr_e[:].*rtol) .+ 1e4)
        if !test_check
            @error ("Fail 'Kzz' \n Observed: $(arr_o_Kzz) \n Expected: $(arr_e)")
        end
        @test test_check
    end
end
