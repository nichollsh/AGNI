"""
Deep Atmospheric Heating test suite (RCE solver)

Runs a baseline (no deep heating) plus a set of deep-heating configurations,
and for each configuration:
1) solves to (approx.) RCE
2) overlays baseline vs heated P–T diagram (saved SVG)
3) checks numerically that the profile changes (ΔT exceeds a threshold)

This suite is *convergence-gated*: results are only compared if the solver
returns `converged=true`. If convergence fails, the case is marked as failed.
"""

println("="^70)
println("TEST SUITE: Deep Atmospheric Heating")
println("="^70)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

ENV["RAD_DIR"] = abspath(joinpath(@__DIR__, "..", "socrates"))
# Headless plotting for GR.
# NOTE: In some headless environments GR writes 0-byte PNGs; SVG is robust.
ENV["GKS_WSTYPE"] = get(ENV, "GKS_WSTYPE", "svg")
ENV["GKSwstype"] = get(ENV, "GKSwstype", "svg")
include(joinpath(@__DIR__, "..", "src", "AGNI.jl"))

using Printf
using Plots
gr()

const ROOT_DIR = abspath(joinpath(@__DIR__, ".."))
const OUT_DIR_BASE = abspath(joinpath(@__DIR__, "..", "out", "test_deep_heating_suite"))
const PLOT_EXT = "svg"
const PLOTS_DIR = joinpath(OUT_DIR_BASE, "plots")
const spectral_file = "res/spectral_files/Dayspring/48/Dayspring.sf"
const star_file = "res/stellar_spectra/sun.txt"

# Planet parameters (warm sub-Neptune like)
const instellation = 3000.0      # W/m²
const tmp_surf_init = 800.0      # K
const gravity = 15.0             # m/s²
const radius = 2.5e7             # m
const p_surf = 100.0             # bar
const p_top = 1e-5               # bar
const nlev = 40                  # reduced for speed

# Solver parameters (convergence-gated; we retry with more robust settings)
const max_steps = 600
const max_runtime = 1200.0       # seconds
const conv_atol = 0.5            # W/m²
const conv_rtol = 5e-3

const SOLVER_ATTEMPTS = [
      (label="LM/backtrack",
       opts=(; method=3, ls_method=2, dx_max=120.0,
             fdw=3.0e-5, fdc=true, fdo=2,
             easy_start=true, detect_plateau=true,
             perturb_all=true, perturb_chem=false,
             ls_increase=1.05, ls_max_steps=25, ls_min_scale=1.0e-5)),
      (label="LM/golden",
       opts=(; method=3, ls_method=1, dx_max=120.0,
             fdw=3.0e-5, fdc=true, fdo=2,
             easy_start=true, detect_plateau=true,
             perturb_all=true, perturb_chem=false,
             ls_increase=1.05, ls_max_steps=25, ls_min_scale=1.0e-5)),
      (label="GN/backtrack",
       opts=(; method=2, ls_method=2, dx_max=80.0,
             fdw=3.0e-5, fdc=true, fdo=2,
             easy_start=true, detect_plateau=true,
             perturb_all=true, perturb_chem=false,
             ls_increase=1.05, ls_max_steps=25, ls_min_scale=1.0e-5)),
]

function get_temp_at_pressure(atmos, p_target::Float64)
      idx = argmin(abs.(atmos.p .- p_target))
      return atmos.tmp[idx], atmos.p[idx]
end

function run_model_rce(out_dir::String;
                                 deep_kwargs::Union{Nothing,NamedTuple}=nothing,
                                 sol_type::Int=1,
                                 verbose::Bool=true)

      rm(out_dir, force=true, recursive=true)
      mkpath(out_dir)

      atmos = AGNI.atmosphere.Atmos_t()
      # Use a dry H2/He composition for solver stability (avoids latent/condensation coupling).
      mf_dict = Dict("H2" => 0.85, "He" => 0.15)

      AGNI.atmosphere.setup!(atmos, ROOT_DIR, out_dir,
                                       spectral_file,
                                       instellation, 1.0, 0.0, 48.19,
                                       tmp_surf_init, gravity, radius,
                                       nlev, p_surf, p_top,
                                       mf_dict, "";
                                       flag_rayleigh=true,
                                       flag_gcontinuum=true,
                                       real_gas=false,
                                       thermo_functions=true)

      AGNI.atmosphere.allocate!(atmos, joinpath(ROOT_DIR, star_file))

      if deep_kwargs !== nothing
            AGNI.atmosphere.set_deep_heating!(atmos; deep_kwargs...)
      end

      AGNI.setpt.isothermal!(atmos, tmp_surf_init)
      AGNI.setpt.dry_adiabat!(atmos)
      AGNI.atmosphere.calc_layer_props!(atmos)

      if verbose
            @info "Starting RCE solver..."
      end

      converged = false
      attempt_label = ""

      for (i, attempt) in enumerate(SOLVER_ATTEMPTS)
            attempt_label = attempt.label
            verbose && @info "Solve attempt $(i)/$(length(SOLVER_ATTEMPTS))" attempt=attempt_label
            converged = AGNI.solver.solve_energy!(atmos,
                                                                  sol_type=sol_type,
                                                                  convect=true,
                                                                  conduct=true,
                                                                  sens_heat=false,
                                                                  latent=false,
                                                                  rainout=false,
                                                                  max_steps=max_steps,
                                                                  max_runtime=max_runtime,
                                                                  conv_atol=conv_atol,
                                                                  conv_rtol=conv_rtol,
                                                                  modplot=0,
                                                                  save_frames=false,
                                                                  plot_jacobian=false,
                                                                  modprint=25;
                                                                  attempt.opts...)
            converged && break
      end

      if verbose
            converged ? @info("Solver converged!", attempt=attempt_label) : @warn("Solver did not fully converge", last_attempt=attempt_label)
      end

      AGNI.energy.calc_fluxes!(atmos, true, false, true, false, true, calc_cf=false)
      AGNI.energy.calc_hrates!(atmos)
      return atmos, converged
end

function plot_pt_overlay(atmos_base, atmos_case;
                                     title_str::String,
                                     outpath::String,
                                     P_dep_pa::Union{Nothing,Float64}=nothing)

      p = plot(xlabel="Temperature [K]", ylabel="Pressure [bar]",
                   yflip=true, yscale=:log10, legend=:bottomleft,
                   title=title_str, size=(650, 500), left_margin=5Plots.mm)
      plot!(p, atmos_base.tmp, atmos_base.p .* 1e-5, label="No deep heating", lw=2.5, color=:blue)
      plot!(p, atmos_case.tmp, atmos_case.p .* 1e-5, label="With deep heating", lw=2.5, color=:red)
      if P_dep_pa !== nothing
            hline!(p, [P_dep_pa * 1e-5], ls=:dot, color=:gray, lw=1.5, label="P_dep")
      end
      savefig(p, outpath)
      if stat(outpath).size == 0
            error("Plot written as 0 bytes: $outpath")
      end
      return nothing
end

function plot_dt(atmos_base, atmos_case;
                         title_str::String,
                         outpath::String,
                         P_dep_pa::Union{Nothing,Float64}=nothing)

      ΔT = atmos_case.tmp .- atmos_base.tmp
      p = plot(xlabel="ΔT [K] (with - without)", ylabel="Pressure [bar]",
                   yflip=true, yscale=:log10, legend=false,
                   title=title_str, size=(650, 500), left_margin=5Plots.mm)
      plot!(p, ΔT, atmos_case.p .* 1e-5, lw=2.5, color=:purple)
      vline!(p, [0.0], ls=:dash, color=:black, lw=1)
      if P_dep_pa !== nothing
            hline!(p, [P_dep_pa * 1e-5], ls=:dot, color=:gray, lw=1.5)
      end
      savefig(p, outpath)
      if stat(outpath).size == 0
            error("Plot written as 0 bytes: $outpath")
      end
      return nothing
end

function plot_flux_overlay(atmos_base, atmos_case;
                                     title_str::String,
                                     outpath::String,
                                     which::Symbol=:net,
                                     P_dep_pa::Union{Nothing,Float64}=nothing)

      flux_base = which === :net ? atmos_base.flux_n : which === :up ? atmos_base.flux_u : atmos_base.flux_d
      flux_case = which === :net ? atmos_case.flux_n : which === :up ? atmos_case.flux_u : atmos_case.flux_d
      label_base = which === :net ? "Net flux (no deep heating)" : which === :up ? "Up flux (no deep heating)" : "Down flux (no deep heating)"
      label_case = which === :net ? "Net flux (with deep heating)" : which === :up ? "Up flux (with deep heating)" : "Down flux (with deep heating)"

      p = plot(xlabel="Flux [W m⁻²]", ylabel="Pressure [bar]",
                   yflip=true, yscale=:log10, legend=:bottomleft,
                   title=title_str, size=(650, 500), left_margin=5Plots.mm)
      plot!(p, flux_base, atmos_base.pl .* 1e-5, label=label_base, lw=2.5, color=:blue)
      plot!(p, flux_case, atmos_case.pl .* 1e-5, label=label_case, lw=2.5, color=:red)
      vline!(p, [0.0], ls=:dash, color=:black, lw=1, label=false)
      if P_dep_pa !== nothing
            hline!(p, [P_dep_pa * 1e-5], ls=:dot, color=:gray, lw=1.5, label="P_dep")
      end
      savefig(p, outpath)
      if stat(outpath).size == 0
            error("Plot written as 0 bytes: $outpath")
      end
      return nothing
end

function plot_dflux(atmos_base, atmos_case;
                              title_str::String,
                              outpath::String,
                              which::Symbol=:net,
                              P_dep_pa::Union{Nothing,Float64}=nothing)

      flux_base = which === :net ? atmos_base.flux_n : which === :up ? atmos_base.flux_u : atmos_base.flux_d
      flux_case = which === :net ? atmos_case.flux_n : which === :up ? atmos_case.flux_u : atmos_case.flux_d
      ΔF = flux_case .- flux_base
      p = plot(xlabel="ΔFlux [W m⁻²] (with - without)", ylabel="Pressure [bar]",
                   yflip=true, yscale=:log10, legend=false,
                   title=title_str, size=(650, 500), left_margin=5Plots.mm)
      plot!(p, ΔF, atmos_case.pl .* 1e-5, lw=2.5, color=:darkgreen)
      vline!(p, [0.0], ls=:dash, color=:black, lw=1)
      if P_dep_pa !== nothing
            hline!(p, [P_dep_pa * 1e-5], ls=:dot, color=:gray, lw=1.5)
      end
      savefig(p, outpath)
      if stat(outpath).size == 0
            error("Plot written as 0 bytes: $outpath")
      end
      return nothing
end

function relerr(obs::Float64, exp::Float64)
      exp == 0.0 && return abs(obs)
      return abs(obs - exp) / abs(exp)
end

mkpath(OUT_DIR_BASE)
rm(PLOTS_DIR, force=true, recursive=true)
mkpath(PLOTS_DIR)

@info "Running BASELINE (no deep heating)..."
base_dir = joinpath(OUT_DIR_BASE, "baseline")
t0 = time()
atmos_base, conv_base = run_model_rce(base_dir; deep_kwargs=nothing, sol_type=1)
@printf("Baseline runtime: %.1f s (converged=%s)\n", time() - t0, conv_base ? "yes" : "no")

baseline_cache = Dict{Int, Tuple{AGNI.atmosphere.Atmos_t, Bool}}(1 => (atmos_base, conv_base))

if !conv_base
      AGNI.atmosphere.deallocate!(atmos_base)
      error("Baseline (sol_type=1) did not converge; aborting deep heating suite.")
end

Teq_est = (instellation * (1.0 - 0.0) * 1.0 / 4.0 / AGNI.phys.σSB)^(0.25)
@printf("Estimated Teq for ohmic proxy: %.2f K\n", Teq_est)

const P_dep_default = 1.0e6
const sigmaP_default = 1.5

cases = [
      (name="generic_eff_pressure",
       label="generic / efficiency / pressure",
      kwargs=(; active=true, P_dep=P_dep_default, sigma_P=sigmaP_default, efficiency=0.01,
                    mechanism=:generic, normalization=:pressure, below_domain=:clamp, power_mode=:efficiency),
      expected_flux=0.01 * instellation,
      dT_min=0.2,
      dF_min=1.0,
      sol_type=1),

      (name="generic_eff_mass",
       label="generic / efficiency / mass",
            kwargs=(; active=true, P_dep=P_dep_default, sigma_P=sigmaP_default, efficiency=0.01,
                    mechanism=:generic, normalization=:mass, below_domain=:clamp, power_mode=:efficiency),
            expected_flux=0.01 * instellation,
            dT_min=0.2,
                dF_min=1.0,
            sol_type=1),

      (name="generic_flux",
       label="generic / flux mode",
       kwargs=(; active=true, P_dep=P_dep_default, sigma_P=sigmaP_default,
                    mechanism=:generic, normalization=:mass, below_domain=:clamp,
                    power_mode=:flux, F_total=100.0),
       expected_flux=100.0,
            dT_min=0.2,
            dF_min=1.0,
            sol_type=1),

      (name="generic_power",
       label="generic / power mode",
       kwargs=(; active=true, P_dep=P_dep_default, sigma_P=sigmaP_default,
                    mechanism=:generic, normalization=:mass, below_domain=:clamp,
                    power_mode=:power, power=5.0e17),
       expected_flux=5.0e17 / (4.0 * π * radius^2),
            dT_min=0.2,
            dF_min=1.0,
            sol_type=1),

      (name="ohmic_eff",
       label="ohmic proxy / efficiency",
            kwargs=(; active=true, P_dep=P_dep_default, sigma_P=sigmaP_default, efficiency=0.01,
                    mechanism=:ohmic, normalization=:mass, below_domain=:clamp, power_mode=:efficiency,
                    ohmic_Tpeak=Teq_est, ohmic_sigmaT=200.0),
            expected_flux=0.01 * instellation,  # Teq matches Tpeak => proxy ~1
            dT_min=0.2,
                dF_min=1.0,
            sol_type=1),

      (name="tidal_power_auto",
       label="tidal / auto power",
       kwargs=(; active=true, P_dep=P_dep_default, sigma_P=sigmaP_default,
                    mechanism=:tidal, normalization=:mass, below_domain=:clamp,
                    power_mode=:power, power=0.0,
                    tidal_e=0.05, tidal_a=4.0e9, tidal_Mstar=1.989e30, tidal_k2=0.3, tidal_Q=1.0e6),
       expected_flux=nothing,
            dT_min=0.2,
            dF_min=1.0,
            sol_type=1),

      (name="boundary_flux_outside",
       label="boundary flux (P_dep below domain)",
       kwargs=(; active=true, P_dep=1.0e12, sigma_P=sigmaP_default,
                    mechanism=:generic, normalization=:mass, below_domain=:boundary_flux,
                        power_mode=:flux, F_total=50.0),
            expected_flux=50.0,
            dT_min=0.2,
                  dF_min=1.0,
            sol_type=3),
]

failures = String[]

for c in cases
      println("\n" * "-"^70)
      @info "Running case: $(c.name) ($(c.label))"

      sol_type = getproperty(c, :sol_type)

      if !haskey(baseline_cache, sol_type)
            @info "Running BASELINE for sol_type=$sol_type (no deep heating)..."
            bdir = joinpath(OUT_DIR_BASE, "baseline_soltype_$(sol_type)")
            b_atm, b_conv = run_model_rce(bdir; deep_kwargs=nothing, sol_type=sol_type)
            baseline_cache[sol_type] = (b_atm, b_conv)
            if !b_conv
                  AGNI.atmosphere.deallocate!(b_atm)
                  error("Baseline (sol_type=$(sol_type)) did not converge; aborting deep heating suite.")
            end
      end

      atmos_base_case, _ = baseline_cache[sol_type]
      case_dir = joinpath(OUT_DIR_BASE, "case_" * c.name)
      t1 = time()
      atmos_case, conv_case = run_model_rce(case_dir; deep_kwargs=c.kwargs, sol_type=sol_type)
      @printf("Case runtime: %.1f s (converged=%s)\n", time() - t1, conv_case ? "yes" : "no")

      if !conv_case
            push!(failures, "$(c.name): solver did not converge")
            AGNI.atmosphere.deallocate!(atmos_case)
            continue
      end

      # Basic numerical checks
      if length(atmos_case.tmp) != length(atmos_base_case.tmp)
            push!(failures, "$(c.name): temperature array length mismatch")
            AGNI.atmosphere.deallocate!(atmos_case)
            continue
      end

      ΔT = atmos_case.tmp .- atmos_base_case.tmp
      max_abs_dT = maximum(abs.(ΔT))
      flux_deep_boa = atmos_case.flux_deep[end]

      ΔFnet = atmos_case.flux_n .- atmos_base_case.flux_n
      max_abs_dFnet = maximum(abs.(ΔFnet))

      @printf("Deep flux at BOA: %.3f W/m²\n", flux_deep_boa)
      @printf("Max |ΔT|: %.3f K\n", max_abs_dT)
      @printf("Max |ΔF_net|: %.3f W/m²\n", max_abs_dFnet)

      # Optional flux expectation check
      if c.expected_flux !== nothing
            r = relerr(flux_deep_boa, c.expected_flux)
            @printf("Expected deep flux: %.3f W/m² (relerr=%.3f)\n", c.expected_flux, r)
            if r > 0.10
                  push!(failures, "$(c.name): deep flux mismatch (obs=$(flux_deep_boa), exp=$(c.expected_flux), relerr=$(r))")
            end
      end

      if !(max_abs_dT >= c.dT_min)
            push!(failures, "$(c.name): ΔT too small (max|ΔT|=$(max_abs_dT) < $(c.dT_min))")
      end

      if !(max_abs_dFnet >= c.dF_min)
            push!(failures, "$(c.name): ΔF_net too small (max|ΔF|=$(max_abs_dFnet) < $(c.dF_min))")
      end

      # Save PT overlays (baseline vs this case)
      Pdep = haskey(c.kwargs, :P_dep) ? getproperty(c.kwargs, :P_dep) : nothing
      plot_pt_overlay(atmos_base_case, atmos_case,
                              title_str="P–T: $(c.name)",
                              outpath=joinpath(PLOTS_DIR, "pt_$(c.name).$(PLOT_EXT)"),
                              P_dep_pa=Pdep === nothing ? nothing : Float64(Pdep))
      plot_dt(atmos_base_case, atmos_case,
                  title_str="ΔT: $(c.name)",
                  outpath=joinpath(PLOTS_DIR, "dT_$(c.name).$(PLOT_EXT)"),
                  P_dep_pa=Pdep === nothing ? nothing : Float64(Pdep))

      plot_flux_overlay(atmos_base_case, atmos_case,
                              title_str="Flux(net): $(c.name)",
                              outpath=joinpath(PLOTS_DIR, "fluxnet_$(c.name).$(PLOT_EXT)"),
                              which=:net,
                              P_dep_pa=Pdep === nothing ? nothing : Float64(Pdep))
      plot_dflux(atmos_base_case, atmos_case,
                  title_str="ΔFlux(net): $(c.name)",
                  outpath=joinpath(PLOTS_DIR, "dFluxnet_$(c.name).$(PLOT_EXT)"),
                  which=:net,
                  P_dep_pa=Pdep === nothing ? nothing : Float64(Pdep))

      # Print a quick comparison at a few pressures
      T1_off, _ = get_temp_at_pressure(atmos_base_case, 1e5)
      T1_on, _ = get_temp_at_pressure(atmos_case, 1e5)
      T10_off, _ = get_temp_at_pressure(atmos_base_case, 1e6)
      T10_on, _ = get_temp_at_pressure(atmos_case, 1e6)
      @printf("ΔT @ 1 bar:  %+.2f K\n", T1_on - T1_off)
      @printf("ΔT @ 10 bar: %+.2f K\n", T10_on - T10_off)

      AGNI.atmosphere.deallocate!(atmos_case)
end

# Write summary
summary_file = joinpath(OUT_DIR_BASE, "summary.txt")
open(summary_file, "w") do f
      println(f, "Deep heating test suite summary")
      println(f, "==============================")
      println(f, "Baseline converged: $(conv_base)")
      println(f, "Cases run: $(length(cases))")
      println(f, "Plots: $(PLOTS_DIR)")
      println(f, "Plot format: $(PLOT_EXT)")
      println(f, "")
      if isempty(failures)
            println(f, "PASS: all cases changed P–T and matched flux expectations.")
      else
            println(f, "FAIL: $(length(failures)) issue(s)")
            for s in failures
                  println(f, "- " * s)
            end
      end
end

for (_, (b_atm, _)) in baseline_cache
      AGNI.atmosphere.deallocate!(b_atm)
end

println("\n" * "="^70)
@info "Plots saved under: $(PLOTS_DIR)"
@info "Summary written to: $summary_file"

if !isempty(failures)
      error("Deep heating test suite FAILED. See summary at $summary_file")
end

@info "Deep heating test suite PASSED"
println("="^70)
