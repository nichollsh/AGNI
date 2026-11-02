# Validation calculations for the Generation III spectral files, as shown on the SOCRATES docs page "Validation of
# spectral files". For each family and band count this writes, as CSV files in OUT_DIR:
#   fix_<family>_<nband>_<composition>_p<psurf>.csv   TOA upward flux in each band on prescribed T(p) profiles
#   runaway_<family>_<nband>.csv                      pure-steam OLR on the saturation curve (families with H2O)
#   lbl_<family>_<nband>.csv, lbl_rfm.csv             H2O-CO2 comparison with the RFM line-by-line model (Breezehome)
#
# Usage, from the AGNI root directory:
#   julia --project=. misc/validation/spectral_files/run_validation.jl Snowline Breezehome Volenfell
# The spectral files are read from SF_DIR/<family>/<nband>/<family>.sf (default res/spectral_files), and the results
# are written to OUT_DIR (default out/spectral_file_validation). The figures are made by plot_validation.py.
using Pkg
const AGNI_DIR = abspath(joinpath(@__DIR__, "..", "..", ".."))
Pkg.activate(AGNI_DIR)
using AGNI, Printf, Logging, DelimitedFiles
import AGNI: atmosphere, energy, setpt, chemistry, rfm
global_logger(ConsoleLogger(stderr, Logging.Error))

const SF_DIR = get(ENV, "SF_DIR", joinpath(AGNI_DIR, "res/spectral_files"))
const OUT    = mkpath(get(ENV, "OUT_DIR", joinpath(AGNI_DIR, "out/spectral_file_validation")))
const STAR   = joinpath(AGNI_DIR, "res/stellar_spectra/sun.txt")
const PAR    = joinpath(AGNI_DIR, "res/parfiles/h2o-co2_4000-5000.par")
const NBS    = (32, 64, 128, 256, 2048)
const TS_RUN = collect(200.0:25.0:2500.0)
const SIGMA  = 5.670374419e-8

# Illustrative compositions (mole fractions), not from a chemistry calculation
const COMPS = Dict(
    "steam"   => Dict("H2O" => 1.0),
    "oxid"    => Dict("H2O" => 0.5, "CO2" => 0.3, "CO" => 0.05, "O2" => 0.05, "SO2" => 0.05, "OCS" => 0.03, "S2" => 0.02),
    "reduced" => Dict("H2" => 0.5, "CH4" => 0.15, "N2" => 0.15, "CO" => 0.1, "NH3" => 0.05, "H2S" => 0.04, "S2" => 0.01),
    "rock"    => Dict("O2" => 0.30, "Na" => 0.30, "SiO" => 0.20, "Fe" => 0.10, "MgO" => 0.05, "K" => 0.03, "CO" => 0.0195, "TiO" => 5e-4),
    "h2oco2"  => Dict("H2O" => 0.7, "CO2" => 0.3),
)
# family => (compositions, surface pressures [bar], surface temperatures [K], stratosphere as a fraction of T_surf)
const FAMS = Dict(
    "Snowline"   => (["steam"], [1.0, 100.0], [700.0, 1500.0, 2500.0], 1/3),
    "Breezehome" => (["steam", "oxid"], [1.0, 100.0], [700.0, 1500.0, 2500.0], 1/3),
    "Cliffshade" => (["reduced"], [1.0, 100.0], [700.0, 1500.0, 2500.0], 1/3),
    "Volenfell"  => (["rock"], [0.1, 10.0], [2000.0, 3000.0, 4000.0], 1/2),
)
# Line-by-line comparison: H2O-CO2 at 1 bar and 700 K, wavenumber range [cm-1] and RFM step [cm-1]. The RFM needs
# line data beyond each end of the range (at least its 25 cm-1 line cut-off); the par file covers 4000-5000 cm-1, and the
# RFM stops with an error for numin = 4025 or numax = 4974 cm-1.
# The RFM nadir radiance is taken at secant 1.66, the diffusivity factor of the Elsasser two-stream scheme that AGNI
# uses in the longwave, so that the RFM flux (pi times radiance) and SOCRATES make the same angular approximation.
const LBL = (comp = "h2oco2", ps = 1.0, ts = 700.0, frac = 1/3, numin = 4026.0, numax = 4973.0, nures = 0.002)

sfpath(fam, n) = joinpath(SF_DIR, fam, string(n), "$fam.sf")

# Uniform composition; dry adiabat from the surface up to an isothermal stratosphere at frac * T_surf
function prescribe!(a, comp, ps, ts, frac; iso=false)
    for g in a.gas_names
        a.gas_ovmr[g] .= get(COMPS[comp], g, 0.0)
    end
    a.p_oboa = ps * 1e5; a.tmp_surf = ts
    fill!(a.tmp, ts); fill!(a.tmpl, ts)
    chemistry.restore_composition!(a)
    iso || (setpt.dry_adiabat!(a); setpt.stratosphere!(a, frac * ts))
    atmosphere.calc_layer_props!(a)
end

function new_atmos(fam, n, gases, ps, ts; kwargs...)
    a = atmosphere.Atmos_t()
    atmosphere.setup!(a, AGNI_DIR, mkpath(joinpath(OUT, "run_$(fam)_$n")), sfpath(fam, n), 1361.0, 0.375, 0.0,
                      get(kwargs, :zenith, 54.74), ts, 9.81, 6.371e6, 40, ps, 1e-5, gases, "";
                      flag_rayleigh=false, overlap_method="ee", filter(p -> p.first != :zenith, kwargs)...) ||
        error("setup! failed")
    atmosphere.allocate!(a, STAR; check_safe_gas=false) || error("allocate! failed")
    return a
end

# TOA upward flux in each band on prescribed profiles; an isothermal atmosphere is checked against sigma T^4
function fixed(fam, n)
    comps, PS, TS, frac = FAMS[fam]
    names = unique(g for c in comps for g in keys(COMPS[c]))
    a = new_atmos(fam, n, Dict(g => 1 / length(names) for g in names), PS[1], TS[1]; flag_gcontinuum=true)
    lo, hi = 1e-2 ./ a.bands_max, 1e-2 ./ a.bands_min
    iso = map(TS) do ts
        prescribe!(a, comps[1], PS[1], ts, frac; iso=true); energy.calc_fluxes!(a; radiative=true)
        a.flux_u_lw[1] / (SIGMA * ts^4) - 1
    end
    for comp in comps, ps in PS
        F = map(TS) do ts
            prescribe!(a, comp, ps, ts, frac); energy.calc_fluxes!(a; radiative=true)
            any(!isfinite, a.band_u_lw) && @warn "non-finite flux" fam n comp ps ts
            copy(a.band_u_lw[1, :])
        end
        writedlm(joinpath(OUT, @sprintf("fix_%s_%d_%s_p%g.csv", fam, n, comp, ps)), hcat(lo, hi, F...), ',')
    end
    @printf("%-10s %5d bands: isothermal OLR/σT⁴ - 1 = %s\n", fam, n, join([@sprintf("%+.1e", x) for x in iso], ", "))
    atmosphere.deallocate!(a)
end

# Pure steam on the saturation curve with a 270 bar water inventory (as docs/notebooks/01_canonical-runaway.ipynb)
function runaway(fam, n)
    a = atmosphere.Atmos_t()
    atmosphere.setup!(a, AGNI_DIR, mkpath(joinpath(OUT, "run_$(fam)_$n")), sfpath(fam, n), 1000.0, 1.0, 0.0, 48.19, 1700.0,
                      9.81, 6.371e6, 100, 270.0, 1e-5, Dict("H2O" => 1.0 - 1e-9, "N2" => 1e-9), "";
                      flag_gcontinuum=true, flag_rayleigh=false, overlap_method="ee",
                      thermo_functions=false, real_gas=false, condensates=["H2O"]) || error("setup! failed")
    atmosphere.allocate!(a, STAR; check_safe_gas=false) || error("allocate! failed")
    olr, ps = Float64[], Float64[]
    for ts in TS_RUN
        a.tmp_surf = ts; fill!(a.tmpl, ts); fill!(a.tmp, ts)
        a.p_oboa = 270.0e5; a.p_boa = a.p_oboa
        atmosphere.generate_pgrid!(a); chemistry.restore_composition!(a); chemistry._sat_surf!(a)
        setpt.dry_adiabat!(a); setpt.saturation!(a, "H2O"); atmosphere.calc_layer_props!(a)
        energy.radtrans!(a, true)
        push!(olr, a.flux_u_lw[1]); push!(ps, a.p_boa / 1e5)
    end
    writedlm(joinpath(OUT, "runaway_$(fam)_$n.csv"), hcat(TS_RUN, ps, olr), ',')
    @printf("%-10s %5d bands: runaway OLR at 300/1000/1500/2000/2200 K = %s\n", fam, n,
            join([@sprintf("%.1f", olr[findfirst(==(t), TS_RUN)]) for t in (300.0, 1000.0, 1500.0, 2000.0, 2200.0)], ", "))
    atmosphere.deallocate!(a)
end

# H2O-CO2 without continua, by SOCRATES with this file and (once) by the RFM. The RFM spectrum is integrated over
# each band of every band count, and also saved as means over 0.1 cm-1 for plotting.
function lbl(fam, n; with_rfm=false)
    a = new_atmos(fam, n, COMPS[LBL.comp], LBL.ps, LBL.ts; flag_gcontinuum=false, zenith=acosd(1 / 1.66),
                  rfm_parfile=PAR)
    prescribe!(a, LBL.comp, LBL.ps, LBL.ts, LBL.frac)
    energy.calc_fluxes!(a; radiative=true)
    lo, hi = 1e-2 ./ a.bands_max, 1e-2 ./ a.bands_min
    writedlm(joinpath(OUT, "lbl_$(fam)_$n.csv"), hcat(lo, hi, a.band_u_lw[1, :]), ',')
    if with_rfm
        t = @elapsed rfm.run_rfm(a, LBL.numin, LBL.numax; nures=LBL.nures) || error("RFM failed")
        wn, fl = a.rfm_wn, a.rfm_fl * 1e-3       # cm-1, W m-2 (cm-1)-1
        nbin = round(Int, 0.1 / LBL.nures)
        m = length(wn) ÷ nbin
        writedlm(joinpath(OUT, "lbl_rfm.csv"),
                 hcat(vec(sum(reshape(wn[1:m*nbin], nbin, m), dims=1)) / nbin,
                      vec(sum(reshape(fl[1:m*nbin], nbin, m), dims=1)) / nbin), ',')
        for nb in NBS
            e = readdlm(joinpath(OUT, "lbl_$(fam)_$nb.csv"), ',')
            band = [sum(fl[(wn .>= l) .& (wn .< h)]) * LBL.nures for (l, h) in zip(e[:, 1], e[:, 2])]
            writedlm(joinpath(OUT, "lbl_$(fam)_$(nb).csv"), hcat(e[:, 1:3], band), ',')
        end
        @printf("RFM: %d points from %.1f to %.1f cm-1 in %.0f s\n", length(wn), wn[1], wn[end], t)
    end
    atmosphere.deallocate!(a)
end

for fam in ARGS
    for n in NBS
        try
            fixed(fam, n)
            fam in ("Snowline", "Breezehome") && runaway(fam, n)
            fam == "Breezehome" && lbl(fam, n; with_rfm=(n == last(NBS)))
        catch e
            @error "failed: $fam $n" exception = (e, catch_backtrace())
        end
        flush(stdout)
    end
end
