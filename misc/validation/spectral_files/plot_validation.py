# Figures and numbers for the SOCRATES docs page "Validation of spectral files", from the CSV files written by
# run_validation.jl: the pure-steam runaway curve, the OLR against band count, a band-resolved spectrum, and the
# comparison with the RFM line-by-line model. Each figure is drawn once for each colour scheme of the docs.
#
# Usage, from the AGNI root directory:
#   python misc/validation/spectral_files/plot_validation.py [figure directory]
# The results are read from OUT_DIR (default out/spectral_file_validation); the figures go to the given directory
# (default OUT_DIR/figures), and the numbers quoted on the docs page are printed.
import os
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

AGNI_DIR = Path(__file__).resolve().parents[3]
OUT = Path(os.environ.get("OUT_DIR", AGNI_DIR / "out" / "spectral_file_validation"))
LIT = AGNI_DIR / "res" / "literature_data" / "runaway"
NBS = [32, 64, 128, 256, 2048]
LBL_RANGE = (4026.0, 4973.0)  # bands compared with the RFM: those within its range (see run_validation.jl)

# Colours from the docs palette (SOCRATES docs/stylesheets/extra.css); the dark scheme has its own steps
THEMES = {
    "light": dict(bg="#F2F5F7", ink="#1C232C", muted="#56606C", grid="#D3DAE1",
                  bands=["#4590C8", "#2A7BBA", "#1B6FA8", "#14406B", "#081A2E"],
                  temps=["#B87A0C", "#E0602A", "#A8261C"],
                  lit=["#C2362B", "#9A6A08", "#593E74", "#56606C"]),
    "dark": dict(bg="#05070B", ink="#E9EEF2", muted="#9AA5B4", grid="#1A2230",
                 bands=["#1B6FA8", "#2E86C8", "#4FA3D9", "#86C1E6", "#D4EAF6"],
                 temps=["#E0A32E", "#F07A3A", "#F04A3A"],
                 lit=["#F0604F", "#E0A32E", "#A98BD0", "#9AA5B4"]),
}
LITS = [("Goldblatt et al. (2013)", "Goldblatt13_data.txt", (0, (5, 2))),
        ("Kopparapu et al. (2013)", "Kopparapu13_data.txt", (0, (1, 1.5))),
        ("Hamano et al. (2015)", "Hamano15_data.txt", (0, (5, 1.5, 1, 1.5))),
        ("Selsis et al. (2023)", "Selsis23_convective.txt", (0, (8, 2, 1, 2, 1, 2)))]
# family, composition, label, surface pressures [bar], surface temperatures [K]
SETS = [("Snowline", "steam", "Snowline: H$_2$O", [1.0, 100.0], [700, 1500, 2500]),
        ("Breezehome", "oxid", "Breezehome: seven gases", [1.0, 100.0], [700, 1500, 2500]),
        ("Volenfell", "rock", "Volenfell: rock vapour", [0.1, 10.0], [2000, 3000, 4000])]


def luminance(c):
    r = np.array([int(c[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    r = np.where(r <= 0.03928, r / 12.92, ((r + 0.055) / 1.055) ** 2.4)
    return 0.2126 * r[0] + 0.7152 * r[1] + 0.0722 * r[2]


def contrast(a, b):
    la, lb = sorted([luminance(a), luminance(b)], reverse=True)
    return (la + 0.05) / (lb + 0.05)


def style(th):
    plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 8.5, "axes.titlelocation": "left", "legend.fontsize": 7.5,
                         "text.color": th["ink"], "axes.labelcolor": th["ink"], "axes.edgecolor": th["muted"],
                         "xtick.color": th["muted"], "ytick.color": th["muted"], "axes.facecolor": "none",
                         "figure.facecolor": "none", "legend.frameon": False, "axes.grid": True,
                         "grid.color": th["grid"], "grid.linewidth": 0.6, "axes.spines.top": False,
                         "axes.spines.right": False, "lines.solid_capstyle": "round"})


def runaway(fam, n):
    d = np.loadtxt(OUT / f"runaway_{fam}_{n}.csv", delimiter=",")
    return d[:, 0], np.where(d[:, 1] > 1e-4, d[:, 2], np.nan)


def lit(f):
    d = np.loadtxt(LIT / f, delimiter=",", skiprows=2)
    o = np.argsort(d[:, 0])
    return d[o, 0], d[o, 1]


def fixed(fam, n, comp, ps):
    d = np.loadtxt(OUT / f"fix_{fam}_{n}_{comp}_p{ps:g}.csv", delimiter=",")
    return d[:, 0], d[:, 1], d[:, 2:]


def lbl(n):
    # band limits [cm-1], SOCRATES and RFM band fluxes [W m-2], for the bands compared with the RFM
    d = np.loadtxt(OUT / f"lbl_Breezehome_{n}.csv", delimiter=",")
    k = (d[:, 0] >= LBL_RANGE[0]) & (d[:, 1] <= LBL_RANGE[1])
    return d[k, 0], d[k, 1], d[k, 2], d[k, 3]


def steps(lo, hi, y):
    o = np.argsort(lo)
    return np.ravel(np.column_stack([lo[o], hi[o]])), np.ravel(np.column_stack([y[o], y[o]]))


def fig_runaway(path, th):
    fig, axs = plt.subplots(1, 2, figsize=(8.2, 3.3))
    for ax in axs:
        for (lab, f, ls), c in zip(LITS, th["lit"]):
            ax.plot(*lit(f), color=c, lw=1.3, ls=ls, label=lab)
        for n, c in zip(NBS, th["bands"]):
            ax.plot(*runaway("Snowline", n), color=c, lw=2.2 if n == 2048 else 1.3, label=f"Snowline, {n} bands",
                    zorder=3 if n == 2048 else 2)
        ax.set_xlabel("Surface temperature [K]")
        ax.set_ylabel("OLR [W m$^{-2}$]")
    axs[0].set(xlim=(250, 1700), ylim=(240, 320), title="Runaway limit")
    axs[1].set(xlim=(200, 2500), ylim=(100, 1500), yscale="log", title="Full range (logarithmic OLR)")
    axs[1].set_yticks([100, 200, 300, 500, 1000])
    axs[1].set_yticklabels(["100", "200", "300", "500", "1000"])
    axs[1].minorticks_off()
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.01), handlelength=3)
    fig.tight_layout(rect=(0, 0.2, 1, 1))
    fig.savefig(path, dpi=160, transparent=True)
    plt.close(fig)


def fig_bands(path, th):
    fig, axs = plt.subplots(2, 3, figsize=(8.2, 4.6), sharex=True)
    for col, (fam, comp, lab, PS, TS) in enumerate(SETS):
        for row, ps in enumerate(PS):
            ax = axs[row, col]
            ref = fixed(fam, 2048, comp, ps)[2].sum(0)
            for k, (t, c) in enumerate(zip(TS, th["temps"])):
                d = [100 * (fixed(fam, n, comp, ps)[2][:, k].sum() / ref[k] - 1) for n in NBS[:-1]]
                ax.plot(NBS[:-1], d, "o-", color=c, lw=1.4, ms=4.5,
                        label=f"{t} K (Volenfell {SETS[2][4][k]} K)")
            ax.axhline(0, color=th["muted"], lw=0.8, ls=(0, (4, 2)), zorder=1)
            ax.set_xscale("log", base=2)
            ax.set_xticks(NBS[:-1])
            ax.set_xticklabels([str(n) for n in NBS[:-1]])
            ax.minorticks_off()
            ax.set_title(f"{lab}, {ps:g} bar" if row == 0 else f"{ps:g} bar")
            if col == 0:
                ax.set_ylabel("OLR difference from\n2048 bands [%]")
            if row == 1:
                ax.set_xlabel("Number of bands")
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, title="Surface temperature", title_fontsize=7.5,
               bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(path, dpi=160, transparent=True)
    plt.close(fig)


def fig_spectrum(path, th, fam="Breezehome", comp="oxid", ps=1.0, k=1):
    fig, ax = plt.subplots(figsize=(8.2, 3.2))
    ymax = 0.0
    for n, c in [(32, th["bands"][0]), (256, th["bands"][3]), (2048, th["muted"])]:
        lo, hi, F = fixed(fam, n, comp, ps)
        wl0, wl1 = 1e4 / hi, 1e4 / lo
        y = F[:, k] / (wl1 - wl0)
        ymax = max(ymax, y.max())
        ax.plot(*steps(wl0, wl1, y), color=c, lw=0.8 if n == 2048 else 1.4, label=f"{n} bands",
                zorder=1 if n == 2048 else 2)
    ax.set(xscale="log", yscale="log", xlim=(0.5, 60), ylim=(ymax * 1e-4, ymax * 2), xlabel="Wavelength [µm]",
           ylabel="Outgoing flux [W m$^{-2}$ µm$^{-1}$]",
           title="Breezehome: seven gases, 1 bar, surface temperature 1500 K")
    ax.set_xticks([0.5, 1, 2, 5, 10, 20, 50])
    ax.set_xticklabels(["0.5", "1", "2", "5", "10", "20", "50"])
    ax.legend(loc="lower center")
    fig.tight_layout()
    fig.savefig(path, dpi=160, transparent=True)
    plt.close(fig)


def fig_lbl(path, th):
    fig, (a, b) = plt.subplots(2, 1, figsize=(8.2, 4.6), sharex=True, gridspec_kw={"height_ratios": [2.4, 1]})
    wn, fl = np.loadtxt(OUT / "lbl_rfm.csv", delimiter=",").T
    a.plot(wn, fl, color=th["muted"], lw=0.4, label="RFM line-by-line (0.1 cm$^{-1}$ means)", zorder=1)
    lo, hi, soc, ref = lbl(2048)
    a.plot(0.5 * (lo + hi), ref / (hi - lo), "o", ms=2.6, mfc="none", mec=th["ink"], mew=0.8,
           label="RFM, mean over each band of the 2048-band file", zorder=3)
    for n, c in [(256, th["bands"][3]), (2048, th["bands"][4])]:
        lo, hi, soc, ref = lbl(n)
        a.plot(*steps(lo, hi, soc / (hi - lo)), color=c, lw=1.4, label=f"Breezehome, {n} bands", zorder=2)
    a.set(yscale="log", ylim=(1e-11, 1), ylabel="Outgoing flux\n[W m$^{-2}$ (cm$^{-1}$)$^{-1}$]",
          title="H$_2$O 0.7, CO$_2$ 0.3, 1 bar, surface temperature 700 K, no continua")
    a.legend(loc="lower center", ncol=2)
    for n, c in zip(NBS[1:], th["bands"][1:]):
        lo, hi, soc, ref = lbl(n)
        b.plot(*steps(lo, hi, 100 * (soc / ref - 1)), color=c, lw=1.3 if n < 2048 else 0.9, label=f"{n} bands")
    b.axhline(0, color=th["muted"], lw=0.8, ls=(0, (4, 2)), zorder=1)
    b.set(xlim=(4000, 5000), xlabel="Wavenumber [cm$^{-1}$]", ylabel="Difference\nfrom RFM [%]")
    b.legend(loc="lower left", ncol=4)
    fig.tight_layout()
    fig.savefig(path, dpi=160, transparent=True)
    plt.close(fig)


def numbers():
    TQ = [400, 1000, 1500, 2000, 2200]
    print("Runaway table, OLR [W m-2] at " + ", ".join(f"{t} K" for t in TQ))
    for n in NBS:
        T, o = runaway("Snowline", n)
        print(f"| Snowline, {n} bands | " + " | ".join(f"{np.interp(t, T, o):.0f}" for t in TQ) + " |")
    for lab, f, _ in LITS:
        T, o = lit(f)
        v = [np.interp(t, T, o, left=np.nan, right=np.nan) for t in TQ]
        print(f"| {lab} | " + " | ".join("" if np.isnan(x) else f"{x:.0f}" for x in v) + " |")
    if (OUT / "runaway_Breezehome_2048.csv").exists():
        d = max(np.nanmax(np.abs(runaway("Breezehome", n)[1] - runaway("Snowline", n)[1])) for n in NBS)
        print(f"largest |Breezehome - Snowline| over all T and band counts: {d:.3g} W m-2")
    T = runaway("Snowline", 2048)[0]
    m = (T >= 300) & (T <= 1500)
    print("Snowline runaway limit, 300-1500 K, all band counts: %.1f-%.1f W m-2"
          % (min(np.nanmin(runaway("Snowline", n)[1][m]) for n in NBS),
             max(np.nanmax(runaway("Snowline", n)[1][m]) for n in NBS)))
    for lab, f, _ in LITS:
        TL, oL = lit(f)
        r = []
        for n in NBS:
            Ts, o = runaway("Snowline", n)
            k = (Ts >= max(400, TL.min())) & (Ts <= min(1500, TL.max()))
            r.append(100 * (o[k] / np.interp(Ts[k], TL, oL) - 1))
        r = np.concatenate(r)
        print(f"Snowline / {lab} - 1, 400-1500 K: {r.min():+.1f} to {r.max():+.1f} %")
    for t in (1700, 1800, 2200):
        v = {n: np.interp(t, *runaway("Snowline", n)) for n in NBS}
        print(f"Snowline OLR at {t} K: " + ", ".join(f"{n}: {x:.0f}" for n, x in v.items()))

    print("\nOLR difference from 2048 bands [%]")
    for fam, comp, lab, PS, TS in SETS:
        for ps in PS:
            ref = fixed(fam, 2048, comp, ps)[2].sum(0)
            print(f"{fam:10s} {ps:5g} bar, T_surf {TS} K")
            for n in NBS[:-1]:
                print(f"    {n:5d} bands " + " ".join(f"{x:+7.2f}" for x in 100 * (fixed(fam, n, comp, ps)[2].sum(0) / ref - 1)))

    print(f"\nBreezehome against the RFM, bands within {LBL_RANGE[0]:g}-{LBL_RANGE[1]:g} cm-1")
    for n in NBS:
        lo, hi, soc, ref = lbl(n)
        if len(lo) == 0:
            print(f"    {n:5d} bands: no band lies within the range")
            continue
        r = 100 * (soc / ref - 1)
        print(f"    {n:5d} bands: {len(lo):3d} bands, {lo.min():.0f}-{hi.max():.0f} cm-1, total flux {soc.sum():.3f} "
              f"against {ref.sum():.3f} W m-2 ({100 * (soc.sum() / ref.sum() - 1):+.1f} %), "
              f"per band {r.min():+.1f} to {r.max():+.1f} %, median |difference| {np.median(np.abs(r)):.1f} %")


if __name__ == "__main__":
    dst = Path(sys.argv[1]) if len(sys.argv) > 1 else OUT / "figures"
    dst.mkdir(parents=True, exist_ok=True)
    for scheme, th in THEMES.items():
        for key in ("bands", "temps", "lit"):
            low = [c for c in th[key] if contrast(c, th["bg"]) < 3]
            if low:
                print(f"warning: {scheme} {key} colours below 3:1 contrast with the page: {low}")
        style(th)
        fig_runaway(dst / f"runaway_{scheme}.png", th)
        fig_bands(dst / f"bands_{scheme}.png", th)
        fig_spectrum(dst / f"spectrum_{scheme}.png", th)
        if (OUT / "lbl_rfm.csv").exists():
            fig_lbl(dst / f"lbl_{scheme}.png", th)
    print("figures saved to", dst)
    numbers()
