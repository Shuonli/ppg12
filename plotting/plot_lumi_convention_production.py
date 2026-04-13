#!/usr/bin/env python3
"""
plot_lumi_convention_production.py — per-pT L*eff plot from the
actual production rerun (not the standalone Python simulation).

Reads results/production_rerun_values.pkl and overlays the new vs old
convention for both periods side by side. Overwrites the ratio figure
produced by plot_lumi_convention.py so the report shows the real rerun.
"""

import os
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 13,
    "axes.labelsize": 14,
    "axes.titlesize": 13,
    "legend.fontsize": 11,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.grid": True,
    "grid.alpha": 0.3,
})

SPHENIX_LABEL = r"$\bf{sPHENIX}$ Internal"
PP_LABEL = r"$p$+$p$ $\sqrt{s}$ = 200 GeV"
FIGDIR = "/sphenix/user/shuhangli/ppg12/plotting/figures/lumi_convention"

LUMI = {
    "1p5mrad": {"L60": 16.8790, "Lnv": 17.1642, "label": "1.5 mrad"},
    "0mrad":   {"L60": 32.7034, "Lnv": 47.2284, "label": "0 mrad"},
}


def main():
    with open("/sphenix/user/shuhangli/ppg12/efficiencytool/results/production_rerun_values.pkl", "rb") as fh:
        data = pickle.load(fh)

    fig = plt.figure(figsize=(14, 8))
    for col, period in enumerate(["1p5mrad", "0mrad"]):
        p = data[period]
        pt_lo = p["pt_lo"]
        pt_hi = p["pt_hi"]
        pt_c = 0.5 * (pt_lo + pt_hi)
        pt_w = 0.5 * (pt_hi - pt_lo)
        L60 = LUMI[period]["L60"]
        Lnv = LUMI[period]["Lnv"]
        prod_old = p["prod_old"]
        prod_new = p["prod_new"]
        ratio = p["ratio"]
        label = LUMI[period]["label"]

        # Top: L*eff
        ax1 = plt.subplot2grid((3, 2), (0, col), rowspan=2)
        ax1.errorbar(pt_c - 0.15, prod_old, xerr=pt_w, fmt="bs", ms=7,
                     label=r"OLD: $L_{60\mathrm{cm}} \times \varepsilon_{\mathrm{old}}$")
        ax1.errorbar(pt_c + 0.15, prod_new, xerr=pt_w, fmt="ro", ms=7,
                     label=r"NEW: $L_{\mathrm{noVtx}} \times \varepsilon_{\mathrm{new}}$")
        ax1.set_ylabel(r"$L \times \varepsilon_\mathrm{vtx+MBD}$  [pb$^{-1}$]")
        ax1.set_title(f"{label} ($L_\\mathrm{{60cm}}$={L60:.2f}, $L_\\mathrm{{noVtx}}$={Lnv:.2f})",
                      fontsize=12)
        ax1.legend(loc="upper right", fontsize=10)
        ax1.set_xticklabels([])
        ax1.set_xlim(7, 37)

        if col == 0:
            ax1.text(0.02, 0.97, SPHENIX_LABEL, transform=ax1.transAxes, fontsize=12, va="top")
            ax1.text(0.02, 0.91, PP_LABEL, transform=ax1.transAxes, fontsize=11, va="top")
            ax1.text(0.02, 0.85, "Production re-run (photon5/10/20 merged)",
                     transform=ax1.transAxes, fontsize=10, va="top")

        # Bottom: ratio
        ax2 = plt.subplot2grid((3, 2), (2, col))
        ax2.errorbar(pt_c, ratio, xerr=pt_w, fmt="ko", ms=7)
        ax2.axhline(1.0, color="gray", ls="--", lw=0.8)
        mean_ratio = ratio.mean() if len(ratio) > 0 else 1.0
        ax2.axhline(mean_ratio, color="red", ls="-", lw=1.5, alpha=0.6,
                    label=f"mean = {mean_ratio:.4f} ({(mean_ratio-1)*100:+.2f}%)")
        ax2.set_xlabel(r"$p_T^{\gamma,\mathrm{truth}}$ [GeV]")
        ax2.set_ylabel("NEW / OLD")
        ax2.set_xlim(7, 37)
        if period == "1p5mrad":
            ax2.set_ylim(0.985, 1.005)
        else:
            ax2.set_ylim(1.10, 1.18)
        ax2.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    out = os.path.join(FIGDIR, "L_times_eff_ratio.pdf")
    fig.savefig(out)
    print(f"Saved {out}")
    plt.close()


if __name__ == "__main__":
    main()
