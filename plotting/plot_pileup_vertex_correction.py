#!/usr/bin/env python3
"""
plot_pileup_vertex_correction.py — Plots for the Section 3 pileup vertex acceptance study.

Reads results/pileup_vertex_acceptance.csv and produces publication-quality figures.
Uses matplotlib with sPHENIX-style labeling.

Figures:
  1. Per-pT vertex_eff: single vs double interaction
  2. Per-pT correction factor for 0 mrad and 1.5 mrad
  3. MBD + vertex acceptance breakdown (single vs double)
  4. MBD hit decomposition stacked bars
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import csv

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 14,
    "legend.fontsize": 12,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.grid": True,
    "grid.alpha": 0.3,
})

SPHENIX_LABEL = r"$\bf{sPHENIX}$ Internal"
PP_LABEL = r"$p$+$p$ $\sqrt{s}$ = 200 GeV"

FIGDIR = "/sphenix/user/shuhangli/ppg12/plotting/figures/pileup_vtx"
CSV_PATH = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/pileup_vertex_acceptance.csv"
PT_EDGES = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)


def load_csv(path):
    """Load CSV into dict of numpy arrays."""
    data = {}
    with open(path) as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    for key in rows[0]:
        try:
            data[key] = np.array([float(r[key]) for r in rows])
        except ValueError:
            data[key] = np.array([r[key] for r in rows])
    return data


def pt_centers_widths(data):
    lo = data["pt_lo"]
    hi = data["pt_hi"]
    return 0.5 * (lo + hi), 0.5 * (hi - lo)


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    d = load_csv(CSV_PATH)
    xc, xw = pt_centers_widths(d)

    # =======================================================================
    # Figure 1: Per-pT vertex_eff single vs double
    # =======================================================================
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.errorbar(xc - 0.3, d["single_eff_mbd_vtx"], xerr=xw, yerr=d["single_err"],
                fmt="s", ms=6, color="blue", label="Single interaction")
    ax.errorbar(xc + 0.3, d["double_eff_mbd_vtx"], xerr=xw, yerr=d["double_err"],
                fmt="o", ms=6, color="red", label="Double interaction")

    ax.set_xlabel(r"$E_T^{\gamma,\,\mathrm{truth}}$ [GeV]")
    ax.set_ylabel(r"$\varepsilon_{\mathrm{vtx+MBD}}$ (MBD N\&S$\geq$1 + $|z_\mathrm{vtx}|<60$ cm)")
    ax.set_xlim(7, 37)
    ax.set_ylim(0.45, 0.95)
    ax.legend(loc="lower left")
    ax.text(0.05, 0.95, SPHENIX_LABEL, transform=ax.transAxes, fontsize=14, va="top")
    ax.text(0.05, 0.89, PP_LABEL, transform=ax.transAxes, fontsize=12, va="top")
    ax.text(0.05, 0.83, "photon MC, trkID truth match", transform=ax.transAxes, fontsize=11, va="top")

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "vertex_eff_single_vs_double.pdf"))
    print(f"Saved vertex_eff_single_vs_double.pdf")
    plt.close()

    # =======================================================================
    # Figure 2: Per-pT correction factor
    # =======================================================================
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.errorbar(xc - 0.3, d["corr_0mrad"], xerr=xw, fmt="s", ms=6,
                color="red", label=r"0 mrad ($f_\mathrm{double}$ = 11.1%)")
    ax.errorbar(xc + 0.3, d["corr_1p5mrad"], xerr=xw, fmt="o", ms=6,
                color="blue", label=r"1.5 mrad ($f_\mathrm{double}$ = 3.9%)")

    ax.axhline(1.0, color="gray", ls="--", lw=1)
    ax.set_xlabel(r"$E_T^{\gamma,\,\mathrm{truth}}$ [GeV]")
    ax.set_ylabel(r"Correction factor ($\varepsilon_\mathrm{blended} / \varepsilon_\mathrm{single}$)")
    ax.set_xlim(7, 37)
    ax.set_ylim(0.995, 1.045)
    ax.legend(loc="upper left")
    ax.text(0.95, 0.95, SPHENIX_LABEL, transform=ax.transAxes, fontsize=14, va="top", ha="right")
    ax.text(0.95, 0.89, "Section 3 pileup correction", transform=ax.transAxes, fontsize=11, va="top", ha="right")

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "correction_factor_vs_pt.pdf"))
    print(f"Saved correction_factor_vs_pt.pdf")
    plt.close()

    # =======================================================================
    # Figure 3: Acceptance breakdown (MBD-only, vtx-only, combined)
    # =======================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for ax, tag, label, color in [(axes[0], "single", "Single interaction", "blue"),
                                   (axes[1], "double", "Double interaction", "red")]:
        ax.errorbar(xc, d[f"{tag}_eff_mbd_vtx"], xerr=xw,
                    fmt="s", ms=6, color="black", label="MBD + vertex (combined)")
        ax.errorbar(xc - 0.4, d[f"{tag}_eff_mbd"], xerr=xw,
                    fmt="^", ms=5, color=color, alpha=0.7, label="MBD N&S≥1 only")
        ax.errorbar(xc + 0.4, d[f"{tag}_eff_vtx"], xerr=xw,
                    fmt="v", ms=5, color="green", alpha=0.7, label="|z_vtx| < 60 cm only")

        ax.set_xlabel(r"$E_T^{\gamma,\,\mathrm{truth}}$ [GeV]")
        ax.set_xlim(7, 37)
        ax.set_title(label, fontsize=14, fontweight="bold")
        ax.legend(loc="lower left", fontsize=10)

    axes[0].set_ylabel("Acceptance")
    axes[0].set_ylim(0.45, 1.0)
    axes[0].text(0.05, 0.97, SPHENIX_LABEL, transform=axes[0].transAxes, fontsize=13, va="top")

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "acceptance_breakdown.pdf"))
    print(f"Saved acceptance_breakdown.pdf")
    plt.close()

    # =======================================================================
    # Figure 4: MBD hit decomposition (stacked bars)
    # =======================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    width = 0.7 * (PT_EDGES[1] - PT_EDGES[0])
    # Use actual bin widths for irregular binning
    widths = (PT_EDGES[1:] - PT_EDGES[:-1]) * 0.7

    for ax, tag, label in [(axes[0], "single", "Single interaction"),
                           (axes[1], "double", "Double interaction")]:
        bottom = np.zeros(len(xc))
        categories = [
            (f"{tag}_frac_both", "Both N & S", "#2ecc71"),
            (f"{tag}_frac_N", "N only", "#3498db"),
            (f"{tag}_frac_S", "S only", "#e67e22"),
            (f"{tag}_frac_neither", "Neither", "#e74c3c"),
        ]
        for key, lbl, col in categories:
            vals = d[key]
            ax.bar(xc, vals, width=widths, bottom=bottom, color=col, edgecolor="black",
                   linewidth=0.5, label=lbl, alpha=0.85)
            bottom += vals

        ax.set_xlabel(r"$E_T^{\gamma,\,\mathrm{truth}}$ [GeV]")
        ax.set_xlim(7, 37)
        ax.set_title(label, fontsize=14, fontweight="bold")
        ax.legend(loc="upper right", fontsize=10)

    axes[0].set_ylabel("Fraction of truth photons")
    axes[0].set_ylim(0, 1.05)
    axes[0].text(0.05, 0.97, SPHENIX_LABEL, transform=axes[0].transAxes, fontsize=13, va="top")
    axes[0].text(0.05, 0.91, r"|$z_\mathrm{vtx}^\mathrm{truth}$| < 60 cm", transform=axes[0].transAxes, fontsize=11, va="top")

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "mbd_hit_decomposition.pdf"))
    print(f"Saved mbd_hit_decomposition.pdf")
    plt.close()

    # =======================================================================
    # Figure 5: Summary — single-panel correction + ratio annotation
    # =======================================================================
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), height_ratios=[2, 1],
                                    sharex=True, gridspec_kw={"hspace": 0.05})

    # Top: vertex_eff
    ax1.errorbar(xc - 0.3, d["single_eff_mbd_vtx"], xerr=xw, yerr=d["single_err"],
                 fmt="s", ms=6, color="blue", label="Single interaction")
    ax1.errorbar(xc + 0.3, d["double_eff_mbd_vtx"], xerr=xw, yerr=d["double_err"],
                 fmt="o", ms=6, color="red", label="Double interaction")
    ax1.set_ylabel(r"$\varepsilon_{\mathrm{vtx+MBD}}$")
    ax1.set_ylim(0.45, 0.95)
    ax1.legend(loc="lower left")
    ax1.text(0.05, 0.95, SPHENIX_LABEL, transform=ax1.transAxes, fontsize=14, va="top")
    ax1.text(0.05, 0.88, PP_LABEL, transform=ax1.transAxes, fontsize=12, va="top")

    # Bottom: ratio (double / single)
    ratio = d["double_eff_mbd_vtx"] / d["single_eff_mbd_vtx"]
    ratio_err = ratio * np.sqrt(
        (d["double_err"] / d["double_eff_mbd_vtx"])**2 +
        (d["single_err"] / d["single_eff_mbd_vtx"])**2
    )
    ax2.errorbar(xc, ratio, xerr=xw, yerr=ratio_err, fmt="ko", ms=5)
    ax2.axhline(1.0, color="gray", ls="--", lw=1)
    mean_ratio = np.average(ratio, weights=1.0/ratio_err**2)
    ax2.axhline(mean_ratio, color="red", ls="-", lw=1.5,
                label=f"mean = {mean_ratio:.3f}")
    ax2.set_xlabel(r"$E_T^{\gamma,\,\mathrm{truth}}$ [GeV]")
    ax2.set_ylabel("Double / Single")
    ax2.set_xlim(7, 37)
    ax2.set_ylim(1.05, 1.35)
    ax2.legend(loc="upper left", fontsize=10)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "vertex_eff_ratio.pdf"))
    print(f"Saved vertex_eff_ratio.pdf")
    plt.close()

    print(f"\nAll figures saved to {FIGDIR}")


if __name__ == "__main__":
    main()
