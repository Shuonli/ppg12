#!/usr/bin/env python3
"""
plot_lumi_convention.py — Plots for the lumi convention validation report.

Creates per-period plots for both 1.5 mrad and 0 mrad:
1. Vertex reweight shape: old vs new
2. Reweighted MC vertex distribution vs data
3. Per-pT L×eff ratio (new/old)
4. Rebin stability check
"""

import os
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


def get_period_data(d, period):
    prefix = f"{period}_"
    return {k[len(prefix):]: d[k] for k in d.files if k.startswith(prefix)}


def plot_reweight_shape_combined(d):
    """Figure: vertex reweight shape for both periods side-by-side."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for ax, period in zip(axes, ["1p5mrad", "0mrad"]):
        p = get_period_data(d, period)
        z = p["z_centers"]
        ax.plot(z, p["old_ratio"], "b-", lw=1, alpha=0.5,
                label="OLD (1 cm, fallback=1.0)")
        ax.plot(z, p["new_ratio_1cm"], "r-", lw=2,
                label="NEW (5 cm rebin + smooth + fallback=0.0)")
        ax.axhline(1.0, color="gray", ls="--", lw=0.8, alpha=0.5)
        ax.set_xlabel(r"$z_\mathrm{vtx}^\mathrm{reco}$ [cm]")
        ax.set_xlim(-100, 100)
        ax.set_ylim(0, 4.5)
        ax.set_title(f"{LUMI[period]['label']}", fontsize=13, fontweight="bold")
        ax.legend(loc="upper right", fontsize=10)
    axes[0].set_ylabel("Reweight (data / sim, normalized)")
    axes[0].text(0.02, 0.97, SPHENIX_LABEL, transform=axes[0].transAxes, fontsize=13, va="top")
    axes[0].text(0.02, 0.91, PP_LABEL, transform=axes[0].transAxes, fontsize=11, va="top")
    axes[0].text(0.02, 0.85, "photon10 MC", transform=axes[0].transAxes, fontsize=10, va="top")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "vertex_reweight_shape.pdf"))
    print("Saved vertex_reweight_shape.pdf")
    plt.close()


def plot_vertex_dist_combined(d):
    """Figure: reweighted MC vs data for both periods, linear + log."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for col, period in enumerate(["1p5mrad", "0mrad"]):
        p = get_period_data(d, period)
        z = p["z_centers"]
        h_data_n = p["h_data_vtx"] / p["h_data_vtx"].sum()
        h_raw_n = p["h_sim_vtx"] / p["h_sim_vtx"].sum()
        h_old_n = p["h_sim_rw_old"] / p["h_sim_rw_old"].sum()
        h_new_n = p["h_sim_rw_new"] / p["h_sim_rw_new"].sum()

        for row, scale in enumerate(["linear", "log"]):
            ax = axes[row, col]
            ax.step(z, h_data_n, "k-", lw=1.5, where="mid", label="Data")
            ax.step(z, h_raw_n, "g--", lw=1, where="mid", alpha=0.7, label="Raw MC")
            ax.step(z, h_old_n, "b-", lw=1, where="mid", alpha=0.8, label="OLD reweighted")
            ax.step(z, h_new_n, "r-", lw=1.5, where="mid", label="NEW reweighted")
            ax.axvline(60, color="gray", ls=":", lw=0.8)
            ax.axvline(-60, color="gray", ls=":", lw=0.8)
            ax.set_xlim(-100, 100)
            ax.set_xlabel(r"$z_\mathrm{vtx}^\mathrm{reco}$ [cm]")
            if scale == "log":
                ax.set_yscale("log")
                ax.set_ylim(1e-6, 1e-1)
            if col == 0:
                ax.set_ylabel("Normalized")
            if row == 0:
                ax.set_title(f"{LUMI[period]['label']}", fontsize=13, fontweight="bold")
                ax.legend(loc="upper right", fontsize=9)

    axes[0, 0].text(0.02, 0.97, SPHENIX_LABEL, transform=axes[0, 0].transAxes, fontsize=12, va="top")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "vertex_dist_comparison.pdf"))
    print("Saved vertex_dist_comparison.pdf")
    plt.close()


def plot_L_times_eff_ratio_combined(d):
    """Figure: per-pT L×eff for both periods side-by-side with ratio panel."""
    fig = plt.figure(figsize=(14, 8))
    for col, period in enumerate(["1p5mrad", "0mrad"]):
        p = get_period_data(d, period)
        pt_lo = p["pt_lo"]
        pt_hi = p["pt_hi"]
        pt_c = 0.5 * (pt_lo + pt_hi)
        pt_w = 0.5 * (pt_hi - pt_lo)
        L60 = LUMI[period]["L60"]
        Lnv = LUMI[period]["Lnv"]
        eff_old = p["eff_old"]
        eff_new = p["eff_new"]
        ratio = p["ratio"]
        label = LUMI[period]["label"]

        # Top: L*eff
        ax1 = plt.subplot2grid((3, 2), (0, col), rowspan=2)
        ax1.errorbar(pt_c - 0.15, L60 * eff_old, xerr=pt_w, fmt="bs", ms=7,
                     label=r"OLD: $L_{60\mathrm{cm}} \times \varepsilon_{\mathrm{old}}$")
        ax1.errorbar(pt_c + 0.15, Lnv * eff_new, xerr=pt_w, fmt="ro", ms=7,
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
        # Adjust y range based on period
        if period == "1p5mrad":
            ax2.set_ylim(0.99, 1.02)
        else:
            ax2.set_ylim(1.08, 1.18)
        ax2.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "L_times_eff_ratio.pdf"))
    print("Saved L_times_eff_ratio.pdf")
    plt.close()


def plot_rebin_stability_combined(d):
    """Figure: rebin stability for both periods side-by-side."""

    def compute_reweight(h_data, h_sim, rebin, smooth_n=1):
        n = len(h_data)
        new_n = n // rebin
        h_d = h_data[:new_n*rebin].reshape(new_n, rebin).sum(axis=1)
        h_s = h_sim[:new_n*rebin].reshape(new_n, rebin).sum(axis=1)
        h_d = h_d / h_d.sum() if h_d.sum() > 0 else h_d
        h_s = h_s / h_s.sum() if h_s.sum() > 0 else h_s
        with np.errstate(divide='ignore', invalid='ignore'):
            r = h_d / h_s
        r[~np.isfinite(r)] = 0
        r[r < 0] = 0
        for _ in range(smooth_n):
            sm = np.copy(r)
            for i in range(1, len(r)-1):
                sm[i] = (r[i-1] + 2*r[i] + r[i+1]) / 4
            r = sm
        return r

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    rebin_factors = [1, 2, 5, 10, 20]
    colors = ["C0", "C1", "C2", "C3", "C4"]
    for ax, period in zip(axes, ["1p5mrad", "0mrad"]):
        p = get_period_data(d, period)
        h_data = p["h_data_vtx"]
        h_sim = p["h_sim_vtx"]
        for rf, col in zip(rebin_factors, colors):
            r = compute_reweight(h_data, h_sim, rf, 1)
            edges = np.linspace(-100, -100 + len(r) * rf, len(r) + 1)
            centers = 0.5 * (edges[:-1] + edges[1:])
            ax.step(centers, r, color=col, where="mid",
                    lw=1.5 if rf == 5 else 1,
                    label=f"rebin={rf} ({rf} cm)")
        ax.axhline(1.0, color="gray", ls="--", lw=0.8)
        ax.set_xlabel(r"$z_\mathrm{vtx}^\mathrm{reco}$ [cm]")
        ax.set_xlim(-100, 100)
        ax.set_ylim(0, 4.5)
        ax.set_title(f"{LUMI[period]['label']}", fontsize=13, fontweight="bold")
        ax.legend(loc="upper right", ncol=2, fontsize=10)
    axes[0].set_ylabel("Reweight (data / sim, normalized)")
    axes[0].text(0.02, 0.97, SPHENIX_LABEL, transform=axes[0].transAxes, fontsize=12, va="top")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "rebin_stability.pdf"))
    print("Saved rebin_stability.pdf")
    plt.close()


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    d = np.load("/sphenix/user/shuhangli/ppg12/efficiencytool/results/lumi_convention_validation.npz")

    print("Available periods:", set(k.split("_")[0] for k in d.files))

    plot_reweight_shape_combined(d)
    plot_vertex_dist_combined(d)
    plot_L_times_eff_ratio_combined(d)
    plot_rebin_stability_combined(d)
    print(f"\nAll plots in {FIGDIR}")


if __name__ == "__main__":
    main()
