#!/usr/bin/env python3
"""
plot_lumi_convention_validation.py — Plots for the Section 3 validation report.

Creates:
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


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    d = np.load("/sphenix/user/shuhangli/ppg12/efficiencytool/results/lumi_convention_validation.npz")

    z = d["z_centers"]
    old_ratio = d["old_ratio"]
    new_ratio_1cm = d["new_ratio_1cm"]
    h_data = d["h_data_vtx"]
    h_sim = d["h_sim_vtx"]
    h_sim_rw_old = d["h_sim_rw_old"]
    h_sim_rw_new = d["h_sim_rw_new"]
    pt_edges = d["pt_edges"]
    eff_old = d["eff_old"]
    eff_new = d["eff_new"]
    ratio = d["ratio"]

    L60 = 16.8790
    Lnv = 17.1642

    # =======================================================================
    # Figure 1: Vertex reweight shape comparison
    # =======================================================================
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(z, old_ratio, "b-", lw=1, alpha=0.5, label="OLD (1 cm bins, fallback=1.0)")
    ax.plot(z, new_ratio_1cm, "r-", lw=2, label="NEW (5 cm rebin + smooth + fallback=0.0)")
    ax.axhline(1.0, color="gray", ls="--", lw=0.8, alpha=0.5)
    ax.set_xlabel("$z_\\mathrm{vtx}^\\mathrm{reco}$ [cm]")
    ax.set_ylabel("Reweight (data / sim, normalized)")
    ax.set_title("Vertex reweight shape: OLD vs NEW", fontsize=13)
    ax.set_xlim(-100, 100)
    ax.set_ylim(0, 4.5)
    ax.legend(loc="upper right")
    ax.text(0.02, 0.97, SPHENIX_LABEL, transform=ax.transAxes, fontsize=13, va="top")
    ax.text(0.02, 0.91, PP_LABEL, transform=ax.transAxes, fontsize=11, va="top")
    ax.text(0.02, 0.85, "photon10 MC (1.5 mrad)", transform=ax.transAxes, fontsize=10, va="top")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "vertex_reweight_shape.pdf"))
    print(f"Saved vertex_reweight_shape.pdf")
    plt.close()

    # =======================================================================
    # Figure 2: Reweighted MC vertex distribution vs data
    # =======================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: linear scale
    ax = axes[0]
    data_norm = h_data / h_data.sum()
    rw_old_norm = h_sim_rw_old / h_sim_rw_old.sum()
    rw_new_norm = h_sim_rw_new / h_sim_rw_new.sum()
    raw_norm = h_sim / h_sim.sum()
    ax.step(z, data_norm, "k-", lw=1.5, where="mid", label="Data")
    ax.step(z, raw_norm, "g--", lw=1, where="mid", alpha=0.7, label="Raw MC (no reweight)")
    ax.step(z, rw_old_norm, "b-", lw=1, where="mid", alpha=0.8, label="OLD reweighted MC")
    ax.step(z, rw_new_norm, "r-", lw=1.5, where="mid", label="NEW reweighted MC")
    ax.axvline(60, color="gray", ls=":", lw=0.8)
    ax.axvline(-60, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("$z_\\mathrm{vtx}^\\mathrm{reco}$ [cm]")
    ax.set_ylabel("Normalized")
    ax.set_xlim(-100, 100)
    ax.set_title("Vertex distribution (linear)", fontsize=13)
    ax.legend(loc="upper right", fontsize=10)
    ax.text(0.02, 0.97, SPHENIX_LABEL, transform=ax.transAxes, fontsize=12, va="top")

    # Right: log scale
    ax = axes[1]
    ax.step(z, data_norm, "k-", lw=1.5, where="mid", label="Data")
    ax.step(z, raw_norm, "g--", lw=1, where="mid", alpha=0.7, label="Raw MC")
    ax.step(z, rw_old_norm, "b-", lw=1, where="mid", alpha=0.8, label="OLD reweighted MC")
    ax.step(z, rw_new_norm, "r-", lw=1.5, where="mid", label="NEW reweighted MC")
    ax.axvline(60, color="gray", ls=":", lw=0.8)
    ax.axvline(-60, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("$z_\\mathrm{vtx}^\\mathrm{reco}$ [cm]")
    ax.set_ylabel("Normalized")
    ax.set_yscale("log")
    ax.set_xlim(-100, 100)
    ax.set_ylim(1e-6, 1e-1)
    ax.set_title("Vertex distribution (log)", fontsize=13)
    ax.legend(loc="lower center", fontsize=10)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "vertex_dist_comparison.pdf"))
    print(f"Saved vertex_dist_comparison.pdf")
    plt.close()

    # =======================================================================
    # Figure 3: Per-pT L×eff ratio (new/old)
    # =======================================================================
    pt_centers = 0.5 * (pt_edges[:-1] + pt_edges[1:])
    pt_widths = 0.5 * (pt_edges[1:] - pt_edges[:-1])
    # Extension bins may be missing (n_bins may differ)
    n_valid = len(ratio)
    if n_valid == 11:
        pt_c_plot = np.array([9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 34])
        pt_w_plot = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2])
    else:
        pt_c_plot = pt_centers[:n_valid]
        pt_w_plot = pt_widths[:n_valid]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                                     gridspec_kw={"height_ratios": [2, 1]})
    # Top: L×eff for both conventions
    ax1.errorbar(pt_c_plot - 0.15, L60 * eff_old, xerr=pt_w_plot,
                 fmt="bs", ms=7, label=r"OLD: $\mathcal{L}_{60\mathrm{cm}} \times \varepsilon_{\mathrm{old}}$")
    ax1.errorbar(pt_c_plot + 0.15, Lnv * eff_new, xerr=pt_w_plot,
                 fmt="ro", ms=7, label=r"NEW: $\mathcal{L}_{\mathrm{noVtx}} \times \varepsilon_{\mathrm{new}}$")
    ax1.set_ylabel(r"$\mathcal{L} \times \varepsilon_\mathrm{vtx+MBD}$  [pb$^{-1}$]")
    ax1.set_title("L × eff product per pT bin (1.5 mrad)", fontsize=13)
    ax1.legend(loc="upper right")
    ax1.set_ylim(6, 13)
    ax1.text(0.02, 0.97, SPHENIX_LABEL, transform=ax1.transAxes, fontsize=13, va="top")
    ax1.text(0.02, 0.90, PP_LABEL, transform=ax1.transAxes, fontsize=11, va="top")

    # Bottom: ratio
    ax2.errorbar(pt_c_plot, ratio, xerr=pt_w_plot, fmt="ko", ms=7)
    ax2.axhline(1.0, color="gray", ls="--", lw=0.8)
    mean_ratio = ratio.mean()
    ax2.axhline(mean_ratio, color="red", ls="-", lw=1.5, alpha=0.6,
                label=f"mean = {mean_ratio:.4f} ({(mean_ratio-1)*100:+.2f}%)")
    ax2.set_xlabel(r"$p_T^{\gamma,\mathrm{truth}}$ [GeV]")
    ax2.set_ylabel("NEW / OLD")
    ax2.set_ylim(0.99, 1.02)
    ax2.legend(loc="upper right", fontsize=10)
    ax2.set_xlim(7, 37)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "L_times_eff_ratio.pdf"))
    print(f"Saved L_times_eff_ratio.pdf")
    plt.close()

    # =======================================================================
    # Figure 4: Rebin stability — test different rebin factors
    # =======================================================================
    def compute_and_integrate(rebin_factor, smooth_n):
        h_d = h_data.copy()
        h_s = h_sim.copy()
        n = len(h_d)
        new_n = n // rebin_factor
        h_d_rb = h_d[:new_n*rebin_factor].reshape(new_n, rebin_factor).sum(axis=1)
        h_s_rb = h_s[:new_n*rebin_factor].reshape(new_n, rebin_factor).sum(axis=1)
        h_d_rb = h_d_rb / h_d_rb.sum() if h_d_rb.sum() > 0 else h_d_rb
        h_s_rb = h_s_rb / h_s_rb.sum() if h_s_rb.sum() > 0 else h_s_rb
        with np.errstate(divide='ignore', invalid='ignore'):
            r = h_d_rb / h_s_rb
        r[~np.isfinite(r)] = 0
        r[r < 0] = 0
        for _ in range(smooth_n):
            smooth = np.copy(r)
            for i in range(1, len(r)-1):
                smooth[i] = (r[i-1] + 2*r[i] + r[i+1])/4
            r = smooth
        return r, new_n*rebin_factor

    fig, ax = plt.subplots(figsize=(9, 6))
    rebin_factors = [1, 2, 5, 10, 20]
    colors = ["C0", "C1", "C2", "C3", "C4"]
    for rf, col in zip(rebin_factors, colors):
        r, n_used = compute_and_integrate(rf, 1)
        # Bin width = 1 cm * rf
        edges_rb = np.linspace(-100, -100 + len(r)*rf, len(r)+1)
        centers_rb = 0.5*(edges_rb[:-1]+edges_rb[1:])
        ax.step(centers_rb, r, color=col, where="mid", lw=1.5 if rf==5 else 1,
                label=f"rebin={rf} ({rf} cm bins)")
    ax.axhline(1.0, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("$z_\\mathrm{vtx}^\\mathrm{reco}$ [cm]")
    ax.set_ylabel("Reweight (data/sim, normalized)")
    ax.set_title("Rebin stability of the fixed reweight", fontsize=13)
    ax.set_xlim(-100, 100)
    ax.set_ylim(0, 4.5)
    ax.legend(loc="upper right", ncol=2)
    ax.text(0.02, 0.97, SPHENIX_LABEL, transform=ax.transAxes, fontsize=13, va="top")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "rebin_stability.pdf"))
    print(f"Saved rebin_stability.pdf")
    plt.close()

    print(f"\nAll plots in {FIGDIR}")


if __name__ == "__main__":
    main()
