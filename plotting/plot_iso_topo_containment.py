#!/usr/bin/env python3
"""
Plot iso_topo distributions (data vs inclusive MC) and photon containment vs R.

Inputs : efficiencytool/results/iso_topo_containment.npz
Outputs:
  plotting/figures/iso_topo_dist_vs_R_data_vs_mc.pdf    (6-panel, 1 per R; data vs inclusive MC, post-common)
  plotting/figures/iso_topo_overlay_per_sample.pdf      (3-panel overlay of 6 R on truth/inclMC/data)
  plotting/figures/photon_containment_vs_R.pdf          (f_contain vs R, pre/post common)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

NPZ = Path("/sphenix/user/shuhangli/ppg12/efficiencytool/results/iso_topo_containment.npz")
OUT_DIR = Path("/sphenix/user/shuhangli/ppg12/plotting/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)
RADII = [("005", 0.05), ("0075", 0.075), ("01", 0.1), ("02", 0.2), ("03", 0.3), ("04", 0.4)]
R_COLORS = ["#000000", "#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#8c564b"]
NPB_CUT = 0.5
ISO_RANGE = (-20, 30)
ISO_BINS  = 100


def sphenix_tag(ax, label_lines):
    ax.text(0.05, 0.92, label_lines[0], transform=ax.transAxes, fontsize=11, fontweight="bold")
    for i, line in enumerate(label_lines[1:], 1):
        ax.text(0.05, 0.92 - 0.06 * i, line, transform=ax.transAxes, fontsize=9)


def hist_norm(values, weights, bins, range_):
    h, edges = np.histogram(values, bins=bins, range=range_, weights=weights, density=False)
    centers = 0.5 * (edges[:-1] + edges[1:])
    integ = h.sum() * (edges[1] - edges[0])
    if integ > 0:
        h = h / integ
    return centers, h


def main():
    d = np.load(NPZ)
    print(f"loaded {NPZ} keys={len(d.files)}")

    # Figure 1: 6-panel, iso_topo_R distribution for each R; overlay data vs MC-inclusive (post-common)
    fig, axs = plt.subplots(2, 3, figsize=(14, 8), sharey=False)
    axs = axs.flatten()
    for i, (rs, r) in enumerate(RADII):
        ax = axs[i]
        # Data post-common
        d_mask = d["d_npb"] > NPB_CUT
        xc, yd = hist_norm(d[f"d_iso_{rs}"][d_mask], np.ones(d_mask.sum()), ISO_BINS, ISO_RANGE)
        # MC inclusive post-common
        m_mask = d["mc_npb"] > NPB_CUT
        _, ym = hist_norm(d[f"mc_iso_{rs}"][m_mask], d["mc_w"][m_mask], ISO_BINS, ISO_RANGE)
        ax.step(xc, ym, where="mid", label="Inclusive MC", color="#1f77b4", lw=1.5)
        ax.errorbar(xc, yd, fmt="o", ms=3, color="black", label="Data", lw=0.8)
        ax.set_xlabel(r"$iso_{topo}^{R}$ [GeV]")
        ax.set_ylabel("Normalised")
        ax.set_title(f"R = {r}")
        ax.legend(fontsize=8, loc="upper right")
        ax.set_yscale("log")
        ax.set_ylim(1e-5, 1)
    sphenix_tag(axs[0], [r"$\bf{sPHENIX}$ Internal",
                         r"$p+p$ $\sqrt{s}=200$ GeV",
                         r"post-common (npb>0.5)"])
    fig.suptitle("cluster_iso_topo distribution, data vs inclusive MC (post-common)", fontsize=11)
    fig.tight_layout()
    out = OUT_DIR / "iso_topo_dist_vs_R_data_vs_mc.pdf"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")

    # Figure 2: overlay of 6 R curves, one panel per sample (truth, inclMC, data); post-common
    fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    sets = [
        ("MC truth-matched photons", "tr_iso_{}", "tr_w", "tr_npb"),
        ("Inclusive MC",             "mc_iso_{}", "mc_w", "mc_npb"),
        ("Data",                     "d_iso_{}",  None,   "d_npb"),
    ]
    for ax, (name, pat, wk, npbk) in zip(axs, sets):
        mask = d[npbk] > NPB_CUT
        for (rs, r), color in zip(RADII, R_COLORS):
            key = pat.format(rs)
            w = d[wk][mask] if wk else np.ones(mask.sum())
            xc, y = hist_norm(d[key][mask], w, ISO_BINS, ISO_RANGE)
            ax.step(xc, y, where="mid", color=color, lw=1.5, label=f"R = {r}")
        ax.set_xlabel(r"$iso_{topo}$ [GeV]")
        ax.set_title(f"{name}  [post-common]")
        ax.set_yscale("log"); ax.set_ylim(1e-5, 1)
    axs[0].set_ylabel("Normalised")
    axs[0].legend(fontsize=8, loc="upper right", ncol=2)
    sphenix_tag(axs[0], [r"$\bf{sPHENIX}$ Internal"])
    fig.tight_layout()
    out = OUT_DIR / "iso_topo_overlay_per_sample.pdf"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")

    # Figure 3: photon containment vs R (MC truth-matched signal)
    fig, ax = plt.subplots(figsize=(7, 5))
    pt_pre  = d["tr_part_Pt"]
    clu_pre = d["tr_clu_Et"]
    w_pre   = d["tr_w"]
    for state_lbl, mask, color in [
        ("pre-common",                            np.ones(len(pt_pre), bool), "#1f77b4"),
        (f"post-common (npb > {NPB_CUT})",        d["tr_npb"] > NPB_CUT,      "#d62728"),
    ]:
        rs_x, fmeans, flow, fhigh = [], [], [], []
        for rs, r in RADII:
            vals = (d[f"tr_iso_{rs}"][mask] + clu_pre[mask]) / pt_pre[mask]
            ww   = w_pre[mask]
            mu   = np.average(vals, weights=ww)
            rms  = np.sqrt(np.average((vals - mu) ** 2, weights=ww))
            rs_x.append(r); fmeans.append(mu)
            flow.append(mu - rms); fhigh.append(mu + rms)
        ax.fill_between(rs_x, flow, fhigh, color=color, alpha=0.2)
        ax.plot(rs_x, fmeans, "o-", color=color, lw=1.5, label=state_lbl)
    ax.axhline(1.0, ls="--", color="gray", lw=1)
    ax.set_xlabel("Inner cone radius R")
    ax.set_ylabel(r"$\langle f_{\rm contain}(R) \rangle = \langle (iso_{topo}^R + E_T^{\rm cluster}) / p_T^{\gamma, truth}\rangle$")
    ax.set_ylim(0.0, 1.2); ax.set_xlim(0, 0.45)
    sphenix_tag(ax, [r"$\bf{sPHENIX}$ Internal",
                     r"$p+p$ $\sqrt{s}=200$ GeV",
                     "MC truth-matched direct photons"])
    ax.legend(fontsize=9, loc="lower right")
    fig.tight_layout()
    out = OUT_DIR / "photon_containment_vs_R.pdf"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
