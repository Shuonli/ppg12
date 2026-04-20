#!/usr/bin/env python3
# Pipeline-side cross-check overlay for PPG12 cross-section weights.
#
# Reads the existing per-sample MC_efficiency_{sample}_bdt_1p5rad.root files
# (produced by RecoEffCalculator_TTreeReader.C with cross_weight * vertex_weight
# already applied and the per-sample truth-pT window already cut) and overlays
# them + the combined sum. Any non-smooth transition at a stitching boundary
# (14 GeV / 30 GeV for photons; 9/14/21/32/42 GeV for jets) indicates a
# weight mismatch between samples.
#
# Complementary to TruthSpectrumOverlay.C which reads raw trees.

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

RESULTS = "/sphenix/user/shuhangli/ppg12/efficiencytool/results"
FIGDIR = "/sphenix/user/shuhangli/ppg12/plotting/figures/truth_spectrum"
VAR_TYPE = "bdt_1p5rad"

os.makedirs(FIGDIR, exist_ok=True)

# (sample, color, label, stitching_boundary_below)
PHOTON_SAMPLES = [
    ("photon5",  "tab:red",    "photon5 [0,14] GeV",  None),
    ("photon10", "tab:blue",   "photon10 [14,30] GeV", 14.0),
    ("photon20", "tab:green",  "photon20 [30,+) GeV",  30.0),
]

JET_SAMPLES = [
    ("jet5",  "tab:red",     "jet5 [7,9] GeV",    None),
    ("jet8",  "tab:orange",  "jet8 [9,14] GeV",   9.0),
    ("jet12", "tab:blue",    "jet12 [14,21] GeV", 14.0),
    ("jet20", "tab:green",   "jet20 [21,32] GeV", 21.0),
    ("jet30", "tab:purple",  "jet30 [32,42] GeV", 32.0),
    ("jet40", "tab:brown",   "jet40 [42,+) GeV",  42.0),
]


def load_hist(path, hname):
    f = uproot.open(path)
    if hname not in f:
        return None
    h = f[hname]
    vals = h.values()
    errs = h.errors()
    edges = h.axes[0].edges()
    centers = 0.5 * (edges[:-1] + edges[1:])
    bw = edges[1] - edges[0]
    return centers, vals, errs, bw, edges


def draw_overlay(samples, hname, title, outpdf, xlim, ylim, boundaries):
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("truth p_{T} [GeV]")
    ax.set_ylabel("weighted entries (xsec_rel * vertex_reweight)")
    ax.set_title(title)

    total_vals = None
    total_bw = None
    total_centers = None

    for (s, color, label, bnd) in samples:
        path = f"{RESULTS}/MC_efficiency_{s}_{VAR_TYPE}.root"
        if not os.path.exists(path):
            print(f"  [miss] {path}")
            continue
        res = load_hist(path, hname)
        if res is None:
            print(f"  [miss hist {hname}] {path}")
            continue
        centers, vals, errs, bw, edges = res
        # Mask zero entries for log plot
        nonzero = vals > 0
        if nonzero.any():
            ax.plot(centers[nonzero], vals[nonzero], color=color,
                    linewidth=2.5, label=label, drawstyle="steps-mid",
                    zorder=5)
        integral = vals.sum()
        print(f"  {s:14s} [{hname}]: integral={integral:.3e}, nonzero={nonzero.sum()}")

        # Initialize / accumulate total (must share binning)
        if total_vals is None:
            total_vals = vals.copy()
            total_centers = centers
            total_bw = bw
        else:
            if len(vals) != len(total_vals):
                print(f"    WARNING: {s} has different binning; skipping from sum")
                continue
            total_vals = total_vals + vals

    if total_vals is not None:
        nonzero = total_vals > 0
        ax.plot(total_centers[nonzero], total_vals[nonzero],
                color="black", linewidth=1.5, linestyle="--",
                marker="o", markersize=3, markerfacecolor="none",
                label="combined sum", alpha=0.8)

    # Draw stitching boundaries
    for b in boundaries:
        ax.axvline(b, color="gray", linestyle="--", alpha=0.5)

    # Ratio of content just above/below each boundary to detect kinks
    print(f"  boundary smoothness check ({title}):")
    if total_vals is not None:
        for b in boundaries:
            # Find bin containing b-0.5 and b+0.5 (bins of width ~1 GeV)
            below_idx = np.argmin(np.abs(total_centers - (b - 0.5)))
            above_idx = np.argmin(np.abs(total_centers - (b + 0.5)))
            if 0 <= below_idx < len(total_vals) and 0 <= above_idx < len(total_vals):
                b_val = total_vals[below_idx]
                a_val = total_vals[above_idx]
                ratio = b_val / a_val if a_val > 0 else 0.0
                flag = "  <-- KINK" if (ratio < 0.5 or ratio > 2.0) else ""
                print(f"    pT={b:.0f} GeV: below={b_val:.3e} above={a_val:.3e} ratio={ratio:.3f}{flag}")

    ax.legend(loc="best", frameon=False, fontsize=10)
    ax.text(0.03, 0.95, "PPG12 pipeline cross-check (weighted by cross_weight * vertex_reweight)",
            transform=ax.transAxes, fontsize=9, verticalalignment="top", style="italic")
    ax.text(0.03, 0.90, f"var_type: {VAR_TYPE}",
            transform=ax.transAxes, fontsize=9, verticalalignment="top")

    fig.tight_layout()
    fig.savefig(outpdf)
    plt.close(fig)
    print(f"  wrote {outpdf}\n")


def main():
    print("=== Photon main pipeline ===")
    draw_overlay(PHOTON_SAMPLES,
                 "h_truth_pT_0",
                 "Truth photon pT per sample + combined (main pipeline, bdt_1p5rad)",
                 f"{FIGDIR}/pipeline_truth_spectrum_photon.pdf",
                 xlim=(5, 45), ylim=(1e2, 1e12),
                 boundaries=[14.0, 30.0])

    print("=== Jet main pipeline ===")
    draw_overlay(JET_SAMPLES,
                 "h_max_truth_jet_pT",
                 "Max truth jet pT per event per sample + combined (main pipeline, bdt_1p5rad)",
                 f"{FIGDIR}/pipeline_truth_spectrum_jet.pdf",
                 xlim=(5, 60), ylim=(1e2, 1e14),
                 boundaries=[9.0, 14.0, 21.0, 32.0, 42.0])


if __name__ == "__main__":
    main()
