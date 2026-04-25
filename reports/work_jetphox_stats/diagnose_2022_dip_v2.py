#!/usr/bin/env python3
"""
Diagnose the [20,22] GeV dip using the 2D histogram h_truth_eta_pT that
MakeJetPHOXhisto.C already writes (100 eta bins [-1,1] x 100 pT bins [0,100]).
Fine pT binning is 1 GeV per bin.

This avoids re-reading the 1e8-entry ntuple.
"""
import os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT_DIR = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/NLO/rootFiles"
OUT_DIR  = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/work_jetphox_stats"

PDFs = {
    "CT14lo":  "jetPHOX",
    "CT14nlo": "jetPHOX_nlo",
    "NNPDF":   "jetPHOX_nnpdf",
}
SCALES = ["05", "10", "20"]

def project_to_y07(h2d):
    """Project h_truth_eta_pT 2D onto 1D pT restricted to |y|<0.7.
       The 2D has y axis = 100 bins from -1 to +1. Cells with center |y| < 0.7
       are kept. Output is 1D in pT."""
    edges_y = h2d.axis(0).edges()
    edges_pt = h2d.axis(1).edges()
    vals = h2d.values()   # shape (ny, npt)
    centers_y = 0.5 * (edges_y[:-1] + edges_y[1:])
    ymask = np.abs(centers_y) < 0.7
    pt_sum = vals[ymask].sum(axis=0)
    return edges_pt, pt_sum

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    fig, (ax_main, ax_res) = plt.subplots(
        2, 1, figsize=(9, 8), sharex=True,
        gridspec_kw={"height_ratios": [3, 2]},
    )
    colors = {"CT14lo": "C0", "CT14nlo": "C1", "NNPDF": "C2"}

    print("=== 1-GeV binned pT distributions (direct+frag combined), mu=1.0 ===\n")
    print(f"{'pT':>6} | {'CT14lo':>10} | {'CT14nlo':>10} | {'NNPDF':>10}")
    print("-"*50)
    store = {}
    for pdf, stem in PDFs.items():
        path = f"{ROOT_DIR}/{stem}_10.root"
        with uproot.open(path) as f:
            h2 = f["h_truth_eta_pT"]
            edges_pt, pt_sum = project_to_y07(h2)
        store[pdf] = (edges_pt, pt_sum)

    edges_pt = store["CT14lo"][0]
    centers = 0.5 * (edges_pt[:-1] + edges_pt[1:])

    # Print values in the interesting range
    for i, c in enumerate(centers):
        if c < 7 or c > 37:
            continue
        row = f"{c:6.1f} |"
        for pdf in ["CT14lo", "CT14nlo", "NNPDF"]:
            row += f" {store[pdf][1][i]:10.4g} |"
        print(row)

    # Plot raw 1-GeV fine binning (arbitrary units — raw weighted counts, not dsigma/dpT yet)
    # For the smoothness test the absolute normalization doesn't matter.
    for pdf, (_, vals) in store.items():
        mask = (centers > 7) & (centers < 38)
        ax_main.step(centers[mask], vals[mask], where="mid", color=colors[pdf], label=pdf)

    ax_main.set_yscale("log")
    ax_main.set_ylabel("raw sum of pdf_weight (arbitrary units)")
    ax_main.set_title("1-GeV fine binning from h_truth_eta_pT (|y|<0.7, mu=1.0)")
    ax_main.axvspan(20, 22, color="red", alpha=0.15, label="[20,22] PPG12 bin")
    ax_main.legend()
    ax_main.grid(True, alpha=0.3)

    # Residual from smooth power-law
    for pdf, (_, vals) in store.items():
        mask_fit = (centers > 10) & (centers < 35) & (vals > 0)
        slope, intercept = np.polyfit(np.log(centers[mask_fit]), np.log(vals[mask_fit]), 1)
        smooth = np.exp(intercept) * centers ** slope
        print(f"  {pdf}: fitted power-law slope = pT^{slope:.2f}")
        ax_res.step(centers[mask_fit], vals[mask_fit] / smooth[mask_fit], where="mid", color=colors[pdf], label=f"{pdf}")

    ax_res.axhline(1.0, color="k", lw=0.5)
    ax_res.axvspan(20, 22, color="red", alpha=0.15)
    ax_res.set_xlabel(r"$p_T$ (GeV)")
    ax_res.set_ylabel("ratio to power-law fit (10-35 GeV)")
    ax_res.set_ylim(0.5, 1.8)
    ax_res.grid(True, alpha=0.3)
    ax_res.legend()

    plt.tight_layout()
    path = f"{OUT_DIR}/diagnose_2022_dip.pdf"
    plt.savefig(path)
    plt.savefig(path.replace(".pdf", ".png"), dpi=150)
    print(f"\nPlot: {path}")

    # Also print the 2-GeV-binned version of the same data for comparison
    # Rebin pairs of 1-GeV bins -> 2-GeV bins
    print("\n=== Rebinned to PPG12 2-GeV bins (same data) ===\n")
    # Keep only pT ∈ [8,36]
    keep = (centers >= 8) & (centers <= 35.5)
    idx = np.where(keep)[0]
    # group by pairs (8-9 with 9-10 -> 8-10, 10-11 with 11-12 -> 10-12, ...)
    # Our 1-GeV centers are 0.5, 1.5, ..., 99.5.  Bin [8,10] corresponds to centers 8.5 and 9.5.
    print(f"{'pT bin':>10} | {'CT14lo':>10} | {'CT14nlo':>10} | {'NNPDF':>10}")
    print("-"*60)
    for lo in [8, 10, 12, 14, 16, 18, 20, 22, 24, 26]:
        hi = lo + 2
        bmask = (centers >= lo - 0.01) & (centers < hi - 0.01)
        row = f"{lo:3d}-{hi:3d}    |"
        for pdf in ["CT14lo", "CT14nlo", "NNPDF"]:
            v = store[pdf][1][bmask].sum()
            row += f" {v:10.4g} |"
        print(row)

if __name__ == "__main__":
    main()
