#!/usr/bin/env python3
"""
3-PDF JETPHOX comparison for PPG12 isolated-photon cross-section.

Reads 9 histograms (3 PDFs × 3 scales) and produces:
  * CSV table of per-bin d^2 sigma/dpT/dy for all 9 combinations
  * Ratio table: CT14nlo/CT14lo, NNPDF/CT14lo
  * PDF ratio plot for the final report

PDFs:
  * CT14lo  — existing Shuonli production (jetPHOX_{05,10,20}.root)
  * CT14nlo — new NLO-PDF production (jetPHOX_nlo_{05,10,20}.root)
  * NNPDF31_nlo_as_0118 — new NLO MC-replica production (jetPHOX_nnpdf_{05,10,20}.root)
"""

import os
import csv
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NLO_DIR = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/NLO/rootFiles"
OUT_DIR = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/work_jetphox_stats"
DETA = 1.4  # |eta|<0.7 -> 2 * 0.7

# Production analysis bins (quoted in the final plot): 8 - 36 GeV
# Truth binning includes overflow: [7,8,10,12,14,16,18,20,22,24,26,28,32,36,45]
# Index 0 = [7,8] (underflow, always 0 in NLO), 13 = [36,45] (overflow)
# Index 1..12 are the 12 PPG12 quoted bins.

PDFs = {
    "CT14lo":  "jetPHOX",
    "CT14nlo": "jetPHOX_nlo",
    "NNPDF":   "jetPHOX_nnpdf",
}
SCALES = ["05", "10", "20"]    # mu = 0.5, 1.0, 2.0 * pT
SCALE_LABEL = {"05": "mu=0.5", "10": "mu=1.0", "20": "mu=2.0"}

def load_all():
    """Return: data[pdf][scale] = (edges, values, errors)."""
    data = {}
    for pdf_label, fname_stem in PDFs.items():
        data[pdf_label] = {}
        for scale in SCALES:
            path = f"{NLO_DIR}/{fname_stem}_{scale}.root"
            with uproot.open(path) as f:
                h = f["h_truth_pT"]
                edges = h.axis().edges()
                values = h.values() / DETA
                errors = h.errors() / DETA
                data[pdf_label][scale] = (edges, values, errors)
    return data

def print_header(cols, widths):
    row = "| " + " | ".join(f"{c:^{widths[i]}}" for i, c in enumerate(cols)) + " |"
    sep = "|" + "|".join("-" * (w+2) for w in widths) + "|"
    print(row)
    print(sep)

def print_row(vals, widths, fmts=None):
    cells = []
    for i, v in enumerate(vals):
        if fmts and fmts[i]:
            cells.append(f"{v:{fmts[i]}}".rjust(widths[i]))
        else:
            cells.append(f"{v}".rjust(widths[i]))
    print("| " + " | ".join(cells) + " |")

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    data = load_all()

    # Reference: CT14lo at mu=1.0 to define the bin grid
    edges_ref, _, _ = data["CT14lo"]["10"]
    nbins = len(edges_ref) - 1

    # Build per-bin table
    rows = []
    for b in range(nbins):
        lo, hi = edges_ref[b], edges_ref[b+1]
        row = {"pt_lo": lo, "pt_hi": hi}
        for pdf in PDFs:
            for scale in SCALES:
                _, vals, errs = data[pdf][scale]
                key = f"{pdf}_{scale}"
                row[key] = vals[b]
                row[key + "_err"] = errs[b]
        # Ratios (central mu=1.0 for both, vs CT14lo mu=1.0 central)
        ct14lo = row["CT14lo_10"] or float("nan")
        row["ratio_CT14nlo_CT14lo"] = row["CT14nlo_10"] / ct14lo if ct14lo else float("nan")
        row["ratio_NNPDF_CT14lo"]   = row["NNPDF_10"]   / ct14lo if ct14lo else float("nan")
        rows.append(row)

    # Write CSV
    csv_path = f"{OUT_DIR}/compare_3pdf_per_bin.csv"
    keys = list(rows[0].keys())
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"Wrote {csv_path}")

    # Pretty table (analysis bins only: 8 - 36 GeV, indices 1..12)
    print("\n--- d^2 sigma / dpT dy  (pb / GeV) at mu=1.0, |y|<0.7 ---")
    cols = ["pT bin", "CT14lo", "CT14nlo", "NNPDF", "CT14nlo/CT14lo", "NNPDF/CT14lo"]
    widths = [12, 11, 11, 11, 15, 13]
    print_header(cols, widths)
    for row in rows:
        if row["pt_lo"] < 8 - 1e-3 or row["pt_hi"] > 36 + 1e-3:
            continue
        bin_str = f"{row['pt_lo']:.0f}-{row['pt_hi']:.0f}"
        print_row(
            [bin_str, row["CT14lo_10"], row["CT14nlo_10"], row["NNPDF_10"],
             row["ratio_CT14nlo_CT14lo"], row["ratio_NNPDF_CT14lo"]],
            widths,
            fmts=[None, ".4g", ".4g", ".4g", ".3f", ".3f"],
        )

    # Plot: 3 bands (CT14lo, CT14nlo, NNPDF) - ratio to CT14lo
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                                   gridspec_kw={"height_ratios": [3, 2]})
    colors = {"CT14lo": "C0", "CT14nlo": "C1", "NNPDF": "C2"}
    labels = {"CT14lo":  r"CT14lo",
              "CT14nlo": r"CT14nlo",
              "NNPDF":   r"NNPDF3.1 $\alpha_s$=0.118"}

    mid = 0.5 * (edges_ref[:-1] + edges_ref[1:])
    analysis_mask = (mid >= 8) & (mid <= 36)

    for pdf in ["CT14lo", "CT14nlo", "NNPDF"]:
        _, v_c, _ = data[pdf]["10"]
        _, v_lo, _ = data[pdf]["20"]  # mu=2  -> scale-down yields smaller xs
        _, v_hi, _ = data[pdf]["05"]  # mu=0.5 -> scale-up yields larger xs
        v_c  = v_c[analysis_mask]
        v_lo = v_lo[analysis_mask]
        v_hi = v_hi[analysis_mask]

        ax1.plot(mid[analysis_mask], v_c, marker='o', ls='-', color=colors[pdf], label=labels[pdf])
        ax1.fill_between(mid[analysis_mask], v_lo, v_hi, color=colors[pdf], alpha=0.25)

        # Ratio to CT14lo central
        _, v_ct14lo_c, _ = data["CT14lo"]["10"]
        v_ct14lo_c = v_ct14lo_c[analysis_mask]
        with np.errstate(divide="ignore", invalid="ignore"):
            r_c  = v_c  / v_ct14lo_c
            r_lo = v_lo / v_ct14lo_c
            r_hi = v_hi / v_ct14lo_c
        ax2.plot(mid[analysis_mask], r_c, marker='o', ls='-', color=colors[pdf], label=labels[pdf])
        ax2.fill_between(mid[analysis_mask], r_lo, r_hi, color=colors[pdf], alpha=0.25)

    ax1.set_yscale("log")
    ax1.set_ylabel(r"d$^2\sigma$/d$p_T$d$y$  [pb/GeV]")
    ax1.set_title(r"JETPHOX NLO: PPG12 fiducial ($\sqrt{s}=200$ GeV, $|y|<0.7$, iso R=0.3, $E_T^{iso}<4$ GeV)")
    ax1.legend(loc="upper right")
    ax1.grid(True, which="both", alpha=0.3)

    ax2.set_xlabel(r"$p_T$  [GeV]")
    ax2.set_ylabel("Ratio to CT14lo (central)")
    ax2.axhline(1.0, color="black", lw=0.8, alpha=0.5)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0.5, 2.0)

    plt.tight_layout()
    plot_path = f"{OUT_DIR}/compare_3pdf.pdf"
    plt.savefig(plot_path)
    plt.savefig(plot_path.replace(".pdf", ".png"), dpi=150)
    print(f"\nWrote {plot_path}")

    # Summary numbers
    print("\n--- Summary (mu=1.0 central values, 8-36 GeV quoted range) ---")
    for pdf in ["CT14lo", "CT14nlo", "NNPDF"]:
        _, v_c, _ = data[pdf]["10"]
        total = np.sum(v_c[analysis_mask] * np.diff(edges_ref)[analysis_mask])
        print(f"  {pdf:8s}: integrated xs 8-36 GeV = {total:.2f} pb")

    # Per-bin ratio summary
    print("\n--- Ratio to CT14lo (per bin), CT14nlo vs NNPDF ---")
    for row in rows:
        if row["pt_lo"] < 8 - 1e-3 or row["pt_hi"] > 36 + 1e-3:
            continue
        print(f"  pT [{row['pt_lo']:.0f},{row['pt_hi']:.0f}] GeV: "
              f"CT14nlo/CT14lo = {row['ratio_CT14nlo_CT14lo']:.3f}  "
              f"NNPDF/CT14lo = {row['ratio_NNPDF_CT14lo']:.3f}")

if __name__ == "__main__":
    main()
