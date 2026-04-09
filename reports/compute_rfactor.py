#!/usr/bin/env python3
"""
Compute ABCD R-factor from TRUE BACKGROUND isoET distributions.

R_true = (A_bkg * D_bkg) / (B_bkg * C_bkg)

where A=tight+iso, B=tight+noniso, C=nontight+iso, D=nontight+noniso.
The isolation cut is ET-dependent: isoET < reco_iso_max_b + reco_iso_max_s * ET.

Uses h_tight_recoisoET_background_0 / h_nontight_recoisoET_background_0 (TH2D)
from the inclusive MC file, which contain ONLY truth-matched non-signal clusters.
"""

import yaml
import uproot
import numpy as np

# --- Config ---
config_path = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/config_bdt_nom.yaml"
root_path = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root"

with open(config_path) as f:
    config = yaml.safe_load(f)

ana = config["analysis"]
iso_b = ana["reco_iso_max_b"]
iso_s = ana["reco_iso_max_s"]
noniso_shift = ana["reco_noniso_min_shift"]

pt_bins = ana["pT_bins"]  # [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
n_bins = len(pt_bins) - 1

print(f"Isolation cut: isoET < {iso_b:.6f} + {iso_s:.7f} * ET")
print(f"Non-isolated region: isoET > iso_cut + {noniso_shift}")
print(f"pT bins: {pt_bins}")
print(f"Number of bins: {n_bins}")
print()

# --- Open ROOT file and extract per-pT-bin background isoET ---
f = uproot.open(root_path)

# TH2D: x = reco pT (12 bins, 8-36 GeV), y = reco isoET (4400 bins, -5 to 50 GeV)
h2_tight = f["h_tight_recoisoET_background_0"]
h2_nontight = f["h_nontight_recoisoET_background_0"]

# 2D arrays: shape (12, 4400)
tight_2d = h2_tight.values()
nontight_2d = h2_nontight.values()
iso_edges = h2_tight.axis(1).edges()  # isoET bin edges

# --- Header ---
header = (
    f"{'pT bin':>12s} | {'ET ctr':>6s} | {'iso cut':>7s} | "
    f"{'A_bkg':>10s} | {'B_bkg':>10s} | {'C_bkg':>10s} | {'D_bkg':>10s} | "
    f"{'R_isoET':>8s} | {'f_iso_T':>8s} | {'f_iso_NT':>8s} | {'ratio':>8s}"
)
print(header)
print("-" * len(header))

results = []

for i in range(n_bins):
    pt_lo = pt_bins[i]
    pt_hi = pt_bins[i + 1]
    et_center = 0.5 * (pt_lo + pt_hi)

    # ET-dependent isolation cut
    iso_cut = iso_b + iso_s * et_center

    # Project pT bin i from TH2D → 1D isoET (background only)
    t_vals = tight_2d[i]
    nt_vals = nontight_2d[i]

    # Find the bin index where the isolation cut falls
    t_cut_bin = np.searchsorted(iso_edges, iso_cut) - 1
    nt_cut_bin = t_cut_bin  # same edges

    # A = tight + isolated (below cut), B = tight + non-isolated (above cut)
    A = np.sum(t_vals[:t_cut_bin])
    B = np.sum(t_vals[t_cut_bin:])

    # C = nontight + isolated (below cut), D = nontight + non-isolated (above cut)
    C = np.sum(nt_vals[:nt_cut_bin])
    D = np.sum(nt_vals[nt_cut_bin:])

    # R-factor
    if B * C > 0:
        R = (A * D) / (B * C)
    else:
        R = float("nan")

    # Isolated fractions
    f_iso_tight = A / (A + B) if (A + B) > 0 else float("nan")
    f_iso_nontight = C / (C + D) if (C + D) > 0 else float("nan")
    ratio = f_iso_tight / f_iso_nontight if f_iso_nontight > 0 else float("nan")

    results.append({
        "pt_lo": pt_lo, "pt_hi": pt_hi, "et_center": et_center,
        "iso_cut": iso_cut, "A": A, "B": B, "C": C, "D": D,
        "R": R, "f_iso_tight": f_iso_tight, "f_iso_nontight": f_iso_nontight,
        "ratio": ratio,
    })

    print(
        f"  [{pt_lo:2d},{pt_hi:2d}] GeV | {et_center:5.1f} | {iso_cut:7.3f} | "
        f"{A:10.1f} | {B:10.1f} | {C:10.1f} | {D:10.1f} | "
        f"{R:8.4f} | {f_iso_tight:8.4f} | {f_iso_nontight:8.4f} | {ratio:8.4f}"
    )

print()

# --- Summary ---
print("=" * 80)
print("SUMMARY")
print("=" * 80)

# Check where R crosses 1.0
prev_R = None
for r in results:
    if prev_R is not None and not np.isnan(r["R"]) and not np.isnan(prev_R):
        if (prev_R < 1.0 and r["R"] >= 1.0) or (prev_R > 1.0 and r["R"] <= 1.0):
            print(f"  R crosses 1.0 between pT [{r['pt_lo']-pt_bins[0]:.0f}] "
                  f"and [{r['pt_lo']},{r['pt_hi']}] GeV "
                  f"(R goes from {prev_R:.4f} to {r['R']:.4f})")
    prev_R = r["R"]

R_vals = [r["R"] for r in results if not np.isnan(r["R"])]
if R_vals:
    print(f"  R range: [{min(R_vals):.4f}, {max(R_vals):.4f}]")
    print(f"  R mean:  {np.mean(R_vals):.4f}")

# Check if R is consistently above or below 1
above = sum(1 for r in R_vals if r > 1.0)
below = sum(1 for r in R_vals if r <= 1.0)
print(f"  Bins with R > 1: {above}/{len(R_vals)}")
print(f"  Bins with R <= 1: {below}/{len(R_vals)}")
print()

# Also print isolated fractions
print("Isolated fraction comparison (tight vs non-tight BDT):")
for r in results:
    marker = " <-- R crosses 1" if abs(r["R"] - 1.0) < 0.05 else ""
    print(f"  pT [{r['pt_lo']:2d},{r['pt_hi']:2d}]: "
          f"f_iso_tight={r['f_iso_tight']:.4f}, "
          f"f_iso_nontight={r['f_iso_nontight']:.4f}, "
          f"ratio={r['ratio']:.4f}{marker}")
