#!/usr/bin/env python3
"""
R-factor vs isolation cut scan for TRUE BACKGROUND.

Reads background-only isoET distributions (tight and non-tight) from inclusive MC,
using TH2D h_tight_recoisoET_background_0 / h_nontight_recoisoET_background_0,
scans the isolation cut threshold, and computes R = (A*D)/(B*C) at each cut.

Key question: Is R ever < 1 for pure background, or is R > 1 everywhere?
"""

import numpy as np
import yaml
import uproot

# ── 1. Read config ──────────────────────────────────────────────────────────
config_path = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/config_bdt_nom.yaml"
with open(config_path) as f:
    config = yaml.safe_load(f)

reco_iso_max_b = config["analysis"]["reco_iso_max_b"]
reco_iso_max_s = config["analysis"]["reco_iso_max_s"]
pT_bins = config["analysis"]["pT_bins"]  # [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]

print(f"Nominal isolation: iso_cut = {reco_iso_max_b:.6f} + {reco_iso_max_s:.7f} * ET")
print(f"pT bins: {pT_bins}")
print(f"Number of pT bins: {len(pT_bins) - 1}")
print()

# ── 2. Read background-only TH2D and project per pT bin ───────────────────
root_file = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root"
f = uproot.open(root_file)

n_bins = len(pT_bins) - 1  # 12

# TH2D: x = reco pT (12 bins), y = reco isoET (4400 bins, -5 to 50 GeV)
h2_tight = f["h_tight_recoisoET_background_0"]
h2_nontight = f["h_nontight_recoisoET_background_0"]

tight_2d = h2_tight.values()      # shape (12, 4400)
nontight_2d = h2_nontight.values()
edges = h2_tight.axis(1).edges()   # isoET bin edges
bin_width = edges[1] - edges[0]

# ── 3. Scan isolation cut for 4 representative pT bins ─────────────────────
representative_bins = {
    0:  "[8, 10]",
    2:  "[12, 14]",
    5:  "[18, 20]",
    10: "[28, 32]",
}

iso_cuts = np.linspace(-2, 15, 100)

print("=" * 100)
print("ISOLATION CUT SCAN: R = (A*D)/(B*C) for pure jet MC background")
print("  A = tight & isolated,  B = tight & non-isolated")
print("  C = nontight & isolated,  D = nontight & non-isolated")
print("=" * 100)

results = {}  # bin_idx -> dict of arrays

for bin_idx, bin_label in representative_bins.items():
    pT_lo = pT_bins[bin_idx]
    pT_hi = pT_bins[bin_idx + 1]
    ET_center = (pT_lo + pT_hi) / 2.0
    nominal_cut = reco_iso_max_b + reco_iso_max_s * ET_center

    tight_vals = tight_2d[bin_idx]
    nontight_vals = nontight_2d[bin_idx]

    R_arr = np.full(len(iso_cuts), np.nan)
    f_iso_tight_arr = np.full(len(iso_cuts), np.nan)
    f_iso_nontight_arr = np.full(len(iso_cuts), np.nan)

    for j, cut in enumerate(iso_cuts):
        # Find the bin index corresponding to this cut value
        # edges[k] <= cut < edges[k+1] means bins 0..k are below, k+1..end are above
        cut_bin = np.searchsorted(edges, cut, side='right') - 1
        # Clamp
        cut_bin = max(0, min(cut_bin, len(tight_vals) - 1))

        # Integrate: A = sum of tight bins with isoET < cut (bins 0..cut_bin-1 inclusive)
        # Note: bin i covers [edges[i], edges[i+1])
        # "below cut" = bins where upper edge <= cut, i.e., edges[i+1] <= cut
        # More precisely: use bin centers. But simpler: integrate up to cut_bin.
        A = np.sum(tight_vals[:cut_bin + 1])      # tight, isoET <= cut (isolated)
        B = np.sum(tight_vals[cut_bin + 1:])       # tight, isoET > cut (non-isolated)
        C = np.sum(nontight_vals[:cut_bin + 1])    # nontight, isoET <= cut
        D = np.sum(nontight_vals[cut_bin + 1:])    # nontight, isoET > cut

        if B > 0 and C > 0:
            R_arr[j] = (A * D) / (B * C)
        if (A + B) > 0:
            f_iso_tight_arr[j] = A / (A + B)
        if (C + D) > 0:
            f_iso_nontight_arr[j] = C / (C + D)

    results[bin_idx] = {
        "iso_cuts": iso_cuts,
        "R": R_arr,
        "f_iso_tight": f_iso_tight_arr,
        "f_iso_nontight": f_iso_nontight_arr,
        "nominal_cut": nominal_cut,
        "ET_center": ET_center,
    }

    # Print table
    print(f"\n{'='*90}")
    print(f"pT bin {bin_idx}: {bin_label} GeV   |   ET_center = {ET_center:.1f} GeV   |   nominal iso cut = {nominal_cut:.3f} GeV")
    print(f"{'='*90}")
    print(f"{'iso_cut':>10s} | {'R':>12s} | {'f_iso_tight':>12s} | {'f_iso_nontight':>14s} | {'note':>10s}")
    print(f"{'-'*10}-+-{'-'*12}-+-{'-'*12}-+-{'-'*14}-+-{'-'*10}")

    # Print every 5th entry plus the one nearest the nominal cut
    nominal_idx = np.argmin(np.abs(iso_cuts - nominal_cut))
    print_indices = set(range(0, len(iso_cuts), 5))
    print_indices.add(nominal_idx)

    for j in sorted(print_indices):
        note = ""
        if j == nominal_idx:
            note = "<-- NOM"
        R_str = f"{R_arr[j]:.6f}" if not np.isnan(R_arr[j]) else "nan"
        ft_str = f"{f_iso_tight_arr[j]:.6f}" if not np.isnan(f_iso_tight_arr[j]) else "nan"
        fn_str = f"{f_iso_nontight_arr[j]:.6f}" if not np.isnan(f_iso_nontight_arr[j]) else "nan"
        print(f"{iso_cuts[j]:10.3f} | {R_str:>12s} | {ft_str:>12s} | {fn_str:>14s} | {note:>10s}")

    # Summary stats for this bin
    valid_R = R_arr[~np.isnan(R_arr)]
    if len(valid_R) > 0:
        print(f"\n  R range: [{valid_R.min():.6f}, {valid_R.max():.6f}]")
        print(f"  R at nominal cut: {R_arr[nominal_idx]:.6f}")
        below_one = np.sum(valid_R < 1.0)
        print(f"  Number of scan points with R < 1: {below_one} / {len(valid_R)}")
        if below_one > 0:
            r_below = iso_cuts[~np.isnan(R_arr) & (R_arr < 1.0)]
            print(f"  iso_cut range where R < 1: [{r_below.min():.3f}, {r_below.max():.3f}]")
    else:
        print(f"\n  No valid R values (likely zero counts everywhere)")


# ── 4. Nominal R for all 12 pT bins ────────────────────────────────────────
print("\n\n")
print("=" * 100)
print("R AT NOMINAL ISOLATION CUT FOR ALL 12 pT BINS (pure jet MC)")
print("=" * 100)
print(f"{'Bin':>4s} | {'pT range':>12s} | {'ET_center':>10s} | {'iso_cut':>10s} | {'A':>14s} | {'B':>14s} | {'C':>14s} | {'D':>14s} | {'R':>12s} | {'f_iso_t':>10s} | {'f_iso_nt':>10s}")
print("-" * 140)

all_R_nominal = []
for i in range(n_bins):
    pT_lo = pT_bins[i]
    pT_hi = pT_bins[i + 1]
    ET_center = (pT_lo + pT_hi) / 2.0
    nominal_cut = reco_iso_max_b + reco_iso_max_s * ET_center

    tight_vals = tight_2d[i]
    nontight_vals = nontight_2d[i]

    cut_bin = np.searchsorted(edges, nominal_cut, side='right') - 1
    cut_bin = max(0, min(cut_bin, len(tight_vals) - 1))

    A = np.sum(tight_vals[:cut_bin + 1])
    B = np.sum(tight_vals[cut_bin + 1:])
    C = np.sum(nontight_vals[:cut_bin + 1])
    D = np.sum(nontight_vals[cut_bin + 1:])

    if B > 0 and C > 0:
        R = (A * D) / (B * C)
    else:
        R = np.nan

    f_iso_t = A / (A + B) if (A + B) > 0 else np.nan
    f_iso_nt = C / (C + D) if (C + D) > 0 else np.nan

    all_R_nominal.append(R)

    R_str = f"{R:.6f}" if not np.isnan(R) else "nan"
    ft_str = f"{f_iso_t:.6f}" if not np.isnan(f_iso_t) else "nan"
    fn_str = f"{f_iso_nt:.6f}" if not np.isnan(f_iso_nt) else "nan"

    print(f"{i:4d} | [{pT_lo:5.0f},{pT_hi:4.0f}] | {ET_center:10.1f} | {nominal_cut:10.4f} | {A:14.1f} | {B:14.1f} | {C:14.1f} | {D:14.1f} | {R_str:>12s} | {ft_str:>10s} | {fn_str:>10s}")

# ── 5. Key question summary ────────────────────────────────────────────────
print("\n\n")
print("=" * 100)
print("KEY FINDINGS")
print("=" * 100)

all_R_nominal = np.array(all_R_nominal)
valid_nominal = all_R_nominal[~np.isnan(all_R_nominal)]

print(f"\nNominal R across all 12 pT bins:")
print(f"  min R = {valid_nominal.min():.6f}")
print(f"  max R = {valid_nominal.max():.6f}")
print(f"  mean R = {valid_nominal.mean():.6f}")
print(f"  All R > 1? {np.all(valid_nominal > 1.0)}")
print(f"  All R < 1? {np.all(valid_nominal < 1.0)}")

print(f"\nIsolation cut scan (4 representative bins, iso_cut from -2 to 15 GeV):")
for bin_idx, bin_label in representative_bins.items():
    res = results[bin_idx]
    valid_R = res["R"][~np.isnan(res["R"])]
    below_one = np.sum(valid_R < 1.0)
    pct = 100 * below_one / len(valid_R) if len(valid_R) > 0 else 0
    rmin = valid_R.min() if len(valid_R) > 0 else np.nan
    rmax = valid_R.max() if len(valid_R) > 0 else np.nan
    print(f"  Bin {bin_idx} {bin_label}: R in [{rmin:.4f}, {rmax:.4f}], "
          f"R < 1 at {below_one}/{len(valid_R)} points ({pct:.1f}%)")

print(f"\nAnswer: Is there ANY isolation cut where R < 1 for pure background?")
any_below = False
for bin_idx in representative_bins:
    res = results[bin_idx]
    valid_R = res["R"][~np.isnan(res["R"])]
    if np.any(valid_R < 1.0):
        any_below = True
        break

if any_below:
    print("  YES -- R < 1 exists at some isolation cut values.")
    print("  This means the ABCD assumption (background factorizes) is violated in a way")
    print("  that makes R cross below 1 for some cut choices.")
else:
    print("  NO -- R >= 1 for ALL scanned cut values in ALL representative bins.")
    print("  The correlation between isolation and shower shape is such that R > 1 always.")

print("\nPhysics interpretation:")
print("  R = (A*D)/(B*C) = 1 means perfect factorization (iso and BDT are independent for bkg)")
print("  R > 1 means tight clusters are MORE likely to be isolated than non-tight ones")
print("  R < 1 means tight clusters are LESS likely to be isolated than non-tight ones")
print("  For pure background, R != 1 indicates correlation between shower shape and isolation")
print("  The ABCD method assumes R = 1; deviation from 1 is a source of systematic bias")
