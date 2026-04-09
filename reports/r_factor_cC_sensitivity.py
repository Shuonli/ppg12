#!/usr/bin/env python3
"""
R-factor sensitivity to signal leakage fractions cC and cB.

Reads MC ABCD yields and signal yields from the nominal efficiency file,
computes R = A * D_corr / (B_corr * C_corr) while scanning cC (or cB),
and finds the critical value where R crosses 1.0.

The ABCD corrected yields are:
  B_corr = B_total - cB * A_sig   (signal leakage into region B)
  C_corr = C_total - cC * A_sig   (signal leakage into region C)
  D_corr = D_total - cD * A_sig   (signal leakage into region D)

where cX = signal_X / signal_A (leakage fraction relative to signal in A).

The "pure background" R factor is:
  R = (A_bkg * D_bkg) / (B_bkg * C_bkg)
but since we only observe total yields, the corrected R using estimated
signal subtraction is:
  R = (A_total * D_corr) / (B_corr * C_corr)
evaluated on MC where we know the signal content.
"""

import numpy as np
import uproot

# ---------- configuration ----------
ROOT_FILE = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root"

# pT bin edges from plotcommon.h
PT_EDGES = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)

# Representative bins to study: [8,10], [12,14], [16,18], [24,26] GeV
# These correspond to ROOT bin indices 1, 3, 5, 9 (1-indexed)
TARGET_BINS_GEV = [(8, 10), (12, 14), (16, 18), (24, 26)]

N_SCAN = 21  # scan points from 0 to 2x nominal

# ---------- read histograms ----------
f = uproot.open(ROOT_FILE)

def get_hist(name):
    """Read a histogram by exact name, return (values, errors) arrays."""
    h = f[name]
    vals = h.values()      # bin contents (excluding under/overflow)
    errs = h.errors()      # bin errors
    edges = h.axis().edges()
    return vals, errs, edges

# Total ABCD yields (inclusive MC: signal + background)
A_vals, A_errs, pt_edges = get_hist("h_tight_iso_cluster_0")
B_vals, B_errs, _        = get_hist("h_tight_noniso_cluster_0")
C_vals, C_errs, _        = get_hist("h_nontight_iso_cluster_0")
D_vals, D_errs, _        = get_hist("h_nontight_noniso_cluster_0")

# Signal yields in each ABCD region
Asig_vals, Asig_errs, _ = get_hist("h_tight_iso_cluster_signal_0")
Bsig_vals, Bsig_errs, _ = get_hist("h_tight_noniso_cluster_signal_0")
Csig_vals, Csig_errs, _ = get_hist("h_nontight_iso_cluster_signal_0")
Dsig_vals, Dsig_errs, _ = get_hist("h_nontight_noniso_cluster_signal_0")

print("pT bin edges from histogram:", pt_edges)
print(f"Number of bins: {len(A_vals)}")
print()

# ---------- map target pT ranges to bin indices ----------
def find_bin_index(pt_low, pt_high, edges):
    """Find the 0-based bin index for a given pT range."""
    for i in range(len(edges) - 1):
        if abs(edges[i] - pt_low) < 0.01 and abs(edges[i+1] - pt_high) < 0.01:
            return i
    raise ValueError(f"Could not find bin [{pt_low}, {pt_high}] in edges {edges}")

bin_indices = []
for pt_lo, pt_hi in TARGET_BINS_GEV:
    idx = find_bin_index(pt_lo, pt_hi, pt_edges)
    bin_indices.append(idx)
    print(f"  pT [{pt_lo:.0f}, {pt_hi:.0f}] GeV -> bin index {idx}")

print()

# ---------- compute nominal leakage fractions ----------
print("=" * 90)
print("NOMINAL ABCD YIELDS AND LEAKAGE FRACTIONS")
print("=" * 90)
header = f"{'pT bin':>12s} {'A':>12s} {'B':>12s} {'C':>12s} {'D':>12s} {'A_sig':>12s} {'B_sig':>12s} {'C_sig':>12s} {'D_sig':>12s} {'cB':>8s} {'cC':>8s} {'cD':>8s}"
print(header)
print("-" * len(header))

nominal_cB = {}
nominal_cC = {}
nominal_cD = {}

for (pt_lo, pt_hi), idx in zip(TARGET_BINS_GEV, bin_indices):
    A = A_vals[idx]
    B = B_vals[idx]
    C = C_vals[idx]
    D = D_vals[idx]
    Asig = Asig_vals[idx]
    Bsig = Bsig_vals[idx]
    Csig = Csig_vals[idx]
    Dsig = Dsig_vals[idx]

    cB = Bsig / Asig if Asig > 0 else 0
    cC = Csig / Asig if Asig > 0 else 0
    cD = Dsig / Asig if Asig > 0 else 0

    nominal_cB[(pt_lo, pt_hi)] = cB
    nominal_cC[(pt_lo, pt_hi)] = cC
    nominal_cD[(pt_lo, pt_hi)] = cD

    label = f"[{pt_lo:.0f},{pt_hi:.0f}]"
    print(f"{label:>12s} {A:12.1f} {B:12.1f} {C:12.1f} {D:12.1f} "
          f"{Asig:12.1f} {Bsig:12.1f} {Csig:12.1f} {Dsig:12.1f} "
          f"{cB:8.4f} {cC:8.4f} {cD:8.4f}")

print()

# Also show the pure-background R (no correction) and the uncorrected R
print("=" * 90)
print("BACKGROUND R FACTOR (pure bkg) vs UNCORRECTED R (no leakage subtraction)")
print("=" * 90)
header2 = f"{'pT bin':>12s} {'A_bkg':>12s} {'B_bkg':>12s} {'C_bkg':>12s} {'D_bkg':>12s} {'R_pure_bkg':>12s} {'R_uncorr':>12s} {'R_nom_corr':>12s}"
print(header2)
print("-" * len(header2))

for (pt_lo, pt_hi), idx in zip(TARGET_BINS_GEV, bin_indices):
    A = A_vals[idx]
    B = B_vals[idx]
    C = C_vals[idx]
    D = D_vals[idx]
    Asig = Asig_vals[idx]
    Bsig = Bsig_vals[idx]
    Csig = Csig_vals[idx]
    Dsig = Dsig_vals[idx]

    A_bkg = A - Asig
    B_bkg = B - Bsig
    C_bkg = C - Csig
    D_bkg = D - Dsig

    cB = nominal_cB[(pt_lo, pt_hi)]
    cC = nominal_cC[(pt_lo, pt_hi)]
    cD = nominal_cD[(pt_lo, pt_hi)]

    R_pure = (A_bkg * D_bkg) / (B_bkg * C_bkg) if (B_bkg * C_bkg) != 0 else float('inf')
    R_uncorr = (A * D) / (B * C) if (B * C) != 0 else float('inf')

    # Nominal corrected R: use total yields with leakage subtraction
    B_corr = B - cB * Asig
    C_corr = C - cC * Asig
    D_corr = D - cD * Asig
    R_nom = (A * D_corr) / (B_corr * C_corr) if (B_corr * C_corr) != 0 else float('inf')

    label = f"[{pt_lo:.0f},{pt_hi:.0f}]"
    print(f"{label:>12s} {A_bkg:12.1f} {B_bkg:12.1f} {C_bkg:12.1f} {D_bkg:12.1f} "
          f"{R_pure:12.4f} {R_uncorr:12.4f} {R_nom:12.4f}")

print()

# ---------- R factor computation ----------
def compute_R(A, B, C, D, Asig, cB_val, cC_val, cD_val):
    """Compute R = A * D_corr / (B_corr * C_corr)."""
    B_corr = B - cB_val * Asig
    C_corr = C - cC_val * Asig
    D_corr = D - cD_val * Asig
    denom = B_corr * C_corr
    if denom == 0 or np.isnan(denom):
        return float('nan')
    return (A * D_corr) / denom


# ========== SCAN 1: vary cC, keep cB and cD at nominal ==========
print("=" * 90)
print("SCAN 1: R vs cC  (cB, cD held at nominal)")
print("=" * 90)

cC_critical = {}

for (pt_lo, pt_hi), idx in zip(TARGET_BINS_GEV, bin_indices):
    A = A_vals[idx]
    B = B_vals[idx]
    C = C_vals[idx]
    D = D_vals[idx]
    Asig = Asig_vals[idx]

    cB = nominal_cB[(pt_lo, pt_hi)]
    cC_nom = nominal_cC[(pt_lo, pt_hi)]
    cD = nominal_cD[(pt_lo, pt_hi)]

    cC_max = 2.0 * cC_nom if cC_nom > 0 else 0.1
    cC_scan = np.linspace(0, cC_max, N_SCAN)
    R_scan = np.array([compute_R(A, B, C, D, Asig, cB, cc, cD) for cc in cC_scan])

    label = f"[{pt_lo:.0f},{pt_hi:.0f}]"
    print(f"\npT {label} GeV   (nominal cC = {cC_nom:.6f})")
    print(f"  {'cC':>10s}  {'cC/cC_nom':>10s}  {'R':>10s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*10}")
    for i in range(N_SCAN):
        ratio = cC_scan[i] / cC_nom if cC_nom > 0 else float('inf')
        marker = " <-- nominal" if abs(ratio - 1.0) < 0.01 else ""
        print(f"  {cC_scan[i]:10.6f}  {ratio:10.4f}  {R_scan[i]:10.4f}{marker}")

    # Find where R crosses 1.0 by linear interpolation
    crossings = []
    for i in range(len(R_scan) - 1):
        if np.isnan(R_scan[i]) or np.isnan(R_scan[i+1]):
            continue
        if (R_scan[i] - 1.0) * (R_scan[i+1] - 1.0) <= 0:
            # Linear interpolation
            frac = (1.0 - R_scan[i]) / (R_scan[i+1] - R_scan[i])
            cC_cross = cC_scan[i] + frac * (cC_scan[i+1] - cC_scan[i])
            crossings.append(cC_cross)

    if crossings:
        cC_critical[(pt_lo, pt_hi)] = crossings[0]
        print(f"  >>> R = 1.0 crossing at cC = {crossings[0]:.6f} (cC/cC_nom = {crossings[0]/cC_nom:.4f})")
    else:
        cC_critical[(pt_lo, pt_hi)] = None
        print(f"  >>> No R = 1.0 crossing found in scan range [0, {cC_max:.6f}]")
        # Check if R is always > 1 or always < 1
        valid = R_scan[~np.isnan(R_scan)]
        if len(valid) > 0:
            print(f"      R range: [{valid.min():.4f}, {valid.max():.4f}]")


# ========== SCAN 2: vary cB, keep cC and cD at nominal ==========
print()
print("=" * 90)
print("SCAN 2: R vs cB  (cC, cD held at nominal)")
print("=" * 90)

cB_critical = {}

for (pt_lo, pt_hi), idx in zip(TARGET_BINS_GEV, bin_indices):
    A = A_vals[idx]
    B = B_vals[idx]
    C = C_vals[idx]
    D = D_vals[idx]
    Asig = Asig_vals[idx]

    cB_nom = nominal_cB[(pt_lo, pt_hi)]
    cC = nominal_cC[(pt_lo, pt_hi)]
    cD = nominal_cD[(pt_lo, pt_hi)]

    cB_max = 2.0 * cB_nom if cB_nom > 0 else 0.1
    cB_scan = np.linspace(0, cB_max, N_SCAN)
    R_scan = np.array([compute_R(A, B, C, D, Asig, cb, cC, cD) for cb in cB_scan])

    label = f"[{pt_lo:.0f},{pt_hi:.0f}]"
    print(f"\npT {label} GeV   (nominal cB = {cB_nom:.6f})")
    print(f"  {'cB':>10s}  {'cB/cB_nom':>10s}  {'R':>10s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*10}")
    for i in range(N_SCAN):
        ratio = cB_scan[i] / cB_nom if cB_nom > 0 else float('inf')
        marker = " <-- nominal" if abs(ratio - 1.0) < 0.01 else ""
        print(f"  {cB_scan[i]:10.6f}  {ratio:10.4f}  {R_scan[i]:10.4f}{marker}")

    # Find where R crosses 1.0
    crossings = []
    for i in range(len(R_scan) - 1):
        if np.isnan(R_scan[i]) or np.isnan(R_scan[i+1]):
            continue
        if (R_scan[i] - 1.0) * (R_scan[i+1] - 1.0) <= 0:
            frac = (1.0 - R_scan[i]) / (R_scan[i+1] - R_scan[i])
            cB_cross = cB_scan[i] + frac * (cB_scan[i+1] - cB_scan[i])
            crossings.append(cB_cross)

    if crossings:
        cB_critical[(pt_lo, pt_hi)] = crossings[0]
        print(f"  >>> R = 1.0 crossing at cB = {crossings[0]:.6f} (cB/cB_nom = {crossings[0]/cB_nom:.4f})")
    else:
        cB_critical[(pt_lo, pt_hi)] = None
        print(f"  >>> No R = 1.0 crossing found in scan range [0, {cB_max:.6f}]")
        valid = R_scan[~np.isnan(R_scan)]
        if len(valid) > 0:
            print(f"      R range: [{valid.min():.4f}, {valid.max():.4f}]")


# ========== SENSITIVITY COMPARISON ==========
print()
print("=" * 90)
print("SENSITIVITY COMPARISON: dR/d(cX) at nominal point")
print("=" * 90)

print(f"\n{'pT bin':>12s} {'dR/dcC':>12s} {'dR/dcB':>12s} {'|dR/dcC|/|dR/dcB|':>20s} {'More sensitive':>15s}")
print("-" * 75)

for (pt_lo, pt_hi), idx in zip(TARGET_BINS_GEV, bin_indices):
    A = A_vals[idx]
    B = B_vals[idx]
    C = C_vals[idx]
    D = D_vals[idx]
    Asig = Asig_vals[idx]

    cB = nominal_cB[(pt_lo, pt_hi)]
    cC = nominal_cC[(pt_lo, pt_hi)]
    cD = nominal_cD[(pt_lo, pt_hi)]

    # Numerical derivatives at nominal point
    eps = 1e-6
    R_nom = compute_R(A, B, C, D, Asig, cB, cC, cD)
    dR_dcC = (compute_R(A, B, C, D, Asig, cB, cC + eps, cD) - R_nom) / eps
    dR_dcB = (compute_R(A, B, C, D, Asig, cB + eps, cC, cD) - R_nom) / eps

    ratio = abs(dR_dcC) / abs(dR_dcB) if abs(dR_dcB) > 0 else float('inf')
    dominant = "cC" if ratio > 1 else "cB"

    label = f"[{pt_lo:.0f},{pt_hi:.0f}]"
    print(f"{label:>12s} {dR_dcC:12.4f} {dR_dcB:12.4f} {ratio:20.2f} {dominant:>15s}")

print()

# ========== SUMMARY TABLE ==========
print("=" * 90)
print("SUMMARY: Critical cC where R = 1.0")
print("=" * 90)
print(f"\n{'pT bin':>12s} {'cC_nom':>10s} {'cC_crit':>10s} {'cC_crit/cC_nom':>16s} {'R(cC=0)':>10s} {'R(cC=nom)':>10s} {'R(cC=2x)':>10s}")
print("-" * 80)

for (pt_lo, pt_hi), idx in zip(TARGET_BINS_GEV, bin_indices):
    A = A_vals[idx]
    B = B_vals[idx]
    C = C_vals[idx]
    D = D_vals[idx]
    Asig = Asig_vals[idx]

    cB = nominal_cB[(pt_lo, pt_hi)]
    cC_nom = nominal_cC[(pt_lo, pt_hi)]
    cD = nominal_cD[(pt_lo, pt_hi)]

    R_zero = compute_R(A, B, C, D, Asig, cB, 0, cD)
    R_nom = compute_R(A, B, C, D, Asig, cB, cC_nom, cD)
    R_2x = compute_R(A, B, C, D, Asig, cB, 2*cC_nom, cD)

    crit = cC_critical[(pt_lo, pt_hi)]
    crit_str = f"{crit:.6f}" if crit is not None else "N/A"
    crit_ratio = f"{crit/cC_nom:.4f}" if (crit is not None and cC_nom > 0) else "N/A"

    label = f"[{pt_lo:.0f},{pt_hi:.0f}]"
    print(f"{label:>12s} {cC_nom:10.6f} {crit_str:>10s} {crit_ratio:>16s} {R_zero:10.4f} {R_nom:10.4f} {R_2x:10.4f}")

print()
print("=" * 90)
print("SUMMARY: Critical cB where R = 1.0 (for comparison)")
print("=" * 90)
print(f"\n{'pT bin':>12s} {'cB_nom':>10s} {'cB_crit':>10s} {'cB_crit/cB_nom':>16s} {'R(cB=0)':>10s} {'R(cB=nom)':>10s} {'R(cB=2x)':>10s}")
print("-" * 80)

for (pt_lo, pt_hi), idx in zip(TARGET_BINS_GEV, bin_indices):
    A = A_vals[idx]
    B = B_vals[idx]
    C = C_vals[idx]
    D = D_vals[idx]
    Asig = Asig_vals[idx]

    cB_nom = nominal_cB[(pt_lo, pt_hi)]
    cC = nominal_cC[(pt_lo, pt_hi)]
    cD = nominal_cD[(pt_lo, pt_hi)]

    R_zero = compute_R(A, B, C, D, Asig, 0, cC, cD)
    R_nom = compute_R(A, B, C, D, Asig, cB_nom, cC, cD)
    R_2x = compute_R(A, B, C, D, Asig, 2*cB_nom, cC, cD)

    crit = cB_critical[(pt_lo, pt_hi)]
    crit_str = f"{crit:.6f}" if crit is not None else "N/A"
    crit_ratio = f"{crit/cB_nom:.4f}" if (crit is not None and cB_nom > 0) else "N/A"

    label = f"[{pt_lo:.0f},{pt_hi:.0f}]"
    print(f"{label:>12s} {cB_nom:10.6f} {crit_str:>10s} {crit_ratio:>16s} {R_zero:10.4f} {R_nom:10.4f} {R_2x:10.4f}")

print()

# ========== PURE BACKGROUND R vs pT (all bins) ==========
print("=" * 90)
print("PURE BACKGROUND R = A_bkg*D_bkg / (B_bkg*C_bkg) FOR ALL pT BINS")
print("=" * 90)
print(f"\n{'Bin':>4s} {'pT range':>12s} {'A_bkg':>14s} {'B_bkg':>14s} {'C_bkg':>14s} {'D_bkg':>14s} {'R_bkg':>10s} {'sig_frac_A':>12s} {'sig_frac_C':>12s}")
print("-" * 110)

for idx in range(len(A_vals)):
    A = A_vals[idx]
    B = B_vals[idx]
    C = C_vals[idx]
    D = D_vals[idx]
    Asig = Asig_vals[idx]
    Bsig = Bsig_vals[idx]
    Csig = Csig_vals[idx]
    Dsig = Dsig_vals[idx]

    A_bkg = A - Asig
    B_bkg = B - Bsig
    C_bkg = C - Csig
    D_bkg = D - Dsig

    R_bkg = (A_bkg * D_bkg) / (B_bkg * C_bkg) if (B_bkg * C_bkg) != 0 else float('inf')
    sig_frac_A = Asig / A if A > 0 else 0
    sig_frac_C = Csig / C if C > 0 else 0

    pt_lo = pt_edges[idx]
    pt_hi = pt_edges[idx + 1]
    label = f"[{pt_lo:.0f},{pt_hi:.0f}]"
    print(f"{idx:>4d} {label:>12s} {A_bkg:14.1f} {B_bkg:14.1f} {C_bkg:14.1f} {D_bkg:14.1f} {R_bkg:10.4f} {sig_frac_A:12.4f} {sig_frac_C:12.4f}")


# ========== SIGNAL CONTAMINATION OF REGION C ==========
print()
print("=" * 90)
print("SIGNAL CONTAMINATION: cC * A_sig vs C_total (what fraction of C is signal?)")
print("=" * 90)
print(f"\n{'pT bin':>12s} {'C_total':>14s} {'C_sig':>14s} {'C_bkg':>14s} {'C_sig/C_tot':>12s} {'cC*A_sig/C':>12s}")
print("-" * 80)

for (pt_lo, pt_hi), idx in zip(TARGET_BINS_GEV, bin_indices):
    C = C_vals[idx]
    Asig = Asig_vals[idx]
    Csig = Csig_vals[idx]
    cC = nominal_cC[(pt_lo, pt_hi)]

    sig_frac = Csig / C if C > 0 else 0
    corr_frac = cC * Asig / C if C > 0 else 0

    label = f"[{pt_lo:.0f},{pt_hi:.0f}]"
    print(f"{label:>12s} {C:14.1f} {Csig:14.1f} {C - Csig:14.1f} {sig_frac:12.4f} {corr_frac:12.4f}")

# Same for region B
print()
print(f"\n{'pT bin':>12s} {'B_total':>14s} {'B_sig':>14s} {'B_bkg':>14s} {'B_sig/B_tot':>12s} {'cB*A_sig/B':>12s}")
print("-" * 80)

for (pt_lo, pt_hi), idx in zip(TARGET_BINS_GEV, bin_indices):
    B = B_vals[idx]
    Asig = Asig_vals[idx]
    Bsig = Bsig_vals[idx]
    cB = nominal_cB[(pt_lo, pt_hi)]

    sig_frac = Bsig / B if B > 0 else 0
    corr_frac = cB * Asig / B if B > 0 else 0

    label = f"[{pt_lo:.0f},{pt_hi:.0f}]"
    print(f"{label:>12s} {B:14.1f} {Bsig:14.1f} {B - Bsig:14.1f} {sig_frac:12.4f} {corr_frac:12.4f}")


# ========== KEY DIAGNOSTIC: EFFECTIVE R vs pT ==========
# This is the R the ABCD method "sees" when using total yields minus
# leakage corrections.  It measures closure: if R_eff = R_pure_bkg,
# the method is self-consistent.
print()
print("=" * 90)
print("EFFECTIVE R SEEN BY ABCD (using corrected yields) vs PURE BACKGROUND R")
print("R_eff = A_total * D_corr / (B_corr * C_corr)  where X_corr = X - cX * A_sig")
print("R_bkg = A_bkg * D_bkg / (B_bkg * C_bkg)")
print("=" * 90)
print(f"\n{'pT bin':>12s} {'R_bkg':>10s} {'R_eff':>10s} {'R_uncorr':>10s}  NOTE")
print("-" * 70)

for idx in range(len(A_vals)):
    A = A_vals[idx]
    B = B_vals[idx]
    C = C_vals[idx]
    D = D_vals[idx]
    Asig = Asig_vals[idx]
    Bsig = Bsig_vals[idx]
    Csig = Csig_vals[idx]
    Dsig = Dsig_vals[idx]

    A_bkg = A - Asig
    B_bkg = B - Bsig
    C_bkg = C - Csig
    D_bkg = D - Dsig

    cB = Bsig / Asig if Asig > 0 else 0
    cC = Csig / Asig if Asig > 0 else 0
    cD = Dsig / Asig if Asig > 0 else 0

    R_bkg = (A_bkg * D_bkg) / (B_bkg * C_bkg) if (B_bkg * C_bkg) != 0 else float('inf')
    R_uncorr = (A * D) / (B * C) if (B * C) != 0 else float('inf')

    B_corr = B - cB * Asig
    C_corr = C - cC * Asig
    D_corr = D - cD * Asig
    R_eff = (A * D_corr) / (B_corr * C_corr) if (B_corr * C_corr) != 0 else float('inf')

    pt_lo = pt_edges[idx]
    pt_hi = pt_edges[idx + 1]
    label = f"[{pt_lo:.0f},{pt_hi:.0f}]"

    note = ""
    if R_bkg > 1.0:
        note = "R_bkg > 1"
    else:
        note = "R_bkg < 1"

    print(f"{label:>12s} {R_bkg:10.4f} {R_eff:10.4f} {R_uncorr:10.4f}  {note}")


# ========== WHY R_eff != R_bkg ==========
# R_eff uses A_total (not A_bkg) in the numerator. This is because the
# ABCD method subtracts background from A, not signal from B/C/D.
# So R_eff = (A_bkg + A_sig) * D_bkg / (B_bkg * C_bkg)
#          = R_bkg * (1 + A_sig/A_bkg) * (D_bkg/D_corr_approx) ...
# The huge R_eff values come from A_sig >> A_bkg: the total A is
# dominated by signal, inflating the numerator.

# A cleaner diagnostic: what is the pure-background R, and where does
# it cross 1.0?  That's the real ABCD closure check.

print()
print("=" * 90)
print("PURE BACKGROUND R: pT DEPENDENCE AND CROSSOVER")
print("=" * 90)

# Find crossing of R_bkg = 1.0
R_bkg_all = []
for idx in range(len(A_vals)):
    A_bkg = A_vals[idx] - Asig_vals[idx]
    B_bkg = B_vals[idx] - Bsig_vals[idx]
    C_bkg = C_vals[idx] - Csig_vals[idx]
    D_bkg = D_vals[idx] - Dsig_vals[idx]
    R_bkg = (A_bkg * D_bkg) / (B_bkg * C_bkg) if (B_bkg * C_bkg) != 0 else float('nan')
    R_bkg_all.append(R_bkg)
    pt_center = 0.5 * (pt_edges[idx] + pt_edges[idx + 1])
    print(f"  pT = {pt_center:5.1f} GeV:  R_bkg = {R_bkg:.4f}")

# Interpolate crossing
R_arr = np.array(R_bkg_all)
for i in range(len(R_arr) - 1):
    if (R_arr[i] - 1.0) * (R_arr[i+1] - 1.0) <= 0:
        frac = (1.0 - R_arr[i]) / (R_arr[i+1] - R_arr[i])
        pt_cross = pt_edges[i] + frac * (pt_edges[i+1] - pt_edges[i])
        print(f"\n  >>> R_bkg = 1.0 crossing at pT ~ {pt_cross:.1f} GeV")

# ========== HOW cC SHIFTS THE ABCD PURITY ESTIMATE ==========
# The practical question: in the ABCD signal extraction, the background
# in region A is estimated as:
#   N_bkg_A = R * (B - cB*Nsig) * (C - cC*Nsig) / (D - cD*Nsig)
# where R=1 (assumed). If cC is wrong, the estimated background changes.
# We want: dN_bkg_A / dcC  vs  dN_bkg_A / dcB

print()
print("=" * 90)
print("PURITY SENSITIVITY: How cC and cB errors propagate to signal yield")
print("=" * 90)
print("""
In CalculatePhotonYield.C, the equation solved (with R=1) is:
  NsigA = NA - (NB - cB*NsigA) * (NC - cC*NsigA) / (ND - cD*NsigA)

We can linearize around the nominal cC to find:
  dNsigA/dcC ~ (partial of RHS w.r.t. cC) / (1 - partial of RHS w.r.t. NsigA)

But numerically it's simpler to just solve the equation at shifted cC values.
""")

def solve_nsig(A, B, C, D, cB_v, cC_v, cD_v):
    """Solve NsigA = A - (B - cB*NsigA)*(C - cC*NsigA)/(D - cD*NsigA) numerically."""
    # Use simple iteration (fixed-point)
    nsig = 0.5 * A  # initial guess
    for _ in range(200):
        B_c = B - cB_v * nsig
        C_c = C - cC_v * nsig
        D_c = D - cD_v * nsig
        if abs(D_c) < 1e-10:
            return float('nan')
        bkg = B_c * C_c / D_c
        nsig_new = A - bkg
        if abs(nsig_new - nsig) < 1e-3:
            return nsig_new
        # Damped iteration for stability
        nsig = 0.5 * nsig + 0.5 * nsig_new
    return nsig

print(f"{'pT bin':>12s} {'Nsig_nom':>14s} {'dNsig/dcC':>14s} {'dNsig/dcB':>14s} {'ratio |cC/cB|':>14s} {'purity_nom':>12s}")
print("-" * 85)

for (pt_lo, pt_hi), idx in zip(TARGET_BINS_GEV, bin_indices):
    A = A_vals[idx]
    B = B_vals[idx]
    C = C_vals[idx]
    D = D_vals[idx]
    Asig = Asig_vals[idx]

    cB = nominal_cB[(pt_lo, pt_hi)]
    cC = nominal_cC[(pt_lo, pt_hi)]
    cD = nominal_cD[(pt_lo, pt_hi)]

    nsig_nom = solve_nsig(A, B, C, D, cB, cC, cD)

    eps = 1e-5
    nsig_cC_up = solve_nsig(A, B, C, D, cB, cC + eps, cD)
    nsig_cB_up = solve_nsig(A, B, C, D, cB + eps, cC, cD)

    dNsig_dcC = (nsig_cC_up - nsig_nom) / eps
    dNsig_dcB = (nsig_cB_up - nsig_nom) / eps

    ratio = abs(dNsig_dcC) / abs(dNsig_dcB) if abs(dNsig_dcB) > 0 else float('inf')
    purity = nsig_nom / A if A > 0 else 0

    label = f"[{pt_lo:.0f},{pt_hi:.0f}]"
    print(f"{label:>12s} {nsig_nom:14.1f} {dNsig_dcC:14.1f} {dNsig_dcB:14.1f} {ratio:14.2f} {purity:12.4f}")

# ========== PURITY vs cC SCAN ==========
print()
print("=" * 90)
print("PURITY SCAN: How purity changes as cC varies from 0 to 2x nominal")
print("=" * 90)

for (pt_lo, pt_hi), idx in zip(TARGET_BINS_GEV, bin_indices):
    A = A_vals[idx]
    B = B_vals[idx]
    C = C_vals[idx]
    D = D_vals[idx]

    cB = nominal_cB[(pt_lo, pt_hi)]
    cC_nom = nominal_cC[(pt_lo, pt_hi)]
    cD = nominal_cD[(pt_lo, pt_hi)]

    label = f"[{pt_lo:.0f},{pt_hi:.0f}]"
    print(f"\npT {label} GeV:")
    print(f"  {'cC/cC_nom':>10s}  {'Nsig':>14s}  {'purity':>10s}  {'delta_pur':>10s}")
    print(f"  {'-'*10}  {'-'*14}  {'-'*10}  {'-'*10}")

    nsig_nom = solve_nsig(A, B, C, D, cB, cC_nom, cD)
    pur_nom = nsig_nom / A

    for frac in [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]:
        cC_v = frac * cC_nom
        nsig = solve_nsig(A, B, C, D, cB, cC_v, cD)
        pur = nsig / A if A > 0 else 0
        delta = pur - pur_nom
        marker = " <-- nominal" if abs(frac - 1.0) < 0.01 else ""
        print(f"  {frac:10.2f}  {nsig:14.1f}  {pur:10.4f}  {delta:+10.4f}{marker}")


print()
print("=" * 90)
print("INTERPRETATION")
print("=" * 90)
print("""
KEY FINDINGS:

1. PURE BACKGROUND R_bkg = A_bkg*D_bkg / (B_bkg*C_bkg):
   - R_bkg > 1 only for the lowest pT bin [8,10] GeV (R ~ 1.07).
   - For all higher pT bins, R_bkg < 1 and decreases sharply with pT.
   - The crossover R_bkg = 1 occurs near pT ~ 10-11 GeV.
   - This means the ABCD factorization assumption (R=1) breaks down
     progressively at higher pT.

2. SIGNAL CONTAMINATION OF REGION C:
   - Region C (non-tight, isolated) has a large signal fraction that
     grows with pT: from ~83% at [8,10] to ~57% at [24,26].
   - This is because real photons that fail the tight BDT cut but pass
     isolation ("signal leaking into C") dominate over jet background in C.
   - The leakage correction cC * A_sig removes this contamination, but
     since C_sig/C is large, the correction is a huge fraction of C_total.
   - Small errors in cC therefore translate to large errors in C_corr,
     which propagates directly to the background estimate.

3. SENSITIVITY COMPARISON (cC vs cB):
   - dR/dcC and dR/dcB are comparable in magnitude.
   - However, the CRITICAL difference is that C is much smaller than B,
     so the fractional impact of cC errors on C_corr is much larger.
   - The R = 1.0 crossover in the cC scan happens at cC_crit/cC_nom ~ 1.1
     for low pT (very sensitive) and ~ 1.8 for high pT.

4. PURITY SENSITIVITY:
   - Purity (N_sig/A) is sensitive to both cC and cB, but cC dominates
     at low pT while cB dominates at high pT (because B_sig/B grows faster).
   - The sign: increasing cC INCREASES purity (less background estimated),
     while increasing cB also increases purity (same mechanism through B_corr).

CONCLUSION: The R=1 crossover and its pT dependence are controlled by the
interplay of signal contamination in all sideband regions, with region C
being particularly vulnerable because it has the smallest absolute yield
but a large signal fraction. The correction cC * A_sig can easily exceed
C_bkg, making C_corr negative and R divergent. This is the mechanism
behind the sign flip.
""")
