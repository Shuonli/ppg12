#!/usr/bin/env python3
"""
R-factor sensitivity to signal leakage fractions cC and cB.

Reads MC ABCD yields from MC_efficiency_bdt_nom.root, computes leakage
fractions, then scans cC (and cB) to find where R crosses 1.0.

Physics:
  Pure background R = (A*D)/(B*C) is always > 1.
  Signal leakage corrections shift R toward (and below) 1.0.
  The corrected R formula is:
      R = A * D_corr / (B_corr * C_corr)
  where B_corr = B - cB * A_sig,  C_corr = C - cC * A_sig,  D_corr = D - cD * A_sig

  We scan cC (and separately cB) to identify which leakage fraction
  controls the crossover R = 1.
"""

import numpy as np
import uproot

# ── Configuration ──────────────────────────────────────────────────────
ROOT_FILE = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root"
PT_EDGES = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
# Bin indices (0-based) for the 4 representative pT bins
# [8,10]=0, [12,14]=2, [16,18]=4, [24,26]=8
REPRESENTATIVE_BINS = {
    "[8,10]":  0,
    "[12,14]": 2,
    "[16,18]": 4,
    "[24,26]": 8,
}
N_SCAN = 21  # 0 to 2x nominal in 21 steps (inclusive)

# ── Read histograms ───────────────────────────────────────────────────
f = uproot.open(ROOT_FILE)

# Total ABCD yields (signal + background MC combined)
A_tot = f["h_tight_iso_cluster_0"].values()
B_tot = f["h_tight_noniso_cluster_0"].values()
C_tot = f["h_nontight_iso_cluster_0"].values()
D_tot = f["h_nontight_noniso_cluster_0"].values()

# Signal-only yields (truth-matched photon)
A_sig = f["h_tight_iso_cluster_signal_0"].values()
B_sig = f["h_tight_noniso_cluster_signal_0"].values()
C_sig = f["h_nontight_iso_cluster_signal_0"].values()
D_sig = f["h_nontight_noniso_cluster_signal_0"].values()

# ── Compute nominal leakage fractions ─────────────────────────────────
# cX = signal_X / signal_A  (fraction of signal in A that leaks into X)
cB_nom = B_sig / A_sig
cC_nom = C_sig / A_sig
cD_nom = D_sig / A_sig

# ── Pure (uncorrected) R ──────────────────────────────────────────────
R_pure = (A_tot * D_tot) / (B_tot * C_tot)

# ── Corrected R at nominal leakage ───────────────────────────────────
def compute_R(A, B, C, D, cB, cC, cD, A_sig):
    """Compute corrected R = A * D_corr / (B_corr * C_corr)."""
    B_corr = B - cB * A_sig
    C_corr = C - cC * A_sig
    D_corr = D - cD * A_sig
    # Guard against division by zero
    denom = B_corr * C_corr
    with np.errstate(divide='ignore', invalid='ignore'):
        R = np.where(np.abs(denom) > 0, A * D_corr / denom, np.nan)
    return R

R_corrected = compute_R(A_tot, B_tot, C_tot, D_tot, cB_nom, cC_nom, cD_nom, A_sig)

# ── Print overview table ──────────────────────────────────────────────
print("=" * 110)
print("NOMINAL ABCD YIELDS AND LEAKAGE FRACTIONS (all 12 pT bins)")
print("=" * 110)
print(f"{'pT bin':>10s} {'A_tot':>12s} {'B_tot':>12s} {'C_tot':>12s} {'D_tot':>12s}"
      f" {'cB':>8s} {'cC':>8s} {'cD':>8s} {'R_pure':>8s} {'R_corr':>8s}")
print("-" * 110)
for i in range(len(A_tot)):
    lo, hi = PT_EDGES[i], PT_EDGES[i + 1]
    print(f"[{lo:.0f},{hi:.0f}]".rjust(10),
          f"{A_tot[i]:12.1f}", f"{B_tot[i]:12.1f}",
          f"{C_tot[i]:12.1f}", f"{D_tot[i]:12.1f}",
          f"{cB_nom[i]:8.4f}", f"{cC_nom[i]:8.4f}", f"{cD_nom[i]:8.4f}",
          f"{R_pure[i]:8.4f}", f"{R_corrected[i]:8.4f}")

# ── Scan functions ────────────────────────────────────────────────────
def scan_leakage(ibin, scan_which="cC"):
    """
    Scan one leakage fraction from 0 to 2x nominal, keeping others at nominal.
    Returns (scan_values, R_values).
    """
    A = A_tot[ibin]
    B = B_tot[ibin]
    C = C_tot[ibin]
    D = D_tot[ibin]
    Asig = A_sig[ibin]
    cb = cB_nom[ibin]
    cc = cC_nom[ibin]
    cd = cD_nom[ibin]

    if scan_which == "cC":
        nom_val = cc
    elif scan_which == "cB":
        nom_val = cb
    else:
        raise ValueError(f"Unknown scan_which={scan_which}")

    scan_vals = np.linspace(0, 2 * nom_val, N_SCAN)
    R_vals = np.empty(N_SCAN)

    for j, sv in enumerate(scan_vals):
        if scan_which == "cC":
            B_corr = B - cb * Asig
            C_corr = C - sv * Asig
            D_corr = D - cd * Asig
        else:  # cB
            B_corr = B - sv * Asig
            C_corr = C - cc * Asig
            D_corr = D - cd * Asig

        denom = B_corr * C_corr
        if abs(denom) > 0:
            R_vals[j] = A * D_corr / denom
        else:
            R_vals[j] = np.nan

    return scan_vals, R_vals


def find_crossing(scan_vals, R_vals, target=1.0):
    """Find the scan value where R crosses target by linear interpolation."""
    for j in range(len(R_vals) - 1):
        r0, r1 = R_vals[j], R_vals[j + 1]
        if np.isnan(r0) or np.isnan(r1):
            continue
        if (r0 - target) * (r1 - target) <= 0:
            # Linear interpolation
            frac = (target - r0) / (r1 - r0) if (r1 - r0) != 0 else 0
            return scan_vals[j] + frac * (scan_vals[j + 1] - scan_vals[j])
    return np.nan


# ── Run scans ─────────────────────────────────────────────────────────
print("\n")
print("=" * 110)
print("SCAN OF cC: R vs cC (cB and cD held at nominal)")
print("=" * 110)

cC_crossings = {}
for label, ibin in REPRESENTATIVE_BINS.items():
    scan_vals, R_vals = scan_leakage(ibin, "cC")
    crossing = find_crossing(scan_vals, R_vals, 1.0)
    cC_crossings[label] = crossing

    print(f"\npT bin {label} GeV  (nominal cC = {cC_nom[ibin]:.6f})")
    print(f"  {'cC':>10s}  {'cC/cC_nom':>10s}  {'R':>10s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*10}")
    for j in range(N_SCAN):
        marker = " <-- R=1 crossing" if not np.isnan(crossing) and j > 0 and \
            (R_vals[j-1] - 1.0) * (R_vals[j] - 1.0) <= 0 else ""
        ratio = scan_vals[j] / cC_nom[ibin] if cC_nom[ibin] > 0 else np.nan
        print(f"  {scan_vals[j]:10.6f}  {ratio:10.4f}  {R_vals[j]:10.6f}{marker}")

    if not np.isnan(crossing):
        print(f"  --> R = 1.0 crossing at cC = {crossing:.6f} "
              f"({crossing / cC_nom[ibin]:.4f}x nominal)")
    else:
        print(f"  --> R never crosses 1.0 in scan range")

print("\n")
print("=" * 110)
print("SCAN OF cB: R vs cB (cC and cD held at nominal)")
print("=" * 110)

cB_crossings = {}
for label, ibin in REPRESENTATIVE_BINS.items():
    scan_vals, R_vals = scan_leakage(ibin, "cB")
    crossing = find_crossing(scan_vals, R_vals, 1.0)
    cB_crossings[label] = crossing

    print(f"\npT bin {label} GeV  (nominal cB = {cB_nom[ibin]:.6f})")
    print(f"  {'cB':>10s}  {'cB/cB_nom':>10s}  {'R':>10s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*10}")
    for j in range(N_SCAN):
        marker = " <-- R=1 crossing" if not np.isnan(crossing) and j > 0 and \
            (R_vals[j-1] - 1.0) * (R_vals[j] - 1.0) <= 0 else ""
        ratio = scan_vals[j] / cB_nom[ibin] if cB_nom[ibin] > 0 else np.nan
        print(f"  {scan_vals[j]:10.6f}  {ratio:10.4f}  {R_vals[j]:10.6f}{marker}")

    if not np.isnan(crossing):
        print(f"  --> R = 1.0 crossing at cB = {crossing:.6f} "
              f"({crossing / cB_nom[ibin]:.4f}x nominal)")
    else:
        print(f"  --> R never crosses 1.0 in scan range")

# ── Sensitivity comparison ────────────────────────────────────────────
print("\n")
print("=" * 110)
print("SENSITIVITY COMPARISON: dR/d(cX) at nominal, and R=1 crossings")
print("=" * 110)

print(f"\n{'pT bin':>10s} {'cB_nom':>8s} {'cC_nom':>8s} {'cD_nom':>8s}"
      f" {'R_pure':>8s} {'R_corr':>8s}"
      f" {'dR/dcC':>10s} {'dR/dcB':>10s}"
      f" {'cC@R=1':>10s} {'cB@R=1':>10s}"
      f" {'cC@R=1/nom':>11s} {'cB@R=1/nom':>11s}")
print("-" * 110)

for label, ibin in REPRESENTATIVE_BINS.items():
    # Numerical derivative dR/dcC at nominal
    eps = cC_nom[ibin] * 0.01 if cC_nom[ibin] > 0 else 1e-6
    _, R_plus = scan_leakage(ibin, "cC")
    # Compute R at nominal +/- eps
    A = A_tot[ibin]; B = B_tot[ibin]; C = C_tot[ibin]; D = D_tot[ibin]; Asig = A_sig[ibin]
    cb = cB_nom[ibin]; cc = cC_nom[ibin]; cd = cD_nom[ibin]

    def R_at(cB_val, cC_val, cD_val):
        B_corr = B - cB_val * Asig
        C_corr = C - cC_val * Asig
        D_corr = D - cD_val * Asig
        d = B_corr * C_corr
        return A * D_corr / d if abs(d) > 0 else np.nan

    dRdcC = (R_at(cb, cc + eps, cd) - R_at(cb, cc - eps, cd)) / (2 * eps)
    dRdcB = (R_at(cb + eps, cc, cd) - R_at(cb - eps, cc, cd)) / (2 * eps)

    cC_cross = cC_crossings[label]
    cB_cross = cB_crossings[label]

    cC_ratio = cC_cross / cc if cc > 0 and not np.isnan(cC_cross) else np.nan
    cB_ratio = cB_cross / cb if cb > 0 and not np.isnan(cB_cross) else np.nan

    print(f"{label:>10s} {cb:8.4f} {cc:8.4f} {cd:8.4f}"
          f" {R_pure[ibin]:8.4f} {R_corrected[ibin]:8.4f}"
          f" {dRdcC:10.4f} {dRdcB:10.4f}",
          end="")

    if not np.isnan(cC_cross):
        print(f" {cC_cross:10.6f}", end="")
    else:
        print(f" {'N/A':>10s}", end="")
    if not np.isnan(cB_cross):
        print(f" {cB_cross:10.6f}", end="")
    else:
        print(f" {'N/A':>10s}", end="")
    if not np.isnan(cC_ratio):
        print(f" {cC_ratio:11.4f}", end="")
    else:
        print(f" {'N/A':>11s}", end="")
    if not np.isnan(cB_ratio):
        print(f" {cB_ratio:11.4f}", end="")
    else:
        print(f" {'N/A':>11s}", end="")
    print()

# ── Interpretation ────────────────────────────────────────────────────
print("\n")
print("=" * 110)
print("INTERPRETATION")
print("=" * 110)

print("""
The R-factor formula (with leakage corrections) is:

    R = A * (D - cD * A_sig) / [(B - cB * A_sig) * (C - cC * A_sig)]

Pure (uncorrected) R = A*D / (B*C) is always > 1 because tight+iso (A)
is enriched in signal, inflating the numerator relative to the denominator.

Signal leakage corrections reduce R toward 1.0:
  - cB subtracts signal from B (denominator increases -> R decreases)
  - cC subtracts signal from C (denominator increases -> R decreases)
  - cD subtracts signal from D (numerator decreases -> R decreases)

The key question: which correction drives R below 1.0?
""")

for label, ibin in REPRESENTATIVE_BINS.items():
    cc_cross = cC_crossings[label]
    cb_cross = cB_crossings[label]
    cc = cC_nom[ibin]
    cb = cB_nom[ibin]

    print(f"  pT {label} GeV:")
    print(f"    R_pure = {R_pure[ibin]:.4f},  R_corrected = {R_corrected[ibin]:.4f}")
    if not np.isnan(cc_cross):
        print(f"    cC scan: R crosses 1.0 at cC = {cc_cross:.6f} = {cc_cross/cc:.3f}x nominal")
    else:
        print(f"    cC scan: R stays {'above' if R_corrected[ibin] > 1 else 'below'} 1.0 for all cC in [0, 2x nominal]")
    if not np.isnan(cb_cross):
        print(f"    cB scan: R crosses 1.0 at cB = {cb_cross:.6f} = {cb_cross/cb:.3f}x nominal")
    else:
        print(f"    cB scan: R stays {'above' if R_corrected[ibin] > 1 else 'below'} 1.0 for all cB in [0, 2x nominal]")
    print()

print("CONCLUSION:")
print("-" * 80)
# Determine which has crossings closer to nominal
has_cC_crossings = sum(1 for v in cC_crossings.values() if not np.isnan(v))
has_cB_crossings = sum(1 for v in cB_crossings.values() if not np.isnan(v))
print(f"  cC scan produces R=1 crossings in {has_cC_crossings}/{len(REPRESENTATIVE_BINS)} bins")
print(f"  cB scan produces R=1 crossings in {has_cB_crossings}/{len(REPRESENTATIVE_BINS)} bins")

# Compare |dR/dcC| vs |dR/dcB| at nominal
for label, ibin in REPRESENTATIVE_BINS.items():
    A = A_tot[ibin]; B = B_tot[ibin]; C = C_tot[ibin]; D = D_tot[ibin]; Asig = A_sig[ibin]
    cb = cB_nom[ibin]; cc = cC_nom[ibin]; cd = cD_nom[ibin]
    eps_c = cc * 0.01 if cc > 0 else 1e-6
    eps_b = cb * 0.01 if cb > 0 else 1e-6

    def R_at_local(cB_val, cC_val, cD_val):
        B_corr = B - cB_val * Asig
        C_corr = C - cC_val * Asig
        D_corr = D - cD_val * Asig
        d = B_corr * C_corr
        return A * D_corr / d if abs(d) > 0 else np.nan

    dRdcC = (R_at_local(cb, cc + eps_c, cd) - R_at_local(cb, cc - eps_c, cd)) / (2 * eps_c)
    dRdcB = (R_at_local(cb + eps_b, cc, cd) - R_at_local(cb - eps_b, cc, cd)) / (2 * eps_b)
    dominant = "cC" if abs(dRdcC) > abs(dRdcB) else "cB"
    ratio = abs(dRdcC) / abs(dRdcB) if abs(dRdcB) > 0 else float('inf')
    print(f"  pT {label}: |dR/dcC|/|dR/dcB| = {ratio:.2f}  -->  {dominant} dominates")
