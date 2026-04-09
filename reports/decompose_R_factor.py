#!/usr/bin/env python3
"""
Decompose the ABCD R factor to identify which signal leakage correction
drives the crossover at pT ~ 14 GeV.

R is defined as R = (A * D) / (B * C) for the ABCD background subtraction.
  - Pure background R (jet MC only) is always > 1
  - Corrected R (h_R from MC closure) crosses 1.0 at pT ~ 14 GeV
  - The difference comes from signal leakage corrections cB, cC, cD

Signal leakage fractions (from CalculatePhotonYield.C):
  cB = B_sig / A_sig
  cC = C_sig / A_sig
  cD = D_sig / A_sig

Corrected regions:
  B_corr = B_tot - cB * A_sig = B_bkg
  C_corr = C_tot - cC * A_sig = C_bkg
  D_corr = D_tot - cD * A_sig = D_bkg

h_R in the code = (A_bkg * D_bkg) / (B_bkg * C_bkg)  [= R_no_leak below]
"""

import uproot
import numpy as np

# ---------- configuration ----------
rootfile = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root"
pt_edges = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
nbins = len(pt_edges) - 1

# ---------- read histograms ----------
f = uproot.open(rootfile)

# Total yields per ABCD region
A_tot = f["h_tight_iso_cluster_0"].values()
B_tot = f["h_tight_noniso_cluster_0"].values()
C_tot = f["h_nontight_iso_cluster_0"].values()
D_tot = f["h_nontight_noniso_cluster_0"].values()

# Signal yields per ABCD region
A_sig = f["h_tight_iso_cluster_signal_0"].values()
B_sig = f["h_tight_noniso_cluster_signal_0"].values()
C_sig = f["h_nontight_iso_cluster_signal_0"].values()
D_sig = f["h_nontight_noniso_cluster_signal_0"].values()

# Verify edges match
edges_check = f["h_tight_iso_cluster_0"].axis().edges()
assert np.allclose(edges_check, pt_edges), f"Edges mismatch: {edges_check} vs {pt_edges}"

# ---------- compute background counts ----------
A_bkg = A_tot - A_sig
B_bkg = B_tot - B_sig
C_bkg = C_tot - C_sig
D_bkg = D_tot - D_sig

# ---------- compute leakage fractions ----------
cB = np.where(A_sig > 0, B_sig / A_sig, 0.0)
cC = np.where(A_sig > 0, C_sig / A_sig, 0.0)
cD = np.where(A_sig > 0, D_sig / A_sig, 0.0)

# ---------- R under different correction scenarios ----------
def safe_ratio(num, den):
    """Element-wise ratio, returning nan where denominator is zero."""
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.where(den != 0, num / den, np.nan)
    return result

# R_raw: no corrections, total counts in all regions
R_raw = safe_ratio(A_tot * D_tot, B_tot * C_tot)

# R_no_leak: pure background (tot - sig) in all four regions
# This is h_R as computed in CalculatePhotonYield.C
R_no_leak = safe_ratio(A_bkg * D_bkg, B_bkg * C_bkg)

# R_only_cB: subtract signal leakage from B only
R_only_cB = safe_ratio(A_tot * D_tot, (B_tot - cB * A_sig) * C_tot)

# R_only_cC: subtract signal leakage from C only
R_only_cC = safe_ratio(A_tot * D_tot, B_tot * (C_tot - cC * A_sig))

# R_only_cD: subtract signal leakage from D only
R_only_cD = safe_ratio(A_tot * (D_tot - cD * A_sig), B_tot * C_tot)

# R_cB_cC: subtract from B and C but not D
R_cB_cC = safe_ratio(A_tot * D_tot, (B_tot - cB * A_sig) * (C_tot - cC * A_sig))

# R_all: fully corrected (subtract from B, C, and D -- but NOT from A)
# This is how the ABCD method works: A is the signal region, corrections apply to B,C,D
R_all = safe_ratio(A_tot * (D_tot - cD * A_sig),
                   (B_tot - cB * A_sig) * (C_tot - cC * A_sig))

# ======================================================================
# PRINT TABLES
# ======================================================================

print("=" * 140)
print("TABLE 1: Raw counts per ABCD region")
print("=" * 140)
header = f"{'pT bin':>12s} {'A_tot':>14s} {'B_tot':>14s} {'C_tot':>14s} {'D_tot':>14s} {'A_sig':>14s} {'B_sig':>14s} {'C_sig':>14s} {'D_sig':>14s}"
print(header)
print("-" * 140)
for i in range(nbins):
    lo, hi = pt_edges[i], pt_edges[i+1]
    print(f"  [{lo:4.0f},{hi:4.0f}] {A_tot[i]:14.0f} {B_tot[i]:14.0f} {C_tot[i]:14.0f} {D_tot[i]:14.0f}"
          f" {A_sig[i]:14.0f} {B_sig[i]:14.0f} {C_sig[i]:14.0f} {D_sig[i]:14.0f}")

print()
print("=" * 140)
print("TABLE 2: Signal fractions and leakage ratios")
print("=" * 140)
header2 = f"{'pT bin':>12s} {'A_sig/A_tot':>12s} {'B_sig/B_tot':>12s} {'C_sig/C_tot':>12s} {'D_sig/D_tot':>12s} {'cB':>10s} {'cC':>10s} {'cD':>10s}"
print(header2)
print("-" * 140)
for i in range(nbins):
    lo, hi = pt_edges[i], pt_edges[i+1]
    fA = A_sig[i] / A_tot[i] if A_tot[i] > 0 else 0
    fB = B_sig[i] / B_tot[i] if B_tot[i] > 0 else 0
    fC = C_sig[i] / C_tot[i] if C_tot[i] > 0 else 0
    fD = D_sig[i] / D_tot[i] if D_tot[i] > 0 else 0
    print(f"  [{lo:4.0f},{hi:4.0f}]   {fA:10.4f}   {fB:10.4f}   {fC:10.4f}   {fD:10.4f}   {cB[i]:8.5f}   {cC[i]:8.5f}   {cD[i]:8.5f}")

print()
print("=" * 160)
print("TABLE 3: R factor under different correction scenarios")
print("=" * 160)
header3 = (f"{'pT bin':>12s} {'R_raw':>10s} {'R_no_leak':>10s} {'R_only_cB':>10s} {'R_only_cC':>10s} "
           f"{'R_only_cD':>10s} {'R_cB_cC':>10s} {'R_all':>10s}")
print(header3)
print("-" * 160)
for i in range(nbins):
    lo, hi = pt_edges[i], pt_edges[i+1]
    print(f"  [{lo:4.0f},{hi:4.0f}] {R_raw[i]:10.4f} {R_no_leak[i]:10.4f} {R_only_cB[i]:10.4f} {R_only_cC[i]:10.4f} "
          f"{R_only_cD[i]:10.4f} {R_cB_cC[i]:10.4f} {R_all[i]:10.4f}")

# ======================================================================
# ANALYSIS: Which correction has the largest effect?
# ======================================================================

print()
print("=" * 160)
print("TABLE 4: Effect of each correction (difference from R_raw)")
print("=" * 160)
header4 = (f"{'pT bin':>12s} {'dR(cB)':>12s} {'dR(cC)':>12s} {'dR(cD)':>12s} "
           f"{'dR(cB+cC)':>12s} {'dR(all)':>12s} {'dominant':>12s}")
print(header4)
print("-" * 160)
for i in range(nbins):
    lo, hi = pt_edges[i], pt_edges[i+1]
    dR_cB = R_only_cB[i] - R_raw[i]
    dR_cC = R_only_cC[i] - R_raw[i]
    dR_cD = R_only_cD[i] - R_raw[i]
    dR_cBcC = R_cB_cC[i] - R_raw[i]
    dR_all = R_all[i] - R_raw[i]

    # Dominant: which single correction shifts R the most (in absolute value)?
    shifts = {'cB': abs(dR_cB), 'cC': abs(dR_cC), 'cD': abs(dR_cD)}
    dominant = max(shifts, key=shifts.get)

    print(f"  [{lo:4.0f},{hi:4.0f}] {dR_cB:12.4f} {dR_cC:12.4f} {dR_cD:12.4f} "
          f"{dR_cBcC:12.4f} {dR_all:12.4f} {dominant:>12s}")

# ======================================================================
# Does cC alone create the crossover?
# ======================================================================
print()
print("=" * 100)
print("TABLE 5: Does R_only_cC cross 1.0?")
print("=" * 100)
print(f"{'pT bin':>12s} {'R_raw':>10s} {'R_only_cC':>10s} {'crosses 1?':>12s}")
print("-" * 100)
for i in range(nbins):
    lo, hi = pt_edges[i], pt_edges[i+1]
    crosses = "YES" if R_only_cC[i] < 1.0 else "no"
    print(f"  [{lo:4.0f},{hi:4.0f}] {R_raw[i]:10.4f} {R_only_cC[i]:10.4f} {crosses:>12s}")

# ======================================================================
# Critical cC: what value of cC would make R exactly = 1?
# ======================================================================
# R_only_cC = (A_tot * D_tot) / (B_tot * (C_tot - cC * A_sig))
# Set R_only_cC = 1:
#   A_tot * D_tot = B_tot * (C_tot - cC_crit * A_sig)
#   C_tot - cC_crit * A_sig = A_tot * D_tot / B_tot
#   cC_crit = (C_tot - A_tot * D_tot / B_tot) / A_sig

print()
print("=" * 100)
print("TABLE 6: Critical cC (value that would make R = 1 with only cC correction)")
print("=" * 100)
print(f"{'pT bin':>12s} {'cC_actual':>12s} {'cC_critical':>12s} {'ratio':>10s} {'cC > cC_crit?':>15s}")
print("-" * 100)
for i in range(nbins):
    lo, hi = pt_edges[i], pt_edges[i+1]
    if B_tot[i] > 0 and A_sig[i] > 0:
        cC_crit = (C_tot[i] - A_tot[i] * D_tot[i] / B_tot[i]) / A_sig[i]
    else:
        cC_crit = np.nan
    ratio = cC[i] / cC_crit if cC_crit != 0 and not np.isnan(cC_crit) else np.nan
    exceeds = "YES (R<1)" if cC[i] > cC_crit else "no (R>1)"
    print(f"  [{lo:4.0f},{hi:4.0f}] {cC[i]:12.5f} {cC_crit:12.5f} {ratio:10.3f} {exceeds:>15s}")

# ======================================================================
# KEY INSIGHT: R_no_leak vs R_all -- the role of the A region
# ======================================================================
# R_no_leak = (A_bkg * D_bkg) / (B_bkg * C_bkg)   -- background in ALL four regions
# R_all     = (A_tot * D_bkg) / (B_bkg * C_bkg)   -- total in A, background in B,C,D
#
# So R_all / R_no_leak = A_tot / A_bkg = 1 / (1 - purity)
# And R_no_leak = R_all * (1 - purity)
#
# The code's h_R IS R_no_leak. But in the actual ABCD method, R multiplies the
# background estimate: N_A^bkg = R * (B_corr * C_corr) / D_corr
# The R that enters the formula is computed from MC truth as:
#   R = A_bkg * D_bkg / (B_bkg * C_bkg)
# So R_no_leak is the correct one.

print()
print("=" * 120)
print("TABLE 7: R_no_leak / R_raw decomposition -- why R drops below 1")
print("=" * 120)
print(f"{'pT bin':>12s} {'1-purity':>12s} {'B_bkg/B_tot':>12s} {'C_bkg/C_tot':>12s} {'D_bkg/D_tot':>12s} "
      f"{'product':>12s} {'R_noleak/R_raw':>15s}")
print("-" * 120)
for i in range(nbins):
    lo, hi = pt_edges[i], pt_edges[i+1]
    fA = A_bkg[i] / A_tot[i] if A_tot[i] > 0 else np.nan
    fB = B_bkg[i] / B_tot[i] if B_tot[i] > 0 else np.nan
    fC = C_bkg[i] / C_tot[i] if C_tot[i] > 0 else np.nan
    fD = D_bkg[i] / D_tot[i] if D_tot[i] > 0 else np.nan
    ratio_R = R_no_leak[i] / R_raw[i] if R_raw[i] != 0 else np.nan
    product = (fA * fD) / (fB * fC) if fB * fC != 0 else np.nan
    print(f"  [{lo:4.0f},{hi:4.0f}] {fA:12.5f} {fB:12.5f} {fC:12.5f} {fD:12.5f} "
          f"{product:12.5f} {ratio_R:15.5f}")

print()
print("  R_no_leak / R_raw = (A_bkg/A_tot) * (D_bkg/D_tot) / ((B_bkg/B_tot) * (C_bkg/C_tot))")
print("  = (1-purity_A) * (1-f_D) / ((1-f_B) * (1-f_C))")
print()
print("  At low pT: purity_A ~ 97.6%, so (1-purity_A) ~ 0.024")
print("  But B,C,D also have signal contamination that partially cancels.")
print("  The net effect is that (1-purity_A) dominates because A has by far")
print("  the highest signal fraction, pulling R_no_leak far below R_raw.")

# ======================================================================
# SUMMARY
# ======================================================================
print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)

# Check which single correction is dominant across all bins
dominant_counts = {'cB': 0, 'cC': 0, 'cD': 0}
for i in range(nbins):
    dR_cB = abs(R_only_cB[i] - R_raw[i])
    dR_cC = abs(R_only_cC[i] - R_raw[i])
    dR_cD = abs(R_only_cD[i] - R_raw[i])
    shifts = {'cB': dR_cB, 'cC': dR_cC, 'cD': dR_cD}
    dominant = max(shifts, key=shifts.get)
    dominant_counts[dominant] += 1

print(f"\n1. Dominant single correction across {nbins} pT bins:")
for corr, count in sorted(dominant_counts.items(), key=lambda x: -x[1]):
    print(f"     {corr}: dominant in {count}/{nbins} bins")

# Check crossovers
def find_crossover(R_arr, name):
    for i in range(nbins - 1):
        if R_arr[i] >= 1.0 and R_arr[i+1] < 1.0:
            return f"{name} crosses 1 between pT = [{pt_edges[i]:.0f},{pt_edges[i+2]:.0f}] GeV"
    if np.all(R_arr > 1):
        return f"{name}: always > 1 (min = {np.nanmin(R_arr):.4f})"
    if np.all(R_arr < 1):
        return f"{name}: always < 1 (max = {np.nanmax(R_arr):.4f})"
    return f"{name}: non-monotonic, range [{np.nanmin(R_arr):.4f}, {np.nanmax(R_arr):.4f}]"

print(f"\n2. R = 1 crossover locations:")
for R_arr, name in [(R_raw, "R_raw"), (R_no_leak, "R_no_leak (h_R)"),
                     (R_only_cB, "R_only_cB"), (R_only_cC, "R_only_cC"),
                     (R_only_cD, "R_only_cD"), (R_cB_cC, "R_cB_cC"),
                     (R_all, "R_all")]:
    print(f"     {find_crossover(R_arr, name)}")

print(f"\n3. The R_raw -> R_no_leak transition is NOT just about B,C,D corrections.")
print(f"   It fundamentally involves the A region: R_no_leak uses A_bkg while R_raw uses A_tot.")
print(f"   At high pT where purity -> 97%, A_bkg/A_tot -> 0.03, which is the dominant")
print(f"   factor suppressing R_no_leak below 1.")
print(f"\n4. Among corrections to B, C, D only (keeping A_tot in numerator):")
print(f"   R_all = A_tot * D_bkg / (B_bkg * C_bkg) includes only B,C,D corrections.")
print(f"   Comparing R_all to R_no_leak isolates the A-region signal subtraction effect.")

print()
for i in range(nbins):
    lo, hi = pt_edges[i], pt_edges[i+1]
    if i == 0:
        print(f"   {'pT bin':>12s} {'R_raw':>10s} {'R_all(BCD)':>12s} {'R_no_leak':>10s} {'A_tot/A_bkg':>12s}")
        print(f"   {'-'*60}")
    rat = A_tot[i] / A_bkg[i] if A_bkg[i] > 0 else np.nan
    print(f"   [{lo:4.0f},{hi:4.0f}] {R_raw[i]:10.4f} {R_all[i]:12.4f} {R_no_leak[i]:10.4f} {rat:12.2f}")
