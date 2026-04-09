#!/usr/bin/env python3
"""
Extract per-pT-bin ABCD yields from MC files to decompose the R factor.
Reads merged (sig+bkg), jet-only, and MC closure files.
"""

import uproot
import numpy as np

# --- File paths ---
merged_file = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root"
jet_file = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_bdt_nom.root"
closure_file = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom_mc.root"

# pT bins
pt_edges = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
n_bins = len(pt_edges) - 1  # 12 bins

# --- Histogram names to read ---
regions = ["tight_iso", "tight_noniso", "nontight_iso", "nontight_noniso"]
region_labels = {"tight_iso": "A", "tight_noniso": "B", "nontight_iso": "C", "nontight_noniso": "D"}

# Build exact histogram name lists
merged_total_names = [f"h_{r}_cluster_0" for r in regions]
merged_signal_names = [f"h_{r}_cluster_signal_0" for r in regions]
merged_bkg_names = [f"h_{r}_cluster_background_0" for r in regions]
jet_names = [f"h_{r}_cluster_0" for r in regions]

# --- Read histograms by exact name ---
def read_hist(f, name):
    """Read a TH1D and return bin contents (excluding under/overflow)."""
    h = f[name]
    vals = h.values()  # bin contents, no overflow
    return vals

# Read merged file
print("Reading merged MC file...")
with uproot.open(merged_file) as f:
    total = {}
    signal = {}
    bkg = {}
    for r in regions:
        label = region_labels[r]
        total[label] = read_hist(f, f"h_{r}_cluster_0")
        signal[label] = read_hist(f, f"h_{r}_cluster_signal_0")
        bkg[label] = read_hist(f, f"h_{r}_cluster_background_0")

# Read jet-only file
print("Reading jet-only MC file...")
with uproot.open(jet_file) as f:
    jet = {}
    for r in regions:
        label = region_labels[r]
        jet[label] = read_hist(f, f"h_{r}_cluster_0")

# Read closure file
print("Reading MC closure file...")
with uproot.open(closure_file) as f:
    h_R = read_hist(f, "h_R")
    h_leak_B = read_hist(f, "h_leak_B")
    h_leak_C = read_hist(f, "h_leak_C")
    h_leak_D = read_hist(f, "h_leak_D")

# --- Verify bin counts ---
for label in "ABCD":
    assert len(total[label]) == n_bins, f"total[{label}] has {len(total[label])} bins, expected {n_bins}"
    assert len(signal[label]) == n_bins, f"signal[{label}] has {len(signal[label])} bins, expected {n_bins}"
    assert len(bkg[label]) == n_bins, f"bkg[{label}] has {len(bkg[label])} bins, expected {n_bins}"
    assert len(jet[label]) == n_bins, f"jet[{label}] has {len(jet[label])} bins, expected {n_bins}"

assert len(h_R) == n_bins, f"h_R has {len(h_R)} bins, expected {n_bins}"
assert len(h_leak_B) == n_bins, f"h_leak_B has {len(h_leak_B)} bins, expected {n_bins}"

# --- Print comprehensive table ---
sep = "-" * 200

print("\n" + "=" * 200)
print("PPG12 R-FACTOR DECOMPOSITION: Per-pT-bin ABCD Yields (MC, nominal BDT)")
print("=" * 200)

# Table 1: Total yields (sig+bkg merged MC)
print(f"\n{'TABLE 1: TOTAL YIELDS (signal + background merged MC)':^200}")
print(sep)
header = f"{'pT bin':>12} | {'A_tot':>12} {'B_tot':>12} {'C_tot':>12} {'D_tot':>12} | {'B/D':>10} {'C/D':>10} {'BC/D':>12} {'A-BC/D':>12}"
print(header)
print(sep)
for i in range(n_bins):
    pt_label = f"[{pt_edges[i]},{pt_edges[i+1]})"
    A, B, C, D = total["A"][i], total["B"][i], total["C"][i], total["D"][i]
    BD = B / D if D > 0 else float('nan')
    CD = C / D if D > 0 else float('nan')
    BCD = (B * C) / D if D > 0 else float('nan')
    A_minus = A - BCD
    print(f"{pt_label:>12} | {A:>12.2f} {B:>12.2f} {C:>12.2f} {D:>12.2f} | {BD:>10.4f} {CD:>10.4f} {BCD:>12.2f} {A_minus:>12.2f}")
print(sep)

# Table 2: Signal yields
print(f"\n{'TABLE 2: SIGNAL-ONLY YIELDS (truth-matched photons in merged MC)':^200}")
print(sep)
header = f"{'pT bin':>12} | {'A_sig':>12} {'B_sig':>12} {'C_sig':>12} {'D_sig':>12} | {'Purity_A':>10} {'Purity_B':>10} {'Purity_C':>10} {'Purity_D':>10}"
print(header)
print(sep)
for i in range(n_bins):
    pt_label = f"[{pt_edges[i]},{pt_edges[i+1]})"
    As, Bs, Cs, Ds = signal["A"][i], signal["B"][i], signal["C"][i], signal["D"][i]
    At, Bt, Ct, Dt = total["A"][i], total["B"][i], total["C"][i], total["D"][i]
    pA = As / At if At > 0 else float('nan')
    pB = Bs / Bt if Bt > 0 else float('nan')
    pC = Cs / Ct if Ct > 0 else float('nan')
    pD = Ds / Dt if Dt > 0 else float('nan')
    print(f"{pt_label:>12} | {As:>12.2f} {Bs:>12.2f} {Cs:>12.2f} {Ds:>12.2f} | {pA:>10.4f} {pB:>10.4f} {pC:>10.4f} {pD:>10.4f}")
print(sep)

# Table 3: Background yields from merged MC
print(f"\n{'TABLE 3: BACKGROUND YIELDS (from merged MC, truth-tagged)':^200}")
print(sep)
header = f"{'pT bin':>12} | {'A_bkg':>12} {'B_bkg':>12} {'C_bkg':>12} {'D_bkg':>12} | {'B_bkg/D_bkg':>12} {'C_bkg/D_bkg':>12} {'BC_bkg/D_bkg':>14}"
print(header)
print(sep)
for i in range(n_bins):
    pt_label = f"[{pt_edges[i]},{pt_edges[i+1]})"
    Ab, Bb, Cb, Db = bkg["A"][i], bkg["B"][i], bkg["C"][i], bkg["D"][i]
    BD = Bb / Db if Db > 0 else float('nan')
    CD = Cb / Db if Db > 0 else float('nan')
    BCD = (Bb * Cb) / Db if Db > 0 else float('nan')
    print(f"{pt_label:>12} | {Ab:>12.2f} {Bb:>12.2f} {Cb:>12.2f} {Db:>12.2f} | {BD:>12.4f} {CD:>12.4f} {BCD:>14.2f}")
print(sep)

# Table 4: Jet-only MC
print(f"\n{'TABLE 4: JET-ONLY MC YIELDS (pure background, no signal photons)':^200}")
print(sep)
header = f"{'pT bin':>12} | {'A_jet':>12} {'B_jet':>12} {'C_jet':>12} {'D_jet':>12} | {'B_jet/D_jet':>12} {'C_jet/D_jet':>12} {'BC_jet/D_jet':>14}"
print(header)
print(sep)
for i in range(n_bins):
    pt_label = f"[{pt_edges[i]},{pt_edges[i+1]})"
    Aj, Bj, Cj, Dj = jet["A"][i], jet["B"][i], jet["C"][i], jet["D"][i]
    BD = Bj / Dj if Dj > 0 else float('nan')
    CD = Cj / Dj if Dj > 0 else float('nan')
    BCD = (Bj * Cj) / Dj if Dj > 0 else float('nan')
    print(f"{pt_label:>12} | {Aj:>12.2f} {Bj:>12.2f} {Cj:>12.2f} {Dj:>12.2f} | {BD:>12.4f} {CD:>12.4f} {BCD:>14.2f}")
print(sep)

# Table 5: Leakage fractions and R
print(f"\n{'TABLE 5: SIGNAL LEAKAGE FRACTIONS AND R FACTOR':^200}")
print(sep)
header = (f"{'pT bin':>12} | {'cB=Bs/As':>10} {'cC=Cs/As':>10} {'cD=Ds/As':>10} | "
          f"{'h_leak_B':>10} {'h_leak_C':>10} {'h_leak_D':>10} | "
          f"{'h_R':>10} | {'R_naive':>10} {'R_corr':>10}")
print(header)
print(sep)
for i in range(n_bins):
    pt_label = f"[{pt_edges[i]},{pt_edges[i+1]})"
    As = signal["A"][i]
    Bs = signal["B"][i]
    Cs = signal["C"][i]
    Ds = signal["D"][i]

    cB = Bs / As if As > 0 else float('nan')
    cC = Cs / As if As > 0 else float('nan')
    cD = Ds / As if As > 0 else float('nan')

    At = total["A"][i]
    Bt = total["B"][i]
    Ct = total["C"][i]
    Dt = total["D"][i]

    # Naive R = (B*C)/(A*D)
    R_naive = (Bt * Ct) / (At * Dt) if (At * Dt) > 0 else float('nan')

    # Corrected R: R = ((B - cB*A)(C - cC*A)) / (A * (D - cD*A))
    B_corr = Bt - cB * At
    C_corr = Ct - cC * At
    D_corr = Dt - cD * At
    R_corr = (B_corr * C_corr) / (At * D_corr) if (At * D_corr) > 0 else float('nan')

    print(f"{pt_label:>12} | {cB:>10.6f} {cC:>10.6f} {cD:>10.6f} | "
          f"{h_leak_B[i]:>10.6f} {h_leak_C[i]:>10.6f} {h_leak_D[i]:>10.6f} | "
          f"{h_R[i]:>10.6f} | {R_naive:>10.6f} {R_corr:>10.6f}")
print(sep)

# Table 6: ABCD closure check
print(f"\n{'TABLE 6: ABCD CLOSURE — Estimated vs True Background in A':^200}")
print(sep)
header = (f"{'pT bin':>12} | {'A_true_bkg':>12} | {'R*BC/D_naive':>14} | {'R*BC/D_corr':>14} | "
          f"{'h_R*BC/D':>14} | {'A_sig_true':>12} | {'A-Rbkg_naive':>14} {'A-Rbkg_corr':>14} {'A-Rbkg_hR':>14}")
print(header)
print(sep)
for i in range(n_bins):
    pt_label = f"[{pt_edges[i]},{pt_edges[i+1]})"
    At = total["A"][i]
    Bt = total["B"][i]
    Ct = total["C"][i]
    Dt = total["D"][i]
    As = signal["A"][i]
    Ab = bkg["A"][i]
    Bs = signal["B"][i]
    Cs = signal["C"][i]
    Ds = signal["D"][i]

    cB = Bs / As if As > 0 else 0
    cC = Cs / As if As > 0 else 0
    cD = Ds / As if As > 0 else 0

    # Naive
    BCD_naive = (Bt * Ct) / Dt if Dt > 0 else 0
    R_naive = (Bt * Ct) / (At * Dt) if (At * Dt) > 0 else 0
    est_bkg_naive = BCD_naive  # This IS R * A ... no. BCD/D * C/D... let me be precise:
    # Actually: estimated bkg in A = R * B * C / D where R=1 for naive, so just BC/D
    # For corrected: bkg_est = (B-cB*A)(C-cC*A)/(D-cD*A)

    B_corr = Bt - cB * At
    C_corr = Ct - cC * At
    D_corr = Dt - cD * At
    est_bkg_corr = (B_corr * C_corr) / D_corr if D_corr > 0 else 0

    # Using h_R from file
    est_bkg_hR = h_R[i] * (Bt * Ct) / Dt if Dt > 0 else 0

    sig_naive = At - BCD_naive
    sig_corr = At - est_bkg_corr
    sig_hR = At - est_bkg_hR

    print(f"{pt_label:>12} | {Ab:>12.2f} | {BCD_naive:>14.2f} | {est_bkg_corr:>14.2f} | "
          f"{est_bkg_hR:>14.2f} | {As:>12.2f} | {sig_naive:>14.2f} {sig_corr:>14.2f} {sig_hR:>14.2f}")
print(sep)

# Table 7: Background ABCD factorization check (bkg-only)
print(f"\n{'TABLE 7: BACKGROUND FACTORIZATION — Does A_bkg = B_bkg*C_bkg/D_bkg?':^200}")
print(sep)
header = (f"{'pT bin':>12} | {'A_bkg':>12} {'BC/D_bkg':>12} {'ratio':>10} | "
          f"{'A_jet':>12} {'BC/D_jet':>12} {'ratio':>10}")
print(header)
print(sep)
for i in range(n_bins):
    pt_label = f"[{pt_edges[i]},{pt_edges[i+1]})"
    Ab, Bb, Cb, Db = bkg["A"][i], bkg["B"][i], bkg["C"][i], bkg["D"][i]
    Aj, Bj, Cj, Dj = jet["A"][i], jet["B"][i], jet["C"][i], jet["D"][i]

    BCD_bkg = (Bb * Cb) / Db if Db > 0 else 0
    ratio_bkg = Ab / BCD_bkg if BCD_bkg > 0 else float('nan')

    BCD_jet = (Bj * Cj) / Dj if Dj > 0 else 0
    ratio_jet = Aj / BCD_jet if BCD_jet > 0 else float('nan')

    print(f"{pt_label:>12} | {Ab:>12.2f} {BCD_bkg:>12.2f} {ratio_bkg:>10.4f} | "
          f"{Aj:>12.2f} {BCD_jet:>12.2f} {ratio_jet:>10.4f}")
print(sep)

# Summary: signal fraction in each region
print(f"\n{'SUMMARY: Signal fraction (purity) in each ABCD region':^200}")
print(sep)
header = f"{'pT bin':>12} | {'f_A (purity)':>12} {'f_B':>12} {'f_C':>12} {'f_D':>12} | {'1-f_A':>10}"
print(header)
print(sep)
for i in range(n_bins):
    pt_label = f"[{pt_edges[i]},{pt_edges[i+1]})"
    fA = signal["A"][i] / total["A"][i] if total["A"][i] > 0 else 0
    fB = signal["B"][i] / total["B"][i] if total["B"][i] > 0 else 0
    fC = signal["C"][i] / total["C"][i] if total["C"][i] > 0 else 0
    fD = signal["D"][i] / total["D"][i] if total["D"][i] > 0 else 0
    print(f"{pt_label:>12} | {fA:>12.4f} {fB:>12.4f} {fC:>12.4f} {fD:>12.4f} | {1-fA:>10.4f}")
print(sep)

print("\nDone.")
