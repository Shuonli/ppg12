#!/usr/bin/env python3
"""
Investigate ABCD leakage fractions and the MC purity correction ratio crossing 1.0.

Reads from:
  - Photon_final_bdt_nom_mc.root: h_leak_B/C/D, gpurity, gpurity_leak, g_purity_truth, g_mc_purity_fit_ratio
  - MC_efficiency_bdt_nom.root: raw signal counts in ABCD regions

Outputs (in plotting/figures/):
  - leakage_fractions_vs_pT.pdf
  - signal_counts_ABCD.pdf
  - purity_comparison.pdf
  - leakage_sensitivity.pdf
"""
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

FINAL_MC = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom_mc.root"
MC_EFF   = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root"
FIGDIR   = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures"

PT_BINS = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
PT_CENTERS = 0.5 * (PT_BINS[:-1] + PT_BINS[1:])
NBINS = 12

f_final = uproot.open(FINAL_MC)
f_eff   = uproot.open(MC_EFF)

def read_th1(f, name):
    h = f[name]
    return h.values(), h.errors()

def read_tgraph(f, name):
    g = f[name]
    x = np.array(g.member("fX"))
    y = np.array(g.member("fY"))
    try:
        ey = np.array(g.member("fEY"))
    except:
        try:
            eyh = np.array(g.member("fEYhigh"))
            eyl = np.array(g.member("fEYlow"))
            ey = 0.5 * (eyh + eyl)
        except:
            ey = np.zeros_like(y)
    return x, y, ey

# Read all data
leak_B_vals, leak_B_errs = read_th1(f_final, "h_leak_B")
leak_C_vals, leak_C_errs = read_th1(f_final, "h_leak_C")
leak_D_vals, leak_D_errs = read_th1(f_final, "h_leak_D")

sig_A, _ = read_th1(f_eff, "h_tight_iso_cluster_signal_0")
sig_B, _ = read_th1(f_eff, "h_tight_noniso_cluster_signal_0")
sig_C, _ = read_th1(f_eff, "h_nontight_iso_cluster_signal_0")
sig_D, _ = read_th1(f_eff, "h_nontight_noniso_cluster_signal_0")

tot_A, _ = read_th1(f_eff, "h_tight_iso_cluster_0")
tot_B, _ = read_th1(f_eff, "h_tight_noniso_cluster_0")
tot_C, _ = read_th1(f_eff, "h_nontight_iso_cluster_0")
tot_D, _ = read_th1(f_eff, "h_nontight_noniso_cluster_0")

bkg_A = tot_A - sig_A
bkg_B = tot_B - sig_B
bkg_C = tot_C - sig_C
bkg_D = tot_D - sig_D

pur_x, pur_y, pur_ey = read_tgraph(f_final, "gpurity")
pur_leak_x, pur_leak_y, pur_leak_ey = read_tgraph(f_final, "gpurity_leak")
truth_x, truth_y, truth_ey = read_tgraph(f_final, "g_purity_truth")
ratio_x, ratio_y, ratio_ey = read_tgraph(f_final, "g_mc_purity_fit_ratio")

# === Print tables ===
print("=" * 80)
print("LEAKAGE FRACTIONS vs pT  (cX = signal_in_X / signal_in_A)")
print("=" * 80)
for i in range(NBINS):
    print(f"  pT={PT_CENTERS[i]:5.1f}: cB={leak_B_vals[i]:.5f}  cC={leak_C_vals[i]:.5f}  cD={leak_D_vals[i]:.5f}")

print("\nSIGNAL CONTAMINATION IN B, C, D:")
for i in range(NBINS):
    sfB = 100*sig_B[i]/tot_B[i] if tot_B[i]>0 else 0
    sfC = 100*sig_C[i]/tot_C[i] if tot_C[i]>0 else 0
    sfD = 100*sig_D[i]/tot_D[i] if tot_D[i]>0 else 0
    print(f"  pT={PT_CENTERS[i]:5.1f}: B={sfB:.1f}%  C={sfC:.1f}%  D={sfD:.1f}%")

print("\nBACKGROUND OVERESTIMATE FACTOR (B*C/D / true_bkg_A):")
for i in range(NBINS):
    tot_bcd = tot_B[i]*tot_C[i]/tot_D[i] if tot_D[i]>0 else 0
    r = tot_bcd / bkg_A[i] if bkg_A[i]>0 else 0
    print(f"  pT={PT_CENTERS[i]:5.1f}: B*C/D={tot_bcd:.0f}  true_bkg={bkg_A[i]:.0f}  ratio={r:.1f}x")

# === Plots ===
plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "axes.titlesize": 15,
    "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 11, "figure.dpi": 150})

# Plot 1: Leakage fractions
fig, axes = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
ax = axes[0]
ax.errorbar(PT_CENTERS, leak_B_vals[:NBINS], yerr=leak_B_errs[:NBINS],
            fmt="o-", label=r"$c_B$ (tight, non-iso)", color="C0", capsize=3, ms=5)
ax.errorbar(PT_CENTERS, leak_C_vals[:NBINS], yerr=leak_C_errs[:NBINS],
            fmt="s-", label=r"$c_C$ (non-tight, iso)", color="C1", capsize=3, ms=5)
ax.errorbar(PT_CENTERS, leak_D_vals[:NBINS], yerr=leak_D_errs[:NBINS],
            fmt="^-", label=r"$c_D$ (non-tight, non-iso)", color="C2", capsize=3, ms=5)
ax.axvspan(13, 15, alpha=0.12, color="red", label="Ratio crosses 1.0")
ax.set_ylabel(r"Leakage fraction $c_X = N_X^{sig}/N_A^{sig}$")
ax.set_title("Signal leakage fractions from region A into B, C, D")
ax.legend(loc="upper right"); ax.set_ylim(bottom=-0.005, top=0.16); ax.grid(True, alpha=0.3)
ax2 = axes[1]
ax2.semilogy(PT_CENTERS, np.maximum(leak_B_vals[:NBINS], 1e-6), "o-", color="C0", ms=4)
ax2.semilogy(PT_CENTERS, np.maximum(leak_C_vals[:NBINS], 1e-6), "s-", color="C1", ms=4)
ax2.semilogy(PT_CENTERS, np.maximum(leak_D_vals[:NBINS], 1e-6), "^-", color="C2", ms=4)
ax2.axvspan(13, 15, alpha=0.12, color="red")
ax2.set_xlabel(r"$p_T$ [GeV]"); ax2.set_ylabel("Log scale")
ax2.set_ylim(1e-4, 0.5); ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FIGDIR}/leakage_fractions_vs_pT.pdf", bbox_inches="tight"); plt.close()

# Plot 2: Signal counts + contamination %
fig, axes = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
ax = axes[0]
ax.plot(PT_CENTERS, sig_A[:NBINS], "o-", label="A (tight+iso)", color="C3", ms=5)
ax.plot(PT_CENTERS, sig_B[:NBINS], "s-", label="B (tight+non-iso)", color="C0", ms=5)
ax.plot(PT_CENTERS, sig_C[:NBINS], "^-", label="C (non-tight+iso)", color="C1", ms=5)
ax.plot(PT_CENTERS, sig_D[:NBINS], "v-", label="D (non-tight+non-iso)", color="C2", ms=5)
ax.axvspan(13, 15, alpha=0.12, color="red")
ax.set_ylabel("Signal count"); ax.set_title("Truth-matched signal in ABCD regions")
ax.legend(); ax.set_yscale("log"); ax.grid(True, alpha=0.3)
ax2 = axes[1]
sfB = np.where(tot_B[:NBINS]>0, 100*sig_B[:NBINS]/tot_B[:NBINS], 0)
sfC = np.where(tot_C[:NBINS]>0, 100*sig_C[:NBINS]/tot_C[:NBINS], 0)
sfD = np.where(tot_D[:NBINS]>0, 100*sig_D[:NBINS]/tot_D[:NBINS], 0)
ax2.plot(PT_CENTERS, sfB, "s-", color="C0", label="Signal % in B", ms=4)
ax2.plot(PT_CENTERS, sfC, "^-", color="C1", label="Signal % in C", ms=4)
ax2.plot(PT_CENTERS, sfD, "v-", color="C2", label="Signal % in D", ms=4)
ax2.axvspan(13, 15, alpha=0.12, color="red")
ax2.set_xlabel(r"$p_T$ [GeV]"); ax2.set_ylabel("Signal fraction [%]")
ax2.legend(fontsize=9); ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FIGDIR}/signal_counts_ABCD.pdf", bbox_inches="tight"); plt.close()

# Plot 3: Purity comparison + ratio
fig, axes = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
ax = axes[0]
n = min(len(pur_x), len(truth_x))
ax.errorbar(truth_x[:n], truth_y[:n], yerr=truth_ey[:n],
            fmt="o-", label="Truth purity", color="black", capsize=3, zorder=5, ms=5)
ax.errorbar(pur_x[:n], pur_y[:n], yerr=pur_ey[:n],
            fmt="s--", label="ABCD (no leakage)", color="C0", capsize=3, ms=5)
ax.errorbar(pur_leak_x[:n], pur_leak_y[:n], yerr=pur_leak_ey[:n],
            fmt="^--", label="ABCD (with leakage)", color="C1", capsize=3, ms=5)
ax.axvspan(13, 15, alpha=0.12, color="red", label="Crossover region")
ax.axhline(1.0, color="gray", ls=":", alpha=0.5)
ax.set_ylabel("Purity"); ax.set_title("MC closure: Truth vs ABCD purity")
ax.legend(loc="lower right", fontsize=10); ax.set_ylim(-0.2, 1.4); ax.grid(True, alpha=0.3)
ax2 = axes[1]
ax2.errorbar(ratio_x[:n], ratio_y[:n], yerr=ratio_ey[:n],
             fmt="o-", color="C3", capsize=3, ms=5, label="Truth / ABCD(leak fit)")
ax2.axhline(1.0, color="gray", ls="--")
ax2.axvspan(13, 15, alpha=0.12, color="red")
ax2.set_xlabel(r"$p_T$ [GeV]"); ax2.set_ylabel("Correction ratio")
ax2.set_ylim(0.5, 2.0); ax2.legend(); ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FIGDIR}/purity_comparison.pdf", bbox_inches="tight"); plt.close()

# Plot 4: Root cause
fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
tot_bcd = np.where(tot_D[:NBINS]>0, tot_B[:NBINS]*tot_C[:NBINS]/tot_D[:NBINS], 0)
ax = axes[0]
ax.plot(PT_CENTERS, tot_bcd, "s-", color="C0", label=r"$B \times C / D$ (ABCD bkg estimate)", ms=5)
ax.plot(PT_CENTERS, bkg_A[:NBINS], "o-", color="C3", label="True background in A", ms=5)
ax.plot(PT_CENTERS, tot_A[:NBINS], "^-", color="black", label="Total in A", ms=5)
ax.axvspan(13, 15, alpha=0.12, color="red")
ax.set_ylabel("Counts"); ax.set_yscale("log")
ax.set_title(r"Root cause: $B\times C/D$ vs true background in A")
ax.legend(fontsize=10); ax.grid(True, alpha=0.3)
ax = axes[1]
ratio_bcd = np.where(bkg_A[:NBINS]>0, tot_bcd/bkg_A[:NBINS], 0)
ax.plot(PT_CENTERS, ratio_bcd, "o-", color="C1", ms=5)
ax.axhline(1.0, color="gray", ls="--", alpha=0.7, label="Perfect ABCD")
ax.axvspan(13, 15, alpha=0.12, color="red")
ax.set_ylabel(r"$B\times C/D$ / true bkg")
ax.set_title("Overestimate factor"); ax.legend(); ax.grid(True, alpha=0.3)
ax = axes[2]
R_bkg = np.where((bkg_B[:NBINS]*bkg_C[:NBINS])>0,
    (bkg_A[:NBINS]*bkg_D[:NBINS])/(bkg_B[:NBINS]*bkg_C[:NBINS]), 0)
ax.plot(PT_CENTERS, R_bkg, "o-", color="C4", ms=6)
ax.axhline(1.0, color="gray", ls="--", alpha=0.7, label="R = 1 (ABCD assumption)")
ax.axvspan(13, 15, alpha=0.12, color="red")
ax.set_xlabel(r"$p_T$ [GeV]")
ax.set_ylabel(r"$R_{bkg} = A_{bkg} D_{bkg} / (B_{bkg} C_{bkg})$")
ax.set_title("Background correlation factor R"); ax.legend(); ax.set_ylim(0, 2.0); ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FIGDIR}/leakage_sensitivity.pdf", bbox_inches="tight"); plt.close()

print("All plots saved.")
