#!/usr/bin/env python3
"""
Plot isoET distributions: tight vs non-tight BDT regions (jet MC),
and tight jet-only vs tight inclusive (signal contamination check).

Physics: The ABCD method assumes BDT and isolation are independent for
background. If the isoET shapes differ between tight and non-tight
background, R != 1 and the ABCD method has systematic bias.

Style: sPHENIX publication style (Helvetica, frameless legends, sPHENIX
header, minor ticks, no grid).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator, MaxNLocator

# ── sPHENIX matplotlib style ─────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'xtick.minor.width': 0.8,
    'ytick.minor.width': 0.8,
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'xtick.minor.size': 3,
    'ytick.minor.size': 3,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'axes.grid': False,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
})

import uproot

# ── Configuration ──────────────────────────────────────────────────
JET_FILE = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_bdt_nom.root"
INC_FILE = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root"
OUTDIR = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures"

PT_EDGES = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
NPTBINS = 12

ISOET_XMIN = -5.0
ISOET_XMAX = 20.0

# sPHENIX-style colors (approximate ROOT kRed+1, kAzure+2)
COL_TIGHT = '#E6342A'      # kRed+1
COL_NONTIGHT = '#3B7DD8'   # kAzure+2
COL_JET = '#E6342A'        # same red for jet-only
COL_INC = '#3B7DD8'        # same blue for inclusive
COL_RATIO = '#222222'


def add_sphenix_header(ax, x=0.04, y=0.96):
    """Add sPHENIX Simulation Internal header to axes."""
    ax.text(x, y, r'$\mathbf{\mathit{sPHENIX}}$ Simulation Internal',
            transform=ax.transAxes, fontsize=11, va='top', ha='left')
    ax.text(x, y - 0.08, r'PYTHIA8 Jet MC',
            transform=ax.transAxes, fontsize=9, va='top', ha='left',
            color='#444444')


def add_sphenix_header_inc(ax, x=0.04, y=0.96):
    """Add sPHENIX header for inclusive MC comparison."""
    ax.text(x, y, r'$\mathbf{\mathit{sPHENIX}}$ Simulation Internal',
            transform=ax.transAxes, fontsize=11, va='top', ha='left')
    ax.text(x, y - 0.08, r'PYTHIA8 Tight BDT',
            transform=ax.transAxes, fontsize=9, va='top', ha='left',
            color='#444444')


def read_hist(root_file, name):
    """Read a TH1D, return (values, edges, centers)."""
    h = root_file[name]
    vals = h.values()
    edges = h.axis().edges()
    centers = 0.5 * (edges[:-1] + edges[1:])
    return vals, edges, centers


def normalize_hist(vals):
    """Normalize to unit area. Returns (normalized, total)."""
    total = vals.sum()
    if total <= 0:
        return vals, 0.0
    return vals / total, total


def ks_test_weighted(vals1, edges, vals2):
    """KS test between two weighted histograms via CDF comparison."""
    if vals1.sum() <= 0 or vals2.sum() <= 0:
        return 0.0, 0.0
    cdf1 = np.cumsum(vals1) / vals1.sum()
    cdf2 = np.cumsum(vals2) / vals2.sum()
    ks_stat = np.max(np.abs(cdf1 - cdf2))
    n1, n2 = vals1.sum(), vals2.sum()
    n_eff = np.sqrt(n1 * n2 / (n1 + n2))
    lam = (n_eff + 0.12 + 0.11 / n_eff) * ks_stat
    if lam < 1e-10:
        return ks_stat, 1.0
    pval = sum(2.0 * ((-1) ** (k + 1)) * np.exp(-2.0 * k * k * lam * lam)
               for k in range(1, 101))
    return ks_stat, max(0.0, min(1.0, pval))


def weighted_mean(vals, centers):
    """Weighted mean of histogram."""
    total = vals.sum()
    if total <= 0:
        return np.nan
    return np.sum(vals * centers) / total


def make_overlay_figure(f_file1, f_file2, hist_prefix1, hist_prefix2,
                        label1, label2, col1, col2, ratio_label,
                        header_func, outpath, summary_header):
    """Generic 4x3 grid isoET overlay with ratio panels."""
    print("=" * 72)
    print(summary_header)
    print("=" * 72)

    fig = plt.figure(figsize=(18, 16))
    outer_gs = gridspec.GridSpec(4, 3, hspace=0.32, wspace=0.26,
                                 left=0.065, right=0.97, top=0.97, bottom=0.04)
    summary_rows = []

    for ipt in range(NPTBINS):
        pt_lo, pt_hi = PT_EDGES[ipt], PT_EDGES[ipt + 1]

        vals1, edges, centers = read_hist(f_file1, f"{hist_prefix1}{ipt}")
        vals2, _, _ = read_hist(f_file2, f"{hist_prefix2}{ipt}")

        if vals1.sum() <= 0 or vals2.sum() <= 0:
            print(f"  pT bin {ipt} [{pt_lo}-{pt_hi}] GeV: EMPTY, skipping")
            summary_rows.append((ipt, pt_lo, pt_hi, np.nan, np.nan,
                                 np.nan, np.nan, np.nan))
            continue

        norm1, _ = normalize_hist(vals1)
        norm2, _ = normalize_hist(vals2)
        mean1 = weighted_mean(vals1, centers)
        mean2 = weighted_mean(vals2, centers)
        ks_stat, ks_pval = ks_test_weighted(vals1, edges, vals2)
        summary_rows.append((ipt, pt_lo, pt_hi, ks_stat, ks_pval,
                             mean1, mean2, mean1 - mean2))

        row, col = ipt // 3, ipt % 3
        inner_gs = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer_gs[row, col],
            height_ratios=[3, 1], hspace=0.0)

        ax_main = fig.add_subplot(inner_gs[0])
        ax_ratio = fig.add_subplot(inner_gs[1], sharex=ax_main)

        mask = (centers >= ISOET_XMIN) & (centers <= ISOET_XMAX)

        # Main panel: overlaid distributions
        ax_main.step(centers[mask], norm1[mask], where='mid',
                     color=col1, linewidth=1.8, label=label1)
        ax_main.step(centers[mask], norm2[mask], where='mid',
                     color=col2, linewidth=1.8, linestyle='--', label=label2)

        ax_main.set_ylabel("Normalized", fontsize=10)
        ax_main.tick_params(labelbottom=False)
        ax_main.set_xlim(ISOET_XMIN, ISOET_XMAX)
        ax_main.xaxis.set_minor_locator(AutoMinorLocator())
        ax_main.yaxis.set_minor_locator(AutoMinorLocator())

        ymax = max(norm1[mask].max(), norm2[mask].max()) * 1.5
        ax_main.set_ylim(0, ymax)

        # pT label and KS test
        ax_main.text(0.96, 0.94,
                     f"${pt_lo} < p_T < {pt_hi}$ GeV",
                     transform=ax_main.transAxes, fontsize=10,
                     ha='right', va='top', fontweight='bold')
        ks_color = '#228B22' if ks_pval > 0.05 else '#CC0000'
        ax_main.text(0.96, 0.80,
                     f"KS stat = {ks_stat:.3f}",
                     transform=ax_main.transAxes, fontsize=8.5,
                     ha='right', va='top', color=ks_color)
        ax_main.text(0.96, 0.69,
                     r"$\langle E_T^{\mathrm{iso}} \rangle$: "
                     f"{mean1:.2f}, {mean2:.2f}",
                     transform=ax_main.transAxes, fontsize=8,
                     ha='right', va='top', color='#555555')

        if ipt == 0:
            header_func(ax_main)
            ax_main.legend(fontsize=9, loc='center left', frameon=False,
                           bbox_to_anchor=(0.0, 0.55))

        # Ratio panel
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(norm2 > 0, norm1 / norm2, np.nan)
        ax_ratio.step(centers[mask], ratio[mask], where='mid',
                      color=COL_RATIO, linewidth=1.2)
        ax_ratio.axhline(1.0, color='#888888', linestyle='--', linewidth=1.0)
        ax_ratio.set_ylabel(ratio_label, fontsize=9)
        ax_ratio.set_xlabel(r"$E_T^{\mathrm{iso}}$ [GeV]", fontsize=10)
        ax_ratio.set_ylim(0.0, 3.0)
        ax_ratio.xaxis.set_minor_locator(AutoMinorLocator())
        ax_ratio.yaxis.set_minor_locator(AutoMinorLocator())
        ax_ratio.yaxis.set_major_locator(MaxNLocator(3))

    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {outpath}")
    plt.close(fig)

    # Print summary table
    print(f"\n{'Bin':>3s}  {'pT range':>12s}  {'KS stat':>10s}  {'KS p-val':>12s}  "
          f"{'<isoET> 1':>10s}  {'<isoET> 2':>11s}  {'Delta':>8s}")
    print("-" * 90)
    for r in summary_rows:
        ipt, pl, ph, ks_s, ks_p, m1, m2, dm = r
        if np.isnan(ks_s):
            print(f"{ipt:3d}  [{pl:4.0f}-{ph:4.0f}] GeV  {'EMPTY':>10s}")
        else:
            print(f"{ipt:3d}  [{pl:4.0f}-{ph:4.0f}] GeV  {ks_s:10.4f}  {ks_p:12.4e}  "
                  f"{m1:10.3f}  {m2:11.3f}  {dm:+8.3f}")
    print("=" * 90)


# ── Main ──────────────────────────────────────────────────────────
f_jet = uproot.open(JET_FILE)
f_inc = uproot.open(INC_FILE)

# Figure 1: Tight vs Non-tight isoET (jet MC — pure background)
make_overlay_figure(
    f_file1=f_jet, f_file2=f_jet,
    hist_prefix1="h_tight_isoET_0_", hist_prefix2="h_nontight_isoET_0_",
    label1="Tight BDT (bkg)", label2="Non-tight BDT (bkg)",
    col1=COL_TIGHT, col2=COL_NONTIGHT,
    ratio_label="T / NT",
    header_func=add_sphenix_header,
    outpath=f"{OUTDIR}/isoET_tight_vs_nontight_jet_bdt_nom.pdf",
    summary_header="Figure 1: Tight vs Non-tight isoET (Jet MC, background only)"
)

# Figure 2: Tight jet-only vs Tight inclusive (signal contamination)
make_overlay_figure(
    f_file1=f_jet, f_file2=f_inc,
    hist_prefix1="h_tight_isoET_0_", hist_prefix2="h_tight_isoET_0_",
    label1="Jet MC (bkg only)", label2="Inclusive MC (sig+bkg)",
    col1=COL_JET, col2=COL_INC,
    ratio_label="Inc / Jet",
    header_func=add_sphenix_header_inc,
    outpath=f"{OUTDIR}/isoET_tight_jetonly_vs_inclusive_bdt_nom.pdf",
    summary_header="Figure 2: Tight isoET — Jet-only vs Inclusive MC"
)

print("\nDone.")
