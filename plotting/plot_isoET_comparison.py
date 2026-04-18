#!/usr/bin/env python3
"""
plot_isoET_comparison.py
Compare three isolation ET definitions: 60 MeV tower threshold, 120 MeV tower
threshold, and topo-cluster R=0.4.

Produces:
  1. Mean isoET vs run number (overlay)
  2. Region A (tight+iso) yield vs run number (overlay)
  3. Region C (nontight+iso) yield vs run number (overlay)
  4. Ratio C/A vs run number (overlay)
  5. Stability summary (RMS/mean for mean isoET and yields)
  6. ROC AUC comparison bar chart per pT bin (from merged ROC file)

Usage:
  python3 plot_isoET_comparison.py [--outdir figures/isoET_comparison]
"""

import argparse
import os
import sys
import numpy as np

import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# sPHENIX-like style
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 10,
    "figure.figsize": (10, 5),
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
})

# Iso definitions: label, ROOT file, color, marker
ISO_DEFS = {
    "60 MeV":   {"file": "efficiencytool/rbrQA_60MeV.root",  "color": "#d62728", "ls": "-"},
    "120 MeV":  {"file": "efficiencytool/rbrQA_120MeV.root", "color": "#2ca02c", "ls": "--"},
    "topo R=0.4": {"file": "efficiencytool/rbrQA_topo.root", "color": "#1f77b4", "ls": "-."},
}

# ROC merged file
ROC_FILE = "efficiencytool/results/roc_isoET_merged_isoroc.root"

# Labels to focus on (matching IsoROC_calculator naming)
ROC_LABELS = {
    "60 MeV":     "iso60",
    "120 MeV":    "iso120",
    "topo R=0.4": "isotopo04",
}
ROC_ALL_LABELS = {
    "70 MeV": "iso70",
    "60 MeV": "iso60",
    "120 MeV": "iso120",
    "topo R=0.3": "isotopo03",
    "topo R=0.4": "isotopo04",
    "topo soft R=0.3": "isotoposoft03",
    "topo soft R=0.4": "isotoposoft04",
}
ROC_COLORS = {
    "70 MeV": "#ff7f0e",
    "60 MeV": "#d62728",
    "120 MeV": "#2ca02c",
    "topo R=0.3": "#9467bd",
    "topo R=0.4": "#1f77b4",
    "topo soft R=0.3": "#8c564b",
    "topo soft R=0.4": "#17becf",
}

SPHENIX_LABEL = r"$\bf{sPHENIX}$ Internal"
PP_LABEL = r"$p$+$p$ $\sqrt{s}$ = 200 GeV"
ETA_LABEL = r"|$\eta^{\gamma}$| < 0.7"

# Crossing angle boundary
RUN_XING_BOUNDARY = 51274


def read_graph(f, name):
    """Read TGraphErrors from uproot file, return (x, y, ey)."""
    g = f[name]
    x = np.array(g.member("fX"))
    y = np.array(g.member("fY"))
    ey = np.array(g.member("fEY")) if hasattr(g, "member") else np.zeros_like(y)
    try:
        ey = np.array(g.member("fEY"))
    except Exception:
        ey = np.zeros_like(y)
    return x, y, ey


def compute_roc_auc(h_sig, h_bg):
    """Compute ROC AUC from signal and background histograms."""
    sig_vals = h_sig.values()
    bg_vals = h_bg.values()
    total_sig = sig_vals.sum()
    total_bg = bg_vals.sum()
    if total_sig <= 0 or total_bg <= 0:
        return 0.0, np.array([]), np.array([])
    cum_sig = np.cumsum(sig_vals) / total_sig
    cum_bg = np.cumsum(bg_vals) / total_bg
    auc = np.trapezoid(cum_sig, cum_bg)
    return auc, cum_bg, cum_sig


def fit_pol1_split(x, y, ey, boundary=RUN_XING_BOUNDARY):
    """Weighted pol1 fits for runs < boundary (0 mrad) and >= boundary (1.5 mrad).

    Returns a 2-tuple (lo, hi) of dicts with keys slope, intercept, x (line
    endpoints), y (line endpoints), n (npts used). Side with <2 points is None.
    """
    out = []
    finite = np.isfinite(y) & np.isfinite(ey) & (y > 0)
    for sel in (x < boundary, x >= boundary):
        m = finite & sel
        if m.sum() < 2:
            out.append(None)
            continue
        xi, yi, eyi = x[m], y[m], ey[m]
        w = np.where(eyi > 0, 1.0 / eyi, 0.0)
        if np.count_nonzero(w) < 2:
            w = np.ones_like(yi)
        try:
            slope, intercept = np.polyfit(xi, yi, 1, w=w)
        except Exception:
            slope, intercept = np.polyfit(xi, yi, 1)
        x_line = np.array([xi.min(), xi.max()])
        y_line = slope * x_line + intercept
        out.append({"slope": float(slope), "intercept": float(intercept),
                    "x": x_line, "y": y_line, "n": int(m.sum())})
    return out


def _slope_str(res):
    return "n/a" if res is None else f"{res['slope']:+.2e}"


def overlay_pol1_fits(ax, x, y, ey, color, boundary=RUN_XING_BOUNDARY):
    """Draw the two regime fits as solid lines, return the fit-result list."""
    fits = fit_pol1_split(x, y, ey, boundary)
    for res in fits:
        if res is None:
            continue
        ax.plot(res["x"], res["y"], color=color, linestyle="-",
                linewidth=1.4, alpha=0.95, zorder=3)
    return fits


def plot_run_overlay(datasets, graph_name, ylabel, title, ax, show_xing=True):
    """Overlay a graph from multiple rbrQA files on given axes with pol1 fits
    computed separately below/above the 0/1.5 mrad crossing-angle boundary."""
    for iso_name, info in datasets.items():
        x, y, ey = info["data"][graph_name]
        # Sort by run number
        idx = np.argsort(x)
        x, y, ey = x[idx], y[idx], ey[idx]
        # Filter out zeros (missing runs)
        mask = y > 0
        xm, ym, eym = x[mask], y[mask], ey[mask]
        fits = overlay_pol1_fits(ax, xm, ym, eym, info["color"])
        label = (f"{iso_name}  b(0mrad)={_slope_str(fits[0])}, "
                 f"b(1.5mrad)={_slope_str(fits[1])}")
        ax.errorbar(xm, ym, yerr=eym,
                     fmt=".", markersize=2, linewidth=0.5, capsize=0,
                     color=info["color"], linestyle="none",
                     label=label, alpha=0.7, zorder=2)
    if show_xing:
        ax.axvline(RUN_XING_BOUNDARY, color="gray", linestyle=":", linewidth=1, alpha=0.6)
        ymin, ymax = ax.get_ylim()
        ax.text(RUN_XING_BOUNDARY + 50, ymax * 0.95, "0/1.5 mrad",
                fontsize=8, color="gray", va="top")
    ax.set_xlabel("Run Number")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="best", framealpha=0.8, fontsize=9)


def compute_stability(x, y, ey):
    """Compute stability metrics for per-run values."""
    mask = y > 0
    y_good = y[mask]
    if len(y_good) < 2:
        return {"mean": 0, "std": 0, "rms_frac": 0, "nruns": 0}
    mean = np.mean(y_good)
    std = np.std(y_good)
    return {
        "mean": mean,
        "std": std,
        "rms_frac": std / mean if mean > 0 else 0,
        "nruns": len(y_good),
    }


def main():
    parser = argparse.ArgumentParser(description="Compare isolation ET definitions")
    parser.add_argument("--outdir", default="figures/isoET_comparison")
    parser.add_argument("--basedir", default="/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12")
    args = parser.parse_args()

    basedir = args.basedir
    outdir = os.path.join(basedir, args.outdir)
    os.makedirs(outdir, exist_ok=True)

    # --- Load rbrQA data ---
    datasets = {}
    graphs_to_read = [
        "gr_avg_isoET", "gr_tight_iso", "gr_tight_noniso",
        "gr_nontight_iso", "gr_nontight_noniso", "gr_common",
        "gr_tight", "gr_nontight",
    ]

    for iso_name, info in ISO_DEFS.items():
        fpath = os.path.join(basedir, info["file"])
        if not os.path.exists(fpath):
            print(f"WARNING: {fpath} not found, skipping {iso_name}")
            continue
        f = uproot.open(fpath)
        data = {}
        for gname in graphs_to_read:
            try:
                data[gname] = read_graph(f, gname)
            except Exception as e:
                print(f"  WARNING: {gname} not found in {fpath}: {e}")
        datasets[iso_name] = {"color": info["color"], "ls": info["ls"], "data": data}

    if not datasets:
        print("ERROR: No rbrQA files found")
        sys.exit(1)

    print(f"Loaded {len(datasets)} isolation definitions: {list(datasets.keys())}")

    # --- Multi-page PDF report ---
    pdf_path = os.path.join(outdir, "isoET_comparison.pdf")
    with PdfPages(pdf_path) as pdf:

        # Page 1: Mean isoET vs run
        fig, ax = plt.subplots(figsize=(12, 5))
        plot_run_overlay(datasets, "gr_avg_isoET",
                         r"$\langle E_T^{\mathrm{iso}} \rangle$ [GeV]",
                         "Mean Isolation ET vs Run Number", ax)
        fig.text(0.15, 0.92, SPHENIX_LABEL, fontsize=12, transform=fig.transFigure)
        fig.text(0.15, 0.88, f"{PP_LABEL}, {ETA_LABEL}", fontsize=10, transform=fig.transFigure)
        pdf.savefig(fig)
        fig.savefig(os.path.join(outdir, "mean_isoET_vs_run.pdf"))
        plt.close(fig)

        # Page 2: Region A (tight+iso) yield vs run
        fig, ax = plt.subplots(figsize=(12, 5))
        plot_run_overlay(datasets, "gr_tight_iso",
                         "Clusters / Lumi [mb$^{-1}$]",
                         "Region A (tight + iso) Yield vs Run Number", ax)
        pdf.savefig(fig)
        fig.savefig(os.path.join(outdir, "tight_iso_yield_vs_run.pdf"))
        plt.close(fig)

        # Page 3: Region C (nontight+iso) yield vs run
        fig, ax = plt.subplots(figsize=(12, 5))
        plot_run_overlay(datasets, "gr_nontight_iso",
                         "Clusters / Lumi [mb$^{-1}$]",
                         "Region C (nontight + iso) Yield vs Run Number", ax)
        pdf.savefig(fig)
        fig.savefig(os.path.join(outdir, "nontight_iso_yield_vs_run.pdf"))
        plt.close(fig)

        # Page 4: Ratio C/A vs run
        fig, ax = plt.subplots(figsize=(12, 5))
        for iso_name, info in datasets.items():
            if "gr_tight_iso" not in info["data"] or "gr_nontight_iso" not in info["data"]:
                continue
            xA, yA, eyA = info["data"]["gr_tight_iso"]
            xC, yC, eyC = info["data"]["gr_nontight_iso"]
            # Match runs (both sorted from same map iteration)
            idx = np.argsort(xA)
            xA, yA, eyA = xA[idx], yA[idx], eyA[idx]
            idx = np.argsort(xC)
            xC, yC, eyC = xC[idx], yC[idx], eyC[idx]
            # Find common runs
            common_runs = np.intersect1d(xA, xC)
            maskA = np.isin(xA, common_runs)
            maskC = np.isin(xC, common_runs)
            yA_c = yA[maskA]
            yC_c = yC[maskC]
            eyA_c = eyA[maskA]
            eyC_c = eyC[maskC]
            x_c = xA[maskA]
            # Compute ratio where A > 0
            valid = yA_c > 0
            ratio = np.where(valid, yC_c / yA_c, 0)
            ratio_err = np.where(valid,
                ratio * np.sqrt((eyA_c / np.where(yA_c > 0, yA_c, 1))**2 +
                                (eyC_c / np.where(yC_c > 0, yC_c, 1))**2), 0)
            x_v, r_v, re_v = x_c[valid], ratio[valid], ratio_err[valid]
            fits = overlay_pol1_fits(ax, x_v, r_v, re_v, info["color"])
            label = (f"{iso_name}  b(0mrad)={_slope_str(fits[0])}, "
                     f"b(1.5mrad)={_slope_str(fits[1])}")
            ax.errorbar(x_v, r_v, yerr=re_v,
                        fmt=".", markersize=2, linewidth=0.5, capsize=0,
                        color=info["color"], label=label, alpha=0.7, zorder=2)
        ax.axvline(RUN_XING_BOUNDARY, color="gray", linestyle=":", linewidth=1, alpha=0.6)
        ax.set_xlabel("Run Number")
        ax.set_ylabel("C / A (nontight iso / tight iso)")
        ax.set_title("ABCD Ratio C/A vs Run Number")
        ax.legend(loc="best", framealpha=0.8, fontsize=9)
        pdf.savefig(fig)
        fig.savefig(os.path.join(outdir, "ratio_CA_vs_run.pdf"))
        plt.close(fig)

        # Page 5: Stability summary table
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.axis("off")
        rows = []
        headers = ["Iso Definition", "Mean <isoET>", "StdDev", "RMS/Mean",
                    "Mean Yield A", "RMS/Mean A", "Mean Yield C", "RMS/Mean C"]
        for iso_name, info in datasets.items():
            s_iso = compute_stability(*info["data"].get("gr_avg_isoET", (np.array([]), np.array([]), np.array([]))))
            s_A = compute_stability(*info["data"].get("gr_tight_iso", (np.array([]), np.array([]), np.array([]))))
            s_C = compute_stability(*info["data"].get("gr_nontight_iso", (np.array([]), np.array([]), np.array([]))))
            rows.append([
                iso_name,
                f"{s_iso['mean']:.3f}",
                f"{s_iso['std']:.3f}",
                f"{s_iso['rms_frac']:.3f}",
                f"{s_A['mean']:.0f}",
                f"{s_A['rms_frac']:.3f}",
                f"{s_C['mean']:.0f}",
                f"{s_C['rms_frac']:.3f}",
            ])
        table = ax.table(cellText=rows, colLabels=headers, loc="center",
                         cellLoc="center", colColours=["#ddeeff"] * len(headers))
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.0, 1.8)
        ax.set_title("Run-by-Run Stability Summary", fontsize=14, pad=20)
        pdf.savefig(fig)
        fig.savefig(os.path.join(outdir, "stability_summary.pdf"))
        plt.close(fig)

        # --- ROC AUC comparison ---
        roc_path = os.path.join(basedir, ROC_FILE)
        if os.path.exists(roc_path):
            froc = uproot.open(roc_path)
            h_pt = froc["h_pT_bins"]
            edges = h_pt.axis().edges()
            n_pt = len(edges) - 1

            # Page 6: ROC AUC bar chart (focused 3 definitions)
            fig, ax = plt.subplots(figsize=(10, 6))
            x_pos = np.arange(n_pt)
            width = 0.25
            for i, (iso_name, roc_label) in enumerate(ROC_LABELS.items()):
                aucs = []
                for ipt in range(n_pt):
                    h_sig = froc.get(f"h_signal_{roc_label}_pt{ipt}", None)
                    h_bg = froc.get(f"h_bg_{roc_label}_pt{ipt}", None)
                    if h_sig and h_bg:
                        auc, _, _ = compute_roc_auc(h_sig, h_bg)
                        aucs.append(auc)
                    else:
                        aucs.append(0)
                bars = ax.bar(x_pos + i * width, aucs, width,
                              label=iso_name, color=ISO_DEFS.get(iso_name, {}).get("color", f"C{i}"),
                              edgecolor="black", linewidth=0.5)
                for bar, auc in zip(bars, aucs):
                    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.002,
                            f"{auc:.3f}", ha="center", va="bottom", fontsize=8)
            ax.set_xticks(x_pos + width)
            ax.set_xticklabels([f"{edges[i]:.0f}-{edges[i+1]:.0f}" for i in range(n_pt)])
            ax.set_xlabel(r"$E_T$ bin [GeV]")
            ax.set_ylabel("ROC AUC")
            ax.set_title("Signal/Background Discrimination: ROC AUC by Isolation Definition")
            ax.set_ylim(0.75, 0.92)
            ax.legend()
            ax.text(0.02, 0.98, SPHENIX_LABEL, fontsize=11, transform=ax.transAxes, va="top")
            ax.text(0.02, 0.93, PP_LABEL, fontsize=9, transform=ax.transAxes, va="top")
            pdf.savefig(fig)
            fig.savefig(os.path.join(outdir, "roc_auc_comparison.pdf"))
            plt.close(fig)

            # Page 7: ROC AUC bar chart (all 7 definitions)
            fig, ax = plt.subplots(figsize=(14, 6))
            n_iso = len(ROC_ALL_LABELS)
            width = 0.8 / n_iso
            for i, (iso_name, roc_label) in enumerate(ROC_ALL_LABELS.items()):
                aucs = []
                for ipt in range(n_pt):
                    h_sig = froc.get(f"h_signal_{roc_label}_pt{ipt}", None)
                    h_bg = froc.get(f"h_bg_{roc_label}_pt{ipt}", None)
                    if h_sig and h_bg:
                        auc, _, _ = compute_roc_auc(h_sig, h_bg)
                        aucs.append(auc)
                    else:
                        aucs.append(0)
                ax.bar(x_pos + i * width, aucs, width,
                       label=iso_name, color=ROC_COLORS.get(iso_name, f"C{i}"),
                       edgecolor="black", linewidth=0.3)
            ax.set_xticks(x_pos + 0.4)
            ax.set_xticklabels([f"{edges[i]:.0f}-{edges[i+1]:.0f}" for i in range(n_pt)])
            ax.set_xlabel(r"$E_T$ bin [GeV]")
            ax.set_ylabel("ROC AUC")
            ax.set_title("Signal/Background Discrimination: All 7 Isolation Definitions")
            ax.set_ylim(0.75, 0.92)
            ax.legend(ncol=2, fontsize=9)
            pdf.savefig(fig)
            fig.savefig(os.path.join(outdir, "roc_auc_all_definitions.pdf"))
            plt.close(fig)

            # Page 8: ROC curves for one representative pT bin (15-20 GeV)
            ipt_rep = 1  # 15-20 GeV
            fig, ax = plt.subplots(figsize=(7, 7))
            for iso_name, roc_label in ROC_ALL_LABELS.items():
                h_sig = froc.get(f"h_signal_{roc_label}_pt{ipt_rep}", None)
                h_bg = froc.get(f"h_bg_{roc_label}_pt{ipt_rep}", None)
                if h_sig and h_bg:
                    auc, bg_eff, sig_eff = compute_roc_auc(h_sig, h_bg)
                    ax.plot(bg_eff, sig_eff, label=f"{iso_name} (AUC={auc:.3f})",
                            color=ROC_COLORS.get(iso_name, "gray"), linewidth=1.5)
            ax.plot([0, 1], [0, 1], "k--", linewidth=0.5, alpha=0.5)
            ax.set_xlabel("Background Efficiency")
            ax.set_ylabel("Signal Efficiency")
            ax.set_title(f"ROC Curves: {edges[ipt_rep]:.0f}-{edges[ipt_rep+1]:.0f} GeV")
            ax.legend(loc="lower right", fontsize=9)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1.05)
            ax.set_aspect("equal")
            ax.grid(True, alpha=0.3)
            ax.text(0.05, 0.95, SPHENIX_LABEL, fontsize=11, transform=ax.transAxes, va="top")
            ax.text(0.05, 0.90, PP_LABEL, fontsize=9, transform=ax.transAxes, va="top")
            pdf.savefig(fig)
            fig.savefig(os.path.join(outdir, "roc_curves_15_20GeV.pdf"))
            plt.close(fig)

            # Page 9: isoET distribution comparison (signal and bg) for one pT bin
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            for ax_i, (cat, cat_title) in enumerate([("signal", "Signal"), ("bg", "Background")]):
                ax = axes[ax_i]
                for iso_name, roc_label in list(ROC_LABELS.items()):
                    h = froc.get(f"h_{cat}_{roc_label}_pt{ipt_rep}", None)
                    if h:
                        vals = h.values()
                        edges_h = h.axis().edges()
                        centers = 0.5 * (edges_h[:-1] + edges_h[1:])
                        # Normalize
                        integral = vals.sum()
                        if integral > 0:
                            vals = vals / integral
                        ax.step(centers, vals, where="mid",
                                color=ISO_DEFS.get(iso_name, {}).get("color", "gray"),
                                label=iso_name, linewidth=1.5)
                ax.set_xlabel(r"$E_T^{\mathrm{iso}}$ [GeV]")
                ax.set_ylabel("Normalized")
                ax.set_title(f"{cat_title} isoET Distribution ({edges[ipt_rep]:.0f}-{edges[ipt_rep+1]:.0f} GeV)")
                ax.set_xlim(-3, 10)
                ax.legend()
                ax.set_yscale("log")
            fig.tight_layout()
            pdf.savefig(fig)
            fig.savefig(os.path.join(outdir, "isoET_distributions_15_20GeV.pdf"))
            plt.close(fig)

        else:
            print(f"WARNING: ROC file not found at {roc_path}")

    print(f"\nAll plots saved to: {outdir}/")
    print(f"Combined PDF: {pdf_path}")

    # --- Print summary table to stdout ---
    print("\n" + "=" * 80)
    print("ISOLATION ET COMPARISON SUMMARY")
    print("=" * 80)
    print(f"\n{'Definition':<16} | {'<isoET>':<10} | {'StdDev':<10} | {'RMS/Mean':<10} | "
          f"{'Yield A':<10} | {'RMS/Mean A':<10} | {'Yield C':<10}")
    print("-" * 100)
    for iso_name, info in datasets.items():
        s_iso = compute_stability(*info["data"].get("gr_avg_isoET", (np.array([]), np.array([]), np.array([]))))
        s_A = compute_stability(*info["data"].get("gr_tight_iso", (np.array([]), np.array([]), np.array([]))))
        s_C = compute_stability(*info["data"].get("gr_nontight_iso", (np.array([]), np.array([]), np.array([]))))
        print(f"{iso_name:<16} | {s_iso['mean']:<10.3f} | {s_iso['std']:<10.3f} | "
              f"{s_iso['rms_frac']:<10.3f} | {s_A['mean']:<10.0f} | {s_A['rms_frac']:<10.3f} | "
              f"{s_C['mean']:<10.0f}")


if __name__ == "__main__":
    main()
