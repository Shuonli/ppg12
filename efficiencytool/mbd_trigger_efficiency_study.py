#!/usr/bin/env python3
"""
MBD trigger efficiency vs truth vertex z — MC study for PPG12.

Uses photon MC samples (photon5/10/20) to emulate the MBD minimum-bias
trigger offline and measure its efficiency as a function of the generated
(truth) vertex z position. The MBD trigger is emulated as a coincidence
of hits on both the north and south MBD arms (N&S >= threshold).

Output:
  - efficiencytool/results/mbd_trigger_eff_study.root  (histograms)
  - efficiencytool/results/mbd_trigger_eff_*.pdf       (plots)
"""

import os
import sys
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# ── Paths ──────────────────────────────────────────────────────────────
BASE = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12"
SAMPLES = {
    "photon5":  {"dir": f"{BASE}/FunWithxgboost/photon5",  "pt_range": (0, 14)},
    "photon10": {"dir": f"{BASE}/FunWithxgboost/photon10", "pt_range": (14, 30)},
    "photon20": {"dir": f"{BASE}/FunWithxgboost/photon20", "pt_range": (30, 100)},
}
OUTDIR = f"{BASE}/efficiencytool/results"
os.makedirs(OUTDIR, exist_ok=True)

# ── Branches to read ──────────────────────────────────────────────────
BRANCHES = [
    "vertexz_truth", "vertexz",
    "mbdnorthhit", "mbdsouthhit",
    "mbdnorthqsum", "mbdsouthqsum",
    "ncluster_CLUSTERINFO_CEMC",
]

# ── Binning ───────────────────────────────────────────────────────────
VTX_BINS_FULL  = np.linspace(-300, 300, 121)   # 5 cm bins, full range
VTX_BINS_ZOOM  = np.linspace(-60, 60, 49)      # 2.5 cm bins, analysis region
VTX_BINS_WIDE  = np.linspace(-150, 150, 61)    # 5 cm bins, intermediate

# MBD trigger emulation thresholds: require hits on BOTH north AND south
THRESHOLDS = [
    (1, 1, r"N$\geq$1 & S$\geq$1"),
    (2, 2, r"N$\geq$2 & S$\geq$2"),
    (1, 0, r"N$\geq$1 (north only)"),
    (0, 1, r"S$\geq$1 (south only)"),
]

# ── Style ─────────────────────────────────────────────────────────────
COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd"]
plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "figure.dpi": 150,
    "savefig.bbox": "tight",
    "savefig.dpi": 150,
})


def load_sample(sample_name, max_files=1):
    """Load one file per sample (10M events is plenty)."""
    info = SAMPLES[sample_name]
    files = sorted(
        f for f in os.listdir(info["dir"])
        if f.endswith(".root") and "bdt" in f
    )[:max_files]

    arrays = []
    for fname in files:
        path = os.path.join(info["dir"], fname)
        print(f"  Reading {path} ...")
        with uproot.open(path) as f:
            tree = f["slimtree"]
            arr = tree.arrays(BRANCHES, library="np")
            arrays.append(arr)
            print(f"    {len(arr['vertexz_truth'])} events")

    # Concatenate
    combined = {}
    for key in BRANCHES:
        combined[key] = np.concatenate([a[key] for a in arrays])
    return combined


def compute_efficiency(vtxz, passed, bins):
    """Compute efficiency and binomial errors in bins of vtxz."""
    h_total, _ = np.histogram(vtxz, bins=bins)
    h_pass, _  = np.histogram(vtxz[passed], bins=bins)

    with np.errstate(divide="ignore", invalid="ignore"):
        eff = np.where(h_total > 0, h_pass / h_total, np.nan)
        # Clopper-Pearson-like binomial error (normal approx for large N)
        err = np.where(
            h_total > 0,
            np.sqrt(eff * (1 - eff) / h_total),
            np.nan,
        )
    centers = 0.5 * (bins[:-1] + bins[1:])
    return centers, eff, err, h_total, h_pass


def plot_efficiency_vs_vtxz(results, bins, title_suffix, fname_suffix):
    """
    Plot trigger efficiency vs truth vtxz for all samples combined,
    with multiple trigger emulation thresholds.
    """
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 8), height_ratios=[3, 1],
        sharex=True, gridspec_kw={"hspace": 0.05}
    )

    # Upper panel: efficiency
    for i, (n_th, s_th, label) in enumerate(THRESHOLDS):
        vtxz_all = results["vertexz_truth"]
        north = results["mbdnorthhit"]
        south = results["mbdsouthhit"]

        passed = (north >= n_th) & (south >= s_th)
        centers, eff, err, h_tot, h_pass = compute_efficiency(
            vtxz_all, passed, bins
        )

        ax1.errorbar(
            centers, eff, yerr=err,
            fmt="o-", ms=3, lw=1.2, capsize=2,
            color=COLORS[i], label=label,
        )

    ax1.set_ylabel("Trigger efficiency")
    ax1.set_ylim(0, 1.05)
    ax1.legend(loc="lower center", ncol=2)
    ax1.set_title(f"MBD trigger efficiency vs truth vertex z {title_suffix}")
    ax1.axhline(1.0, color="gray", ls="--", lw=0.5)
    ax1.grid(axis="y", alpha=0.3)

    # Lower panel: event count
    h_tot_ref, _ = np.histogram(results["vertexz_truth"], bins=bins)
    centers_ref = 0.5 * (bins[:-1] + bins[1:])
    ax2.bar(
        centers_ref, h_tot_ref,
        width=np.diff(bins), color="gray", alpha=0.5, label="All events"
    )
    ax2.set_ylabel("Events / bin")
    ax2.set_xlabel("Truth vertex z [cm]")
    ax2.set_yscale("log")
    ax2.legend(fontsize=10)

    fig.align_ylabels()
    outpath = os.path.join(OUTDIR, f"mbd_trigger_eff_{fname_suffix}.pdf")
    fig.savefig(outpath)
    plt.close(fig)
    print(f"  Saved: {outpath}")
    return outpath


def plot_hit_asymmetry(results, fname_suffix):
    """2D histogram: north hits vs south hits, and asymmetry vs vtxz."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    vtxz = results["vertexz_truth"]
    north = results["mbdnorthhit"]
    south = results["mbdsouthhit"]

    # Panel 1: 2D north vs south hits (all vtxz)
    ax = axes[0]
    h, xedges, yedges, im = ax.hist2d(
        north, south,
        bins=[np.arange(-0.5, 25.5, 1), np.arange(-0.5, 25.5, 1)],
        cmin=1, norm=matplotlib.colors.LogNorm(),
    )
    fig.colorbar(im, ax=ax, label="Events")
    ax.set_xlabel("MBD North hits")
    ax.set_ylabel("MBD South hits")
    ax.set_title("All events")
    ax.plot([0, 24], [0, 24], "r--", lw=0.8)

    # Panel 2: Mean north and south hits vs vtxz
    ax = axes[1]
    bins_vz = np.linspace(-200, 200, 81)
    centers_vz = 0.5 * (bins_vz[:-1] + bins_vz[1:])
    idx = np.digitize(vtxz, bins_vz) - 1
    valid = (idx >= 0) & (idx < len(centers_vz))

    mean_n = np.full(len(centers_vz), np.nan)
    mean_s = np.full(len(centers_vz), np.nan)
    for ib in range(len(centers_vz)):
        mask = (idx == ib) & valid
        if mask.sum() > 10:
            mean_n[ib] = north[mask].mean()
            mean_s[ib] = south[mask].mean()

    ax.plot(centers_vz, mean_n, "b-", lw=1.5, label="North")
    ax.plot(centers_vz, mean_s, "r-", lw=1.5, label="South")
    ax.set_xlabel("Truth vertex z [cm]")
    ax.set_ylabel("Mean MBD hits")
    ax.legend()
    ax.set_title("Mean hits vs vertex z")
    ax.axvline(0, color="gray", ls=":", lw=0.5)
    ax.grid(alpha=0.3)

    # Panel 3: Hit asymmetry (N-S)/(N+S) vs vtxz
    ax = axes[2]
    total = north.astype(float) + south.astype(float)
    asym = np.where(total > 0, (north - south) / total, np.nan)

    mean_asym = np.full(len(centers_vz), np.nan)
    for ib in range(len(centers_vz)):
        mask = (idx == ib) & valid & np.isfinite(asym)
        if mask.sum() > 10:
            mean_asym[ib] = asym[mask].mean()

    ax.plot(centers_vz, mean_asym, "k-", lw=1.5)
    ax.set_xlabel("Truth vertex z [cm]")
    ax.set_ylabel(r"$\langle(N_{\rm north} - N_{\rm south})/(N_{\rm north} + N_{\rm south})\rangle$")
    ax.set_title("Hit asymmetry vs vertex z")
    ax.axhline(0, color="gray", ls="--", lw=0.8)
    ax.axvline(0, color="gray", ls=":", lw=0.5)
    ax.set_ylim(-1, 1)
    ax.grid(alpha=0.3)

    fig.suptitle("MBD hit response vs truth vertex z (photon MC)", fontsize=15, y=1.02)
    fig.tight_layout()
    outpath = os.path.join(OUTDIR, f"mbd_hit_asymmetry_{fname_suffix}.pdf")
    fig.savefig(outpath)
    plt.close(fig)
    print(f"  Saved: {outpath}")
    return outpath


def plot_efficiency_vs_vtxz_per_sample(all_data, bins, fname_suffix):
    """Show efficiency for each photon sample separately (N&S>=1)."""
    fig, ax = plt.subplots(figsize=(10, 6))

    for i, (sample, data) in enumerate(all_data.items()):
        pt_lo, pt_hi = SAMPLES[sample]["pt_range"]
        vtxz = data["vertexz_truth"]
        passed = (data["mbdnorthhit"] >= 1) & (data["mbdsouthhit"] >= 1)
        centers, eff, err, _, _ = compute_efficiency(vtxz, passed, bins)

        ax.errorbar(
            centers, eff, yerr=err,
            fmt="o-", ms=3, lw=1.2, capsize=2,
            color=COLORS[i],
            label=f"{sample} ({pt_lo}–{pt_hi} GeV)",
        )

    ax.set_xlabel("Truth vertex z [cm]")
    ax.set_ylabel("Trigger efficiency (N&S$\\geq$1)")
    ax.set_ylim(0, 1.05)
    ax.legend()
    ax.set_title("MBD trigger efficiency by photon sample")
    ax.axhline(1.0, color="gray", ls="--", lw=0.5)
    ax.grid(axis="y", alpha=0.3)

    outpath = os.path.join(OUTDIR, f"mbd_trigger_eff_per_sample_{fname_suffix}.pdf")
    fig.savefig(outpath)
    plt.close(fig)
    print(f"  Saved: {outpath}")
    return outpath


def plot_efficiency_with_cluster_cut(results, bins, fname_suffix):
    """
    Compare trigger efficiency for:
    - All events (inclusive denominator)
    - Events with >= 1 EMCal cluster (analysis-relevant denominator)
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    vtxz = results["vertexz_truth"]
    north = results["mbdnorthhit"]
    south = results["mbdsouthhit"]
    ncl = results["ncluster_CLUSTERINFO_CEMC"]

    passed_trig = (north >= 1) & (south >= 1)

    # All events
    centers, eff_all, err_all, _, _ = compute_efficiency(vtxz, passed_trig, bins)
    ax.errorbar(
        centers, eff_all, yerr=err_all,
        fmt="o-", ms=3, lw=1.2, capsize=2, color=COLORS[0],
        label="All events",
    )

    # Events with >= 1 cluster
    has_cl = ncl >= 1
    vtxz_cl = vtxz[has_cl]
    passed_cl = passed_trig[has_cl]
    centers_cl, eff_cl, err_cl, _, _ = compute_efficiency(vtxz_cl, passed_cl, bins)
    ax.errorbar(
        centers_cl, eff_cl, yerr=err_cl,
        fmt="s-", ms=3, lw=1.2, capsize=2, color=COLORS[1],
        label=r"Events with $\geq$1 EMCal cluster",
    )

    ax.set_xlabel("Truth vertex z [cm]")
    ax.set_ylabel("Trigger efficiency (N&S$\\geq$1)")
    ax.set_ylim(0, 1.05)
    ax.legend()
    ax.set_title("MBD trigger efficiency: inclusive vs cluster-selected events")
    ax.axhline(1.0, color="gray", ls="--", lw=0.5)
    ax.grid(axis="y", alpha=0.3)

    outpath = os.path.join(OUTDIR, f"mbd_trigger_eff_cluster_cut_{fname_suffix}.pdf")
    fig.savefig(outpath)
    plt.close(fig)
    print(f"  Saved: {outpath}")
    return outpath


def plot_vtxz_distributions(results, fname_suffix):
    """
    Show truth vertex z distribution with and without trigger,
    and the reco vertex z for comparison.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    vtxz_truth = results["vertexz_truth"]
    vtxz_reco = results["vertexz"]
    north = results["mbdnorthhit"]
    south = results["mbdsouthhit"]

    passed = (north >= 1) & (south >= 1)
    bins = np.linspace(-200, 200, 81)

    # Panel 1: truth vtxz with and without trigger
    ax1.hist(vtxz_truth, bins=bins, histtype="step", lw=1.5,
             color="black", label="All events")
    ax1.hist(vtxz_truth[passed], bins=bins, histtype="stepfilled",
             alpha=0.4, color=COLORS[0], label="Pass N&S$\\geq$1")
    ax1.hist(vtxz_truth[~passed], bins=bins, histtype="stepfilled",
             alpha=0.4, color=COLORS[1], label="Fail N&S$\\geq$1")
    ax1.set_xlabel("Truth vertex z [cm]")
    ax1.set_ylabel("Events / 5 cm")
    ax1.legend()
    ax1.set_title("Truth vertex z distribution")
    ax1.set_yscale("log")

    # Panel 2: reco vtxz for passing events (valid reco)
    valid_reco = passed & (np.abs(vtxz_reco) < 200) & (vtxz_reco != -9999)
    bins_reco = np.linspace(-60, 60, 49)
    ax2.hist(vtxz_truth[valid_reco], bins=bins_reco, histtype="step",
             lw=1.5, color="black", label="Truth z (trigger pass)")
    ax2.hist(vtxz_reco[valid_reco], bins=bins_reco, histtype="step",
             lw=1.5, color=COLORS[0], ls="--", label="Reco z (MBD vertex)")
    ax2.set_xlabel("Vertex z [cm]")
    ax2.set_ylabel("Events / 2.5 cm")
    ax2.legend()
    ax2.set_title("Truth vs reco vertex (trigger-passing, |z|<60)")

    fig.tight_layout()
    outpath = os.path.join(OUTDIR, f"mbd_vtxz_distributions_{fname_suffix}.pdf")
    fig.savefig(outpath)
    plt.close(fig)
    print(f"  Saved: {outpath}")
    return outpath


def print_summary(results, label):
    """Print efficiency summary for key vertex windows."""
    vtxz = results["vertexz_truth"]
    north = results["mbdnorthhit"]
    south = results["mbdsouthhit"]

    passed = (north >= 1) & (south >= 1)

    print(f"\n{'='*70}")
    print(f"  MBD Trigger Efficiency Summary — {label}")
    print(f"{'='*70}")
    print(f"  Total events: {len(vtxz):,}")
    print(f"  Pass N&S>=1:  {passed.sum():,} ({passed.mean()*100:.2f}%)")

    windows = [
        ("Full range", -300, 300),
        ("|z| < 200 cm", -200, 200),
        ("|z| < 100 cm", -100, 100),
        ("|z| < 60 cm", -60, 60),
        ("|z| < 30 cm", -30, 30),
        ("|z| < 10 cm", -10, 10),
    ]
    print(f"\n  {'Window':<20s}  {'Events':>10s}  {'Pass':>10s}  {'Eff':>8s}  {'Eff (N>=2&S>=2)':>16s}")
    print(f"  {'-'*20}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*16}")
    for wlabel, zlo, zhi in windows:
        mask = (vtxz >= zlo) & (vtxz < zhi)
        n_tot = mask.sum()
        n_pass = (mask & passed).sum()
        eff = n_pass / n_tot if n_tot > 0 else 0

        passed2 = (north >= 2) & (south >= 2)
        n_pass2 = (mask & passed2).sum()
        eff2 = n_pass2 / n_tot if n_tot > 0 else 0

        print(f"  {wlabel:<20s}  {n_tot:>10,}  {n_pass:>10,}  {eff:>8.4f}  {eff2:>16.4f}")

    # Efficiency variation within |z| < 30 cm
    mask30 = (vtxz > -30) & (vtxz < 30)
    bins_fine = np.linspace(-30, 30, 13)
    h_tot, _ = np.histogram(vtxz[mask30], bins=bins_fine)
    h_pass, _ = np.histogram(vtxz[mask30 & passed], bins=bins_fine)
    with np.errstate(divide="ignore", invalid="ignore"):
        eff_bins = np.where(h_tot > 0, h_pass / h_tot, np.nan)
    valid = np.isfinite(eff_bins)
    if valid.sum() > 0:
        print(f"\n  Within |z| < 30 cm (5 cm bins):")
        print(f"    Max efficiency: {np.nanmax(eff_bins):.4f}")
        print(f"    Min efficiency: {np.nanmin(eff_bins):.4f}")
        print(f"    Variation:      {np.nanmax(eff_bins) - np.nanmin(eff_bins):.4f}")
        print(f"    Relative var:   {(np.nanmax(eff_bins) - np.nanmin(eff_bins)) / np.nanmean(eff_bins) * 100:.2f}%")


def save_to_root(results, all_data, outpath):
    """Save key histograms to a ROOT file via uproot."""
    vtxz = results["vertexz_truth"]
    north = results["mbdnorthhit"]
    south = results["mbdsouthhit"]

    bins_full = VTX_BINS_WIDE
    passed_ns1 = (north >= 1) & (south >= 1)
    passed_ns2 = (north >= 2) & (south >= 2)

    h_total, _ = np.histogram(vtxz, bins=bins_full)
    h_pass1, _ = np.histogram(vtxz[passed_ns1], bins=bins_full)
    h_pass2, _ = np.histogram(vtxz[passed_ns2], bins=bins_full)

    with uproot.recreate(outpath) as f:
        f["h_vtxz_all"]       = (h_total.astype(np.float64), bins_full)
        f["h_vtxz_pass_ns1"]  = (h_pass1.astype(np.float64), bins_full)
        f["h_vtxz_pass_ns2"]  = (h_pass2.astype(np.float64), bins_full)

        # Per-sample
        for sample, data in all_data.items():
            vz = data["vertexz_truth"]
            p1 = (data["mbdnorthhit"] >= 1) & (data["mbdsouthhit"] >= 1)
            h_t, _ = np.histogram(vz, bins=bins_full)
            h_p, _ = np.histogram(vz[p1], bins=bins_full)
            f[f"h_vtxz_all_{sample}"]      = (h_t.astype(np.float64), bins_full)
            f[f"h_vtxz_pass_ns1_{sample}"] = (h_p.astype(np.float64), bins_full)

    print(f"  Saved ROOT file: {outpath}")


# ── Main ──────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("  MBD Trigger Efficiency Study — PPG12")
    print("=" * 60)

    # Load samples
    all_data = {}
    for sample in SAMPLES:
        print(f"\nLoading {sample} ...")
        all_data[sample] = load_sample(sample, max_files=1)

    # Combine all samples (unweighted — each has comparable statistics)
    combined = {}
    for key in BRANCHES:
        combined[key] = np.concatenate([all_data[s][key] for s in SAMPLES])
    print(f"\nCombined: {len(combined['vertexz_truth']):,} events")

    # Print summary tables
    print_summary(combined, "All photon samples combined")
    for sample in SAMPLES:
        print_summary(all_data[sample], sample)

    # Generate plots
    print("\n--- Generating plots ---")

    # 1. Efficiency vs vtxz — full range
    plot_efficiency_vs_vtxz(
        combined, VTX_BINS_FULL,
        "(all photon samples)", "full_range"
    )

    # 2. Efficiency vs vtxz — zoom to analysis region
    plot_efficiency_vs_vtxz(
        combined, VTX_BINS_ZOOM,
        r"(|z| < 60 cm, analysis region)", "zoom"
    )

    # 3. Efficiency vs vtxz — intermediate range
    plot_efficiency_vs_vtxz(
        combined, VTX_BINS_WIDE,
        r"(|z| < 150 cm)", "wide"
    )

    # 4. Per-sample comparison
    plot_efficiency_vs_vtxz_per_sample(all_data, VTX_BINS_WIDE, "wide")
    plot_efficiency_vs_vtxz_per_sample(all_data, VTX_BINS_ZOOM, "zoom")

    # 5. Hit asymmetry
    plot_hit_asymmetry(combined, "combined")

    # 6. Cluster cut comparison
    plot_efficiency_with_cluster_cut(combined, VTX_BINS_WIDE, "wide")
    plot_efficiency_with_cluster_cut(combined, VTX_BINS_ZOOM, "zoom")

    # 7. Vertex z distributions
    plot_vtxz_distributions(combined, "combined")

    # 8. Save histograms to ROOT
    save_to_root(
        combined, all_data,
        os.path.join(OUTDIR, "mbd_trigger_eff_study.root")
    )

    print("\n--- Done ---")


if __name__ == "__main__":
    main()
