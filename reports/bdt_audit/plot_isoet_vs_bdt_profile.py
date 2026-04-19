#!/usr/bin/env python3
"""Profile plots of recoIsoET vs BDT score for every trained variant.

Swaps axes of the original isoET-vs-BDT scatter so BDT score is on X, recoIsoET on Y,
and overlays a binned profile of the mean recoIsoET at each BDT score bin. One plot
per (variant, cluster-node) combination (22 total), plus a per-cluster-node grid
summary.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import uproot

SIM_ROOT = Path("/sphenix/user/shuhangli/ppg12/FunWithxgboost")

SIGNAL_SAMPLES = ["photon5", "photon10", "photon20"]
BACKGROUND_SAMPLES = ["jet10", "jet15", "jet20", "jet30"]

VARIANTS = [
    "base", "base_vr", "base_v0", "base_v1", "base_v2", "base_v3",
    "base_E", "base_v0E", "base_v1E", "base_v2E", "base_v3E",
]

CLUSTER_NODES = {
    "split":    "CLUSTERINFO_CEMC",
    "nosplit":  "CLUSTERINFO_CEMC_NO_SPLIT",
}

# Matches training range
ET_MIN, ET_MAX = 6.0, 40.0
ETA_ABS_MAX = 0.7


def load_sample(sample: str, entry_stop: int, node: str, variants: list[str]) -> dict:
    """Return flat numpy arrays for one MC sample, selecting ET/|eta| cuts."""
    path = SIM_ROOT / sample / "bdt_split.root"
    if not path.exists():
        return None
    tree = uproot.open(str(path))[uproot.open(str(path)).keys()[0]]

    bdt_branches = [f"cluster_bdt_{node}_{v}" for v in variants]
    want = [
        f"cluster_Et_{node}",
        f"cluster_Eta_{node}",
        f"cluster_iso_04_{node}",
    ] + bdt_branches

    missing = [b for b in want if b not in tree.keys()]
    if missing:
        print(f"  [{sample}/{node}] missing branches: {missing[:3]} ... skipping")
        return None

    t0 = time.time()
    arr = tree.arrays(want, entry_stop=entry_stop, library="ak")
    et = ak.flatten(arr[f"cluster_Et_{node}"]).to_numpy()
    eta = ak.flatten(arr[f"cluster_Eta_{node}"]).to_numpy()
    iso = ak.flatten(arr[f"cluster_iso_04_{node}"]).to_numpy()
    mask = (et > ET_MIN) & (et < ET_MAX) & (np.abs(eta) < ETA_ABS_MAX) & np.isfinite(iso)

    out = {"et": et[mask], "eta": eta[mask], "iso": iso[mask]}
    for v in variants:
        b = ak.flatten(arr[f"cluster_bdt_{node}_{v}"]).to_numpy()
        out[v] = b[mask]

    print(f"  [{sample}/{node}] {mask.sum():>7d} clusters "
          f"({len(et):>7d} raw, {time.time()-t0:.1f}s)")
    return out


def concat_samples(results: list[dict], variants: list[str]) -> dict:
    combined = {"iso": np.concatenate([r["iso"] for r in results])}
    for v in variants:
        combined[v] = np.concatenate([r[v] for r in results])
    return combined


def profile(bdt: np.ndarray, iso: np.ndarray, nbins: int = 30,
            bdt_range=(0.0, 1.0)) -> tuple:
    """Mean + SEM of iso per BDT bin."""
    edges = np.linspace(*bdt_range, nbins + 1)
    centers = 0.5 * (edges[1:] + edges[:-1])
    means = np.full(nbins, np.nan)
    sems = np.full(nbins, np.nan)
    counts = np.zeros(nbins, dtype=int)
    idx = np.clip(np.digitize(bdt, edges) - 1, 0, nbins - 1)
    for i in range(nbins):
        sel = idx == i
        n = sel.sum()
        counts[i] = n
        if n >= 10:
            means[i] = iso[sel].mean()
            sems[i] = iso[sel].std(ddof=1) / np.sqrt(n)
    return centers, means, sems, counts


def plot_single_variant(ax, sig: dict, bkg: dict, variant: str,
                        iso_ylim=(-2, 35), show_legend=True, title_prefix=""):
    """Draw scatter + profile for one variant on an existing Axes."""
    # sparse scatter
    rng = np.random.default_rng(42)
    n_scatter = 4000
    for data, label, color, alpha in [
        (bkg, "Background", "red", 0.25),
        (sig, "Signal", "blue", 0.35),
    ]:
        x = data[variant]
        y = data["iso"]
        n = len(x)
        if n > n_scatter:
            sel = rng.choice(n, size=n_scatter, replace=False)
            x, y = x[sel], y[sel]
        ax.scatter(x, y, s=2, alpha=alpha, color=color, label=label, rasterized=True)

    for data, label, color, marker in [
        (sig, "Signal (profile)", "navy", "o"),
        (bkg, "Background (profile)", "darkred", "s"),
    ]:
        c, m, e, n = profile(data[variant], data["iso"])
        ok = ~np.isnan(m)
        ax.errorbar(c[ok], m[ok], yerr=e[ok], fmt=marker + "-",
                    color=color, lw=1.5, ms=4, capsize=2, label=label)

    # Pearson r on combined (to reproduce original annotation)
    all_bdt = np.concatenate([sig[variant], bkg[variant]])
    all_iso = np.concatenate([sig["iso"], bkg["iso"]])
    r = np.corrcoef(all_bdt, all_iso)[0, 1] if len(all_bdt) > 10 else np.nan

    ax.set_xlim(0, 1)
    ax.set_ylim(*iso_ylim)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("BDT Score")
    ax.set_ylabel("recoisoET (GeV)")
    ax.set_title(f"{title_prefix}{variant}")
    ax.text(0.05, 0.95, f"Pearson r = {r:+.3f}", transform=ax.transAxes,
            fontsize=9, verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))
    if show_legend:
        ax.legend(loc="upper right", fontsize=8)


def add_sphenix_header(fig, text="sPHENIX Simulation Internal"):
    fig.text(0.01, 0.99, text, style="italic", fontsize=11,
             ha="left", va="top", weight="bold")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--entries", type=int, default=150_000,
                    help="events per MC file (default 150k)")
    ap.add_argument("--outdir", default="reports/figures/bdt_audit")
    ap.add_argument("--only", default=None, help="restrict to one cluster node: split|nosplit")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    nodes = [args.only] if args.only else list(CLUSTER_NODES.keys())

    for node_key in nodes:
        node = CLUSTER_NODES[node_key]
        print(f"\n=== cluster node: {node_key} ({node}) ===")

        print("Loading signal samples …")
        sig_parts = []
        for s in SIGNAL_SAMPLES:
            r = load_sample(s, args.entries, node, VARIANTS)
            if r is not None:
                sig_parts.append(r)

        print("Loading background samples …")
        bkg_parts = []
        for s in BACKGROUND_SAMPLES:
            r = load_sample(s, args.entries, node, VARIANTS)
            if r is not None:
                bkg_parts.append(r)

        if not sig_parts or not bkg_parts:
            print(f"  {node_key}: not enough samples loaded, skipping")
            continue

        sig = concat_samples(sig_parts, VARIANTS)
        bkg = concat_samples(bkg_parts, VARIANTS)
        print(f"  combined: signal={len(sig['iso'])}, bkg={len(bkg['iso'])}")

        # Per-variant individual plots
        for v in VARIANTS:
            fig, ax = plt.subplots(figsize=(6.5, 5.5))
            plot_single_variant(ax, sig, bkg, v, title_prefix=f"{node_key}: ")
            add_sphenix_header(fig)
            fig.tight_layout(rect=(0, 0, 1, 0.96))
            out = outdir / f"profile_isoET_vs_BDT_{v}_{node_key}.pdf"
            fig.savefig(out, dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"  wrote {out.name}")

        # Grid summary (4x3, 11 panels + 1 spare)
        fig, axes = plt.subplots(3, 4, figsize=(18, 13))
        axes = axes.flatten()
        for i, v in enumerate(VARIANTS):
            plot_single_variant(axes[i], sig, bkg, v, show_legend=(i == 0),
                                title_prefix="")
        # hide spare
        for j in range(len(VARIANTS), len(axes)):
            axes[j].axis("off")
        fig.suptitle(f"recoIsoET vs BDT score — {node_key} "
                     f"({CLUSTER_NODES[node_key]})", fontsize=14)
        add_sphenix_header(fig)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        out = outdir / f"profile_isoET_vs_BDT_grid_{node_key}.pdf"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {out.name}")


if __name__ == "__main__":
    sys.exit(main())
