#!/usr/bin/env python3
"""Per-variant ROC and BDT-score grid plots for post_preliminary_updates.

Produces, for each cluster node (split / nosplit):
  - roc_grid_{node}.pdf            : 4x3 grid, one panel per variant, 5 ROC curves
                                     (one per ET training bin)
  - bdt_score_grid_{node}.pdf      : 4x3 grid, one panel per variant, signal/bkg BDT
                                     score histograms in the 15-20 GeV ET bin

Same data pipeline as plot_isoet_vs_bdt_profile.py: streams MC samples from
/sphenix/user/shuhangli/ppg12/FunWithxgboost/{sample}/bdt_split.root, applies
training fiducial cuts.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import uproot
from sklearn.metrics import roc_auc_score, roc_curve

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

ET_BINS = [(6, 10), (10, 15), (15, 20), (20, 25), (25, 35)]
ET_BIN_COLORS = ["C0", "C1", "C2", "C3", "C4"]
BDTSCORE_ETBIN = (15, 20)  # histogram ET-bin for Figure 2 grid

ET_MIN_GLOBAL, ET_MAX_GLOBAL = 6.0, 40.0
ETA_ABS_MAX = 0.7


def load_sample(sample: str, entry_stop: int, node: str, variants: list[str]) -> dict:
    path = SIM_ROOT / sample / "bdt_split.root"
    if not path.exists():
        return None
    tree = uproot.open(str(path))[uproot.open(str(path)).keys()[0]]

    bdt_branches = [f"cluster_bdt_{node}_{v}" for v in variants]
    want = [f"cluster_Et_{node}", f"cluster_Eta_{node}"] + bdt_branches
    missing = [b for b in want if b not in tree.keys()]
    if missing:
        print(f"  [{sample}/{node}] missing: {missing[:2]}  skipping")
        return None

    t0 = time.time()
    arr = tree.arrays(want, entry_stop=entry_stop, library="ak")
    et = ak.flatten(arr[f"cluster_Et_{node}"]).to_numpy()
    eta = ak.flatten(arr[f"cluster_Eta_{node}"]).to_numpy()
    mask = ((et > ET_MIN_GLOBAL) & (et < ET_MAX_GLOBAL)
            & (np.abs(eta) < ETA_ABS_MAX))
    out = {"et": et[mask]}
    for v in variants:
        out[v] = ak.flatten(arr[f"cluster_bdt_{node}_{v}"]).to_numpy()[mask]
    print(f"  [{sample}/{node}] {mask.sum():>7d} clusters ({time.time()-t0:.1f}s)")
    return out


def concat(parts, variants):
    combined = {"et": np.concatenate([p["et"] for p in parts])}
    for v in variants:
        combined[v] = np.concatenate([p[v] for p in parts])
    return combined


def add_sphenix_header(fig, text="sPHENIX Simulation Internal"):
    fig.text(0.01, 0.995, text, style="italic", fontsize=11,
             ha="left", va="top", weight="bold")


def draw_variant_roc(ax, sig: dict, bkg: dict, variant: str,
                     show_legend: bool = False):
    for (lo, hi), col in zip(ET_BINS, ET_BIN_COLORS):
        ss = (sig["et"] >= lo) & (sig["et"] < hi)
        bb = (bkg["et"] >= lo) & (bkg["et"] < hi)
        if ss.sum() < 50 or bb.sum() < 50:
            continue
        y = np.concatenate([np.ones(ss.sum()), np.zeros(bb.sum())])
        s = np.concatenate([sig[variant][ss], bkg[variant][bb]])
        fpr, tpr, _ = roc_curve(y, s)
        auc = roc_auc_score(y, s)
        ax.plot(fpr, tpr, "-", color=col, lw=1.3,
                label=f"{lo}-{hi} (AUC={auc:.3f})")
    ax.plot([0, 1], [0, 1], "k--", lw=0.7, alpha=0.5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.01)
    ax.set_xlabel("FPR"); ax.set_ylabel("TPR")
    ax.grid(True, alpha=0.3)
    ax.set_title(variant, fontsize=10)
    if show_legend:
        ax.legend(fontsize=7, loc="lower right")


def draw_variant_bdtscore(ax, sig: dict, bkg: dict, variant: str,
                          et_lo: float, et_hi: float, show_legend: bool = False):
    ss = (sig["et"] >= et_lo) & (sig["et"] < et_hi)
    bb = (bkg["et"] >= et_lo) & (bkg["et"] < et_hi)
    if ss.sum() < 20 or bb.sum() < 20:
        ax.text(0.5, 0.5, "low stats", ha="center", va="center",
                transform=ax.transAxes, fontsize=9)
        ax.set_title(variant, fontsize=10)
        return
    bins = np.linspace(0, 1, 41)
    ax.hist(bkg[variant][bb], bins=bins, histtype="stepfilled",
            color="steelblue", alpha=0.55, density=True, label="Background")
    ax.hist(sig[variant][ss], bins=bins, histtype="step", color="darkorange",
            lw=1.8, density=True, label="Signal")
    ax.set_xlim(0, 1)
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("BDT score"); ax.set_ylabel("a.u.")
    ax.set_title(variant, fontsize=10)
    if show_legend:
        ax.legend(fontsize=8, loc="upper center")


def grid_figure(title: str, draw_fn, sig, bkg, outpath: Path, **kwargs):
    fig, axes = plt.subplots(3, 4, figsize=(18, 12))
    axes = axes.flatten()
    for i, v in enumerate(VARIANTS):
        draw_fn(axes[i], sig, bkg, v, show_legend=(i == 0), **kwargs)
    for j in range(len(VARIANTS), len(axes)):
        axes[j].axis("off")
    fig.suptitle(title, fontsize=14)
    add_sphenix_header(fig)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {outpath.name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--entries", type=int, default=400_000)
    ap.add_argument("--outdir", default="reports/figures/bdt_audit")
    ap.add_argument("--only", default=None, help="split|nosplit")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    nodes = [args.only] if args.only else list(CLUSTER_NODES.keys())

    for node_key in nodes:
        node = CLUSTER_NODES[node_key]
        print(f"\n=== {node_key} ({node}) ===")

        print("Loading signal …")
        sig_parts = [r for s in SIGNAL_SAMPLES
                     if (r := load_sample(s, args.entries, node, VARIANTS)) is not None]
        print("Loading background …")
        bkg_parts = [r for s in BACKGROUND_SAMPLES
                     if (r := load_sample(s, args.entries, node, VARIANTS)) is not None]
        if not sig_parts or not bkg_parts:
            print(f"  insufficient samples, skipping {node_key}")
            continue
        sig = concat(sig_parts, VARIANTS); bkg = concat(bkg_parts, VARIANTS)
        print(f"  combined sig={len(sig['et'])}, bkg={len(bkg['et'])}")

        grid_figure(
            f"ROC curves per variant — {node_key} ({node})",
            draw_variant_roc, sig, bkg,
            outdir / f"roc_grid_{node_key}.pdf",
        )
        lo, hi = BDTSCORE_ETBIN
        grid_figure(
            f"BDT score distributions ({lo} < ET < {hi} GeV) — {node_key} ({node})",
            draw_variant_bdtscore, sig, bkg,
            outdir / f"bdt_score_grid_{node_key}.pdf",
            et_lo=lo, et_hi=hi,
        )


if __name__ == "__main__":
    sys.exit(main())
