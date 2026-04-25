"""
Reweighting verification for PPG12 BDT training.

Loads the split- and nosplit-variant shapes files, reproduces the DataLoader
preprocessing (background subsampling) and KinematicReweighter pipeline as
applied during training, and produces multi-panel QA PDFs:

  reports/bdt_audit/reweight_qa_split.pdf
  reports/bdt_audit/reweight_qa_nosplit.pdf

Also prints flatness metrics to stdout.

Usage:
  python verify_reweighting.py [--variant split|nosplit|both] [--frac 1.0]
"""
from __future__ import annotations

import argparse
import os
import sys
from copy import deepcopy
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

REPO = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12"
FWX = os.path.join(REPO, "FunWithxgboost")
sys.path.insert(0, FWX)

from data_loader import DataLoader  # noqa: E402
from reweighting import KinematicReweighter  # noqa: E402

OUT_DIR = os.path.join(REPO, "reports", "bdt_audit")
os.makedirs(OUT_DIR, exist_ok=True)

JET_BASES = ["jet5", "jet12", "jet15", "jet20", "jet30", "jet40", "jet50"]
PHO_BASES = ["photon5", "photon10", "photon20"]


def existing_files(variant: str) -> Tuple[List[str], List[str]]:
    """Return (signal_paths, background_paths) that exist with > 1 line."""
    sig, bkg = [], []
    for s in PHO_BASES:
        p = os.path.join(FWX, f"shapes_{variant}_{s}.txt")
        if os.path.exists(p) and _nonempty(p):
            sig.append(p)
    for j in JET_BASES:
        p = os.path.join(FWX, f"shapes_{variant}_{j}.txt")
        if os.path.exists(p) and _nonempty(p):
            bkg.append(p)
    return sig, bkg


def _nonempty(path: str) -> bool:
    try:
        with open(path) as f:
            f.readline()  # header
            return f.readline().strip() != ""
    except OSError:
        return False


def build_config(base_cfg: Dict, sig: List[str], bkg: List[str]) -> Dict:
    cfg = deepcopy(base_cfg)
    cfg["data"]["use_single_file_set"] = True
    cfg["data"]["single_file_set"] = {"signal": sig, "background": bkg}
    return cfg


def load_with_pipeline(cfg: Dict, frac: float) -> pd.DataFrame:
    loader = DataLoader(cfg)
    df = loader._load_single_file_set_once()
    if frac < 1.0:
        df = df.sample(frac=frac, random_state=42).reset_index(drop=True)
    rew = KinematicReweighter(cfg)
    df_w = rew.apply_all_reweighting(df)
    return df_w


def _flat_residuals(values: np.ndarray, weights: np.ndarray, edges: np.ndarray) -> np.ndarray:
    counts, _ = np.histogram(values, bins=edges, weights=weights)
    nonzero = counts[counts > 0]
    if nonzero.size == 0:
        return np.zeros_like(counts)
    mean = nonzero.mean()
    return (counts - mean) / mean


def _spline_check(df: pd.DataFrame, label: int, n_bins: int = 20) -> Dict[str, float]:
    """Re-fit the same UnivariateSpline used inside KinematicReweighter and
    sample it on a fine grid to flag negative / oscillatory PDF values."""
    from scipy.interpolate import UnivariateSpline

    out = {}
    mask = df["label"] == label
    et = df.loc[mask, "cluster_Et"].values
    if len(et) == 0:
        return out
    edges = np.linspace(et.min(), et.max(), n_bins + 1)
    hist, _ = np.histogram(et, bins=edges, density=True)
    hist *= n_bins
    centres = 0.5 * (edges[:-1] + edges[1:])
    sp = UnivariateSpline(centres, hist, s=0.0)
    grid = np.linspace(centres[0], centres[-1], 2000)
    pdf = sp(grid)
    out["et_spline_min_pdf"] = float(pdf.min())
    out["et_spline_max_pdf"] = float(pdf.max())
    out["et_spline_neg_frac"] = float(np.mean(pdf < 0))
    out["et_spline_below_clip_frac"] = float(np.mean(pdf < 1e-3))
    return out


def make_qa_pdf(df: pd.DataFrame, variant: str, out_pdf: str) -> Dict:
    """Produce the multi-panel QA pdf and return numeric summary."""
    sig = df[df["label"] == 1]
    bkg = df[df["label"] == 0]

    et_edges = np.linspace(min(sig["cluster_Et"].min(), bkg["cluster_Et"].min()),
                           max(sig["cluster_Et"].max(), bkg["cluster_Et"].max()),
                           41)
    eta_edges = np.linspace(-0.7, 0.7, 29)
    vz_edges = np.linspace(min(sig["vertexz"].min(), bkg["vertexz"].min()),
                           max(sig["vertexz"].max(), bkg["vertexz"].max()),
                           41)
    res_edges = np.linspace(et_edges[0], et_edges[-1], 21)

    fig, axes = plt.subplots(3, 3, figsize=(18, 14))

    # Row 1: ET raw / ET weighted / ET ratio (sig/bkg)
    ax = axes[0, 0]
    ax.hist(sig["cluster_Et"], bins=et_edges, histtype="step", color="C0",
            density=True, label=f"signal raw ({len(sig)})")
    ax.hist(bkg["cluster_Et"], bins=et_edges, histtype="step", color="C3",
            density=True, label=f"background raw ({len(bkg)})")
    ax.set_xlabel("ET [GeV]"); ax.set_ylabel("a.u."); ax.legend(fontsize=8); ax.set_title("ET, raw")

    ax = axes[0, 1]
    ax.hist(sig["cluster_Et"], bins=et_edges, histtype="step", color="C0",
            weights=sig["weight"], density=True, label="signal weighted")
    ax.hist(bkg["cluster_Et"], bins=et_edges, histtype="step", color="C3",
            weights=bkg["weight"], density=True, label="background weighted")
    ax.set_xlabel("ET [GeV]"); ax.set_ylabel("a.u."); ax.legend(fontsize=8); ax.set_title("ET, weighted")

    ax = axes[0, 2]
    sig_h, _ = np.histogram(sig["cluster_Et"], bins=et_edges, weights=sig["weight"])
    bkg_h, _ = np.histogram(bkg["cluster_Et"], bins=et_edges, weights=bkg["weight"])
    ratio = np.divide(sig_h, bkg_h, out=np.full_like(sig_h, np.nan, dtype=float),
                      where=bkg_h > 0)
    centres = 0.5 * (et_edges[:-1] + et_edges[1:])
    ax.plot(centres, ratio, "o-")
    ax.axhline(1.0, color="grey", ls="--")
    ax.set_xlabel("ET [GeV]"); ax.set_ylabel("sig/bkg weighted")
    ax.set_title("class-balance per ET bin (target = 1)")

    # Row 2: eta raw / eta weighted / vertex z (raw + weighted)
    ax = axes[1, 0]
    ax.hist(sig["cluster_Eta"], bins=eta_edges, histtype="step", color="C0",
            density=True, label="signal raw")
    ax.hist(bkg["cluster_Eta"], bins=eta_edges, histtype="step", color="C3",
            density=True, label="background raw")
    ax.set_xlabel("eta"); ax.legend(fontsize=8); ax.set_title("eta, raw")

    ax = axes[1, 1]
    ax.hist(sig["cluster_Eta"], bins=eta_edges, histtype="step", color="C0",
            weights=sig["weight"], density=True, label="signal weighted")
    ax.hist(bkg["cluster_Eta"], bins=eta_edges, histtype="step", color="C3",
            weights=bkg["weight"], density=True, label="background weighted")
    ax.set_xlabel("eta"); ax.legend(fontsize=8); ax.set_title("eta, weighted")

    ax = axes[1, 2]
    ax.hist(sig["vertexz"], bins=vz_edges, histtype="step", color="C0",
            density=True, label="signal raw", linestyle=":")
    ax.hist(bkg["vertexz"], bins=vz_edges, histtype="step", color="C3",
            density=True, label="background raw", linestyle=":")
    ax.hist(sig["vertexz"], bins=vz_edges, histtype="step", color="C0",
            weights=sig["weight"], density=True, label="signal weighted")
    ax.hist(bkg["vertexz"], bins=vz_edges, histtype="step", color="C3",
            weights=bkg["weight"], density=True, label="background weighted")
    ax.set_xlabel("vertexz [cm]"); ax.legend(fontsize=7)
    ax.set_title("vertex z (vertex_reweight = OFF)")

    # Row 3: weights log-y for sig + bkg, and flatness residuals
    ax = axes[2, 0]
    ax.hist(sig["weight"], bins=80, histtype="step", color="C0",
            label=f"signal w (n={len(sig)})")
    ax.hist(bkg["weight"], bins=80, histtype="step", color="C3",
            label=f"background w (n={len(bkg)})")
    ax.set_yscale("log")
    ax.set_xlabel("weight"); ax.legend(fontsize=8); ax.set_title("weight distribution")

    ax = axes[2, 1]
    ax.bar(np.arange(len(res_edges) - 1) - 0.2,
           _flat_residuals(sig["cluster_Et"].values, sig["weight"].values, res_edges),
           width=0.4, color="C0", label="signal")
    ax.bar(np.arange(len(res_edges) - 1) + 0.2,
           _flat_residuals(bkg["cluster_Et"].values, bkg["weight"].values, res_edges),
           width=0.4, color="C3", label="background")
    ax.axhline(0, color="grey", lw=0.7)
    ax.axhline(0.1, color="grey", lw=0.5, ls="--"); ax.axhline(-0.1, color="grey", lw=0.5, ls="--")
    ax.set_xlabel("ET bin index"); ax.set_ylabel("(N - mean) / mean")
    ax.legend(fontsize=8); ax.set_title("ET flatness residuals (target |.| < 0.1)")

    ax = axes[2, 2]
    eta_res_edges = np.linspace(-0.7, 0.7, 21)
    ax.bar(np.arange(len(eta_res_edges) - 1) - 0.2,
           _flat_residuals(sig["cluster_Eta"].values, sig["weight"].values, eta_res_edges),
           width=0.4, color="C0", label="signal")
    ax.bar(np.arange(len(eta_res_edges) - 1) + 0.2,
           _flat_residuals(bkg["cluster_Eta"].values, bkg["weight"].values, eta_res_edges),
           width=0.4, color="C3", label="background")
    ax.axhline(0, color="grey", lw=0.7)
    ax.axhline(0.1, color="grey", lw=0.5, ls="--"); ax.axhline(-0.1, color="grey", lw=0.5, ls="--")
    ax.set_xlabel("eta bin index"); ax.set_ylabel("(N - mean) / mean")
    ax.legend(fontsize=8); ax.set_title("eta flatness residuals")

    fig.suptitle(f"Reweighting QA — variant = {variant}", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_pdf)
    plt.close(fig)

    summary = _summarize(df, variant, et_edges, eta_edges)
    summary["et_flat_res_max_sig"] = float(np.max(np.abs(
        _flat_residuals(sig["cluster_Et"].values, sig["weight"].values, res_edges))))
    summary["et_flat_res_max_bkg"] = float(np.max(np.abs(
        _flat_residuals(bkg["cluster_Et"].values, bkg["weight"].values, res_edges))))
    summary["eta_flat_res_max_sig"] = float(np.max(np.abs(
        _flat_residuals(sig["cluster_Eta"].values, sig["weight"].values,
                        np.linspace(-0.7, 0.7, 21)))))
    summary["eta_flat_res_max_bkg"] = float(np.max(np.abs(
        _flat_residuals(bkg["cluster_Eta"].values, bkg["weight"].values,
                        np.linspace(-0.7, 0.7, 21)))))
    return summary


def _summarize(df: pd.DataFrame, variant: str, et_edges, eta_edges) -> Dict:
    sig = df[df["label"] == 1]
    bkg = df[df["label"] == 0]
    s = {
        "variant": variant,
        "n_sig": int(len(sig)),
        "n_bkg": int(len(bkg)),
    }
    for name, sub in [("sig", sig), ("bkg", bkg)]:
        w = sub["weight"].values
        s.update({
            f"w_{name}_min": float(w.min()),
            f"w_{name}_p50": float(np.median(w)),
            f"w_{name}_mean": float(w.mean()),
            f"w_{name}_p99": float(np.percentile(w, 99)),
            f"w_{name}_max": float(w.max()),
            f"w_{name}_sum": float(w.sum()),
        })
    s["spline_sig"] = _spline_check(df, 1)
    s["spline_bkg"] = _spline_check(df, 0)
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", choices=["split", "nosplit", "both"], default="both")
    ap.add_argument("--frac", type=float, default=1.0,
                    help="row fraction sub-sample for speed/memory")
    args = ap.parse_args()

    with open(os.path.join(FWX, "config.yaml")) as f:
        base_cfg = yaml.safe_load(f)

    variants = ["split", "nosplit"] if args.variant == "both" else [args.variant]
    summaries = {}
    for v in variants:
        sig, bkg = existing_files(v)
        print(f"\n=== variant {v} ===")
        print(f"  signal files     : {[os.path.basename(p) for p in sig]}")
        print(f"  background files : {[os.path.basename(p) for p in bkg]}")
        cfg = build_config(base_cfg, sig, bkg)
        df = load_with_pipeline(cfg, args.frac)
        out_pdf = os.path.join(OUT_DIR, f"reweight_qa_{v}.pdf")
        s = make_qa_pdf(df, v, out_pdf)
        summaries[v] = s
        print(f"  pdf written      : {out_pdf}")
        print(f"  N signal / bkg   : {s['n_sig']} / {s['n_bkg']}")
        print(f"  weight (sig)     min={s['w_sig_min']:.3g} mean={s['w_sig_mean']:.3g} "
              f"p99={s['w_sig_p99']:.3g} max={s['w_sig_max']:.3g}")
        print(f"  weight (bkg)     min={s['w_bkg_min']:.3g} mean={s['w_bkg_mean']:.3g} "
              f"p99={s['w_bkg_p99']:.3g} max={s['w_bkg_max']:.3g}")
        print(f"  ET flat |res|max sig={s['et_flat_res_max_sig']:.3f} "
              f"bkg={s['et_flat_res_max_bkg']:.3f}")
        print(f"  eta flat |res|max sig={s['eta_flat_res_max_sig']:.3f} "
              f"bkg={s['eta_flat_res_max_bkg']:.3f}")
        print(f"  ET spline (sig)  : {s['spline_sig']}")
        print(f"  ET spline (bkg)  : {s['spline_bkg']}")

    return summaries


if __name__ == "__main__":
    main()
