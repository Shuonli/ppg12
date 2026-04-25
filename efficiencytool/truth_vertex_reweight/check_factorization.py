#!/usr/bin/env python3
"""
Stage 0 factorization sanity check.

Verifies, on a real double-interaction MC sample (and a nominal
single-interaction sample for cross-reference), that the two truth
vertices can be treated as independent and identically distributed.

Three gates:

  G1: |corr(z_h, z_mb)| < 0.02
  G2: |std(z_h) - std(z_mb)| < 2 σ_stat
  G3: |std(z_h)_double - std(z_h)_nominal| < 2 σ_stat

If any gate fails, the single-function w(z) assumption is broken and
the factorized fit cannot proceed.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import yaml

import fit_core as fc


def _load_both(double_path: str, single_path: str, branches: dict, n_events: int | None):
    double = fc.load_sample(
        fc.SampleSpec(name="double", path=double_path, xs=1.0, n_gen=1.0, kind="double_mc"),
        branches, n_events=n_events, require_mb=True,
    )
    single = fc.load_sample(
        fc.SampleSpec(name="single", path=single_path, xs=1.0, n_gen=1.0, kind="single_mc"),
        branches, n_events=n_events,
    )
    return double, single


def _std_err(arr: np.ndarray) -> tuple[float, float]:
    n = len(arr)
    if n < 2:
        return (float("nan"), float("nan"))
    s = float(arr.std(ddof=1))
    # standard error of sample std for Gaussian data: σ / sqrt(2(N-1))
    err = s / np.sqrt(2.0 * (n - 1))
    return s, err


def _corr_err(r: float, n: int) -> float:
    # Fisher-z approximation for the standard error of Pearson's r
    if n < 4:
        return float("nan")
    return (1.0 - r * r) / np.sqrt(n - 1)


def _save_scatter(zh: np.ndarray, zmb: np.ndarray, outpath: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    # 2D scatter (or hexbin for large N)
    ax = axes[0]
    if len(zh) > 50000:
        ax.hexbin(zh, zmb, gridsize=60, mincnt=1, cmap="viridis")
    else:
        ax.scatter(zh, zmb, s=3, alpha=0.3, color="C0")
    ax.axhline(0, color="k", lw=0.5, ls="--", alpha=0.5)
    ax.axvline(0, color="k", lw=0.5, ls="--", alpha=0.5)
    ax.set_xlabel(r"$z_{\rm truth}^{\rm hard}$ (cm)")
    ax.set_ylabel(r"$z_{\rm truth}^{\rm MB}$ (cm)")
    ax.set_title("Factorization scatter")
    ax.set_aspect("equal")
    lim = max(abs(zh).max(), abs(zmb).max()) * 1.05
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    # 1D overlay
    ax = axes[1]
    rng = (-200, 200)
    ax.hist(zh, bins=80, range=rng, histtype="step", color="C0", label="hard")
    ax.hist(zmb, bins=80, range=rng, histtype="step", color="C3", label="MB")
    ax.set_xlabel(r"$z_{\rm truth}$ (cm)")
    ax.set_ylabel("events")
    ax.set_title("1D truth distributions")
    ax.legend()

    fig.suptitle("Double-interaction truth vertex factorization check")
    fig.tight_layout()
    fig.savefig(outpath, dpi=120)
    plt.close(fig)
    print(f"  scatter + 1D plot saved to {outpath}")


def check(config: dict, double_name: str, single_name: str, out_pdf: Path, n_events: int | None) -> bool:
    branches = config["branches"]
    double_spec = config["inputs"]["double_mc"][double_name]
    single_spec = config["inputs"]["single_mc"][single_name]

    print("=== Stage 0: factorization sanity check ===\n")
    print(f"double MC: {double_name}  ({double_spec['path']})")
    print(f"single MC: {single_name}  ({single_spec['path']})\n")

    double, single = _load_both(double_spec["path"], single_spec["path"], branches, n_events)
    print(f"loaded {len(double.z_h)} double-MC events, {len(single.z_h)} single-MC events")

    # Filter sentinels — double MC events that did not get a second
    # truth vertex (shouldn't happen after reprocess, but be defensive).
    zh = double.z_h
    zmb = double.z_mb
    mask = (zh > fc.SENTINEL_F + 1) & (zmb > fc.SENTINEL_F + 1) & np.isfinite(zh) & np.isfinite(zmb)
    zh = zh[mask]
    zmb = zmb[mask]
    dropped = len(double.z_h) - len(zh)
    print(f"  {dropped} double-MC events dropped (missing z_h or z_mb)")

    zh_nom_mask = np.isfinite(single.z_h) & (single.z_h > fc.SENTINEL_F + 1)
    zh_nom = single.z_h[zh_nom_mask]

    # Stats
    std_zh, err_zh = _std_err(zh)
    std_zmb, err_zmb = _std_err(zmb)
    std_zh_nom, err_zh_nom = _std_err(zh_nom)
    mean_zh = float(zh.mean())
    mean_zmb = float(zmb.mean())

    print("\n--- Widths ---")
    print(f"  std(z_h, double)   = {std_zh:.3f} ± {err_zh:.3f} cm  (mean {mean_zh:+.3f})")
    print(f"  std(z_mb, double)  = {std_zmb:.3f} ± {err_zmb:.3f} cm  (mean {mean_zmb:+.3f})")
    print(f"  std(z_h, single)   = {std_zh_nom:.3f} ± {err_zh_nom:.3f} cm")

    r = float(np.corrcoef(zh, zmb)[0, 1])
    err_r = _corr_err(r, len(zh))
    print(f"\n--- Correlation ---")
    print(f"  corr(z_h, z_mb)    = {r:+.5f} ± {err_r:.5f}")

    # Gates
    print("\n--- Gates ---")

    # G1: |corr| < 0.02
    G1 = abs(r) < 0.02
    print(f"  [{'PASS' if G1 else 'FAIL'}]  G1: |corr| = {abs(r):.5f} < 0.02")

    # G2: |std(z_h) - std(z_mb)| < 2σ
    diff_hmb = abs(std_zh - std_zmb)
    err_diff = float(np.hypot(err_zh, err_zmb))
    G2 = diff_hmb < 2.0 * err_diff
    print(f"  [{'PASS' if G2 else 'FAIL'}]  G2: |Δstd(h, mb)| = {diff_hmb:.3f} < 2σ = {2*err_diff:.3f} cm")

    # G3: |std(z_h)_double - std(z_h)_single| < 2σ
    diff_ds = abs(std_zh - std_zh_nom)
    err_diff_ds = float(np.hypot(err_zh, err_zh_nom))
    G3 = diff_ds < 2.0 * err_diff_ds
    print(f"  [{'PASS' if G3 else 'FAIL'}]  G3: |Δstd(double, single)| = {diff_ds:.3f} < 2σ = {2*err_diff_ds:.3f} cm")

    all_pass = G1 and G2 and G3

    # Scatter plot — always produced for eyeball check
    try:
        _save_scatter(zh, zmb, out_pdf)
    except Exception as e:
        print(f"WARNING: failed to save scatter plot: {e}")

    print()
    if all_pass:
        print("ALL GATES PASS — factorized fit is valid; proceed to fit_truth_vertex_reweight.py")
    else:
        print("ONE OR MORE GATES FAILED — the single-function w(z) assumption is broken.")
        print("Escalate to the non-parametric 2D approach (not implemented in this directory).")

    return all_pass


def main(argv: list[str] | None = None) -> int:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(here / "config.yaml"),
                    help="path to config.yaml")
    ap.add_argument("--double-sample", default="jet12_double",
                    help="key under inputs.double_mc to use")
    ap.add_argument("--single-sample", default="jet12",
                    help="key under inputs.single_mc to use for cross-check")
    ap.add_argument("--n-events", type=int, default=None,
                    help="limit number of events per file (for quick tests)")
    ap.add_argument("--out-pdf", default=None,
                    help="output PDF for scatter plot (default: output/factorization_<sample>.pdf)")
    args = ap.parse_args(argv)

    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Resolve "TODO" placeholder paths for a quick failure
    dm = config["inputs"]["double_mc"][args.double_sample]["path"]
    sm = config["inputs"]["single_mc"][args.single_sample]["path"]
    for label, p in (("double", dm), ("single", sm)):
        if "TODO" in p or not p:
            print(f"ERROR: {label} sample path in config is still a placeholder: {p}")
            print("Fill in config.yaml with the reprocessed slimtree paths first.")
            return 2

    out_base = config["output"]["base_dir"]
    out_base_path = (here / out_base) if not Path(out_base).is_absolute() else Path(out_base)
    if args.out_pdf is None:
        out_pdf = out_base_path / f"factorization_{args.double_sample}.pdf"
    else:
        out_pdf = Path(args.out_pdf)

    ok = check(config, args.double_sample, args.single_sample, out_pdf, args.n_events)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
