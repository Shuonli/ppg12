#!/usr/bin/env python3
"""
6-panel closure plot for the truth vertex reweight fit.

Panels:
  1. Reco vertex: data vs mixed MC (linear, within 2× closure window)
  2. Reco vertex: same, log y, full binning range
  3. Reweight functions w(z) (log)
  4. Residual D/M per bin, with Poisson error bars and ±1σ stats band
  5. Convergence trajectory: χ² and max|r−1| vs iteration
  6. Pull distribution (D−M)/σ inside the closure window
"""
from __future__ import annotations

import argparse
import logging
import math
import sys
from pathlib import Path
from typing import Tuple

import numpy as np
import yaml

import fit_core as fc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sphenix_style import (
    set_sphenix_rcparams, sphenix_label, sphenix_period_text,
)
set_sphenix_rcparams()


LOG = logging.getLogger("plot_closure")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _gauss(z: np.ndarray, sigma: float, norm: float = 1.0) -> np.ndarray:
    return norm * np.exp(-0.5 * (z / sigma) ** 2) / (sigma * math.sqrt(2.0 * math.pi))


def _step_hist(ax, centers, values, **kw):
    edges = np.zeros(len(centers) + 1)
    if len(centers) > 1:
        dz = centers[1] - centers[0]
    else:
        dz = 1.0
    edges[:-1] = centers - dz / 2.0
    edges[-1] = centers[-1] + dz / 2.0
    ax.step(edges, np.concatenate([values, [values[-1]]]), where="post", **kw)


def _ratio_and_sigma(
    D: np.ndarray, M: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (ratio, σ_ratio) bin-by-bin for D/M with Poisson errors.
    σ(D/M) = (D/M) · sqrt(1/D + 1/M). Bins with D=0 or M=0 get ratio=NaN."""
    ratio = np.full_like(D, np.nan, dtype=np.float64)
    sigma = np.full_like(D, np.nan, dtype=np.float64)
    valid = (D > 0) & (M > 0)
    ratio[valid] = D[valid] / M[valid]
    sigma[valid] = ratio[valid] * np.sqrt(1.0 / D[valid] + 1.0 / M[valid])
    return ratio, sigma


def _summary_stat(resid: np.ndarray, mask: np.ndarray) -> Tuple[float, float, float]:
    """Return (max_abs_dev, rms_dev, mean_dev) inside the mask."""
    r = resid[mask]
    r = r[np.isfinite(r)]
    if len(r) == 0:
        return float("nan"), float("nan"), float("nan")
    dev = r - 1.0
    return float(np.max(np.abs(dev))), float(np.sqrt(np.mean(dev ** 2))), float(np.mean(dev))


# ---------------------------------------------------------------------------
# Main plot
# ---------------------------------------------------------------------------

def plot(out_dir: Path, config: dict, period_name: str, cache_dir: Path,
         sample_match=None, z_cut_override=None) -> int:
    yaml_path = out_dir / "fit_result.yaml"
    if not yaml_path.exists():
        LOG.error("missing fit_result.yaml at %s", yaml_path)
        return 2
    with open(yaml_path) as f:
        doc = yaml.safe_load(f)

    bnode = doc["binning"]
    centers = np.asarray(bnode["centers"])
    z_min, z_max = bnode["z_min"], bnode["z_max"]
    n_bins = bnode["n_bins"]
    binning = fc.Binning(z_min=z_min, z_max=z_max, n_bins=n_bins)
    bin_width = float(centers[1] - centers[0]) if len(centers) > 1 else 1.0

    s1 = doc["stage1_parametric"]
    s2 = doc["stage2_residual"]
    sigma_gen = float(s1["sigma_gen"])
    sigma_tgt = float(s1["sigma_tgt"])
    w_param = np.asarray(s1["w_values"])
    c_raw = np.asarray(s2["c_raw"])
    c_sm = np.asarray(s2["c_smoothed"])
    w_total = np.asarray(s2["w_total"])
    resid_final = np.asarray(s2["residual_final"])

    iter_doc = doc.get("iterative")
    has_iter = iter_doc is not None
    w_iter = w_seed = resid_iter = None
    iter_history = None
    sigma_seed = None
    iter_converged = False
    iter_max_res = float("nan")
    if has_iter:
        w_iter = np.asarray(iter_doc["w_iterative"])
        w_seed = np.asarray(iter_doc["w_seed"])
        resid_iter = np.asarray(iter_doc["residual_iterative"])
        sigma_seed = float(iter_doc["sigma_seed_tgt"])
        iter_history = iter_doc["history"]
        iter_converged = bool(iter_doc["converged"])
        iter_max_res = float(iter_doc["max_residual_in_window"])

    # Reload inputs so we can rebuild M histograms and compute proper errors
    from fit_truth_vertex_reweight import load_period_inputs
    LOG.info("reloading inputs for plotting ...")
    inputs = load_period_inputs(
        config, period_name, cache_dir=cache_dir,
        sample_match=sample_match, n_events=None,
        z_cut_override=z_cut_override,
    )
    D = inputs.D.astype(np.float64)
    f_s, f_d = inputs.period.f_single, inputs.period.f_double
    z_cut = float(inputs.period.z_cut)
    closure_mask = np.abs(centers) < z_cut

    def _w_unit(z):
        return np.ones_like(z, dtype=np.float64)

    def _w_param_fn(z):
        return fc.parametric_weight(z, sigma_gen, sigma_tgt)

    c_func = fc.c_interpolator(centers, c_sm)

    M_noRW = fc.mixed_mc_histogram(
        inputs.single_samples, inputs.double_samples, binning, f_s, f_d, _w_unit
    )
    M_noRW *= D.sum() / max(M_noRW.sum(), 1e-30)

    M_s1 = fc.mixed_mc_histogram(
        inputs.single_samples, inputs.double_samples, binning, f_s, f_d, _w_param_fn
    )
    M_s1 *= D.sum() / max(M_s1.sum(), 1e-30)

    # Stage 1+2 is not shown in the nominal plots; skip the forward
    # model to save time. Kept as zeros for downstream compatibility.
    M_s2 = np.zeros_like(M_s1)

    M_seed = None
    M_iter = None
    if has_iter:
        def _w_seed_fn(z):
            return fc.parametric_weight(z, sigma_gen, sigma_seed)
        M_seed = fc.mixed_mc_histogram(
            inputs.single_samples, inputs.double_samples, binning, f_s, f_d, _w_seed_fn,
        )
        M_seed *= D.sum() / max(M_seed.sum(), 1e-30)
        w_iter_func = fc.c_interpolator(centers, w_iter)
        M_iter = fc.mixed_mc_histogram(
            inputs.single_samples, inputs.double_samples, binning, f_s, f_d, w_iter_func,
        )
        M_iter *= D.sum() / max(M_iter.sum(), 1e-30)

    # Bin-by-bin ratios with Poisson errors
    r_s1, e_s1 = _ratio_and_sigma(D, M_s1)
    r_s2, e_s2 = _ratio_and_sigma(D, M_s2)
    r_seed, e_seed = (_ratio_and_sigma(D, M_seed) if M_seed is not None else (None, None))
    if has_iter:
        r_iter, e_iter = _ratio_and_sigma(D, M_iter)
    else:
        r_iter, e_iter = r_s2, e_s2

    # Summary stats inside the closure window
    max_s1, rms_s1, mean_s1 = _summary_stat(r_s1, closure_mask)
    max_s2, rms_s2, mean_s2 = _summary_stat(r_s2, closure_mask)
    max_seed = _summary_stat(r_seed, closure_mask)[0] if r_seed is not None else float("nan")
    max_it, rms_it, mean_it = _summary_stat(r_iter, closure_mask)

    # ------------------------------------------------------------------
    # Figure: 2×3 layout
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.3)
    ax1 = fig.add_subplot(gs[0, 0])   # log reco closure
    ax2 = fig.add_subplot(gs[0, 1])   # linear reco closure (zoom)
    ax3 = fig.add_subplot(gs[0, 2])   # reweight functions
    ax4 = fig.add_subplot(gs[1, 0:2]) # residual D/M with error bars (span 2)
    ax5 = fig.add_subplot(gs[1, 2])   # convergence trajectory

    # =============================================================
    # Panel 1 — reco closure, log y, full range
    # =============================================================
    _step_hist(ax1, centers, D, color="k", lw=1.5, label="data")
    _step_hist(ax1, centers, M_noRW, color="C7", lw=1.0, ls=":", label="MC, no reweight")
    if M_seed is not None:
        _step_hist(ax1, centers, M_seed, color="C1", lw=1.0, ls="-.",
                   label=fr"MC, seed Gaussian ($\sigma_{{\rm seed}}={sigma_seed:.1f}$)")
    _step_hist(ax1, centers, M_s1, color="C0", lw=1.0, ls="--",
               label=fr"MC, parametric Gaussian ($\sigma_{{\rm tgt}}={sigma_tgt:.1f}$)")
    if has_iter:
        _step_hist(ax1, centers, M_iter, color="C3", lw=1.6, label="MC, iterative")
    ax1.axvspan(-z_cut, +z_cut, color="C2", alpha=0.08, zorder=0)
    ax1.axvline(-z_cut, color="C2", lw=0.8, ls="-.", alpha=0.6)
    ax1.axvline(+z_cut, color="C2", lw=0.8, ls="-.", alpha=0.6,
                label=fr"closure ±{z_cut:.0f} cm")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$z_{\rm reco}$ (cm)")
    ax1.set_ylabel("events")
    ax1.set_title("Reco vertex (log)")
    ax1.legend(fontsize=7.5, loc="lower center")
    ax1.grid(alpha=0.3)
    ax1.set_xlim(z_min, z_max)

    # =============================================================
    # Panel 2 — reco closure, linear y, zoomed to ±2 z_cut
    # =============================================================
    _step_hist(ax2, centers, D, color="k", lw=1.5, label="data")
    if has_iter:
        _step_hist(ax2, centers, M_iter, color="C3", lw=1.6, label="MC, iterative")
    _step_hist(ax2, centers, M_s1, color="C0", lw=1.0, ls="--",
               label="MC, parametric Gaussian")
    if M_seed is not None:
        _step_hist(ax2, centers, M_seed, color="C1", lw=1.0, ls="-.",
                   label="MC, seed Gaussian")
    ax2.axvspan(-z_cut, +z_cut, color="C2", alpha=0.08, zorder=0)
    ax2.set_xlabel(r"$z_{\rm reco}$ (cm)")
    ax2.set_ylabel("events")
    ax2.set_title("Reco vertex (linear, zoom)")
    ax2.legend(fontsize=8, loc="upper right")
    ax2.grid(alpha=0.3)
    zoom_lim = max(1.5 * z_cut, 1.0)
    zoom_lim = min(zoom_lim, z_max)
    ax2.set_xlim(-zoom_lim, +zoom_lim)

    # =============================================================
    # Panel 3 — reweight functions
    # =============================================================
    _step_hist(ax3, centers, w_param, color="C0", lw=1.0, ls="--",
               label=fr"$w_{{\rm param}}$ ($\sigma_{{\rm tgt}}={sigma_tgt:.1f}$)")
    if has_iter:
        _step_hist(ax3, centers, w_seed, color="C7", lw=1.0, ls=":",
                   label=fr"seed Gaussian ($\sigma={sigma_seed:.1f}$)")
        _step_hist(ax3, centers, w_iter, color="C3", lw=1.6, label=r"$w_{\rm iterative}$")
    ax3.axvspan(-z_cut, +z_cut, color="C2", alpha=0.08, zorder=0)
    ax3.axhline(1.0, color="k", lw=0.5, ls=":", alpha=0.5)
    ax3.set_yscale("log")
    ax3.set_xlabel(r"$z$ (cm)")
    ax3.set_ylabel("weight")
    ax3.set_title(fr"Reweight functions  ($\sigma_{{\rm gen}}={sigma_gen:.1f}$ cm)")
    ax3.legend(fontsize=8, loc="lower center")
    ax3.grid(alpha=0.3, which="both")
    ax3.set_xlim(z_min, z_max)

    # =============================================================
    # Panel 4 — residual D/M with Poisson error bars
    # =============================================================
    # Show the σ=1% / 2% reference bands as a proxy for "acceptable closure"
    ax4.fill_between(
        [z_min, z_max], 0.98, 1.02,
        color="gray", alpha=0.15, zorder=0, label="±2% band",
    )
    ax4.fill_between(
        [z_min, z_max], 0.99, 1.01,
        color="gray", alpha=0.25, zorder=0, label="±1% band",
    )
    ax4.axhline(1.0, color="k", lw=0.6, ls=":", alpha=0.6)
    # Closure window shading
    ax4.axvspan(-z_cut, +z_cut, color="C2", alpha=0.08, zorder=0)

    # Parametric Gaussian (stage 1, best-fit σ_tgt from χ² scan) and
    # seed Gaussian (σ_seed = std(z_r^data), starting point of iteration)
    # shown as thin step lines for context (no error bars).
    # Stage 1+2 (parametric + residual correction) is intentionally not
    # shown: it is never used in the nominal pipeline.
    if r_seed is not None:
        _step_hist(ax4, centers, r_seed, color="C1", lw=1.0, ls="-.",
                   label=f"seed Gaussian (max|r−1|={max_seed:.3f})")
    _step_hist(ax4, centers, r_s1, color="C0", lw=1.0, ls="--",
               label=f"parametric Gaussian (max|r−1|={max_s1:.3f})")

    # Iterative — points with stat error bars (this is the main result)
    if has_iter:
        valid = np.isfinite(r_iter) & np.isfinite(e_iter)
        ax4.errorbar(
            centers[valid], r_iter[valid], yerr=e_iter[valid],
            fmt="o", color="C3", ecolor="C3", markersize=3,
            elinewidth=1.0, capsize=0,
            label=f"iterative (max|r−1|={max_it:.3f}, rms={rms_it:.3f})",
        )

    ax4.axvline(-z_cut, color="C2", lw=0.8, ls="-.", alpha=0.7)
    ax4.axvline(+z_cut, color="C2", lw=0.8, ls="-.", alpha=0.7)
    ax4.set_xlabel(r"$z_{\rm reco}$ (cm)")
    ax4.set_ylabel(r"$D / M$")

    title_stats = (
        f"iter: max|r−1|={max_it:.3f}, rms={rms_it:.3f}, mean={mean_it:+.3f}"
    )
    ax4.set_title(f"Residual D/M in closure window   [{title_stats}]", fontsize=10)

    # Adaptive y-range based on iterative residuals in the closure window
    r_in = r_iter[closure_mask] if has_iter else r_s1[closure_mask]
    r_in = r_in[np.isfinite(r_in)]
    if len(r_in) > 0:
        lo = min(0.9, float(np.min(r_in)) - 0.02)
        hi = max(1.1, float(np.max(r_in)) + 0.02)
        ax4.set_ylim(lo, hi)
    else:
        ax4.set_ylim(0.5, 1.5)
    ax4.set_xlim(-1.5 * z_cut, +1.5 * z_cut)
    ax4.legend(fontsize=8, loc="lower center", ncol=3)
    ax4.grid(alpha=0.3)

    # =============================================================
    # Panel 5 — convergence trajectory (χ² and max|r−1| vs iteration)
    # =============================================================
    if has_iter and iter_history:
        iters = [h["iteration"] for h in iter_history]
        chi2_hist = [h["chi2_in_window"] for h in iter_history]
        max_res_hist = [h["max_residual"] for h in iter_history]
        ax5_r = ax5.twinx()
        l1 = ax5.plot(iters, chi2_hist, "o-", color="C0", markersize=3,
                      label=r"$\chi^2$ in window")
        ax5.set_yscale("log")
        ax5.set_xlabel("iteration")
        ax5.set_ylabel(r"$\chi^2$ in closure window", color="C0")
        ax5.tick_params(axis="y", labelcolor="C0")
        ax5.grid(alpha=0.3, axis="x")

        l2 = ax5_r.plot(iters, max_res_hist, "s-", color="C3", markersize=3,
                        label=r"max $|r-1|$")
        ax5_r.set_ylabel(r"max $|D/M - 1|$ in window", color="C3")
        ax5_r.tick_params(axis="y", labelcolor="C3")
        ax5_r.set_ylim(0.0, max(0.1, max(max_res_hist) * 1.1))
        lines = l1 + l2
        labels = [l.get_label() for l in lines]
        ax5.legend(lines, labels, fontsize=8, loc="upper right")
        ax5.set_title(
            f"Iteration trajectory ({len(iters)} iters, "
            f"{'converged' if iter_converged else 'not converged'})",
            fontsize=10,
        )
    else:
        ax5.text(0.5, 0.5, "no iterative result", ha="center", va="center",
                 transform=ax5.transAxes)
        ax5.set_xticks([])
        ax5.set_yticks([])

    # =============================================================
    # sPHENIX Internal + crossing-angle label on every panel, plus
    # the per-period fit metadata in the suptitle.
    # =============================================================
    for _ax in (ax1, ax2, ax3, ax4, ax5):
        sphenix_label(_ax, period=period_name, loc="upper left", size=9.5)
    period_tag = sphenix_period_text(period_name)
    match_str = f" [photon10 mix]" if sample_match == "photon10" else (
        f" [match={sample_match!r}]" if sample_match else ""
    )
    fig.suptitle(
        rf"Truth vertex reweight — $p+p$ {period_tag}{match_str}   "
        rf"($f_{{\rm double}}={inputs.period.f_double:.3f}$,  "
        rf"closure $|z|<{z_cut:.0f}$ cm, {int(closure_mask.sum())} bins)",
        fontsize=11, y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    out_path = out_dir / "closure.pdf"
    fig.savefig(out_path)
    plt.close(fig)
    LOG.info("wrote %s", out_path)
    return 0


def main(argv=None):
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(here / "config.yaml"))
    ap.add_argument("--period", required=True, choices=("1p5mrad", "0mrad"))
    ap.add_argument("--match", default=None)
    ap.add_argument("--z-cut", type=float, default=None,
                    help="override z_cut (must match fit run)")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    with open(args.config) as f:
        config = yaml.safe_load(f)

    out_base = config["output"]["base_dir"]
    out_base_path = (here / out_base) if not Path(out_base).is_absolute() else Path(out_base)
    out_dir = Path(args.out_dir) if args.out_dir else (out_base_path / args.period)
    cache_dir = Path(args.cache_dir) if args.cache_dir else (out_base_path / "cache")
    return plot(out_dir, config, args.period, cache_dir,
                sample_match=args.match, z_cut_override=args.z_cut)


if __name__ == "__main__":
    sys.exit(main())
