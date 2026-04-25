#!/usr/bin/env python3
"""
Cross-period and cross-method comparison plots for the truth vertex
reweight fit. Supplements per-period closure.pdf with:

  1. method_comparison_<period>.pdf
        Three-curve reco overlay (no reweight / seed Gaussian / iterative)
        + zoom panel + residual with error bars, for each period.
        Stage 1+2 (parametric + residual correction) is intentionally
        not shown: it is never used in the nominal pipeline.
  2. convergence_<period>.pdf
        χ² and max|r−1| trajectory, with annotations on converged state.
  3. pull_distribution_<period>.pdf
        Histogram of (D−M)/σ inside the closure window, plus cumulative.
  4. two_period_overlay.pdf
        Side-by-side 1.5 mrad vs 0 mrad for w(z), M(z_r), residual D/M.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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


LOG = logging.getLogger("plot_comparisons")


def _load_period(out_dir: Path, config: dict, period: str, cache_dir: Path,
                 sample_match: Optional[str] = None) -> Dict:
    yaml_path = out_dir / "fit_result.yaml"
    if not yaml_path.exists():
        raise FileNotFoundError(f"missing {yaml_path}")
    with open(yaml_path) as f:
        doc = yaml.safe_load(f)

    bnode = doc["binning"]
    centers = np.asarray(bnode["centers"])
    z_min, z_max = bnode["z_min"], bnode["z_max"]
    n_bins = bnode["n_bins"]
    binning = fc.Binning(z_min=z_min, z_max=z_max, n_bins=n_bins)

    s1 = doc["stage1_parametric"]
    s2 = doc["stage2_residual"]
    iter_doc = doc.get("iterative")

    from fit_truth_vertex_reweight import load_period_inputs
    inputs = load_period_inputs(
        config, period, cache_dir=cache_dir,
        sample_match=sample_match, n_events=None, z_cut_override=None,
    )
    D = inputs.D.astype(np.float64)
    f_s, f_d = inputs.period.f_single, inputs.period.f_double
    z_cut = float(inputs.period.z_cut)
    closure_mask = np.abs(centers) < z_cut

    sigma_gen = float(s1["sigma_gen"])
    sigma_tgt = float(s1["sigma_tgt"])

    def _w_unit(z):
        return np.ones_like(z, dtype=np.float64)

    def _w_param(z):
        return fc.parametric_weight(z, sigma_gen, sigma_tgt)

    c_func = fc.c_interpolator(centers, np.asarray(s2["c_smoothed"]))

    M_noRW = fc.mixed_mc_histogram(inputs.single_samples, inputs.double_samples,
                                    binning, f_s, f_d, _w_unit)
    M_noRW *= D.sum() / max(M_noRW.sum(), 1e-30)
    M_s1 = fc.mixed_mc_histogram(inputs.single_samples, inputs.double_samples,
                                  binning, f_s, f_d, _w_param)
    M_s1 *= D.sum() / max(M_s1.sum(), 1e-30)
    # Stage 1+2 is not shown in the nominal method comparison; skip the
    # forward-model call. Keep the key for downstream compatibility.
    M_s2 = np.zeros_like(M_s1)

    M_iter = None
    M_seed = None
    w_iter = None
    w_seed = None
    sigma_seed = None
    if iter_doc is not None:
        w_iter = np.asarray(iter_doc["w_iterative"])
        w_seed = np.asarray(iter_doc["w_seed"])
        sigma_seed = float(iter_doc["sigma_seed_tgt"])
        def _w_seed_fn(z, g=sigma_gen, s=sigma_seed):
            return fc.parametric_weight(z, g, s)
        M_seed = fc.mixed_mc_histogram(inputs.single_samples, inputs.double_samples,
                                        binning, f_s, f_d, _w_seed_fn)
        M_seed *= D.sum() / max(M_seed.sum(), 1e-30)
        w_iter_func = fc.c_interpolator(centers, w_iter)
        M_iter = fc.mixed_mc_histogram(inputs.single_samples, inputs.double_samples,
                                        binning, f_s, f_d, w_iter_func)
        M_iter *= D.sum() / max(M_iter.sum(), 1e-30)

    return {
        "period": period,
        "doc": doc,
        "centers": centers,
        "binning": binning,
        "closure_mask": closure_mask,
        "z_cut": z_cut,
        "z_min": z_min,
        "z_max": z_max,
        "inputs": inputs,
        "D": D,
        "M_noRW": M_noRW,
        "M_s1": M_s1,
        "M_s2": M_s2,
        "M_seed": M_seed,
        "M_iter": M_iter,
        "w_param": np.asarray(s1["w_values"]),
        "w_seed": w_seed,
        "w_iter": w_iter,
        "sigma_gen": sigma_gen,
        "sigma_tgt": sigma_tgt,
        "sigma_seed": sigma_seed,
    }


def _ratio_sigma(D: np.ndarray, M: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    r = np.full_like(D, np.nan, dtype=np.float64)
    s = np.full_like(D, np.nan, dtype=np.float64)
    valid = (D > 0) & (M > 0)
    r[valid] = D[valid] / M[valid]
    s[valid] = r[valid] * np.sqrt(1.0 / D[valid] + 1.0 / M[valid])
    return r, s


def _step_hist(ax, centers, values, **kw):
    edges = np.zeros(len(centers) + 1)
    if len(centers) > 1:
        dz = centers[1] - centers[0]
    else:
        dz = 1.0
    edges[:-1] = centers - dz / 2.0
    edges[-1] = centers[-1] + dz / 2.0
    ax.step(edges, np.concatenate([values, [values[-1]]]), where="post", **kw)


# ---------------------------------------------------------------------------
# Plot 1: method comparison for one period
# ---------------------------------------------------------------------------

def plot_method_comparison(pd: Dict, out_path: Path) -> None:
    centers = pd["centers"]
    D = pd["D"]
    z_cut = pd["z_cut"]
    closure_mask = pd["closure_mask"]
    z_min = pd["z_min"]
    z_max = pd["z_max"]

    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.22, height_ratios=[2, 1])
    ax_log = fig.add_subplot(gs[0, 0])
    ax_lin = fig.add_subplot(gs[0, 1])
    ax_res = fig.add_subplot(gs[1, :])

    sigma_seed = pd.get("sigma_seed")
    sigma_tgt = pd.get("sigma_tgt")
    param_label = (
        fr"MC (parametric Gaussian, $\sigma_{{\rm tgt}}={sigma_tgt:.1f}$)"
        if sigma_tgt is not None else "MC (parametric Gaussian)"
    )
    seed_label = (
        fr"MC (seed Gaussian, $\sigma_{{\rm seed}}={sigma_seed:.1f}$)"
        if sigma_seed is not None else "MC (seed Gaussian)"
    )
    methods = [
        ("MC (no reweight)", pd["M_noRW"], "C7", ":", 1.0),
    ]
    if pd.get("M_seed") is not None:
        methods.append((seed_label, pd["M_seed"], "C1", "-.", 1.0))
    methods.append((param_label, pd["M_s1"], "C0", "--", 1.0))
    if pd["M_iter"] is not None:
        methods.append(("MC (iterative)", pd["M_iter"], "C3", "-", 1.8))

    # Top-left: log
    _step_hist(ax_log, centers, D, color="k", lw=1.5, label="data")
    for name, M, color, ls, lw in methods:
        _step_hist(ax_log, centers, M, color=color, lw=lw, ls=ls, label=name)
    ax_log.axvspan(-z_cut, +z_cut, color="C2", alpha=0.08)
    ax_log.set_yscale("log")
    ax_log.set_xlabel(r"$z_{\rm reco}$ (cm)")
    ax_log.set_ylabel("events")
    ax_log.set_title(f"{pd['period']}: reco vertex (log)")
    ax_log.legend(fontsize=8, loc="lower center")
    ax_log.grid(alpha=0.3)
    ax_log.set_xlim(z_min, z_max)

    # Top-right: linear, zoomed
    _step_hist(ax_lin, centers, D, color="k", lw=1.5, label="data")
    for name, M, color, ls, lw in methods:
        if "no reweight" in name:
            continue
        _step_hist(ax_lin, centers, M, color=color, lw=lw, ls=ls, label=name)
    ax_lin.axvspan(-z_cut, +z_cut, color="C2", alpha=0.08)
    ax_lin.set_xlabel(r"$z_{\rm reco}$ (cm)")
    ax_lin.set_ylabel("events")
    ax_lin.set_title(f"{pd['period']}: reco vertex (linear, zoom)")
    ax_lin.legend(fontsize=8, loc="upper right")
    ax_lin.grid(alpha=0.3)
    zoom = min(z_max, max(1.5 * z_cut, 1.0))
    ax_lin.set_xlim(-zoom, zoom)

    # Bottom: residual ratios with error bars
    ax_res.fill_between([z_min, z_max], 0.98, 1.02, color="gray", alpha=0.2, label="±2%")
    ax_res.fill_between([z_min, z_max], 0.99, 1.01, color="gray", alpha=0.35, label="±1%")
    ax_res.axhline(1.0, color="k", lw=0.5, ls=":", alpha=0.5)
    ax_res.axvspan(-z_cut, +z_cut, color="C2", alpha=0.08)

    for name, M, color, ls, lw in methods:
        if "no reweight" in name:
            continue
        r, e = _ratio_sigma(D, M)
        if "iterative" in name:
            valid = np.isfinite(r) & np.isfinite(e)
            ax_res.errorbar(centers[valid], r[valid], yerr=e[valid],
                            fmt="o", color=color, markersize=3, elinewidth=1.0,
                            capsize=0, label=name)
        else:
            _step_hist(ax_res, centers, r, color=color, lw=lw, ls=ls, label=name)

    ax_res.set_xlabel(r"$z_{\rm reco}$ (cm)")
    ax_res.set_ylabel(r"$D / M$")
    ax_res.set_xlim(-1.5 * z_cut, 1.5 * z_cut)
    # Adaptive y range
    if pd["M_iter"] is not None:
        r_it, _ = _ratio_sigma(D, pd["M_iter"])
        r_in = r_it[closure_mask]
        r_in = r_in[np.isfinite(r_in)]
        if len(r_in) > 0:
            lo = min(0.9, float(np.min(r_in)) - 0.02)
            hi = max(1.1, float(np.max(r_in)) + 0.02)
            ax_res.set_ylim(lo, hi)
    ax_res.legend(fontsize=8, loc="lower center", ncol=3)
    ax_res.grid(alpha=0.3)
    ax_res.set_title(f"{pd['period']}: closure D/M (in window ±{z_cut:.0f} cm)",
                      fontsize=10)

    for _ax in (ax_log, ax_lin, ax_res):
        sphenix_label(_ax, period=pd["period"], loc="upper left", size=9.5)
    fig.suptitle(
        f"Method comparison — $p+p$ {sphenix_period_text(pd['period'])} "
        f"(photon10 mix)",
        fontsize=11, y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(out_path)
    plt.close(fig)
    LOG.info("wrote %s", out_path)


# ---------------------------------------------------------------------------
# Plot 2: convergence trajectory
# ---------------------------------------------------------------------------

def plot_convergence(pd: Dict, out_path: Path) -> None:
    iter_doc = pd["doc"].get("iterative")
    if iter_doc is None or not iter_doc.get("history"):
        return
    history = iter_doc["history"]
    iters = np.array([h["iteration"] for h in history])
    chi2 = np.array([h["chi2_in_window"] for h in history])
    max_r = np.array([h["max_residual"] for h in history])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    axes[0].plot(iters, chi2, "o-", color="C0", markersize=4)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("iteration")
    axes[0].set_ylabel(r"$\chi^2$ in closure window")
    axes[0].set_title(fr"{sphenix_period_text(pd['period'])}: $\chi^2$ trajectory")
    axes[0].grid(alpha=0.3, which="both")
    sphenix_label(axes[0], period=pd["period"], loc="upper right", size=9.0)
    sphenix_label(axes[1], period=pd["period"], loc="upper right", size=9.0)
    # mark the convergence state
    axes[0].annotate(
        f"final = {chi2[-1]:.1f}\nreduction = {chi2[0]/chi2[-1]:.1f}×",
        xy=(iters[-1], chi2[-1]), xycoords="data",
        xytext=(-60, 30), textcoords="offset points",
        fontsize=9, color="C0",
        arrowprops=dict(arrowstyle="->", color="C0", alpha=0.5),
    )

    axes[1].plot(iters, max_r, "s-", color="C3", markersize=4)
    axes[1].axhline(0.01, color="k", lw=0.5, ls=":", alpha=0.5, label="1% target")
    axes[1].set_xlabel("iteration")
    axes[1].set_ylabel(r"$\max |D/M - 1|$ in window")
    axes[1].set_title(f"{sphenix_period_text(pd['period'])}: convergence")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    axes[1].annotate(
        f"final = {max_r[-1]:.3f}\n"
        f"{'converged' if iter_doc['converged'] else 'not converged'}",
        xy=(iters[-1], max_r[-1]), xycoords="data",
        xytext=(-80, 20), textcoords="offset points",
        fontsize=9, color="C3",
        arrowprops=dict(arrowstyle="->", color="C3", alpha=0.5),
    )

    fig.suptitle(
        fr"Iterative convergence trajectory — $p+p$ {sphenix_period_text(pd['period'])}",
        fontsize=11, y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path)
    plt.close(fig)
    LOG.info("wrote %s", out_path)


# ---------------------------------------------------------------------------
# Plot 3: pull distribution
# ---------------------------------------------------------------------------

def plot_pull(pd: Dict, out_path: Path) -> None:
    D = pd["D"]
    M_iter = pd["M_iter"]
    if M_iter is None:
        return
    closure_mask = pd["closure_mask"]
    with np.errstate(divide="ignore", invalid="ignore"):
        err2 = np.where(D > 0, D, 1.0) + np.where(M_iter > 0, M_iter, 1.0)
        pull = (D - M_iter) / np.sqrt(err2)
    pull_in = pull[closure_mask]
    pull_in = pull_in[np.isfinite(pull_in)]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left: pull histogram
    axes[0].hist(pull_in, bins=max(10, len(pull_in) // 2), color="C3",
                  edgecolor="k", alpha=0.6)
    axes[0].axvline(0, color="k", lw=0.5, ls="--")
    axes[0].set_xlabel(r"pull = $(D - M) / \sqrt{D+M}$")
    axes[0].set_ylabel("bins")
    mu = float(np.mean(pull_in))
    sigma = float(np.std(pull_in, ddof=1)) if len(pull_in) > 1 else float("nan")
    axes[0].set_title(
        f"{pd['period']}: pull distribution "
        f"($N={len(pull_in)}$, $\\mu={mu:+.2f}$, $\\sigma={sigma:.2f}$)",
        fontsize=10,
    )
    axes[0].grid(alpha=0.3)
    # Reference N(0,1)
    if len(pull_in) > 0:
        x = np.linspace(min(-4, pull_in.min() - 0.5),
                         max(4, pull_in.max() + 0.5), 200)
        from math import pi, sqrt
        # Expected N(0,1) scaled to bin count
        n_bins_plot = max(10, len(pull_in) // 2)
        x_range = float(max(x) - min(x))
        n_expected = len(pull_in) / n_bins_plot
        axes[0].plot(
            x,
            n_expected * x_range / n_bins_plot * (1 / (sqrt(2 * pi))) * np.exp(-0.5 * x ** 2),
            color="C0", lw=1.2, ls="--", alpha=0.7, label=r"$\mathcal{N}(0,1)$",
        )
        axes[0].legend(fontsize=8)

    # Right: per-bin pull vs z
    centers = pd["centers"]
    mask = closure_mask
    axes[1].fill_between([-pd["z_cut"], pd["z_cut"]], -1, 1,
                          color="gray", alpha=0.2, label="±1")
    axes[1].fill_between([-pd["z_cut"], pd["z_cut"]], -2, 2,
                          color="gray", alpha=0.1, label="±2")
    axes[1].axhline(0, color="k", lw=0.5, ls="--")
    axes[1].plot(centers[mask], pull[mask], "o", color="C3", markersize=4)
    axes[1].set_xlabel(r"$z_{\rm reco}$ (cm)")
    axes[1].set_ylabel("pull per bin")
    axes[1].set_title(f"{pd['period']}: per-bin pull in closure window", fontsize=10)
    axes[1].legend(fontsize=8, loc="upper right")
    axes[1].grid(alpha=0.3)

    for _ax in axes:
        sphenix_label(_ax, period=pd["period"], loc="upper left", size=9.0)
    fig.suptitle(
        f"Pull diagnostics — $p+p$ {sphenix_period_text(pd['period'])}",
        fontsize=11, y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path)
    plt.close(fig)
    LOG.info("wrote %s", out_path)


# ---------------------------------------------------------------------------
# Plot 4: two-period overlay
# ---------------------------------------------------------------------------

def plot_two_period_overlay(pds: List[Dict], out_path: Path) -> None:
    if len(pds) < 2:
        return
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # Top-left: reco MC iterative vs data, both periods
    for pd, color in zip(pds, ["C0", "C3"]):
        _step_hist(axes[0, 0], pd["centers"], pd["D"], color=color, lw=1.5,
                    ls="-", label=f"{pd['period']} data")
        if pd["M_iter"] is not None:
            _step_hist(axes[0, 0], pd["centers"], pd["M_iter"], color=color,
                        lw=1.0, ls="--", label=f"{pd['period']} MC iter")
        axes[0, 0].axvline(-pd["z_cut"], color=color, lw=0.5, ls=":", alpha=0.5)
        axes[0, 0].axvline(+pd["z_cut"], color=color, lw=0.5, ls=":", alpha=0.5)
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_xlabel(r"$z_{\rm reco}$ (cm)")
    axes[0, 0].set_ylabel("events")
    axes[0, 0].set_title("Reco vertex: data vs iterative MC (both periods)")
    axes[0, 0].legend(fontsize=8, loc="lower center")
    axes[0, 0].grid(alpha=0.3)

    # Top-right: residual D/M iterative, both periods
    for pd, color in zip(pds, ["C0", "C3"]):
        if pd["M_iter"] is None:
            continue
        r, e = _ratio_sigma(pd["D"], pd["M_iter"])
        valid = np.isfinite(r) & np.isfinite(e)
        axes[0, 1].errorbar(pd["centers"][valid], r[valid], yerr=e[valid],
                             fmt="o", markersize=3, color=color,
                             elinewidth=1.0, capsize=0,
                             label=f"{pd['period']}")
        axes[0, 1].axvline(-pd["z_cut"], color=color, lw=0.5, ls=":", alpha=0.5)
        axes[0, 1].axvline(+pd["z_cut"], color=color, lw=0.5, ls=":", alpha=0.5)
    axes[0, 1].fill_between([-200, 200], 0.98, 1.02, color="gray", alpha=0.2)
    axes[0, 1].fill_between([-200, 200], 0.99, 1.01, color="gray", alpha=0.35)
    axes[0, 1].axhline(1.0, color="k", lw=0.5, ls=":", alpha=0.5)
    axes[0, 1].set_xlabel(r"$z_{\rm reco}$ (cm)")
    axes[0, 1].set_ylabel(r"$D / M$")
    axes[0, 1].set_title("Iterative residual in closure window")
    axes[0, 1].set_xlim(-200, 200)
    axes[0, 1].set_ylim(0.9, 1.1)
    axes[0, 1].legend(fontsize=8, loc="lower center")
    axes[0, 1].grid(alpha=0.3)

    # Bottom-left: w_iterative, both periods (log)
    for pd, color in zip(pds, ["C0", "C3"]):
        if pd["w_iter"] is not None:
            _step_hist(axes[1, 0], pd["centers"], pd["w_iter"], color=color,
                        lw=1.5, label=f"{pd['period']}")
        axes[1, 0].axvline(-pd["z_cut"], color=color, lw=0.5, ls=":", alpha=0.5)
        axes[1, 0].axvline(+pd["z_cut"], color=color, lw=0.5, ls=":", alpha=0.5)
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_xlabel(r"$z$ (cm)")
    axes[1, 0].set_ylabel(r"$w_{\rm iter}(z)$")
    axes[1, 0].set_title("Iterative reweight function, both periods")
    axes[1, 0].legend(fontsize=8, loc="upper right")
    axes[1, 0].grid(alpha=0.3, which="both")

    # Bottom-right: summary table as text
    axes[1, 1].axis("off")
    lines = ["Period summary:"]
    for pd in pds:
        it = pd["doc"].get("iterative")
        if it:
            lines.append(
                f"  {pd['period']:8s}: z_cut={pd['z_cut']:.1f} cm, "
                f"σ_seed={it['sigma_seed_tgt']:.1f} cm"
            )
            lines.append(
                f"    f_double={pd['inputs'].period.f_double:.3f}, "
                f"iter={it['n_iterations']}"
            )
            lines.append(
                f"    max|r−1|={it['max_residual_in_window']:.3f}, "
                f"χ²/ndf={it['chi2_in_window']/max(it['ndf_in_window'],1):.2f}"
            )
            lines.append(
                f"    converged={bool(it['converged'])}"
            )
            lines.append("")
    axes[1, 1].text(0.02, 0.98, "\n".join(lines), transform=axes[1, 1].transAxes,
                     fontsize=10, va="top", ha="left", family="monospace")

    for _ax in (axes[0, 0], axes[0, 1], axes[1, 0]):
        sphenix_label(_ax, period=None, loc="upper left", size=9.5)
    fig.suptitle("Two-period overlay — $p+p$ (photon10 mix)",
                 fontsize=11, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(out_path)
    plt.close(fig)
    LOG.info("wrote %s", out_path)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main(argv=None):
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(here / "config.yaml"))
    ap.add_argument("--periods", nargs="+", default=("1p5mrad", "0mrad"))
    ap.add_argument("--match", default=None)
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
    cache_dir = Path(args.cache_dir) if args.cache_dir else (out_base_path / "cache")

    pds: List[Dict] = []
    for period in args.periods:
        period_dir = Path(args.out_dir) / period if args.out_dir else (out_base_path / period)
        try:
            pd = _load_period(period_dir, config, period, cache_dir, sample_match=args.match)
            pds.append(pd)
            plot_method_comparison(pd, period_dir / "method_comparison.pdf")
            plot_convergence(pd, period_dir / "convergence.pdf")
            plot_pull(pd, period_dir / "pull_distribution.pdf")
        except Exception as e:
            LOG.error("period %s: %s", period, e)

    if len(pds) >= 2:
        two_period_out = out_base_path / "two_period_overlay.pdf"
        plot_two_period_overlay(pds, two_period_out)

    return 0


if __name__ == "__main__":
    sys.exit(main())
