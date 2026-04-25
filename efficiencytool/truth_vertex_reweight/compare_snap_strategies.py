#!/usr/bin/env python3
"""
Compare four post-iteration snap strategies on the truth-vertex reweight.

Strategies:
  (0) current     : snap outside |z|<z_cut to the seed Gaussian (production)
  (a) none        : no snap; use iteration's w(z) everywhere
  (b) wide        : snap only where data statistics truly vanish, i.e. at a
                    wider boundary z_wide derived from 99.9% of |z_r^data|
  (c) tail_fit    : fit log(w) = a + b·z² to the outer inner-window bins and
                    analytically extrapolate that Gaussian tail beyond z_cut

For each strategy and each period, the script computes the mixed-MC reco
histogram M and the closure residual D/M, and produces a comparison plot
plus a metrics table.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fit_core as fc
from fit_truth_vertex_reweight import (
    load_period_inputs, _data_reco_rms, _normalize_to_total,
    _gaussian_smooth_1bin,
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sphenix_style import set_sphenix_rcparams, sphenix_label, sphenix_period_text
set_sphenix_rcparams()

LOG = logging.getLogger("compare_snap_strategies")


def run_iteration_no_snap(inputs, sigma_gen, config) -> Tuple[np.ndarray, np.ndarray, List[dict], float]:
    """Re-run the iteration loop (identical to stage_iterative_fit) but STOP
    right before the post-iteration snap, returning the pre-snap w_raw."""
    iter_cfg = config["fit"].get("iterative", {})
    centers = inputs.binning.centers
    D = inputs.D
    f_s = inputs.period.f_single
    f_d = inputs.period.f_double
    z_cut = inputs.period.z_cut

    closure_mask = np.abs(centers) < z_cut
    update_mask = np.ones_like(centers, dtype=bool)

    sigma_seed = _data_reco_rms(D, centers)
    w_seed = fc.parametric_weight(centers, sigma_gen, sigma_seed)
    w = w_seed.copy()

    max_iter = int(iter_cfg.get("max_iterations", 30))
    alpha = float(iter_cfg.get("damping", 0.5))
    eps = float(iter_cfg.get("convergence_max_resid", 0.01))
    smooth_n = int(iter_cfg.get("residual_smooth_bins", 1))
    clip_lo, clip_hi = iter_cfg.get("update_clip", [0.5, 2.0])

    _all_zh_single = np.concatenate([
        s.z_h[s.z_h > -9000] for s in inputs.single_samples if s.z_h is not None
    ])

    history: List[dict] = []
    D_total = float(D.sum())

    for it in range(max_iter):
        w_func = fc.c_interpolator(centers, w)
        M = fc.mixed_mc_histogram(
            inputs.single_samples, inputs.double_samples,
            inputs.binning, f_s, f_d, w_func,
        )
        M = _normalize_to_total(M, D_total)

        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.ones_like(D, dtype=np.float64)
            valid = update_mask & (M > 0)
            r[valid] = D[valid] / M[valid]

        max_res = float(np.max(np.abs(r[closure_mask] - 1.0))) if closure_mask.any() else 0.0
        chi2_inside = fc.chi2(D, M, mask=closure_mask)
        history.append({"iteration": it, "max_residual": max_res, "chi2_in_window": chi2_inside})

        if max_res < eps:
            LOG.info("[%s] converged at iter %d", inputs.period.name, it)
            break

        r_smooth = r.copy()
        for _ in range(smooth_n):
            inner = r_smooth[update_mask]
            if len(inner) >= 3:
                r_smooth[update_mask] = _gaussian_smooth_1bin(inner)

        update = np.power(r_smooth[update_mask], alpha)
        update = np.clip(update, float(clip_lo), float(clip_hi))
        w[update_mask] *= update

        if int(update_mask.sum()) >= 3:
            idx_in = np.where(update_mask)[0]
            w[idx_in[0]] = w[idx_in[1]]
            w[idx_in[-1]] = w[idx_in[-2]]

        w_func_tmp = fc.c_interpolator(centers, w)
        mean_w = float(np.mean(w_func_tmp(_all_zh_single)))
        if mean_w > 1e-12:
            w /= mean_w

    return w, w_seed, history, sigma_seed


def strategy_current(w_raw, w_seed, centers, z_cut, **kw):
    w = w_raw.copy()
    mask_out = np.abs(centers) >= z_cut
    w[mask_out] = w_seed[mask_out]
    return w


def strategy_none(w_raw, w_seed, centers, z_cut, **kw):
    return w_raw.copy()


def strategy_wide(w_raw, w_seed, centers, z_cut, z_wide, **kw):
    w = w_raw.copy()
    mask_out = np.abs(centers) >= z_wide
    w[mask_out] = w_seed[mask_out]
    return w


def strategy_tail_fit(w_raw, w_seed, centers, z_cut, n_tail=6, **kw):
    """Fit log(w_raw) = a + b·z² over the last `n_tail` inner-window bins on
    each side of z_cut, then extrapolate that Gaussian functional form to
    bins outside the closure window."""
    w = w_raw.copy()
    closure_mask = np.abs(centers) < z_cut

    idx_inside = np.where(closure_mask)[0]
    if len(idx_inside) < 2 * n_tail + 1:
        LOG.warning("tail_fit: not enough inside bins; falling back to seed snap")
        return strategy_current(w_raw, w_seed, centers, z_cut)

    # Fit region: outer n_tail bins on each side of z_cut (still inside)
    fit_idx = np.concatenate([idx_inside[:n_tail], idx_inside[-n_tail:]])
    z_fit = centers[fit_idx]
    w_fit = w_raw[fit_idx]
    mask_pos = w_fit > 0
    if mask_pos.sum() < 4:
        LOG.warning("tail_fit: too few positive w_raw points in fit region; using seed")
        return strategy_current(w_raw, w_seed, centers, z_cut)

    y = np.log(w_fit[mask_pos])
    x = z_fit[mask_pos] ** 2
    # Linear fit y = a + b·x
    b, a = np.polyfit(x, y, 1)
    # Extrapolate outside closure window
    mask_out = ~closure_mask
    w[mask_out] = np.exp(a + b * centers[mask_out] ** 2)
    return w


def build_M(inputs, w, D_total):
    w_func = fc.c_interpolator(inputs.binning.centers, w)
    M = fc.mixed_mc_histogram(
        inputs.single_samples, inputs.double_samples, inputs.binning,
        inputs.period.f_single, inputs.period.f_double, w_func,
    )
    return _normalize_to_total(M, D_total)


def metrics(D, M, closure_mask):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(M > 0, D / M, 1.0)
    r_in = r[closure_mask]
    r_in = r_in[np.isfinite(r_in)]
    max_res = float(np.max(np.abs(r_in - 1.0))) if len(r_in) else float("nan")
    rms_res = float(np.sqrt(np.mean((r_in - 1.0) ** 2))) if len(r_in) else float("nan")
    chi2 = fc.chi2(D, M, mask=closure_mask)
    return r, max_res, rms_res, chi2


def _step(ax, centers, values, **kw):
    edges = np.zeros(len(centers) + 1)
    dz = centers[1] - centers[0]
    edges[:-1] = centers - dz / 2.0
    edges[-1] = centers[-1] + dz / 2.0
    ax.step(edges, np.concatenate([values, [values[-1]]]), where="post", **kw)


def _truth_rms_z_h(single_samples, double_samples, w, binning, f_s, f_d):
    """Return the RMS of z_h in the mixed-MC truth distribution under w(z)."""
    w_func = fc.c_interpolator(binning.centers, w)
    zh_all = []
    wt_all = []
    for s in single_samples:
        if s.z_h is None:
            continue
        m = s.z_h > -9000
        zh_all.append(s.z_h[m])
        wt_all.append(f_s * s.spec.event_weight * w_func(s.z_h[m]))
    for s in double_samples:
        if s.z_h is None:
            continue
        m = s.z_h > -9000
        if s.z_mb is not None:
            m &= (s.z_mb > -9000)
        zh = s.z_h[m]
        zh_all.append(zh)
        wdbl = f_d * s.spec.event_weight * w_func(zh)
        if s.z_mb is not None:
            wdbl *= w_func(s.z_mb[m])
        wt_all.append(wdbl)
    z = np.concatenate(zh_all)
    w_evt = np.concatenate(wt_all)
    w_sum = w_evt.sum()
    if w_sum <= 0:
        return float("nan")
    mean = float((z * w_evt).sum() / w_sum)
    var = float(((z - mean) ** 2 * w_evt).sum() / w_sum)
    return float(np.sqrt(max(var, 0.0)))


def compare_period(inputs, sigma_gen, config, z_wide) -> Dict:
    w_raw, w_seed, hist, sigma_seed = run_iteration_no_snap(inputs, sigma_gen, config)
    centers = inputs.binning.centers
    z_cut = inputs.period.z_cut
    closure_mask = np.abs(centers) < z_cut

    strategies = {
        "current (snap→seed)": strategy_current,
        "none (no snap)":       strategy_none,
        "wide (snap→seed @99.9%)": strategy_wide,
        "tail_fit (Gauss extrap)":  strategy_tail_fit,
    }
    results = {}
    D = inputs.D
    D_total = float(D.sum())
    for name, fn in strategies.items():
        w = fn(w_raw, w_seed, centers, z_cut, z_wide=z_wide)
        M = build_M(inputs, w, D_total)
        r, max_r, rms_r, chi2 = metrics(D, M, closure_mask)
        # Truth z_h RMS under this w (key for efficiency)
        rms_zh = _truth_rms_z_h(inputs.single_samples, inputs.double_samples,
                                w, inputs.binning, inputs.period.f_single,
                                inputs.period.f_double)
        # Fraction of MC events with |w| beyond 10× seed peak (sanity check)
        w_at_zero = np.interp(0.0, centers, w)
        results[name] = dict(w=w, M=M, r=r, max_r=max_r, rms_r=rms_r, chi2=chi2,
                             rms_zh=rms_zh, w0=w_at_zero)
    return {
        "period": inputs.period.name,
        "centers": centers,
        "D": D,
        "z_cut": z_cut,
        "z_wide": z_wide,
        "w_seed": w_seed,
        "w_raw": w_raw,
        "sigma_seed": sigma_seed,
        "results": results,
        "closure_mask": closure_mask,
    }


def plot_period(pd: Dict, out_path: Path):
    centers = pd["centers"]
    z_cut = pd["z_cut"]
    z_wide = pd["z_wide"]
    res = pd["results"]

    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.22, height_ratios=[1, 1])
    ax_w  = fig.add_subplot(gs[0, 0])
    ax_M  = fig.add_subplot(gs[0, 1])
    ax_r  = fig.add_subplot(gs[1, 0])
    ax_tx = fig.add_subplot(gs[1, 1]); ax_tx.axis("off")

    colors = {
        "current (snap→seed)": "C0",
        "none (no snap)":        "C3",
        "wide (snap→seed @99.9%)":  "C2",
        "tail_fit (Gauss extrap)":  "C4",
    }
    styles = {
        "current (snap→seed)": "--",
        "none (no snap)":        "-",
        "wide (snap→seed @99.9%)":  "-.",
        "tail_fit (Gauss extrap)":  ":",
    }

    # Panel: w(z) log
    _step(ax_w, centers, pd["w_seed"], color="0.5", lw=0.8, ls=":", label="seed Gauss (pre-iter)")
    for name, d in res.items():
        _step(ax_w, centers, d["w"], color=colors[name], lw=1.4, ls=styles[name],
              label=name)
    ax_w.axvspan(-z_cut, +z_cut, color="C1", alpha=0.10, zorder=0)
    ax_w.axvline(-z_cut, color="C1", lw=0.8, ls="-.", alpha=0.6)
    ax_w.axvline(+z_cut, color="C1", lw=0.8, ls="-.", alpha=0.6)
    if z_wide != z_cut:
        ax_w.axvline(-z_wide, color="C2", lw=0.8, ls=":", alpha=0.6)
        ax_w.axvline(+z_wide, color="C2", lw=0.8, ls=":", alpha=0.6)
    ax_w.set_yscale("log")
    ax_w.set_xlabel(r"$z$ (cm)")
    ax_w.set_ylabel(r"$w(z)$")
    ax_w.set_title("Reweight function per strategy")
    ax_w.legend(fontsize=8, loc="lower center")
    ax_w.grid(alpha=0.3, which="both")

    # Panel: mixed-MC reco vs data (log)
    _step(ax_M, centers, pd["D"], color="k", lw=1.5, label="data")
    for name, d in res.items():
        _step(ax_M, centers, d["M"], color=colors[name], lw=1.2, ls=styles[name],
              label=name)
    ax_M.axvspan(-z_cut, +z_cut, color="C1", alpha=0.10, zorder=0)
    ax_M.set_yscale("log")
    ax_M.set_xlabel(r"$z_{\rm reco}$ (cm)")
    ax_M.set_ylabel("events")
    ax_M.set_title("Mixed-MC reco per strategy")
    ax_M.legend(fontsize=8, loc="lower center")
    ax_M.grid(alpha=0.3, which="both")

    # Panel: residual D/M
    ax_r.fill_between([-200, 200], 0.98, 1.02, color="gray", alpha=0.2)
    ax_r.fill_between([-200, 200], 0.99, 1.01, color="gray", alpha=0.35)
    ax_r.axhline(1.0, color="k", lw=0.5, ls=":", alpha=0.5)
    ax_r.axvspan(-z_cut, +z_cut, color="C1", alpha=0.10, zorder=0)
    for name, d in res.items():
        _step(ax_r, centers, d["r"], color=colors[name], lw=1.2, ls=styles[name],
              label=name)
    ax_r.set_xlabel(r"$z_{\rm reco}$ (cm)")
    ax_r.set_ylabel(r"$D/M$")
    ax_r.set_xlim(-1.5 * z_cut, 1.5 * z_cut)
    # adaptive y range
    r_all = np.concatenate([
        d["r"][pd["closure_mask"]] for d in res.values()
    ])
    r_all = r_all[np.isfinite(r_all)]
    lo = max(0.5, float(np.min(r_all)) - 0.02)
    hi = min(2.0, float(np.max(r_all)) + 0.02)
    ax_r.set_ylim(lo, hi)
    ax_r.set_title("Closure D/M per strategy")
    ax_r.legend(fontsize=8, loc="lower center")
    ax_r.grid(alpha=0.3)

    # Text panel: metrics table
    lines = [f"Period: {pd['period']}   z_cut={z_cut:.0f} cm",
             f"z_wide (99.9%-tile): {z_wide:.0f} cm",
             f"σ_seed = {pd['sigma_seed']:.2f} cm",
             ""]
    hdr = f"{'strategy':<24s} {'χ²':>9s} {'max|r-1|':>10s} {'rms':>9s} {'w(0)':>7s} {'σ(zh)':>8s}"
    lines.append(hdr)
    lines.append("-" * len(hdr))
    for name, d in res.items():
        lines.append(
            f"{name:<24s} {d['chi2']:>9.1f} {d['max_r']:>10.4f} {d['rms_r']:>9.4f} "
            f"{d['w0']:>7.2f} {d['rms_zh']:>8.2f}"
        )
    ax_tx.text(0.0, 1.0, "\n".join(lines), transform=ax_tx.transAxes,
               va="top", ha="left", family="monospace", fontsize=9)

    for _a in (ax_w, ax_M, ax_r):
        sphenix_label(_a, period=pd["period"], loc="upper left", size=9.5)
    fig.suptitle(f"Snap-strategy comparison — $p+p$ {sphenix_period_text(pd['period'])}  (photon10 mix)",
                 fontsize=11, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    fig.savefig(out_path)
    plt.close(fig)
    LOG.info("wrote %s", out_path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(HERE / "config.yaml"))
    ap.add_argument("--periods", nargs="+", default=("1p5mrad", "0mrad"))
    ap.add_argument("--match", default="photon10")
    ap.add_argument("-v", "--verbose", action="store_true")
    ap.add_argument("--out-dir", default=str(HERE / "output"))
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    with open(args.config) as f:
        config = yaml.safe_load(f)

    out_base = Path(args.out_dir)
    cache_dir = out_base / "cache"

    summary_lines: List[str] = []
    for period in args.periods:
        LOG.info("=== period %s ===", period)
        inputs = load_period_inputs(
            config, period, cache_dir=cache_dir,
            sample_match=args.match, n_events=None, z_cut_override=None,
        )
        # σ_gen from single-MC truth z_h (same as production)
        sigma_gen, _ = fc.fit_sigma_gen(inputs.single_samples)
        # z_wide: 99.9 percentile of |z_r^data|
        from fit_truth_vertex_reweight import _filtered_data_z  # noqa
        # reload D non-period-filtered? We already have inputs.D; derive from centers.
        # simpler: compute from raw data sample available in fit
        # Use inputs.period to recover data; or fall back to a wider percentile
        # of D itself.
        z_wide = _data_percentile_from_hist(inputs.D, inputs.binning.centers, 99.9)
        LOG.info("  z_wide (99.9 %%-tile from data hist) = %.1f cm", z_wide)

        pd = compare_period(inputs, sigma_gen, config, z_wide=z_wide)

        out_path = out_base / f"snap_compare_{period}.pdf"
        plot_period(pd, out_path)
        dst = Path("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/figures") / f"snap_compare_{period}.pdf"
        import shutil
        shutil.copy2(out_path, dst)
        LOG.info("staged to %s", dst)

        summary_lines.append(f"\n=== {period} (z_cut={pd['z_cut']:.0f}, z_wide={z_wide:.0f}, σ_seed={pd['sigma_seed']:.2f}) ===")
        for name, d in pd["results"].items():
            summary_lines.append(
                f"  {name:<26s}  chi2={d['chi2']:>9.1f}  max|r-1|={d['max_r']:.4f}  "
                f"rms={d['rms_r']:.4f}  w(0)={d['w0']:.2f}  σ(zh)={d['rms_zh']:.2f}"
            )

    print("\n".join(summary_lines))


def _data_percentile_from_hist(D, centers, pct):
    """Return the |z| percentile that corresponds to cumulative fraction pct
    of D (two-sided, assuming symmetric shape around 0)."""
    D = np.asarray(D, dtype=float)
    total = D.sum()
    if total <= 0:
        return float(np.max(np.abs(centers)))
    # Sort bins by |z|, accumulate
    idx = np.argsort(np.abs(centers))
    cum = np.cumsum(D[idx]) / total * 100.0
    # find first index where cum >= pct
    hit = np.searchsorted(cum, pct)
    hit = min(hit, len(idx) - 1)
    return float(np.abs(centers[idx[hit]]))


if __name__ == "__main__":
    sys.exit(main())
