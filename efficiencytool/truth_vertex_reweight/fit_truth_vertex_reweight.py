#!/usr/bin/env python3
"""
Two-stage truth vertex reweight fit (event-level).

Stage 1 — Parametric Gaussian.  For each trial σ_tgt, reweight every MC
event by the Gaussian ratio w(z) evaluated at its truth vertex (once
for single MC, twice for double MC), histogram z_reco, mix single+double
at the physical fractions, and minimize χ² against data.

Stage 2 — Residual correction.  Compute c(z_r) = D/M_stage1 and apply
it as an extra per-vertex factor alongside w(z). The stage-2 integral-
preserving scale is solved analytically (quadratic in α because c is
applied twice for double events).

Outputs (under output/<period>/):

  fit_result.yaml  — σ_gen, σ_tgt, χ² per stage, residuals, w / c arrays
  reweight.root    — TH1 histograms + TNamed metadata
  fit.log          — textual log
"""
from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np
import yaml

import fit_core as fc

try:
    import ROOT                                    # noqa: F401
    HAS_ROOT = True
except ImportError:
    HAS_ROOT = False


LOG = logging.getLogger("fit_truth_vertex_reweight")


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class PeriodConfig:
    name: str
    run_min: int
    run_max: int
    lumi: float
    f_double: float
    z_cut: float = 60.0
    f_single: float = field(init=False)

    def __post_init__(self):
        self.f_single = 1.0 - self.f_double


@dataclass
class FitInputs:
    binning: fc.Binning
    period: PeriodConfig
    single_samples: List[fc.LoadedSample]
    double_samples: List[fc.LoadedSample]
    data_sample: fc.LoadedSample
    D: np.ndarray                                  # data reco histogram


def load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def build_sample_spec(name: str, node: dict, kind: str) -> fc.SampleSpec:
    return fc.SampleSpec(
        name=name,
        path=node["path"],
        xs=float(node.get("xs", 1.0)),
        n_gen=float(node.get("n_gen", 1.0)),
        kind=kind,
    )


def _sample_match(name: str, pattern: Optional[str]) -> bool:
    if pattern is None or pattern == "":
        return True
    return pattern in name


def _filtered_data_z(data_sample: fc.LoadedSample, run_min: int, run_max: int) -> np.ndarray:
    """Return the data reco vertex array for events inside the period and
    with valid z_r (excluding the −9999 sentinel). Also applies the
    trigger-passed mask if attached, so percentile diagnostics match what
    the fit actually uses downstream."""
    z = data_sample.z_r
    mask = np.isfinite(z) & (z > fc.SENTINEL_F + 1)
    if data_sample.runnumber is not None:
        mask = mask & (data_sample.runnumber >= run_min) & (data_sample.runnumber < run_max)
    if data_sample.trigger_passed is not None:
        mask = mask & data_sample.trigger_passed.astype(bool)
    return z[mask]


def _resolve_z_cut(p_node: dict, z_data_period: np.ndarray) -> float:
    """Return z_cut: from z_cut_percentile (preferred) or explicit z_cut."""
    pct = p_node.get("z_cut_percentile")
    if pct is not None and len(z_data_period) > 0:
        z_cut = float(np.percentile(np.abs(z_data_period), float(pct)))
        z_cut = float(np.ceil(z_cut))
        LOG.info("z_cut derived from %.1f-th percentile of |z_data_reco|: %.1f cm",
                 float(pct), z_cut)
        return z_cut
    z_cut = float(p_node.get("z_cut", 60.0))
    LOG.info("z_cut explicit: %.1f cm", z_cut)
    return z_cut


def load_period_inputs(
    config: dict,
    period_name: str,
    cache_dir: Path,
    sample_match: Optional[str] = None,
    n_events: Optional[int] = None,
    z_cut_override: Optional[float] = None,
) -> FitInputs:
    b = config["branches"]
    binning = fc.Binning(**config["binning"])
    p_node = config["periods"][period_name]
    run_min = int(p_node["run_min"])
    run_max = int(p_node["run_max"])

    # 1. Load data first (cached) — needed for z_cut percentile derivation.
    data_node = config["inputs"]["data"]
    data_path = data_node["path"] if isinstance(data_node, dict) else data_node
    if "TODO" in str(data_path):
        raise FileNotFoundError("inputs.data.path is a placeholder")
    data_spec = fc.SampleSpec(name="data", path=data_path, xs=1.0, n_gen=1.0, kind="data")
    LOG.info("loading data from %s", data_path)
    # Pull trigger filter from config (mirrors RecoEffCalculator's `trigger_used`).
    # If null/absent, no trigger filter is applied (legacy behaviour).
    trig_node = config.get("data_trigger_bit")
    trigger_bit = int(trig_node) if trig_node is not None else None
    if trigger_bit is not None:
        LOG.info("data trigger filter: scaledtrigger[%d] != 0", trigger_bit)
    data_sample = fc.load_sample_cached(data_spec, b, cache_dir, n_events=n_events,
                                        trigger_bit=trigger_bit)

    # 2. Period-filter and percentile diagnostics.
    z_data_period = _filtered_data_z(data_sample, run_min, run_max)
    if len(z_data_period) == 0:
        raise RuntimeError(f"no data events in period {period_name} (run [{run_min}, {run_max}))")
    LOG.info("data in period %s: %d events, |z| percentiles "
             "50=%.1f / 90=%.1f / 95=%.1f / 99=%.1f / 99.5=%.1f cm",
             period_name, len(z_data_period),
             float(np.percentile(np.abs(z_data_period), 50)),
             float(np.percentile(np.abs(z_data_period), 90)),
             float(np.percentile(np.abs(z_data_period), 95)),
             float(np.percentile(np.abs(z_data_period), 99)),
             float(np.percentile(np.abs(z_data_period), 99.5)))

    # 3. Derive z_cut (closure window) from data percentile or fallback.
    if z_cut_override is not None:
        z_cut = float(z_cut_override)
        LOG.info("z_cut override from CLI: %.1f cm", z_cut)
    else:
        z_cut = _resolve_z_cut(p_node, z_data_period)
    if z_cut > binning.z_max:
        LOG.warning("z_cut (%.1f) > binning z_max (%.1f); clipping z_cut to z_max",
                    z_cut, binning.z_max)
        z_cut = binning.z_max

    period = PeriodConfig(
        name=period_name,
        run_min=run_min,
        run_max=run_max,
        lumi=float(p_node["lumi"]),
        f_double=float(p_node["f_double"]),
        z_cut=z_cut,
    )
    LOG.info("period %s: runs [%d, %d), lumi=%.3f pb^-1, f_double=%.4f, z_cut=%.1f cm",
             period.name, period.run_min, period.run_max, period.lumi, period.f_double, period.z_cut)
    LOG.info("binning (derivation window): z=[%.1f, %.1f] cm, %d bins (%.2f cm/bin)",
             binning.z_min, binning.z_max, binning.n_bins,
             (binning.z_max - binning.z_min) / binning.n_bins)

    # 4. Load MC samples (cached, no dependence on binning).
    single_samples: List[fc.LoadedSample] = []
    for name, node in config["inputs"]["single_mc"].items():
        if not _sample_match(name, sample_match):
            LOG.info("skipping single MC %s (not matching filter %r)", name, sample_match)
            continue
        if "TODO" in str(node["path"]):
            raise FileNotFoundError(f"single_mc.{name}.path is a placeholder")
        spec = build_sample_spec(name, node, kind="single_mc")
        LOG.info("loading single MC %s (xs=%.3e, n_gen=%.3e)", name, spec.xs, spec.n_gen)
        single_samples.append(fc.load_sample_cached(spec, b, cache_dir, n_events=n_events))

    double_samples: List[fc.LoadedSample] = []
    for name, node in config["inputs"]["double_mc"].items():
        if not _sample_match(name, sample_match):
            LOG.info("skipping double MC %s (not matching filter %r)", name, sample_match)
            continue
        if "TODO" in str(node["path"]):
            raise FileNotFoundError(f"double_mc.{name}.path is a placeholder")
        spec = build_sample_spec(name, node, kind="double_mc")
        LOG.info("loading double MC %s (xs=%.3e, n_gen=%.3e)", name, spec.xs, spec.n_gen)
        double_samples.append(
            fc.load_sample_cached(spec, b, cache_dir, n_events=n_events, require_mb=True)
        )

    if not single_samples or not double_samples:
        raise RuntimeError("filter matched 0 single or double samples; cannot fit")

    # 5. Build the data histogram with the binning.
    LOG.info("building D[%d] for period %s ...", binning.n_bins, period.name)
    D = fc.data_reco_histogram(
        data_sample, binning, period=(run_min, run_max)
    )
    LOG.info("  D sum = %.3e (entries: %d)", D.sum(), int(D.sum()))

    return FitInputs(
        binning=binning, period=period,
        single_samples=single_samples, double_samples=double_samples,
        data_sample=data_sample, D=D,
    )


# ---------------------------------------------------------------------------
# Stage 1: parametric Gaussian scan
# ---------------------------------------------------------------------------

@dataclass
class Stage1Result:
    sigma_gen: float
    sigma_gen_err: float
    sigma_tgt: float
    chi2: float
    ndf: int
    w_parametric: np.ndarray              # w(centers) for saving
    grid_sigma: np.ndarray
    grid_chi2: np.ndarray
    M_stage1: np.ndarray                  # mixed MC at the best σ_tgt, rescaled to D.sum


def _normalize_to_total(h: np.ndarray, total: float) -> np.ndarray:
    s = h.sum()
    if s <= 0:
        return h
    return h * (total / s)


def stage1_fit(inputs: FitInputs, config: dict) -> Stage1Result:
    fit_cfg = config["fit"]

    if fit_cfg.get("sigma_gen") is not None:
        sigma_gen = float(fit_cfg["sigma_gen"])
        sigma_gen_err = 0.0
        LOG.info("using hardcoded sigma_gen = %.3f cm", sigma_gen)
    else:
        sigma_gen, sigma_gen_err = fc.fit_sigma_gen(inputs.single_samples)
        LOG.info("measured sigma_gen = %.3f ± %.3f cm (from nominal single MC)",
                 sigma_gen, sigma_gen_err)

    sr = fit_cfg["sigma_tgt_range"]
    ss = float(fit_cfg["sigma_tgt_step"])
    grid = np.arange(sr[0], sr[1] + 0.5 * ss, ss)
    LOG.info("stage 1 grid: sigma_tgt in [%.1f, %.1f] cm step %.2f (%d trials)",
             sr[0], sr[1], ss, len(grid))

    f_s = inputs.period.f_single
    f_d = inputs.period.f_double
    D = inputs.D
    D_total = float(D.sum())

    # χ² is restricted to bins inside the period's |z| < z_cut analysis window
    fit_mask = np.abs(inputs.binning.centers) < inputs.period.z_cut
    LOG.info("stage 1 χ² mask: %d / %d bins inside |z|<%.1f cm",
             int(fit_mask.sum()), len(fit_mask), inputs.period.z_cut)

    chi2_values = np.empty(len(grid))
    for i, sig in enumerate(grid):
        w_func = (lambda z, g=sigma_gen, s=sig: fc.parametric_weight(z, g, s))
        M = fc.mixed_mc_histogram(
            inputs.single_samples, inputs.double_samples,
            inputs.binning, f_s, f_d, w_func,
        )
        M = _normalize_to_total(M, D_total)
        chi2_values[i] = fc.chi2(D, M, mask=fit_mask)

    best_idx = int(np.argmin(chi2_values))
    best_sigma = float(grid[best_idx])
    best_chi2 = float(chi2_values[best_idx])

    # Parabolic refinement
    if 0 < best_idx < len(grid) - 1:
        x = grid[best_idx - 1: best_idx + 2]
        y = chi2_values[best_idx - 1: best_idx + 2]
        denom = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
        if abs(denom) > 1e-9:
            a = (x[2] * (y[1] - y[0]) + x[1] * (y[0] - y[2]) + x[0] * (y[2] - y[1])) / denom
            b = (x[2] ** 2 * (y[0] - y[1]) + x[1] ** 2 * (y[2] - y[0]) + x[0] ** 2 * (y[1] - y[2])) / denom
            if abs(a) > 1e-12:
                x_min = -b / (2 * a)
                if sr[0] <= x_min <= sr[1]:
                    best_sigma = float(x_min)

    w_best = (lambda z, g=sigma_gen, s=best_sigma: fc.parametric_weight(z, g, s))
    M_best = fc.mixed_mc_histogram(
        inputs.single_samples, inputs.double_samples,
        inputs.binning, f_s, f_d, w_best,
    )
    M_best = _normalize_to_total(M_best, D_total)

    ndf = int(fit_mask.sum()) - 1
    LOG.info("stage 1 best sigma_tgt = %.3f cm, chi2 = %.3f, ndf = %d (chi2/ndf = %.3f)",
             best_sigma, best_chi2, ndf, best_chi2 / max(ndf, 1))

    # Per-bin w(z) snapshot for output
    w_centers = fc.parametric_weight(inputs.binning.centers, sigma_gen, best_sigma)

    return Stage1Result(
        sigma_gen=sigma_gen, sigma_gen_err=sigma_gen_err,
        sigma_tgt=best_sigma, chi2=best_chi2, ndf=ndf,
        w_parametric=w_centers,
        grid_sigma=grid, grid_chi2=chi2_values,
        M_stage1=M_best,
    )


# ---------------------------------------------------------------------------
# Stage 2: residual correction
# ---------------------------------------------------------------------------

@dataclass
class Stage2Result:
    c_raw: np.ndarray
    c_smoothed: np.ndarray
    w_total: np.ndarray                   # w_parametric * c at bin centers
    M_stage2: np.ndarray
    residual_stage2: np.ndarray           # D / M_stage2
    max_residual: float
    alpha: float                          # integral-preserving rescale on c


def _gaussian_smooth_1bin(x: np.ndarray) -> np.ndarray:
    if len(x) < 3:
        return x.copy()
    pad = np.concatenate([[x[0]], x, [x[-1]]])
    return 0.25 * pad[:-2] + 0.5 * pad[1:-1] + 0.25 * pad[2:]


def stage2_fit(inputs: FitInputs, stage1: Stage1Result, config: dict) -> Stage2Result:
    fit_cfg = config["fit"]

    D = inputs.D
    M1 = stage1.M_stage1

    # Stage 2 residual is only defined inside the analysis window.
    # Outside |z|<z_cut we set c(z) = 1 (no correction) — those bins
    # are dropped from the χ² and the analysis cuts them anyway.
    fit_mask = np.abs(inputs.binning.centers) < inputs.period.z_cut

    with np.errstate(divide="ignore", invalid="ignore"):
        c_raw = np.ones_like(D, dtype=np.float64)
        valid = fit_mask & (M1 > 0)
        c_raw[valid] = D[valid] / M1[valid]
    in_window_max = float(np.max(np.abs(c_raw[fit_mask] - 1.0))) if fit_mask.any() else 0.0
    LOG.info("stage 2 raw residual (inside |z|<%.1f): min=%.4f max=%.4f max|c-1|=%.4f",
             inputs.period.z_cut, float(c_raw[fit_mask].min()),
             float(c_raw[fit_mask].max()), in_window_max)

    c = c_raw.copy()
    n_smooth = int(fit_cfg.get("residual_smooth_bins", 1))
    for _ in range(n_smooth):
        # smooth only the in-window slice, keep c=1 outside
        in_slice = c[fit_mask]
        if len(in_slice) >= 3:
            c[fit_mask] = _gaussian_smooth_1bin(in_slice)

    clip_lo, clip_hi = fit_cfg.get("residual_clip", [0.1, 10.0])
    c[fit_mask] = np.clip(c[fit_mask], float(clip_lo), float(clip_hi))

    if fit_cfg.get("residual_edge_clamp", True) and int(fit_mask.sum()) >= 3:
        idx_in = np.where(fit_mask)[0]
        c[idx_in[0]] = c[idx_in[1]]
        c[idx_in[-1]] = c[idx_in[-2]]

    # The legacy "integral-preserving alpha rescale" is removed: it was
    # both buggy (the rescale target was the stage-1 histogram after
    # renormalization to D.sum, but T_s/T_d were natural integrals) and
    # useless (M_stage2 is renormalized to D.sum right after, so any
    # constant rescale of c cancels). Keep c at its natural magnitude
    # (typically O(1)), which is also what downstream analysis wants
    # when it applies the per-vertex correction.
    alpha = 1.0

    # Build the final M_stage2 with the rescaled c
    sigma_gen, sigma_tgt = stage1.sigma_gen, stage1.sigma_tgt
    w_func = (lambda z, g=sigma_gen, s=sigma_tgt: fc.parametric_weight(z, g, s))
    c_func = fc.c_interpolator(inputs.binning.centers, c)
    f_s, f_d = inputs.period.f_single, inputs.period.f_double
    M2 = fc.mixed_mc_histogram(
        inputs.single_samples, inputs.double_samples,
        inputs.binning, f_s, f_d, w_func, c_func,
    )
    M2 = _normalize_to_total(M2, float(D.sum()))

    with np.errstate(divide="ignore", invalid="ignore"):
        resid2 = np.where(M2 > 0, D / M2, 1.0)
    max_res2_global = float(np.max(np.abs(resid2 - 1.0)))
    max_res2 = float(np.max(np.abs(resid2[fit_mask] - 1.0)))
    LOG.info("stage 2 final residual: max|D/M-1| = %.4f inside |z|<%.1f, %.4f global",
             max_res2, inputs.period.z_cut, max_res2_global)

    w_total = stage1.w_parametric * c
    return Stage2Result(
        c_raw=c_raw, c_smoothed=c, w_total=w_total,
        M_stage2=M2, residual_stage2=resid2, max_residual=max_res2,
        alpha=float(alpha),
    )


# ---------------------------------------------------------------------------
# Iterative data-driven stage (Lucy–Richardson-like multiplicative update)
# ---------------------------------------------------------------------------

@dataclass
class IterativeResult:
    sigma_seed_tgt: float                  # σ_tgt used for the parametric seed
    w_seed: np.ndarray                     # seed Gaussian weights at bin centers
    w_iterative: np.ndarray                # final w(z) at bin centers
    n_iterations: int
    converged: bool
    history: List[dict]                    # per-iter: max_residual, chi2_in_window
    M_iterative: np.ndarray                # final mixed MC histogram (renormalized)
    residual_iterative: np.ndarray         # final D / M
    max_residual_in_window: float
    chi2_in_window: float
    ndf_in_window: int


def _data_reco_rms(D: np.ndarray, centers: np.ndarray) -> float:
    """Weighted std of the data reco vertex distribution."""
    s = float(D.sum())
    if s <= 0:
        return 0.0
    mean = float((centers * D).sum() / s)
    var = float(((centers - mean) ** 2 * D).sum() / s)
    return float(np.sqrt(max(var, 0.0)))


def stage_iterative_fit(inputs: FitInputs, sigma_gen: float, config: dict) -> IterativeResult:
    fit_cfg = config["fit"]
    iter_cfg = fit_cfg.get("iterative", {})

    centers = inputs.binning.centers
    D = inputs.D
    f_s = inputs.period.f_single
    f_d = inputs.period.f_double
    z_cut = inputs.period.z_cut

    # Two masks:
    #   closure_mask: bins inside |z| < z_cut. Convergence and χ² are
    #                 evaluated only on this window.
    #   update_mask:  bins where w(z) is iteratively modified. Full
    #                 binning range. No post-iteration snap: the final
    #                 w(z) is data-driven everywhere (comparison study
    #                 showed no-snap improves 0 mrad closure by ~27% in
    #                 χ² with no degradation at 1.5 mrad).
    closure_mask = np.abs(centers) < z_cut
    update_mask = np.ones_like(centers, dtype=bool)

    # ----- Seed: Gaussian where σ_tgt = std(data reco) -----
    sigma_data = _data_reco_rms(D, centers)
    sigma_seed = sigma_data
    LOG.info("iterative seed: σ_data_reco = %.3f cm → using as σ_tgt for Gaussian seed",
             sigma_seed)
    w_seed = fc.parametric_weight(centers, sigma_gen, sigma_seed)
    w = w_seed.copy()

    # ----- Iteration parameters -----
    max_iter = int(iter_cfg.get("max_iterations", 30))
    alpha = float(iter_cfg.get("damping", 0.5))
    eps = float(iter_cfg.get("convergence_max_resid", 0.01))
    smooth_n = int(iter_cfg.get("residual_smooth_bins", 1))
    clip_lo, clip_hi = iter_cfg.get("update_clip", [0.5, 2.0])

    LOG.info("iterative fit: max_iter=%d damping=%.2f conv_eps=%.4f smooth=%d update_clip=[%.2f,%.2f]",
             max_iter, alpha, eps, smooth_n, clip_lo, clip_hi)
    LOG.info("iterative update window: %d bins (full binning ±%.1f cm)",
             int(update_mask.sum()), inputs.binning.z_max)
    LOG.info("iterative closure window: %d bins inside |z|<%.1f cm",
             int(closure_mask.sum()), z_cut)
    LOG.info("post-iteration: no snap; w(z) is data-driven over the full range")

    # Precompute single-MC truth-z for per-iteration normalization.
    # After each multiplicative update we rescale w so that <w(z)>_MC = 1,
    # which preserves the physics DI mix fraction in the forward model.
    _all_zh_single = np.concatenate([
        s.z_h[s.z_h > -9000] for s in inputs.single_samples if s.z_h is not None
    ])

    history: List[dict] = []
    D_total = float(D.sum())
    converged = False
    last_M = None
    last_resid = None
    last_max_res = float("inf")
    n_done = 0

    for it in range(max_iter):
        n_done = it + 1

        # 1. Build mixed-MC reco histogram with current w
        w_func = fc.c_interpolator(centers, w)
        M = fc.mixed_mc_histogram(
            inputs.single_samples, inputs.double_samples,
            inputs.binning, f_s, f_d, w_func,
        )
        M = _normalize_to_total(M, D_total)
        last_M = M

        # 2. Residual everywhere in the update window
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.ones_like(D, dtype=np.float64)
            valid = update_mask & (M > 0)
            r[valid] = D[valid] / M[valid]
        last_resid = r

        # Closure metrics: only on the narrow closure window
        max_res = float(np.max(np.abs(r[closure_mask] - 1.0))) if closure_mask.any() else 0.0
        chi2_inside = fc.chi2(D, M, mask=closure_mask)
        history.append({
            "iteration": it,
            "max_residual": max_res,
            "chi2_in_window": chi2_inside,
        })
        last_max_res = max_res

        LOG.info("iter %2d: max|r-1|=%.4f (closure window) chi2_in_window=%.2f",
                 it, max_res, chi2_inside)

        if max_res < eps:
            converged = True
            LOG.info("iterative: converged at iter %d (tolerance %.4f reached)", it, eps)
            break

        # 3. Smooth the residual over the update window
        r_smooth = r.copy()
        for _ in range(smooth_n):
            inner = r_smooth[update_mask]
            if len(inner) >= 3:
                r_smooth[update_mask] = _gaussian_smooth_1bin(inner)

        # 4. Damping + clip
        update = np.power(r_smooth[update_mask], alpha)
        update = np.clip(update, float(clip_lo), float(clip_hi))

        # 5. Multiplicative update — over the full update window
        w[update_mask] *= update

        # Edge clamp at the boundaries of the update window
        if int(update_mask.sum()) >= 3:
            idx_in = np.where(update_mask)[0]
            w[idx_in[0]] = w[idx_in[1]]
            w[idx_in[-1]] = w[idx_in[-2]]

        # 6. Renormalize w so <w(z)>_MC = 1 over the single-MC truth-z
        #    distribution. Without this, <w> drifts below 1 as the target
        #    is narrower than MC, silently reducing the effective DI fraction
        #    in the forward model (singles scale as <w>, doubles as <w>²).
        w_func_tmp = fc.c_interpolator(centers, w)
        mean_w = float(np.mean(w_func_tmp(_all_zh_single)))
        if mean_w > 1e-12:
            w /= mean_w

    if not converged:
        LOG.warning("iterative: did not converge after %d iterations (final max|r-1|=%.4f)",
                    max_iter, last_max_res)

    # No post-iteration snap: the iterative w(z) is kept over the full
    # binning range. See `compare_snap_strategies.py` for the study that
    # motivated this choice.

    ndf_in = int(closure_mask.sum()) - 1
    return IterativeResult(
        sigma_seed_tgt=sigma_seed,
        w_seed=w_seed,
        w_iterative=w,
        n_iterations=n_done,
        converged=converged,
        history=history,
        M_iterative=last_M if last_M is not None else np.zeros_like(D),
        residual_iterative=last_resid if last_resid is not None else np.ones_like(D),
        max_residual_in_window=last_max_res,
        chi2_in_window=float(history[-1]["chi2_in_window"]) if history else float("nan"),
        ndf_in_window=ndf_in,
    )


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _git_hash() -> str:
    try:
        h = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
            cwd=Path(__file__).resolve().parent,
        ).decode().strip()
        return h
    except Exception:
        return "unknown"


def save_yaml(path: Path, inputs: FitInputs, stage1: Stage1Result, stage2: Stage2Result,
              config: dict, sample_match: Optional[str] = None,
              iterative: Optional[IterativeResult] = None) -> None:
    doc = {
        "period": inputs.period.name,
        "run_range": [inputs.period.run_min, inputs.period.run_max],
        "lumi_pb_inv": inputs.period.lumi,
        "f_double": inputs.period.f_double,
        "git_hash": _git_hash(),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "sample_match": sample_match or "",
        "binning": {
            "z_min": inputs.binning.z_min,
            "z_max": inputs.binning.z_max,
            "n_bins": inputs.binning.n_bins,
            "centers": inputs.binning.centers.tolist(),
        },
        "stage1_parametric": {
            "sigma_gen": float(stage1.sigma_gen),
            "sigma_gen_err": float(stage1.sigma_gen_err),
            "sigma_tgt": float(stage1.sigma_tgt),
            "chi2": float(stage1.chi2),
            "ndf": int(stage1.ndf),
            "chi2_per_ndf": float(stage1.chi2 / max(stage1.ndf, 1)),
            "w_values": stage1.w_parametric.tolist(),
            "grid_sigma": stage1.grid_sigma.tolist(),
            "grid_chi2": stage1.grid_chi2.tolist(),
        },
        "stage2_residual": {
            "alpha": stage2.alpha,
            "c_raw": stage2.c_raw.tolist(),
            "c_smoothed": stage2.c_smoothed.tolist(),
            "w_total": stage2.w_total.tolist(),
            "residual_max": stage2.max_residual,
            "residual_final": stage2.residual_stage2.tolist(),
        },
        "sample_list": {
            "single_mc": [s.spec.name for s in inputs.single_samples],
            "double_mc": [s.spec.name for s in inputs.double_samples],
            "data": inputs.data_sample.spec.path,
        },
    }
    if iterative is not None:
        doc["iterative"] = {
            "sigma_seed_tgt": float(iterative.sigma_seed_tgt),
            "n_iterations": int(iterative.n_iterations),
            "converged": bool(iterative.converged),
            "max_residual_in_window": float(iterative.max_residual_in_window),
            "chi2_in_window": float(iterative.chi2_in_window),
            "ndf_in_window": int(iterative.ndf_in_window),
            "chi2_per_ndf": float(iterative.chi2_in_window / max(iterative.ndf_in_window, 1)),
            "w_seed": iterative.w_seed.tolist(),
            "w_iterative": iterative.w_iterative.tolist(),
            "residual_iterative": iterative.residual_iterative.tolist(),
            "history": iterative.history,
        }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.dump(doc, f, default_flow_style=False, sort_keys=False)
    LOG.info("wrote %s", path)


def save_root(path: Path, inputs: FitInputs, stage1: Stage1Result, stage2: Stage2Result,
              iterative: Optional[IterativeResult] = None) -> None:
    if not HAS_ROOT:
        LOG.warning("PyROOT not available; skipping %s", path)
        return
    import ROOT
    ROOT.gROOT.SetBatch(True)
    path.parent.mkdir(parents=True, exist_ok=True)
    f = ROOT.TFile(str(path), "RECREATE")

    n = inputs.binning.n_bins
    edges = inputs.binning.edges

    def _make_th1(name: str, title: str, values: np.ndarray):
        h = ROOT.TH1D(name, title, n, edges[0], edges[-1])
        for i, v in enumerate(values):
            h.SetBinContent(i + 1, float(v))
            h.SetBinError(i + 1, 0.0)
        h.SetDirectory(0)
        return h

    th1_list = [
        ("h_c_truth", "truth-vertex residual correction c(z)", stage2.c_smoothed),
        ("h_c_truth_raw", "raw c(z) before smoothing", stage2.c_raw),
        ("h_w_parametric", "parametric stage-1 w(z)", stage1.w_parametric),
        ("h_w_total", "total stage-1 * stage-2 w(z) * c(z)", stage2.w_total),
        ("h_D", "data reco vertex", inputs.D),
        ("h_M_stage1", "mixed MC after stage 1", stage1.M_stage1),
        ("h_M_stage2", "mixed MC after stage 2", stage2.M_stage2),
    ]
    if iterative is not None:
        th1_list += [
            ("h_w_seed", "iterative seed w(z) (Gaussian σ_tgt = std(D))", iterative.w_seed),
            ("h_w_iterative", "iterative data-driven w(z)", iterative.w_iterative),
            ("h_M_iterative", "mixed MC after iterative reweight", iterative.M_iterative),
        ]
    for name, title, values in th1_list:
        h = _make_th1(name, title, values)
        h.Write()

    meta_pairs = [
        ("sigma_gen", stage1.sigma_gen),
        ("sigma_tgt", stage1.sigma_tgt),
        ("alpha_stage2", stage2.alpha),
        ("period", inputs.period.name),
        ("run_min", inputs.period.run_min),
        ("run_max", inputs.period.run_max),
        ("f_double", inputs.period.f_double),
        ("z_cut", inputs.period.z_cut),
        ("lumi_pb_inv", inputs.period.lumi),
        ("git_hash", _git_hash()),
        ("timestamp", time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())),
        ("sample_list_single", ",".join(s.spec.name for s in inputs.single_samples)),
        ("sample_list_double", ",".join(s.spec.name for s in inputs.double_samples)),
    ]
    if iterative is not None:
        meta_pairs += [
            ("iter_sigma_seed", iterative.sigma_seed_tgt),
            ("iter_n_iterations", iterative.n_iterations),
            ("iter_converged", int(iterative.converged)),
            ("iter_max_residual_in_window", iterative.max_residual_in_window),
        ]
    for key, value in meta_pairs:
        ROOT.TNamed(key, str(value)).Write()

    f.Close()
    LOG.info("wrote %s", path)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run(config: dict, period_name: str, out_dir: Path, cache_dir: Path,
        sample_match: Optional[str], n_events: Optional[int],
        z_cut_override: Optional[float] = None) -> int:
    inputs = load_period_inputs(
        config, period_name, cache_dir=cache_dir,
        sample_match=sample_match, n_events=n_events,
        z_cut_override=z_cut_override,
    )

    stage1 = stage1_fit(inputs, config)
    stage2 = stage2_fit(inputs, stage1, config)

    if stage2.max_residual > 0.02:
        LOG.warning(
            "stage 2 post-fit residual is %.1f%% > 2%% threshold",
            100 * stage2.max_residual,
        )

    iterative: Optional[IterativeResult] = None
    iter_cfg = config["fit"].get("iterative", {})
    if iter_cfg.get("enabled", False):
        LOG.info("=== iterative data-driven stage ===")
        iterative = stage_iterative_fit(inputs, stage1.sigma_gen, config)
        if iterative.converged:
            LOG.info("iterative converged in %d iterations, final max|r-1|=%.4f",
                     iterative.n_iterations, iterative.max_residual_in_window)
        else:
            LOG.warning("iterative did not fully converge: %d iters, final max|r-1|=%.4f",
                        iterative.n_iterations, iterative.max_residual_in_window)

    save_yaml(out_dir / "fit_result.yaml", inputs, stage1, stage2, config, sample_match,
              iterative=iterative)
    save_root(out_dir / "reweight.root", inputs, stage1, stage2, iterative=iterative)

    # Warn but do not exit non-zero on large χ²/ndf — we still want the
    # closure plot and YAML so we can diagnose the mismatch.
    if not np.isfinite(stage1.chi2):
        LOG.error("stage 1 chi2 is not finite — fit failed")
        return 10
    if stage1.chi2 / max(stage1.ndf, 1) > 5:
        LOG.warning("stage 1 chi2/ndf = %.1f > 5 — parametric model does not close",
                    stage1.chi2 / max(stage1.ndf, 1))
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(here / "config.yaml"))
    ap.add_argument("--period", required=True, choices=("1p5mrad", "0mrad"))
    ap.add_argument("--match", default=None,
                    help="only load MC samples whose name contains this substring")
    ap.add_argument("--n-events", type=int, default=None,
                    help="limit events per sample (for quick test runs)")
    ap.add_argument("--z-cut", type=float, default=None,
                    help="override the z_cut closure window (cm); bypasses "
                         "z_cut / z_cut_percentile in the config")
    ap.add_argument("--out-dir", default=None,
                    help="output directory (default: output/<period>/)")
    ap.add_argument("--cache-dir", default=None,
                    help="cache directory for .npz event arrays "
                         "(default: output/cache/)")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    fc.LOG.setLevel(logging.INFO)

    config = load_config(Path(args.config))

    out_base = config["output"]["base_dir"]
    out_base_path = (here / out_base) if not Path(out_base).is_absolute() else Path(out_base)
    out_dir = Path(args.out_dir) if args.out_dir else (out_base_path / args.period)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(args.cache_dir) if args.cache_dir else (out_base_path / "cache")
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Tee log
    log_path = out_dir / "fit.log"
    fh = logging.FileHandler(log_path, mode="w")
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%H:%M:%S"))
    logging.getLogger().addHandler(fh)

    LOG.info("=== fit_truth_vertex_reweight %s ===", args.period)
    LOG.info("config: %s", args.config)
    LOG.info("output: %s", out_dir)
    LOG.info("cache:  %s", cache_dir)
    if args.match:
        LOG.info("sample filter: substring %r", args.match)

    try:
        return run(config, args.period, out_dir, cache_dir,
                   sample_match=args.match, n_events=args.n_events,
                   z_cut_override=args.z_cut)
    except FileNotFoundError as e:
        LOG.error("%s", e)
        return 2
    except Exception as e:
        LOG.exception("fit failed: %s", e)
        return 3


if __name__ == "__main__":
    sys.exit(main())
