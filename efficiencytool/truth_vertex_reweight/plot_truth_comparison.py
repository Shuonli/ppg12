#!/usr/bin/env python3
"""
Compare truth vertex z_h distributions under three reweighting methods:
  1. No reweight (raw MC)
  2. Old production reco-level f(z_r) = D(z_r)/M(z_r)
  3. New iterative w(z_h) [single: w(z_h), double: w(z_h)*w(z_mb)]

Produces a 2-panel figure (1.5 mrad | 0 mrad) showing how the iterative
method properly narrows the truth vertex distribution, while the reco-level
method leaves double-interaction MC at the "half-narrow, half-wide" compromise.
"""
from __future__ import annotations

import argparse
import logging
import math
import shutil
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_core as fc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sphenix_style import set_sphenix_rcparams, sphenix_label, sphenix_period_text
set_sphenix_rcparams()

LOG = logging.getLogger("plot_truth_comparison")


def _valid(z):
    return z > -9000.0


def _weighted_rms(centers, h):
    s = h.sum()
    if s <= 0:
        return 0.0
    return float(np.sqrt(np.average(centers ** 2, weights=h)))


def _step_hist(ax, centers, values, **kw):
    if len(centers) < 2:
        return
    dz = centers[1] - centers[0]
    edges = np.zeros(len(centers) + 1)
    edges[:-1] = centers - dz / 2.0
    edges[-1] = centers[-1] + dz / 2.0
    ax.step(edges, np.concatenate([values, [values[-1]]]), where="post", **kw)


def _build_truth_hist(
    single_samples, double_samples, binning, f_s, f_d,
    truth_w_func=None, reco_w_func=None,
):
    """Histogram z_h (primary truth vertex) for mixed MC under a given
    weighting scheme. Events without a valid reco vertex are included
    only when no reco_w_func is applied (matching the analysis treatment)."""
    edges = np.linspace(binning.z_min, binning.z_max, binning.n_bins + 1)
    h = np.zeros(binning.n_bins, dtype=np.float64)

    for s in single_samples:
        zh = s.z_h
        zr = s.z_r
        mask = _valid(zh)
        if reco_w_func is not None:
            mask &= _valid(zr)
        zh_m, zr_m = zh[mask], zr[mask]
        w = np.full(len(zh_m), s.spec.event_weight * f_s, dtype=np.float64)
        if truth_w_func is not None:
            w *= truth_w_func(zh_m)
        if reco_w_func is not None:
            w *= reco_w_func(zr_m)
        h += np.histogram(zh_m, bins=edges, weights=w)[0]

    for s in double_samples:
        zh = s.z_h
        zr = s.z_r
        zmb = s.z_mb
        mask = _valid(zh) & _valid(zr)
        if zmb is not None:
            mask &= _valid(zmb)
        zh_m, zr_m = zh[mask], zr[mask]
        zmb_m = zmb[mask] if zmb is not None else None
        w = np.full(len(zh_m), s.spec.event_weight * f_d, dtype=np.float64)
        if truth_w_func is not None:
            w *= truth_w_func(zh_m)
            if zmb_m is not None:
                w *= truth_w_func(zmb_m)
        if reco_w_func is not None:
            w *= reco_w_func(zr_m)
        h += np.histogram(zh_m, bins=edges, weights=w)[0]

    return h


def plot(config, periods, cache_dir, out_base, sample_match):
    from fit_truth_vertex_reweight import load_period_inputs

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))

    for ax, period in zip(axes, periods):
        LOG.info("processing %s ...", period)

        inputs = load_period_inputs(
            config, period, cache_dir=cache_dir,
            sample_match=sample_match, n_events=None, z_cut_override=None,
        )

        yaml_path = out_base / period / "fit_result.yaml"
        with open(yaml_path) as f:
            doc = yaml.safe_load(f)

        bnode = doc["binning"]
        centers = np.asarray(bnode["centers"])
        binning = fc.Binning(z_min=bnode["z_min"], z_max=bnode["z_max"],
                             n_bins=bnode["n_bins"])

        iter_doc = doc.get("iterative")
        if iter_doc is None:
            LOG.warning("no iterative result for %s, skipping", period)
            continue

        w_iter = np.asarray(iter_doc["w_iterative"])
        w_iter_func = fc.c_interpolator(centers, w_iter)
        sigma_seed = float(iter_doc["sigma_seed_tgt"])

        D = inputs.D.astype(np.float64)
        f_s, f_d = inputs.period.f_single, inputs.period.f_double

        def _w_unit(z):
            return np.ones_like(z, dtype=np.float64)

        M_noRW = fc.mixed_mc_histogram(
            inputs.single_samples, inputs.double_samples,
            binning, f_s, f_d, _w_unit,
        )

        D_n = D / max(D.sum(), 1e-30)
        M_n = M_noRW / max(M_noRW.sum(), 1e-30)
        f_reco_arr = np.where(M_n > 0, D_n / M_n, 1.0)
        f_reco_func = fc.c_interpolator(centers, f_reco_arr)

        h_noRW = _build_truth_hist(
            inputs.single_samples, inputs.double_samples,
            binning, f_s, f_d,
        )
        h_reco = _build_truth_hist(
            inputs.single_samples, inputs.double_samples,
            binning, f_s, f_d, reco_w_func=f_reco_func,
        )
        h_iter = _build_truth_hist(
            inputs.single_samples, inputs.double_samples,
            binning, f_s, f_d, truth_w_func=w_iter_func,
        )

        for h in (h_noRW, h_reco, h_iter):
            s = h.sum()
            if s > 0:
                h /= s

        rms_noRW = _weighted_rms(centers, h_noRW)
        rms_reco = _weighted_rms(centers, h_reco)
        rms_iter = _weighted_rms(centers, h_iter)

        _step_hist(ax, centers, h_noRW, color="C7", lw=1.0, ls=":",
                   label=fr"no reweight ($\sigma\approx${rms_noRW:.0f} cm)")
        _step_hist(ax, centers, h_reco, color="C0", lw=1.5, ls="--",
                   label=fr"reco $f(z_r)$ ($\sigma\approx${rms_reco:.0f} cm)")
        _step_hist(ax, centers, h_iter, color="C3", lw=1.8,
                   label=fr"iterative $w(z_h)$ ($\sigma\approx${rms_iter:.0f} cm)")

        ax.set_xlabel(r"truth vertex $z_h$ (cm)")
        ax.set_ylabel("area-normalized")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(alpha=0.3)
        ax.set_xlim(-200, 200)

        sphenix_label(ax, period=period, loc="upper left", size=9.5)

    fig.suptitle(
        r"Truth vertex $z_h$ — reco-level $f(z_r)$ vs iterative $w(z_h)$ reweight"
        "  (photon10 mix)",
        fontsize=11, y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.955])

    out_path = out_base / "truth_vertex_comparison.pdf"
    fig.savefig(out_path)
    plt.close(fig)
    LOG.info("wrote %s", out_path)

    dst = Path(
        "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/"
        "reports/figures/vertex_reweight_truth_comparison.pdf"
    )
    shutil.copy2(out_path, dst)
    LOG.info("staged to %s", dst)
    return 0


def main(argv=None):
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(here / "config.yaml"))
    ap.add_argument("--match", default="photon10")
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

    out_base_str = config["output"]["base_dir"]
    out_base = (here / out_base_str) if not Path(out_base_str).is_absolute() else Path(out_base_str)
    cache_dir = Path(args.cache_dir) if args.cache_dir else (out_base / "cache")
    return plot(config, ["1p5mrad", "0mrad"], cache_dir, out_base, args.match)


if __name__ == "__main__":
    sys.exit(main())
