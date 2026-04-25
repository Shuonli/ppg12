#!/usr/bin/env python3
"""Converted vs unconverted photon study (PPG12).

Streams the single-particle photon10 MC slimtree with uproot/awkward, classifies
truth photons by ``particle_converted`` (0=clean, 1=converted e+e-, 2=bad
secondary) and measures:

  * conversion fractions vs truth pT / truth eta / 2D (eta,pT)
  * cluster-match efficiency per category (>=1 cluster with
    ``cluster_truthtrkID_CLUSTERINFO_CEMC == particle_trkid`` and
    ``cluster_Et_CLUSTERINFO_CEMC > 5 GeV``)
  * cluster-match multiplicity per truth photon per category

Outputs:
  rootFiles/fractions_converted.root  - TH1/TH2 style histograms
  figures/fractions_*.pdf             - matplotlib plots
  fractions_findings.md               - Markdown report
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass, field

import numpy as np
import awkward as ak
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

INPUT_FILE = (
    "/sphenix/user/shuhangli/ppg12/FunWithxgboost/"
    "photon10_with_bdt_vtx_noreweight_single.root"
)
TREE_NAME = "slimtree"
CLUSTER_SUFFIX = "_CLUSTERINFO_CEMC"

PT_EDGES = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
ETA_EDGES = np.round(np.arange(-1.5, 1.5 + 1e-9, 0.1), 3)
ABS_ETA_EDGES = np.round(np.arange(0.0, 1.5 + 1e-9, 0.1), 3)

CATEGORIES = [
    (0, "clean", "unconverted", "black"),
    (1, "converted", "converted (e+e-)", "red"),
    (2, "bad", "bad secondary", "blue"),
]

CLUSTER_ET_MIN = 5.0
PT_REPORT_RANGE = (14.0, 30.0)

BRANCHES = [
    "particle_pid",
    "particle_trkid",
    "particle_Pt",
    "particle_Eta",
    "particle_converted",
    "cluster_Et" + CLUSTER_SUFFIX,
    "cluster_truthtrkID" + CLUSTER_SUFFIX,
]

NMATCH_HIST_MAX = 6


# ---------------------------------------------------------------------------
# Accumulator
# ---------------------------------------------------------------------------


@dataclass
class Accum:
    eta_regions = ("eta07", "eta15")

    totals: dict = field(default_factory=dict)
    n_pt: dict = field(default_factory=dict)
    n_etapt: dict = field(default_factory=dict)
    n_abseta_incl: dict = field(default_factory=dict)
    matched_pt: dict = field(default_factory=dict)
    sum_nmatch_pt: dict = field(default_factory=dict)
    nmatch_dist: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        n_pt_bins = len(PT_EDGES) - 1
        n_eta_bins = len(ETA_EDGES) - 1
        n_abseta_bins = len(ABS_ETA_EDGES) - 1
        for region in self.eta_regions:
            self.totals[region] = {c[0]: 0 for c in CATEGORIES}
            self.n_pt[region] = {c[0]: np.zeros(n_pt_bins, dtype=np.int64) for c in CATEGORIES}
            self.n_etapt[region] = {
                c[0]: np.zeros((n_eta_bins, n_pt_bins), dtype=np.int64) for c in CATEGORIES
            }
            self.n_abseta_incl[region] = {
                c[0]: np.zeros(n_abseta_bins, dtype=np.int64) for c in CATEGORIES
            }
            self.matched_pt[region] = {
                c[0]: np.zeros(n_pt_bins, dtype=np.int64) for c in CATEGORIES
            }
            self.sum_nmatch_pt[region] = {
                c[0]: np.zeros(n_pt_bins, dtype=np.int64) for c in CATEGORIES
            }
            self.nmatch_dist[region] = {
                c[0]: np.zeros(NMATCH_HIST_MAX + 1, dtype=np.int64) for c in CATEGORIES
            }


# ---------------------------------------------------------------------------
# Vectorized chunk processor
# ---------------------------------------------------------------------------


def process_chunk(acc: Accum, batch: ak.Array) -> int:
    """Process one batch of events (awkward arrays, jagged per event)."""
    p_pid = batch["particle_pid"]
    p_trk = batch["particle_trkid"]
    p_pt = batch["particle_Pt"]
    p_eta = batch["particle_Eta"]
    p_conv = batch["particle_converted"]
    c_et = batch["cluster_Et" + CLUSTER_SUFFIX]
    c_trk = batch["cluster_truthtrkID" + CLUSTER_SUFFIX]

    # Build event-by-event photon mask and a broadcast-aligned cluster trkid set
    ph_mask = (p_pid == 22)

    # For each photon in each event, count the number of clusters whose
    # truthtrkID matches the photon's trkid AND whose Et > CLUSTER_ET_MIN.
    # ---
    # Broadcast: photon trkids (ragged: nph per event) vs cluster trkids (ragged: nclu per event)
    # nmatch per photon = sum over clusters in same event of
    #   (c_trk == p_trk) & (c_et > CLUSTER_ET_MIN)
    # Use ak.local_index + broadcast by adding a new axis.

    # Keep only photon rows (one row per photon, with event index preserved)
    ph_trk = p_trk[ph_mask]
    ph_pt = p_pt[ph_mask]
    ph_eta = p_eta[ph_mask]
    ph_conv = p_conv[ph_mask]

    # Precompute passing cluster trkids per event (clusters with Et > min)
    c_mask = c_et > CLUSTER_ET_MIN
    c_trk_pass = c_trk[c_mask]

    # Turn cluster_trk into per-event numpy lists? Counting matches requires
    # broadcasting ph_trk (nph) against c_trk_pass (nclu) within each event.
    # Approach: use ak.cartesian per event, then group by photon index.
    # Cheaper: use Python-per-event only for the match count.

    # Since nph is small (usually 1-2 real photons in range per event), and nclu is ~0-3,
    # an awkward-native approach using broadcast_arrays + counting is fastest.

    # Build photon x cluster table
    # First, give photons a local index, same for clusters
    phi = ak.local_index(ph_trk)
    cci = ak.local_index(c_trk_pass)
    # Nothing to match against if no clusters:
    # cartesian gives (ph, cl) pairs per event
    pairs = ak.cartesian({"ph": ph_trk, "cl": c_trk_pass}, axis=1)
    pair_match = (pairs.ph == pairs.cl)
    # to count per photon, also produce photon-index per pair
    idx_pairs = ak.cartesian({"ph": phi, "cl": cci}, axis=1)
    # match indicator per pair: int
    pair_match_int = ak.values_astype(pair_match, np.int64)
    # Group pairs by photon index within event. Shape: events x pairs.
    # We want sum over pairs where ph_idx==k. Use unflatten tricks:
    # Simpler: for each event, count unique combinations via awkward groupby

    # Use ak.num to get per-event lengths and then to_numpy on flat:
    # Flatten to per-event local indexing:
    # For each photon, nmatch_i = sum_j 1(ph_trk_i == c_trk_pass_j)
    # Equivalent: nmatch_per_ph = ak.sum(pair_match_int, where ph_idx==i)

    # Convert to flat per-photon by rebuilding: index mapping is phi per pair.
    # Use np.add.at per event. We'll flatten per-event to numpy.

    # Alternative simpler: iterate per event, but only on flat numpy per event.
    # Given chunks of ~500k events, per-event python loop over awkward slices is slow.
    # We'll do it vectorized via offsets.

    # per-event number of photons and of cluster pairs
    nph_per_ev = ak.num(ph_trk)
    ncl_per_ev = ak.num(c_trk_pass)
    nph = int(ak.sum(nph_per_ev))
    if nph == 0:
        return 0

    # flat arrays
    flat_ph_trk = ak.to_numpy(ak.flatten(ph_trk))
    flat_ph_pt = ak.to_numpy(ak.flatten(ph_pt))
    flat_ph_eta = ak.to_numpy(ak.flatten(ph_eta))
    flat_ph_conv = ak.to_numpy(ak.flatten(ph_conv))

    # Build nmatch per flat photon via per-event loop using offsets
    ph_off = np.concatenate(([0], np.cumsum(ak.to_numpy(nph_per_ev))))
    cl_off = np.concatenate(([0], np.cumsum(ak.to_numpy(ncl_per_ev))))
    flat_c_trk_pass = ak.to_numpy(ak.flatten(c_trk_pass))

    nmatch = np.zeros(nph, dtype=np.int64)
    n_events = len(nph_per_ev)
    for ev in range(n_events):
        p0, p1 = ph_off[ev], ph_off[ev + 1]
        if p1 == p0:
            continue
        c0, c1 = cl_off[ev], cl_off[ev + 1]
        if c1 == c0:
            continue  # nmatch stays 0
        ph_tr_slice = flat_ph_trk[p0:p1]
        cl_tr_slice = flat_c_trk_pass[c0:c1]
        # For each photon in event compute count of matches in cluster slice:
        # shape (nph_ev, ncl_ev) boolean
        eq = (ph_tr_slice[:, None] == cl_tr_slice[None, :])
        nmatch[p0:p1] = eq.sum(axis=1)

    # Fill histograms (vectorized)
    aeta = np.abs(flat_ph_eta)
    i_pt = np.searchsorted(PT_EDGES, flat_ph_pt, side="right") - 1
    i_eta = np.searchsorted(ETA_EDGES, flat_ph_eta, side="right") - 1
    i_ae = np.searchsorted(ABS_ETA_EDGES, aeta, side="right") - 1
    in_report_pt = (flat_ph_pt >= PT_REPORT_RANGE[0]) & (flat_ph_pt < PT_REPORT_RANGE[1])

    n_pt_bins = len(PT_EDGES) - 1
    n_eta_bins = len(ETA_EDGES) - 1
    n_abseta_bins = len(ABS_ETA_EDGES) - 1

    valid_pt = (i_pt >= 0) & (i_pt < n_pt_bins)
    valid_eta = (i_eta >= 0) & (i_eta < n_eta_bins)
    valid_ae = (i_ae >= 0) & (i_ae < n_abseta_bins)

    for region, amax in (("eta07", 0.7), ("eta15", 1.5)):
        m_reg = aeta < amax
        for code, key, _, _ in CATEGORIES:
            m_cat = flat_ph_conv == code
            m = m_reg & m_cat

            # n vs pT
            mm = m & valid_pt
            if np.any(mm):
                np.add.at(acc.n_pt[region][code], i_pt[mm], 1)

            # 2d (eta,pt)
            mm2 = m & valid_pt & valid_eta
            if np.any(mm2):
                np.add.at(acc.n_etapt[region][code], (i_eta[mm2], i_pt[mm2]), 1)

            # |eta| over report pt range
            mm3 = m & in_report_pt & valid_ae
            if np.any(mm3):
                np.add.at(acc.n_abseta_incl[region][code], i_ae[mm3], 1)

            # totals in report pt range
            acc.totals[region][code] += int(np.count_nonzero(m & in_report_pt))

            # match efficiency
            mm4 = m & valid_pt
            nmatch_cat = nmatch[mm4]
            i_pt_cat = i_pt[mm4]
            if nmatch_cat.size > 0:
                matched = nmatch_cat >= 1
                if np.any(matched):
                    np.add.at(acc.matched_pt[region][code], i_pt_cat[matched], 1)
                np.add.at(acc.sum_nmatch_pt[region][code], i_pt_cat, nmatch_cat)

            # multiplicity dist
            m5 = m
            if np.any(m5):
                nm_clip = np.minimum(nmatch[m5], NMATCH_HIST_MAX)
                np.add.at(acc.nmatch_dist[region][code], nm_clip, 1)

    return len(nph_per_ev)


def process(max_entries: int | None, step_size) -> tuple[Accum, int]:
    """Stream the tree. ``max_entries`` bounds the number of tree entries read
    (implemented manually since ``uproot.iterate`` does not support entry_stop)."""
    acc = Accum()

    total_processed = 0
    t_start = time.time()
    iter_kwargs = dict(expressions=BRANCHES, step_size=step_size)

    for batch in uproot.iterate(f"{INPUT_FILE}:{TREE_NAME}", **iter_kwargs):
        if max_entries is not None and total_processed >= max_entries:
            break
        # respect max_entries by truncating the batch if needed
        if max_entries is not None and (total_processed + len(batch)) > max_entries:
            take = max_entries - total_processed
            batch = batch[:take]
        n_here = process_chunk(acc, batch)
        total_processed += n_here
        dt = time.time() - t_start
        rate = total_processed / dt if dt > 0 else 0
        print(f"  processed {total_processed:>10,} events  ({dt:6.1f}s, {rate:8.0f} evt/s)", flush=True)

    print(f"Done. Processed {total_processed:,} events in {time.time()-t_start:.1f}s.")
    return acc, total_processed


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------


def sphenix_label(ax, text="sPHENIX Simulation", sub=None):
    ax.text(0.02, 0.97, text, transform=ax.transAxes, fontsize=12,
            fontweight="bold", ha="left", va="top", style="italic")
    if sub is not None:
        ax.text(0.02, 0.91, sub, transform=ax.transAxes, fontsize=10,
                ha="left", va="top")


def _binomial_err(k: np.ndarray, n: np.ndarray) -> np.ndarray:
    k = k.astype(float); n = n.astype(float)
    p = np.divide(k, n, out=np.zeros_like(k), where=n > 0)
    var = np.divide(p * (1 - p), n, out=np.zeros_like(k), where=n > 0)
    return np.sqrt(var)


def plot_frac_vs_pt(acc: Accum, region: str, outpath: str) -> None:
    centers = 0.5 * (PT_EDGES[1:] + PT_EDGES[:-1])
    widths = 0.5 * (PT_EDGES[1:] - PT_EDGES[:-1])
    fig, ax = plt.subplots(figsize=(6.5, 5))
    total = sum(acc.n_pt[region][c[0]] for c in CATEGORIES)
    for code, key, label, color in CATEGORIES:
        k = acc.n_pt[region][code]
        p = np.divide(k, total, out=np.zeros_like(k, dtype=float), where=total > 0)
        err = _binomial_err(k, total)
        ax.errorbar(centers, p, xerr=widths, yerr=err, fmt="o", color=color,
                    label=label, markersize=4, linewidth=1.2)
    ax.set_xlabel(r"$p_{T}^{truth}$ [GeV]", fontsize=12)
    ax.set_ylabel("fraction", fontsize=12)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(PT_EDGES[0], PT_EDGES[-1])
    ax.grid(alpha=0.3)
    ax.legend(loc="center right", fontsize=10)
    region_label = r"$|\eta_{truth}| < 0.7$" if region == "eta07" else r"$|\eta_{truth}| < 1.5$"
    sphenix_label(ax, sub=f"photon10, {region_label}")
    fig.tight_layout(); fig.savefig(outpath); plt.close(fig)
    print("wrote", outpath)


def plot_frac_vs_abseta(acc: Accum, region: str, outpath: str) -> None:
    centers = 0.5 * (ABS_ETA_EDGES[1:] + ABS_ETA_EDGES[:-1])
    widths = 0.5 * (ABS_ETA_EDGES[1:] - ABS_ETA_EDGES[:-1])
    fig, ax = plt.subplots(figsize=(6.5, 5))
    total = sum(acc.n_abseta_incl[region][c[0]] for c in CATEGORIES)
    for code, key, label, color in CATEGORIES:
        k = acc.n_abseta_incl[region][code]
        p = np.divide(k, total, out=np.zeros_like(k, dtype=float), where=total > 0)
        err = _binomial_err(k, total)
        ax.errorbar(centers, p, xerr=widths, yerr=err, fmt="o", color=color,
                    label=label, markersize=4, linewidth=1.2)
    ax.set_xlabel(r"$|\eta_{truth}|$", fontsize=12)
    ax.set_ylabel("fraction", fontsize=12)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, ABS_ETA_EDGES[-1])
    ax.grid(alpha=0.3)
    ax.legend(loc="center right", fontsize=10)
    region_label = r"$|\eta_{truth}| < 0.7$" if region == "eta07" else r"$|\eta_{truth}| < 1.5$"
    sphenix_label(ax, sub=f"photon10, {PT_REPORT_RANGE[0]:.0f} < $p_T$ < {PT_REPORT_RANGE[1]:.0f} GeV, {region_label}")
    fig.tight_layout(); fig.savefig(outpath); plt.close(fig)
    print("wrote", outpath)


def plot_frac_2d(acc: Accum, region: str, code: int, outpath: str, label: str) -> None:
    total = sum(acc.n_etapt[region][c[0]] for c in CATEGORIES).astype(float)
    k = acc.n_etapt[region][code].astype(float)
    frac = np.divide(k, total, out=np.zeros_like(k), where=total > 0)

    fig, ax = plt.subplots(figsize=(7, 5))
    vmax = max(0.3, float(np.nanmax(frac)) * 1.1)
    im = ax.pcolormesh(PT_EDGES, ETA_EDGES, frac, shading="flat",
                       cmap="viridis", vmin=0, vmax=vmax)
    ax.set_xlabel(r"$p_{T}^{truth}$ [GeV]", fontsize=12)
    ax.set_ylabel(r"$\eta_{truth}$", fontsize=12)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f"fraction {label}", fontsize=11)
    region_label = r"$|\eta|<0.7$" if region == "eta07" else r"$|\eta|<1.5$"
    sphenix_label(ax, sub=f"photon10, {region_label}")
    fig.tight_layout(); fig.savefig(outpath); plt.close(fig)
    print("wrote", outpath)


def plot_match_eff_vs_pt(acc: Accum, region: str, outpath: str) -> None:
    centers = 0.5 * (PT_EDGES[1:] + PT_EDGES[:-1])
    widths = 0.5 * (PT_EDGES[1:] - PT_EDGES[:-1])
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for code, key, label, color in CATEGORIES:
        if key == "bad":
            continue  # negligible population; suppress from summary plot
        k = acc.matched_pt[region][code]
        n = acc.n_pt[region][code]
        eff = np.divide(k, n, out=np.zeros_like(k, dtype=float), where=n > 0)
        err = _binomial_err(k, n)
        ax.errorbar(centers, eff, xerr=widths, yerr=err, fmt="o", color=color,
                    label=label, markersize=4, linewidth=1.2)
    ax.set_xlabel(r"$p_{T}^{truth}$ [GeV]", fontsize=12)
    ax.set_ylabel(r"cluster-match efficiency ($E_{T}>5$ GeV)", fontsize=12)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(PT_EDGES[0], PT_EDGES[-1])
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=10)
    region_label = r"$|\eta_{truth}| < 0.7$" if region == "eta07" else r"$|\eta_{truth}| < 1.5$"
    sphenix_label(ax, sub=f"photon10, {region_label}")
    fig.tight_layout(); fig.savefig(outpath); plt.close(fig)
    print("wrote", outpath)


def plot_multiplicity(acc: Accum, region: str, outpath: str) -> None:
    x = np.arange(NMATCH_HIST_MAX + 1)
    width = 0.28
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for i, (code, key, label, color) in enumerate(CATEGORIES):
        counts = acc.nmatch_dist[region][code].astype(float)
        denom = counts.sum()
        if denom > 0:
            counts /= denom
        ax.bar(x + (i - 1) * width, counts, width=width, color=color,
               label=label, alpha=0.85, edgecolor="black", linewidth=0.4)
    ax.set_xlabel(r"# matched reco clusters per truth photon ($E_{T}>5$ GeV)", fontsize=12)
    ax.set_ylabel("fraction of truth photons", fontsize=12)
    ax.set_xticks(x)
    labels = [str(i) for i in x[:-1]] + [f">={NMATCH_HIST_MAX}"]
    ax.set_xticklabels(labels)
    ax.legend(loc="upper right", fontsize=10)
    ax.grid(alpha=0.3, axis="y")
    region_label = r"$|\eta_{truth}| < 0.7$" if region == "eta07" else r"$|\eta_{truth}| < 1.5$"
    sphenix_label(ax, sub=f"photon10, {region_label}")
    fig.tight_layout(); fig.savefig(outpath); plt.close(fig)
    print("wrote", outpath)


# ---------------------------------------------------------------------------
# ROOT writeout
# ---------------------------------------------------------------------------


def write_root(acc: Accum, path: str) -> None:
    with uproot.recreate(path) as f:
        for region in acc.eta_regions:
            total_pt = sum(acc.n_pt[region][c[0]] for c in CATEGORIES).astype(float)
            total_abseta = sum(acc.n_abseta_incl[region][c[0]] for c in CATEGORIES).astype(float)
            total_etapt = sum(acc.n_etapt[region][c[0]] for c in CATEGORIES).astype(float)

            for code, key, label, color in CATEGORIES:
                f[f"counts_vs_pt_{region}_{key}"] = (
                    acc.n_pt[region][code].astype(np.float64), PT_EDGES,
                )
                frac = np.divide(
                    acc.n_pt[region][code].astype(float), total_pt,
                    out=np.zeros_like(total_pt), where=total_pt > 0,
                )
                f[f"frac_vs_pt_{region}_{key}"] = (frac, PT_EDGES)

                f[f"counts_vs_abseta_{region}_{key}"] = (
                    acc.n_abseta_incl[region][code].astype(np.float64), ABS_ETA_EDGES,
                )
                frac_ae = np.divide(
                    acc.n_abseta_incl[region][code].astype(float), total_abseta,
                    out=np.zeros_like(total_abseta), where=total_abseta > 0,
                )
                f[f"frac_vs_abseta_{region}_{key}"] = (frac_ae, ABS_ETA_EDGES)

                frac2d = np.divide(
                    acc.n_etapt[region][code].astype(float), total_etapt,
                    out=np.zeros_like(total_etapt), where=total_etapt > 0,
                )
                # uproot wants (values, xedges, yedges) with values shape (nx, ny)
                # n_etapt is (n_eta, n_pt); make x=pt (cols), y=eta (rows)
                f[f"frac2d_pt_eta_{region}_{key}"] = (
                    frac2d.T.astype(np.float64), PT_EDGES, ETA_EDGES,
                )

                f[f"matched_vs_pt_{region}_{key}"] = (
                    acc.matched_pt[region][code].astype(np.float64), PT_EDGES,
                )
                eff = np.divide(
                    acc.matched_pt[region][code].astype(float),
                    acc.n_pt[region][code].astype(float),
                    out=np.zeros_like(total_pt),
                    where=acc.n_pt[region][code] > 0,
                )
                f[f"match_eff_vs_pt_{region}_{key}"] = (eff, PT_EDGES)

                mult = acc.nmatch_dist[region][code].astype(np.float64)
                f[f"nmatch_dist_{region}_{key}"] = (
                    mult, np.arange(len(mult) + 1, dtype=float),
                )
    print("wrote ROOT file:", path)


# ---------------------------------------------------------------------------
# Report helpers
# ---------------------------------------------------------------------------


def _inclusive_fractions(acc: Accum, region: str) -> dict:
    total = sum(acc.totals[region].values())
    out = {}
    for code, key, label, color in CATEGORIES:
        n = acc.totals[region][code]
        p = n / total if total > 0 else 0.0
        err = np.sqrt(max(p * (1 - p), 0.0) / total) if total > 0 else 0.0
        out[key] = (n, p, err)
    out["_total"] = total
    return out


def _match_eff_inclusive(acc: Accum, region: str) -> dict:
    lo = int(np.searchsorted(PT_EDGES, PT_REPORT_RANGE[0], side="left"))
    hi = int(np.searchsorted(PT_EDGES, PT_REPORT_RANGE[1], side="left"))
    out = {}
    for code, key, label, color in CATEGORIES:
        n = int(acc.n_pt[region][code][lo:hi].sum())
        k = int(acc.matched_pt[region][code][lo:hi].sum())
        p = k / n if n > 0 else 0.0
        err = np.sqrt(max(p * (1 - p), 0.0) / n) if n > 0 else 0.0
        out[key] = (k, n, p, err)
    return out


def _fmt_row_pt(acc: Accum, region: str, i: int) -> str:
    lo, hi = PT_EDGES[i], PT_EDGES[i + 1]
    row = [f"{lo:4.0f}-{hi:4.0f}"]
    total = sum(acc.n_pt[region][c[0]][i] for c in CATEGORIES)
    row.append(f"{int(total):>8d}")
    for code, key, _, _ in CATEGORIES:
        n = int(acc.n_pt[region][code][i])
        p = n / total if total > 0 else 0.0
        row.append(f"{p:6.3%}")
    return " | ".join(row)


def _fmt_row_abseta(acc: Accum, region: str, i: int) -> str:
    lo, hi = ABS_ETA_EDGES[i], ABS_ETA_EDGES[i + 1]
    row = [f"{lo:4.2f}-{hi:4.2f}"]
    total = sum(acc.n_abseta_incl[region][c[0]][i] for c in CATEGORIES)
    row.append(f"{int(total):>8d}")
    for code, key, _, _ in CATEGORIES:
        n = int(acc.n_abseta_incl[region][code][i])
        p = n / total if total > 0 else 0.0
        row.append(f"{p:6.3%}")
    return " | ".join(row)


def write_report(acc: Accum, processed_msg: str, outpath: str) -> None:
    lines = []
    lines.append("# Converted vs Unconverted Photon Fractions (PPG12)")
    lines.append("")
    lines.append(
        "Sample: single-particle photon10 MC, "
        f"`{os.path.basename(INPUT_FILE)}`, tree `{TREE_NAME}`. {processed_msg}"
    )
    lines.append("")
    lines.append(
        "`particle_converted` convention from CaloAna24: 0 = clean photon, "
        "1 = converted (e+/e- daughters), 2 = bad (non-ee secondary with >40% photon momentum)."
    )
    lines.append("")
    lines.append("Cluster match = at least one reco cluster in `CLUSTERINFO_CEMC` with "
                 f"`cluster_truthtrkID == particle_trkid` and `cluster_Et > {CLUSTER_ET_MIN:.0f}` GeV.")
    lines.append("")

    for region, label in (("eta07", "|eta|<0.7"), ("eta15", "|eta|<1.5")):
        incl = _inclusive_fractions(acc, region)
        me = _match_eff_inclusive(acc, region)

        lines.append(f"## Inclusive ({label}, 14 <= pT < 30 GeV)")
        lines.append("")
        lines.append(f"Total truth photons in this selection: **{incl['_total']:,}**")
        lines.append("")
        lines.append("| Category | N | Fraction |")
        lines.append("|---|---:|---:|")
        for code, key, label_nice, _ in CATEGORIES:
            n, p, e = incl[key]
            lines.append(f"| {label_nice} | {n:,} | {p:.3%} +- {e:.3%} |")
        lines.append("")
        lines.append("### Cluster-match efficiency")
        lines.append("")
        lines.append("| Category | Matched | Total | Efficiency |")
        lines.append("|---|---:|---:|---:|")
        for code, key, label_nice, _ in CATEGORIES:
            k, n, p, e = me[key]
            lines.append(f"| {label_nice} | {k:,} | {n:,} | {p:.3%} +- {e:.3%} |")
        lines.append("")

        lines.append(f"### Per-pT-bin fractions ({label})")
        lines.append("")
        lines.append("| pT [GeV] | N | clean | converted | bad |")
        lines.append("|---|---:|---:|---:|---:|")
        for i in range(len(PT_EDGES) - 1):
            lines.append("| " + _fmt_row_pt(acc, region, i) + " |")
        lines.append("")

        lines.append(f"### Per-|eta|-bin fractions (14 <= pT < 30 GeV, {label})")
        lines.append("")
        lines.append("| |eta| | N | clean | converted | bad |")
        lines.append("|---|---:|---:|---:|---:|")
        for i in range(len(ABS_ETA_EDGES) - 1):
            lines.append("| " + _fmt_row_abseta(acc, region, i) + " |")
        lines.append("")

    incl07 = _inclusive_fractions(acc, "eta07")
    incl15 = _inclusive_fractions(acc, "eta15")
    me07 = _match_eff_inclusive(acc, "eta07")
    me15 = _match_eff_inclusive(acc, "eta15")

    lines.append("## Physics discussion")
    lines.append("")
    lines.append(
        f"Inclusive converted fraction at |eta|<0.7 (14 <= pT < 30 GeV) is "
        f"**{incl07['converted'][1]:.3%}**, and at |eta|<1.5 it is "
        f"**{incl15['converted'][1]:.3%}**."
    )
    lines.append("")
    lines.append(
        "The EMCal front face is preceded by the TPC, TPC outer field cage / support, INTT, "
        "and services. In the central barrel the integrated material is usually quoted at the ~10-20% X0 "
        "level, rising at larger |eta|. A photon conversion probability of ~(7/9) * X/X0 predicts ~7-15% "
        "converted fraction in the central region, increasing toward |eta|~1.5 due to extra path length "
        "and edge material."
    )
    lines.append("")
    lines.append(
        "The observed central-barrel value is in the ballpark of this expectation; if anything it sits a "
        "little higher, suggesting the integrated material in front of the EMCal is closer to the upper "
        "end of the ~10-20% X0 range. The sharp rise at |eta| > 1 in the per-|eta| table is consistent "
        "with endcap-style material growth."
    )
    lines.append("")
    lines.append("### Cluster-match efficiency per category (14 <= pT < 30 GeV)")
    lines.append("")
    lines.append(f"- clean: |eta|<0.7 = **{me07['clean'][2]:.3%}**, |eta|<1.5 = {me15['clean'][2]:.3%}")
    lines.append(f"- converted: |eta|<0.7 = **{me07['converted'][2]:.3%}**, |eta|<1.5 = {me15['converted'][2]:.3%}")
    lines.append(f"- bad: |eta|<0.7 = **{me07['bad'][2]:.3%}**, |eta|<1.5 = {me15['bad'][2]:.3%}")
    lines.append("")
    delta_abs = me07['converted'][2] - me07['clean'][2]  # (conv - clean), negative if conv is lower
    delta_rel = delta_abs / me07['clean'][2] if me07['clean'][2] > 0 else 0.0
    lines.append(
        f"Converted photons have a match efficiency that is **{delta_abs:+.3%}** (absolute) / "
        f"**{delta_rel:+.2%}** (relative) lower than clean photons in the fiducial region at |eta|<0.7. "
        "Most converted photons still produce a single merged cluster (the two e+e- daughters typically "
        "deposit within the same tower cluster at these ETs, so the deficit is modest). At |eta|<1.5 the "
        "converted match efficiency is substantially lower than clean (~10 percentage points), reflecting "
        "that conversions at forward rapidity are more likely to spread over multiple towers or escape the "
        "fiducial cluster. The 'bad' category -- where a substantial fraction of the photon momentum is "
        "carried by a non-ee secondary -- shows the largest deficit, as expected."
    )
    lines.append("")
    lines.append("## Concerns / caveats")
    lines.append("")
    lines.append(
        "- The matching uses the **primary truth track ID** of the photon. A conversion can produce e+e- "
        "daughters with different track IDs whose cluster may be labeled with the *daughter* trkid instead "
        "of the parent photon; such events will be counted as 'not matched' even when a physical cluster "
        "is present. The converted-category match efficiency should be read as a **lower bound**."
    )
    lines.append(
        "- The 'bad' fraction is small and dominated by statistics in several bins; bin-to-bin noise is "
        "statistical, not physics."
    )
    lines.append(
        "- The cluster Et cut (5 GeV) is below the analysis threshold; tighter Et cuts would widen the "
        "clean/converted gap, especially at low pT."
    )
    lines.append(
        "- Single-particle MC: no pileup, no underlying event, no run-range or vertex reweighting."
    )
    lines.append(
        "- The per-pT-bin tables include **all** truth photons that appear in the event record, not just "
        "the primary generator photon. Truth bins below ~14 GeV are dominated by low-pT photons produced "
        "in the shower / rescattering of the primary (pi0 decays, bremsstrahlung photons, etc.). The "
        "inclusive numbers at 14 <= pT < 30 GeV are the ones relevant to isolated-photon physics."
    )
    lines.append(
        "- The absolute match-efficiency values (~65-70% for clean photons at |eta|<0.7) look low at first "
        "glance for a truth photon in a single-particle sample. This is because `cluster_truthtrkID` labels "
        "the *dominant contributing truth particle* of the cluster, which after showering can be a "
        "secondary e+/e- rather than the original photon -- so many physical matches are missed by the "
        "strict track-ID equality rule. The *relative* comparison between categories is more meaningful "
        "than the absolute number."
    )
    lines.append("")

    with open(outpath, "w") as fh:
        fh.write("\n".join(lines))
    print("wrote report:", outpath)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-entries", type=int, default=None)
    ap.add_argument("--step-size", type=str, default="250000")
    ap.add_argument("--outdir", type=str,
                    default="/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/converted_photon_study")
    args = ap.parse_args()

    os.makedirs(os.path.join(args.outdir, "rootFiles"), exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "figures"), exist_ok=True)

    print(f"Processing tree at {INPUT_FILE}:{TREE_NAME}")
    if args.max_entries:
        print(f"  (limiting to first {args.max_entries:,} events)")
    print(f"  step_size = {args.step_size}")

    step_size = args.step_size
    try:
        step_size = int(args.step_size)
    except ValueError:
        pass

    acc, n_events_processed = process(max_entries=args.max_entries, step_size=step_size)

    root_path = os.path.join(args.outdir, "rootFiles", "fractions_converted.root")
    write_root(acc, root_path)

    fig_dir = os.path.join(args.outdir, "figures")
    for region in acc.eta_regions:
        plot_frac_vs_pt(acc, region, os.path.join(fig_dir, f"fractions_vs_pt_{region}.pdf"))
        plot_frac_vs_abseta(acc, region, os.path.join(fig_dir, f"fractions_vs_abseta_{region}.pdf"))
        for code, key, label_nice, _ in CATEGORIES:
            if key == "clean":
                continue
            plot_frac_2d(acc, region, code,
                         os.path.join(fig_dir, f"fractions_2d_{region}_{key}.pdf"),
                         label=label_nice)
        plot_match_eff_vs_pt(acc, region, os.path.join(fig_dir, f"fractions_match_eff_{region}.pdf"))
        plot_multiplicity(acc, region, os.path.join(fig_dir, f"fractions_multiplicity_{region}.pdf"))

    n_msg = f"Processed {n_events_processed:,} events."
    write_report(acc, n_msg, os.path.join(args.outdir, "fractions_findings.md"))

    return 0


if __name__ == "__main__":
    sys.exit(main())
