#!/usr/bin/env python3
"""Converted vs unconverted photon shower-shape study.

Reads the photon10 MC file, matches reco clusters to truth photons by trkID,
and fills 1D histograms of shower-shape observables per (variable, pt_bin,
converted_category) category.

Inputs:
  /sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_with_bdt_vtx_noreweight_single.root
  tree = slimtree
  NODE_CEMC    = CLUSTERINFO_CEMC          (shower-shape + CNN_prob)
  NODE_NOSPLIT = CLUSTERINFO_CEMC_NO_SPLIT (BDT score branch)

Selection:
  truth: particle_pid==22, |particle_Eta|<0.7, particle_truth_iso_03<4,
         particle_Pt in [8, 36]
  reco: cluster_truthtrkID == particle_trkid AND cluster_Et > 5 GeV
        (pick highest-ET match per truth photon)

Categories:
  particle_converted  0 -> unconverted
                      1 -> converted (e+e-)
                      2 -> bad (non-e+e- secondary carrying >40% momentum)

pT binning (matched-cluster Et as x-axis):
  ptRanges = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]  (12 bins)

Run:
  python3 analyze_showershape.py [--max-events N]
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

HERE = Path(__file__).resolve().parent
OUTDIR = HERE / "rootFiles"
FIGDIR = HERE / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

INPUT_FILE = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_with_bdt_vtx_noreweight_single.root"
TREE = "slimtree"

NODE_CEMC = "CLUSTERINFO_CEMC"
NODE_NOSPLIT = "CLUSTERINFO_CEMC_NO_SPLIT"

# pT bins in analysis (using cluster_Et of matched cluster as x-axis)
PT_BINS = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
N_PT = len(PT_BINS) - 1  # 12

# Truth-level selection window
PT_LO_TRUTH, PT_HI_TRUTH = 8.0, 36.0
ETA_ABS_MAX = 0.7
ISO03_MAX = 4.0
PID_PHOTON = 22
MATCH_ET_MIN = 5.0  # GeV, minimum cluster ET to count as a match

# Categories: (name, flag value, display label)
CATEGORIES = [("unconv", 0), ("conv", 1), ("bad", 2)]
CAT_NAMES = [c[0] for c in CATEGORIES]

# ---------------------------------------------------------------------------
# Variables to study.  Each entry:
#   name    : key for histogram storage (ASCII, safe)
#   branch  : the branch on NODE_CEMC (without the _<NODE> suffix)
#   title   : x-axis title in LaTeX / ROOT style
#   nbins, xlo, xhi
#   derived : None for straight branch, else a tuple (num_branch, den_branch)
#             for simple ratios
# ---------------------------------------------------------------------------
VARIABLES = [
    # ----- shower widths -----
    {"name": "weta_cogx", "branch": "cluster_weta_cogx",
     "title": r"$w_{\eta,\,\mathrm{CoG\text{-}x}}$", "nb": 60, "lo": 0.0, "hi": 1.5},
    {"name": "wphi_cogx", "branch": "cluster_wphi_cogx",
     "title": r"$w_{\phi,\,\mathrm{CoG\text{-}x}}$", "nb": 60, "lo": 0.0, "hi": 1.5},
    # ----- ring ET sums (individual) -----
    {"name": "et1", "branch": "cluster_et1",
     "title": r"$E_{T,1}$ / cluster ET",  "nb": 60, "lo": 0.0, "hi": 1.2},
    {"name": "et2", "branch": "cluster_et2",
     "title": r"$E_{T,2}$ / cluster ET",  "nb": 60, "lo": 0.0, "hi": 1.2},
    {"name": "et3", "branch": "cluster_et3",
     "title": r"$E_{T,3}$ / cluster ET",  "nb": 60, "lo": 0.0, "hi": 1.2},
    {"name": "et4", "branch": "cluster_et4",
     "title": r"$E_{T,4}$ / cluster ET",  "nb": 60, "lo": 0.0, "hi": 1.2},
    # ----- CNN / topological probs -----
    {"name": "CNN_prob", "branch": "cluster_CNN_prob",
     "title": r"CNN photon prob.",        "nb": 60, "lo": 0.0, "hi": 1.0},
    {"name": "cluster_prob", "branch": "cluster_prob",
     "title": r"cluster prob. ($\chi^2$)",  "nb": 60, "lo": 0.0, "hi": 1.0},
    {"name": "merged_prob", "branch": "cluster_merged_prob",
     "title": r"cluster merged prob.",    "nb": 60, "lo": 0.0, "hi": 1.0},
]

# Derived ratios (computed from the loaded e## / et## branches)
RATIOS = [
    {"name": "e11_over_e33", "num": "cluster_e11", "den": "cluster_e33",
     "title": r"$E_{1\times1}\,/\,E_{3\times3}$", "nb": 60, "lo": 0.0, "hi": 1.1},
    {"name": "e32_over_e35", "num": "cluster_e32", "den": "cluster_e35",
     "title": r"$E_{3\times2}\,/\,E_{3\times5}$", "nb": 60, "lo": 0.0, "hi": 1.2},
    {"name": "et2_over_et1", "num": "cluster_et2", "den": "cluster_et1",
     "title": r"$E_{T,2}\,/\,E_{T,1}$", "nb": 60, "lo": 0.0, "hi": 1.5},
    {"name": "et3_over_et1", "num": "cluster_et3", "den": "cluster_et1",
     "title": r"$E_{T,3}\,/\,E_{T,1}$", "nb": 60, "lo": 0.0, "hi": 1.5},
    {"name": "et4_over_et1", "num": "cluster_et4", "den": "cluster_et1",
     "title": r"$E_{T,4}\,/\,E_{T,1}$", "nb": 60, "lo": 0.0, "hi": 1.5},
]

# BDT on NO_SPLIT node
BDT_VAR = {"name": "bdt_score_nosplit", "branch": "cluster_bdt",
           "title": r"BDT photon score (NO_SPLIT node)",
           "nb": 60, "lo": 0.0, "hi": 1.0}

# ---------------------------------------------------------------------------

def suf(b: str, node: str) -> str:
    return f"{b}_{node}"


# All CEMC branches we need to load (scalar features + ratio numerators/denoms)
SHOWERSHAPE_BRANCHES = sorted(set(
    [v["branch"] for v in VARIABLES] +
    [r["num"] for r in RATIOS] +
    [r["den"] for r in RATIOS]
))

BRANCHES_CEMC = [suf(b, NODE_CEMC) for b in SHOWERSHAPE_BRANCHES] + [
    suf("cluster_Et", NODE_CEMC),
    suf("cluster_Eta", NODE_CEMC),
    suf("cluster_truthtrkID", NODE_CEMC),
]

BRANCHES_NOSPLIT = [
    suf("cluster_Et", NODE_NOSPLIT),
    suf("cluster_bdt", NODE_NOSPLIT),
    suf("cluster_truthtrkID", NODE_NOSPLIT),
]

BRANCHES_TRUTH = [
    "nparticles",
    "particle_pid",
    "particle_trkid",
    "particle_Pt",
    "particle_Eta",
    "particle_truth_iso_03",
    "particle_converted",
]

ALL_BRANCHES = BRANCHES_TRUTH + BRANCHES_CEMC + BRANCHES_NOSPLIT


# ---------------------------------------------------------------------------

class HistStore:
    """1D histograms per (variable, pt_bin, category) with inclusive key."""

    def __init__(self):
        # Flat dict: key="{var}__{cat}__pt{N}" or "{var}__{cat}__incl"
        self.hist = {}
        self.edges = {}
        self.counts = {c: 0 for c in CAT_NAMES}  # total entries in selection per cat
        # Per-pT, per-cat counts for reporting
        self.count_bin = {c: np.zeros(N_PT, dtype=np.int64) for c in CAT_NAMES}
        # Sum / Sum-of-squares for means per variable per pT bin per category
        self.sum_v = {}
        self.sum_v2 = {}
        self.nentries_v = {}

    def ensure(self, var_name: str, nb: int, lo: float, hi: float):
        edges = np.linspace(lo, hi, nb + 1)
        for c in CAT_NAMES:
            for i in range(N_PT):
                k = f"{var_name}__{c}__pt{i}"
                if k not in self.hist:
                    self.hist[k] = np.zeros(nb, dtype=np.float64)
                    self.edges[k] = edges
            k = f"{var_name}__{c}__incl"
            if k not in self.hist:
                self.hist[k] = np.zeros(nb, dtype=np.float64)
                self.edges[k] = edges

        # stats
        if var_name not in self.sum_v:
            self.sum_v[var_name] = {c: np.zeros(N_PT, dtype=np.float64) for c in CAT_NAMES}
            self.sum_v2[var_name] = {c: np.zeros(N_PT, dtype=np.float64) for c in CAT_NAMES}
            self.nentries_v[var_name] = {c: np.zeros(N_PT, dtype=np.int64) for c in CAT_NAMES}

    def fill_per_var(self, var_name: str, nb: int, lo: float, hi: float,
                     values: np.ndarray, pt_idx: np.ndarray, conv: np.ndarray):
        """Fill (var, cat, ptbin) for all matched photons.

        values[i]: observable for photon i (nan or inf -> skipped)
        pt_idx[i]: index in [0, N_PT-1] or -1 if out of range
        conv[i]:   0/1/2
        """
        self.ensure(var_name, nb, lo, hi)

        valid = np.isfinite(values) & (pt_idx >= 0)
        if not np.any(valid):
            return

        v = values[valid]
        p = pt_idx[valid]
        c = conv[valid]

        # bin edges already stored
        edges = self.edges[f"{var_name}__{CAT_NAMES[0]}__pt0"]
        bin_idx = np.clip(((v - lo) / (hi - lo) * nb).astype(np.int64), 0, nb - 1)

        for cat_name, cat_val in CATEGORIES:
            mcat = c == cat_val
            if not np.any(mcat):
                continue
            v_c = v[mcat]
            b_c = bin_idx[mcat]
            p_c = p[mcat]

            # Fill per-pT
            for i_pt in range(N_PT):
                mpt = p_c == i_pt
                if not np.any(mpt):
                    continue
                bb = b_c[mpt]
                # np.add.at is ok here; or use bincount
                counts = np.bincount(bb, minlength=nb)
                self.hist[f"{var_name}__{cat_name}__pt{i_pt}"] += counts

                # stats
                vv = v_c[mpt]
                self.sum_v[var_name][cat_name][i_pt] += vv.sum()
                self.sum_v2[var_name][cat_name][i_pt] += (vv * vv).sum()
                self.nentries_v[var_name][cat_name][i_pt] += len(vv)

            # Fill inclusive
            counts_incl = np.bincount(b_c, minlength=nb)
            self.hist[f"{var_name}__{cat_name}__incl"] += counts_incl


# ---------------------------------------------------------------------------

def match_clusters_vec(c_trk_ak, c_ET_ak, photon_event_idx, photon_trkid,
                        match_et_min=MATCH_ET_MIN):
    """Fully vectorized cluster matcher.

    For each truth photon (same event as a set of clusters), we want the cluster
    with the largest ET whose cluster_truthtrkID == particle_trkid AND ET>threshold.

    Implementation: build one "pair" row per (photon, cluster-in-same-event)
    combination, mask by trk match & ET>thresh, then groupby-argmax on photon idx.

    c_trk_ak, c_ET_ak : awkward.Array (one sub-list per event)
    photon_event_idx  : np.int64 array of shape (n_photons,), batch-local event idx
    photon_trkid      : np.int64 array of shape (n_photons,)

    Returns
    -------
    matched_abs_idx : np.int64 array of shape (n_photons,) with index into the
        flattened cluster arrays; -1 if no match or ET<=match_et_min.
    offsets : cumulative cluster offsets, length = n_events+1 (useful for later gathers)
    """
    cnt = ak.to_numpy(ak.num(c_trk_ak)).astype(np.int64)
    offsets = np.concatenate(([0], np.cumsum(cnt))).astype(np.int64)
    trk_flat = ak.to_numpy(ak.flatten(c_trk_ak)).astype(np.int64)
    ET_flat = ak.to_numpy(ak.flatten(c_ET_ak)).astype(np.float64)

    n_phot = len(photon_trkid)
    matched_abs_idx = np.full(n_phot, -1, dtype=np.int64)
    if n_phot == 0 or trk_flat.size == 0:
        return matched_abs_idx, offsets

    # For each photon i, figure out how many clusters are in its event
    n_clus_per_photon = cnt[photon_event_idx]  # shape (n_phot,)
    # Start offset per photon into flat cluster arrays
    start_per_photon = offsets[photon_event_idx]

    # Build a "pair" array by repeating photon info for each cluster in event
    pair_photon_idx = np.repeat(np.arange(n_phot, dtype=np.int64), n_clus_per_photon)
    # Cluster-local index 0..cnt-1 within photon's event:
    # per-photon local offsets using np.concatenate(arange(n) for n in counts)
    # Fast way: cumulative minus first-element trick
    total_pairs = int(n_clus_per_photon.sum())
    if total_pairs == 0:
        return matched_abs_idx, offsets

    # Build per-pair local cluster index
    # pair_local_cluster = arange per photon
    # use np.arange(total_pairs) - repeated_start_per_pair
    starts_per_pair = np.repeat(start_per_photon, n_clus_per_photon)
    # global absolute cluster index for each pair
    # arange within each photon: obtain via np.arange(total_pairs) - cum_start
    # simpler: pair_abs_idx = starts_per_pair + (arange(total_pairs) - cum_start_per_photon)
    cum_start = np.concatenate(([0], np.cumsum(n_clus_per_photon))).astype(np.int64)
    local_idx = np.arange(total_pairs, dtype=np.int64) - np.repeat(cum_start[:-1], n_clus_per_photon)
    pair_abs_cluster_idx = starts_per_pair + local_idx

    # Get truth trkID per photon, expanded to each pair
    pair_photon_trkid = np.repeat(photon_trkid, n_clus_per_photon)
    pair_cluster_trk = trk_flat[pair_abs_cluster_idx]
    pair_cluster_ET = ET_flat[pair_abs_cluster_idx]

    # Pair mask: trk matches AND ET > threshold
    pair_mask = (pair_cluster_trk == pair_photon_trkid) & (pair_cluster_ET > match_et_min)
    if not pair_mask.any():
        return matched_abs_idx, offsets

    # For each photon, pick pair with max ET among pair_mask
    pair_ET_masked = np.where(pair_mask, pair_cluster_ET, -np.inf)
    # Sort by (photon_idx, -ET) and take first for each photon: vectorized via
    # segment-argmax. Use np.maximum.reduceat after sorting:
    # Easier: use bincount+argmax emulation via pandas-free groupby.
    # Approach: for each photon i, find index with max ET in its pairs.
    # Convert pair_photon_idx to contiguous 0..n_phot-1 already.
    # Segment boundaries are cum_start (size n_phot+1).

    # Use np.maximum.reduceat on the ET array sliced per photon:
    # If a photon has 0 pairs, cum_start[i] == cum_start[i+1]; reduceat would fail.
    # Handle that by iterating only over photons that have at least one pair mask true.
    # Build best pair ET per photon using reduceat on pair_ET_masked
    # reduceat needs non-empty segments. Replace zero-length with sentinel to skip.
    has_clus = cum_start[1:] > cum_start[:-1]  # shape (n_phot,)
    # For photons with zero clusters in their event, skip entirely (matched stays -1).
    # np.maximum.reduceat with segment boundaries cum_start[:-1] works, but for empty segments
    # it returns the boundary value. Easier: compute segment max via sorted approach.
    # Alternative simple approach: argmax per-segment via np.add.reduceat over a trick.

    # Simpler approach: since pairs are grouped contiguously by photon (we built them that way),
    # we can iterate only over photons that have any True in pair_mask.
    # Vectorized: for each photon, best cluster abs idx = argmax of pair_ET_masked within segment.
    # Use a fast per-segment argmax:
    # pair_ET_masked is (total_pairs,) and segments are cum_start[:-1], cum_start[1:].
    # Compute running max idx via np.maximum.accumulate + reset at segment boundaries.
    # That's complex. For simplicity and speed, since total_pairs is usually small
    # (n_phot * avg clusters ~ n_phot * 4-10), the per-photon loop over segments is fast.

    # Do per-photon loop over segments (total loop iter = n_phot only)
    for i in range(n_phot):
        a = cum_start[i]
        b = cum_start[i + 1]
        if b == a:
            continue
        seg = pair_ET_masked[a:b]
        k = int(np.argmax(seg))
        if not np.isfinite(seg[k]):
            continue
        matched_abs_idx[i] = pair_abs_cluster_idx[a + k]

    return matched_abs_idx, offsets


def process(max_events: int | None = None, step_size: int = 500_000):
    store = HistStore()

    # Ensure all variables (including BDT and ratios) have storage
    for v in VARIABLES:
        store.ensure(v["name"], v["nb"], v["lo"], v["hi"])
    for r in RATIOS:
        store.ensure(r["name"], r["nb"], r["lo"], r["hi"])
    store.ensure(BDT_VAR["name"], BDT_VAR["nb"], BDT_VAR["lo"], BDT_VAR["hi"])

    total_photons = 0  # truth photons entering selection
    total_matched = 0
    events_seen = 0

    with uproot.open(INPUT_FILE) as f:
        t = f[TREE]
        entries = t.num_entries
        if max_events is not None:
            entries = min(entries, max_events)
        print(f"[analyze] processing {entries:,} / {t.num_entries:,} events, step={step_size:,}")
        for batch in t.iterate(
            ALL_BRANCHES,
            step_size=step_size,
            entry_stop=entries,
        ):
            n_batch = len(batch)
            events_seen += n_batch

            pid = batch["particle_pid"]
            pt_truth = batch["particle_Pt"]
            eta_truth = batch["particle_Eta"]
            iso = batch["particle_truth_iso_03"]
            conv = batch["particle_converted"]
            trkid = batch["particle_trkid"]

            sel = (
                (pid == PID_PHOTON)
                & (pt_truth >= PT_LO_TRUTH)
                & (pt_truth <= PT_HI_TRUTH)
                & (np.abs(eta_truth) < ETA_ABS_MAX)
                & (iso < ISO03_MAX)
            )
            n_sel = int(ak.sum(sel))
            if n_sel == 0:
                continue
            total_photons += n_sel

            photon_pt = ak.flatten(pt_truth[sel]).to_numpy()
            photon_conv = ak.flatten(conv[sel]).to_numpy()
            photon_trkid = ak.flatten(trkid[sel]).to_numpy()
            counts_per_event = ak.to_numpy(ak.sum(sel, axis=1))
            photon_event_idx = np.repeat(np.arange(n_batch), counts_per_event)

            # --- CEMC node matching (vectorized via flattened arrays) ---
            c_ET_cemc = batch[suf("cluster_Et", NODE_CEMC)]
            c_trk_cemc = batch[suf("cluster_truthtrkID", NODE_CEMC)]
            matched_cemc, offsets_cemc = match_clusters_vec(
                c_trk_cemc, c_ET_cemc, photon_event_idx, photon_trkid
            )

            matched_mask = matched_cemc >= 0
            total_matched += int(matched_mask.sum())
            if not np.any(matched_mask):
                continue

            # Gather: take flattened branch + matched_abs_idx -> per-photon value
            def gather(branch: str):
                arr_flat = ak.to_numpy(ak.flatten(batch[suf(branch, NODE_CEMC)]))
                out = np.full(len(matched_cemc), np.nan, dtype=np.float64)
                m = matched_mask
                out[m] = arr_flat[matched_cemc[m]]
                return out

            # Compute matched cluster ET to drive pT binning
            cluster_ET = gather("cluster_Et")
            # pt_idx based on cluster_ET (BDT convention)
            pt_idx = np.digitize(cluster_ET, PT_BINS) - 1
            pt_idx = np.where((pt_idx < 0) | (pt_idx >= N_PT), -1, pt_idx).astype(np.int64)

            # per-cat counts for logging
            for cat_name, cat_val in CATEGORIES:
                mm = (photon_conv == cat_val) & matched_mask
                store.counts[cat_name] += int(mm.sum())
                # per pT bin
                mmpt = mm & (pt_idx >= 0)
                for i in range(N_PT):
                    store.count_bin[cat_name][i] += int(((pt_idx == i) & mm).sum())

            # ----- Straight-branch variables -----
            # Use masked fills: nan when not matched
            for v in VARIABLES:
                values = gather(v["branch"])
                # Only matched ones are valid
                values[~matched_mask] = np.nan
                store.fill_per_var(v["name"], v["nb"], v["lo"], v["hi"],
                                   values, pt_idx, photon_conv)

            # ----- Ratios -----
            num_cache = {}
            for r in RATIOS:
                if r["num"] not in num_cache:
                    num_cache[r["num"]] = gather(r["num"])
                if r["den"] not in num_cache:
                    num_cache[r["den"]] = gather(r["den"])
                num = num_cache[r["num"]]
                den = num_cache[r["den"]]
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = np.where(den > 0, num / den, np.nan)
                ratio[~matched_mask] = np.nan
                store.fill_per_var(r["name"], r["nb"], r["lo"], r["hi"],
                                   ratio, pt_idx, photon_conv)

            # ----- BDT on NO_SPLIT node (rematch on that node, vectorized) -----
            c_ET_ns = batch[suf("cluster_Et", NODE_NOSPLIT)]
            c_trk_ns = batch[suf("cluster_truthtrkID", NODE_NOSPLIT)]
            matched_ns, _ = match_clusters_vec(
                c_trk_ns, c_ET_ns, photon_event_idx, photon_trkid
            )

            bdt_flat = ak.to_numpy(ak.flatten(batch[suf("cluster_bdt", NODE_NOSPLIT)]))
            ET_ns_flat = ak.to_numpy(ak.flatten(c_ET_ns))

            bdt_vals = np.full(len(matched_ns), np.nan, dtype=np.float64)
            cluster_ET_ns = np.full(len(matched_ns), np.nan, dtype=np.float64)
            m_ns = matched_ns >= 0
            if np.any(m_ns):
                bdt_vals[m_ns] = bdt_flat[matched_ns[m_ns]]
                cluster_ET_ns[m_ns] = ET_ns_flat[matched_ns[m_ns]]

            pt_idx_ns = np.digitize(cluster_ET_ns, PT_BINS) - 1
            pt_idx_ns = np.where((pt_idx_ns < 0) | (pt_idx_ns >= N_PT), -1, pt_idx_ns).astype(np.int64)

            store.fill_per_var(BDT_VAR["name"], BDT_VAR["nb"], BDT_VAR["lo"], BDT_VAR["hi"],
                               bdt_vals, pt_idx_ns, photon_conv)

            if events_seen % (5 * step_size) < step_size:
                print(f"  events={events_seen:,}  truth={total_photons:,}  matched={total_matched:,}")

    print(f"[analyze] DONE  events={events_seen:,}  truth={total_photons:,}  matched={total_matched:,}")
    for cat in CAT_NAMES:
        print(f"    cat={cat:>6s}: matched in selection = {store.counts[cat]:,}")
    return store, events_seen, total_photons, total_matched


# ---------------------------------------------------------------------------
def save_root(store: HistStore, events_seen: int, total_photons: int,
              total_matched: int, outpath: Path):
    """Write TH1D per (variable, cat, pt_bin) to ROOT using uproot.writing."""
    out_trees = {}
    # TH1 per key
    import uproot as up
    with up.recreate(outpath) as fout:
        for k, h in store.hist.items():
            edges = store.edges[k]
            fout[k] = (h, edges)
        # Stats
        meta_keys = []
        for v in VARIABLES + RATIOS + [BDT_VAR]:
            name = v["name"]
            if name not in store.sum_v:
                continue
            for cat in CAT_NAMES:
                means = np.where(store.nentries_v[name][cat] > 0,
                                 store.sum_v[name][cat] / np.maximum(store.nentries_v[name][cat], 1),
                                 np.nan)
                var = np.where(store.nentries_v[name][cat] > 0,
                               store.sum_v2[name][cat] / np.maximum(store.nentries_v[name][cat], 1) - means ** 2,
                               np.nan)
                rms = np.sqrt(np.maximum(var, 0))
                # Write TH1D-like of means with dummy edges = PT_BINS
                fout[f"mean__{name}__{cat}"] = (means, PT_BINS)
                fout[f"rms__{name}__{cat}"] = (rms, PT_BINS)
                fout[f"nent__{name}__{cat}"] = (store.nentries_v[name][cat].astype(np.float64), PT_BINS)
                meta_keys.append(f"mean__{name}__{cat}")
    print(f"[save] wrote {outpath}  (with {len(store.hist)} hists + {len(meta_keys)} mean graphs)")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def _legend_counts(store: HistStore, cat: str, pt_bin: int | None) -> int:
    if pt_bin is None:
        return int(sum(store.count_bin[cat]))
    return int(store.count_bin[cat][pt_bin])


MIN_CAT_ENTRIES = 200  # suppress a category on a given panel below this count
                        # (bad photons are ~280 inclusive across all pT bins,
                        # and tiny per pT bin, so they get dropped from panels;
                        # this avoids noisy spikes. Counts in findings.md document
                        # the actual bad-photon yield).


def plot_variable(store: HistStore, var_info: dict, outdir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    name = var_info["name"]
    nb, lo, hi = var_info["nb"], var_info["lo"], var_info["hi"]
    edges = np.linspace(lo, hi, nb + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = edges[1:] - edges[:-1]

    style = {
        "unconv": dict(color="black", marker="o", fillstyle="full", label="Unconverted"),
        "conv":   dict(color="red",   marker="s", fillstyle="none", label="Converted"),
        "bad":    dict(color="blue",  marker="^", fillstyle="none", label="Bad photon"),
    }

    # Multi-panel: 4x3 = 12 pT bins
    fig, axes = plt.subplots(4, 3, figsize=(12, 14), sharex=False, sharey=False)
    axes = axes.flatten()
    for i_pt in range(N_PT):
        ax = axes[i_pt]
        ymax_guide = 0.0
        for cat_name, _ in CATEGORIES:
            h = store.hist.get(f"{name}__{cat_name}__pt{i_pt}")
            if h is None or h.sum() < MIN_CAT_ENTRIES:
                continue
            norm = h.sum()
            y = h / norm / widths
            n_in_bin = _legend_counts(store, cat_name, i_pt)
            label = f"{style[cat_name]['label']} (N={n_in_bin:,})"
            ax.step(centers, y, where="mid", color=style[cat_name]["color"], label=label,
                    linewidth=1.3)
            ax.errorbar(centers, y,
                        yerr=np.sqrt(np.maximum(h, 0)) / norm / widths,
                        fmt=style[cat_name]["marker"],
                        color=style[cat_name]["color"],
                        mfc=style[cat_name]["color"] if style[cat_name]["fillstyle"] == "full" else "white",
                        mec=style[cat_name]["color"],
                        markersize=3.5, linestyle="", capsize=0, linewidth=0)
            # Track y-scale from unconv+conv only (bad tends to spike)
            if cat_name != "bad":
                ymax_guide = max(ymax_guide, float(y.max()))
        ax.set_xlabel(var_info["title"])
        ax.set_ylabel("Normalized")
        pt_lo, pt_hi = PT_BINS[i_pt], PT_BINS[i_pt + 1]
        ax.set_title(fr"${pt_lo:.0f} \leq p_T < {pt_hi:.0f}$ GeV", fontsize=10)
        ax.grid(alpha=0.3)
        if ymax_guide > 0:
            ax.set_ylim(0, 1.3 * ymax_guide)
            ax.legend(loc="best", fontsize=8)
        # sPHENIX label (top-left small)
        ax.text(0.03, 0.97, r"$\bf{\it{sPHENIX}}$ Simulation",
                transform=ax.transAxes, ha="left", va="top", fontsize=8)
    fig.suptitle(f"{var_info['title']}  (matched cluster, converted-photon categories)",
                 fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    out = outdir / f"showershape_{name}.pdf"
    fig.savefig(out)
    plt.close(fig)

    # Inclusive canvas (all pT in [8,36])
    fig2, ax = plt.subplots(figsize=(7, 5))
    ymax_guide = 0.0
    for cat_name, _ in CATEGORIES:
        h = store.hist.get(f"{name}__{cat_name}__incl")
        if h is None or h.sum() < MIN_CAT_ENTRIES:
            continue
        norm = h.sum()
        y = h / norm / widths
        n_in = _legend_counts(store, cat_name, None)
        label = f"{style[cat_name]['label']} (N={n_in:,})"
        ax.step(centers, y, where="mid", color=style[cat_name]["color"], label=label,
                linewidth=1.3)
        ax.errorbar(centers, y,
                    yerr=np.sqrt(np.maximum(h, 0)) / norm / widths,
                    fmt=style[cat_name]["marker"],
                    color=style[cat_name]["color"],
                    mfc=style[cat_name]["color"] if style[cat_name]["fillstyle"] == "full" else "white",
                    mec=style[cat_name]["color"],
                    markersize=4, linestyle="", capsize=0)
        if cat_name != "bad":
            ymax_guide = max(ymax_guide, float(y.max()))
    ax.set_xlabel(var_info["title"])
    ax.set_ylabel("Normalized")
    ax.set_title(fr"Photon10 MC, $8 \leq p_T < 36$ GeV (matched clusters)")
    ax.grid(alpha=0.3)
    use_log_y_incl = any(key in name for key in ("bdt_score", "CNN_prob", "cluster_prob", "merged_prob"))
    if use_log_y_incl and ymax_guide > 0:
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 5 * ymax_guide)
        ax.legend(loc="best", fontsize=9)
    elif ymax_guide > 0:
        ax.set_ylim(0, 1.3 * ymax_guide)
        ax.legend(loc="best", fontsize=9)
    ax.text(0.03, 0.97, r"$\bf{\it{sPHENIX}}$ Simulation",
            transform=ax.transAxes, ha="left", va="top", fontsize=9)
    out2 = outdir / f"showershape_{name}_inclusive.pdf"
    fig2.tight_layout()
    fig2.savefig(out2)
    plt.close(fig2)
    print(f"  [plot] wrote {out} and {out2}")


def plot_mean_vs_pt(store: HistStore, var_info: dict, outdir: Path):
    """Plot mean vs pT for a given variable, one line per category."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    name = var_info["name"]
    fig, ax = plt.subplots(figsize=(7, 5))
    pt_centers = 0.5 * (PT_BINS[:-1] + PT_BINS[1:])
    style = {
        "unconv": dict(color="black", marker="o", fillstyle="full", label="Unconverted"),
        "conv":   dict(color="red",   marker="s", fillstyle="none", label="Converted"),
        "bad":    dict(color="blue",  marker="^", fillstyle="none", label="Bad photon"),
    }
    for cat_name, _ in CATEGORIES:
        if cat_name == "bad":
            continue  # negligible population; suppress from summary plots
        n = store.nentries_v[name][cat_name]
        if n.sum() == 0:
            continue
        mean = np.where(n > 0, store.sum_v[name][cat_name] / np.maximum(n, 1), np.nan)
        mean2 = np.where(n > 0, store.sum_v2[name][cat_name] / np.maximum(n, 1), np.nan)
        var = np.maximum(mean2 - mean ** 2, 0)
        err = np.where(n > 0, np.sqrt(var) / np.sqrt(np.maximum(n, 1)), np.nan)
        total = int(n.sum())
        ax.errorbar(pt_centers, mean, yerr=err,
                    fmt=style[cat_name]["marker"],
                    color=style[cat_name]["color"],
                    mfc=style[cat_name]["color"] if style[cat_name]["fillstyle"] == "full" else "white",
                    mec=style[cat_name]["color"],
                    label=f"{style[cat_name]['label']} (N={total:,})",
                    markersize=6, capsize=2, linewidth=1.3)
    ax.set_xlabel(r"matched cluster $E_T$ [GeV]")
    ax.set_ylabel(f"mean {var_info['title']}")
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    ax.text(0.03, 0.97, r"$\bf{\it{sPHENIX}}$ Simulation",
            transform=ax.transAxes, ha="left", va="top", fontsize=9)
    ax.set_title(f"Mean of {var_info['title']} vs $p_T$")
    fig.tight_layout()
    out = outdir / f"showershape_{name}_mean.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"  [plot] wrote {out}")


# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-events", type=int, default=None,
                        help="Cap on number of events to process (default: full file)")
    parser.add_argument("--step-size", type=int, default=500_000)
    parser.add_argument("--plots-only", action="store_true",
                        help="Skip processing; read histograms from ROOT and replot")
    args = parser.parse_args()

    outroot = OUTDIR / "showershape_converted.root"

    if args.plots_only:
        print("[main] --plots-only: reading back histograms from", outroot)
        import uproot as up
        store = HistStore()
        # Rebuild skeleton
        for v in VARIABLES + RATIOS + [BDT_VAR]:
            store.ensure(v["name"], v["nb"], v["lo"], v["hi"])
        with up.open(outroot) as fin:
            keys = fin.keys(filter_classname="TH1*", cycle=False)
            for k in keys:
                kn = k.rsplit(";", 1)[0]
                arr = fin[k].values()
                if kn in store.hist:
                    store.hist[kn] = arr.astype(np.float64)
            for v in VARIABLES + RATIOS + [BDT_VAR]:
                name = v["name"]
                for cat in CAT_NAMES:
                    try:
                        mean = fin[f"mean__{name}__{cat}"].values()
                        nent = fin[f"nent__{name}__{cat}"].values()
                        # Reconstruct sum_v = mean * nent (approx)
                        store.sum_v[name][cat] = mean * nent
                        store.nentries_v[name][cat] = nent.astype(np.int64)
                        rms = fin[f"rms__{name}__{cat}"].values()
                        store.sum_v2[name][cat] = (rms ** 2 + mean ** 2) * nent
                    except Exception:
                        pass
        # Don't know per-bin counts precisely from ROOT; estimate from a reference var (first VAR)
        ref = VARIABLES[0]["name"]
        for cat in CAT_NAMES:
            for i in range(N_PT):
                k = f"{ref}__{cat}__pt{i}"
                if k in store.hist:
                    store.count_bin[cat][i] = int(store.hist[k].sum())
            store.counts[cat] = int(sum(store.count_bin[cat]))
    else:
        store, events_seen, total_photons, total_matched = process(
            max_events=args.max_events, step_size=args.step_size
        )
        save_root(store, events_seen, total_photons, total_matched, outroot)

    # Plots
    print("[main] making plots...")
    all_vars = VARIABLES + RATIOS + [BDT_VAR]
    for v in all_vars:
        plot_variable(store, v, FIGDIR)
        plot_mean_vs_pt(store, v, FIGDIR)

    # Final counts summary (nice for terminal)
    print("\n[summary] matched photons per category:")
    for cat_name, _ in CATEGORIES:
        print(f"  {cat_name}: total={store.counts[cat_name]:,}")
        for i in range(N_PT):
            print(f"     pT[{PT_BINS[i]:.0f}-{PT_BINS[i+1]:.0f}]: {store.count_bin[cat_name][i]:,}")


if __name__ == "__main__":
    main()
