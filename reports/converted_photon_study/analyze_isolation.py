#!/usr/bin/env python3
"""Analyze isolation-energy distributions for converted vs unconverted photons
in the PPG12 single-particle photon10 MC sample.

Reads photon10_with_bdt_vtx_noreweight_single.root with uproot.iterate,
matches truth photons to clusters, categorizes by particle_converted, and
fills 1D iso histograms in pT bins. Saves to ROOT + matplotlib PDFs.
"""

from __future__ import annotations

import os
import sys
import json
import argparse
from collections import defaultdict

import numpy as np
import uproot

# Headless ROOT / matplotlib
os.environ.setdefault("DISPLAY", "")
import ROOT  # noqa: E402

ROOT.gROOT.SetBatch(True)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
FILE_PATH = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_with_bdt_vtx_noreweight_single.root"
TREE_NAME = "slimtree"
OUT_DIR = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/converted_photon_study"
ROOT_OUT = os.path.join(OUT_DIR, "rootFiles", "isolation_converted.root")
FIG_DIR = os.path.join(OUT_DIR, "figures")

# Canonical analysis pT edges (GeV)
PT_EDGES = np.array(
    [8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 32.0, 36.0],
    dtype=np.float64,
)
INCLUSIVE_LABEL = "8to36"

# Parametric iso cut from efficiencytool/config_bdt_nom.yaml
RECO_ISO_MAX_B = 0.502095
RECO_ISO_MAX_S = 0.0433036

# Conversion categories
CATEGORIES = {0: "unconverted", 1: "converted", 2: "bad"}
CAT_ORDER = [0, 1, 2]

# Iso branch list (base names; full branch name has _CLUSTERINFO_CEMC suffix)
ISO_VARS_R03 = [
    "cluster_iso_03",
    "cluster_iso_03_emcal",
    "cluster_iso_03_hcalin",
    "cluster_iso_03_hcalout",
]
ISO_VARS_R04 = [
    "cluster_iso_04",
    "cluster_iso_04_emcal",
    "cluster_iso_04_hcalin",
    "cluster_iso_04_hcalout",
]
ISO_VARS = ISO_VARS_R03 + ISO_VARS_R04
ISO_XLABEL = {
    "cluster_iso_03": "E_{T}^{iso}(R=0.3) [GeV]",
    "cluster_iso_03_emcal": "E_{T}^{iso, EMCal}(R=0.3) [GeV]",
    "cluster_iso_03_hcalin": "E_{T}^{iso, HCalIn}(R=0.3) [GeV]",
    "cluster_iso_03_hcalout": "E_{T}^{iso, HCalOut}(R=0.3) [GeV]",
    "cluster_iso_04": "E_{T}^{iso}(R=0.4) [GeV]",
    "cluster_iso_04_emcal": "E_{T}^{iso, EMCal}(R=0.4) [GeV]",
    "cluster_iso_04_hcalin": "E_{T}^{iso, HCalIn}(R=0.4) [GeV]",
    "cluster_iso_04_hcalout": "E_{T}^{iso, HCalOut}(R=0.4) [GeV]",
}

# Histogram ranges
HIST_RANGE_TOTAL = (-5.0, 15.0)
HIST_RANGE_SUB = (-5.0, 10.0)
HIST_NBINS = 100


def pt_bin_label(lo: float, hi: float) -> str:
    """Canonical label for pT bins (e.g. '8to10')."""
    return f"{int(round(lo))}to{int(round(hi))}"


def make_hist_name(var: str, cat_name: str, pt_label: str) -> str:
    return f"h_{var}_{cat_name}_pt{pt_label}"


# ---------------------------------------------------------------------------
# Selection kernel on a single batch
# ---------------------------------------------------------------------------
def process_batch(arr: dict, hists: dict, stats: dict) -> int:
    """Process one batch of arrays from uproot.iterate.

    Args:
        arr: dict of jagged arrays from uproot (awkward-like)
        hists: nested dict hists[var][cat][pt_label] -> (values accumulator)
        stats: running stats per (cat, pt_label)

    Returns:
        number of selected truth-photon/matched-cluster pairs in this batch
    """
    import awkward as ak

    particle_pid = arr["particle_pid"]
    particle_trkid = arr["particle_trkid"]
    particle_Eta = arr["particle_Eta"]
    particle_Pt = arr["particle_Pt"]  # unused but kept for completeness
    particle_truth_iso_03 = arr["particle_truth_iso_03"]
    particle_converted = arr["particle_converted"]

    cluster_E = arr["cluster_E_CLUSTERINFO_CEMC"]  # unused
    cluster_Et = arr["cluster_Et_CLUSTERINFO_CEMC"]
    cluster_Eta = arr["cluster_Eta_CLUSTERINFO_CEMC"]  # unused
    cluster_truthtrkID = arr["cluster_truthtrkID_CLUSTERINFO_CEMC"]

    iso_arrays = {
        var: arr[f"{var}_CLUSTERINFO_CEMC"] for var in ISO_VARS
    }

    # Truth photon selection
    truth_sel = (
        (particle_pid == 22)
        & (np.abs(particle_Eta) < 0.7)
        & (particle_truth_iso_03 < 4.0)
    )
    # No direct truth pT in selection logic (per task: categorize on cluster_Et)

    # Flatten into per-(event, particle) selection
    sel_particles = truth_sel
    n_sel_truth = int(ak.sum(sel_particles))
    if n_sel_truth == 0:
        return 0

    # For each selected truth photon, find the matched cluster with highest cluster_Et
    # whose cluster_truthtrkID matches particle_trkid.
    # Because events have >=1 truth photon and several clusters, we need a
    # per-event cross-product. We handle this by iterating event-by-event in a
    # vectorized fashion per event using awkward.

    # Keep only events with >=1 selected truth photon
    sel_trkid = particle_trkid[sel_particles]
    sel_conv = particle_converted[sel_particles]

    # For matching we need to broadcast: for each sel-particle, find best cluster in event.
    # Approach: build (event, cluster) x (event, particle) cross match using cartesian.

    # Keep per-event parallel lists of selected (trkid, converted).
    # Now cross-match against clusters in that event.

    # Cartesian per event between selected truth trkids and clusters
    # axis=1 means within each event
    pairs = ak.cartesian(
        {
            "trkid": sel_trkid,
            "conv": sel_conv,
            "cluster_Et": cluster_Et,
            "cluster_id": cluster_truthtrkID,
            **{f"iso_{v}": iso_arrays[v] for v in ISO_VARS},
        },
        axis=1,
        nested=False,
    )

    # We need to match trkid vs cluster_id and pick highest cluster_Et.
    # But ak.cartesian above gives a flat cross product per event; we need a
    # different construction so each selected truth has ALL clusters for the
    # same event. Let me use ak.argcartesian instead.
    # Simpler: iterate using ak.broadcast_arrays per event.

    # Better approach: for each (event, selected_particle), filter clusters with
    # cluster_truthtrkID == trkid, then pick argmax of cluster_Et.

    # Build: for each event, for each selected particle, list of matched clusters.
    # Use ak.cartesian with nested=True? Then the pair record will have arrays.

    # Construct the match-mask using a cross-join:
    # shape: per event, list[particle] of list[cluster]
    cross = ak.cartesian(
        [sel_trkid, cluster_truthtrkID], axis=1, nested=True
    )
    # cross[event][ipart][icluster].0 is sel trkid, .1 is cluster trkid
    match_mask = cross["0"] == cross["1"]

    # Build per-(event, particle) list of cluster indices that match.
    # We'll take a local index array for clusters and mask it.
    cluster_idx = ak.local_index(cluster_truthtrkID, axis=1)
    # Broadcast cluster_idx so it aligns with the nested cartesian shape
    # (per event, per particle, per cluster).
    _, cluster_idx_nested = ak.broadcast_arrays(sel_trkid[..., None], cluster_idx[:, None, :])

    # Broadcast cluster_Et similarly to pick the argmax among matched clusters.
    _, cluster_Et_nested = ak.broadcast_arrays(sel_trkid[..., None], cluster_Et[:, None, :])

    # Apply match mask; set non-matching to -inf so they lose argmax
    cluster_Et_matched = ak.where(match_mask, cluster_Et_nested, -np.inf)

    # For each (event, particle), argmax along the last axis
    # But if NO cluster matched, argmax returns None (option type). Fill with -1
    # and track with any_match to filter those out.
    any_match = ak.any(match_mask, axis=-1)  # shape: (event, particle)
    best_idx_local = ak.argmax(cluster_Et_matched, axis=-1, keepdims=False)  # (event, particle)
    best_idx_safe = ak.fill_none(best_idx_local, -1)

    # Flatten per-event lists for simpler indexing
    nev = len(cluster_truthtrkID)

    # Flatten events x particles
    n_part_per_event = ak.num(sel_trkid, axis=1)
    event_idx_per_part = np.repeat(np.arange(nev), np.asarray(n_part_per_event))
    flat_any_match = np.asarray(ak.flatten(any_match, axis=1))
    flat_best_idx = np.asarray(ak.flatten(best_idx_safe, axis=1)).astype(np.int64)
    flat_conv = np.asarray(ak.flatten(sel_conv, axis=1)).astype(np.int64)

    # Convert cluster arrays to awkward -> numpy via flat + offsets
    cluster_Et_arr = ak.to_numpy(ak.flatten(cluster_Et, axis=1))
    offsets = np.asarray(ak.num(cluster_Et, axis=1))
    cluster_starts = np.concatenate([[0], np.cumsum(offsets)])[:-1]

    # For each selected particle, flat cluster index = event_start + local_idx
    mask_valid = flat_any_match
    evt_idx_v = event_idx_per_part[mask_valid]
    loc_idx_v = flat_best_idx[mask_valid]
    flat_cluster_pos = cluster_starts[evt_idx_v] + loc_idx_v

    best_cluster_Et_v = cluster_Et_arr[flat_cluster_pos]

    # Build a flat dict of iso values for the matched cluster
    iso_flat_arrays = {}
    for var in ISO_VARS:
        arr_np = ak.to_numpy(ak.flatten(iso_arrays[var], axis=1))
        iso_flat_arrays[var] = arr_np[flat_cluster_pos]

    conv_v = flat_conv[mask_valid]

    # Require cluster_Et > 5 GeV
    et_cut = best_cluster_Et_v > 5.0
    evt_idx_v = evt_idx_v[et_cut]
    best_Et = best_cluster_Et_v[et_cut]
    conv_v = conv_v[et_cut]
    for var in ISO_VARS:
        iso_flat_arrays[var] = iso_flat_arrays[var][et_cut]

    n_accepted = int(best_Et.size)

    if n_accepted == 0:
        return 0

    # Fill histograms (per category x pT-bin)
    pt_bin_idx = np.digitize(best_Et, PT_EDGES) - 1  # 0..Nbins-1; outside -> -1 or Nbins
    n_pt_bins = len(PT_EDGES) - 1
    in_range = (pt_bin_idx >= 0) & (pt_bin_idx < n_pt_bins) & (best_Et < PT_EDGES[-1])

    for cat in CAT_ORDER:
        cat_name = CATEGORIES[cat]
        cat_mask = conv_v == cat

        # Inclusive (8<=Et<36)
        incl_mask = cat_mask & in_range
        stats[(cat, INCLUSIVE_LABEL)]["n"] += int(incl_mask.sum())
        for var in ISO_VARS:
            vals = iso_flat_arrays[var][incl_mask]
            if vals.size:
                hists[var][cat_name][INCLUSIVE_LABEL].append(vals)
                stats[(cat, INCLUSIVE_LABEL)][f"sum_{var}"] += float(vals.sum())

        # Per pT bin
        for ib in range(n_pt_bins):
            bin_mask = cat_mask & (pt_bin_idx == ib) & in_range
            if not bin_mask.any():
                continue
            lbl = pt_bin_label(PT_EDGES[ib], PT_EDGES[ib + 1])
            stats[(cat, lbl)]["n"] += int(bin_mask.sum())
            for var in ISO_VARS:
                vals = iso_flat_arrays[var][bin_mask]
                hists[var][cat_name][lbl].append(vals)
                stats[(cat, lbl)][f"sum_{var}"] += float(vals.sum())

    return n_accepted


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------
def run_analysis(max_events: int | None) -> tuple[dict, dict]:
    """Iterate over the file and return (hists, stats).

    hists[var][cat_name][pt_label] = list of np.ndarray chunks
    stats[(cat, pt_label)] = dict with 'n', 'sum_<var>'
    """
    # Initialise structures
    hists = {
        var: {CATEGORIES[c]: defaultdict(list) for c in CAT_ORDER}
        for var in ISO_VARS
    }
    pt_labels = [pt_bin_label(PT_EDGES[i], PT_EDGES[i + 1]) for i in range(len(PT_EDGES) - 1)]
    pt_labels.append(INCLUSIVE_LABEL)
    stats = {
        (c, lbl): {"n": 0, **{f"sum_{v}": 0.0 for v in ISO_VARS}}
        for c in CAT_ORDER for lbl in pt_labels
    }

    # Branches to read (keep minimal)
    branches = [
        "particle_pid",
        "particle_trkid",
        "particle_Pt",
        "particle_Eta",
        "particle_truth_iso_03",
        "particle_converted",
        "cluster_E_CLUSTERINFO_CEMC",
        "cluster_Et_CLUSTERINFO_CEMC",
        "cluster_Eta_CLUSTERINFO_CEMC",
        "cluster_truthtrkID_CLUSTERINFO_CEMC",
    ]
    for v in ISO_VARS:
        branches.append(f"{v}_CLUSTERINFO_CEMC")

    total_events = 0
    total_accepted = 0
    step = 500_000

    f = uproot.open(FILE_PATH)
    t = f[TREE_NAME]
    n_total = t.num_entries
    print(f"[info] tree has {n_total:,} events; max_events={max_events}", flush=True)

    entry_stop = n_total if max_events is None else min(max_events, n_total)

    for arr in t.iterate(
        branches,
        step_size=step,
        library="ak",
        entry_stop=entry_stop,
    ):
        n_batch = int(len(arr))
        n_acc = process_batch(arr, hists, stats)
        total_events += n_batch
        total_accepted += n_acc
        print(
            f"  [batch] events={n_batch:,} accepted={n_acc:,} "
            f"running_events={total_events:,} running_accepted={total_accepted:,}",
            flush=True,
        )

    print(f"[info] total processed events = {total_events:,}", flush=True)
    print(f"[info] total accepted truth-photon+cluster matches = {total_accepted:,}", flush=True)
    return hists, stats


# ---------------------------------------------------------------------------
# Histogram writing (ROOT)
# ---------------------------------------------------------------------------
def write_root_file(hists: dict) -> None:
    fout = ROOT.TFile(ROOT_OUT, "RECREATE")
    fout.cd()
    n_pt_bins = len(PT_EDGES) - 1
    labels = [pt_bin_label(PT_EDGES[i], PT_EDGES[i + 1]) for i in range(n_pt_bins)]
    labels.append(INCLUSIVE_LABEL)

    for var in ISO_VARS:
        is_sub = not var.endswith("iso_03") and not var.endswith("iso_04")
        rng = HIST_RANGE_SUB if is_sub else HIST_RANGE_TOTAL
        for cat in CAT_ORDER:
            cat_name = CATEGORIES[cat]
            for lbl in labels:
                name = make_hist_name(var, cat_name, lbl)
                h = ROOT.TH1D(name, name, HIST_NBINS, rng[0], rng[1])
                chunks = hists[var][cat_name].get(lbl, [])
                for ch in chunks:
                    for v in ch:
                        h.Fill(float(v))
                h.SetDirectory(fout)
                h.Write()
    fout.Close()
    print(f"[info] wrote ROOT histograms to {ROOT_OUT}", flush=True)


# ---------------------------------------------------------------------------
# Per-category summary statistics
# ---------------------------------------------------------------------------
def compute_summary(hists: dict, stats: dict) -> dict:
    """Compute per-category means and parametric iso pass fraction.

    Returns a nested dict suitable for JSON serialization.
    """
    summary = {}
    n_pt_bins = len(PT_EDGES) - 1
    pt_labels = [pt_bin_label(PT_EDGES[i], PT_EDGES[i + 1]) for i in range(n_pt_bins)]
    pt_labels.append(INCLUSIVE_LABEL)

    for cat in CAT_ORDER:
        cat_name = CATEGORIES[cat]
        summary[cat_name] = {}
        for lbl in pt_labels:
            # Concatenate cluster_iso_03 chunks for this cat, pt-bin
            chunks_iso03 = hists["cluster_iso_03"][cat_name].get(lbl, [])
            vals_iso03 = (
                np.concatenate(chunks_iso03) if chunks_iso03 else np.array([], dtype=np.float32)
            )
            n = int(vals_iso03.size)
            # For the cut, we need cluster_Et per entry — but we didn't keep it!
            # We only kept flat iso arrays. Workaround: store per-pT-bin — then
            # use bin centre for the ET value to compute cut threshold.
            # Better: use the stored pT-bin label midpoint OR compute an
            # ET-integrated cut: cut = b + s * <Et_bin>.
            if lbl == INCLUSIVE_LABEL:
                # inclusive: approximate <ET> via histogram? we didn't save ET.
                # Fall back to pT bin centres weighted by counts.
                # Here, just use a single representative ET = 14 GeV (median).
                # Better: store average ET per category during the loop.
                et_rep = 14.0
            else:
                lo, hi = [float(x) for x in lbl.split("to")]
                et_rep = 0.5 * (lo + hi)
            iso_cut = RECO_ISO_MAX_B + RECO_ISO_MAX_S * et_rep
            n_pass = int(np.sum(vals_iso03 < iso_cut)) if n > 0 else 0

            entry = {
                "n": n,
                "iso_cut": iso_cut,
                "n_pass_iso03": n_pass,
                "pass_frac_iso03": (n_pass / n) if n > 0 else np.nan,
            }
            # Means of each iso variable
            for var in ISO_VARS:
                chunks = hists[var][cat_name].get(lbl, [])
                if not chunks:
                    entry[f"mean_{var}"] = np.nan
                    entry[f"median_{var}"] = np.nan
                    continue
                arr = np.concatenate(chunks)
                entry[f"mean_{var}"] = float(np.mean(arr))
                entry[f"median_{var}"] = float(np.median(arr))
            summary[cat_name][lbl] = entry
    return summary


# ---------------------------------------------------------------------------
# Pass-fraction using per-entry cluster_Et (2nd pass, needed for accuracy)
# ---------------------------------------------------------------------------
def compute_pass_fraction_exact(max_events: int | None) -> dict:
    """Second pass to compute the parametric iso pass fraction using per-entry
    cluster_Et (instead of bin midpoint) for accuracy.
    """
    import awkward as ak

    branches = [
        "particle_pid",
        "particle_trkid",
        "particle_Pt",
        "particle_Eta",
        "particle_truth_iso_03",
        "particle_converted",
        "cluster_Et_CLUSTERINFO_CEMC",
        "cluster_truthtrkID_CLUSTERINFO_CEMC",
        "cluster_iso_03_CLUSTERINFO_CEMC",
    ]

    n_pt_bins = len(PT_EDGES) - 1
    pt_labels = [pt_bin_label(PT_EDGES[i], PT_EDGES[i + 1]) for i in range(n_pt_bins)]
    pt_labels.append(INCLUSIVE_LABEL)

    counts = {(c, lbl): {"n": 0, "n_pass": 0, "sum_Et": 0.0}
              for c in CAT_ORDER for lbl in pt_labels}

    f = uproot.open(FILE_PATH)
    t = f[TREE_NAME]
    n_total = t.num_entries
    entry_stop = n_total if max_events is None else min(max_events, n_total)

    for arr in t.iterate(branches, step_size=500_000, library="ak", entry_stop=entry_stop):
        particle_pid = arr["particle_pid"]
        particle_trkid = arr["particle_trkid"]
        particle_Eta = arr["particle_Eta"]
        particle_truth_iso_03 = arr["particle_truth_iso_03"]
        particle_converted = arr["particle_converted"]
        cluster_Et = arr["cluster_Et_CLUSTERINFO_CEMC"]
        cluster_truthtrkID = arr["cluster_truthtrkID_CLUSTERINFO_CEMC"]
        cluster_iso03 = arr["cluster_iso_03_CLUSTERINFO_CEMC"]

        truth_sel = (
            (particle_pid == 22)
            & (np.abs(particle_Eta) < 0.7)
            & (particle_truth_iso_03 < 4.0)
        )

        sel_trkid = particle_trkid[truth_sel]
        sel_conv = particle_converted[truth_sel]

        cross = ak.cartesian([sel_trkid, cluster_truthtrkID], axis=1, nested=True)
        match_mask = cross["0"] == cross["1"]

        # Broadcast cluster arrays against the nested (event, particle, cluster) shape
        _, cluster_Et_nested = ak.broadcast_arrays(sel_trkid[..., None], cluster_Et[:, None, :])
        _, cluster_iso_nested = ak.broadcast_arrays(sel_trkid[..., None], cluster_iso03[:, None, :])
        cluster_Et_matched = ak.where(match_mask, cluster_Et_nested, -np.inf)

        any_match = ak.any(match_mask, axis=-1)
        best_idx = ak.argmax(cluster_Et_matched, axis=-1, keepdims=False)
        best_idx_safe = ak.fill_none(best_idx, -1)

        # Flatten to per-particle entries
        nev = len(cluster_Et)
        n_part_per_event = ak.num(sel_trkid, axis=1)
        event_idx_per_part = np.repeat(np.arange(nev), np.asarray(n_part_per_event))
        flat_any = np.asarray(ak.flatten(any_match, axis=1))
        flat_best_idx = np.asarray(ak.flatten(best_idx_safe, axis=1)).astype(np.int64)
        flat_conv = np.asarray(ak.flatten(sel_conv, axis=1)).astype(np.int64)

        cluster_Et_np = ak.to_numpy(ak.flatten(cluster_Et, axis=1))
        cluster_iso_np = ak.to_numpy(ak.flatten(cluster_iso03, axis=1))
        offsets = np.asarray(ak.num(cluster_Et, axis=1))
        cluster_starts = np.concatenate([[0], np.cumsum(offsets)])[:-1]

        mask_valid = flat_any
        evt_idx_v = event_idx_per_part[mask_valid]
        loc_idx_v = flat_best_idx[mask_valid]
        flat_pos = cluster_starts[evt_idx_v] + loc_idx_v

        best_Et = cluster_Et_np[flat_pos]
        best_iso = cluster_iso_np[flat_pos]
        conv_v = flat_conv[mask_valid]

        et_cut = best_Et > 5.0
        best_Et = best_Et[et_cut]
        best_iso = best_iso[et_cut]
        conv_v = conv_v[et_cut]

        # parametric cut per entry
        iso_max_per = RECO_ISO_MAX_B + RECO_ISO_MAX_S * best_Et
        passed = best_iso < iso_max_per

        pt_bin_idx = np.digitize(best_Et, PT_EDGES) - 1
        in_range = (pt_bin_idx >= 0) & (pt_bin_idx < n_pt_bins) & (best_Et < PT_EDGES[-1])

        for cat in CAT_ORDER:
            cat_mask = conv_v == cat

            incl_mask = cat_mask & in_range
            counts[(cat, INCLUSIVE_LABEL)]["n"] += int(incl_mask.sum())
            counts[(cat, INCLUSIVE_LABEL)]["n_pass"] += int((incl_mask & passed).sum())
            counts[(cat, INCLUSIVE_LABEL)]["sum_Et"] += float(best_Et[incl_mask].sum())

            for ib in range(n_pt_bins):
                bin_mask = cat_mask & (pt_bin_idx == ib) & in_range
                if not bin_mask.any():
                    continue
                lbl = pt_bin_label(PT_EDGES[ib], PT_EDGES[ib + 1])
                counts[(cat, lbl)]["n"] += int(bin_mask.sum())
                counts[(cat, lbl)]["n_pass"] += int((bin_mask & passed).sum())
                counts[(cat, lbl)]["sum_Et"] += float(best_Et[bin_mask].sum())

    return counts


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_distributions(hists: dict) -> list[str]:
    """Make 1D overlay plots per (variable, pT bin) and cumulative plots.

    Returns list of PDF paths created.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    cat_colors = {"unconverted": "black", "converted": "red", "bad": "blue"}
    cat_markers = {"unconverted": "o", "converted": "s", "bad": "^"}
    cat_fill = {"unconverted": "full", "converted": "none", "bad": "none"}

    n_pt_bins = len(PT_EDGES) - 1
    pt_labels = [pt_bin_label(PT_EDGES[i], PT_EDGES[i + 1]) for i in range(n_pt_bins)]
    pt_labels.append(INCLUSIVE_LABEL)

    os.makedirs(FIG_DIR, exist_ok=True)
    pdf_paths = []

    # One multi-panel PDF per variable, panels = pT bins (grid) + inclusive
    for var in ISO_VARS:
        pdf_path = os.path.join(FIG_DIR, f"isolation_{var}.pdf")
        is_sub = not var.endswith("iso_03") and not var.endswith("iso_04")
        rng = HIST_RANGE_SUB if is_sub else HIST_RANGE_TOTAL
        xlabel = ISO_XLABEL[var]

        with PdfPages(pdf_path) as pdf:
            # grid: 3 x 5 = 15 panels, use 13 (12 pt bins + inclusive)
            ncols = 4
            nrows = int(np.ceil((n_pt_bins + 1) / ncols))
            fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
            axes = axes.ravel()

            for i, lbl in enumerate(pt_labels):
                ax = axes[i]
                bins = np.linspace(rng[0], rng[1], HIST_NBINS + 1)

                for cat in CAT_ORDER:
                    cat_name = CATEGORIES[cat]
                    chunks = hists[var][cat_name].get(lbl, [])
                    if not chunks:
                        continue
                    vals = np.concatenate(chunks)
                    n = vals.size
                    if n == 0:
                        continue
                    counts, edges = np.histogram(vals, bins=bins)
                    widths = np.diff(edges)
                    # normalise to area 1
                    norm = counts / (counts.sum() * widths)
                    centers = 0.5 * (edges[:-1] + edges[1:])
                    lw = 1.6
                    ls = "-" if cat_name == "unconverted" else ("--" if cat_name == "converted" else ":")
                    ax.step(
                        edges,
                        np.concatenate([[norm[0]], norm]),
                        where="pre",
                        color=cat_colors[cat_name],
                        linewidth=lw,
                        linestyle=ls,
                        label=f"{cat_name} (N={n:,})",
                    )

                ax.set_xlabel(xlabel)
                ax.set_ylabel("Normalized")
                title = f"Inclusive 8-36 GeV" if lbl == INCLUSIVE_LABEL else f"cluster E_T in {lbl.replace('to',' - ')} GeV"
                ax.set_title(title, fontsize=9)
                ax.set_xlim(rng)
                ax.set_yscale("log")
                ax.legend(fontsize=7, loc="upper right")
                ax.grid(True, alpha=0.3)

                if i == 0:
                    ax.text(
                        0.02, 0.98,
                        r"$\bf{\it{sPHENIX}}$ Simulation",
                        transform=ax.transAxes,
                        ha="left", va="top",
                        fontsize=9,
                    )
            # Hide unused axes
            for j in range(len(pt_labels), len(axes)):
                axes[j].set_visible(False)

            fig.suptitle(f"Isolation: {var}  (converted vs unconverted vs bad)")
            fig.tight_layout(rect=[0, 0, 1, 0.97])
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        pdf_paths.append(pdf_path)
        print(f"[info] wrote {pdf_path}", flush=True)

    # Cumulative distribution plot, one per iso_03 totals only (most useful)
    for var in ["cluster_iso_03", "cluster_iso_04"]:
        pdf_path = os.path.join(FIG_DIR, f"isolation_{var}_cumulative.pdf")
        rng = HIST_RANGE_TOTAL
        xlabel = ISO_XLABEL[var]

        with PdfPages(pdf_path) as pdf:
            ncols = 4
            nrows = int(np.ceil((n_pt_bins + 1) / ncols))
            fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
            axes = axes.ravel()

            for i, lbl in enumerate(pt_labels):
                ax = axes[i]
                for cat in CAT_ORDER:
                    cat_name = CATEGORIES[cat]
                    chunks = hists[var][cat_name].get(lbl, [])
                    if not chunks:
                        continue
                    vals = np.concatenate(chunks)
                    if vals.size == 0:
                        continue
                    vals_sorted = np.sort(vals)
                    cdf = np.arange(1, vals_sorted.size + 1) / vals_sorted.size
                    ls = "-" if cat_name == "unconverted" else ("--" if cat_name == "converted" else ":")
                    ax.plot(
                        vals_sorted,
                        cdf,
                        color=cat_colors[cat_name],
                        linewidth=1.5,
                        linestyle=ls,
                        label=f"{cat_name} (N={vals.size:,})",
                    )

                # Mark parametric iso cut at bin midpoint
                if lbl != INCLUSIVE_LABEL:
                    lo, hi = [float(x) for x in lbl.split("to")]
                    et_mid = 0.5 * (lo + hi)
                    iso_cut = RECO_ISO_MAX_B + RECO_ISO_MAX_S * et_mid
                    ax.axvline(iso_cut, color="grey", linestyle=":", linewidth=1.2,
                               label=f"iso cut ({iso_cut:.2f})")

                ax.set_xlabel(xlabel)
                ax.set_ylabel("Cumulative fraction")
                title = f"Inclusive 8-36 GeV" if lbl == INCLUSIVE_LABEL else f"cluster E_T in {lbl.replace('to',' - ')} GeV"
                ax.set_title(title, fontsize=9)
                ax.set_xlim(rng)
                ax.set_ylim(0, 1.05)
                ax.legend(fontsize=7, loc="lower right")
                ax.grid(True, alpha=0.3)

                if i == 0:
                    ax.text(
                        0.02, 0.98,
                        r"$\bf{\it{sPHENIX}}$ Simulation",
                        transform=ax.transAxes,
                        ha="left", va="top",
                        fontsize=9,
                    )
            for j in range(len(pt_labels), len(axes)):
                axes[j].set_visible(False)

            fig.suptitle(f"Isolation cumulative: {var}")
            fig.tight_layout(rect=[0, 0, 1, 0.97])
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        pdf_paths.append(pdf_path)
        print(f"[info] wrote {pdf_path}", flush=True)

    return pdf_paths


def plot_means_vs_pt(hists: dict) -> list[str]:
    """Plot mean iso vs pT for 4 sub-components (iso_03, iso_04)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cat_colors = {"unconverted": "black", "converted": "red", "bad": "blue"}
    cat_markers = {"unconverted": "o", "converted": "s", "bad": "^"}

    n_pt_bins = len(PT_EDGES) - 1
    pt_centres = 0.5 * (PT_EDGES[:-1] + PT_EDGES[1:])
    pt_labels = [pt_bin_label(PT_EDGES[i], PT_EDGES[i + 1]) for i in range(n_pt_bins)]

    paths = []
    for radius, vars_for_radius in [("03", ISO_VARS_R03), ("04", ISO_VARS_R04)]:
        fig, axes = plt.subplots(2, 2, figsize=(11, 8))
        axes = axes.ravel()
        labels_short = {
            f"cluster_iso_{radius}": "Total",
            f"cluster_iso_{radius}_emcal": "EMCal",
            f"cluster_iso_{radius}_hcalin": "HCalIn",
            f"cluster_iso_{radius}_hcalout": "HCalOut",
        }
        for ax, var in zip(axes, vars_for_radius):
            for cat in CAT_ORDER:
                cat_name = CATEGORIES[cat]
                if cat_name == "bad":
                    continue  # negligible population; suppress from summary plot
                means = []
                errs = []
                for lbl in pt_labels:
                    chunks = hists[var][cat_name].get(lbl, [])
                    if not chunks:
                        means.append(np.nan)
                        errs.append(np.nan)
                        continue
                    arr = np.concatenate(chunks)
                    if arr.size == 0:
                        means.append(np.nan)
                        errs.append(np.nan)
                        continue
                    means.append(float(np.mean(arr)))
                    errs.append(float(np.std(arr, ddof=1) / np.sqrt(arr.size)) if arr.size > 1 else 0.0)
                means = np.array(means)
                errs = np.array(errs)
                ax.errorbar(
                    pt_centres, means, yerr=errs,
                    color=cat_colors[cat_name],
                    marker=cat_markers[cat_name],
                    markerfacecolor=(cat_colors[cat_name] if cat_name == "unconverted" else "white"),
                    linestyle="-" if cat_name == "unconverted" else ("--" if cat_name == "converted" else ":"),
                    label=cat_name,
                    markersize=7, linewidth=1.3,
                )
            ax.set_xlabel(r"cluster $E_T$ [GeV]")
            ax.set_ylabel(rf"$\langle$iso_{radius}{' '+labels_short[var].lower() if labels_short[var]!='Total' else ''}$\rangle$ [GeV]")
            ax.set_title(f"R=0.{int(radius)} — {labels_short[var]}", fontsize=10)
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=9)
            ax.axhline(0, color="grey", linewidth=0.5)
            ax.text(
                0.02, 0.98,
                r"$\bf{\it{sPHENIX}}$ Simulation",
                transform=ax.transAxes,
                ha="left", va="top",
                fontsize=9,
            )
            ax.text(
                0.02, 0.90,
                r"photon10, $|\eta^{\gamma}|<0.7$",
                transform=ax.transAxes,
                ha="left", va="top",
                fontsize=8, color="gray",
            )
        fig.suptitle(f"Mean isolation vs cluster $E_T$ (R=0.{int(radius)})")
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        out = os.path.join(FIG_DIR, f"isolation_mean_vs_pt_R0{radius}.pdf")
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
        paths.append(out)
        print(f"[info] wrote {out}", flush=True)
    return paths


def plot_pass_fraction(counts: dict) -> str:
    """Plot parametric iso pass-fraction vs pT per category."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cat_colors = {"unconverted": "black", "converted": "red", "bad": "blue"}
    cat_markers = {"unconverted": "o", "converted": "s", "bad": "^"}

    n_pt_bins = len(PT_EDGES) - 1
    pt_centres = 0.5 * (PT_EDGES[:-1] + PT_EDGES[1:])
    pt_labels = [pt_bin_label(PT_EDGES[i], PT_EDGES[i + 1]) for i in range(n_pt_bins)]

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for cat in CAT_ORDER:
        cat_name = CATEGORIES[cat]
        if cat_name == "bad":
            continue  # negligible population; suppress from summary plot
        fracs = []
        errs = []
        for lbl in pt_labels:
            c = counts.get((cat, lbl), {"n": 0, "n_pass": 0})
            n = c["n"]
            n_pass = c["n_pass"]
            if n == 0:
                fracs.append(np.nan)
                errs.append(0)
                continue
            p = n_pass / n
            fracs.append(p)
            errs.append(np.sqrt(max(p * (1 - p), 0) / n))
        fracs = np.array(fracs)
        errs = np.array(errs)
        ax.errorbar(
            pt_centres, fracs, yerr=errs,
            color=cat_colors[cat_name],
            marker=cat_markers[cat_name],
            markerfacecolor=(cat_colors[cat_name] if cat_name == "unconverted" else "white"),
            linestyle="-" if cat_name == "unconverted" else ("--" if cat_name == "converted" else ":"),
            label=cat_name,
            markersize=8, linewidth=1.5,
        )
    ax.set_xlabel(r"cluster $E_T$ [GeV]")
    ax.set_ylabel("Fraction passing parametric iso cut (iso_03 < b + s·$E_T$)")
    ax.set_ylim(0.0, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11, loc="lower right")
    ax.text(
        0.02, 0.98,
        r"$\bf{\it{sPHENIX}}$ Simulation",
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=11,
    )
    ax.text(
        0.02, 0.92,
        r"photon10, $|\eta^{\gamma}|<0.7$",
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=9, color="gray",
    )
    ax.text(
        0.02, 0.86,
        rf"iso cut: $\max = {RECO_ISO_MAX_B:.4f} + {RECO_ISO_MAX_S:.5f} \cdot E_T$ [GeV]",
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=9,
    )
    out = os.path.join(FIG_DIR, "isolation_pass_fraction_vs_pt.pdf")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"[info] wrote {out}", flush=True)
    return out


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--max-events", type=int, default=None,
                   help="Process at most N events (default: all)")
    p.add_argument("--json-out", default=os.path.join(OUT_DIR, "summary.json"))
    args = p.parse_args()

    os.makedirs(os.path.dirname(ROOT_OUT), exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)

    print("=" * 70)
    print("[Pass 1] Fill per-cat/pT iso histograms")
    print("=" * 70)
    hists, stats = run_analysis(args.max_events)

    print("[info] writing ROOT output")
    write_root_file(hists)

    print("=" * 70)
    print("[Pass 2] Exact parametric iso pass-fraction using per-entry cluster_Et")
    print("=" * 70)
    counts = compute_pass_fraction_exact(args.max_events)

    # Combine into summary
    summary = compute_summary(hists, stats)
    # Replace pass-fraction fields with the exact ones
    for (cat, lbl), c in counts.items():
        cat_name = CATEGORIES[cat]
        entry = summary.get(cat_name, {}).get(lbl, {})
        entry["n_exact"] = c["n"]
        entry["n_pass_exact"] = c["n_pass"]
        entry["pass_frac_exact"] = (c["n_pass"] / c["n"]) if c["n"] > 0 else np.nan
        entry["mean_cluster_Et"] = (c["sum_Et"] / c["n"]) if c["n"] > 0 else np.nan
        summary[cat_name][lbl] = entry

    with open(args.json_out, "w") as fp:
        json.dump(summary, fp, indent=2, default=float)
    print(f"[info] wrote summary to {args.json_out}", flush=True)

    print("=" * 70)
    print("[Plot] distributions / cumulative")
    print("=" * 70)
    plot_distributions(hists)

    print("=" * 70)
    print("[Plot] mean iso vs pT")
    print("=" * 70)
    plot_means_vs_pt(hists)

    print("=" * 70)
    print("[Plot] pass fraction vs pT")
    print("=" * 70)
    plot_pass_fraction(counts)

    print("Done.")


if __name__ == "__main__":
    main()
