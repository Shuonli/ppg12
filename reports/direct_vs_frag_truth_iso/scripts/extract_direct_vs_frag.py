#!/usr/bin/env python3
"""
PPG12 - Direct vs Fragmentation photon fraction extraction from Pythia truth.

Uses the CaloAna24 classification stored in the slimtree branch
``particle_photonclass``:

   1 = direct (2->2 hard scatter with two colored incomers/outgoers, e.g.
       qg -> gamma q, q qbar -> gamma g)
   2 = fragmentation (1->2 parton splitting q -> q gamma inside the shower)
   3 = decay photon
   0 = unclassified

See ``anatreemaker/source/CaloAna24.cc::photon_type`` (lines 2034-2168).

Weighting mirrors ``efficiencytool/RecoEffCalculator.C`` and
``efficiencytool/CrossSectionWeights.h``:

   weight = {photon5,photon10,photon20}cross / photon20cross
   photon5cross  = 146359.3
   photon10cross = 6944.675
   photon20cross = 130.4461

Overlap handling: for each event use max truth photon pT (over ALL photons
with pid==22) and keep only events with that max in the sample's truth pT
window (photon5: [0,14], photon10: [14,30], photon20: [30,200]).  Mirrors
RecoEffCalculator.C lines 1289-1294.  Prevents double counting.

Implementation: vectorized with ``awkward`` arrays; no per-event Python loop.

Outputs
-------
  reports/direct_vs_frag_truth_iso/direct_vs_frag_truth.root
  reports/direct_vs_frag_truth_iso/fraction_scan.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass

import awkward as ak
import numpy as np
import uproot

# ---------------------------------------------------------------------------
# Cross-section weights (in sync with efficiencytool/CrossSectionWeights.h)
# ---------------------------------------------------------------------------
PHOTON5_CROSS = 146359.3
PHOTON10_CROSS = 6944.675
PHOTON20_CROSS = 130.4461


@dataclass
class Sample:
    name: str
    path: str
    pt_lower: float
    pt_upper: float
    weight: float  # cross / photon20cross  (relative)


SAMPLES = [
    Sample("photon5",
           "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon5/bdt_split.root",
           0.0, 14.0, PHOTON5_CROSS / PHOTON20_CROSS),
    Sample("photon10",
           "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root",
           14.0, 30.0, PHOTON10_CROSS / PHOTON20_CROSS),
    Sample("photon20",
           "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon20/bdt_split.root",
           30.0, 200.0, 1.0),
]

# ---------------------------------------------------------------------------
# Analysis binning
# ---------------------------------------------------------------------------
TRUTH_ETA_MAX = 0.7
TRUTH_PT_MIN = 8.0
TRUTH_PT_MAX = 36.0

# Coarse truth pT bins (GeV)
PT_BINS = [(8.0, 14.0), (14.0, 22.0), (22.0, 36.0)]

# Iso cut scan (GeV).  math.inf -> no iso cut.
ISO_CUTS = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, math.inf]

# Fine iso binning: 0.25 GeV edges from 0 to 15 GeV (60 bins)
ISO_EDGES = np.arange(0.0, 15.0 + 1e-9, 0.25)  # [0, 0.25, ..., 15.0]
N_ISO_BINS = len(ISO_EDGES) - 1  # 60

BRANCHES = [
    "particle_pid",
    "particle_Pt",
    "particle_Eta",
    "particle_photonclass",
    "particle_truth_iso_03",
]


def process_sample(sample: Sample, step_size="1 GB", verbose: bool = True):
    """Stream over the file in chunks; use fully vectorized selections."""
    result = {
        "hist": {
            ("direct", k): np.zeros(N_ISO_BINS, dtype=np.float64)
            for k in range(len(PT_BINS))
        },
        "counts": defaultdict(float),
        "n_events_total": 0,
        "n_events_kept": 0,
        "n_photons_fiducial": 0,
        "n_decay_photons_fiducial": 0,
        "n_unclassified_photons_fiducial": 0,
    }
    for k in range(len(PT_BINS)):
        result["hist"][("frag", k)] = np.zeros(N_ISO_BINS, dtype=np.float64)

    if not os.path.exists(sample.path):
        raise FileNotFoundError(f"Missing sample file: {sample.path}")

    if verbose:
        print(f"[{sample.name}] opening {sample.path}", flush=True)
        print(f"[{sample.name}] truth pT window [{sample.pt_lower}, {sample.pt_upper}] GeV, "
              f"weight={sample.weight:.6g}", flush=True)

    f = uproot.open(sample.path)
    tree = f["slimtree"]
    n_total = tree.num_entries
    if verbose:
        print(f"[{sample.name}] total events: {n_total:,}", flush=True)
    result["n_events_total"] = n_total

    w_sample = sample.weight

    # Stream over the file in chunks to keep memory under control.
    chunk_index = 0
    for chunk in tree.iterate(BRANCHES, step_size=step_size, library="ak"):
        chunk_index += 1
        n_ev = len(chunk)
        if verbose:
            print(f"[{sample.name}]   chunk {chunk_index}: {n_ev:,} events", flush=True)

        pids = chunk["particle_pid"]
        pts = chunk["particle_Pt"]
        etas = chunk["particle_Eta"]
        clses = chunk["particle_photonclass"]
        isos = chunk["particle_truth_iso_03"]

        # Per-particle photon mask (inside each event)
        is_photon = (pids == 22)

        # Max photon pT per event (for overlap-window cut)
        # ak.fill_none in case an event has no photons (shouldn't happen in
        # signal MC but be safe).
        ph_pts = ak.mask(pts, is_photon)  # particles_per_event, masked
        max_ph_pt_per_ev = ak.fill_none(ak.max(ph_pts, axis=1), -1.0)

        # Event-level window cut
        ev_in_window = (max_ph_pt_per_ev >= sample.pt_lower) & \
                       (max_ph_pt_per_ev < sample.pt_upper)
        # note: use half-open interval; for photon20 pt_upper=200 covers tail
        result["n_events_kept"] += int(ak.sum(ev_in_window))

        # Broadcast event mask down to particle level
        ev_mask_broad = ev_in_window[:, None]  # shape (n_ev, 1) but awkward
        # We'll use it via "ak.broadcast_arrays" implicitly via &

        # Build per-photon selection (flattened later):
        fid_mask = is_photon & \
                   (abs(etas) < TRUTH_ETA_MAX) & \
                   (pts >= TRUTH_PT_MIN) & (pts < TRUTH_PT_MAX)

        # Combine event-level cut
        sel_mask = fid_mask & ev_mask_broad

        # Flatten arrays, keep only selected photons
        flat_pt = ak.to_numpy(ak.flatten(pts[sel_mask]))
        flat_cls = ak.to_numpy(ak.flatten(clses[sel_mask]))
        flat_iso = ak.to_numpy(ak.flatten(isos[sel_mask]))

        n_fid = len(flat_pt)
        result["n_photons_fiducial"] += n_fid
        if n_fid == 0:
            continue

        # Tally decay/unclassified
        is_decay = (flat_cls == 3)
        is_unkn = (flat_cls == 0) | (flat_cls < 0) | \
                  (~np.isin(flat_cls, [1, 2, 3]))
        result["n_decay_photons_fiducial"] += int(is_decay.sum())
        result["n_unclassified_photons_fiducial"] += int(is_unkn.sum())

        # Keep only direct (1) / frag (2)
        keep = (flat_cls == 1) | (flat_cls == 2)
        if not keep.any():
            continue
        flat_pt = flat_pt[keep]
        flat_cls = flat_cls[keep]
        flat_iso = flat_iso[keep]

        # Truth pT bin index (vectorized)
        pt_bin_idx = np.full(len(flat_pt), -1, dtype=np.int8)
        for k, (lo, hi) in enumerate(PT_BINS):
            pt_bin_idx[(flat_pt >= lo) & (flat_pt < hi)] = k
        in_ptbin = pt_bin_idx >= 0

        flat_pt = flat_pt[in_ptbin]
        flat_cls = flat_cls[in_ptbin]
        flat_iso = flat_iso[in_ptbin]
        pt_bin_idx = pt_bin_idx[in_ptbin]

        # Clip iso for the fine histogram (overflow to last bin)
        iso_for_hist = np.clip(flat_iso, 0.0, ISO_EDGES[-1] - 1e-9)

        # Loop over the 6 (class, pt_bin) cells (very cheap)
        for cls_val, cls_name in [(1, "direct"), (2, "frag")]:
            for k in range(len(PT_BINS)):
                cell_mask = (flat_cls == cls_val) & (pt_bin_idx == k)
                if not cell_mask.any():
                    continue
                cell_iso = iso_for_hist[cell_mask]
                # Fine histogram (weighted, constant sample weight)
                hvals, _ = np.histogram(cell_iso, bins=ISO_EDGES)
                result["hist"][(cls_name, k)] += hvals * w_sample

                raw_iso = flat_iso[cell_mask]
                # Iso cut scan (vectorized)
                for iso_cut in ISO_CUTS:
                    if math.isinf(iso_cut):
                        result["counts"][(cls_name, k, iso_cut)] += \
                            len(raw_iso) * w_sample
                    else:
                        result["counts"][(cls_name, k, iso_cut)] += \
                            int((raw_iso < iso_cut).sum()) * w_sample

    if verbose:
        print(f"[{sample.name}] events kept: {result['n_events_kept']:,} / {n_total:,}", flush=True)
        print(f"[{sample.name}] fiducial photons: {result['n_photons_fiducial']:,}  "
              f"(decay: {result['n_decay_photons_fiducial']:,}, "
              f"unclassified: {result['n_unclassified_photons_fiducial']:,})", flush=True)
    return result


def merge_results(per_sample):
    combined = {
        "hist": {
            (cls, k): np.zeros(N_ISO_BINS, dtype=np.float64)
            for cls in ("direct", "frag")
            for k in range(len(PT_BINS))
        },
        "counts": defaultdict(float),
    }
    for res in per_sample.values():
        for key, hist in res["hist"].items():
            combined["hist"][key] += hist
        for key, val in res["counts"].items():
            combined["counts"][key] += val
    return combined


def write_outputs(combined, per_sample, out_root, out_csv):
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "pt_bin_low", "pt_bin_high", "iso_cut",
            "N_direct", "N_frag", "N_total",
            "f_direct", "f_frag",
        ])
        for k, (lo, hi) in enumerate(PT_BINS):
            for iso_cut in ISO_CUTS:
                nd = combined["counts"].get(("direct", k, iso_cut), 0.0)
                nf = combined["counts"].get(("frag", k, iso_cut), 0.0)
                tot = nd + nf
                fd = nd / tot if tot > 0 else 0.0
                ff = nf / tot if tot > 0 else 0.0
                iso_str = "inf" if math.isinf(iso_cut) else f"{iso_cut}"
                w.writerow([lo, hi, iso_str,
                            f"{nd:.6g}", f"{nf:.6g}", f"{tot:.6g}",
                            f"{fd:.6f}", f"{ff:.6f}"])
    print(f"[write] wrote {out_csv}", flush=True)

    with uproot.recreate(out_root) as fout:
        for (cls_name, k), hist in combined["hist"].items():
            lo, hi = PT_BINS[k]
            hname = f"h_iso_{cls_name}_pt{int(lo)}_{int(hi)}"
            fout[hname] = (hist, ISO_EDGES)
        for sname, res in per_sample.items():
            for (cls_name, k), hist in res["hist"].items():
                lo, hi = PT_BINS[k]
                hname = f"h_iso_{cls_name}_pt{int(lo)}_{int(hi)}_{sname}"
                fout[hname] = (hist, ISO_EDGES)
        ptlo_arr, pthi_arr, iso_arr = [], [], []
        nd_arr, nf_arr, tot_arr, fd_arr, ff_arr = [], [], [], [], []
        for k, (lo, hi) in enumerate(PT_BINS):
            for iso_cut in ISO_CUTS:
                nd = combined["counts"].get(("direct", k, iso_cut), 0.0)
                nf = combined["counts"].get(("frag", k, iso_cut), 0.0)
                tot = nd + nf
                fd = nd / tot if tot > 0 else 0.0
                ff = nf / tot if tot > 0 else 0.0
                iso_val = 1e9 if math.isinf(iso_cut) else float(iso_cut)
                ptlo_arr.append(lo); pthi_arr.append(hi); iso_arr.append(iso_val)
                nd_arr.append(nd); nf_arr.append(nf); tot_arr.append(tot)
                fd_arr.append(fd); ff_arr.append(ff)
        fout["scan_tree"] = {
            "pt_bin_low":  np.asarray(ptlo_arr, dtype=np.float32),
            "pt_bin_high": np.asarray(pthi_arr, dtype=np.float32),
            "iso_cut":     np.asarray(iso_arr, dtype=np.float32),
            "N_direct":    np.asarray(nd_arr, dtype=np.float64),
            "N_frag":      np.asarray(nf_arr, dtype=np.float64),
            "N_total":     np.asarray(tot_arr, dtype=np.float64),
            "f_direct":    np.asarray(fd_arr, dtype=np.float64),
            "f_frag":      np.asarray(ff_arr, dtype=np.float64),
        }
    print(f"[write] wrote {out_root}", flush=True)


def sanity_print(combined, per_sample):
    print("\n================================================================", flush=True)
    print("Sanity: f_direct at iso<4 GeV (PPG12 nominal fiducial)", flush=True)
    print("================================================================", flush=True)
    print(f"  {'pT bin [GeV]':<16}{'N_direct':>14}{'N_frag':>14}"
          f"{'f_direct':>12}{'f_frag':>12}", flush=True)
    for k, (lo, hi) in enumerate(PT_BINS):
        nd = combined["counts"].get(("direct", k, 4.0), 0.0)
        nf = combined["counts"].get(("frag", k, 4.0), 0.0)
        tot = nd + nf
        fd = nd / tot if tot > 0 else 0.0
        print(f"  {f'[{lo:.0f},{hi:.0f})':<16}{nd:>14.4g}{nf:>14.4g}"
              f"{fd:>12.4f}{1-fd:>12.4f}", flush=True)
    print("\n  Per-sample contribution at iso<4 GeV (sum over all pT bins):", flush=True)
    print(f"  {'sample':<12}{'N_direct':>14}{'N_frag':>14}"
          f"{'direct/(d+f)':>16}", flush=True)
    for sname, res in per_sample.items():
        nd = sum(res["counts"].get(("direct", k, 4.0), 0.0)
                 for k in range(len(PT_BINS)))
        nf = sum(res["counts"].get(("frag", k, 4.0), 0.0)
                 for k in range(len(PT_BINS)))
        tot = nd + nf
        fd = nd / tot if tot > 0 else 0.0
        print(f"  {sname:<12}{nd:>14.4g}{nf:>14.4g}{fd:>16.4f}", flush=True)
    print("================================================================", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir",
                   default=os.path.join(os.path.dirname(__file__), ".."))
    p.add_argument("--samples", nargs="+", default=None)
    p.add_argument("--step-size", default="1 GB",
                   help="uproot iterate step_size (e.g. '500 MB', '1 GB')")
    args = p.parse_args()

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    out_root = os.path.join(outdir, "direct_vs_frag_truth.root")
    out_csv = os.path.join(outdir, "fraction_scan.csv")

    samples = SAMPLES
    if args.samples:
        samples = [s for s in SAMPLES if s.name in args.samples]
        if not samples:
            print("No matching samples; aborting.", file=sys.stderr)
            sys.exit(1)

    per_sample = {}
    for s in samples:
        per_sample[s.name] = process_sample(s, step_size=args.step_size)

    combined = merge_results(per_sample)
    write_outputs(combined, per_sample, out_root, out_csv)
    sanity_print(combined, per_sample)


if __name__ == "__main__":
    main()
