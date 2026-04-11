#!/usr/bin/env python3
"""
extract_data_vertex.py — Step B2 of the Section 3 luminosity correction.

Extracts the z-vertex distribution from data per run. This is the beam vertex
profile that will be convolved with the MC-derived MBD response function.

The vertex distribution is measured from ALL events (no trigger or vertex cut),
using the MBD-reconstructed vertex. This matches the denominator of the
no-vertex-cut luminosity L(noVtx).

Output: results/data_vertex_distributions.root
  - h_vtxz_run_{runnumber}: per-run vertex distribution
  - h_vtxz_period_0mrad:    combined 0 mrad vertex distribution
  - h_vtxz_period_1p5mrad:  combined 1.5 mrad vertex distribution
  - h_vtxz_all:             all runs combined

Usage:
    python3 extract_data_vertex.py --config config_bdt_nom.yaml
"""

import argparse
import glob
import os
from collections import defaultdict

import numpy as np
import uproot
import awkward as ak
import yaml

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
Z_EDGES = np.linspace(-200, 200, 201)  # must match extract_mbd_response.py
N_Z_BINS = len(Z_EDGES) - 1
XING_RUN = 51274  # 0 mrad / 1.5 mrad boundary


def main():
    parser = argparse.ArgumentParser(description="Extract data vertex distributions per run")
    parser.add_argument("--config", default="config_bdt_nom.yaml")
    parser.add_argument("--results", default="results")
    parser.add_argument("--max-events-per-run", type=int, default=-1,
                        help="Limit events per run for testing")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    data_glob = config["input"]["data_file"]
    treename = config["input"]["tree"]

    files = sorted(glob.glob(data_glob))
    print(f"Processing {len(files)} data files for vertex extraction...")

    # Per-run vertex histograms
    run_vtxz = defaultdict(lambda: np.zeros(N_Z_BINS))
    run_nevents = defaultdict(int)

    for fi, fname in enumerate(files):
        try:
            f = uproot.open(fname)
            tree = f[treename]
        except Exception as e:
            print(f"  SKIP {fname}: {e}")
            continue

        for batch in tree.iterate(["runnumber", "vertexz"], step_size=500000, library="np"):
            rns = batch["runnumber"]
            vtxz = batch["vertexz"]

            for rn in np.unique(rns):
                mask = rns == rn
                z_vals = vtxz[mask]

                # Histogram the vertex distribution
                hist, _ = np.histogram(z_vals, bins=Z_EDGES)
                run_vtxz[rn] += hist
                run_nevents[rn] += int(mask.sum())

        if (fi + 1) % 10 == 0 or fi == len(files) - 1:
            print(f"  [{fi+1}/{len(files)}] {len(run_nevents)} runs extracted")

    # --- Combine into period distributions ---
    vtxz_0mrad = np.zeros(N_Z_BINS)
    vtxz_1p5mrad = np.zeros(N_Z_BINS)
    vtxz_all = np.zeros(N_Z_BINS)
    n_runs_0 = 0
    n_runs_1p5 = 0

    for rn in sorted(run_nevents.keys()):
        vtxz_all += run_vtxz[rn]
        if rn < XING_RUN:
            vtxz_0mrad += run_vtxz[rn]
            n_runs_0 += 1
        else:
            vtxz_1p5mrad += run_vtxz[rn]
            n_runs_1p5 += 1

    # --- Summary ---
    print(f"\n{'='*60}")
    print("DATA VERTEX DISTRIBUTION SUMMARY")
    print(f"{'='*60}")
    print(f"Total runs: {len(run_nevents)} ({n_runs_0} at 0 mrad, {n_runs_1p5} at 1.5 mrad)")
    print(f"Total events: {sum(run_nevents.values()):,}")

    for label, hist in [("0 mrad", vtxz_0mrad), ("1.5 mrad", vtxz_1p5mrad), ("All", vtxz_all)]:
        total = hist.sum()
        if total == 0:
            continue
        centers = 0.5 * (Z_EDGES[:-1] + Z_EDGES[1:])
        mean = np.average(centers, weights=hist)
        rms = np.sqrt(np.average((centers - mean)**2, weights=hist))
        within_60 = hist[(Z_EDGES[:-1] >= -60) & (Z_EDGES[:-1] < 60)].sum()
        print(f"\n  {label}: {total:.0f} events")
        print(f"    Mean:         {mean:.1f} cm")
        print(f"    RMS:          {rms:.1f} cm")
        print(f"    Within ±60:   {within_60:.0f} ({100*within_60/total:.1f}%)")

    # --- Save ---
    outpath = os.path.join(args.results, "data_vertex_distributions.root")
    with uproot.recreate(outpath) as fout:
        fout["h_vtxz_period_0mrad"] = (vtxz_0mrad, Z_EDGES)
        fout["h_vtxz_period_1p5mrad"] = (vtxz_1p5mrad, Z_EDGES)
        fout["h_vtxz_all"] = (vtxz_all, Z_EDGES)

        # Save per-run (only runs with >1000 events to avoid clutter)
        n_saved = 0
        for rn in sorted(run_nevents.keys()):
            if run_nevents[rn] > 1000:
                fout[f"h_vtxz_run_{rn}"] = (run_vtxz[rn], Z_EDGES)
                n_saved += 1

        # z_edges are embedded in the histogram axis — no need to save separately

    print(f"\nSaved: {outpath} ({n_saved} per-run histograms)")


if __name__ == "__main__":
    main()
