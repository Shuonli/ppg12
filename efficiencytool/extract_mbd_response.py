#!/usr/bin/env python3
"""
extract_mbd_response.py — Step B1 of the Section 3 luminosity correction.

Extracts the MBD response function P(MBD fires AND |z_reco|<60 | z_truth, pT)
from single-interaction and double-interaction MC.

This is a DETECTOR property, independent of the beam vertex profile.
It will be convolved with the data vertex distribution to compute f_corr.

Also extracts Pass₂ from double-interaction MC for the pileup correction.

Output: results/mbd_response_function.root
  - h_response_single[pT_bin]:  P(MBD+vtx | z_truth) per pT bin, single interaction
  - h_response_double[pT_bin]:  P(MBD+vtx | z_truth) per pT bin, double interaction
  - h_denom_single[pT_bin]:     denominator (truth photons per z_truth bin)
  - h_numer_single[pT_bin]:     numerator (passing MBD + vtx)
  - h_response_single_2d:       2D (z_truth, pT) response

Usage:
    python3 extract_mbd_response.py
    python3 extract_mbd_response.py --max-events 2000000
"""

import argparse
import os
import numpy as np
import uproot
import awkward as ak

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
ETA_MAX = 0.7
PT_MIN = 8.0
PT_MAX = 36.0
ISO_MAX = 4.0
VERTEX_CUT_RECO = 60.0  # reco vertex cut applied in data

PT_EDGES = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
N_PT_BINS = len(PT_EDGES) - 1

# z_truth binning — wide enough to cover the full MC truth vertex range
Z_EDGES = np.linspace(-200, 200, 201)  # 2 cm bins from -200 to +200 cm
N_Z_BINS = len(Z_EDGES) - 1

XSEC_WEIGHTS = {
    "photon5":  146359.3 / 130.4461,
    "photon10": 6944.675 / 130.4461,
    "photon20": 1.0,
    "photon10_double": 6944.675 / 130.4461,
}

SAMPLES_SINGLE = {
    "photon5":  "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon5/bdt_split.root",
    "photon10": "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root",
    "photon20": "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon20/bdt_split.root",
}

SAMPLES_DOUBLE = {
    "photon10_double": "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_double/bdt_split.root",
}

BRANCHES = [
    "vertexz", "vertexz_truth", "mbdnorthhit", "mbdsouthhit",
    "particle_Pt", "particle_Eta", "particle_pid",
    "particle_photonclass", "particle_truth_iso_03",
]


def process_sample(path, weight, max_events, label):
    """Extract per-(z_truth, pT) acceptance from one sample."""
    f = uproot.open(path)
    tree = f["slimtree"]

    # 2D accumulators: [z_bin, pT_bin]
    denom = np.zeros((N_Z_BINS, N_PT_BINS))
    numer = np.zeros((N_Z_BINS, N_PT_BINS))

    n_events = 0
    n_truth = 0

    for batch in tree.iterate(BRANCHES, step_size=200000):
        vtxz = ak.to_numpy(batch["vertexz"])
        vtxz_truth = ak.to_numpy(batch["vertexz_truth"])
        mbd_n = ak.to_numpy(batch["mbdnorthhit"])
        mbd_s = ak.to_numpy(batch["mbdsouthhit"])

        p_pt = batch["particle_Pt"]
        p_eta = batch["particle_Eta"]
        p_pid = batch["particle_pid"]
        p_class = batch["particle_photonclass"]
        p_iso = batch["particle_truth_iso_03"]

        nevt = len(vtxz)
        n_events += nevt
        if max_events > 0 and n_events > max_events:
            break

        # Truth photon selection
        pmask = ((p_pid == 22) & (p_class < 3) & (p_iso < ISO_MAX)
                 & (abs(p_eta) < ETA_MAX) & (p_pt >= PT_MIN) & (p_pt < PT_MAX))

        selected_pt = p_pt[pmask]
        mbd_fires = (mbd_n >= 1) & (mbd_s >= 1)
        reco_vtx_ok = np.abs(vtxz) < VERTEX_CUT_RECO

        for ievt in range(nevt):
            pts = ak.to_numpy(selected_pt[ievt])
            if len(pts) == 0:
                continue

            z_t = vtxz_truth[ievt]
            iz = np.searchsorted(Z_EDGES, z_t, side="right") - 1
            if iz < 0 or iz >= N_Z_BINS:
                continue

            for pt in pts:
                ipt = np.searchsorted(PT_EDGES, pt, side="right") - 1
                if ipt < 0 or ipt >= N_PT_BINS:
                    continue

                n_truth += 1
                denom[iz, ipt] += weight

                if mbd_fires[ievt] and reco_vtx_ok[ievt]:
                    numer[iz, ipt] += weight

        if max_events > 0 and n_events >= max_events:
            break

    print(f"  [{label}] {n_events} events, {n_truth} truth photons")
    return denom, numer


def main():
    parser = argparse.ArgumentParser(description="Extract MBD response function")
    parser.add_argument("--max-events", type=int, default=2_000_000)
    parser.add_argument("--results", default="results")
    args = parser.parse_args()

    os.makedirs(args.results, exist_ok=True)

    # --- Single interaction ---
    print("\n=== Single-interaction MC ===")
    denom_single = np.zeros((N_Z_BINS, N_PT_BINS))
    numer_single = np.zeros((N_Z_BINS, N_PT_BINS))

    for name, path in SAMPLES_SINGLE.items():
        w = XSEC_WEIGHTS[name]
        d, n = process_sample(path, w, args.max_events, name)
        denom_single += d
        numer_single += n

    # --- Double interaction ---
    print("\n=== Double-interaction MC ===")
    denom_double = np.zeros((N_Z_BINS, N_PT_BINS))
    numer_double = np.zeros((N_Z_BINS, N_PT_BINS))

    for name, path in SAMPLES_DOUBLE.items():
        w = XSEC_WEIGHTS[name]
        d, n = process_sample(path, w, args.max_events, name)
        denom_double += d
        numer_double += n

    # --- Compute response functions ---
    with np.errstate(divide="ignore", invalid="ignore"):
        response_single = np.where(denom_single > 0, numer_single / denom_single, 0)
        response_double = np.where(denom_double > 0, numer_double / denom_double, 0)

    # --- Validation: integrate response over flat z to get pT-dependent vertex_eff ---
    print(f"\n{'='*70}")
    print("VALIDATION: Integrate response over z_truth (flat) to recover vertex_eff")
    print(f"{'='*70}")
    print(f"\n{'pT bin':<14} {'vertex_eff':>12} {'(from resp)':>12} {'Pass2':>10}")
    print("-" * 50)

    for ipt in range(N_PT_BINS):
        # weighted average of response over z (weighted by denominator = truth photon density)
        d_z = denom_single[:, ipt]
        n_z = numer_single[:, ipt]
        total_d = d_z.sum()
        total_n = n_z.sum()
        vtx_eff = total_n / total_d if total_d > 0 else 0

        # Same for double
        d_z2 = denom_double[:, ipt]
        n_z2 = numer_double[:, ipt]
        total_d2 = d_z2.sum()
        total_n2 = n_z2.sum()
        pass2 = total_n2 / total_d2 if total_d2 > 0 else 0

        print(f"{PT_EDGES[ipt]:4.0f}-{PT_EDGES[ipt+1]:4.0f} GeV  {vtx_eff:12.4f} {'(MC avg)':>12} {pass2:10.4f}")

    # --- Save to ROOT file ---
    outpath = os.path.join(args.results, "mbd_response_function.root")
    with uproot.recreate(outpath) as fout:
        # Save 2D response as histograms
        fout["h_response_single_2d"] = (response_single, Z_EDGES, PT_EDGES)
        fout["h_response_double_2d"] = (response_double, Z_EDGES, PT_EDGES)
        fout["h_denom_single_2d"] = (denom_single, Z_EDGES, PT_EDGES)
        fout["h_numer_single_2d"] = (numer_single, Z_EDGES, PT_EDGES)
        fout["h_denom_double_2d"] = (denom_double, Z_EDGES, PT_EDGES)
        fout["h_numer_double_2d"] = (numer_double, Z_EDGES, PT_EDGES)

        # Save per-pT 1D response vs z_truth
        for ipt in range(N_PT_BINS):
            tag = f"{PT_EDGES[ipt]:.0f}_{PT_EDGES[ipt+1]:.0f}"
            fout[f"h_response_single_pt{tag}"] = (response_single[:, ipt], Z_EDGES)
            fout[f"h_response_double_pt{tag}"] = (response_double[:, ipt], Z_EDGES)
            fout[f"h_denom_single_pt{tag}"] = (denom_single[:, ipt], Z_EDGES)
            fout[f"h_denom_double_pt{tag}"] = (denom_double[:, ipt], Z_EDGES)

        # pT and z edges are embedded in the histogram axes — no separate save needed

    print(f"\nSaved: {outpath}")
    print(f"  {N_Z_BINS} z bins × {N_PT_BINS} pT bins")
    print(f"  z range: [{Z_EDGES[0]:.0f}, {Z_EDGES[-1]:.0f}] cm, bin width {Z_EDGES[1]-Z_EDGES[0]:.0f} cm")


if __name__ == "__main__":
    main()
