#!/usr/bin/env python3
"""
study_pileup_vertex_acceptance.py — Section 3 pileup vertex+MBD acceptance study.

Measures the MBD coincidence + z-vertex acceptance (vertex_eff) separately for
single-interaction and double-interaction MC.  Computes per-pT blended correction
factors for both crossing-angle periods (0 mrad, 1.5 mrad).

This quantifies the "Section 3 correction" from the sPHENIX luminosity note:
pileup overlay increases the effective vertex+MBD acceptance because the overlay
collision provides MBD hits and shifts the averaged vertex toward the fiducial cut.

Usage:
    python3 study_pileup_vertex_acceptance.py
    python3 study_pileup_vertex_acceptance.py --max-events 500000
"""

import argparse
import os
import sys

import numpy as np
import uproot
import awkward as ak

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
VERTEX_CUT = 60.0
ETA_MAX = 0.7
PT_MIN = 8.0
PT_MAX = 36.0
ISO_MAX = 4.0

PT_EDGES = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
N_PT_BINS = len(PT_EDGES) - 1

# Event-level double-interaction fractions (from calc_pileup_range.C)
# These are P(n=2|n>=1) + P(n>=3|n>=1) folded together
PILEUP_FRACTIONS = {
    "0mrad":   {"f_double": 0.111, "label": "0 mrad"},
    "1p5mrad": {"f_double": 0.039, "label": "1.5 mrad"},
}

# Cross-section weights (from CrossSectionWeights.h, normalized to photon20)
XSEC_WEIGHTS = {
    "photon5":  146359.3 / 130.4461,
    "photon10": 6944.675 / 130.4461,
    "photon20": 1.0,
}

# Sample file paths and pT truth windows
SAMPLES_SINGLE = {
    "photon5":  {
        "path": "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon5/bdt_split.root",
        "pt_lo": 0, "pt_hi": 14, "weight": XSEC_WEIGHTS["photon5"],
    },
    "photon10": {
        "path": "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root",
        "pt_lo": 14, "pt_hi": 30, "weight": XSEC_WEIGHTS["photon10"],
    },
    "photon20": {
        "path": "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon20/bdt_split.root",
        "pt_lo": 30, "pt_hi": 200, "weight": XSEC_WEIGHTS["photon20"],
    },
}

SAMPLES_DOUBLE = {
    "photon10_double": {
        "path": "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_double/bdt_split.root",
        "pt_lo": 10, "pt_hi": 100, "weight": XSEC_WEIGHTS["photon10"],
    },
}

BRANCHES = [
    "vertexz", "vertexz_truth", "mbdnorthhit", "mbdsouthhit",
    "particle_Pt", "particle_Eta", "particle_pid",
    "particle_photonclass", "particle_truth_iso_03",
]


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
def process_sample(path, weight, max_events, label):
    """Extract vertex+MBD acceptance from one sample."""
    f = uproot.open(path)
    tree = f["slimtree"]

    # Per-pT accumulators (weighted)
    denom = np.zeros(N_PT_BINS)
    numer_mbd_vtx = np.zeros(N_PT_BINS)
    numer_mbd_only = np.zeros(N_PT_BINS)
    numer_vtx_only = np.zeros(N_PT_BINS)
    numer_neither = np.zeros(N_PT_BINS)

    # MBD hit decomposition (weighted)
    mbd_both = np.zeros(N_PT_BINS)
    mbd_north_only = np.zeros(N_PT_BINS)
    mbd_south_only = np.zeros(N_PT_BINS)
    mbd_neither = np.zeros(N_PT_BINS)

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
        truth_vtx_ok = np.abs(vtxz_truth) < VERTEX_CUT
        reco_vtx_ok = np.abs(vtxz) < VERTEX_CUT
        mbd_fires = (mbd_n >= 1) & (mbd_s >= 1)
        has_north = mbd_n >= 1
        has_south = mbd_s >= 1

        for ievt in range(nevt):
            pts = ak.to_numpy(selected_pt[ievt])
            if len(pts) == 0 or not truth_vtx_ok[ievt]:
                continue

            for pt in pts:
                ibin = np.searchsorted(PT_EDGES, pt, side="right") - 1
                if ibin < 0 or ibin >= N_PT_BINS:
                    continue

                n_truth += 1
                denom[ibin] += weight

                # MBD hit decomposition
                if has_north[ievt] and has_south[ievt]:
                    mbd_both[ibin] += weight
                elif has_north[ievt]:
                    mbd_north_only[ibin] += weight
                elif has_south[ievt]:
                    mbd_south_only[ibin] += weight
                else:
                    mbd_neither[ibin] += weight

                # Acceptance categories
                if mbd_fires[ievt] and reco_vtx_ok[ievt]:
                    numer_mbd_vtx[ibin] += weight
                if mbd_fires[ievt]:
                    numer_mbd_only[ibin] += weight
                if reco_vtx_ok[ievt]:
                    numer_vtx_only[ibin] += weight

        if max_events > 0 and n_events >= max_events:
            break

    print(f"  [{label}] {n_events} events, {n_truth} truth photons")

    return {
        "denom": denom,
        "numer_mbd_vtx": numer_mbd_vtx,
        "numer_mbd_only": numer_mbd_only,
        "numer_vtx_only": numer_vtx_only,
        "mbd_both": mbd_both,
        "mbd_north_only": mbd_north_only,
        "mbd_south_only": mbd_south_only,
        "mbd_neither": mbd_neither,
        "n_events": n_events,
        "n_truth": n_truth,
    }


def merge_results(results_list):
    """Merge results from multiple samples (weighted addition)."""
    merged = {}
    for key in results_list[0]:
        if isinstance(results_list[0][key], np.ndarray):
            merged[key] = sum(r[key] for r in results_list)
        else:
            merged[key] = sum(r[key] for r in results_list)
    return merged


def compute_efficiencies(result):
    """Compute efficiency arrays from raw counts."""
    d = result["denom"]
    with np.errstate(divide="ignore", invalid="ignore"):
        eff_mbd_vtx = np.where(d > 0, result["numer_mbd_vtx"] / d, 0)
        eff_mbd_only = np.where(d > 0, result["numer_mbd_only"] / d, 0)
        eff_vtx_only = np.where(d > 0, result["numer_vtx_only"] / d, 0)

        frac_both = np.where(d > 0, result["mbd_both"] / d, 0)
        frac_north = np.where(d > 0, result["mbd_north_only"] / d, 0)
        frac_south = np.where(d > 0, result["mbd_south_only"] / d, 0)
        frac_neither = np.where(d > 0, result["mbd_neither"] / d, 0)

        # Binomial errors
        err_mbd_vtx = np.where(d > 0, np.sqrt(eff_mbd_vtx * (1 - eff_mbd_vtx) / d), 0)

    return {
        "eff_mbd_vtx": eff_mbd_vtx,
        "eff_mbd_only": eff_mbd_only,
        "eff_vtx_only": eff_vtx_only,
        "err_mbd_vtx": err_mbd_vtx,
        "frac_both": frac_both,
        "frac_north_only": frac_north,
        "frac_south_only": frac_south,
        "frac_neither": frac_neither,
    }


def save_root(single_eff, double_eff, corrections, outpath):
    """Save results to ROOT file using uproot."""
    pt_centers = 0.5 * (PT_EDGES[:-1] + PT_EDGES[1:])
    pt_widths = 0.5 * (PT_EDGES[1:] - PT_EDGES[:-1])

    data = {}

    # Single-interaction efficiencies
    for key in ["eff_mbd_vtx", "eff_mbd_only", "eff_vtx_only", "err_mbd_vtx",
                "frac_both", "frac_north_only", "frac_south_only", "frac_neither"]:
        data[f"single_{key}"] = single_eff[key].astype(np.float64)
    # Double
    for key in ["eff_mbd_vtx", "eff_mbd_only", "eff_vtx_only", "err_mbd_vtx",
                "frac_both", "frac_north_only", "frac_south_only", "frac_neither"]:
        data[f"double_{key}"] = double_eff[key].astype(np.float64)
    # Corrections
    for period, corr in corrections.items():
        data[f"correction_{period}"] = corr["correction"].astype(np.float64)
        data[f"eff_blended_{period}"] = corr["eff_blended"].astype(np.float64)
    # pT info
    data["pt_centers"] = pt_centers.astype(np.float64)
    data["pt_lo"] = PT_EDGES[:-1].astype(np.float64)
    data["pt_hi"] = PT_EDGES[1:].astype(np.float64)

    with uproot.recreate(outpath) as fout:
        # Save as a TTree for easy access
        fout["results"] = data
        # Also save as individual histograms
        for key in ["eff_mbd_vtx", "eff_mbd_only", "eff_vtx_only"]:
            for tag in ["single", "double"]:
                vals = data[f"{tag}_{key}"]
                # Create histogram-like arrays
                fout[f"h_{tag}_{key}"] = (vals, PT_EDGES)

        for period in corrections:
            fout[f"h_correction_{period}"] = (corrections[period]["correction"], PT_EDGES)

    print(f"Saved ROOT: {outpath}")


def save_csv(single_eff, double_eff, corrections, outpath):
    """Save per-pT results table to CSV."""
    with open(outpath, "w") as fh:
        header = "pt_lo,pt_hi"
        header += ",single_eff_mbd_vtx,single_eff_mbd,single_eff_vtx,single_err"
        header += ",double_eff_mbd_vtx,double_eff_mbd,double_eff_vtx,double_err"
        header += ",single_frac_both,single_frac_N,single_frac_S,single_frac_neither"
        header += ",double_frac_both,double_frac_N,double_frac_S,double_frac_neither"
        for period in corrections:
            header += f",corr_{period},eff_blend_{period}"
        fh.write(header + "\n")

        for i in range(N_PT_BINS):
            row = f"{PT_EDGES[i]:.0f},{PT_EDGES[i+1]:.0f}"
            for tag, eff in [("single", single_eff), ("double", double_eff)]:
                row += f",{eff['eff_mbd_vtx'][i]:.6f},{eff['eff_mbd_only'][i]:.6f}"
                row += f",{eff['eff_vtx_only'][i]:.6f},{eff['err_mbd_vtx'][i]:.6f}"
            for tag, eff in [("single", single_eff), ("double", double_eff)]:
                row += f",{eff['frac_both'][i]:.6f},{eff['frac_north_only'][i]:.6f}"
                row += f",{eff['frac_south_only'][i]:.6f},{eff['frac_neither'][i]:.6f}"
            for period, corr in corrections.items():
                row += f",{corr['correction'][i]:.6f},{corr['eff_blended'][i]:.6f}"
            fh.write(row + "\n")

    print(f"Saved CSV: {outpath}")


def print_summary(single_eff, double_eff, corrections, single_raw, double_raw):
    """Print formatted summary to stdout."""
    print(f"\n{'='*70}")
    print("PILEUP VERTEX+MBD ACCEPTANCE STUDY — RESULTS")
    print(f"{'='*70}")

    print(f"\n--- Integrated vertex_eff (MBD+vtx) ---")
    for tag, eff, raw in [("Single", single_eff, single_raw),
                          ("Double", double_eff, double_raw)]:
        d = raw["denom"].sum()
        n = raw["numer_mbd_vtx"].sum()
        e = n / d if d > 0 else 0
        err = np.sqrt(e * (1 - e) / d) if d > 0 else 0
        print(f"  {tag:8s}: {e:.4f} ± {err:.4f}  (N={d:.0f})")

    print(f"\n  Ratio (double/single): {double_eff['eff_mbd_vtx'].mean() / single_eff['eff_mbd_vtx'].mean():.3f}")

    print(f"\n--- MBD hit decomposition (integrated, %) ---")
    print(f"  {'Category':<20s} {'Single':>10s} {'Double':>10s}")
    for key, label in [("frac_both", "Both N&S"),
                       ("frac_north_only", "N only"),
                       ("frac_south_only", "S only"),
                       ("frac_neither", "Neither")]:
        s = np.average(single_eff[key], weights=single_raw["denom"])
        d = np.average(double_eff[key], weights=double_raw["denom"])
        print(f"  {label:<20s} {100*s:10.1f}% {100*d:10.1f}%")

    print(f"\n--- Per-pT correction factors ---")
    for period, corr in corrections.items():
        f_d = corr["f_double"]
        label = corr["label"]
        print(f"\n  {label} (f_double = {f_d:.1%}):")
        print(f"  {'pT bin':<14s} {'ε_single':>10s} {'ε_double':>10s} {'ε_blend':>10s} {'Correction':>12s}")
        print(f"  {'-'*58}")
        for i in range(N_PT_BINS):
            print(f"  {PT_EDGES[i]:4.0f}-{PT_EDGES[i+1]:4.0f} GeV   "
                  f"{single_eff['eff_mbd_vtx'][i]:10.4f} "
                  f"{double_eff['eff_mbd_vtx'][i]:10.4f} "
                  f"{corr['eff_blended'][i]:10.4f} "
                  f"{corr['correction'][i]:11.4f} ({(corr['correction'][i]-1)*100:+.1f}%)")

        # Integrated correction
        int_blend = np.average(corr["eff_blended"], weights=single_raw["denom"])
        int_single = np.average(single_eff["eff_mbd_vtx"], weights=single_raw["denom"])
        int_corr = int_blend / int_single if int_single > 0 else 1
        print(f"\n  Integrated correction: {int_corr:.4f} ({(int_corr-1)*100:+.1f}%)")
        print(f"  → Cross-section change: {(1-int_corr)*100:+.1f}%")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Section 3 pileup vertex acceptance study")
    parser.add_argument("--max-events", type=int, default=2_000_000)
    parser.add_argument("--results", default="results")
    args = parser.parse_args()

    os.makedirs(args.results, exist_ok=True)

    # --- Process single-interaction samples ---
    print("\n=== Single-interaction MC (cross-section weighted) ===")
    single_results = []
    for name, info in SAMPLES_SINGLE.items():
        r = process_sample(info["path"], info["weight"], args.max_events, name)
        single_results.append(r)
    single_raw = merge_results(single_results)
    single_eff = compute_efficiencies(single_raw)

    # --- Process double-interaction samples ---
    print("\n=== Double-interaction MC ===")
    double_results = []
    for name, info in SAMPLES_DOUBLE.items():
        r = process_sample(info["path"], info["weight"], args.max_events, name)
        double_results.append(r)
    double_raw = merge_results(double_results)
    double_eff = compute_efficiencies(double_raw)

    # --- Compute blended corrections per period ---
    corrections = {}
    for period, pinfo in PILEUP_FRACTIONS.items():
        f_d = pinfo["f_double"]
        f_s = 1.0 - f_d
        eff_blend = f_s * single_eff["eff_mbd_vtx"] + f_d * double_eff["eff_mbd_vtx"]
        correction = np.where(single_eff["eff_mbd_vtx"] > 0,
                              eff_blend / single_eff["eff_mbd_vtx"], 1.0)
        corrections[period] = {
            "f_double": f_d,
            "label": pinfo["label"],
            "eff_blended": eff_blend,
            "correction": correction,
        }

    # --- Output ---
    print_summary(single_eff, double_eff, corrections, single_raw, double_raw)

    root_path = os.path.join(args.results, "pileup_vertex_acceptance.root")
    csv_path = os.path.join(args.results, "pileup_vertex_acceptance.csv")
    save_root(single_eff, double_eff, corrections, root_path)
    save_csv(single_eff, double_eff, corrections, csv_path)


if __name__ == "__main__":
    main()
