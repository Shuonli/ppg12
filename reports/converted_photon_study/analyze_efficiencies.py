#!/usr/bin/env python3
"""Three-stage selection efficiency for converted vs unconverted truth photons.

For each truth photon passing the fiducial selection (|eta|<0.7, truth_iso_03<4 GeV,
8 <= pT < 36 GeV) this script records, per conversion category:
  - reco:  at least one cluster on CLUSTERINFO_CEMC with cluster_truthtrkID == trkid
           and cluster_Et > MATCH_ET_MIN
  - ID:    at least one cluster on CLUSTERINFO_CEMC_NO_SPLIT with matching trkid
           whose BDT exceeds 0.80 - 0.015 * cluster_Et (evaluated per its own ET)
  - iso:   CEMC-matched cluster above the ET threshold with
           cluster_iso_03 < 0.502095 + 0.0433036 * cluster_Et
  - total: reco AND ID AND iso (same truth photon passing all three cuts)

Binned in truth pT with plotcommon.h ptRanges = [8,10,12,14,16,18,20,22,24,26,28,32,36].

Outputs:
  rootFiles/efficiencies_converted.root  TH1D per (cat, stage) plus total N_truth
  rootFiles/efficiencies_summary.json    per-pT + inclusive efficiencies with
                                          binomial errors
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import awkward as ak
import numpy as np
import uproot
from tqdm import tqdm

HERE = Path(__file__).resolve().parent
ROOT_OUT = HERE / "rootFiles" / "efficiencies_converted.root"
JSON_OUT = HERE / "rootFiles" / "efficiencies_summary.json"

INPUT = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_with_bdt_vtx_noreweight_single.root"
TREE = "slimtree"

PT_EDGES = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
N_PT = len(PT_EDGES) - 1

ETA_MAX = 0.7
TRUTH_ISO_MAX = 4.0
MATCH_ET_MIN = 8.0  # reco-stage threshold (matches analysis pT range start)

# Parametric cuts (efficiencytool/config_bdt_nom.yaml)
BDT_INTERCEPT = 0.80
BDT_SLOPE = -0.015
ISO_INTERCEPT = 0.502095
ISO_SLOPE = 0.0433036

# Conversion-category codes: 0=unconv, 1=conv, 2=bad
CATS = {0: "unconv", 1: "conv", 2: "bad"}

TRUTH_BRANCHES = [
    "particle_pid", "particle_Pt", "particle_Eta",
    "particle_E", "particle_trkid", "particle_converted",
    "particle_truth_iso_03",
]
CEMC_BRANCHES = [
    "cluster_Et_CLUSTERINFO_CEMC",
    "cluster_truthtrkID_CLUSTERINFO_CEMC",
    "cluster_iso_03_CLUSTERINFO_CEMC",
]
NOSPLIT_BRANCHES = [
    "cluster_Et_CLUSTERINFO_CEMC_NO_SPLIT",
    "cluster_truthtrkID_CLUSTERINFO_CEMC_NO_SPLIT",
    "cluster_bdt_CLUSTERINFO_CEMC_NO_SPLIT",
]


def pt_index(pt: np.ndarray) -> np.ndarray:
    """Return pT-bin index [0..N_PT-1] for each entry, or -1 if out of range."""
    idx = np.digitize(pt, PT_EDGES) - 1
    idx = np.where((idx >= 0) & (idx < N_PT), idx, -1)
    return idx


def process_file(max_events: int | None = None, step_size: int = 200_000):
    # Counters per (category, pt_bin):
    #   n_truth       = fiducial truth photons
    #   n_reco        = CEMC-matched (iso denominator)
    #   n_reco_ns     = NO_SPLIT-matched (BDT denominator — different node)
    #   n_reco_both   = BOTH CEMC and NO_SPLIT matched (joint denominator for ID+iso)
    #   n_id          = NO_SPLIT-matched and BDT passes
    #   n_iso         = CEMC-matched and iso passes
    #   n_all         = both matched AND BDT AND iso
    counters = {
        c: {k: np.zeros(N_PT, dtype=np.int64)
            for k in ("n_truth", "n_reco", "n_reco_ns", "n_reco_both",
                       "n_id", "n_iso", "n_all")}
        for c in CATS.values()
    }

    total_processed = 0
    branches = (
        TRUTH_BRANCHES + CEMC_BRANCHES + NOSPLIT_BRANCHES
    )

    with uproot.open(INPUT) as f:
        tree = f[TREE]
        n_total = tree.num_entries
        cap = n_total if max_events is None else min(n_total, max_events)
        pbar = tqdm(total=cap, unit="evt")
        for arrays in tree.iterate(branches, step_size=step_size, library="ak"):
            if total_processed >= cap:
                break
            n_batch = len(arrays)
            if total_processed + n_batch > cap:
                arrays = arrays[: cap - total_processed]
                n_batch = len(arrays)

            # Flatten per-cluster arrays to numpy for matching
            t_pid = arrays["particle_pid"]
            t_pt = arrays["particle_Pt"]
            t_eta = arrays["particle_Eta"]
            t_iso = arrays["particle_truth_iso_03"]
            t_trk = arrays["particle_trkid"]
            t_conv = arrays["particle_converted"]

            # fiducial mask at truth level (still nested per event)
            fid = (
                (t_pid == 22)
                & (np.abs(t_eta) < ETA_MAX)
                & (t_iso < TRUTH_ISO_MAX)
                & (t_pt >= PT_EDGES[0])
                & (t_pt < PT_EDGES[-1])
            )

            # Cluster collections (per event, flat arrays)
            cemc_et = arrays["cluster_Et_CLUSTERINFO_CEMC"]
            cemc_trk = arrays["cluster_truthtrkID_CLUSTERINFO_CEMC"]
            cemc_iso03 = arrays["cluster_iso_03_CLUSTERINFO_CEMC"]

            nos_et = arrays["cluster_Et_CLUSTERINFO_CEMC_NO_SPLIT"]
            nos_trk = arrays["cluster_truthtrkID_CLUSTERINFO_CEMC_NO_SPLIT"]
            nos_bdt = arrays["cluster_bdt_CLUSTERINFO_CEMC_NO_SPLIT"]

            n_evt = len(t_pt)
            for ievt in range(n_evt):
                fid_mask = np.asarray(fid[ievt])
                if not np.any(fid_mask):
                    continue
                trk_arr = np.asarray(t_trk[ievt])
                pt_arr = np.asarray(t_pt[ievt])
                conv_arr = np.asarray(t_conv[ievt])

                # Per-event cluster arrays
                ce_et = np.asarray(cemc_et[ievt])
                ce_trk = np.asarray(cemc_trk[ievt])
                ce_iso = np.asarray(cemc_iso03[ievt])
                ns_et = np.asarray(nos_et[ievt])
                ns_trk = np.asarray(nos_trk[ievt])
                ns_bdt = np.asarray(nos_bdt[ievt])

                fid_idx = np.where(fid_mask)[0]
                for ip in fid_idx:
                    pt = pt_arr[ip]
                    ipt = int(np.digitize(pt, PT_EDGES) - 1)
                    if ipt < 0 or ipt >= N_PT:
                        continue
                    code = int(conv_arr[ip])
                    cat = CATS.get(code)
                    if cat is None:
                        continue
                    trk = int(trk_arr[ip])

                    # CEMC match: highest-ET cluster with matching trkid and ET>MATCH
                    ce_match_mask = (ce_trk == trk) & (ce_et > MATCH_ET_MIN)
                    cemc_ok = bool(np.any(ce_match_mask))
                    iso_ok = False
                    if cemc_ok:
                        ce_cand = np.where(ce_match_mask)[0]
                        best = ce_cand[np.argmax(ce_et[ce_cand])]
                        iso_thr = ISO_INTERCEPT + ISO_SLOPE * ce_et[best]
                        iso_ok = bool(ce_iso[best] < iso_thr)

                    # NO_SPLIT match: highest-ET cluster with matching trkid and ET>MATCH
                    ns_match_mask = (ns_trk == trk) & (ns_et > MATCH_ET_MIN)
                    ns_ok = bool(np.any(ns_match_mask))
                    id_ok = False
                    if ns_ok:
                        ns_cand = np.where(ns_match_mask)[0]
                        best_ns = ns_cand[np.argmax(ns_et[ns_cand])]
                        bdt_thr = BDT_INTERCEPT + BDT_SLOPE * ns_et[best_ns]
                        id_ok = bool(ns_bdt[best_ns] > bdt_thr)

                    counters[cat]["n_truth"][ipt] += 1
                    if cemc_ok:
                        counters[cat]["n_reco"][ipt] += 1
                    if ns_ok:
                        counters[cat]["n_reco_ns"][ipt] += 1
                    if cemc_ok and ns_ok:
                        counters[cat]["n_reco_both"][ipt] += 1
                    if id_ok:
                        counters[cat]["n_id"][ipt] += 1
                    if iso_ok:
                        counters[cat]["n_iso"][ipt] += 1
                    if cemc_ok and id_ok and iso_ok:
                        counters[cat]["n_all"][ipt] += 1

            total_processed += n_batch
            pbar.update(n_batch)
        pbar.close()

    return counters, total_processed


def binomial_err(k: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Standard clopper-like binomial approximation with n=0 guard."""
    n_safe = np.where(n > 0, n, 1)
    p = k / n_safe
    err = np.sqrt(np.maximum(p * (1 - p), 0.0) / n_safe)
    return np.where(n > 0, err, 0.0)


def build_summary(counters):
    out = {"pt_edges": PT_EDGES.tolist()}
    inclusive = {}
    per_pt = {}
    for cat, c in counters.items():
        n_truth = c["n_truth"]
        n_reco = c["n_reco"]
        n_reco_ns = c["n_reco_ns"]
        n_id = c["n_id"]
        n_iso = c["n_iso"]
        n_all = c["n_all"]

        with np.errstate(invalid="ignore", divide="ignore"):
            eps_reco = np.where(n_truth > 0, n_reco / n_truth, np.nan)
            eps_id = np.where(n_truth > 0, n_id / n_truth, np.nan)
            eps_iso = np.where(n_truth > 0, n_iso / n_truth, np.nan)
            eps_all = np.where(n_truth > 0, n_all / n_truth, np.nan)
            # Conditional: use matching denominator (same cluster node)
            eps_id_cond = np.where(n_reco_ns > 0, n_id / n_reco_ns, np.nan)
            eps_iso_cond = np.where(n_reco > 0, n_iso / n_reco, np.nan)

        per_pt[cat] = {
            "n_truth": n_truth.tolist(),
            "n_reco": n_reco.tolist(),
            "n_reco_ns": n_reco_ns.tolist(),
            "n_id": n_id.tolist(),
            "n_iso": n_iso.tolist(),
            "n_all": n_all.tolist(),
            "eps_reco": eps_reco.tolist(),
            "eps_reco_err": binomial_err(n_reco, n_truth).tolist(),
            "eps_id": eps_id.tolist(),
            "eps_id_err": binomial_err(n_id, n_truth).tolist(),
            "eps_iso": eps_iso.tolist(),
            "eps_iso_err": binomial_err(n_iso, n_truth).tolist(),
            "eps_all": eps_all.tolist(),
            "eps_all_err": binomial_err(n_all, n_truth).tolist(),
            "eps_id_cond": eps_id_cond.tolist(),
            "eps_id_cond_err": binomial_err(n_id, n_reco_ns).tolist(),
            "eps_iso_cond": eps_iso_cond.tolist(),
            "eps_iso_cond_err": binomial_err(n_iso, n_reco).tolist(),
        }

        N = n_truth.sum()
        R = n_reco.sum()
        Rns = n_reco_ns.sum()
        I = n_id.sum()
        Is = n_iso.sum()
        A = n_all.sum()
        safe = N if N > 0 else 1
        inclusive[cat] = {
            "n_truth": int(N), "n_reco": int(R), "n_reco_ns": int(Rns),
            "n_id": int(I), "n_iso": int(Is), "n_all": int(A),
            "eps_reco": R / safe, "eps_id": I / safe,
            "eps_iso": Is / safe, "eps_all": A / safe,
            "eps_id_cond": (I / Rns) if Rns > 0 else None,
            "eps_iso_cond": (Is / R) if R > 0 else None,
        }
    out["per_pt"] = per_pt
    out["inclusive"] = inclusive
    return out


def write_root(counters):
    # Write TH1D-like arrays so downstream tools can read back
    ROOT_OUT.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(ROOT_OUT) as f:
        for cat, c in counters.items():
            for k, arr in c.items():
                f[f"{k}_{cat}"] = (arr.astype(np.float64),
                                   PT_EDGES.astype(np.float64))


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-events", type=int, default=None,
                        help="Cap number of events processed (default: full file)")
    parser.add_argument("--step-size", type=int, default=200_000)
    args = parser.parse_args()

    counters, n_processed = process_file(
        max_events=args.max_events, step_size=args.step_size)
    print(f"[efficiencies] processed {n_processed:,} events")

    write_root(counters)
    print(f"[efficiencies] wrote {ROOT_OUT}")

    summary = build_summary(counters)
    with open(JSON_OUT, "w") as jf:
        json.dump(summary, jf, indent=2)
    print(f"[efficiencies] wrote {JSON_OUT}")

    # Short text summary
    for cat in ("unconv", "conv"):
        inc = summary["inclusive"][cat]
        print(f"  {cat}: N={inc['n_truth']:,} "
              f"eps_reco={inc['eps_reco']:.4f} "
              f"eps_id={inc['eps_id']:.4f} "
              f"eps_iso={inc['eps_iso']:.4f} "
              f"eps_all={inc['eps_all']:.4f} "
              f"eps_id|reco={inc['eps_id_cond']:.4f} "
              f"eps_iso|reco={inc['eps_iso_cond']:.4f}")


if __name__ == "__main__":
    main()
