#!/usr/bin/env python3
"""Converted vs unconverted photon response study.

Reads photon10 single-particle MC, matches reco clusters to truth photons
by truthtrkID, and fills per-(pT, eta, conversion-category) response and
angular histograms.

Inputs:
  /sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_with_bdt_vtx_noreweight_single.root
  tree=slimtree

Output:
  rootFiles/response_converted.root

Run as:
  python3 analyze_response.py [--max-events N]
"""
from __future__ import annotations

import argparse
import os
import pickle
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

PT_BINS = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
ETA_BINS = np.array([-0.7, -0.35, 0.0, 0.35, 0.7], dtype=float)
R_MIN, R_MAX, R_NBIN = 0.0, 1.3, 130

CATEGORIES = [("unconv", 0), ("conv", 1), ("bad", 2)]
CAT_NAMES = [c[0] for c in CATEGORIES]

# Selection window (slightly inside sample range)
PT_LO, PT_HI = 10.0, 34.0
ETA_ABS_MAX = 0.7
PID_PHOTON = 22
MATCH_ET_MIN = 3.0  # GeV, threshold for "found" match

CLUSTER_SUFFIX = "_CLUSTERINFO_CEMC"

# Branches to read
BRANCHES = [
    "nparticles",
    "particle_pid",
    "particle_trkid",
    "particle_Pt",
    "particle_E",
    "particle_Eta",
    "particle_Phi",
    "particle_converted",
    f"cluster_E{CLUSTER_SUFFIX}",
    f"cluster_Et{CLUSTER_SUFFIX}",
    f"cluster_Eta{CLUSTER_SUFFIX}",
    f"cluster_Phi{CLUSTER_SUFFIX}",
    f"cluster_truthtrkID{CLUSTER_SUFFIX}",
]


def dphi_wrap(dphi: np.ndarray) -> np.ndarray:
    """Wrap phi difference to [-pi, pi]."""
    return (dphi + np.pi) % (2.0 * np.pi) - np.pi


class HistStore:
    """Minimal numpy-based histograms we'll convert to TH1F at the end.

    Keys:
      R_E[cat][pT_idx][eta_idx]   # per-(pT,eta) response (E)
      R_ET[cat][pT_idx][eta_idx]
      R_E_inclusive[cat][eta_idx]  # inclusive pT
      R_ET_inclusive[cat][eta_idx]
      dEta[cat], dPhi[cat], dR[cat]  # inclusive
      match_count_num[cat][pT_idx][eta_idx]  # matched truth photons
      match_count_den[cat][pT_idx][eta_idx]  # all truth photons in selection
    """

    def __init__(self):
        self.nPt = len(PT_BINS) - 1
        self.nEta = len(ETA_BINS) - 1
        self.R_E = {c: np.zeros((self.nPt, self.nEta, R_NBIN)) for c in CAT_NAMES}
        self.R_ET = {c: np.zeros((self.nPt, self.nEta, R_NBIN)) for c in CAT_NAMES}
        self.R_E_pt = {c: np.zeros((self.nPt, R_NBIN)) for c in CAT_NAMES}  # eta-inclusive per-pT
        self.R_ET_pt = {c: np.zeros((self.nPt, R_NBIN)) for c in CAT_NAMES}
        # Inclusive (all pT, all eta)
        self.R_E_inc = {c: np.zeros(R_NBIN) for c in CAT_NAMES}
        self.R_ET_inc = {c: np.zeros(R_NBIN) for c in CAT_NAMES}
        # Per-eta inclusive pT
        self.R_E_eta = {c: np.zeros((self.nEta, R_NBIN)) for c in CAT_NAMES}

        # Angular
        self.dEta = {c: np.zeros(200) for c in CAT_NAMES}
        self.dPhi = {c: np.zeros(200) for c in CAT_NAMES}
        self.dR = {c: np.zeros(200) for c in CAT_NAMES}

        # Match counts
        self.num = {c: np.zeros((self.nPt, self.nEta)) for c in CAT_NAMES}
        self.den = {c: np.zeros((self.nPt, self.nEta)) for c in CAT_NAMES}

        # For later Gaussian fits, stash raw (R, pt, eta, cat) samples up to cap
        self.samples = {c: [] for c in CAT_NAMES}
        self.MAX_SAMPLES_PER_CAT = 800_000

    def fill_truth(self, cat_name: str, pt_idx: np.ndarray, eta_idx: np.ndarray):
        """Fill denominator (all truth photons in selection) per bin."""
        valid = (pt_idx >= 0) & (eta_idx >= 0)
        if not np.any(valid):
            return
        np.add.at(self.den[cat_name], (pt_idx[valid], eta_idx[valid]), 1)

    def fill_matched(
        self,
        cat_name: str,
        pt_idx: np.ndarray,
        eta_idx: np.ndarray,
        R_E: np.ndarray,
        R_ET: np.ndarray,
        dEta: np.ndarray,
        dPhi: np.ndarray,
    ):
        valid = (pt_idx >= 0) & (eta_idx >= 0)
        pt_idx = pt_idx[valid]
        eta_idx = eta_idx[valid]
        R_E = R_E[valid]
        R_ET = R_ET[valid]
        dEta = dEta[valid]
        dPhi = dPhi[valid]
        if len(pt_idx) == 0:
            return

        # Match-count numerator
        np.add.at(self.num[cat_name], (pt_idx, eta_idx), 1)

        # Bin indices for response histograms
        b_E = np.clip(((R_E - R_MIN) / (R_MAX - R_MIN) * R_NBIN).astype(int), 0, R_NBIN - 1)
        b_ET = np.clip(((R_ET - R_MIN) / (R_MAX - R_MIN) * R_NBIN).astype(int), 0, R_NBIN - 1)
        np.add.at(self.R_E[cat_name], (pt_idx, eta_idx, b_E), 1)
        np.add.at(self.R_ET[cat_name], (pt_idx, eta_idx, b_ET), 1)
        np.add.at(self.R_E_pt[cat_name], (pt_idx, b_E), 1)
        np.add.at(self.R_ET_pt[cat_name], (pt_idx, b_ET), 1)
        np.add.at(self.R_E_eta[cat_name], (eta_idx, b_E), 1)
        self.R_E_inc[cat_name] += np.bincount(b_E, minlength=R_NBIN)
        self.R_ET_inc[cat_name] += np.bincount(b_ET, minlength=R_NBIN)

        # Angular
        de_bin = np.clip(((dEta + 0.1) / 0.2 * 200).astype(int), 0, 199)
        dp_bin = np.clip(((dPhi + 0.1) / 0.2 * 200).astype(int), 0, 199)
        dr_bin = np.clip((np.hypot(dEta, dPhi) / 0.2 * 200).astype(int), 0, 199)
        self.dEta[cat_name] += np.bincount(de_bin, minlength=200)
        self.dPhi[cat_name] += np.bincount(dp_bin, minlength=200)
        self.dR[cat_name] += np.bincount(dr_bin, minlength=200)

        # Stash samples for fits
        cur = (
            np.sum([a.shape[0] for a in self.samples[cat_name]])
            if isinstance(self.samples[cat_name], list)
            else 0
        )
        remaining = self.MAX_SAMPLES_PER_CAT - cur
        if remaining > 0:
            n_add = min(remaining, len(pt_idx))
            self.samples[cat_name].append(
                np.stack(
                    [R_E[:n_add], R_ET[:n_add], pt_idx[:n_add].astype(float), eta_idx[:n_add].astype(float)],
                    axis=1,
                )
            )

    def finalize_samples(self):
        for c in CAT_NAMES:
            if self.samples[c]:
                self.samples[c] = np.vstack(self.samples[c])
            else:
                self.samples[c] = np.zeros((0, 4))


def process(max_events: int | None = None):
    store = HistStore()
    total = 0
    events_seen = 0
    with uproot.open(INPUT_FILE) as f:
        t = f[TREE]
        entries = t.num_entries
        if max_events is not None:
            entries = min(entries, max_events)
        print(f"[analyze] processing {entries} entries (total file {t.num_entries})")
        for batch in t.iterate(
            BRANCHES,
            step_size=500_000,
            entry_stop=entries,
        ):
            n_batch = len(batch)
            events_seen += n_batch
            # Truth photon selection
            pid = batch["particle_pid"]
            pt = batch["particle_Pt"]
            eta = batch["particle_Eta"]
            phi = batch["particle_Phi"]
            energy = batch["particle_E"]
            conv = batch["particle_converted"]
            trkid = batch["particle_trkid"]

            sel = (
                (pid == PID_PHOTON)
                & (pt >= PT_LO)
                & (pt < PT_HI)
                & (np.abs(eta) < ETA_ABS_MAX)
            )
            n_truth_photons = ak.sum(sel)
            if n_truth_photons == 0:
                continue
            total += int(n_truth_photons)

            # Vectorized matching via ak.cartesian: for each truth photon in each event,
            # find the highest-ET cluster whose truthtrkID equals the photon's trkid.
            ptrk_evt = trkid[sel]  # jagged (evt, photon)
            ppt_evt = pt[sel]
            peta_evt = eta[sel]
            pphi_evt = phi[sel]
            pE_evt = energy[sel]
            pconv_evt = conv[sel]

            c_E = batch[f"cluster_E{CLUSTER_SUFFIX}"]
            c_ET = batch[f"cluster_Et{CLUSTER_SUFFIX}"]
            c_eta = batch[f"cluster_Eta{CLUSTER_SUFFIX}"]
            c_phi = batch[f"cluster_Phi{CLUSTER_SUFFIX}"]
            c_trk = batch[f"cluster_truthtrkID{CLUSTER_SUFFIX}"]

            # pairs: evt -> photon -> cluster (nested cartesian)
            pairs = ak.cartesian({"p": ptrk_evt, "c": c_trk}, axis=1, nested=True)
            match_mask = pairs["p"] == pairs["c"]
            # Broadcast cluster observables to (evt, photon, cluster)
            c_ET_b = ak.broadcast_arrays(ptrk_evt[..., None], c_ET[:, None, ...])[1]
            c_E_b = ak.broadcast_arrays(ptrk_evt[..., None], c_E[:, None, ...])[1]
            c_eta_b = ak.broadcast_arrays(ptrk_evt[..., None], c_eta[:, None, ...])[1]
            c_phi_b = ak.broadcast_arrays(ptrk_evt[..., None], c_phi[:, None, ...])[1]

            NEG_INF = -1.0e30
            c_ET_masked = ak.where(match_mask, c_ET_b, NEG_INF)
            best_idx = ak.argmax(c_ET_masked, axis=2, keepdims=True)

            best_ET = ak.firsts(c_ET_masked[best_idx], axis=2)
            best_E = ak.firsts(ak.where(match_mask, c_E_b, 0.0)[best_idx], axis=2)
            best_eta = ak.firsts(ak.where(match_mask, c_eta_b, 0.0)[best_idx], axis=2)
            best_phi = ak.firsts(ak.where(match_mask, c_phi_b, 0.0)[best_idx], axis=2)

            # Flatten to per-photon 1D arrays
            photon_pt = ak.flatten(ppt_evt).to_numpy()
            photon_eta = ak.flatten(peta_evt).to_numpy()
            photon_phi = ak.flatten(pphi_evt).to_numpy()
            photon_E = ak.flatten(pE_evt).to_numpy()
            photon_conv = ak.flatten(pconv_evt).to_numpy()

            matched_ET = ak.to_numpy(ak.fill_none(ak.flatten(best_ET), NEG_INF))
            matched_E = ak.to_numpy(ak.fill_none(ak.flatten(best_E), 0.0))
            matched_eta = ak.to_numpy(ak.fill_none(ak.flatten(best_eta), 0.0))
            matched_phi = ak.to_numpy(ak.fill_none(ak.flatten(best_phi), 0.0))
            # Apply ET threshold
            match_valid = matched_ET > MATCH_ET_MIN
            matched_E = np.where(match_valid, matched_E, np.nan)
            matched_ET = np.where(match_valid, matched_ET, np.nan)
            matched_eta = np.where(match_valid, matched_eta, np.nan)
            matched_phi = np.where(match_valid, matched_phi, np.nan)

            # Digitize pT and eta into bins
            pt_idx_all = np.digitize(photon_pt, PT_BINS) - 1
            eta_idx_all = np.digitize(photon_eta, ETA_BINS) - 1
            # Clip outside bins to -1 for "invalid"
            pt_idx_all = np.where((pt_idx_all < 0) | (pt_idx_all >= len(PT_BINS) - 1), -1, pt_idx_all)
            eta_idx_all = np.where((eta_idx_all < 0) | (eta_idx_all >= len(ETA_BINS) - 1), -1, eta_idx_all)

            # Fill denominator per category
            for cat_name, cat_val in CATEGORIES:
                msk = photon_conv == cat_val
                if not np.any(msk):
                    continue
                store.fill_truth(cat_name, pt_idx_all[msk], eta_idx_all[msk])

            # Compute responses (nan safe)
            R_E = matched_E / photon_E
            R_ET = matched_ET / photon_pt
            dEta = matched_eta - photon_eta
            dPhi = dphi_wrap(matched_phi - photon_phi)

            # Matched mask
            matched_mask = ~np.isnan(matched_E)

            for cat_name, cat_val in CATEGORIES:
                cmsk = (photon_conv == cat_val) & matched_mask
                if not np.any(cmsk):
                    continue
                store.fill_matched(
                    cat_name,
                    pt_idx_all[cmsk],
                    eta_idx_all[cmsk],
                    R_E[cmsk],
                    R_ET[cmsk],
                    dEta[cmsk],
                    dPhi[cmsk],
                )
            if events_seen % 1_000_000 < 500_000:
                print(f"  events_seen={events_seen}, total truth photons so far={total}")

    store.finalize_samples()
    print(f"[analyze] done. events_seen={events_seen}, truth photons in selection={total}")
    return store, events_seen


def save_store(store: HistStore, events_seen: int, outpath: Path):
    import uproot

    # Prepare TH1 objects manually via uproot.writing
    # We'll write TH1D using uproot 5's histogram interface.
    # Use edges for R histograms.
    R_edges = np.linspace(R_MIN, R_MAX, R_NBIN + 1)
    de_edges = np.linspace(-0.1, 0.1, 201)
    dp_edges = np.linspace(-0.1, 0.1, 201)
    dr_edges = np.linspace(0.0, 0.2, 201)

    out = {}

    # Per-(pT, eta, cat)
    for c in CAT_NAMES:
        # Per-pT, eta-inclusive
        for ip in range(store.nPt):
            counts = store.R_E_pt[c][ip]
            out[f"R_E_{c}_pt{ip}"] = (counts, R_edges)
            out[f"R_ET_{c}_pt{ip}"] = (store.R_ET_pt[c][ip], R_edges)
        # Per-eta, pT-inclusive
        for ie in range(store.nEta):
            out[f"R_E_{c}_eta{ie}"] = (store.R_E_eta[c][ie], R_edges)
        # 2D (pT x eta)
        for ip in range(store.nPt):
            for ie in range(store.nEta):
                out[f"R_E_{c}_pt{ip}_eta{ie}"] = (store.R_E[c][ip, ie], R_edges)
                out[f"R_ET_{c}_pt{ip}_eta{ie}"] = (store.R_ET[c][ip, ie], R_edges)
        # inclusive
        out[f"R_E_{c}_inc"] = (store.R_E_inc[c], R_edges)
        out[f"R_ET_{c}_inc"] = (store.R_ET_inc[c], R_edges)
        out[f"dEta_{c}"] = (store.dEta[c], de_edges)
        out[f"dPhi_{c}"] = (store.dPhi[c], dp_edges)
        out[f"dR_{c}"] = (store.dR[c], dr_edges)

    # Match counts as TH2 (pT, eta)
    pt_edges = PT_BINS
    eta_edges = ETA_BINS
    for c in CAT_NAMES:
        out[f"num_{c}"] = (store.num[c], pt_edges, eta_edges)
        out[f"den_{c}"] = (store.den[c], pt_edges, eta_edges)

    with uproot.recreate(outpath) as fout:
        for name, payload in out.items():
            if len(payload) == 2:
                counts, edges = payload
                fout[name] = (counts.astype(np.float64), edges)
            else:
                counts, xe, ye = payload
                fout[name] = (counts.astype(np.float64), xe, ye)

    # Also save samples pickle for fitting
    sample_path = outpath.parent / "samples.pkl"
    with open(sample_path, "wb") as sf:
        pickle.dump({"samples": store.samples, "events_seen": events_seen}, sf)
    print(f"[analyze] wrote {outpath} and {sample_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-events", type=int, default=None, help="cap on entries read")
    args = ap.parse_args()
    store, events_seen = process(max_events=args.max_events)
    save_store(store, events_seen, OUTDIR / "response_converted.root")


if __name__ == "__main__":
    main()
