#!/usr/bin/env python3
"""
Independent re-derivation of f_direct = N_direct / (N_direct + N_frag) from raw
signal MC for PPG12 direct-vs-frag investigation.

This script is written from scratch and does NOT re-use any logic from
extract_direct_vs_frag.py. It reads the three photon-signal MC files,
applies fiducial / truth-class / truth-iso / event-pT-window cuts and
per-sample cross-section weighting, and prints f_direct per pT bin at
iso < 4 GeV for comparison against claimed values.

Claimed values (iso < 4 GeV):
  [ 8,14) -> 0.7456
  [14,22) -> 0.7933
  [22,36) -> 0.8181
"""

import os
import sys
import numpy as np
import uproot
import awkward as ak

# --- Inputs -----------------------------------------------------------------
BASE = "/sphenix/user/shuhangli/ppg12/FunWithxgboost"

# (sample_name, file_path, pT window on max photon pT in event, cross-section)
SAMPLES = [
    ("photon5",  f"{BASE}/photon5/bdt_split.root",  (0.0,  14.0), 146359.3),
    ("photon10", f"{BASE}/photon10/bdt_split.root", (14.0, 30.0),   6944.675),
    ("photon20", f"{BASE}/photon20/bdt_split.root", (30.0, 200.0),   130.4461),
]
REF_XSEC = 130.4461  # photon20cross; weights are per-sample xsec / this

TREE = "slimtree"

# --- Selection / binning ----------------------------------------------------
ETA_MAX = 0.7
PT_MIN, PT_MAX = 8.0, 36.0
CLASS_DIRECT, CLASS_FRAG = 1, 2
ISO_CUT = 4.0  # GeV

PT_BINS = [(8.0, 14.0), (14.0, 22.0), (22.0, 36.0)]

BRANCHES = [
    "particle_pid",
    "particle_Pt",
    "particle_Eta",
    "particle_photonclass",
    "particle_truth_iso_03",
]


def process_sample(path, pt_window, weight):
    """Return (pt_flat, cls_flat, iso_flat, w_flat) for this sample."""
    pt_lo, pt_hi = pt_window
    pt_list, cls_list, iso_list, w_list = [], [], [], []

    with uproot.open(path) as f:
        t = f[TREE]
        nev_total = t.num_entries
        nev_keep = 0
        for arrays in t.iterate(BRANCHES, step_size="200 MB", library="ak"):
            pid   = arrays["particle_pid"]
            pt    = arrays["particle_Pt"]
            eta   = arrays["particle_Eta"]
            cls   = arrays["particle_photonclass"]
            iso   = arrays["particle_truth_iso_03"]

            # Event-level max photon pT (among pid==22). Events with no photon
            # get -inf so they fail the window, except pt_lo=0 for photon5
            # which still fails because no photon can pass fiducial anyway.
            is_photon = (pid == 22)
            pt_photon = ak.mask(pt, is_photon)  # non-photons become None
            max_pt = ak.fill_none(ak.max(pt_photon, axis=1), -1.0)
            ev_keep = (max_pt >= pt_lo) & (max_pt < pt_hi)
            nev_keep += int(ak.sum(ev_keep))

            # Apply event mask by subsetting all arrays
            pid   = pid[ev_keep]
            pt    = pt[ev_keep]
            eta   = eta[ev_keep]
            cls   = cls[ev_keep]
            iso   = iso[ev_keep]

            # Particle-level selection: pid==22, |eta|<0.7, 8<=pt<36,
            # class in {1,2}. iso cut applied later when we bin f_direct.
            m = (pid == 22) & (abs(eta) < ETA_MAX) \
                & (pt >= PT_MIN) & (pt < PT_MAX) \
                & ((cls == CLASS_DIRECT) | (cls == CLASS_FRAG))

            pt_sel  = ak.flatten(pt[m])
            cls_sel = ak.flatten(cls[m])
            iso_sel = ak.flatten(iso[m])
            pt_list.append(ak.to_numpy(pt_sel))
            cls_list.append(ak.to_numpy(cls_sel))
            iso_list.append(ak.to_numpy(iso_sel))
            w_list.append(np.full(len(pt_sel), weight, dtype=np.float64))

        print(f"  [{os.path.basename(os.path.dirname(path))}] "
              f"events total={nev_total}, kept in window={nev_keep}, "
              f"selected particles={sum(len(x) for x in pt_list)}")

    if not pt_list:
        return (np.array([]),) * 4
    return (np.concatenate(pt_list),
            np.concatenate(cls_list),
            np.concatenate(iso_list),
            np.concatenate(w_list))


def main():
    all_pt, all_cls, all_iso, all_w = [], [], [], []
    for name, path, win, xsec in SAMPLES:
        w = xsec / REF_XSEC
        print(f"Processing {name} (xsec={xsec}, rel weight={w:.6g}, "
              f"window=[{win[0]},{win[1]}))")
        pt, cls, iso, ww = process_sample(path, win, w)
        all_pt.append(pt); all_cls.append(cls); all_iso.append(iso); all_w.append(ww)

    pt  = np.concatenate(all_pt)
    cls = np.concatenate(all_cls)
    iso = np.concatenate(all_iso)
    w   = np.concatenate(all_w)

    print(f"\nTotal selected particles across samples: {len(pt)}")

    iso_mask = iso < ISO_CUT
    claimed = {(8.0,14.0): 0.7456, (14.0,22.0): 0.7933, (22.0,36.0): 0.8181}

    print("\npT bin       N_direct(w)     N_frag(w)     f_direct    claimed   diff   status")
    print("-" * 88)
    any_fail = False
    for lo, hi in PT_BINS:
        m_bin = iso_mask & (pt >= lo) & (pt < hi)
        w_bin  = w[m_bin]
        cls_bin = cls[m_bin]
        nd = float(np.sum(w_bin[cls_bin == CLASS_DIRECT]))
        nf = float(np.sum(w_bin[cls_bin == CLASS_FRAG]))
        denom = nd + nf
        f = nd / denom if denom > 0 else float("nan")
        c = claimed[(lo, hi)]
        diff = abs(f - c)
        status = "PASS" if diff < 0.01 * c else "FAIL"
        if status == "FAIL":
            any_fail = True
        print(f"[{lo:>4.0f},{hi:>4.0f})  {nd:12.2f}  {nf:12.2f}   "
              f"{f:8.4f}   {c:7.4f}  {diff:7.4f}  {status}")

    print("-" * 88)
    print("Overall:", "PASS" if not any_fail else "FAIL")


if __name__ == "__main__":
    main()
