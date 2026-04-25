#!/usr/bin/env python3
"""
Topo-cluster isolation + photon self-energy containment study.

Study A: isolation distributions vs cone radius.
  - Truth-matched isolated direct photons (signal MC only)      [for containment physics]
  - Inclusive MC (signal + jets, cross-section weighted)        [matches what is compared to data]
  - Data (part_*_with_bdt_split.root, no truth match)
  Each with pre-common and post-common (npb_score > 0.5) selections.

Study A': tower-based iso scan — same samples, plain tower branches
  cluster_iso_{R}_{NODE} for R in {0.05, 0.075, 0.20, 0.30, 0.40}.

Study A'': NEW tower iso variants (added Apr 2026)
  - per-calo split at R=0.3, R=0.4 (no threshold)
  - per-tower threshold scan at R=0.3 (60 / 70 / 120 MeV, per calo)
  - EMCal-only 70 MeV at R in {0.05, 0.075, 0.10, 0.20, 0.30}
  - SUB1-subtracted iso at R=0.3, R=0.4 (per calo)
  - ownership-excluded iso (cluster_iso_excl_R) at R in {0.05, 0.075, 0.10, 0.20, 0.30, 0.40}

Study B: photon self-energy containment — TRUTH-MATCHED SIGNAL MC ONLY
  f_contain(R) = (cluster_iso_topo_R + cluster_ET) / particle_truth_Pt
  Both pre-common and post-common.

Output: results/iso_topo_containment.npz
"""
import argparse
import numpy as np
import uproot
from pathlib import Path

NODE = "CLUSTERINFO_CEMC"
RADII      = [("005", 0.05), ("0075", 0.075), ("01", 0.1), ("02", 0.2), ("03", 0.3), ("04", 0.4)]
# Plain tower cluster_iso_R branches: no R=0.1 in the current slimtree
RADII_TOWER = [("005", 0.05), ("0075", 0.075),                 ("02", 0.2), ("03", 0.3), ("04", 0.4)]

# NEW tower iso suffixes — flat list; each element is the string between
# "cluster_iso_" and "_{NODE}" in the branch name. Every entry is saved raw to npz
# as `iso_x_{suffix}` for all three samples; the ROOT-file maker sums them as needed.
TOWER_ISO_EXTRA = [
    # (A) per-calo split at R=0.3, R=0.4 (no threshold)  → sum reproduces the plain branch
    "03_emcal", "03_hcalin", "03_hcalout",
    "04_emcal", "04_hcalin", "04_hcalout",
    # (B) threshold scan at R=0.3 (60 / 70 / 120 MeV, per calo)
    "03_60_emcal",  "03_60_hcalin",  "03_60_hcalout",
    "03_70_emcal",  "03_70_hcalin",  "03_70_hcalout",
    "03_120_emcal", "03_120_hcalin", "03_120_hcalout",
    # (C) EMCal-only 70 MeV at scanning R (03_70_emcal already in (B))
    "005_70_emcal", "0075_70_emcal", "01_70_emcal", "02_70_emcal",
    # (D) SUB1-subtracted at R=0.3, R=0.4 (per calo)
    "03_sub1_emcal", "03_sub1_hcalin", "03_sub1_hcalout",
    "04_sub1_emcal", "04_sub1_hcalin", "04_sub1_hcalout",
    # (E) ownership-excluded iso — branch name is cluster_iso_excl_{R}_{NODE}
    "excl_005", "excl_0075", "excl_01", "excl_02", "excl_03", "excl_04",
]

NPB_CUT = 0.5

# Cross-sections (pb) from CrossSectionWeights.h. Per-event weight = xsec / n_events.
# Relative normalization is sufficient for shape comparison.
XSEC = {
    "photon5":  146359.3,
    "photon10": 6944.675,
    "photon20": 130.4461,
    "jet5":     1.3878e+08,
    "jet8":     1.15e+07,
    "jet12":    1.4903e+06,
    "jet20":    6.2623e+04,
    "jet30":    2.5298e+03,
    "jet40":    1.3553e+02,
}
SIGNAL_SAMPLES = ["photon5", "photon10", "photon20"]
JET_SAMPLES    = ["jet5", "jet8", "jet12", "jet20", "jet30", "jet40"]
N_EVENTS_NOM   = 1e7  # approx for all samples

MC_PATH = lambda s: f"/sphenix/user/shuhangli/ppg12/FunWithxgboost/{s}/bdt_split.root"
DATA_GLOB = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/condorout/part_*_with_bdt_split.root"


def _flat_iso(chunk, suffix, nclu, present):
    """Flatten a per-event iso branch to a per-cluster 1D array; NaN if absent."""
    key = f"cluster_iso_{suffix}_{NODE}"
    if key in present and key in chunk:
        return np.concatenate([chunk[key][i][:nclu[i]] for i in range(len(nclu))])
    n = int(np.sum(nclu))
    return np.full(n, np.nan, dtype=np.float32)


def _present_branches(path_glob: str):
    """Open the first ROOT matching path_glob and return the set of branch names."""
    import glob as _glob
    paths = sorted(_glob.glob(path_glob)) if any(c in path_glob for c in "*?[") else [path_glob]
    with uproot.open(f"{paths[0]}:slimtree") as t:
        return set(t.keys())


def read_inclusive(path_glob: str, weight: float, max_events=None):
    """Iso (topo + tower + new tower variants) + npb + cluster_Et, no truth match."""
    iso      = {rs: [] for rs, _ in RADII}
    iso_tw   = {rs: [] for rs, _ in RADII_TOWER}
    iso_xtra = {sfx: [] for sfx in TOWER_ISO_EXTRA}
    clu_Ets, npbs, ws = [], [], []
    present = _present_branches(path_glob)
    branches = [
        f"ncluster_{NODE}", f"cluster_Et_{NODE}", f"cluster_npb_score_{NODE}",
    ] + [f"cluster_iso_topo_{rs}_{NODE}" for rs, _ in RADII] \
      + [f"cluster_iso_{rs}_{NODE}"       for rs, _ in RADII_TOWER] \
      + [f"cluster_iso_{sfx}_{NODE}"      for sfx in TOWER_ISO_EXTRA
         if f"cluster_iso_{sfx}_{NODE}" in present]
    missing_xtra = [sfx for sfx in TOWER_ISO_EXTRA
                    if f"cluster_iso_{sfx}_{NODE}" not in present]
    if missing_xtra:
        print(f"  WARNING: {path_glob} missing {len(missing_xtra)} extra iso branches "
              f"(will fill NaN): {missing_xtra}")
    total = 0
    for chunk in uproot.iterate(f"{path_glob}:slimtree", branches,
                                 step_size=50_000, library="np"):
        nclu = chunk[f"ncluster_{NODE}"]
        clu_Et  = np.concatenate([chunk[f"cluster_Et_{NODE}"][i][:nclu[i]]        for i in range(len(nclu))])
        clu_npb = np.concatenate([chunk[f"cluster_npb_score_{NODE}"][i][:nclu[i]] for i in range(len(nclu))])
        clu_Ets.append(clu_Et); npbs.append(clu_npb)
        ws.append(np.full(len(clu_Et), weight, dtype=np.float32))
        for rs, _ in RADII:
            iso[rs].append(np.concatenate([chunk[f"cluster_iso_topo_{rs}_{NODE}"][i][:nclu[i]]
                                           for i in range(len(nclu))]))
        for rs, _ in RADII_TOWER:
            iso_tw[rs].append(np.concatenate([chunk[f"cluster_iso_{rs}_{NODE}"][i][:nclu[i]]
                                              for i in range(len(nclu))]))
        for sfx in TOWER_ISO_EXTRA:
            iso_xtra[sfx].append(_flat_iso(chunk, sfx, nclu, present))
        total += len(clu_Et)
        if max_events and total >= max_events:
            break
    return ({rs: np.concatenate(iso[rs])     for rs, _ in RADII},
            {rs: np.concatenate(iso_tw[rs])  for rs, _ in RADII_TOWER},
            {sfx: np.concatenate(iso_xtra[sfx]) for sfx in TOWER_ISO_EXTRA},
            np.concatenate(clu_Ets),
            np.concatenate(npbs),
            np.concatenate(ws))


def process_mc_truth(path: str, weight: float, max_events=None):
    """Truth-matched direct-iso photon pass — topo + tower + new tower variants."""
    iso    = {rs: [] for rs, _ in RADII}
    iso_tw = {rs: [] for rs, _ in RADII_TOWER}
    iso_xtra = {sfx: [] for sfx in TOWER_ISO_EXTRA}
    clu_Ets, part_Pts, npbs, ws = [], [], [], []
    present = _present_branches(path)
    branches = [
        f"ncluster_{NODE}", f"cluster_Et_{NODE}", f"cluster_Eta_{NODE}",
        f"cluster_truthtrkID_{NODE}", f"cluster_npb_score_{NODE}",
    ] + [f"cluster_iso_topo_{rs}_{NODE}" for rs, _ in RADII] \
      + [f"cluster_iso_{rs}_{NODE}"       for rs, _ in RADII_TOWER] \
      + [f"cluster_iso_{sfx}_{NODE}"      for sfx in TOWER_ISO_EXTRA
         if f"cluster_iso_{sfx}_{NODE}" in present] + [
        "nparticles", "particle_trkid", "particle_Pt",
        "particle_pid", "particle_photonclass", "particle_truth_iso_03",
    ]
    missing_xtra = [sfx for sfx in TOWER_ISO_EXTRA
                    if f"cluster_iso_{sfx}_{NODE}" not in present]
    if missing_xtra:
        print(f"  WARNING: {path} missing {len(missing_xtra)} extra iso branches "
              f"(will fill NaN): {missing_xtra}")
    total = 0
    with uproot.open(f"{path}:slimtree") as t:
        print(f"  truth-match on {path}: {t.num_entries:,} events")
        for chunk in t.iterate(branches, step_size=50_000, library="np"):
            nclu = chunk[f"ncluster_{NODE}"]
            clu_Et  = np.concatenate([chunk[f"cluster_Et_{NODE}"][i][:nclu[i]]        for i in range(len(nclu))])
            clu_tID = np.concatenate([chunk[f"cluster_truthtrkID_{NODE}"][i][:nclu[i]] for i in range(len(nclu))])
            clu_npb = np.concatenate([chunk[f"cluster_npb_score_{NODE}"][i][:nclu[i]] for i in range(len(nclu))])
            npart = chunk["nparticles"]
            p_trkid = chunk["particle_trkid"]; p_pt = chunk["particle_Pt"]
            p_pid = chunk["particle_pid"]; p_pcl = chunk["particle_photonclass"]
            p_iso3 = chunk["particle_truth_iso_03"]

            matched_pt = np.full(len(clu_Et), -1.0, dtype=np.float32)
            mask_direct = np.zeros(len(clu_Et), dtype=bool)
            cursor = 0
            for ev in range(len(nclu)):
                nc = nclu[ev]
                if nc == 0: continue
                this_trk = p_trkid[ev][:npart[ev]]
                for ci in range(nc):
                    tid = clu_tID[cursor + ci]
                    if tid < 0: continue
                    idx = np.where(this_trk == tid)[0]
                    if len(idx) == 0: continue
                    p = idx[0]
                    if p_pid[ev][p] != 22 or p_pcl[ev][p] != 1 or p_iso3[ev][p] >= 4.0:
                        continue
                    matched_pt[cursor + ci] = p_pt[ev][p]
                    mask_direct[cursor + ci] = True
                cursor += nc

            m = mask_direct
            clu_Ets.append(clu_Et[m]); part_Pts.append(matched_pt[m])
            npbs.append(clu_npb[m])
            ws.append(np.full(m.sum(), weight, dtype=np.float32))
            for rs, _ in RADII:
                iso[rs].append(_flat_iso_topo(chunk, rs, nclu)[m])
            for rs, _ in RADII_TOWER:
                iso_tw[rs].append(np.concatenate([chunk[f"cluster_iso_{rs}_{NODE}"][i][:nclu[i]]
                                                  for i in range(len(nclu))])[m])
            for sfx in TOWER_ISO_EXTRA:
                iso_xtra[sfx].append(_flat_iso(chunk, sfx, nclu, present)[m])
            total += len(clu_Et)
            if max_events and total >= max_events:
                break
    return ({rs: np.concatenate(iso[rs])    for rs, _ in RADII},
            {rs: np.concatenate(iso_tw[rs]) for rs, _ in RADII_TOWER},
            {sfx: np.concatenate(iso_xtra[sfx]) for sfx in TOWER_ISO_EXTRA},
            np.concatenate(clu_Ets),
            np.concatenate(part_Pts),
            np.concatenate(npbs),
            np.concatenate(ws))


def _flat_iso_topo(chunk, rs, nclu):
    key = f"cluster_iso_topo_{rs}_{NODE}"
    return np.concatenate([chunk[key][i][:nclu[i]] for i in range(len(nclu))])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-mc",   type=int, default=150_000,
                    help="max clusters per inclusive-MC sample (default 150k for speed; "
                         "~13x faster than the original 2M per-sample cap)")
    ap.add_argument("--max-data", type=int, default=300_000,
                    help="max data clusters (default 300k, ~10x faster than 3M)")
    ap.add_argument("--max-signal-mc", type=int, default=700_000,
                    help="max clusters scanned per truth-signal-MC sample (default 700k). "
                         "Only direct-iso photons kept, so final truth-matched count is much smaller.")
    ap.add_argument("--out", default="results/iso_topo_containment.npz")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # ---------- MC: truth-matched signal (containment) ----------
    t_iso = {rs: [] for rs, _ in RADII}
    t_iso_tw = {rs: [] for rs, _ in RADII_TOWER}
    t_iso_xtra = {sfx: [] for sfx in TOWER_ISO_EXTRA}
    t_clu_Et, t_part_Pt, t_npb, t_w = [], [], [], []
    for s in SIGNAL_SAMPLES:
        w = XSEC[s] / N_EVENTS_NOM
        print(f"MC truth-match {s}  (weight {w:.3e})")
        iso, iso_tw, iso_xtra, cEt, pPt, np_, ww = process_mc_truth(
            MC_PATH(s), w, max_events=args.max_signal_mc)
        for rs, _ in RADII:       t_iso[rs].append(iso[rs])
        for rs, _ in RADII_TOWER: t_iso_tw[rs].append(iso_tw[rs])
        for sfx in TOWER_ISO_EXTRA: t_iso_xtra[sfx].append(iso_xtra[sfx])
        t_clu_Et.append(cEt); t_part_Pt.append(pPt); t_npb.append(np_); t_w.append(ww)
    T       = {rs:  np.concatenate(t_iso[rs])     for rs, _ in RADII}
    T_tw    = {rs:  np.concatenate(t_iso_tw[rs])  for rs, _ in RADII_TOWER}
    T_xtra  = {sfx: np.concatenate(t_iso_xtra[sfx]) for sfx in TOWER_ISO_EXTRA}
    t_clu_Et = np.concatenate(t_clu_Et); t_part_Pt = np.concatenate(t_part_Pt)
    t_npb    = np.concatenate(t_npb);    t_w       = np.concatenate(t_w)
    print(f"  truth-matched direct-iso photons: {len(t_clu_Et):,}")

    # ---------- MC: inclusive (signal + jets, cross-section weighted) ----------
    mi_iso    = {rs: [] for rs, _ in RADII}
    mi_iso_tw = {rs: [] for rs, _ in RADII_TOWER}
    mi_iso_xtra = {sfx: [] for sfx in TOWER_ISO_EXTRA}
    mi_clu_Et, mi_npb, mi_w = [], [], []
    for s in SIGNAL_SAMPLES + JET_SAMPLES:
        w = XSEC[s] / N_EVENTS_NOM
        print(f"MC inclusive {s} (weight {w:.3e})")
        iso, iso_tw, iso_xtra, cEt, np_, ww = read_inclusive(
            MC_PATH(s), w, max_events=args.max_mc)
        for rs, _ in RADII:       mi_iso[rs].append(iso[rs])
        for rs, _ in RADII_TOWER: mi_iso_tw[rs].append(iso_tw[rs])
        for sfx in TOWER_ISO_EXTRA: mi_iso_xtra[sfx].append(iso_xtra[sfx])
        mi_clu_Et.append(cEt); mi_npb.append(np_); mi_w.append(ww)
    M      = {rs: np.concatenate(mi_iso[rs])    for rs, _ in RADII}
    M_tw   = {rs: np.concatenate(mi_iso_tw[rs]) for rs, _ in RADII_TOWER}
    M_xtra = {sfx: np.concatenate(mi_iso_xtra[sfx]) for sfx in TOWER_ISO_EXTRA}
    mi_clu_Et = np.concatenate(mi_clu_Et); mi_npb = np.concatenate(mi_npb); mi_w = np.concatenate(mi_w)
    print(f"  inclusive MC clusters: {len(mi_clu_Et):,}")

    # ---------- DATA ----------
    print("DATA (77 parts)")
    D_iso, D_iso_tw, D_iso_xtra, d_clu_Et, d_npb, _ = read_inclusive(
        DATA_GLOB, 1.0, max_events=args.max_data)
    print(f"  data clusters: {len(d_clu_Et):,}")

    # Use compressed npz to stay under GPFS user quota on higher-stats runs
    np.savez_compressed(out_path,
        **{f"tr_iso_{rs}":    T[rs]    for rs, _ in RADII},       tr_clu_Et=t_clu_Et,  tr_part_Pt=t_part_Pt, tr_npb=t_npb, tr_w=t_w,
        **{f"tr_iso_tw_{rs}": T_tw[rs] for rs, _ in RADII_TOWER},
        **{f"tr_iso_x_{sfx}": T_xtra[sfx] for sfx in TOWER_ISO_EXTRA},
        **{f"mc_iso_{rs}":    M[rs]    for rs, _ in RADII},       mc_clu_Et=mi_clu_Et, mc_npb=mi_npb, mc_w=mi_w,
        **{f"mc_iso_tw_{rs}": M_tw[rs] for rs, _ in RADII_TOWER},
        **{f"mc_iso_x_{sfx}": M_xtra[sfx] for sfx in TOWER_ISO_EXTRA},
        **{f"d_iso_{rs}":     D_iso[rs]    for rs, _ in RADII},   d_clu_Et=d_clu_Et, d_npb=d_npb,
        **{f"d_iso_tw_{rs}":  D_iso_tw[rs] for rs, _ in RADII_TOWER},
        **{f"d_iso_x_{sfx}":  D_iso_xtra[sfx] for sfx in TOWER_ISO_EXTRA},
    )
    print(f"Saved {out_path}")

    # ---------- Summary tables ----------
    def stats(vals, w, mask):
        vv = vals[mask]; ww = w[mask] if w is not None else None
        if len(vv) == 0: return None
        if ww is None:
            return vv.mean(), vv.std(), np.median(vv), len(vv)
        mu = np.average(vv, weights=ww)
        rms = np.sqrt(np.average((vv - mu) ** 2, weights=ww))
        order = np.argsort(vv); cw = np.cumsum(ww[order])
        med = vv[order][np.searchsorted(cw, cw[-1]/2.0)]
        return mu, rms, med, len(vv)

    for name, iso_map, clu_Et, npb, w in [
        ("MC truth-matched photons", T, t_clu_Et, t_npb, t_w),
        ("MC inclusive",             M, mi_clu_Et, mi_npb, mi_w),
        ("DATA",                     D_iso, d_clu_Et, d_npb, None),
    ]:
        for state_lbl, m in [("pre-common", np.ones(len(clu_Et), dtype=bool)),
                             (f"post-common (npb>{NPB_CUT})", npb > NPB_CUT)]:
            n_sel = int(m.sum())
            print(f"\n[Study A] {name}  [{state_lbl}]  N={n_sel:,}")
            print(f"  {'R':<8}{'mean':>10}{'RMS':>10}{'median':>10}")
            for rs, r in RADII:
                s = stats(iso_map[rs], w, m)
                if s is None: print(f"  {r:<8.3f}  (empty)"); continue
                mu, rms, med, n = s
                print(f"  {r:<8.3f}{mu:>10.4f}{rms:>10.4f}{med:>10.4f}")

    # Containment (MC truth-matched only)
    for state_lbl, m in [("pre-common", np.ones(len(t_clu_Et), dtype=bool)),
                         (f"post-common (npb>{NPB_CUT})", t_npb > NPB_CUT)]:
        print(f"\n[Study B] Containment f_contain(R) = (iso_R + cluster_ET) / truth_Pt  [{state_lbl}]")
        print(f"  {'R':<8}{'<f_contain>':>15}{'<leakage>':>15}")
        for rs, r in RADII:
            f = (T[rs][m] + t_clu_Et[m]) / t_part_Pt[m]
            if len(f) == 0: print(f"  {r:<8.3f} (empty)"); continue
            ww = t_w[m]
            mu = np.average(f, weights=ww)
            print(f"  {r:<8.3f}{mu:>15.4f}{1 - mu:>15.4f}")


if __name__ == "__main__":
    main()
