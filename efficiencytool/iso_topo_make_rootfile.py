#!/usr/bin/env python3
"""
Read iso_topo_containment.npz and write a ROOT file with:
 - 1D distributions h_iso_{sample}_{R}_{cut}  (topo)
 - 1D distributions h_iso_tw_{sample}_{R}_{cut}  (plain tower)
 - New-variant distributions h_iso_tw_{variant}_{sample}_{...}_{cut}  (tower)
 - 1D containment f_contain per R per truth-pT bin (TH1D)
 - Summary numbers as a TTree

Output: efficiencytool/results/iso_topo_containment.root

Sample  in {tr, mc, d}    (truth-matched sig / inclusive MC / data)
Cut     in {pre, post}    (pre-common / post-common npb>0.5)
R       in {005,0075,01,02,03,04} for topo, {005,0075,02,03,04} for plain tower

New tower variants (added Apr 2026):
 - thr_nothr_{03,04}: per-calo split summed, no threshold  (= plain total, sanity)
 - thr_60_03, thr_70_03, thr_120_03: sum of 3 calos with per-tower threshold 60/70/120 MeV
 - emcal70_{005,0075,01,02,03}: EMCal-only, 70 MeV threshold, scanning R
 - sub1_{03,04}: sum of 3 calos, SUB1 (background) subtracted
 - excl_{005,0075,01,02,03,04}: ownership-excluded iso (cluster_iso_excl_R)

Truth-pT bins for (b): [8,14], [14,20], [20,28], [28,40]  (only makes sense for tr)
"""
import numpy as np
import uproot
from pathlib import Path

NPZ = Path("/sphenix/user/shuhangli/ppg12/efficiencytool/results/iso_topo_containment.npz")
OUT = Path("/sphenix/user/shuhangli/ppg12/efficiencytool/results/iso_topo_containment.root")
NPB = 0.5
RADII       = [("005", 0.05), ("0075", 0.075), ("01", 0.1), ("02", 0.2), ("03", 0.3), ("04", 0.4)]
RADII_TOWER = [("005", 0.05), ("0075", 0.075),               ("02", 0.2), ("03", 0.3), ("04", 0.4)]
PT_BIN_EDGES = [8, 14, 20, 28, 40]

ISO_BINS = 100
ISO_LO, ISO_HI = -20, 30
CLU_ET_BINS = [(8, 14), (14, 20), (20, 28), (28, 40)]

# New tower variants — each builds a single per-cluster value from one or more npz keys.
# For sums (threshold scan / SUB1), lookup returns a list of suffixes whose npz arrays are added.
# For direct branches (emcal70, excl), list has one suffix.
TOWER_VARIANTS = {
    # (A) per-calo-split sum with no threshold — sanity check vs plain
    "thr_nothr_03": ["03_emcal", "03_hcalin", "03_hcalout"],
    "thr_nothr_04": ["04_emcal", "04_hcalin", "04_hcalout"],
    # (B) threshold scan at R=0.3 (EMCal+HCalIn+HCalOut)
    "thr_60_03":    ["03_60_emcal",  "03_60_hcalin",  "03_60_hcalout"],
    "thr_70_03":    ["03_70_emcal",  "03_70_hcalin",  "03_70_hcalout"],
    "thr_120_03":   ["03_120_emcal", "03_120_hcalin", "03_120_hcalout"],
    # (C) EMCal-only, 70 MeV threshold, scanning R (single branch each)
    "emcal70_005":  ["005_70_emcal"],
    "emcal70_0075": ["0075_70_emcal"],
    "emcal70_01":   ["01_70_emcal"],
    "emcal70_02":   ["02_70_emcal"],
    "emcal70_03":   ["03_70_emcal"],
    # (D) SUB1 subtracted (EMCal+HCalIn+HCalOut)
    "sub1_03":      ["03_sub1_emcal", "03_sub1_hcalin", "03_sub1_hcalout"],
    "sub1_04":      ["04_sub1_emcal", "04_sub1_hcalin", "04_sub1_hcalout"],
    # (E) ownership-excluded (single branch each)
    "excl_005":  ["excl_005"],
    "excl_0075": ["excl_0075"],
    "excl_01":   ["excl_01"],
    "excl_02":   ["excl_02"],
    "excl_03":   ["excl_03"],
    "excl_04":   ["excl_04"],
}


def build_variant(d, sample, suffixes):
    """Return summed per-cluster iso array for this variant; NaN if any component is NaN."""
    arrs = [d[f"{sample}_iso_x_{s}"] for s in suffixes]
    out = np.zeros_like(arrs[0], dtype=np.float64)
    any_nan = np.zeros_like(arrs[0], dtype=bool)
    for a in arrs:
        out += np.nan_to_num(a, nan=0.0)
        any_nan |= np.isnan(a)
    out[any_nan] = np.nan
    return out


def _hist_from_vals(vals, w, bins, range_):
    m = ~np.isnan(vals)
    vv = vals[m]
    ww = w[m] if w is not None else None
    h, edges = np.histogram(vv, bins=bins, range=range_, weights=ww)
    integ = h.sum() * (edges[1] - edges[0])
    if integ > 0:
        h = h / integ
    return h, edges


def main():
    d = np.load(NPZ)
    print(f"loaded {NPZ}  ({len(d.files)} keys)")

    with uproot.recreate(OUT) as f:
        # Study A — topo iso_topo_R distributions (inclusive + per cluster_Et bin)
        for sample, iso_pat, ck, npbk, wk in [
            ("tr", "tr_iso_{}", "tr_clu_Et", "tr_npb", "tr_w"),
            ("mc", "mc_iso_{}", "mc_clu_Et", "mc_npb", "mc_w"),
            ("d",  "d_iso_{}",  "d_clu_Et",  "d_npb",  None),
        ]:
            npb = d[npbk]
            clu_Et = d[ck]
            w = d[wk] if wk else None
            for cut_name, base_mask in [("pre",  np.ones(len(npb), bool)),
                                        ("post", npb > NPB)]:
                bin_variants = [("incl", base_mask)]
                for lo, hi in CLU_ET_BINS:
                    bin_variants.append((f"clEt_{lo}_{hi}",
                                          base_mask & (clu_Et >= lo) & (clu_Et < hi)))
                for bin_label, m in bin_variants:
                    ww = w[m] if w is not None else None
                    for rs, r in RADII:
                        vals = d[iso_pat.format(rs)][m]
                        if len(vals) == 0: continue
                        h, edges = _hist_from_vals(vals, ww, ISO_BINS, (ISO_LO, ISO_HI))
                        if bin_label == "incl":
                            f[f"h_iso_{sample}_{rs}_{cut_name}"] = (h, edges)
                        else:
                            f[f"h_iso_{sample}_{rs}_{cut_name}_{bin_label}"] = (h, edges)
        print("  Study A (topo) distributions written")

        # Study A' — plain tower iso_R distributions
        for sample, iso_pat, ck, npbk, wk in [
            ("tr", "tr_iso_tw_{}", "tr_clu_Et", "tr_npb", "tr_w"),
            ("mc", "mc_iso_tw_{}", "mc_clu_Et", "mc_npb", "mc_w"),
            ("d",  "d_iso_tw_{}",  "d_clu_Et",  "d_npb",  None),
        ]:
            if iso_pat.format(RADII_TOWER[0][0]) not in d.files:
                print(f"    tower branches missing for sample {sample}, skipping")
                continue
            npb = d[npbk]
            clu_Et = d[ck]
            w = d[wk] if wk else None
            for cut_name, base_mask in [("pre",  np.ones(len(npb), bool)),
                                        ("post", npb > NPB)]:
                bin_variants = [("incl", base_mask)]
                for lo, hi in CLU_ET_BINS:
                    bin_variants.append((f"clEt_{lo}_{hi}",
                                          base_mask & (clu_Et >= lo) & (clu_Et < hi)))
                for bin_label, m in bin_variants:
                    ww = w[m] if w is not None else None
                    for rs, r in RADII_TOWER:
                        vals = d[iso_pat.format(rs)][m]
                        if len(vals) == 0: continue
                        h, edges = _hist_from_vals(vals, ww, ISO_BINS, (ISO_LO, ISO_HI))
                        if bin_label == "incl":
                            f[f"h_iso_tw_{sample}_{rs}_{cut_name}"] = (h, edges)
                        else:
                            f[f"h_iso_tw_{sample}_{rs}_{cut_name}_{bin_label}"] = (h, edges)
        print("  Study A' (plain tower) distributions written")

        # Study A'' — NEW tower iso variants
        # Check that the npz has these keys
        sample_map = {
            "tr": ("tr_clu_Et", "tr_npb", "tr_w"),
            "mc": ("mc_clu_Et", "mc_npb", "mc_w"),
            "d":  ("d_clu_Et",  "d_npb",  None),
        }
        nvar_written = 0
        for vname, suffixes in TOWER_VARIANTS.items():
            # sanity: need every suffix available in npz
            needed_keys = [f"tr_iso_x_{s}" for s in suffixes] + \
                          [f"mc_iso_x_{s}" for s in suffixes] + \
                          [f"d_iso_x_{s}"  for s in suffixes]
            if not all(k in d.files for k in needed_keys):
                print(f"    variant {vname} missing extra-iso keys in npz, skipping")
                continue
            for sample, (ck, npbk, wk) in sample_map.items():
                vals_all = build_variant(d, sample, suffixes)
                npb = d[npbk]
                clu_Et = d[ck]
                w = d[wk] if wk else None
                for cut_name, base_mask in [("pre",  np.ones(len(npb), bool)),
                                            ("post", npb > NPB)]:
                    bin_variants = [("incl", base_mask)]
                    for lo, hi in CLU_ET_BINS:
                        bin_variants.append((f"clEt_{lo}_{hi}",
                                              base_mask & (clu_Et >= lo) & (clu_Et < hi)))
                    for bin_label, m in bin_variants:
                        ww = w[m] if w is not None else None
                        vals = vals_all[m]
                        if len(vals) == 0: continue
                        # Skip if all-NaN (i.e., variant not available for any sample-component)
                        if np.all(np.isnan(vals)): continue
                        h, edges = _hist_from_vals(vals, ww, ISO_BINS, (ISO_LO, ISO_HI))
                        if bin_label == "incl":
                            f[f"h_iso_tw_{vname}_{sample}_{cut_name}"] = (h, edges)
                        else:
                            f[f"h_iso_tw_{vname}_{sample}_{cut_name}_{bin_label}"] = (h, edges)
            nvar_written += 1
        print(f"  Study A'' ({nvar_written} new tower variants) distributions written")

        # Study B — photon containment per truth-pT bin, topo (truth-matched only)
        tr_pt = d["tr_part_Pt"]
        tr_clu = d["tr_clu_Et"]
        tr_w  = d["tr_w"]
        tr_npb = d["tr_npb"]

        for cut_name, cut_mask in [("pre",  np.ones(len(tr_pt), bool)),
                                   ("post", tr_npb > NPB)]:
            for rs, r in RADII:
                vals = (d[f"tr_iso_{rs}"][cut_mask] + tr_clu[cut_mask]) / tr_pt[cut_mask]
                h, edges = np.histogram(vals, bins=100, range=(0, 1.5),
                                        weights=tr_w[cut_mask])
                integ = h.sum() * (edges[1] - edges[0])
                if integ > 0: h = h / integ
                f[f"h_fcont_{cut_name}_inclusive_{rs}"] = (h, edges)
            for ipt in range(len(PT_BIN_EDGES) - 1):
                lo, hi = PT_BIN_EDGES[ipt], PT_BIN_EDGES[ipt+1]
                ptmask = cut_mask & (tr_pt >= lo) & (tr_pt < hi)
                for rs, r in RADII:
                    vals = (d[f"tr_iso_{rs}"][ptmask] + tr_clu[ptmask]) / tr_pt[ptmask]
                    if len(vals) == 0: continue
                    h, edges = np.histogram(vals, bins=100, range=(0, 1.5),
                                            weights=tr_w[ptmask])
                    integ = h.sum() * (edges[1] - edges[0])
                    if integ > 0: h = h / integ
                    f[f"h_fcont_{cut_name}_pt_{lo}_{hi}_{rs}"] = (h, edges)
        print("  Study B containment histograms written (pT-binned)")

        # Summary TTree for containment
        arrs = {"cut": [], "pt_lo": [], "pt_hi": [], "r": [], "mean": [], "rms": []}
        for cut_name, cut_mask in [("pre",  np.ones(len(tr_pt), bool)),
                                   ("post", tr_npb > NPB)]:
            bin_list = [("inclusive", 8, 40)] + [(f"pt_{PT_BIN_EDGES[i]}_{PT_BIN_EDGES[i+1]}",
                                                   PT_BIN_EDGES[i], PT_BIN_EDGES[i+1])
                                                  for i in range(len(PT_BIN_EDGES)-1)]
            for lbl, lo, hi in bin_list:
                if lbl == "inclusive":
                    mask = cut_mask
                else:
                    mask = cut_mask & (tr_pt >= lo) & (tr_pt < hi)
                if mask.sum() == 0: continue
                for rs, r in RADII:
                    vals = (d[f"tr_iso_{rs}"][mask] + tr_clu[mask]) / tr_pt[mask]
                    ww = tr_w[mask]
                    mu = np.average(vals, weights=ww)
                    rms = np.sqrt(np.average((vals - mu)**2, weights=ww))
                    arrs["cut"].append(cut_name)
                    arrs["pt_lo"].append(lo)
                    arrs["pt_hi"].append(hi)
                    arrs["r"].append(r)
                    arrs["mean"].append(mu)
                    arrs["rms"].append(rms)
        f["fcont_summary"] = {k: np.array(arrs[k]) for k in arrs}
        print(f"  fcont_summary tree: {len(arrs['r'])} rows")

    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
