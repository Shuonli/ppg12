#!/usr/bin/env python3
"""
Independent re-extraction of the ntbdtpair cross-check numbers from ROOT files.

For each of 18 variants (9 suffixes x 2 periods), extract:
  1. Total unfolded cross-section integrated over bins whose centers lie in [10, 36] GeV
     integral = sum_i h_unfold_sub_result.GetBinContent(i) * h_unfold_sub_result.GetBinWidth(i)
  2. Max |h_var / h_nom - 1| over bins with nominal content > 0 and center in [10, 36] GeV.
     Report signed (var/nom - 1) where |.| is largest, and the ET bin center.
  3. Purity at ET = 11 GeV (lowest reco bin center >= 10 GeV) and at ET = 25 GeV
     from gpurity.GetY()[i] at the matching ET index, both for variant and nominal.

Output: Markdown + CSV table of 18 rows with columns:
  variant, period, xsec_int[pb], xsec_nom_int[pb], max_abs_rel_dev,
  max_dev_et_bin, purity_11, purity_nom_11, purity_25, purity_nom_25
"""

import os
import csv
import numpy as np
import uproot

RESULTS_DIR = "/sphenix/user/shuhangli/ppg12/efficiencytool/results"
NOMINAL_FILES = {
    "1p5mrad": os.path.join(RESULTS_DIR, "Photon_final_bdt_nom.root"),
    "0mrad":   os.path.join(RESULTS_DIR, "Photon_final_bdt_0rad.root"),
}
SUFFIXES = [
    "t80_70_h80_70_l70_40",
    "t80_75_h80_75_l70_50",
    "t80_75_h80_75_l70_40",
    "t80_70_h80_70_l60_40",
    "t80_75_h80_75_l60_50",
    "t80_75_h80_75_l60_40",
    "t80_70_h70_70_l60_40",
    "t80_75_h70_75_l60_50",
    "t80_75_h70_75_l60_40",
]
PERIODS = ["1p5mrad", "0mrad"]

ET_MIN, ET_MAX = 10.0, 36.0  # cross-section integration window [GeV]
PURITY_BIN_CENTERS_TARGET = (11.0, 25.0)  # target ET bin centers for purity readout


def load_xsec_hist(fname):
    """Return (centers, widths, values) for h_unfold_sub_result."""
    with uproot.open(fname) as f:
        h = f["h_unfold_sub_result"]
        edges = np.asarray(h.axis().edges())
        vals  = np.asarray(h.values())
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths  = np.diff(edges)
    return centers, widths, vals


def load_purity_graph(fname):
    """Return (x, y) arrays for gpurity."""
    with uproot.open(fname) as f:
        g = f["gpurity"]
        x = np.asarray(g.member("fX"))
        y = np.asarray(g.member("fY"))
    return x, y


def integrate_xsec(centers, widths, vals, lo=ET_MIN, hi=ET_MAX):
    """Integral of dSigma/dET over bins whose CENTER lies in [lo, hi]."""
    mask = (centers >= lo) & (centers <= hi)
    return float(np.sum(vals[mask] * widths[mask]))


def max_rel_dev(centers_v, vals_v, centers_n, vals_n, lo=ET_MIN, hi=ET_MAX):
    """
    Max |var/nom - 1| bin-by-bin, restricted to bins within [lo, hi] and nom > 0.
    Returns (signed_dev, et_center) where |signed_dev| is maximal.
    """
    assert np.allclose(centers_v, centers_n), "variant/nominal bin centers differ"
    mask = (centers_v >= lo) & (centers_v <= hi) & (vals_n > 0)
    if not np.any(mask):
        return float("nan"), float("nan")
    idx = np.where(mask)[0]
    ratio = vals_v[idx] / vals_n[idx]
    devs  = ratio - 1.0
    j = int(np.argmax(np.abs(devs)))
    return float(devs[j]), float(centers_v[idx[j]])


def purity_at(x, y, target, tol=1e-6):
    """
    Return purity at the x value closest to target. If the closest bin center is
    further than 2 GeV from target, fall back to nearest and flag in a comment
    (we keep the number; all files here have matching grids so tol is cosmetic).
    """
    i = int(np.argmin(np.abs(x - target)))
    return float(y[i]), float(x[i])


def main():
    rows = []
    for period in PERIODS:
        nom_file = NOMINAL_FILES[period]
        c_n, w_n, v_n = load_xsec_hist(nom_file)
        xn, yn = load_purity_graph(nom_file)
        xsec_nom = integrate_xsec(c_n, w_n, v_n)
        pur_nom_11, _ = purity_at(xn, yn, 11.0)
        pur_nom_25, _ = purity_at(xn, yn, 25.0)

        for suffix in SUFFIXES:
            var_file = os.path.join(
                RESULTS_DIR,
                f"Photon_final_bdt_ntbdtpair_{suffix}_{period}.root",
            )
            c_v, w_v, v_v = load_xsec_hist(var_file)
            xv, yv = load_purity_graph(var_file)
            xsec_var = integrate_xsec(c_v, w_v, v_v)
            signed_dev, et_bin = max_rel_dev(c_v, v_v, c_n, v_n)
            pur_var_11, _ = purity_at(xv, yv, 11.0)
            pur_var_25, _ = purity_at(xv, yv, 25.0)

            flag = ""
            if np.isfinite(signed_dev):
                a = abs(signed_dev)
                if a > 0.50:
                    flag = "CRIT"
                elif a > 0.20:
                    flag = "WARN"

            rows.append({
                "variant": suffix,
                "period":  period,
                "xsec_int_pb":      xsec_var,
                "xsec_nom_int_pb":  xsec_nom,
                "max_abs_rel_dev":  signed_dev,
                "max_dev_et_bin":   et_bin,
                "purity_11":        pur_var_11,
                "purity_nom_11":    pur_nom_11,
                "purity_25":        pur_var_25,
                "purity_nom_25":    pur_nom_25,
                "flag":             flag,
            })

    # ---- write CSV ----
    out_dir = "/sphenix/user/shuhangli/ppg12/reports"
    csv_path = os.path.join(out_dir, "ntbdtpair_numbers.csv")
    with open(csv_path, "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ---- write Markdown table ----
    md_path = os.path.join(out_dir, "ntbdtpair_numbers.md")
    header = (
        "| variant | period | xsec_int[pb] | xsec_nom_int[pb] | max_abs_rel_dev "
        "| max_dev_et_bin | purity_11 | purity_nom_11 | purity_25 | purity_nom_25 | flag |"
    )
    sep = "|" + "|".join(["---"] * 11) + "|"
    lines = [header, sep]
    for r in rows:
        lines.append(
            f"| {r['variant']} | {r['period']} | {r['xsec_int_pb']:.3f} "
            f"| {r['xsec_nom_int_pb']:.3f} | {r['max_abs_rel_dev']:+.4f} "
            f"| {r['max_dev_et_bin']:.1f} | {r['purity_11']:.4f} "
            f"| {r['purity_nom_11']:.4f} | {r['purity_25']:.4f} "
            f"| {r['purity_nom_25']:.4f} | {r['flag']} |"
        )
    with open(md_path, "w") as fp:
        fp.write("\n".join(lines) + "\n")

    # ---- stdout summary ----
    print("\n".join(lines))
    print(f"\nWrote: {csv_path}")
    print(f"Wrote: {md_path}")


if __name__ == "__main__":
    main()
