#!/usr/bin/env python3
"""
extract_chi2_nominal_vs_mixed.py

Compare data vs nominal MC and data vs mixed (double-interaction blended) MC
by computing chi2/ndf for shower-shape variable projections.

Replicates the chi2 computation from plot_showershapes_variations.C:
  - Rebin the TH2 x-axis
  - Project X (sum over all Y bins)
  - Normalize to unit area
  - chi2 = sum_{bins in range} (d - m)^2 / (ed^2 + em^2),  ndf = nbins_contributing - 1
  - Skip bins where e2 <= 0 or m <= 0 (matching ROOT macro logic)

Output: formatted table to stdout + CSV to plotting/figures/di_comparison/chi2_nominal_vs_mixed.csv

Usage:
    python3 extract_chi2_nominal_vs_mixed.py                           # default suffixes
    python3 extract_chi2_nominal_vs_mixed.py --suffixes showershape    # single suffix
    python3 extract_chi2_nominal_vs_mixed.py --suffixes showershape showershape_0rad showershape_1p5mrad
"""

import argparse
import csv
import os
import sys

import numpy as np
import uproot

# ── Configuration ──────────────────────────────────────────────────────────────

BASE = "/sphenix/user/shuhangli/ppg12"
RESULTS = os.path.join(BASE, "efficiencytool", "results")

# Default suffix list for the DI comparison. Adjust with --suffixes.
# Display labels are auto-derived from the suffix (showershape_0rad → "0mrad",
# showershape_1p5mrad → "1.5mrad", showershape → "all").
DEFAULT_SUFFIXES = ["showershape", "showershape_0rad", "showershape_1p5mrad"]


def build_file_dict(suffix):
    return {
        "data":    os.path.join(RESULTS, f"data_histoshower_shape_{suffix}.root"),
        "nominal": os.path.join(RESULTS, f"MC_efficiencyshower_shape_jet_inclusive_{suffix}.root"),
        "mixed":   os.path.join(RESULTS, f"MC_efficiencyshower_shape_jet_inclusive_combined_{suffix}.root"),
    }


def display_label(suffix):
    if suffix.endswith("_0rad"):      return "0mrad"
    if suffix.endswith("_1p5mrad"):   return "1.5mrad"
    if suffix.endswith("_1p5rad"):    return "1.5mrad"
    if suffix == "showershape":       return "all"
    return suffix


VARIABLES = ["h2d_weta_cogx", "h2d_et1", "h2d_e11_to_e33", "h2d_bdt"]

CUTS = [0, 1, 2, 3]

PT_BINS = 5  # pt0..pt4
PT_EDGES = [10, 14, 18, 22, 28, 30]  # GeV

ETA_BIN = 0

# Axis settings matching getAxisSettings in plot_showershapes_variations.C
# Key is the variable name after stripping "h2d_" prefix
AXIS_SETTINGS = {
    "weta_cogx":  {"xmin": 0.0, "xmax": 2.0, "rebin": 4},
    "et1":        {"xmin": 0.3, "xmax": 1.0, "rebin": 1},
    "e11_to_e33": {"xmin": 0.0, "xmax": 1.0, "rebin": 4},
    "bdt":        {"xmin": 0.0, "xmax": 1.0, "rebin": 2},
}


# ── Helper functions ───────────────────────────────────────────────────────────

def get_projection_and_errors(th2, rebin):
    """
    From a TH2F (uproot object):
      1. Get values and variances (2D)
      2. Project onto X (sum over Y axis)
      3. Rebin the projection
      4. Return (rebinned_values, rebinned_errors, rebinned_edges)
    """
    vals, xedges, yedges = th2.to_numpy()

    # Get variances: from fSumw2 if available, otherwise Poisson (vals)
    sumw2_arr = np.array(th2.member("fSumw2"))
    if len(sumw2_arr) > 0:
        nx = len(xedges) - 1
        ny = len(yedges) - 1
        sumw2_2d = sumw2_arr.reshape(nx + 2, ny + 2)
        variances = sumw2_2d[1:-1, 1:-1]  # strip overflow
    else:
        variances = np.abs(vals)  # Poisson: variance = |N|

    # Project X: sum over Y (axis=1)
    proj_vals = vals.sum(axis=1)
    proj_vars = variances.sum(axis=1)

    # Rebin
    if rebin > 1:
        norigbins = len(proj_vals)
        nnewbins = norigbins // rebin
        # Truncate to exact multiple of rebin
        proj_vals = proj_vals[: nnewbins * rebin].reshape(nnewbins, rebin).sum(axis=1)
        proj_vars = proj_vars[: nnewbins * rebin].reshape(nnewbins, rebin).sum(axis=1)
        new_xedges = xedges[:: rebin][: nnewbins + 1]
    else:
        new_xedges = xedges

    proj_errs = np.sqrt(np.maximum(proj_vars, 0.0))
    return proj_vals, proj_errs, new_xedges


def compute_chi2(data_vals, data_errs, mc_vals, mc_errs, xedges, xmin, xmax):
    """
    Compute chi2/ndf between normalized data and MC projections.
    """
    centers = (xedges[:-1] + xedges[1:]) / 2.0
    mask = (centers >= xmin) & (centers <= xmax)

    d_sel = data_vals[mask].copy()
    m_sel = mc_vals[mask].copy()
    de_sel = data_errs[mask].copy()
    me_sel = mc_errs[mask].copy()

    # Normalize to unit area
    d_integral = d_sel.sum()
    m_integral = m_sel.sum()

    if d_integral <= 1e-12 or m_integral <= 1e-12:
        return np.nan, 0

    # Scale values and errors
    d_sel /= d_integral
    de_sel /= d_integral
    m_sel /= m_integral
    me_sel /= m_integral

    chi2 = 0.0
    ndf = 0
    for i in range(len(d_sel)):
        d = d_sel[i]
        m = m_sel[i]
        ed = de_sel[i]
        em = me_sel[i]
        e2 = ed * ed + em * em
        if e2 <= 0 or m <= 0:
            continue
        chi2 += (d - m) ** 2 / e2
        ndf += 1

    if ndf > 1:
        ndf -= 1
    else:
        return np.nan, 0

    return chi2 / ndf, ndf


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--suffixes", nargs="+", default=DEFAULT_SUFFIXES,
                        help=(f"Config suffix list (default: {DEFAULT_SUFFIXES}). "
                              "Each suffix S is expanded to the file triplet "
                              "(data_histoshower_shape_S.root, MC_..._jet_inclusive_S.root, "
                              "MC_..._jet_inclusive_combined_S.root)."))
    args = parser.parse_args()

    suffixes = args.suffixes

    # Open all files (skip missing files with a warning rather than hard-failing)
    open_files = {}
    for s in suffixes:
        paths = build_file_dict(s)
        missing = [p for p in paths.values() if not os.path.isfile(p)]
        if missing:
            for p in missing:
                print(f"[WARN] suffix={s}: file not found: {p}", file=sys.stderr)
            print(f"[WARN] suffix={s}: skipping due to missing file(s)", file=sys.stderr)
            continue
        open_files[s] = {tag: uproot.open(p) for tag, p in paths.items()}

    if not open_files:
        print("ERROR: no usable suffixes", file=sys.stderr)
        sys.exit(1)

    results = []

    for var in VARIABLES:
        varname = var.replace("h2d_", "")
        settings = AXIS_SETTINGS[varname]
        xmin = settings["xmin"]
        xmax = settings["xmax"]
        rebin = settings["rebin"]

        for icut in CUTS:
            for ipt in range(PT_BINS):
                ptlo = PT_EDGES[ipt]
                pthi = PT_EDGES[ipt + 1]

                for s in open_files:
                    histname = f"{var}_eta{ETA_BIN}_pt{ipt}_cut{icut}"

                    f_data = open_files[s]["data"]
                    f_nom = open_files[s]["nominal"]
                    f_mix = open_files[s]["mixed"]

                    # Check histogram exists
                    if histname not in f_data:
                        print(f"WARNING: {histname} not in data ({s}), skipping")
                        continue
                    if histname not in f_nom:
                        print(f"WARNING: {histname} not in nominal MC ({s}), skipping")
                        continue
                    if histname not in f_mix:
                        print(f"WARNING: {histname} not in mixed MC ({s}), skipping")
                        continue

                    h_data = f_data[histname]
                    h_nom = f_nom[histname]
                    h_mix = f_mix[histname]

                    d_vals, d_errs, d_edges = get_projection_and_errors(h_data, rebin)
                    n_vals, n_errs, n_edges = get_projection_and_errors(h_nom, rebin)
                    m_vals, m_errs, m_edges = get_projection_and_errors(h_mix, rebin)

                    chi2_nom, ndf_nom = compute_chi2(
                        d_vals, d_errs, n_vals, n_errs, d_edges, xmin, xmax
                    )
                    chi2_mix, ndf_mix = compute_chi2(
                        d_vals, d_errs, m_vals, m_errs, d_edges, xmin, xmax
                    )

                    # Percent change: (mixed - nominal) / nominal * 100
                    if not np.isnan(chi2_nom) and chi2_nom > 0:
                        pct = (chi2_mix - chi2_nom) / chi2_nom * 100.0
                    else:
                        pct = np.nan

                    results.append({
                        "variable": varname,
                        "cut": f"cut{icut}",
                        "pt_bin": f"pt{ipt} [{ptlo}-{pthi}]",
                        "crossing_angle": display_label(s),
                        "suffix": s,
                        "chi2_ndf_nominal": chi2_nom,
                        "ndf_nominal": ndf_nom,
                        "chi2_ndf_mixed": chi2_mix,
                        "ndf_mixed": ndf_mix,
                        "pct_change": pct,
                    })

    # ── Print formatted table ──────────────────────────────────────────────────
    header = (
        f"{'Variable':<14} {'Cut':<5} {'pT bin':<16} {'Angle':<9} "
        f"{'chi2/ndf nom':>14} {'chi2/ndf mix':>14} {'% change':>10}"
    )
    print()
    print("=" * len(header))
    print("Chi2/ndf comparison: Data vs Nominal MC  vs  Data vs Mixed (double-interaction) MC")
    print("=" * len(header))
    print(header)
    print("-" * len(header))

    for r in results:
        nom_str = f"{r['chi2_ndf_nominal']:.3f}" if not np.isnan(r["chi2_ndf_nominal"]) else "N/A"
        mix_str = f"{r['chi2_ndf_mixed']:.3f}" if not np.isnan(r["chi2_ndf_mixed"]) else "N/A"
        pct_str = f"{r['pct_change']:+.1f}%" if not np.isnan(r["pct_change"]) else "N/A"

        print(
            f"{r['variable']:<14} {r['cut']:<5} {r['pt_bin']:<16} {r['crossing_angle']:<9} "
            f"{nom_str:>14} {mix_str:>14} {pct_str:>10}"
        )

    print("-" * len(header))
    print(f"Total entries: {len(results)}")
    print()

    # ── Write CSV ──────────────────────────────────────────────────────────────
    csv_path = os.path.join(BASE, "plotting", "figures", "di_comparison", "chi2_nominal_vs_mixed.csv")
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "variable", "cut", "pt_bin", "crossing_angle", "suffix",
                "chi2_ndf_nominal", "ndf_nominal",
                "chi2_ndf_mixed", "ndf_mixed",
                "pct_change",
            ],
        )
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    print(f"CSV written to: {csv_path}")


if __name__ == "__main__":
    main()
