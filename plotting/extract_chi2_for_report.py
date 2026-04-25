#!/usr/bin/env python3
"""
Extract chi2/ndf and p-values for the double interaction shower shape report.

Computes chi2/ndf between unit-normalized data and inclusive MC histograms,
matching the procedure in plot_showershapes_variations.C:
  - TH2F -> RebinX -> SetRangeUser -> ProjectionX -> Sumw2 -> Scale(1/integral)
  - chi2 = sum_bins (d - m)^2 / (ed^2 + em^2), skip bins with e2<=0 or m<=0
  - ndf = count - 1
  - p-value = scipy.stats.chi2.sf(chi2, ndf)

Usage:
    python3 extract_chi2_for_report.py                                 # default suffixes
    python3 extract_chi2_for_report.py --suffixes showershape          # single suffix
    python3 extract_chi2_for_report.py --suffixes showershape showershape_0rad showershape_1p5mrad
"""

import argparse

import numpy as np
import uproot
from scipy.stats import chi2 as chi2_dist

# ---- Configuration ----

RESULTS = "/sphenix/user/shuhangli/ppg12/efficiencytool/results"

# Default suffix list. "showershape" = all-range bare; "_0rad"/"_1p5mrad" = per-angle feeders.
DEFAULT_SUFFIXES = ["showershape", "showershape_0rad", "showershape_1p5mrad"]


def build_files(suffix):
    return {
        "data":     f"{RESULTS}/data_histoshower_shape_{suffix}.root",
        "nom_incl": f"{RESULTS}/MC_efficiencyshower_shape_jet_inclusive_{suffix}.root",
        "mix_incl": f"{RESULTS}/MC_efficiencyshower_shape_jet_inclusive_combined_{suffix}.root",
    }


def display_label(suffix):
    if suffix.endswith("_0rad"):      return "0mrad"
    if suffix.endswith("_1p5mrad"):   return "1.5mrad"
    if suffix.endswith("_1p5rad"):    return "1.5mrad"
    if suffix == "showershape":       return "all"
    return suffix


VARIABLES = {
    "weta_cogx":   {"xmin": 0.0, "xmax": 2.0, "rebin": 4},
    "bdt":         {"xmin": 0.0, "xmax": 1.0, "rebin": 2},
    "et1":         {"xmin": 0.3, "xmax": 1.0, "rebin": 1},
    "e11_to_e33":  {"xmin": 0.0, "xmax": 1.0, "rebin": 4},
}

PT_BINS = [0, 1, 3]
CUT = 1

# From plotcommon.h
PT_EDGES = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]


def project_and_prepare(rootfile, histname, rebin, xmin, xmax):
    """Read TH2F, project to X axis, rebin, restrict range, normalize to unit area."""
    h2 = rootfile[histname]
    vals2d = h2.values()       # shape (nx, ny)
    vars2d = h2.variances()    # shape (nx, ny), sum-of-weights-squared per bin
    xedges = h2.axes[0].edges()

    # ProjectionX: sum over Y axis
    proj_vals = np.sum(vals2d, axis=1)
    proj_vars = np.sum(vars2d, axis=1)

    # Rebin X
    if rebin > 1:
        n = len(proj_vals)
        n_new = n // rebin
        new_vals = np.zeros(n_new)
        new_vars = np.zeros(n_new)
        new_edges = np.zeros(n_new + 1)
        for i in range(n_new):
            sl = slice(i * rebin, (i + 1) * rebin)
            new_vals[i] = np.sum(proj_vals[sl])
            new_vars[i] = np.sum(proj_vars[sl])
            new_edges[i] = xedges[i * rebin]
        new_edges[n_new] = xedges[n_new * rebin]
        proj_vals = new_vals
        proj_vars = new_vars
        xedges = new_edges

    centers = 0.5 * (xedges[:-1] + xedges[1:])

    # Restrict to [xmin, xmax] by bin center
    mask = (centers >= xmin) & (centers <= xmax)
    proj_vals = proj_vals[mask]
    proj_vars = proj_vars[mask]
    centers = centers[mask]

    # Unit normalization
    integral = np.sum(proj_vals)
    if integral > 1e-12:
        proj_vars = proj_vars / (integral * integral)
        proj_vals = proj_vals / integral

    errors = np.sqrt(np.maximum(proj_vars, 0.0))
    return proj_vals, errors, centers


def compute_chi2(data_vals, data_errs, mc_vals, mc_errs):
    """Compute chi2 and ndf matching plot_showershapes_variations.C."""
    chi2_sum = 0.0
    count = 0
    for d, m, ed, em in zip(data_vals, mc_vals, data_errs, mc_errs):
        e2 = ed * ed + em * em
        if e2 <= 0 or m <= 0:
            continue
        chi2_sum += (d - m) ** 2 / e2
        count += 1
    if count > 1:
        ndf = count - 1
        pval = float(chi2_dist.sf(chi2_sum, ndf))
        return chi2_sum, ndf, pval
    return None, None, None


def main():
    import os
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--suffixes", nargs="+", default=DEFAULT_SUFFIXES,
                        help=f"Config suffix list (default: {DEFAULT_SUFFIXES}).")
    args = parser.parse_args()

    # Open all ROOT files; skip suffixes whose files are missing
    open_files = {}
    for s in args.suffixes:
        files = build_files(s)
        missing = [p for p in files.values() if not os.path.isfile(p)]
        if missing:
            for p in missing:
                print(f"[WARN] suffix={s}: missing {p}")
            continue
        open_files[s] = {k: uproot.open(p) for k, p in files.items()}

    if not open_files:
        print("ERROR: no usable suffixes")
        return

    # Print header
    header = (
        f"{'variable':<14s} {'xangle':>8s} {'pt':>4s} "
        f"| {'chi2_nom':>9s} {'ndf':>4s} {'c2/n_nom':>9s} {'pval_nom':>9s} "
        f"| {'chi2_mix':>9s} {'ndf':>4s} {'c2/n_mix':>9s} {'pval_mix':>9s}"
    )
    sep = "-" * len(header)
    print(header)
    print(sep)

    results = []

    for var, settings in VARIABLES.items():
        xmin = settings["xmin"]
        xmax = settings["xmax"]
        rebin = settings["rebin"]

        for s in open_files:
            for ipt in PT_BINS:
                histname = f"h2d_{var}_eta0_pt{ipt}_cut{CUT}"
                pt_label = f"pt{ipt}"
                pt_range = f"{PT_EDGES[ipt]}-{PT_EDGES[ipt+1]}"

                d_vals, d_errs, _ = project_and_prepare(
                    open_files[s]["data"], histname, rebin, xmin, xmax)
                n_vals, n_errs, _ = project_and_prepare(
                    open_files[s]["nom_incl"], histname, rebin, xmin, xmax)
                m_vals, m_errs, _ = project_and_prepare(
                    open_files[s]["mix_incl"], histname, rebin, xmin, xmax)

                c2_nom, ndf_nom, pv_nom = compute_chi2(d_vals, d_errs, n_vals, n_errs)
                c2_mix, ndf_mix, pv_mix = compute_chi2(d_vals, d_errs, m_vals, m_errs)

                # Format output
                def fmtrow(c2, ndf, pv):
                    if c2 is None:
                        return "      ---", " ---", "      ---", "      ---"
                    return (f"{c2:9.1f}", f"{ndf:4d}",
                            f"{c2/ndf:9.2f}", f"{pv:9.4f}")

                cn, nn, rn, pn = fmtrow(c2_nom, ndf_nom, pv_nom)
                cm, nm, rm, pm = fmtrow(c2_mix, ndf_mix, pv_mix)

                xangle_disp = display_label(s)
                print(f"{var:<14s} {xangle_disp:>8s} {pt_label:>4s} "
                      f"| {cn} {nn} {rn} {pn} "
                      f"| {cm} {nm} {rm} {pm}")

                results.append({
                    "variable": var, "xangle": xangle_disp, "suffix": s,
                    "pt_bin": pt_label, "pt_range": pt_range,
                    "chi2_nom": c2_nom, "ndf_nom": ndf_nom, "pval_nom": pv_nom,
                    "chi2_mix": c2_mix, "ndf_mix": ndf_mix, "pval_mix": pv_mix,
                })

    # Close files
    for s in open_files:
        for key in open_files[s]:
            open_files[s][key].close()

    # LaTeX table
    print("\n\n=== LaTeX table rows ===\n")
    print(r"\begin{tabular}{llc|rrc|rrc}")
    print(r"\toprule")
    print(r"Variable & Crossing & $p_T$ [GeV] "
          r"& $\chi^2$/ndf (nominal) & $p$-value "
          r"& $\chi^2$/ndf (mixed) & $p$-value \\")
    print(r"\midrule")

    for r in results:
        if r["chi2_nom"] is None:
            continue
        c2n, nn, pn = r["chi2_nom"], r["ndf_nom"], r["pval_nom"]
        c2m, nm, pm = r["chi2_mix"], r["ndf_mix"], r["pval_mix"]
        var_tex = r["variable"].replace("_", r"\_")
        xa = r["xangle"]
        print(f"  {var_tex} & {xa} & {r['pt_range']} "
              f"& {c2n:.1f}/{nn} = {c2n/nn:.2f} & {pn:.4f} "
              f"& {c2m:.1f}/{nm} = {c2m/nm:.2f} & {pm:.4f} \\\\")

    print(r"\bottomrule")
    print(r"\end{tabular}")


if __name__ == "__main__":
    main()
