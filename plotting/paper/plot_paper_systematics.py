"""plot_paper_systematics.py -- paper Fig.~\\ref{fig:syst_sum_rel}

Lifts the plot_breakdown() canvas from plotting/calc_syst_bdt.py and writes
it directly into PPG12-Paper/figures/syst_sum_rel.pdf. Aggregation logic is
imported from the parent script so the paper figure stays in sync with the
analysis-note one as long as the input ROOT files are unchanged. Cosmetic
edits live in this file and do NOT touch calc_syst_bdt.py.

Run:
    cd plotting/paper
    python plot_paper_systematics.py            # uses default results dir
    python plot_paper_systematics.py --results /path/to/results
"""
from __future__ import annotations

import argparse
import os
import sys

import ROOT

# Allow `from calc_syst_bdt import ...` even when run from plotting/paper/.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
sys.path.insert(0, _PARENT_DIR)
# Also expose efficiencytool/ for the make_bdt_variations import inside calc_syst_bdt.
sys.path.insert(0, os.path.join(os.path.dirname(_PARENT_DIR), "efficiencytool"))

# Load plotcommon.h via the interpreter (no ACLiC compile -- the header
# defines globals at file scope which are not compatible with +g).
ROOT.gROOT.ProcessLine(f'.L {os.path.join(_PARENT_DIR, "plotcommon.h")}')

from calc_syst_bdt import (  # noqa: E402
    SYST_TYPES,
    SYST_GROUPS,
    FINAL_SYSTS,
    FLAT_SYSTS,
    GROUP_COLORS,
    aggregate_type,
    add_flat,
    build_variant_map,
    load_spectrum,
    quadrature_sum,
    _build_stat_rel,
)

# Paper-specific output location and filename.
_PAPER_FIGDIR = os.path.normpath(
    os.path.join(_PARENT_DIR, "..", "PPG12-Paper", "figures"))
_PAPER_NAME = "syst_sum_rel.pdf"


# ----------------------------------------------------------------------
# plot_breakdown copy with paper-only cosmetics. Edit freely.
# ----------------------------------------------------------------------
def plot_breakdown_paper(group_results: dict, total: tuple,
                         args: argparse.Namespace) -> None:
    ROOT.init_plot()

    c = ROOT.TCanvas("c_breakdown_paper", "", 900, 600)

    ROOT.frame_et_truth.GetYaxis().SetRangeUser(-0.58, 0.5)
    # Paper x-range: 12 < ETg < 32 GeV (the reported analysis range).
    ROOT.frame_et_truth.GetXaxis().SetRangeUser(12.0, 32.0)
    ROOT.frame_et_truth.SetXTitle("#it{E}_{T}^{#gamma} [GeV]")
    ROOT.frame_et_truth.SetYTitle("Relative difference")
    ROOT.frame_et_truth.Draw("axis")

    keep_alive = []  # PyROOT GC guard

    # Stat-uncertainty band intentionally NOT drawn on the paper figure --
    # the journal draft shows the systematic envelope only.

    h_rl_tot, h_rh_tot = total[2], total[3]
    h_rl_tot_plot = h_rl_tot.Clone("h_rl_tot_plot")
    h_rl_tot_plot.SetDirectory(0)
    h_rl_tot_plot.Scale(-1)
    keep_alive.append(h_rl_tot_plot)
    for h in (h_rl_tot_plot, h_rh_tot):
        h.SetMarkerStyle(20)
        h.SetLineColor(ROOT.kBlack)
        h.SetLineWidth(4)
        h.SetMarkerColor(ROOT.kBlack)
    h_rl_tot_plot.Draw("same ][ HIST")
    h_rh_tot.Draw("same ][ HIST")

    ROOT.myMarkerLineText(0.25, 0.36, 0, ROOT.kBlack, 20, ROOT.kBlack, 0,
                          "Total", 0.05, True)

    # Paper legend order and human-readable labels. Drawing order is
    # decoupled from legend order so that escale stays on top of the
    # other bands (drawn last) while the legend reads in the order
    # the user wants.
    GROUP_ORDER = ["escale", "eres", "purity", "efficiency",
                   "di_fraction", "npb", "unfolding"]
    GROUP_LABELS = {
        "escale":      "Energy scale",
        "eres":        "Energy resolution",
        "purity":      "Purity",
        "efficiency":  "Efficiency",
        "di_fraction": "Pileup",
        "npb":         "Non-collisional background",
        "unfolding":   "Unfolding",
    }
    legend_index = {grp: i for i, grp in enumerate(GROUP_ORDER)}

    present = [g for g in GROUP_ORDER if g in group_results]
    # escale drawn last so its band sits on top of the others.
    draw_order = [g for g in present if g != "escale"] + \
                 [g for g in present if g == "escale"]

    switchover = 3
    for grp in draw_order:
        h_al, h_ah, h_rl, h_rh = group_results[grp]
        col = GROUP_COLORS.get(grp, ROOT.kGray)

        h_rl_plot = h_rl.Clone(f"h_rl_plot_{grp}")
        h_rl_plot.SetDirectory(0)
        h_rl_plot.Scale(-1)
        keep_alive.append(h_rl_plot)

        h_rl_plot.SetMarkerStyle(24)
        h_rl_plot.SetLineStyle(1)
        h_rl_plot.SetLineColor(col)
        h_rl_plot.SetMarkerColor(col)
        h_rl_plot.SetLineWidth(3 if grp == "escale" else 2)

        h_rh.SetMarkerStyle(20)
        h_rh.SetLineStyle(1)
        h_rh.SetLineColor(col)
        h_rh.SetMarkerColor(col)
        h_rh.SetLineWidth(3 if grp == "escale" else 2)

        h_rh.Draw("same ][ HIST")
        h_rl_plot.Draw("same ][ HIST")

        # Text position uses the user-requested legend order, NOT draw
        # order, so escale stays at the top of the legend.
        idx = legend_index[grp]
        xshift = 0.33 if idx >= switchover else 0.0
        yshift = 0.05 * (switchover + 1) if idx >= switchover else 0.0
        label  = GROUP_LABELS.get(grp, grp.replace("_", " "))
        ROOT.myMarkerLineText(0.25 + xshift, 0.31 - 0.05 * idx + yshift,
                              0, col, 20, col, 0,
                              label, 0.05, True)

    ROOT.myText(0.20, 0.88, 1, ROOT.strleg1.c_str(), 0.05)
    ROOT.myText(0.20, 0.83, 1, ROOT.strleg2.c_str(), 0.05)

    out = os.path.join(_PAPER_FIGDIR, _PAPER_NAME)
    c.SaveAs(out)
    c.Close()
    print(f"[paper-syst] wrote {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results",
                   default="/sphenix/user/shuhangli/ppg12/efficiencytool/results")
    p.add_argument("--outdir",
                   default=os.path.join(_PARENT_DIR, "rootFiles"),
                   help="Working dir for syst ROOT files (intermediate).")
    p.add_argument("--nom", default="bdt_nom")
    p.add_argument("--histogram", default="h_unfold_sub_result")
    p.add_argument("--skip-missing", action="store_true", default=True)
    args = p.parse_args()
    # plot_breakdown / _build_stat_rel expect args.figdir to exist for the
    # stat-only intermediate it writes; map it to the paper output dir.
    args.figdir = _PAPER_FIGDIR
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    print(f"[paper-syst] loading nominal: Photon_final_{args.nom}.root")
    h_nom = load_spectrum(args.nom, args.results, args.histogram)

    vmap = build_variant_map()

    type_results = {}
    for type_name in SYST_TYPES:
        try:
            result = aggregate_type(type_name, vmap, h_nom,
                                    args.results, args.histogram,
                                    args.skip_missing)
        except (FileNotFoundError, KeyError, OSError) as e:
            print(f"  [ERROR] {type_name}: {e}")
            continue
        if result is None:
            continue
        type_results[type_name] = result

    group_results = {}
    for grp, members in SYST_GROUPS.items():
        comps = [type_results[t] for t in members if t in type_results]
        if not comps:
            continue
        group_results[grp] = quadrature_sum(comps)

    total_comps = [group_results[g] for g in FINAL_SYSTS if g in group_results]
    if not total_comps:
        print("[paper-syst] ERROR: no group results -- nothing to plot")
        return 1
    total = add_flat(*quadrature_sum(total_comps), h_nom, FLAT_SYSTS)

    plot_breakdown_paper(group_results, total, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
