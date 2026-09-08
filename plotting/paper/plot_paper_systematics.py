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


def _legend_line_text(x, y, color, width, text, tsize=0.05, keep=None, lstyle=1):
    """Legend swatch: a horizontal line of `width` plus a black text label.

    Mirrors the msize==0 branch of BlairUtils myMarkerLineText (swatch line
    to the left of the label) but exposes the line width -- myMarkerLineText
    hardcodes it to 2, which is why the legend could not match the per-curve
    histogram line widths. BlairUtils.C is shared, so the override lives here.
    """
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextAlign(12)
    lat.SetTextSize(tsize)
    lat.DrawLatex(x, y, text)
    ln = ROOT.TLine()
    ln.SetLineColor(color)
    ln.SetLineStyle(lstyle)
    ln.SetLineWidth(width)
    ln.DrawLineNDC(x - 0.95 * tsize, y, x - 0.15 * tsize, y)
    if keep is not None:
        keep.extend((lat, ln))


def _flatten_left_edge(h, x_below=11.0, x_first=13.0):
    """Copy the first-visible (12-14 GeV) bin value into the bin just left of
    the 12 GeV axis start (the 10-12 GeV bin).

    The figure x-range begins exactly at 12 GeV, so the step histogram would
    otherwise draw a vertical riser on the y-axis between the (off-frame)
    10-12 bin and the 12-14 bin. Equalising the two makes that riser
    zero-height, so each curve enters the frame flat with no spurious edge
    line. Cosmetic only -- the off-frame 10-12 bin is never reported.
    """
    b_below = h.FindBin(x_below)
    b_first = h.FindBin(x_first)
    if b_below >= 1 and b_first >= 1:
        h.SetBinContent(b_below, h.GetBinContent(b_first))


# ----------------------------------------------------------------------
# plot_breakdown copy with paper-only cosmetics. Edit freely.
# ----------------------------------------------------------------------
def plot_breakdown_paper(group_results: dict, total: tuple,
                         args: argparse.Namespace) -> None:
    ROOT.init_plot()

    # Paper-only colour overrides (cosmetic): Pileup (di_fraction, was teal
    # kCyan+2) read green-ish next to Energy scale (green), and Unfolding
    # (was magenta kMagenta+2) read violet next to Purity (violet). Recolour
    # them to neutral gray / brown so all seven groups are separable. The
    # shared calc_syst_bdt.GROUP_COLORS is left untouched.
    paper_colors = dict(GROUP_COLORS)
    paper_colors["efficiency"] = ROOT.kRed -4
    paper_colors["di_fraction"] = ROOT.kMagenta +1
    # paper_colors["di_fraction"] = ROOT.TColor.GetColor("#7f7f7f")  # gray  (Pileup)
    paper_colors["unfolding"]   = ROOT.TColor.GetColor("#8c564b")  # brown (Unfolding)

    c = ROOT.TCanvas("c_breakdown_paper", "", 900, 600)

    ROOT.frame_et_truth.GetYaxis().SetRangeUser(-0.48, 0.36)
    # Paper x-range: 12 < ETg < 32 GeV (the reported analysis range).
    # Slight inset (12.1, 31.9) to suppress bin-edge vertical-line artifacts
    # flagged by Jamie on Fig. 6.
    ROOT.frame_et_truth.GetXaxis().SetRangeUser(12.0, 31.9)
    ROOT.frame_et_truth.SetXTitle("#it{E}_{T}^{#gamma} [GeV]")
    ROOT.frame_et_truth.SetYTitle("Relative Uncertainty")
    ROOT.frame_et_truth.GetYaxis().SetTitleOffset(1.1)
    ROOT.frame_et_truth.GetXaxis().SetTitleOffset(1.0)
    ROOT.frame_et_truth.GetYaxis().SetTitleSize(0.058)
    ROOT.frame_et_truth.GetXaxis().SetTitleSize(0.057)
    ROOT.frame_et_truth.Draw("axis")

    keep_alive = []  # PyROOT GC guard

    # Stat-uncertainty band intentionally NOT drawn on the paper figure --
    # the journal draft shows the systematic envelope only.

    h_rl_tot, h_rh_tot = total[2], total[3]
    h_rl_tot_plot = h_rl_tot.Clone("h_rl_tot_plot")
    h_rl_tot_plot.SetDirectory(0)
    h_rl_tot_plot.Scale(-1)
    keep_alive.append(h_rl_tot_plot)
    LW_TOTAL = 3  # Total histogram line width; legend swatch matches it
    for h in (h_rl_tot_plot, h_rh_tot):
        h.SetMarkerStyle(20)
        h.SetLineColor(ROOT.kGray+3)
        h.SetLineWidth(LW_TOTAL)
        h.SetMarkerColor(ROOT.kGray+3)
        _flatten_left_edge(h)
    h_rl_tot_plot.Draw("same ][ HIST")
    h_rh_tot.Draw("same ][ HIST")

    _legend_line_text(0.25, 0.36, ROOT.kGray+3, LW_TOTAL, "Total", 0.05, keep_alive)

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
    # Per-group ROOT line style -- applied to BOTH the histogram curve and its
    # legend swatch so they stay in sync (same idea as the per-group width).
    # escale kept solid (it is the dominant band); the rest get distinct
    # dashed/dotted patterns so the groups are separable in B/W too.
    # 1=solid 2=dashed 3=dotted 4=dash-dot 7=long-dash 8=dash-dot-dot
    # 10=long-dash-dot. (Style 9 "80 20" avoided -- its 80-unit dash fills the
    # short legend swatch solid; 8 keeps a visible gap there.)
    GROUP_STYLES = {
        "escale":      1,
        "eres":        1,
        "purity":      1,
        "efficiency":  4,
        "di_fraction": 7,
        "npb":         3,
        "unfolding":   2,
    }
    legend_index = {grp: i for i, grp in enumerate(GROUP_ORDER)}

    present = [g for g in GROUP_ORDER if g in group_results]
    # escale drawn last so its band sits on top of the others.
    draw_order = [g for g in present if g != "escale"] + \
                 [g for g in present if g == "escale"]

    switchover = 3
    for grp in draw_order:
        h_al, h_ah, h_rl, h_rh = group_results[grp]
        col = paper_colors.get(grp, ROOT.kGray)

        h_rl_plot = h_rl.Clone(f"h_rl_plot_{grp}")
        h_rl_plot.SetDirectory(0)
        h_rl_plot.Scale(-1)
        keep_alive.append(h_rl_plot)

        # Group line width + style -- used for BOTH the histogram and its
        # legend swatch (below) so the two always match. Bump escale width
        # here to emphasise it; the legend tracks automatically.
        lw = 2
        # lw = 3 if grp == "escale" else 2
        ls = GROUP_STYLES.get(grp, 1)
        h_rl_plot.SetMarkerStyle(24)
        h_rl_plot.SetLineStyle(ls)
        h_rl_plot.SetLineColor(col)
        h_rl_plot.SetMarkerColor(col)
        h_rl_plot.SetLineWidth(lw)

        h_rh.SetMarkerStyle(20)
        h_rh.SetLineStyle(ls)
        h_rh.SetLineColor(col)
        h_rh.SetMarkerColor(col)
        h_rh.SetLineWidth(lw)

        _flatten_left_edge(h_rl_plot)
        _flatten_left_edge(h_rh)

        h_rh.Draw("same ][ HIST")
        h_rl_plot.Draw("same ][ HIST")

        # Text position uses the user-requested legend order, NOT draw
        # order, so escale stays at the top of the legend.
        idx = legend_index[grp]
        xshift = 0.31 if idx >= switchover else 0.0
        yshift = 0.05 * (switchover + 1) if idx >= switchover else 0.0
        label  = GROUP_LABELS.get(grp, grp.replace("_", " "))
        _legend_line_text(0.25 + xshift, 0.31 - 0.05 * idx + yshift,
                          col, lw, label, 0.05, keep_alive, lstyle=ls)

    ROOT.myText(0.20, 0.86, 1, "#bf{#it{sPHENIX}}", 0.062)
    # ROOT.myText(0.20, 0.88, 1, ROOT.strleg1.c_str(), 0.05)
    ROOT.myText(0.91, 0.88, 1, ROOT.strleg2.c_str(), 0.05, 1)
    # ROOT.myText(0.20, 0.83, 1, ROOT.strleg2.c_str(), 0.05)

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
