#!/usr/bin/env python3
"""Plot the relative difference between two result variants.

Usage:
    python plot_var_comparison.py VAR1 VAR2 [options]

    VAR1, VAR2  — var_type strings, e.g. bdt_nom, bdt_tightbdt50
                  These must correspond to Photon_final_{VAR}.root files
                  in the results directory.

Examples:
    python plot_var_comparison.py bdt_nom bdt_tightbdt50
    python plot_var_comparison.py bdt_nom bdt_energyscale26 --label1 "nominal" --label2 "escale +2.6%"
    python plot_var_comparison.py bdt_nom bdt_0rad --outdir rootFiles --figdir figures
"""

import argparse
import os
import sys

import ROOT

# ---------------------------------------------------------------------------
# Load sPHENIX style from plotcommon.h (same dir as this script)
# ---------------------------------------------------------------------------
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(_THIS_DIR, "plotcommon.h"))


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("var1", help="Reference var_type (denominator), e.g. bdt_nom")
    p.add_argument("var2", help="Comparison var_type (numerator), e.g. bdt_tightbdt50")
    p.add_argument("--results", default="/sphenix/user/shuhangli/ppg12/efficiencytool/results",
                   help="Directory containing Photon_final_*.root files")
    p.add_argument("--histogram", default="h_unfold_sub_result",
                   help="Histogram name inside each ROOT file (default: h_unfold_sub_result)")
    p.add_argument("--figdir", default="figures",
                   help="Output directory for the PDF plot (default: figures/)")
    p.add_argument("--label1", default=None,
                   help="Legend label for VAR1 (default: VAR1 string)")
    p.add_argument("--label2", default=None,
                   help="Legend label for VAR2 (default: VAR2 string)")
    p.add_argument("--output", default=None,
                   help="Output PDF filename (default: figures/compare_{VAR1}_vs_{VAR2}.pdf)")
    return p.parse_args()


def load_spectrum(var_type, results_dir, histogram):
    fname = os.path.join(results_dir, f"Photon_final_{var_type}.root")
    f = ROOT.TFile.Open(fname, "READ")
    if not f or f.IsZombie():
        sys.exit(f"ERROR: cannot open {fname}")
    h = f.Get(histogram)
    if not h:
        f.Close()
        sys.exit(f"ERROR: histogram '{histogram}' not found in {fname}")
    h = h.Clone(f"h_{var_type}")
    h.SetDirectory(0)
    f.Close()
    return h


def yrange_from_hist(h, padding=0.3):
    """Auto Y range: scans all bins, anchors at 0, adds padding."""
    vmax, vmin = 0.0, 0.0
    for i in range(1, h.GetNbinsX() + 1):
        v = h.GetBinContent(i)
        if v > vmax:
            vmax = v
        if v < vmin:
            vmin = v
    span = max(vmax - vmin, 0.04)
    margin = span * padding
    return vmin - margin, vmax + margin


def main():
    args = parse_args()
    os.makedirs(args.figdir, exist_ok=True)

    label1 = args.label1 if args.label1 else args.var1
    label2 = args.label2 if args.label2 else args.var2

    outfile = args.output if args.output else \
        os.path.join(args.figdir, f"compare_{args.var1}_vs_{args.var2}.pdf")

    ROOT.init_plot()

    h1 = load_spectrum(args.var1, args.results, args.histogram)
    h2 = load_spectrum(args.var2, args.results, args.histogram)

    # relative difference: (h2 - h1) / h1, bin-by-bin
    h_rel = h2.Clone("h_rel")
    h_rel.SetDirectory(0)
    nbins = h1.GetNbinsX()
    for i in range(1, nbins + 1):
        denom = h1.GetBinContent(i)
        if denom != 0:
            h_rel.SetBinContent(i, (h2.GetBinContent(i) - denom) / denom)
        else:
            h_rel.SetBinContent(i, 0.0)
        h_rel.SetBinError(i, 0.0)

    ylo, yhi = yrange_from_hist(h_rel)

    c = ROOT.TCanvas("c_compare", "", 800, 600)
    c.cd()
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.08)
    c.SetTopMargin(0.07)
    c.SetBottomMargin(0.13)

    ROOT.frame_et_truth.SetYTitle(f"({label2} #minus {label1}) / {label1}")
    ROOT.frame_et_truth.GetYaxis().SetNdivisions(506)
    ROOT.frame_et_truth.GetYaxis().SetRangeUser(ylo, yhi)
    ROOT.frame_et_truth.GetXaxis().SetRangeUser(10, 30)
    ROOT.frame_et_truth.GetXaxis().SetNdivisions(505)
    ROOT.frame_et_truth.SetXTitle("#it{E}_{T}^{#gamma} [GeV]")
    ROOT.frame_et_truth.Draw("axis")

    h_rel.SetMarkerStyle(0)
    h_rel.SetLineColor(ROOT.kBlack)
    h_rel.SetLineWidth(2)
    h_rel.Draw("same ][ HIST")

    ROOT.linezero.Draw("L")

    ROOT.myText(0.55, 0.88, 1, ROOT.strleg1.c_str(), 0.045)
    ROOT.myText(0.55, 0.83, 1, ROOT.strleg2.c_str(), 0.045)
    ROOT.myText(0.18, 0.88, 1, f"ref: {label1}", 0.045)
    ROOT.myText(0.18, 0.83, 1, f"cmp: {label2}", 0.045)

    c.SaveAs(outfile)
    c.Close()
    print(f"Saved: {outfile}")


if __name__ == "__main__":
    main()
