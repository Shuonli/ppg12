#!/usr/bin/env python3
"""
Compute the per-PDF 68% CL uncertainty band from the per-member histograms
produced by MakeJetPHOXhisto_pdfmem.C.

Inputs (one per PDF):
    rootFiles/jetPHOX_<pdf>_10_pdfmem.root
        h_truth_pT_m0, h_truth_pT_m1, ... (one per LHAPDF member)

Output:
    rootFiles/jetPHOX_pdfunc.root
        For each PDF tag <pdf>:
            h_pdf_central_<pdf>     TH1F  central member (or replica mean)
            g_pdf_band_<pdf>        TGraphAsymmErrors  (central, +PDF, -PDF)

Recipes (PDF4LHC15):

    Hessian:
        Δσ^+ = (1/C) * sqrt( Σ_pairs [ max(σ_{2i-1} - σ_0, σ_{2i} - σ_0, 0) ]^2 )
        Δσ^- = (1/C) * sqrt( Σ_pairs [ max(σ_0 - σ_{2i-1}, σ_0 - σ_{2i}, 0) ]^2 )
        Pairs are members 1..2*N_pairs ; member 0 is central.
        For CT/CTEQ at 90% CL: C = 1.645  (convert to 68% CL).
        For MSHT (already 68% CL): C = 1.

    MC replicas (NNPDF):
        σ_0 (reported)  = (1/N) * Σ_r σ_r       (replica mean)
        Δσ              = sqrt( (1/(N-1)) * Σ_r (σ_r - <σ>)^2 )
        N = 100, replicas are members 1..100 ; member 0 is the
        published "central" replica (we use the mean, per NNPDF
        recommendation).

PDF-specific config (member layout, CL conversion factor) is encoded in
PDF_INFO below. CT18NLO has 28 eigenvector pairs (members 1..56) plus 2
alpha_s variants (members 57, 58); the alpha_s pair is excluded from the
PDF eigenvector band by default and would be combined separately if
needed.
"""

import argparse
import math
import sys
from array import array

import uproot
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)

PDF_INFO = {
    # tag       (kind,           n_pairs_or_dirs,  C,      n_total_members)
    "ct18":    ("hessian",       28,               1.645,  59),  # 1 + 56 pairs + 2 alpha_s
    "nlo":     ("hessian",       28,               1.645,  57),  # CT14NLO: 1 + 56 pairs
    "cteq":    ("hessian",       22,               1.645,  45),  # 1 + 44 pairs
    "msht":    ("hessian",       32,               1.0,    65),  # 1 + 64 pairs (68% CL native)
    "nnpdf":   ("replicas",      100,              None,   101), # 1 central + 100 replicas
    "nnpdf4":  ("replicas",      100,              None,   101), # 1 central + 100 replicas
    "nnpdfh":  ("symmhessian",   100,              1.0,    101), # 1 central + 100 eigendirs (68% CL)
    "nnpdf4h": ("symmhessian",   50,               1.0,    51),  # NNPDF4.0 NNLO Hessian: 1 + 50 eigendirs (68% CL)
    "ct18nnlo": ("hessian",      28,               1.645,  59),  # CT18NNLO: 1 + 56 pairs + 2 alpha_s
    "pdf4lhc21": ("symmhessian", 40,               1.0,    41),  # PDF4LHC21_40 (NNLO): 1 + 40 eigendirs (68% CL)
}

ROOT_FILE_FMT = "rootFiles/jetPHOX_{tag}_10_pdfmem.root"
OUT_FILE = "rootFiles/jetPHOX_pdfunc.root"


def read_members(path, n_members):
    """Load all member histograms from an _pdfmem.root file.
    Returns: edges (n+1,), values (n_members, n_bins), errors (n_members, n_bins)
    where errors are the Sumw2 MC uncertainties per bin."""
    f = uproot.open(path)
    edges = None
    vals = []
    errs = []
    for m in range(n_members):
        key = f"h_truth_pT_m{m}"
        if key not in f:
            raise KeyError(f"missing {key} in {path}")
        h = f[key]
        if edges is None:
            edges = h.axis().edges()
        vals.append(h.values())
        errs.append(h.errors())
    return edges, np.array(vals), np.array(errs)


def compute_hessian_band(vals, n_pairs, C):
    """Return (central, dplus, dminus) per bin for a Hessian PDF.
    vals shape: (1 + 2*n_pairs + extras, n_bins)."""
    central = vals[0]
    n_bins = vals.shape[1]
    dplus = np.zeros(n_bins)
    dminus = np.zeros(n_bins)
    for p in range(n_pairs):
        sp = vals[2 * p + 1]
        sm = vals[2 * p + 2]
        dp = np.maximum.reduce([sp - central, sm - central, np.zeros(n_bins)])
        dm = np.maximum.reduce([central - sp, central - sm, np.zeros(n_bins)])
        dplus += dp ** 2
        dminus += dm ** 2
    dplus = np.sqrt(dplus) / C
    dminus = np.sqrt(dminus) / C
    return central, dplus, dminus


def compute_symmhessian_band(vals, n_dirs, C):
    """Return (central, dplus, dminus) per bin for a symmetric Hessian PDF
    (NNPDF Hessian-converted variant). Each member 1..n_dirs is one
    eigendirection; the band is symmetric.
        Δσ = (1/C) * sqrt( Σ_i (σ_i - σ_0)^2 )
    For NNPDF Hessian-converted sets, C = 1 (already 68% CL)."""
    central = vals[0]
    n_bins = vals.shape[1]
    delta_sq = np.zeros(n_bins)
    for i in range(1, n_dirs + 1):
        delta_sq += (vals[i] - central) ** 2
    delta = np.sqrt(delta_sq) / C
    return central, delta, delta


def compute_replica_band(vals, errs, n_replicas):
    """Return (central, dplus, dminus) per bin for a replica PDF (symmetric).
    vals shape: (1 + n_replicas, n_bins).  Use member 0 as central; PDF
    uncertainty is the naive std-dev across replicas at 68% CL.

    Caveat: for replica PDFs with wide per-event weight variation (e.g.
    NNPDF3.1 high-pT), individual replicas can carry large MC statistical
    noise. The reported band then includes PDF + MC noise together. The
    cleaner solution for such cases is to use the Hessian-converted variant
    (NNPDF31_nlo_as_0118_hessian) and the Hessian master formula, which
    avoids per-replica MC contamination entirely."""
    rep_vals = vals[1 : 1 + n_replicas]
    central = vals[0]
    delta = rep_vals.std(axis=0, ddof=1)
    return central, delta, delta


def make_tgraph(edges, central, dplus, dminus, name):
    """Build a TGraphAsymmErrors with x=bin centres, half-bin-width x errors."""
    n = len(central)
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    bin_half = 0.5 * (edges[1:] - edges[:-1])
    g = ROOT.TGraphAsymmErrors(n)
    g.SetName(name)
    for i in range(n):
        g.SetPoint(i, float(bin_centers[i]), float(central[i]))
        g.SetPointError(i, float(bin_half[i]), float(bin_half[i]),
                        float(dminus[i]), float(dplus[i]))
    return g


def make_central_hist(edges, central, name):
    """Return a TH1F with the central values (no errors set)."""
    n = len(central)
    h = ROOT.TH1F(name, name, n, array("d", edges))
    for i in range(n):
        h.SetBinContent(i + 1, float(central[i]))
        h.SetBinError(i + 1, 0.0)
    return h


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=OUT_FILE,
                        help="Output ROOT file path")
    args = parser.parse_args()

    out_root = ROOT.TFile(args.out, "RECREATE")

    print(f"{'tag':>8s}  {'kind':>9s}  {'n_members':>9s}  {'central(20-22)':>15s}  {'+/- (20-22)':>20s}")
    for tag, (kind, n_p, C, n_total) in PDF_INFO.items():
        path = ROOT_FILE_FMT.format(tag=tag)
        try:
            edges, vals, errs = read_members(path, n_total)
        except (FileNotFoundError, KeyError) as e:
            print(f"  SKIP {tag}: {e}")
            continue

        if kind == "hessian":
            central, dplus, dminus = compute_hessian_band(vals, n_p, C)
        elif kind == "symmhessian":
            central, dplus, dminus = compute_symmhessian_band(vals, n_p, C)
        elif kind == "replicas":
            central, dplus, dminus = compute_replica_band(vals, errs, n_p)
        else:
            raise ValueError(f"unknown kind {kind}")

        # Find the 20-22 GeV bin for the printout
        idx_2022 = None
        for i in range(len(edges) - 1):
            if abs(edges[i] - 20.0) < 0.01 and abs(edges[i + 1] - 22.0) < 0.01:
                idx_2022 = i
                break
        if idx_2022 is not None:
            c = central[idx_2022]
            dp = dplus[idx_2022]
            dm = dminus[idx_2022]
            band_str = f"+{100*dp/c:5.2f}% / -{100*dm/c:5.2f}%"
        else:
            c = float("nan")
            band_str = "n/a"

        print(f"  {tag:>8s}  {kind:>9s}  {n_total:>9d}  {c:15.4f}  {band_str:>20s}")

        h_cent = make_central_hist(edges, central, f"h_pdf_central_{tag}")
        g_band = make_tgraph(edges, central, dplus, dminus, f"g_pdf_band_{tag}")
        out_root.cd()
        h_cent.Write()
        g_band.Write()

    out_root.Close()
    print(f"\nWritten to {args.out}")


if __name__ == "__main__":
    main()
