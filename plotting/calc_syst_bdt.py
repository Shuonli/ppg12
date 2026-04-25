#!/usr/bin/env python3
"""Automated BDT systematic uncertainty pipeline.

Reads systematic type metadata from make_bdt_variations.py, computes per-type
deviations from nominal, aggregates into groups and total, writes ROOT files,
and produces plots following the existing sPHENIX syst_*.C style.

Usage:
    python calc_syst_bdt.py [options]

Output ROOT files (rootFiles/):
    syst_bdt_{type}.root   — per syst_type  (h_dev_low/high, h_dev_rel_low/high)
    syst_bdt_{group}.root  — per group      (quadrature sum of member types)
    syst_bdt_total.root    — total + lumi   (quadrature sum of all FINAL_SYSTS)

Output figures (figures/):
    syst_bdt_rel_{type}.pdf   — per-type two-pad plot (spectrum + rel. deviation)
    syst_bdt_breakdown.pdf    — all groups + total on one canvas
"""

import argparse
import os
import sys

import ROOT

# ---------------------------------------------------------------------------
# Load sPHENIX style + helpers from plotcommon.h
# ---------------------------------------------------------------------------
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(_THIS_DIR, "plotcommon.h"))  # interpreted mode (no +)

# ---------------------------------------------------------------------------
# Import systematic metadata from make_bdt_variations.py
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(_THIS_DIR, "..", "efficiencytool"))
from make_bdt_variations import (  # noqa: E402
    VARIANTS, SYST_TYPES, SYST_GROUPS, FINAL_SYSTS, FLAT_SYSTS, LUMI_SYST,
)

# ---------------------------------------------------------------------------
# Colors for groups in the breakdown plot
# ---------------------------------------------------------------------------
GROUP_COLORS = {
    # Phase-2 (2026-04-25) restructure: 9 Resp groups + flat sources.
    # Colors chosen for distinguishability in the breakdown overlay
    # (avoid kGray/kBlack which collide with the total band).
    "photon_id":      ROOT.kAzure + 7,    # blue — BDT cuts + NPB
    "abcd_region":    ROOT.kCyan + 2,     # teal — noniso boundary
    "purity_method":  ROOT.kViolet + 1,   # violet — fit form + non-closure
    "iso_resolution": ROOT.kRed + 1,      # red   — mciso pedestal
    "acceptance":     ROOT.kPink + 7,     # pink  — phi-symm tower mask
    "unfolding":      ROOT.kMagenta + 2,  # magenta — reweight + iter scan
    "escale":         ROOT.kGreen - 2,    # green — energy scale
    "eres":           ROOT.kOrange + 1,   # orange — energy resolution
    "di_fraction":    ROOT.kSpring - 6,   # yellow-green — DI blending (Phase-3)
    # Legacy keys (kept for backward compat with old breakdown plots):
    "purity":         ROOT.kViolet + 1,
    "eff":            ROOT.kAzure + 7,
    "mbd":            ROOT.kViolet - 1,
}


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", default="/sphenix/user/shuhangli/ppg12/efficiencytool/results",
                   help="Directory containing Photon_final_*.root files")
    p.add_argument("--outdir",  default="rootFiles",
                   help="Output directory for ROOT files (default: rootFiles/)")
    p.add_argument("--figdir",  default="figures",
                   help="Output directory for plots (default: figures/)")
    p.add_argument("--nom",     default="bdt_nom",
                   help="var_type string for the nominal file (default: bdt_nom)")
    p.add_argument("--histogram", default="h_unfold_sub_result",
                   help="Histogram name inside each ROOT file (default: h_unfold_sub_result)")
    p.add_argument("--skip-missing", action="store_true",
                   help="Skip syst types whose result files are missing")
    return p.parse_args()


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
def load_spectrum(var_type: str, results_dir: str, histogram: str) -> ROOT.TH1F:
    """Open result file, clone histogram (detached from file), return TH1F."""
    fname = os.path.join(results_dir, f"Photon_final_{var_type}.root")
    f = ROOT.TFile.Open(fname, "READ")
    if not f or f.IsZombie():
        raise FileNotFoundError(f"Cannot open: {fname}")
    h = f.Get(histogram)
    if not h:
        f.Close()
        raise KeyError(f"Histogram '{histogram}' not found in {fname}")
    h = h.Clone(f"h_{var_type}")
    h.SetDirectory(0)
    f.Close()
    return h


def write_syst_root(outdir: str, tag: str,
                    h_al, h_ah, h_rl, h_rh) -> None:
    """Write four deviation histograms to rootFiles/syst_bdt_{tag}.root."""
    fname = os.path.join(outdir, f"syst_bdt_{tag}.root")
    fout = ROOT.TFile(fname, "RECREATE")
    for h, name in [(h_al, "h_dev_low"),  (h_ah, "h_dev_high"),
                    (h_rl, "h_dev_rel_low"), (h_rh, "h_dev_rel_high")]:
        hout = h.Clone(name)
        hout.SetDirectory(fout)
        hout.Write()
    fout.Close()


def write_syst_sum(outdir: str, h_al, h_ah, h_rl, h_rh) -> None:
    """Write total syst in the format expected by plot_final_selection.C.

    File: rootFiles/syst_sum.root
    Histograms: h_sum_low, h_sum_high, h_sum_rel_low, h_sum_rel_high
    """
    fname = os.path.join(outdir, "syst_sum.root")
    fout = ROOT.TFile(fname, "RECREATE")
    for h, name in [(h_al, "h_sum_low"),  (h_ah, "h_sum_high"),
                    (h_rl, "h_sum_rel_low"), (h_rh, "h_sum_rel_high")]:
        hout = h.Clone(name)
        hout.SetDirectory(fout)
        hout.Write()
    fout.Close()


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------
def calc_delta(h_var: ROOT.TH1F, h_nom: ROOT.TH1F, tag: str):
    """Return (h_dev_abs, h_dev_rel) — per-bin signed |var - nom|.

    No stat-subtract. The earlier stat-subtract path
    (sqrt(max(0, dev^2 - sigma_stat^2))) was OVER-aggressive: ROOT
    propagates sigma_dev = sqrt(sigma_var^2 + sigma_nom^2) which assumes
    uncorrelated samples, but PPG12 systematic variants share the same
    underlying MC events with the nominal (only the cut decision differs
    per event), so the true sigma_dev is much smaller than the propagated
    value. The over-subtracted shift was silently zeroing real systematic
    deviations, especially at high pT where MC stat is small.

    Standard treatment matches most published xsec analyses: take |dev|
    directly. MC-stat noise contamination is small (<<dev) for ABCD-
    propagated systematics and only matters for high-pT bins with very
    low statistics, where it is treated as a known conservative
    over-estimate rather than subtracted.
    """
    h_dev = h_var.Clone(f"h_dev_{tag}")
    h_dev.SetDirectory(0)
    h_dev.Add(h_nom, -1.0)
    h_dev_rel = h_dev.Clone(f"h_dev_rel_{tag}")
    h_dev_rel.SetDirectory(0)
    for i in range(1, h_nom.GetNbinsX() + 1):
        nom = h_nom.GetBinContent(i)
        if nom != 0:
            h_dev_rel.SetBinContent(i, h_dev.GetBinContent(i) / nom)
        else:
            h_dev_rel.SetBinContent(i, 0.0)
        h_dev_rel.SetBinError(i, 0.0)
    return h_dev, h_dev_rel


def _zero_clone(h_ref: ROOT.TH1F, name: str) -> ROOT.TH1F:
    h = h_ref.Clone(name)
    h.Reset()
    h.SetDirectory(0)
    return h


def build_variant_map() -> dict:
    """Return {syst_type: {role: [variant_name, ...]}} from VARIANTS metadata."""
    vmap = {t: {"up": [], "down": [], "one_sided": [], "max": []}
            for t in SYST_TYPES}
    for v in VARIANTS:
        st = v.get("syst_type")
        sr = v.get("syst_role")
        if st is None or sr is None or st not in vmap:
            continue
        vmap[st][sr].append(v["name"])
    return vmap


def aggregate_type(syst_type: str, vmap: dict, h_nom: ROOT.TH1F,
                   results_dir: str, histogram: str,
                   skip_missing: bool):
    """
    Compute (h_al, h_ah, h_rl, h_rh) for one syst_type.
    Returns None for placeholder types or missing files (when skip_missing=True).

    Mode logic:
      two_sided  : low=|delta(down)|, high=|delta(up)| per bin; mirrors if one direction absent
      one_sided  : low=high=|delta(var)| (symmetric)
      max        : low=high=bin-wise max of |delta| over all "max" variants
      placeholder: returns None with a warning
    """
    mode = SYST_TYPES[syst_type]["mode"]
    if mode == "placeholder":
        print(f"  [PLACEHOLDER] {syst_type} — skipped (no config variant defined yet)")
        return None

    roles = vmap[syst_type]
    nbins = h_nom.GetNbinsX()
    h_al = _zero_clone(h_nom, f"h_al_{syst_type}")
    h_ah = _zero_clone(h_nom, f"h_ah_{syst_type}")
    h_rl = _zero_clone(h_nom, f"h_rl_{syst_type}")
    h_rh = _zero_clone(h_nom, f"h_rh_{syst_type}")

    def load_var(name):
        var_type = f"bdt_{name}"
        try:
            return load_spectrum(var_type, results_dir, histogram)
        except (FileNotFoundError, KeyError, OSError) as e:
            if skip_missing:
                print(f"  [SKIP] {var_type}: {e}")
                return None
            raise

    def fill_abs(h_abs, h_rel, h_var):
        hd, hr = calc_delta(h_var, h_nom, f"{syst_type}_{h_var.GetName()}")
        for i in range(1, nbins + 1):
            h_abs.SetBinContent(i, abs(hd.GetBinContent(i)))
            h_rel.SetBinContent(i, abs(hr.GetBinContent(i)))

    if mode == "two_sided":
        up_ok = False
        if roles["up"]:
            hv = load_var(roles["up"][0])
            if hv is not None:
                fill_abs(h_ah, h_rh, hv)
                up_ok = True
        down_ok = False
        if roles["down"]:
            hv = load_var(roles["down"][0])
            if hv is not None:
                fill_abs(h_al, h_rl, hv)
                down_ok = True
        # Mirror if only one direction was available
        if not up_ok:
            h_ah = h_al.Clone(f"h_ah_{syst_type}"); h_ah.SetDirectory(0)
            h_rh = h_rl.Clone(f"h_rh_{syst_type}"); h_rh.SetDirectory(0)
        if not down_ok:
            h_al = h_ah.Clone(f"h_al_{syst_type}"); h_al.SetDirectory(0)
            h_rl = h_rh.Clone(f"h_rl_{syst_type}"); h_rl.SetDirectory(0)

    elif mode == "one_sided":
        if not roles["one_sided"]:
            return None
        hv = load_var(roles["one_sided"][0])
        if hv is None:
            return None
        hd, hr = calc_delta(hv, h_nom, syst_type)
        for i in range(1, nbins + 1):
            v  = abs(hd.GetBinContent(i))
            vr = abs(hr.GetBinContent(i))
            h_al.SetBinContent(i, v);  h_ah.SetBinContent(i, v)
            h_rl.SetBinContent(i, vr); h_rh.SetBinContent(i, vr)

    elif mode == "max":
        for vname in roles["max"]:
            hv = load_var(vname)
            if hv is None:
                continue
            hd, hr = calc_delta(hv, h_nom, f"{syst_type}_{vname}")
            for i in range(1, nbins + 1):
                v  = abs(hd.GetBinContent(i))
                vr = abs(hr.GetBinContent(i))
                if v  > h_al.GetBinContent(i): h_al.SetBinContent(i, v)
                if vr > h_rl.GetBinContent(i): h_rl.SetBinContent(i, vr)
        # For max mode, low == high
        h_ah = h_al.Clone(f"h_ah_{syst_type}"); h_ah.SetDirectory(0)
        h_rh = h_rl.Clone(f"h_rh_{syst_type}"); h_rh.SetDirectory(0)

    return h_al, h_ah, h_rl, h_rh


def quadrature_sum(components: list) -> tuple:
    """Add a list of (h_al, h_ah, h_rl, h_rh) tuples in quadrature bin-by-bin."""
    ref = components[0][0]
    nbins = ref.GetNbinsX()
    h_al = _zero_clone(ref, "h_qs_al")
    h_ah = _zero_clone(ref, "h_qs_ah")
    h_rl = _zero_clone(ref, "h_qs_rl")
    h_rh = _zero_clone(ref, "h_qs_rh")
    for i in range(1, nbins + 1):
        al2 = sum(c[0].GetBinContent(i) ** 2 for c in components)
        ah2 = sum(c[1].GetBinContent(i) ** 2 for c in components)
        rl2 = sum(c[2].GetBinContent(i) ** 2 for c in components)
        rh2 = sum(c[3].GetBinContent(i) ** 2 for c in components)
        h_al.SetBinContent(i, al2 ** 0.5)
        h_ah.SetBinContent(i, ah2 ** 0.5)
        h_rl.SetBinContent(i, rl2 ** 0.5)
        h_rh.SetBinContent(i, rh2 ** 0.5)
    return h_al, h_ah, h_rl, h_rh


def add_flat(h_al, h_ah, h_rl, h_rh, h_nom, flat_systs=None) -> tuple:
    """Add a dict of flat fractional uncertainties in quadrature.

    flat_systs: {name: {"down": frac, "up": frac}}. Default = FLAT_SYSTS
    from make_bdt_variations (lumi + l1_plateau). Each source contributes
    a multiplicative-flat uncertainty applied post-unfold; this is the
    correct treatment for sources independent of bin migration (lumi from
    MBD xsec, L1 photon plateau, MBD trigger eff if added).
    """
    if flat_systs is None:
        flat_systs = FLAT_SYSTS
    nbins = h_nom.GetNbinsX()
    for src_name, sym in flat_systs.items():
        f_lo = sym.get("down", 0.0)
        f_hi = sym.get("up",   0.0)
        if f_lo == 0.0 and f_hi == 0.0:
            continue
        for i in range(1, nbins + 1):
            nom = h_nom.GetBinContent(i)
            h_al.SetBinContent(i, (h_al.GetBinContent(i) ** 2 + (f_lo * nom) ** 2) ** 0.5)
            h_ah.SetBinContent(i, (h_ah.GetBinContent(i) ** 2 + (f_hi * nom) ** 2) ** 0.5)
            h_rl.SetBinContent(i, (h_rl.GetBinContent(i) ** 2 + f_lo ** 2) ** 0.5)
            h_rh.SetBinContent(i, (h_rh.GetBinContent(i) ** 2 + f_hi ** 2) ** 0.5)
    return h_al, h_ah, h_rl, h_rh


# Backward-compat alias — will be removed once all callers migrate to add_flat.
def add_lumi(h_al, h_ah, h_rl, h_rh, h_nom) -> tuple:
    """Deprecated: use add_flat(...) instead. Adds only the lumi contribution."""
    return add_flat(h_al, h_ah, h_rl, h_rh, h_nom, {"lumi": LUMI_SYST})


# ---------------------------------------------------------------------------
# Plotting helpers (PyROOT, following syst_*.C / calcSyst.C style)
# ---------------------------------------------------------------------------
def _yrange_from_hists(hists, padding=0.3):
    """
    Return (ylo, yhi) for a ratio plot that comfortably fits all histograms.
    padding is the fractional margin added above the max and below the min.
    Always includes zero; minimum full-range is ±0.02.
    """
    vmax = 0.0
    vmin = 0.0
    for h in hists:
        for i in range(1, h.GetNbinsX() + 1):
            v = h.GetBinContent(i)
            if v > vmax:
                vmax = v
            if v < vmin:
                vmin = v
    span = max(vmax - vmin, 0.04)
    margin = span * padding
    return vmin - margin, vmax + margin


def plot_syst_type(figdir: str, type_name: str, result: tuple,
                   h_nom: ROOT.TH1F, results_dir: str, histogram: str,
                   vmap: dict, args) -> None:
    """
    Single-pad ratio-only canvas per syst_type.
    Y axis range is auto-scaled to the actual relative deviation values.
    """
    h_al, h_ah, h_rl, h_rh = result

    ROOT.init_plot()

    c = ROOT.TCanvas(f"c_{type_name}", "", 800, 600)
    c.cd()
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.08)
    c.SetTopMargin(0.07)
    c.SetBottomMargin(0.13)

    # Build the two plotted histograms: -h_rl (down, blue) and h_rh (up, red)
    h_rl_plot = h_rl.Clone(f"h_rl_plot_{type_name}")
    h_rl_plot.SetDirectory(0)
    h_rl_plot.Scale(-1)

    # Auto-scale Y axis to fit both bands with padding
    ylo, yhi = _yrange_from_hists([h_rl_plot, h_rh])

    ROOT.frame_et_truth.SetYTitle("Relative difference")
    ROOT.frame_et_truth.GetYaxis().SetNdivisions(506)
    ROOT.frame_et_truth.GetYaxis().SetRangeUser(ylo, yhi)
    ROOT.frame_et_truth.GetXaxis().SetRangeUser(10, 32)
    ROOT.frame_et_truth.GetXaxis().SetNdivisions(505)
    ROOT.frame_et_truth.SetXTitle("#it{E}_{T}^{#gamma} [GeV]")
    ROOT.frame_et_truth.Draw("axis")

    h_rl_plot.SetMarkerStyle(20)
    h_rl_plot.SetMarkerColor(ROOT.kBlue)
    h_rl_plot.SetLineColor(ROOT.kBlue)
    h_rl_plot.SetLineWidth(2)
    h_rl_plot.Draw("same ][ HIST")

    h_rh.SetMarkerStyle(20)
    h_rh.SetMarkerColor(ROOT.kRed)
    h_rh.SetLineColor(ROOT.kRed)
    h_rh.SetLineWidth(2)
    h_rh.Draw("same ][ HIST")

    ROOT.linezero.Draw("L")

    ROOT.myText(0.55, 0.88, 1, ROOT.strleg1.c_str(), 0.045)
    ROOT.myText(0.55, 0.83, 1, ROOT.strleg2.c_str(), 0.045)
    ROOT.myMarkerLineText(0.18, 0.88, 1, ROOT.kRed,  20, ROOT.kRed,  1, "up var.",   0.05, True)
    ROOT.myMarkerLineText(0.18, 0.82, 1, ROOT.kBlue, 20, ROOT.kBlue, 1, "down var.", 0.05, True)
    ROOT.myText(0.18, 0.76, 1, type_name.replace("_", " "), 0.05)

    out = os.path.join(figdir, f"syst_bdt_rel_{type_name}.pdf")
    c.SaveAs(out)
    c.Close()


def plot_breakdown(figdir: str, group_results: dict,
                   total: tuple, args) -> None:
    """
    Single-canvas breakdown plot showing all group bands + total.
    Mirrors the style of calcSyst.C.
    """
    ROOT.init_plot()

    c = ROOT.TCanvas("c_breakdown", "", 900, 600)

    ROOT.frame_et_truth.GetYaxis().SetRangeUser(-0.58, 0.5)
    ROOT.frame_et_truth.GetXaxis().SetRangeUser(ROOT.pTmin, ROOT.pTmax)
    ROOT.frame_et_truth.SetXTitle("#it{E}_{T}^{#gamma} [GeV]")
    ROOT.frame_et_truth.SetYTitle("Relative difference")
    ROOT.frame_et_truth.Draw("axis")

    # Keep explicit Python references to all drawn clones so PyROOT does not
    # garbage-collect one side of the bands before SaveAs().
    keep_alive = []

    # Draw total (black, thick)
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

    switchover = 3
    legend_entries = list(group_results.items())   # ordered as in SYST_GROUPS
    # Draw escale last so its lower band does not get hidden by overlays.
    draw_order = [e for e in legend_entries if e[0] != "escale"] + \
                 [e for e in legend_entries if e[0] == "escale"]
    for isys, (grp, (h_al, h_ah, h_rl, h_rh)) in enumerate(draw_order):
        col = GROUP_COLORS.get(grp, ROOT.kGray)

        h_rl_plot = h_rl.Clone(f"h_rl_plot_{grp}")
        h_rl_plot.SetDirectory(0)
        h_rl_plot.Scale(-1)
        keep_alive.append(h_rl_plot)

        # Keep both sides solid for consistency in the breakdown panel.
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

        xshift = 0.33 if isys >= switchover else 0.0
        yshift = 0.05 * (switchover + 1) if isys >= switchover else 0.0
        ROOT.myMarkerLineText(0.25 + xshift, 0.31 - 0.05 * isys + yshift,
                              0, col, 20, col, 0,
                              grp.replace("_", " "), 0.05, True)

    ROOT.myText(0.20, 0.88, 1, ROOT.strleg1.c_str(), 0.05)
    ROOT.myText(0.20, 0.83, 1, ROOT.strleg2.c_str(), 0.05)

    out = os.path.join(figdir, "syst_bdt_breakdown.pdf")
    c.SaveAs(out)
    c.Close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    ROOT.init_plot()
    print(f"Loading nominal: {args.nom}")
    h_nom = load_spectrum(args.nom, args.results, args.histogram)

    vmap = build_variant_map()

    # ------------------------------------------------------------------
    # 1. Per-type systematics
    # ------------------------------------------------------------------
    type_results = {}
    print("\n--- Per-type systematics ---")
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
        write_syst_root(args.outdir, type_name, *result)
        plot_syst_type(args.figdir, type_name, result, h_nom,
                       args.results, args.histogram, vmap, args)
        print(f"  [OK]  {type_name}")

    # ------------------------------------------------------------------
    # 2. Group-level quadrature sums
    # ------------------------------------------------------------------
    group_results = {}
    print("\n--- Group quadrature sums ---")
    for grp, members in SYST_GROUPS.items():
        comps = [type_results[t] for t in members if t in type_results]
        if not comps:
            print(f"  [SKIP] group:{grp} — no member types computed")
            continue
        group_results[grp] = quadrature_sum(comps)
        write_syst_root(args.outdir, grp, *group_results[grp])
        # Per-group syst plot (the analysis note references these by group
        # rather than by type — e.g. "Photon Identification" combines tight
        # + nontight + npb_cut into a single photon_id band).
        plot_syst_type(args.figdir, f"group_{grp}", group_results[grp], h_nom,
                       args.results, args.histogram, vmap, args)
        print(f"  [OK]  group:{grp}  (members: {[t for t in members if t in type_results]})")

    # ------------------------------------------------------------------
    # 3. Total (quadrature sum of FINAL_SYSTS + lumi)
    # ------------------------------------------------------------------
    print("\n--- Total systematic ---")
    total_comps = [group_results[g] for g in FINAL_SYSTS if g in group_results]
    if not total_comps:
        print("  [ERROR] No group results available — cannot compute total.")
        return
    total = add_flat(*quadrature_sum(total_comps), h_nom, FLAT_SYSTS)
    write_syst_root(args.outdir, "total", *total)
    write_syst_sum(args.outdir, *total)
    included = [g for g in FINAL_SYSTS if g in group_results]
    flat_names = list(FLAT_SYSTS.keys())
    print(f"  [OK]  total (groups: {included} + flat: {flat_names})")
    print(f"  [OK]  syst_sum.root  (h_sum_low/high, h_sum_rel_low/high — for plot_final_selection.C)")

    # ------------------------------------------------------------------
    # 4. Summary breakdown plot
    # ------------------------------------------------------------------
    print("\n--- Generating breakdown plot ---")
    plot_breakdown(args.figdir, group_results, total, args)
    print(f"  [OK]  {os.path.join(args.figdir, 'syst_bdt_breakdown.pdf')}")

    print(f"\nDone.")
    print(f"  ROOT files : {args.outdir}/")
    print(f"  Figures    : {args.figdir}/")


if __name__ == "__main__":
    main()
