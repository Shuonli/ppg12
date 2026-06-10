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
    # Note-aligned grouping (2026-04-28): purity merges tight/nontight/noniso/
    # fit/mc-closure; efficiency = iso_resolution; npb is its own group.
    "escale":      ROOT.kGreen - 2,    # green — energy scale
    "eres":        ROOT.kOrange + 1,   # orange — energy resolution
    "purity":      ROOT.kViolet + 1,   # violet — tight+nontight+noniso+fit+mc-closure
    "efficiency":  ROOT.kRed + 1,      # red   — iso pedestal
    "unfolding":   ROOT.kMagenta + 2,  # magenta — reweight + iter scan
    "di_fraction": ROOT.kCyan + 2,     # teal-cyan — DI blending (distinct from escale green)
    "npb":         ROOT.kAzure + 7,    # blue — NPB cut
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
    from make_bdt_variations (currently {lumi}). Each source contributes
    a multiplicative-flat uncertainty applied post-unfold; this is the
    correct treatment for sources independent of bin migration (lumi from
    MBD xsec). The L1-plateau term was retired on 2026-04-28 once the
    bit-30 turn-on is corrected per-bin in CalculatePhotonYield.
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
def _yrange_from_hists(hists, padding=0.3, x_min=12.0, x_max=32.0):
    """
    Return (ylo, yhi) for a ratio plot that comfortably fits all histograms,
    looking only at bins whose center falls in (x_min, x_max] — i.e. the
    visible analysis range. Bins outside this window (truth overflow/under-
    flow [8,10] / [32,45], etc.) are ignored so a single huge edge bin can't
    blow up the y axis.

    padding is the fractional margin added above the max and below the min.
    Always includes zero; minimum full-range is ~0.04.
    """
    vmax = 0.0
    vmin = 0.0
    for h in hists:
        for i in range(1, h.GetNbinsX() + 1):
            xc = h.GetBinCenter(i)
            if xc < x_min or xc > x_max:
                continue
            v = h.GetBinContent(i)
            if v > vmax:
                vmax = v
            if v < vmin:
                vmin = v
    span = max(vmax - vmin, 0.04)
    margin = span * padding
    return vmin - margin, vmax + margin


# ---------------------------------------------------------------------------
# Physical labels for variants — used in the per-syst-type plots so the
# legend shows what was actually varied (e.g. "ET scale +2.6%") instead of
# the generic "up var." / "down var.". Variants not in this dict fall back
# to their config name with underscores replaced by spaces.
# ---------------------------------------------------------------------------
VARIANT_LABELS = {
    # Energy scale / resolution
    "energyscale148up":   "#it{E}_{T} scale #times 1.0148",
    "energyscale148down": "#it{E}_{T} scale #times 0.9852",
    "energyscale26up":   "#it{E}_{T} scale #times 1.026",
    "energyscale26down": "#it{E}_{T} scale #times 0.974",
    "energyresolution0": "no extra smearing (#sigma_{E}/E = 0%)",
    "energyresolution8": "extra smearing #sigma_{E}/E = 8%",
    # ABCD region
    "noniso04": "non-iso gap +0.4 GeV",
    "noniso10": "non-iso gap +1.0 GeV",
    # Photon ID
    "npb03":      "NPB cut = 0.3 (looser)",
    "npb07":      "NPB cut = 0.7 (tighter)",
    "tightup_p05": "tight-BDT intercept +0.05",
    "ntdown_m10":  "non-tight BDT intercept -0.10",
    # Purity
    "purity_pade":          "purity erf fit",
    "mc_purity_correction": "MC closure correction",
    # Iso resolution
    "mciso_no_shift": "MC iso pedestal shift = 0",
    # Unfolding
    "no_unfolding_reweighting": "unweighted response matrix",
    "unfold_iter1": "Bayes iter = 1",
    "unfold_iter3": "Bayes iter = 3",
    "unfold_iter4": "Bayes iter = 4",
    # Acceptance
    "mask_phisymm_or": "#phi-symm tower-mask OR",
    # Double interaction
    "di_frac_p20_0rad":     "f_{DI}(0 mrad) #times 1.20",
    "di_frac_m20_0rad":     "f_{DI}(0 mrad) #times 0.80",
    "di_frac_p20_1p5mrad":  "f_{DI}(1.5 mrad) #times 1.20",
    "di_frac_m20_1p5mrad":  "f_{DI}(1.5 mrad) #times 0.80",
    "di_frac_fit":          "f_{DI} from data BDT-score #chi^{2} fit",
}


def _variant_label(name: str) -> str:
    """Look up a physical legend label for a variant; fall back to the cleaned
    config name."""
    return VARIANT_LABELS.get(name, name.replace("_", " "))


def _draw_line_legend(x, y, color, text, tsize=0.045, dashed=False):
    """Draw a short coloured line + label at NDC (x, y)."""
    ln = ROOT.TLine()
    ln.SetLineColor(color)
    ln.SetLineWidth(2)
    ln.SetLineStyle(2 if dashed else 1)
    ln.DrawLineNDC(x - 0.045, y, x - 0.010, y)
    lab = ROOT.TLatex()
    lab.SetNDC()
    lab.SetTextAlign(12)
    lab.SetTextSize(tsize)
    lab.DrawLatex(x, y, text)


def plot_syst_type(figdir: str, type_name: str, result: tuple,
                   h_nom: ROOT.TH1F, results_dir: str, histogram: str,
                   vmap: dict, args) -> None:
    """Per-syst-type ratio-only canvas.

    Rendering depends on the syst_type's mode:
      two_sided  : two SIGNED relative-deviation lines, one per variant, with
                   physical legend labels (e.g. "ET scale x 1.026").
      one_sided  : a single SIGNED line for the lone variant. Caption notes
                   that the budget symmetrizes |delta|.
      max        : one SIGNED line per "max" variant. Caption notes that the
                   budget takes the per-bin envelope of |delta|.
      group_*    : (called via type_name="group_<grp>") -- legacy mirror
                   envelope is plotted because there is no single variant.

    The deviation is *signed* (var - nom)/nom and is NOT flipped to ensure
    the band stays on one side of zero — letting it cross zero is the honest
    visualization that motivates the symmetrization step.
    """
    h_al, h_ah, h_rl, h_rh = result

    ROOT.init_plot()

    c = ROOT.TCanvas(f"c_{type_name}", "", 800, 600)
    c.cd()
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.08)
    c.SetTopMargin(0.07)
    c.SetBottomMargin(0.13)

    is_group = type_name.startswith("group_") or type_name == "total"
    mode = SYST_TYPES.get(type_name, {}).get("mode", "")
    roles = vmap.get(type_name, {}) if not is_group else {}

    # Load signed (var - nom)/nom for each variant we want to draw.
    def signed_rel(vname):
        var_type = f"bdt_{vname}"
        try:
            hv = load_spectrum(var_type, results_dir, histogram)
        except (FileNotFoundError, KeyError, OSError):
            return None
        _, hr = calc_delta(hv, h_nom, f"plot_{type_name}_{vname}")
        return hr

    # Decide what to draw and gather (h, color, label) entries.
    entries = []  # (h, color, label)
    palette = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2,
               ROOT.kViolet + 1, ROOT.kOrange + 7, ROOT.kCyan + 2]
    caption_note = ""

    if is_group:
        # Fall back to the legacy envelope plot (mirrors |delta|): h_rh up
        # in red, -h_rl down in blue. Useful for group_*/total summaries.
        h_rl_plot = h_rl.Clone(f"h_rl_plot_{type_name}")
        h_rl_plot.SetDirectory(0)
        h_rl_plot.Scale(-1)
        entries.append((h_rh,        ROOT.kRed,  "envelope (up)"))
        entries.append((h_rl_plot,   ROOT.kBlue, "envelope (down)"))
    elif mode == "two_sided":
        for i, role in enumerate(("up", "down")):
            for vname in roles.get(role, []):
                hr = signed_rel(vname)
                if hr is None: continue
                entries.append((hr, palette[i], _variant_label(vname)))
    elif mode == "one_sided":
        for vname in roles.get("one_sided", []):
            hr = signed_rel(vname)
            if hr is None: continue
            entries.append((hr, palette[0], _variant_label(vname)))
        caption_note = "Single-sided variation; |#it{#delta}| symmetrized in the systematic budget."
    elif mode == "max":
        for i, vname in enumerate(roles.get("max", [])):
            hr = signed_rel(vname)
            if hr is None: continue
            entries.append((hr, palette[i % len(palette)], _variant_label(vname)))
        caption_note = "Per-bin |#it{#delta}|-envelope of these variations is taken in the budget."

    if not entries:
        # Nothing to draw — fall back to legacy behaviour so the figure still
        # appears (mirrored abs envelope).
        h_rl_plot = h_rl.Clone(f"h_rl_plot_{type_name}")
        h_rl_plot.SetDirectory(0)
        h_rl_plot.Scale(-1)
        entries = [(h_rh,      ROOT.kRed,  "up envelope"),
                   (h_rl_plot, ROOT.kBlue, "down envelope")]

    # Auto-scale Y axis to fit all drawn lines with padding.
    ylo, yhi = _yrange_from_hists([h for h, _, _ in entries])

    ROOT.frame_et_truth.SetYTitle("Relative difference")
    ROOT.frame_et_truth.GetYaxis().SetNdivisions(506)
    ROOT.frame_et_truth.GetYaxis().SetRangeUser(ylo, yhi)
    ROOT.frame_et_truth.GetXaxis().SetRangeUser(12, 32)
    ROOT.frame_et_truth.GetXaxis().SetNdivisions(505)
    ROOT.frame_et_truth.SetXTitle("#it{E}_{T}^{#gamma} [GeV]")
    ROOT.frame_et_truth.Draw("axis")

    for h, color, _ in entries:
        h.SetMarkerStyle(20)
        h.SetMarkerColor(color)
        h.SetLineColor(color)
        h.SetLineWidth(2)
        h.Draw("same ][ HIST")

    ROOT.linezero.Draw("L")

    ROOT.myText(0.55, 0.88, 1, ROOT.strleg1.c_str(), 0.045)
    ROOT.myText(0.55, 0.83, 1, ROOT.strleg2.c_str(), 0.045)

    # Per-variant legend entries (top-left).
    leg_x = 0.22
    leg_y = 0.88
    leg_dy = 0.06
    for i, (_, color, label) in enumerate(entries):
        _draw_line_legend(leg_x, leg_y - i * leg_dy, color, label, tsize=0.04)

    # Type label below the legend block.
    type_y = leg_y - len(entries) * leg_dy - 0.02
    ROOT.myText(0.18, type_y, 1, type_name.replace("_", " "), 0.045)

    # Caption note at the bottom, if any.
    if caption_note:
        cap = ROOT.TLatex()
        cap.SetNDC()
        cap.SetTextAlign(13)
        cap.SetTextSize(0.035)
        cap.DrawLatex(0.15, 0.20, caption_note)

    out = os.path.join(figdir, f"syst_bdt_rel_{type_name}.pdf")
    c.SaveAs(out)
    c.Close()


def _build_stat_rel(args) -> "ROOT.TH1F | None":
    """Read the nominal final-yield histogram and build a per-bin
    relative statistical uncertainty histogram, on the same binning as
    the syst breakdown.

    Returns None if the nominal file or histogram cannot be opened.
    """
    nom_path = os.path.join(args.results,
                            f"Photon_final_{args.nom}.root")
    f = ROOT.TFile.Open(nom_path, "READ")
    if not f or f.IsZombie():
        print(f"[stat overlay] Cannot open {nom_path}; stat curve skipped.")
        return None
    h_nom = f.Get(args.histogram)
    if not h_nom:
        print(f"[stat overlay] Missing {args.histogram} in {nom_path}; "
              f"stat curve skipped.")
        f.Close()
        return None
    h_stat = h_nom.Clone("h_stat_rel")
    h_stat.SetDirectory(0)
    h_stat.Reset()
    for i in range(1, h_nom.GetNbinsX() + 1):
        y = h_nom.GetBinContent(i)
        e = h_nom.GetBinError(i)
        h_stat.SetBinContent(i, (e / y) if y > 0 else 0.0)
        h_stat.SetBinError(i, 0.0)
    f.Close()
    return h_stat


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

    # Statistical uncertainty curve (per-bin relative). Drawn first so the
    # syst bands sit on top.
    h_stat_rel = _build_stat_rel(args)
    h_stat_rel_neg = None
    if h_stat_rel is not None:
        h_stat_rel.SetLineColor(ROOT.kGray + 2)
        h_stat_rel.SetLineStyle(2)
        h_stat_rel.SetLineWidth(2)
        h_stat_rel.SetMarkerSize(0)
        keep_alive.append(h_stat_rel)
        h_stat_rel_neg = h_stat_rel.Clone("h_stat_rel_neg")
        h_stat_rel_neg.SetDirectory(0)
        h_stat_rel_neg.Scale(-1.0)
        keep_alive.append(h_stat_rel_neg)
        h_stat_rel.Draw("same ][ HIST")
        h_stat_rel_neg.Draw("same ][ HIST")
        ROOT.myMarkerLineText(0.25, 0.41, 0, ROOT.kGray + 2, 0,
                              ROOT.kGray + 2, 2,
                              "Stat (#pm#sigma)", 0.05, True)

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
