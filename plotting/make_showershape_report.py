#!/usr/bin/env python3
"""Generate a LaTeX report with 9-panel shower shape figures for each config_showershape variant.

Usage:
    python3 make_showershape_report.py                            # generate .tex files only
    python3 make_showershape_report.py --compile                  # generate + compile PDF
    python3 make_showershape_report.py --suffix showershape_nom   # single variant
    python3 make_showershape_report.py --suffix showershape_nom --compile

Outputs:
    showershape_report/showershape_report.tex          (master document)
    showershape_report/{suffix}/showershape_{suffix}.tex  (per-variant fragment)
    showershape_report/{suffix}/figures/               (PDFs copied from figures/{suffix}/)
"""

import os
import glob
import difflib
import shutil
import argparse
import subprocess

import yaml

CONFIG_DIR  = "../efficiencytool"
FIGURES_SRC = "figures"           # where plot_showershapes_variations.C writes PDFs
OUTPUT_DIR  = "showershape_report"

# 9 shower shape variables shown in a 3x3 grid per figure
FIGURE_VARS = [
    "e11_to_e33", "e32_to_e35", "bdt",
    "weta_cogx",  "wphi_cogx",  "et1",
    "et2",        "et3",        "et4",
]

ALL_CUTS = ["cut0", "cut1", "cut2", "cut3"]
CUT_LABELS = {
    "cut0": "no nbkg cut",
    "cut1": "nbkg cut",
    "cut2": "tight photon ID",
    "cut3": "non-tight (sideband)",
}

ALL_PFXS  = ["dis"]
MIXED_PFXS = ["dis_mixed"]

# Human-readable labels per figure prefix, used in subsection headers and captions
PFX_LABELS = {
    "dis":       "Nominal",
    "dis_mixed": "Double-interaction mixed (0.776 \u00d7 single + 0.224 \u00d7 double)",
}


# ---------------------------------------------------------------------------
# Discovery helpers
# ---------------------------------------------------------------------------

NOM_SUFFIX = "showershape"   # config_showershape.yaml is the reference


def get_suffixes():
    """Return list of showershape config suffixes with showershape (nom) first."""
    pattern = os.path.join(CONFIG_DIR, "config_showershape*.yaml")
    configs = sorted(glob.glob(pattern))
    suffixes = []
    for path in configs:
        basename = os.path.basename(path)           # config_showershape_nom.yaml
        stem = basename[len("config_"):]            # showershape_nom.yaml
        stem = stem.removesuffix(".yaml")           # showershape_nom
        suffixes.append(stem)

    # Put the bare showershape (reference config) first
    if NOM_SUFFIX in suffixes:
        suffixes.remove(NOM_SUFFIX)
        suffixes = [NOM_SUFFIX] + sorted(suffixes)
    else:
        suffixes = sorted(suffixes)
    return suffixes


def get_config_diff(suffix):
    """Return lines to display for the config block in a section.

    For the reference config (showershape): return the full file contents.
    For variants: return the unified diff against config_showershape.yaml.
    Returns an empty list if files are missing.
    """
    nom_path = os.path.join(CONFIG_DIR, f"config_{NOM_SUFFIX}.yaml")
    var_path = os.path.join(CONFIG_DIR, f"config_{suffix}.yaml")

    if not os.path.exists(nom_path):
        return []

    with open(nom_path) as fh:
        nom_lines = fh.readlines()

    if suffix == NOM_SUFFIX:
        return nom_lines

    if not os.path.exists(var_path):
        return []

    with open(var_path) as fh:
        var_lines = fh.readlines()

    diff = list(difflib.unified_diff(
        nom_lines, var_lines,
        fromfile=f"config_{NOM_SUFFIX}.yaml",
        tofile=f"config_{suffix}.yaml",
        n=2,
    ))
    return diff


def get_pt_bins(config_path):
    """Load YAML config and return pT_bins list."""
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    return cfg["analysis"]["pT_bins"]


# ---------------------------------------------------------------------------
# Plot collection and copying
# ---------------------------------------------------------------------------

def collect_plots(suffix, pfxs=None, cuts=None, pt_indices=None, ieta=0):
    """Return list of (src_path, dest_basename) for all existing PDFs."""
    if pfxs is None:
        pfxs = ALL_PFXS
    if cuts is None:
        cuts = ALL_CUTS
    if pt_indices is None:
        pt_indices = [1, 2, 3]

    result = []
    seen = set()
    for pfx in pfxs:
        for cut in cuts:
            for ipt in pt_indices:
                for var in FIGURE_VARS:
                    fname = f"{pfx}_{var}_eta{ieta}_pt{ipt}_{cut}.pdf"
                    src = os.path.join(FIGURES_SRC, suffix, fname)
                    if os.path.exists(src) and src not in seen:
                        result.append((src, fname))
                        seen.add(src)
    return result


def copy_figures(src_dest_pairs, dest_dir):
    """Copy source PDFs into dest_dir, skipping already-current files."""
    os.makedirs(dest_dir, exist_ok=True)
    copied = 0
    for src, basename in src_dest_pairs:
        dst = os.path.join(dest_dir, basename)
        if not os.path.exists(dst) or os.path.getmtime(src) > os.path.getmtime(dst):
            shutil.copy2(src, dst)
            copied += 1
    skipped = len(src_dest_pairs) - copied
    print(f"  Copied {copied} figure(s) to {dest_dir}/  ({skipped} already up to date)")


# ---------------------------------------------------------------------------
# LaTeX helpers
# ---------------------------------------------------------------------------

def escape_latex(s):
    return s.replace("_", r"\_").replace("%", r"\%")


def latex_escape_alltt(s):
    """Escape characters special in alltt."""
    s = s.replace("\\", "\\textbackslash{}")
    s = s.replace("{",  "\\{")
    s = s.replace("}",  "\\}")
    s = s.replace("$",  "\\$")
    s = s.replace("%",  "\\%")
    s = s.replace("#",  "\\#")
    s = s.replace("_",  "\\_")
    s = s.replace("^",  "\\^{}")
    s = s.replace("&",  "\\&")
    s = s.replace("~",  "\\textasciitilde{}")
    return s


def write_diff_block(f, lines):
    """Render diff/config lines in alltt with +/- coloring."""
    f.write("\\begin{alltt}\\small\n")
    for raw in lines:
        line = raw.rstrip("\n")
        escaped = latex_escape_alltt(line)
        if line.startswith("+++") or line.startswith("---"):
            f.write(f"\\textcolor{{difffile}}{{{escaped}}}\n")
        elif line.startswith("+"):
            f.write(f"\\textcolor{{diffadd}}{{{escaped}}}\n")
        elif line.startswith("-"):
            f.write(f"\\textcolor{{diffdel}}{{{escaped}}}\n")
        elif line.startswith("@@"):
            f.write(f"\\textcolor{{diffhunk}}{{{escaped}}}\n")
        else:
            f.write(escaped + "\n")
    f.write("\\end{alltt}\n\n")


# ---------------------------------------------------------------------------
# Fragment writer
# ---------------------------------------------------------------------------

def write_variant_fragment(suffix, pt_bins, out_tex_path,
                           pfxs=None, cuts=None, pt_indices=None, ieta=0):
    """Write a per-variant .tex fragment (no \\documentclass).

    \\includegraphics paths are relative to the master's compilation directory
    (showershape_report/), so they use the form:
        {suffix}/figures/{fname}
    """
    if pfxs is None:
        pfxs = ALL_PFXS
    if cuts is None:
        cuts = ALL_CUTS
    if pt_indices is None:
        pt_indices = [1, 2, 3]

    # Only include pfxs that have at least one figure for this suffix — lets
    # --mixed be passed globally without emitting MISSING-box sections for
    # configs that don't have double-interaction blending figures.
    def _pfx_has_any_figure(pfx):
        for cut in cuts:
            for ipt in pt_indices:
                for var in FIGURE_VARS:
                    fname = f"{pfx}_{var}_eta{ieta}_pt{ipt}_{cut}.pdf"
                    if os.path.exists(os.path.join(FIGURES_SRC, suffix, fname)):
                        return True
        return False
    pfxs = [p for p in pfxs if _pfx_has_any_figure(p)]

    # Directory where this fragment's figures live
    figs_rel = f"{suffix}/figures"   # relative to showershape_report/

    with open(out_tex_path, "w") as f:
        # Config diff block
        config_lines = get_config_diff(suffix)
        if config_lines:
            if suffix == NOM_SUFFIX:
                f.write("\\subsection*{Nominal configuration}\n\n")
            else:
                f.write("\\subsection*{Changes from nominal}\n\n")
            write_diff_block(f, config_lines)
        elif suffix != NOM_SUFFIX:
            f.write("\\subsection*{Changes from nominal}\n\n"
                    "\\textit{(config file not found)}\n\n")

        for pfx in pfxs:
            pfx_label = PFX_LABELS.get(pfx, pfx)
            if len(pfxs) > 1:
                f.write(f"\\subsection{{{escape_latex(pfx_label)}}}\n\n")
            for icut in cuts:
                cut_label = CUT_LABELS.get(icut, icut)
                heading = "\\subsubsection" if len(pfxs) > 1 else "\\subsection"
                f.write(f"{heading}{{{escape_latex(cut_label)}}}\n\n")

                for ipt in pt_indices:
                    # Build figure block
                    pt_lo = pt_bins[ipt] if ipt < len(pt_bins) else "?"
                    pt_hi = pt_bins[ipt + 1] if (ipt + 1) < len(pt_bins) else "?"

                    # Annotation macros (self-documenting in source)
                    f.write(f"\\renewcommand{{\\ptval}}{{{pfx}_pt{ipt}}}\n")
                    f.write(f"\\renewcommand{{\\mycut}}{{{icut}}}\n")
                    f.write(f"\\renewcommand{{\\mypfx}}{{{pfx}}}\n\n")

                    f.write("\\begin{figure}[htbp]\n")
                    f.write("  \\centering\n")
                    for var in FIGURE_VARS:
                        fname = f"{pfx}_{var}_eta{ieta}_pt{ipt}_{icut}.pdf"
                        full_path = os.path.join(FIGURES_SRC, suffix, fname)
                        tex_path = f"{figs_rel}/{fname}"
                        if os.path.exists(full_path):
                            f.write(f"  \\includegraphics[width=0.32\\linewidth]{{{tex_path}}}\n")
                        else:
                            f.write(f"  \\fbox{{\\parbox{{0.30\\linewidth}}{{\\centering\\small MISSING: {escape_latex(var)}}}}}\n")
                    f.write(f"  \\caption{{Shower shape distributions ({escape_latex(pfx_label)}, Data vs MC), "
                            f"${pt_lo} < E_T < {pt_hi}$~GeV, {escape_latex(cut_label)}.}}\n")
                    label = f"fig:ss_{suffix}_pt{ipt}_{icut}_{pfx}"
                    f.write(f"  \\label{{{label}}}\n")
                    f.write("\\end{figure}\n\n")

                f.write("\\clearpage\n\n")

    print(f"  Written fragment: {out_tex_path}")


# ---------------------------------------------------------------------------
# Master document writer
# ---------------------------------------------------------------------------

def write_master_tex(suffixes, output_dir):
    """Write the master showershape_report.tex."""
    tex_path = os.path.join(output_dir, "showershape_report.tex")
    with open(tex_path, "w") as f:
        f.write(r"""\documentclass[11pt,letterpaper]{article}
\usepackage{graphicx}
\usepackage[margin=0.75in]{geometry}
\usepackage{hyperref}
\usepackage{fancyhdr}
\usepackage{float}
\usepackage{alltt}
\usepackage{xcolor}
\definecolor{diffadd}{rgb}{0.0, 0.50, 0.15}
\definecolor{diffdel}{rgb}{0.75, 0.10, 0.10}
\definecolor{diffhunk}{rgb}{0.10, 0.30, 0.70}
\definecolor{difffile}{rgb}{0.40, 0.40, 0.40}

% Annotation macros used inside fragments
\newcommand{\ptval}{pt1}
\newcommand{\mycut}{cut0}
\newcommand{\mypfx}{dis}

\pagestyle{fancy}
\fancyhf{}
\rhead{Shower Shape Variations}
\lhead{\nouppercase{\leftmark}}
\rfoot{\thepage}

\begin{document}

\title{Shower Shape Distributions: Configuration Variations}
\date{\today}
\maketitle
\tableofcontents
\newpage

""")
        for suffix in suffixes:
            f.write(f"\\section{{\\texttt{{{escape_latex(suffix)}}}}}\n\n")
            fragment = f"{suffix}/showershape_{suffix}.tex"
            f.write(f"\\input{{{fragment}}}\n")
            f.write("\\clearpage\n\n")

        f.write("\\end{document}\n")

    print(f"Written master: {tex_path}")
    return tex_path


# ---------------------------------------------------------------------------
# PDF compilation
# ---------------------------------------------------------------------------

def compile_pdf(tex_file):
    """Run pdflatex twice (for TOC) from the directory containing tex_file."""
    tex_dir  = os.path.dirname(os.path.abspath(tex_file))
    tex_name = os.path.basename(tex_file)
    for run in range(2):
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", tex_name],
            cwd=tex_dir, capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f"pdflatex failed on run {run + 1}")
            print(result.stdout[-3000:])
            print(result.stderr[-1000:])
            return False
    pdf = os.path.join(tex_dir, tex_name.removesuffix(".tex") + ".pdf")
    if os.path.exists(pdf):
        print(f"PDF produced: {pdf}")
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--compile", action="store_true",
                        help="Run pdflatex after generating the .tex files")
    parser.add_argument("--output-dir", default=OUTPUT_DIR,
                        help=f"Output directory (default: {OUTPUT_DIR})")
    parser.add_argument("--suffix", default=None,
                        help="Process only this suffix (e.g. showershape_nom)")
    parser.add_argument("--cuts", nargs="+", default=ALL_CUTS,
                        choices=ALL_CUTS, metavar="CUT",
                        help="Cut levels to include (default: all four)")
    parser.add_argument("--pfx", nargs="+", default=ALL_PFXS,
                        metavar="PFX",
                        help="Plot filename prefixes (default: dis)")
    parser.add_argument("--mixed", action="store_true",
                        help="Append 'dis_mixed' prefix to include double-interaction "
                             "mixed shower shape plots (requires plot_showershapes_variations.C "
                             "to be run with use_mixed=true first)")
    parser.add_argument("--pt-indices", nargs="+", type=int, default=[1, 2, 3],
                        metavar="N",
                        help="pT bin indices to include (default: 1 2 3)")
    args = parser.parse_args()

    # Discover suffixes
    all_suffixes = get_suffixes()
    if not all_suffixes:
        print(f"No config_showershape*.yaml files found in {CONFIG_DIR}/")
        print("Run from the plotting/ directory, or check that configs exist.")
        return

    if args.suffix:
        if args.suffix not in all_suffixes:
            print(f"Suffix '{args.suffix}' not found. Available:")
            for s in all_suffixes:
                print(f"  {s}")
            return
        suffixes = [args.suffix]
    else:
        suffixes = all_suffixes

    print(f"Processing {len(suffixes)} config variant(s):")
    for s in suffixes:
        print(f"  {s}")

    if args.mixed and "dis_mixed" not in args.pfx:
        args.pfx = args.pfx + ["dis_mixed"]
        print("--mixed: appending 'dis_mixed' prefix")

    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)

    for suffix in suffixes:
        print(f"\n--- {suffix} ---")
        config_path = os.path.join(CONFIG_DIR, f"config_{suffix}.yaml")
        if not os.path.exists(config_path):
            print(f"  [SKIP] Config file not found: {config_path}")
            continue

        pt_bins = get_pt_bins(config_path)

        variant_dir  = os.path.join(out_dir, suffix)
        dest_figures = os.path.join(variant_dir, "figures")
        os.makedirs(variant_dir, exist_ok=True)

        # Collect and copy figures
        src_dest = collect_plots(suffix, args.pfx, args.cuts, args.pt_indices)
        copy_figures(src_dest, dest_figures)

        # Write fragment
        fragment_path = os.path.join(variant_dir, f"showershape_{suffix}.tex")
        write_variant_fragment(
            suffix, pt_bins, fragment_path,
            pfxs=args.pfx, cuts=args.cuts, pt_indices=args.pt_indices,
        )

    # Write master always from all discovered suffixes, regardless of whether
    # --suffix was used to rebuild a subset. Previous behavior pinned the master
    # to the single --suffix, silently dropping the other 9 fragments from the PDF.
    master_tex = write_master_tex(all_suffixes, out_dir)

    if args.compile:
        print("\nCompiling PDF...")
        compile_pdf(master_tex)


if __name__ == "__main__":
    main()
