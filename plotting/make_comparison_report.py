#!/usr/bin/env python3
"""Generate a LaTeX report comparing pairs of BDT config variants side-by-side.

For each pair (A, B), every plot type is shown with A on the left and B on the
right, making it easy to see the effect of a config change across all plots.

Usage:
    python3 make_comparison_report.py \\
        --pair bdt_nom bdt_tightbdt50 \\
        --pair bdt_nom bdt_npb03 \\
        [--compile] [--output-dir comparison_report]

Outputs:
    comparison_report/comparison_report.tex
    comparison_report/figures/<all referenced PDFs copied here>

Run after make_all_selection.sh has produced plots in figures/.
"""

import os
import glob
import re
import shutil
import difflib
import argparse
import subprocess

FIGURES_DIR = "figures"
CONFIG_DIR  = "../efficiencytool"
OUTPUT_DIR  = "selection_report"
TEX_NAME    = "comparison_report.tex"

# Number of pT bins (matches NptBins in plotcommon.h)
N_PT_BINS = 12

# Ordered plot groups with display label and expected filenames (use {s} for suffix).
PLOT_GROUPS = [
    ("Sideband (Data)", [
        "et_sbs_{s}.pdf",
        "et_sbs_ratio_{s}.pdf",
        "leakage_fraction_et_{s}.pdf",
    ]),
    ("Sideband (MC)", [
        "et_sbs_sim_{s}.pdf",
        "et_sbs_ratio_sim_{s}.pdf",
        "leakage_fraction_et_sim_{s}.pdf",
    ]),
    ("Purity", [
        "purity_{s}.pdf",
        "purity_sim_{s}.pdf",
    ]),
    ("Efficiency", [
        "eff_reco_{s}.pdf",
        "eff_iso_{s}.pdf",
        "eff_id_{s}.pdf",
        "eff_total_{s}.pdf",
        "eff_mbd_{s}.pdf",
        "eff_photon_{s}.pdf",
    ]),
    ("Final result", [
        "final_{s}.pdf",
        "final_phenix_{s}.pdf",
        "final_common_cluster_{s}.pdf",
        "final_tight_iso_cluster_{s}.pdf",
        "final_all_{s}.pdf",
    ]),
    ("Isolation fit (CONF)", [
        "h1D_iso_{s}_0_2.pdf",
        "h1D_iso_fit_{s}_0_2.pdf",
        "h1D_iso_{s}_3_5.pdf",
        "h1D_iso_fit_{s}_3_5.pdf",
        "h1D_iso_{s}_6_8.pdf",
        "h1D_iso_fit_{s}_6_8.pdf",
        "h1D_iso_{s}_9_9.pdf",
        "h1D_iso_fit_{s}_9_9.pdf",
    ]),
]


# ---------------------------------------------------------------------------
# Helpers (shared logic with make_selection_report.py)
# ---------------------------------------------------------------------------

def escape_latex(s):
    """Escape underscores and other special characters for LaTeX text mode."""
    return s.replace("_", r"\_")


def latex_escape_alltt(s):
    """Escape characters special in alltt: \\, {, } and other LaTeX metacharacters."""
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


def make_caption(fpath, suffix):
    """Generate a descriptive caption from the filename."""
    name = os.path.basename(fpath).removesuffix(".pdf")
    if name.endswith("_" + suffix):
        plot_type = name[: -(len(suffix) + 1)]
    else:
        plot_type = name
    labels = {
        "et_sbs":                  "Sideband ET spectrum (data)",
        "et_sbs_ratio":            "Sideband ET ratio (data)",
        "leakage_fraction_et":     "Leakage fraction vs ET (data)",
        "et_sbs_sim":              "Sideband ET spectrum (MC)",
        "et_sbs_ratio_sim":        "Sideband ET ratio (MC)",
        "leakage_fraction_et_sim": "Leakage fraction vs ET (MC)",
        "purity":                  "Purity (data)",
        "purity_sim":              "Purity (MC)",
        "eff_reco":                "Reconstruction efficiency",
        "eff_iso":                 "Isolation efficiency",
        "eff_id":                  "Identification efficiency",
        "eff_total":               "Total efficiency",
        "eff_mbd":                 "MBD trigger efficiency",
        "eff_photon":              "Combined photon efficiency",
        "final":                   "Final cross section result",
        "final_phenix":            "Comparison with PHENIX",
        "final_common_cluster":    "Common cluster yield (data vs MC)",
        "final_tight_iso_cluster": "Tight iso cluster yield (data vs MC)",
        "final_all":               "All correction steps",
    }
    if plot_type in labels:
        return labels[plot_type]
    m = re.match(r"iso_ET_(tight|nontight|comb)_pt(\d+)$", plot_type)
    if m:
        kind_labels = {"tight": "Tight", "nontight": "Non-tight", "comb": "Combined"}
        return f"IsoET {kind_labels[m.group(1)]}, pT bin {m.group(2)}"
    m = re.match(r"h1D_iso(_fit)?_\w+_(\d+)_(\d+)$", plot_type)
    if m:
        kind = "fit " if m.group(1) else ""
        return f"IsoET {kind}template (pT bins {m.group(2)}\\textendash {m.group(3)})"
    return escape_latex(plot_type)


def copy_figures(all_src_paths, dest_figures_dir):
    """Copy source PDFs into dest_figures_dir, skipping already-current copies."""
    os.makedirs(dest_figures_dir, exist_ok=True)
    copied = 0
    for src in all_src_paths:
        dst = os.path.join(dest_figures_dir, os.path.basename(src))
        if not os.path.exists(dst) or os.path.getmtime(src) > os.path.getmtime(dst):
            shutil.copy2(src, dst)
            copied += 1
    print(f"Copied {copied} figure(s) to {dest_figures_dir}/  "
          f"({len(all_src_paths) - copied} already up to date)")


def compile_pdf(tex_file):
    """Run pdflatex twice (for TOC) from the directory containing the tex file."""
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
# Comparison-specific logic
# ---------------------------------------------------------------------------

def get_pair_diff(suffix_a, suffix_b):
    """Return unified diff lines between config_A.yaml and config_B.yaml."""
    path_a = os.path.join(CONFIG_DIR, f"config_{suffix_a}.yaml")
    path_b = os.path.join(CONFIG_DIR, f"config_{suffix_b}.yaml")

    lines_a, lines_b = [], []
    if os.path.exists(path_a):
        with open(path_a) as fh:
            lines_a = fh.readlines()
    if os.path.exists(path_b):
        with open(path_b) as fh:
            lines_b = fh.readlines()

    if not lines_a and not lines_b:
        return []

    return list(difflib.unified_diff(
        lines_a, lines_b,
        fromfile=f"config_{suffix_a}.yaml",
        tofile=f"config_{suffix_b}.yaml",
        n=2,
    ))


def _resolve_plot(suffix):
    """Return a dict mapping pattern-key -> path (or None) for one suffix.

    Handles both PLOT_GROUPS patterns and IsoET patterns.
    """
    mapping = {}

    for _, patterns in PLOT_GROUPS:
        for pat in patterns:
            fname = pat.replace("{s}", suffix)
            fpath = os.path.join(FIGURES_DIR, fname)
            mapping[pat] = fpath if os.path.exists(fpath) else None

    for ipt in range(N_PT_BINS):
        for plot_type in ("tight", "nontight", "comb"):
            fname = f"iso_ET_{plot_type}_pt{ipt}_{suffix}.pdf"
            pat   = f"iso_ET_{plot_type}_pt{{ipt}}_{{}}"  # just a unique key
            fpath = os.path.join(FIGURES_DIR, fname)
            mapping[fname] = fpath if os.path.exists(fpath) else None

    return mapping


def collect_pair_plots(suffix_a, suffix_b):
    """Return ordered list of (group_label, [(plot_key, path_a|None, path_b|None)]).

    Only includes groups where at least one plot exists in either config.
    """
    result = []

    for label, patterns in PLOT_GROUPS:
        entries = []
        for pat in patterns:
            fname_a = pat.replace("{s}", suffix_a)
            fname_b = pat.replace("{s}", suffix_b)
            path_a = os.path.join(FIGURES_DIR, fname_a)
            path_b = os.path.join(FIGURES_DIR, fname_b)
            pa = path_a if os.path.exists(path_a) else None
            pb = path_b if os.path.exists(path_b) else None
            if pa is not None or pb is not None:
                entries.append((pat, pa, pb))
        if entries:
            result.append((label, entries))

    # IsoET: tight only shown; still scan nontight/comb to skip them
    isoet_entries = []
    for ipt in range(N_PT_BINS):
        for plot_type in ("tight", "nontight", "comb"):
            fname_a = f"iso_ET_{plot_type}_pt{ipt}_{suffix_a}.pdf"
            fname_b = f"iso_ET_{plot_type}_pt{ipt}_{suffix_b}.pdf"
            path_a = os.path.join(FIGURES_DIR, fname_a)
            path_b = os.path.join(FIGURES_DIR, fname_b)
            pa = path_a if os.path.exists(path_a) else None
            pb = path_b if os.path.exists(path_b) else None
            if plot_type == "tight" and (pa is not None or pb is not None):
                isoet_entries.append((fname_a, pa, pb))

    if isoet_entries:
        result.append(("Isolation ET", isoet_entries))

    return result


def write_pair_figure_block(f, plot_entries, suffix_a, suffix_b, width=0.48):
    """Write one figure environment per plot type, with A left and B right.

    plot_entries: list of (plot_key, path_a|None, path_b|None)
    Each subfloat either includes the graphic or shows a "not found" note.
    Entire figure is skipped only if both sides are missing.
    """
    esc_a = escape_latex(suffix_a)
    esc_b = escape_latex(suffix_b)

    for _key, path_a, path_b in plot_entries:
        if path_a is None and path_b is None:
            continue  # shouldn't happen after collect_pair_plots filtering

        # Derive a representative path for caption purposes (prefer whichever exists)
        cap_ref = path_a if path_a is not None else path_b
        cap_ref_suffix = suffix_a if path_a is not None else suffix_b

        cap_a = make_caption(cap_ref, cap_ref_suffix) + f" --- \\texttt{{{esc_a}}}"
        cap_b = make_caption(cap_ref, cap_ref_suffix) + f" --- \\texttt{{{esc_b}}}"

        f.write("\\begin{figure}[htbp]\n")
        f.write("  \\centering\n")

        # Left subfloat (config A)
        if path_a is not None:
            tex_a = "figures/" + os.path.basename(path_a)
            f.write(f"  \\subfloat[{cap_a}]{{\\includegraphics[width={width}\\linewidth]{{{tex_a}}}}}\n")
        else:
            f.write(f"  \\subfloat[{cap_a}]{{\\textit{{not found}}}}\n")

        f.write("  \\hfill\n")

        # Right subfloat (config B)
        if path_b is not None:
            tex_b = "figures/" + os.path.basename(path_b)
            f.write(f"  \\subfloat[{cap_b}]{{\\includegraphics[width={width}\\linewidth]{{{tex_b}}}}}\n")
        else:
            f.write(f"  \\subfloat[{cap_b}]{{\\textit{{not found}}}}\n")

        f.write("\\end{figure}\n\n")


def generate_tex(pairs, output_tex):
    """Write the LaTeX comparison document; returns list of all src plot paths."""
    all_src_paths = []

    with open(output_tex, "w") as f:
        f.write(r"""\documentclass[11pt,letterpaper]{article}
\usepackage{graphicx}
\usepackage{subfig}
\usepackage[margin=1in]{geometry}
\usepackage{hyperref}
\usepackage{booktabs}
\usepackage{fancyhdr}
\usepackage{alltt}
\usepackage{xcolor}
\definecolor{diffadd}{rgb}{0.0, 0.50, 0.15}
\definecolor{diffdel}{rgb}{0.75, 0.10, 0.10}
\definecolor{diffhunk}{rgb}{0.10, 0.30, 0.70}
\definecolor{difffile}{rgb}{0.40, 0.40, 0.40}

\pagestyle{fancy}
\fancyhf{}
\rhead{BDT Config Comparisons}
\lhead{\nouppercase{\leftmark}}
\rfoot{\thepage}

\begin{document}

\title{Selection Plots: Config Comparisons}
\date{\today}
\maketitle
\tableofcontents
\newpage

""")

        for suffix_a, suffix_b in pairs:
            esc_a = escape_latex(suffix_a)
            esc_b = escape_latex(suffix_b)
            section_title = f"\\texttt{{{esc_a}}} vs \\texttt{{{esc_b}}}"
            f.write(f"\\section{{{section_title}}}\n\n")

            # Config diff block
            diff_lines = get_pair_diff(suffix_a, suffix_b)
            if diff_lines:
                f.write("\\subsection*{Config diff (A $\\to$ B)}\n\n")
                write_diff_block(f, diff_lines)
            else:
                f.write("\\subsection*{Config diff}\n\n"
                        "\\textit{(config files not found or identical)}\n\n")

            plots_by_group = collect_pair_plots(suffix_a, suffix_b)

            if not plots_by_group:
                f.write("\\textit{No plots found for either configuration.}\n\n")
                f.write("\\clearpage\n\n")
                continue

            for group_label, entries in plots_by_group:
                f.write(f"\\subsection{{{group_label}}}\n\n")

                # Collect existing paths for copying
                for _key, pa, pb in entries:
                    if pa is not None:
                        all_src_paths.append(pa)
                    if pb is not None:
                        all_src_paths.append(pb)

                write_pair_figure_block(f, entries, suffix_a, suffix_b)
                f.write("\\clearpage\n\n")

        f.write("\\end{document}\n")

    total_figs = len(all_src_paths)
    print(f"Written: {output_tex}  ({len(pairs)} pair(s), {total_figs} figure(s))")
    return all_src_paths


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--pair", nargs=2, metavar=("A", "B"), action="append", default=[],
        help="Pair of config suffixes to compare (repeat for multiple pairs)",
    )
    parser.add_argument("--compile", action="store_true",
                        help="Run pdflatex after generating the .tex file")
    parser.add_argument("--output-dir", default=OUTPUT_DIR,
                        help=f"Output directory (default: {OUTPUT_DIR})")
    args = parser.parse_args()

    if not args.pair:
        parser.error("Provide at least one --pair A B")

    pairs = [tuple(p) for p in args.pair]

    # Validate configs exist; warn but proceed if missing
    for suffix_a, suffix_b in pairs:
        for suffix in (suffix_a, suffix_b):
            cfg = os.path.join(CONFIG_DIR, f"config_{suffix}.yaml")
            if not os.path.exists(cfg):
                print(f"Warning: config file not found: {cfg}")

    out_dir    = args.output_dir
    dest_fig   = os.path.join(out_dir, "figures")
    output_tex = os.path.join(out_dir, TEX_NAME)

    os.makedirs(out_dir, exist_ok=True)

    all_src_paths = generate_tex(pairs, output_tex)
    # Deduplicate while preserving order
    seen = set()
    unique_paths = []
    for p in all_src_paths:
        if p not in seen:
            seen.add(p)
            unique_paths.append(p)
    copy_figures(unique_paths, dest_fig)

    if args.compile:
        compile_pdf(output_tex)


if __name__ == "__main__":
    main()
