#!/usr/bin/env python3
"""Generate a LaTeX report collecting selection plots from all BDT config variants.

Usage:
    python3 make_selection_report.py [--compile] [--output-dir selection_report]

Outputs:
    selection_report/selection_report.tex
    selection_report/figures/<all referenced PDFs copied here>

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
TEX_NAME    = "selection_report.tex"

# Number of pT bins (matches NptBins in plotcommon.h)
N_PT_BINS = 11  # noqa: matches NptBins in plotcommon.h

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
    ("MC purity correction", [
        "mc_purity_correction_{s}.pdf",
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


def get_suffixes():
    """Extract sorted list of suffixes from config_bdt*.yaml files.

    The suffix is the filename stem after the first underscore, matching the
    shell logic in make_all_selection.sh:
        input="${config#*_}"   # strip up to first underscore (removes 'config_')
        input="${input%.yaml}" # strip .yaml
    """
    pattern = os.path.join(CONFIG_DIR, "config_bdt*.yaml")
    configs = sorted(glob.glob(pattern))
    suffixes = []
    for path in configs:
        basename = os.path.basename(path)           # e.g. config_bdt_nom.yaml
        stem = basename.split("_", 1)[1]            # strip 'config_' -> bdt_nom.yaml
        stem = stem.removesuffix(".yaml")           # -> bdt_nom
        # Only include per-period (_0rad / _1p5mrad) sections for the nominal
        # variant. For all other variants, keep only the all-range bare config
        # to avoid bloating the selection report.
        if stem.endswith("_0rad") or stem.endswith("_1p5mrad"):
            if stem not in ("bdt_nom_0rad", "bdt_nom_1p5mrad"):
                continue
        suffixes.append(stem)
    return suffixes


def get_config_diff(suffix):
    """Return lines to display for the config block in a section.

    For bdt_nom: returns the full file contents (serves as reference).
    For variants: returns the unified diff against config_bdt_nom.yaml.
    Returns an empty list if files are missing.
    """
    nom_path = os.path.join(CONFIG_DIR, "config_bdt_nom.yaml")
    var_path = os.path.join(CONFIG_DIR, f"config_{suffix}.yaml")

    if not os.path.exists(nom_path):
        return []

    with open(nom_path) as fh:
        nom_lines = fh.readlines()

    if suffix == "bdt_nom":
        return nom_lines

    if not os.path.exists(var_path):
        return []

    with open(var_path) as fh:
        var_lines = fh.readlines()

    diff = list(difflib.unified_diff(
        nom_lines, var_lines,
        fromfile="config_bdt_nom.yaml",
        tofile=f"config_{suffix}.yaml",
        n=2,       # lines of context around each change
    ))
    return diff


def collect_plots(suffix):
    """Return ordered list of (group_label, [existing_src_pdf_paths]) for a suffix."""
    seen = set()
    result = []

    for label, patterns in PLOT_GROUPS:
        found = []
        for pat in patterns:
            fname = pat.replace("{s}", suffix)
            fpath = os.path.join(FIGURES_DIR, fname)
            if os.path.exists(fpath) and fpath not in seen:
                found.append(fpath)
                seen.add(fpath)
        result.append((label, found))

    # IsoET tight/nontight/comb plots: not displayed in the report (they
    # multiply per-config page count by ~12). Pre-mark as seen to protect
    # against any future catch-all picking them up.
    for plot_type in ("tight", "nontight", "comb"):
        for ipt in range(N_PT_BINS):
            fname = f"iso_ET_{plot_type}_pt{ipt}_{suffix}.pdf"
            fpath = os.path.join(FIGURES_DIR, fname)
            if os.path.exists(fpath):
                seen.add(fpath)

    return result


def copy_figures(all_src_paths, dest_figures_dir):
    """Copy source PDFs into dest_figures_dir, skipping already-current copies."""
    os.makedirs(dest_figures_dir, exist_ok=True)
    copied = 0
    for src in all_src_paths:
        dst = os.path.join(dest_figures_dir, os.path.basename(src))
        # Only copy if source is newer or destination doesn't exist
        if not os.path.exists(dst) or os.path.getmtime(src) > os.path.getmtime(dst):
            shutil.copy2(src, dst)
            copied += 1
    print(f"Copied {copied} figure(s) to {dest_figures_dir}/  "
          f"({len(all_src_paths) - copied} already up to date)")


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
    # Strip trailing suffix to get the plot type key
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
        "mc_purity_correction":    "MC purity correction factor (g_{truth}^{MC}/f_{fit})",
    }
    if plot_type in labels:
        return labels[plot_type]
    # IsoET pattern: iso_ET_{tight|nontight|comb}_pt{N}
    m = re.match(r"iso_ET_(tight|nontight|comb)_pt(\d+)$", plot_type)
    if m:
        kind_labels = {"tight": "Tight", "nontight": "Non-tight", "comb": "Combined"}
        return f"IsoET {kind_labels[m.group(1)]}, pT bin {m.group(2)}"
    # CONF_plots pattern: h1D_iso[_fit]_{suffix}_{lowbin}_{highbin}
    m = re.match(r"h1D_iso(_fit)?_\w+_(\d+)_(\d+)$", plot_type)
    if m:
        kind = "fit " if m.group(1) else ""
        return f"IsoET {kind}template (pT bins {m.group(2)}\\textendash {m.group(3)})"
    return escape_latex(plot_type)


def figures_in_pairs(paths):
    """Yield successive pairs (or a single for the last item if odd count)."""
    it = iter(paths)
    for a in it:
        b = next(it, None)
        yield (a, b)


def write_figure_block(f, src_paths, suffix, width=0.48):
    """Write a LaTeX figure environment with up to two subfloats side-by-side.

    src_paths are the original source paths; includegraphics uses just
    figures/<basename> since the tex file lives in the output directory.
    """
    for a, b in figures_in_pairs(src_paths):
        tex_a = "figures/" + os.path.basename(a)
        f.write("\\begin{figure}[htbp]\n")
        f.write("  \\centering\n")
        cap_a = make_caption(a, suffix)
        f.write(f"  \\subfloat[{cap_a}]{{\\includegraphics[width={width}\\linewidth]{{{tex_a}}}}}\n")
        if b is not None:
            tex_b = "figures/" + os.path.basename(b)
            cap_b = make_caption(b, suffix)
            f.write("  \\hfill\n")
            f.write(f"  \\subfloat[{cap_b}]{{\\includegraphics[width={width}\\linewidth]{{{tex_b}}}}}\n")
        f.write("\\end{figure}\n\n")


def generate_tex(suffixes, output_tex):
    """Write the LaTeX document; returns list of all src plot paths referenced."""
    # Put bdt_nom first, then sort the rest alphabetically
    ordered = []
    if "bdt_nom" in suffixes:
        ordered.append("bdt_nom")
    ordered += sorted(s for s in suffixes if s != "bdt_nom")

    all_src_paths = []

    with open(output_tex, "w") as f:
        f.write(r"""\documentclass[11pt,letterpaper]{article}
\usepackage{graphicx}
\usepackage{subfig}
\usepackage{chngcntr}
\counterwithin{figure}{section}
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
\rhead{BDT Selection Plots}
\lhead{\nouppercase{\leftmark}}
\rfoot{\thepage}

\begin{document}

\title{Selection Plots: BDT Variations}
\date{\today}
\maketitle
\tableofcontents
\newpage

""")

        for suffix in ordered:
            plots_by_group = collect_plots(suffix)
            total_plots = sum(len(ps) for _, ps in plots_by_group)

            section_name = r"\texttt{" + escape_latex(suffix) + "}"
            # \clearpage flushes deferred floats from the previous section so
            # \phantomsection's hyperref anchor lands on THIS section's page,
            # not stuck on the prior section's overflow.
            f.write("\\clearpage\n")
            f.write("\\phantomsection\n")
            f.write(f"\\section{{{section_name}}}\n\n")

            # Config block: full file for nominal, unified diff for variants
            config_lines = get_config_diff(suffix)
            if config_lines:
                if suffix == "bdt_nom":
                    f.write("\\subsection*{Nominal configuration}\n\n")
                else:
                    f.write("\\subsection*{Changes from nominal}\n\n")
                write_diff_block(f, config_lines)
            elif suffix != "bdt_nom":
                f.write("\\subsection*{Changes from nominal}\n\n"
                        "\\textit{(config file not found)}\n\n")

            if total_plots == 0:
                f.write("\\textit{No plots found for this configuration.}\n\n")
                continue

            for group_label, paths in plots_by_group:
                if not paths:
                    continue
                all_src_paths.extend(paths)
                f.write(f"\\subsection{{{group_label}}}\n\n")
                write_figure_block(f, paths, suffix)
                f.write("\\clearpage\n\n")

        f.write("\\end{document}\n")

    print(f"Written: {output_tex}  ({len(ordered)} sections, {len(all_src_paths)} figures)")
    return all_src_paths


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


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--compile", action="store_true",
                        help="Run pdflatex after generating the .tex file")
    parser.add_argument("--output-dir", default=OUTPUT_DIR,
                        help=f"Output directory (default: {OUTPUT_DIR})")
    args = parser.parse_args()

    suffixes = get_suffixes()
    if not suffixes:
        print(f"No config_bdt*.yaml files found in {CONFIG_DIR}/")
        print("Run from the plotting/ directory, or check that configs exist.")
        return

    print(f"Found {len(suffixes)} BDT config variants:")
    for s in suffixes:
        print(f"  {s}")

    out_dir     = args.output_dir
    dest_fig    = os.path.join(out_dir, "figures")
    output_tex  = os.path.join(out_dir, TEX_NAME)

    os.makedirs(out_dir, exist_ok=True)

    all_src_paths = generate_tex(suffixes, output_tex)
    copy_figures(all_src_paths, dest_fig)

    if args.compile:
        compile_pdf(output_tex)


if __name__ == "__main__":
    main()
