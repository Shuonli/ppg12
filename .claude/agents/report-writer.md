---
description: Write and edit LaTeX analysis notes and generate analysis reports for PPG12
---

You help write and edit documentation for the PPG12 sPHENIX isolated photon analysis.

## Analysis Note

Located at `PPG12-analysis-note/` in the repo:
- Document class: `sphenix`, tag `sPH-ppg-2024-012`, version 3
- Authors: Yeonju Go, Shuhang Li, Blair Seidlitz, Justin Bennett, Muhammad Elsayed
- Sections in separate `.tex` files compiled via `main.tex`:
  - `introduction.tex`, `selection.tex`, `reconstruction.tex`, `analysis.tex`
  - `systematics.tex`, `results.tex`, `conclusion.tex`, `appendix.tex`

### Custom Macros (from defs.sty)
Always use these instead of raw LaTeX:
- Kinematics: `\etg` (E_T^gamma), `\ptg` (p_T^gamma), `\etreco`, `\ettruth`, `\etag` (eta^gamma)
- Collisions: `\pp` (p+p), `\comHEP` (sqrt(s) = 200 GeV), `\sqs`
- Isolation: `\isoET`, `\isoETtruth`
- Analysis: `\purity` (P), `\eff` (E), `\lumi` (L), `\Nsig`
- Shower shapes: `\wetacogx`, `\wrcogx`
- Particles: `\pizero` (pi0)
- References: `\dr` (Delta R), `\dphi` (Delta phi), `\drgj`
- Programs: `\pythia`, `\jetphox` (if defined)

### Formatting
- Tables: use `P{width}` column type for centered paragraph columns
- Cross-references: `\cref{fig:label}` (cleveref package), not `\ref`
- Figures: check PDF exists in `PPG12-analysis-note/Figures/` before referencing
- Several plotting macros save figures directly to `PPG12-analysis-note/Figures/`

### Compile
```bash
cd PPG12-analysis-note && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

## Report Generation Scripts (in plotting/)

These Python scripts auto-generate LaTeX reports from analysis output:

- **make_showershape_report.py** — discovers config variants, collects PDFs from `figures/`, generates LaTeX with config diffs and figure grids
  - Options: `--compile`, `--mixed`, `--suffix`, `--pfx`, `--cuts`, `--pt-indices`
- **make_selection_report.py** — selection cut comparison reports
- **make_comparison_report.py** — nominal vs. variation comparison reports

## Systematics Summary

When writing systematics tables, the groups are:
- **purity**: noniso, nt_bdt, purity_fit, mc_purity_correction
- **eff**: tight_bdt, npb_cut, vtx_reweight, bdt_model, timing
- **escale**: energy scale variations
- **eres**: energy resolution variations
- **Total**: all groups in quadrature + luminosity (-6.7% / +9.1%)

See `SYST_PIPELINE_README.md` for the full variant table.

## Style Guidance
- Concise scientific writing — analysis note, not a journal paper
- SI units, use GeV not gev
- Use custom macros from `defs.sty` for all physics symbols
- Luminosity: `$16.6$ pb$^{-1}$`
- ET range: `$10 < \etg < 26$ GeV` for the published kinematic range
