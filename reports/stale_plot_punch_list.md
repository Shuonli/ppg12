# Stale-Plot Punch List

**Date:** 2026-04-23
**Scope:** Figures in `reports/*.tex` that are stale relative to the current pipeline state.

## Executive summary

- **Figures are NOT missing** — Agent W's "129 missing / 26%" initial inventory was a false positive because it didn't resolve `\graphicspath` / `\figdir` LaTeX macros that each report uses to find its figures. The figures exist under `plotting/figures/` and subdirs.
- **Figures ARE stale** — today's (2026-04-23) pipeline ROOT files (`efficiencytool/results/Photon_final_bdt_*.root`) were regenerated under the **allz nominal + post-Apr-20 full-run-range + tower-mask** pipeline. Most figures in reports date 2026-04-19 to 2026-04-22, so they were produced with pre-allz-nominal or pre-tower-mask inputs.
- **All required input ROOT files exist** for the plotting scripts — no pipeline rerun needed before the plotting step.
- **Action:** re-run 6 plotting scripts. The resulting PDFs overwrite the stale ones in `plotting/figures/` (and subdirs). Reports then re-render the new plots automatically on next pdflatex.

## Path-resolution convention per report

| Report | `\graphicspath` or `\figdir` | Resolved base dir |
|--------|------------------------------|-------------------|
| `tower_acceptance_2026-04-22.tex` | `\newcommand{\figdir}{../plotting/figures}` | `plotting/figures/` |
| `donut_iso_variants.tex` | `\graphicspath{{../plotting/figures/}}` | `plotting/figures/` |
| `purity_nonclosure_ntbdt.tex` | `\graphicspath{{figures/}{/gpfs/.../plotting/figures/}}` | `plotting/figures/` (fallback) |
| `sample_combining_check.tex` | `\graphicspath{{../plotting/figures/sample_combining_check/}}` | `plotting/figures/sample_combining_check/` |

## Per-report punch list

### 1. `tower_acceptance_2026-04-22.tex` (14 figures)

| figure | current mtime | plotting code | input ROOT | regen command (from repo root) |
|--------|--------------|---------------|------------|------------------------------|
| `official_tower_{preselect,common,tight}_{mc,data}.pdf` | 2026-04-22 | `plotting/quantify_tower_mismatch_official.C` | `efficiencytool/results/Photon_final_bdt_nom.root` (today 15:07) | `cd plotting && root -l -b -q quantify_tower_mismatch_official.C` |
| `alt_check_phisymm.pdf`, `phisymm_variant_bybin.pdf` | 2026-04-22 | `plotting/plot_phisymm_variant_bybin.C` | `Photon_final_bdt_mask_phisymm_{preselect,common,tight,or}.root` (today 12:26–28) | `cd plotting && root -l -b -q plot_phisymm_variant_bybin.C` |
| `final_bdt_*.pdf` (4 mask variants) | 2026-04-22 | part of `plot_phisymm_variant_bybin.C` | same | covered by above |

**Text impact:** LOW — +2.5% envelope already confirmed across all 4 mask variants; numbers unlikely to shift materially. Still re-run to reflect the latest `Photon_final_bdt_nom.root`.

### 2. `donut_iso_variants.tex` (73 figures)

| figure family | current mtime | plotting code | input ROOT | regen command |
|---------------|--------------|---------------|------------|---------------|
| `purity_sim_donutFull_scan.pdf`, `purity_sim_donutExcl_scan.pdf` | 2026-04-19 | `plotting/plot_purity_sim_donut_scan.C` | `Photon_final_bdt_donut{Full,Excl}_{005,0075,01,02}_mc.root` (today) | `cd plotting && root -l -b -q plot_purity_sim_donut_scan.C` |
| `isoET_donutFull_*.pdf`, `isoET_donutExcl_*.pdf` | 2026-04-19 | `plotting/plot_isoET_donut_mc_scan.C`, `plot_isoET_donut_scan.C` | same family | `cd plotting && root -l -b -q plot_isoET_donut_scan.C && root -l -b -q plot_isoET_donut_mc_scan.C` |
| `final_bdt_donut*.pdf`, `h1D_iso_bdt_*.pdf` | 2026-04-19 | scripts TBD (grep `plotting/` for specific names if needed) | `Photon_final_bdt_donut*.root` | grep-then-run |

**Text impact:** MED — donut variants are presented as cross-checks; numeric-value shifts under new iso params (0.502095 / 0.0433036) may invalidate the "purity recovery" claims slightly. Review numbers after regen.

### 3. `purity_nonclosure_ntbdt.tex` (18 figures)

| figure | current mtime | plotting code | input ROOT | regen command |
|--------|--------------|---------------|------------|---------------|
| `purity_sim_ntbdt_scan.pdf`, `purity_sim_flat_scan.pdf` | 2026-04-19 | `plotting/plot_purity_sim_ntbdt_flat_scan.C` | `Photon_final_bdt_nom_mc.root`, `Photon_final_bdt_ntbdtmin{05,10}_mc.root`, `Photon_final_bdt_flat_t{85,90,95}_nt{40,50,70}_mc.root` (all today) | `cd plotting && root -l -b -q plot_purity_sim_ntbdt_flat_scan.C` |
| `purity_closure_ntbdt_scan.pdf` | 2026-04-19 | `plotting/plot_purity_closure_ntbdt_flat_scan.C` | same family | `cd plotting && root -l -b -q plot_purity_closure_ntbdt_flat_scan.C` |
| `final_bdt_*.pdf` cross-section panels | 2026-04-19 | script TBD | Photon_final_bdt_flat_t*.root family | grep-then-run |

**Text impact:** HIGH — user explicitly DROPPED this as obsolete (per "drop #15 this is obsolete"). The report's findings (mean |NC| 0.062–0.073) may not reflect current pipeline. Consider whether to regen at all or just retire the report. **Recommendation:** mark report as historical/retired in `reports/`.

### 4. `sample_combining_check.tex` (12 figures)

| figure | current mtime | plotting code | input ROOT | regen command |
|--------|--------------|---------------|------------|---------------|
| `signal_truth_pT_overlay_{0rad,1p5mrad}.pdf`, `signal_ratio_per_sample_*.pdf`, `signal_boundary_zoom_*.pdf` | 2026-04-21 | `plotting/plot_sample_combining_signal.C` | `MC_efficiency_photon{5,10,20}_bdt_nom_{0rad,1p5mrad}.root` (today 14:36) | `cd plotting && root -l -b -q plot_sample_combining_signal.C` |
| `background_*.pdf` (6 plots) | 2026-04-21 | `plotting/plot_sample_combining_background.C` | `MC_efficiency_jet{8,12,20,30,40}_bdt_nom_{0rad,1p5mrad}.root` (today 14:36) | `cd plotting && root -l -b -q plot_sample_combining_background.C` |

**Text impact:** LOW — the report establishes combining arithmetic correctness at machine precision. Regen confirms the new `_nom_0rad`/`_nom_1p5mrad` feeder files still combine cleanly.

## Rest-of-reports survey

Agent W's data (adjusted for graphicspath resolution) suggests the remaining reports are mostly OK on disk but many figures are also 1–5 days old. Priority re-run list:

1. **High priority** (reports cited by the Phase 2 analysis note updates):
   - `post_preliminary_updates.tex` (if referenced in §5 prose)
   - `iso_cone_fix.tex`
   - `truth_vertex_closure_nominal.tex`
   - `double_interaction_note.tex` (DI consolidated — compiles to the appendix PDF)

2. **Medium priority** (diagnostic / audit reports):
   - `bdt_training_audit.tex`
   - `condor_pipeline_review_run28.md`
   - `xsec_weight_audit.md`

3. **Deprecate candidates** (per Shuonli's drops):
   - `purity_nonclosure_ntbdt.tex` (item #15 dropped as obsolete)

## Recommendation

Two-step execution:
1. **Step 1 (batch regen):** run the 6 scripts identified above via a shell loop. Log output to `reports/figures/regen_2026-04-23.log`. ~5–10 min total.
2. **Step 2 (audit):** re-run Agent W's mtime scan after Step 1 to confirm all 4 worst-hit reports now have figure mtimes > 2026-04-23 00:00. Compare selected figures against previous PDFs for unexpected numerical shifts (spot-check 3 figures per report).

No report text rewriting needed unless Step 2 reveals large numerical drifts.

## Cross-cutting note on `tower_acceptance_2026-04-22.tex`

The report's `\figdir{../plotting/figures}` macro is well-formed — the figures compile correctly. Agent W flagged this as "compilation issue" but that was wrong. If the report previously failed to compile, check for other causes (missing packages, broken `\ref`).

## Limitations of this audit

- Agents X and Y (plotting-code map + input-ROOT semantic audit) were killed before returning full output. This punch list is based on Agent W's figure inventory + direct grep + manual inspection for the 4 worst-case reports.
- Not every figure in every report was checked — only the ones Agent W flagged as "missing" (all of which were actually present but stale). A complete audit of the other 17 reports is deferred.
- "Text impact" assessments are heuristic; actual numerical drift must be verified post-regen.
