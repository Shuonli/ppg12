# Stage 5: Plotting and Systematics

## Key Files

| File | Purpose |
|------|---------|
| `plotting/plot_final_selection.C` | Primary final cross-section plot (816 lines) |
| `plotting/calc_syst_bdt.py` | Automated systematic aggregation (534 lines) |
| `plotting/plotcommon.h` | Shared pT bins, frames, legend strings (92 lines) |
| `plotting/BlairUtils.C` | Error-band construction, TGraph helpers (543 lines) |
| `plotting/plot_efficiency.C` | Efficiency plots (258 lines) |
| `plotting/plot_purity_selection.C` | Purity with fit bands (114 lines) |
| `plotting/plot_sideband_selection.C` | ABCD yield plots (235 lines) |
| `plotting/plot_showershapes_variations.C` | Shower shapes per config (670 lines) |
| `plotting/plot_unfold_iter.C` | Unfolding convergence (356 lines) |
| `plotting/make_showershape_report.py` | LaTeX report for shower shape figures (447 lines) |
| `plotting/make_selection_report.py` | LaTeX report for selection/efficiency (413 lines) |
| `plotting/make_comparison_report.py` | Side-by-side comparison report (456 lines) |

## Shared Utilities

### plotcommon.h

```cpp
const int NptBins = 12;
const float ptRanges[NptBins + 1] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36};
```

Legend strings:
- `strleg1 = "#bf{#it{sPHENIX}} Internal"`
- `strleg2_1 = "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, 16.6 pb^{-1}"`

Helper: `calcDelta(h1, h2, name)` returns `(h1-h2, (h1-h2)/h2)`.

### BlairUtils.C

Key functions: `myTGraphErrorsDivide`, `myMakeBand`, `myAddtoBand`, `TH1TOTGraph`, `myText`, `myMarkerText`, `legStyle`.

## Final Cross-Section Plot (plot_final_selection.C)

```cpp
void plot_final_selection(string tune = "bdt_nom")
```

**Input:** `Photon_final_{tune}.root`, `syst_sum.root`, JETPHOX ROOT files, Werner Vogelsang text files.

**Output figures (all in `figures/`):**
1. `final_{tune}.pdf` -- Cross-section + Theory/Data ratio with systematic bands
2. `final_common_cluster_{tune}.pdf` -- Common cluster yield
3. `final_tight_iso_cluster_{tune}.pdf` -- Tight-iso yield
4. `final_all_{tune}.pdf` -- Analysis pipeline stages
5. `final_phenix_{tune}.pdf` -- PHENIX comparison

Reads luminosity from YAML config for legend text. Divides by `deta = 1.4` for d2sigma/deta/dET presentation.

## Systematic Aggregation (calc_syst_bdt.py)

```bash
python calc_syst_bdt.py [--results DIR] [--outdir rootFiles] [--figdir figures] [--nom bdt_nom]
```

Imports `VARIANTS`, `SYST_TYPES`, `SYST_GROUPS`, `FINAL_SYSTS`, `LUMI_SYST` from `efficiencytool/make_bdt_variations.py`.

### Three-Tier Aggregation

1. **Per-type:** Load `Photon_final_bdt_{variant}.root`, compute delta from nominal. Modes: `two_sided`, `one_sided`, `max`.
2. **Per-group:** Quadrature sum of member types. Groups: purity, eff, escale, eres, mbd, unfolding.
3. **Total:** Quadrature sum of all groups + flat luminosity uncertainty.

### Output

- `rootFiles/syst_bdt_{type}.root` -- per-type (h_dev_low, h_dev_high, h_dev_rel_low, h_dev_rel_high)
- `rootFiles/syst_bdt_{group}.root` -- per-group
- `rootFiles/syst_sum.root` -- total (consumed by plot_final_selection.C)
- `figures/syst_bdt_rel_{type}.pdf` -- per-type plots
- `figures/syst_bdt_breakdown.pdf` -- all groups on one canvas

## Other Important Plotting Macros

### Efficiency (plot_efficiency.C)

Input: `MC_efficiency_{suffix}.root`, `Photon_final_{suffix}.root`.
Output: `eff_reco_`, `eff_iso_`, `eff_id_`, `eff_total_`, `eff_mbd_`, `eff_photon_{suffix}.pdf`.

### Purity (plot_purity_selection.C)

Input: `Photon_final_{suffix}.root`, `Photon_final_{suffix}_mc.root`.
Output: `purity_{suffix}.pdf` with fit confidence bands.

### ABCD Sidebands (plot_sideband_selection.C)

Input: `data_histo_{suffix}.root`, `MC_efficiency_{suffix}.root`, `MC_efficiency_jet_{suffix}.root`.
Output: `et_sbs_{suffix}.pdf`, `et_sbs_ratio_{suffix}.pdf`, `leakage_fraction_et_{suffix}.pdf`.

### Shower Shapes (plot_showershapes_variations.C)

Config-driven: reads pT_bins, bdt_bins from YAML. Produces ratio plots with chi2/ndf. Output in `figures/{suffix}/`.

### Unfolding (plot_unfold_iter.C)

Shows Bayesian unfolding convergence over 10 iterations and closure tests.

## Report Generators

### make_showershape_report.py

Discovers `config_showershape*.yaml` variants, collects shower shape figures, generates LaTeX with 3x3 variable grids per (cut, pT bin).

### make_selection_report.py

Collects all selection/efficiency/purity/final plots per BDT config variant. Shows unified diff against nominal for non-nominal variants.

### make_comparison_report.py

Side-by-side two-column layout comparing pairs of config variants.

```bash
python3 make_comparison_report.py --pair bdt_nom bdt_tightbdt50 --pair bdt_nom bdt_npb03
```

## Output Conventions

- All figures: PDF format, saved to `plotting/figures/`
- Exception: per-variant shower shapes in `figures/{config_suffix}/`
- Naming: `{plot_type}_{suffix}.pdf` (e.g., `final_bdt_nom.pdf`)

## Batch Execution

```bash
# All selection plots for one config
bash make_selection_plots.sh bdt_nom

# All configs
bash make_all_bdt_selection.sh

# Legacy C++ systematics (alternative to calc_syst_bdt.py)
bash make_syst.sh
```

## Gotchas

- plotcommon.h legend says "R=0.3" (`strleg4`) but nominal reco isolation uses topo R=0.4 (`use_topo_iso: 2` in config)
- plotcommon.h legend says "16.6 pb^-1" but nominal config lumi is 16.2735 pb^-1

---

## Plotting Framework Overview

### Macro Categories

The plotting macros are organized by physics topic. Each macro includes `plotcommon.h` and is executed through the ROOT interpreter.

**Common utilities:**
- `plotcommon.h` -- shared pT bins (8-36 GeV, 12 bins), frame histograms, legend strings, `init_plot()`, `calcDelta()` helper
- `BlairUtils.C` -- error-band construction (`myMakeBand`, `myAddtoBand`), TGraphErrors division with error propagation (`myTGraphErrorsDivide`), TH1-to-TGraph conversion (`TH1TOTGraph`), text and marker drawing helpers
- `draw.C` -- generic drawing functions: `draw_1D_single_plot`, `draw_1D_multiple_plot_ratio` with configurable axes, legends, and ratio pads
- `sPhenixStyle.C` -- sPHENIX publication style (loaded via `plotcommon.h`)

**Final results:** `plot_final.C`, `plot_final_selection.C` -- cross-section vs ET with JETPHOX NLO comparison, ratio panels, systematic bands.

**Efficiency:** `plot_efficiency.C` -- reco, isolation, ID, total, MBD, and combined efficiencies vs truth ET.

**Purity:** `plot_purity.C`, `plot_purity_sim.C`, `plot_purity_selection.C` -- data and MC purity with fit confidence bands.

**ABCD sidebands:** `plot_sideband.C`, `plot_sideband_selection.C` -- ET spectra in ABCD regions, leakage fractions, B/A C/A D/A ratios.

**Shower shapes and BDT:** `plot_showershapes.C`, `plot_showershapes_selections.C`, `plot_showershapes_variations.C` -- distributions of shower shape variables (prob, CNN_prob, weta_cogx, wphi_cogx, et1-4, e11/e33, e17/e77, e32/e35, bdt) comparing data vs signal MC vs jet MC.

**Unfolding:** `plot_unfold_iter.C`, `Closure.C`, `plot_response.C` -- convergence, closure tests, response matrix visualization.

**Cluster and event timing:** `plot_cluster_timing.C`, `plot_cluster_mbd_time.C`, `plot_tower_timing.C`, `plot_leading_tower_mean_time_2d3d.C`, `plot_npb_time_bkgsub.C` -- timing distributions, cluster-MBD correlations, tower-level performance.

**Special studies:** `plot_reweight.C`, `plot_combine.C`, `plot_double_interaction.C`, `plot_rbrQA.C`, `plot_vertex.C`, `plot_xjg.C`, `plot_truth_iso.C`, `plot_satdiff.C`, `plot_mbd_sigma_efficiency.C`, `CONF_plots.C`.

**Systematic uncertainty macros:** Individual `syst_*.C` files (iso, nt, ntf, fit, escale, eres, eff, tight, fudge, mbd, nor, vertex, CNNvar) each produce `rootFiles/syst_*.root` with deviation histograms. `calcSyst.C` sums all sources in quadrature.

### Naming Conventions

All output figures are PDF format saved to `plotting/figures/`, with the pattern `{plot_type}_{suffix}.pdf`. Per-variant shower shape plots go under `figures/{config_suffix}/`. NLO comparison macros live in the `sphenix_nlo/` subdirectory.

## Shower Shape Plotting Updates

### plot_showershapes_selections.C Structure

This macro consolidates all shower shape plots into a single file, organized in four parts:

**Part 0 -- Standard data/MC comparison (original plots).** Normalized 1D distributions of all shower shape variables, overlaying data, signal MC, and background MC. Produced for four cut levels: cut0 (no nbkg weighting), cut1 (with nbkg weighting), cut2 (tight selection), cut3 (non-tight selection). Also produces background profile plots showing mean isolation ET vs each shower shape variable with the correlation factor displayed. Output filenames: `dis_{variable}_eta{N}_pt{N}_cut{N}.pdf` and `pfx_{variable}_eta{N}_pt{N}_cut{N}.pdf`.

**Part 1 -- BDT overlay plots (new).** Shows how the isolation ET distribution changes across three BDT score bins: 0.0-0.3, 0.3-0.7, and 0.7-1.0. This reveals the correlation between BDT score and isolation behavior. Output: `iso_BDT_overlay_eta{N}_pt{N}.pdf`. Requires BDT-binned histograms from `ShowerShapeCheck.C` (the `h2d_bdt` histograms). The macro handles missing histograms gracefully if `ShowerShapeCheck.C` has not been rerun with BDT binning enabled.

**Part 2 -- Isolation region overlay plots (new).** Compares shower shape distributions between isolated (cut2, tight) and non-isolated (cut3, non-tight) regions for background MC only. Shows how shower shape variables differ between the signal-like and sideband regions. Output: `{variable}_iso_overlay_eta{N}_pt{N}.pdf`.

**Part 3 -- Data vs MC isolation ET comparison (new).** Compares isolation ET distributions between data and inclusive jet MC for the tight selection. Used to validate the background model in the signal region. Output: `isoET_datamc_eta{N}_pt{N}.pdf`.

Total output: approximately 500+ PDF files across all eta bins, pT bins, variables, and cut levels.

### Inclusive File Handling

The plotting macro was updated to use `_inclusive` file suffixes for combined MC files. This resolved an issue where non-inclusive files produced by older code did not contain the `h2d_bdt` histograms needed by Part 1. The `hadd` commands in `run_showershape.sh` were also updated to merge inclusive files:

```bash
hadd -f MC_efficiencyshower_shape_signal_inclusive.root \
    MC_efficiencyshower_shape_photon5_inclusive.root \
    MC_efficiencyshower_shape_photon10_inclusive.root \
    MC_efficiencyshower_shape_photon20_inclusive.root

hadd -f MC_efficiencyshower_shape_jet_inclusive.root \
    MC_efficiencyshower_shape_jet30_inclusive.root \
    MC_efficiencyshower_shape_jet10_inclusive.root \
    MC_efficiencyshower_shape_jet15_inclusive.root \
    MC_efficiencyshower_shape_jet20_inclusive.root
```

### Shower Shape Variation Workflow

The full workflow for producing shower shape comparison plots across systematic variations has two steps.

**Step 1 -- Produce figures with `plot_showershapes_variations.C`.** This macro loops over all `config_showershape*.yaml` files found in `../efficiencytool/` and writes PDFs to `figures/{suffix}/`.

```bash
# All variants, nominal MC
root -l -b -q 'plot_showershapes_variations.C'

# Single variant
root -l -b -q 'plot_showershapes_variations.C("config_showershape_nom.yaml")'
```

The `use_mixed` parameter controls whether to read combined single+double interaction MC files instead of per-suffix files:

```bash
# All variants, combined single+double interaction MC
root -l -b -q 'plot_showershapes_variations.C("config_showershape.yaml", true)'
```

When `use_mixed=true`:
- Signal MC is read from `MC_efficiencyshower_shape_photon10_combined_showershape.root`
- Background MC is read from `MC_efficiencyshower_shape_jet12_combined_inclusive_showershape.root` (inclusive only; no non-inclusive combined file exists -- the code falls back gracefully)
- Figure filename prefix becomes `dis_mixed_` instead of `dis_`

**Step 2 -- Build the LaTeX report with `make_showershape_report.py`.** Discovers all `config_showershape*.yaml` variants, collects the PDFs from `figures/`, and generates a multi-section LaTeX document with config diffs and figure grids.

Key CLI flags:

| Flag | Effect |
|------|--------|
| `--compile` | Run pdflatex after generating .tex |
| `--suffix SUFFIX` | Process only one variant |
| `--mixed` | Include double-interaction mixed plots (prefix `dis_mixed`) |
| `--pfx PFX [PFX ...]` | Explicit prefix list (default: `dis`) |
| `--cuts CUT [CUT ...]` | Subset of cut levels: cut0, cut1, cut2, cut3 |
| `--pt-indices N [N ...]` | pT bin indices to include (default: 1 2 3) |
| `--output-dir DIR` | Override output directory (default: `showershape_report/`) |

Example usage:

```bash
# Nominal plots only
python3 make_showershape_report.py --compile

# Nominal + double-interaction mixed side by side
python3 make_showershape_report.py --mixed --compile

# Mixed only
python3 make_showershape_report.py --pfx dis_mixed --compile

# Single variant, both prefixes
python3 make_showershape_report.py --suffix showershape --mixed --compile
```

When multiple prefixes are active (e.g. `dis` + `dis_mixed`), the report adds a subsection per prefix and demotes cut sections to subsubsections. Figure captions include the prefix label (Nominal or Double-interaction mixed).

Output structure:

```
showershape_report/
  showershape_report.tex
  showershape_report.pdf
  {suffix}/
    showershape_{suffix}.tex
    figures/
```

### Generating BDT-Binned Histograms

The BDT overlay plots (Part 1) require BDT-binned histograms that are produced by `ShowerShapeCheck.C`. To generate them:

```bash
cd efficiencytool
root -l -b -q 'ShowerShapeCheck.C("config_nom.yaml", "photon20")'
root -l -b -q 'ShowerShapeCheck.C("config_nom.yaml", "jet10")'
root -l -b -q 'ShowerShapeCheck.C("config_nom.yaml", "data")'
```

After running `ShowerShapeCheck.C` and merging with `run_showershape.sh`, the combined inclusive files will contain the `h2d_bdt` histograms needed by `plot_showershapes_selections.C`.
