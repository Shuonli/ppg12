# Plotting Directory Report -- PPG12 Isolated Photon Cross-Section

## 1. Purpose

The `plotting/` directory produces all figures for the PPG12 isolated photon cross-section analysis: final cross-section vs JETPHOX NLO, PHENIX comparison, ABCD sideband yields, purity, efficiency, unfolding, response matrices, shower-shape distributions, systematic uncertainty breakdowns, BDT model comparisons, isolation ET fits, MBD timing/pileup studies, and automated LaTeX reports collecting these figures by config variant.

---

## 2. Key Files (every .C, .py, .h with approximate line count)

### Shared Utilities
| File | Lines | Description |
|------|-------|-------------|
| `plotcommon.h` | 93 | Shared pT bins, frames, legend strings, `init_plot()`, `calcDelta()` |
| `BlairUtils.C` | 544 | TGraph error-band construction, text/legend/marker drawing helpers |
| `sPhenixStyle.C` | (external) | sPHENIX publication style, loaded from CVMFS |
| `draw.C` | 330 | Legacy drawing helper (used only by `Closure.C`) |

### Final Cross-Section Plots
| File | Lines | Description |
|------|-------|-------------|
| `plot_final_selection.C` | 816 | **Primary**: final cross-section, JETPHOX+Werner NLO comparison, PHENIX overlay, analysis pipeline steps |
| `plot_final.C` | 0 (empty) | Placeholder (content is now in `plot_final_selection.C`) |
| `plot_final.C~` | 205 | Emacs backup of an older version |
| `plot_final_backup.C` | 783 | Earlier version of the final plot (hardcoded lumi 49.562 pb-1) |
| `plot_final_backup250327_v1.C` | 934 | Older backup (March 2025) |

### Systematics
| File | Lines | Description |
|------|-------|-------------|
| `calc_syst_bdt.py` | 534 | **Primary**: automated BDT systematic pipeline (per-type, group, total + lumi) |
| `calcSyst.C` | 208 | Legacy C++ systematics aggregation (reads syst_*.root from rootFiles/) |
| `syst_escale.C` | 161 | Energy scale systematic |
| `syst_eres.C` | 161 | Energy resolution systematic |
| `syst_iso.C` | 161 | Isolation cut systematic |
| `syst_tight.C` | 161 | Tight BDT cut systematic |
| `syst_mbd.C` | 161 | MBD efficiency systematic |
| `syst_ntf.C` | 164 | Non-tight fit systematic |
| `syst_fit.C` | 187 | Purity fit systematic |
| `syst_nor.C` | 144 | Unfolding (normalization) systematic |
| `syst_fudge.C` | 144 | Fudge factor systematic |
| `syst_eff.C` | 90 | Efficiency group (quadrature of tight + fudge) |
| `syst_purity.C` | 90 | Purity group (quadrature of iso + ntf + fit) |
| `syst_vertex.C` | 144 | Vertex systematic |
| `syst_nt.C` | 143 | Non-tight systematic (legacy) |
| `syst_CNNvar.C` | 121 | CNN variation systematic (legacy) |
| `examplecalcSyst.C` | 1295 | Example/template for systematic calculation |

### Efficiency and Purity
| File | Lines | Description |
|------|-------|-------------|
| `plot_efficiency.C` | 258 | Reco/iso/ID/total/MBD efficiencies, MBD fraction breakdown |
| `plot_purity.C` | 64 | Data purity with/without signal leakage correction |
| `plot_purity_sim.C` | 50 | MC purity (truth, with/without leakage correction) |
| `plot_purity_selection.C` | 114 | Purity with fit confidence bands (parameterized suffix) |
| `plot_purity_sim_selection.C` | 70 | MC purity with parameterized suffix |
| `plot_mc_purity_correction.C` | 53 | MC truth/purity-fit ratio vs ET |

### ABCD Sideband
| File | Lines | Description |
|------|-------|-------------|
| `plot_sideband.C` | 275 | Data ABCD yields and B/A, C/A, D/A ratios (hardcoded paths) |
| `plot_sideband_selection.C` | 235 | Data ABCD yields (parameterized suffix) |
| `plot_sideband_sim.C` | 177 | MC jet ABCD yields (hardcoded paths) |
| `plot_sideband_sim_selection.C` | 183 | MC jet ABCD yields (parameterized suffix) |
| `plot_sideband_compare_left_right.C` | 356 | Left/right sideband comparison |

### Shower Shape and BDT
| File | Lines | Description |
|------|-------|-------------|
| `plot_showershapes.C` | 379 | Shower-shape 2D profiles (data/sig/bkg), saves to analysis note |
| `plot_showershapes_variations.C` | 670 | Per-config-variant shower shapes with ratio plots and chi2/ndf |
| `plot_showershapes_selections.C` | 1027 | Shower shapes with selection cuts, double-interaction mixing |
| `plot_bdt_variations.C` | 197 | Compare cross-section across BDT model variants |
| `chi2_weight_scan.C` | 835 | Chi2 scan of double-interaction weight for data-MC agreement |

### Isolation ET
| File | Lines | Description |
|------|-------|-------------|
| `plot_isoET.C` | 310 | Isolation ET distributions (tight/nontight), data vs MC, per pT bin |
| `CONF_plots.C` | 228 | Isolation ET template fits in summed pT ranges |
| `plot_background_recoisoET_overlay.C` | 148 | Reco isoET overlay for background studies |

### Unfolding
| File | Lines | Description |
|------|-------|-------------|
| `plot_unfold_iter.C` | 356 | Bayesian unfolding convergence (10 iterations), closure tests |
| `Closure.C` | 110 | Full/half closure test using RooUnfold |
| `plot_response.C` | 58 | Response matrix visualization (reweighted and unweighted) |
| `plot_reweight.C` | 114 | Unfolding prior comparison (data unfolded vs MC prior) |

### Vertex and MC Combination
| File | Lines | Description |
|------|-------|-------------|
| `plot_vertex.C` | 39 | Data vs MC vertex z distributions and ratio |
| `plot_vertex_check.C` | 283 | Detailed vertex studies, saves to analysis note |
| `plot_combine.C` | 182 | Photon5/10/20 MC sample combination and fit residuals |
| `plot_SB.C` | 63 | Signal/background ratio vs ET from shower shape MC |

### Timing and QA
| File | Lines | Description |
|------|-------|-------------|
| `plot_tower_timing.C` | 1695 | Tower-level timing studies (largest file in directory) |
| `plot_cluster_timing.C` | 509 | Cluster timing distributions |
| `plot_clustertime.C` | 160 | Cluster time plots (legacy version) |
| `plot_cluster_timing_compare_leading.C` | 147 | Leading tower timing comparison |
| `plot_cluster_mbd_time.C` | 283 | Cluster-MBD timing correlations |
| `plot_time_energy_corr.C` | 129 | Time-energy correlation |
| `plot_time_energy_corr_fixtime.C` | 104 | Fixed-time energy correlation |
| `plot_leading_tower_mean_time_2d3d.C` | 331 | 2D/3D tower timing maps |
| `plot_cluster_rbr_QA.C` | 397 | Run-by-run QA for clusters |
| `plot_rbrQA.C` | 126 | Run-by-run QA (compact version) |
| `plot_saturation.C` | 491 | Tower saturation studies |
| `plot_satdiff.C` | 91 | Saturation difference plots |
| `plot_mbd_sigma_efficiency.C` | 164 | MBD sigma selection efficiency |

### Double Interaction / Pileup
| File | Lines | Description |
|------|-------|-------------|
| `plot_double_interaction.C` | 630 | Full toy double-interaction visualization |
| `plot_npb_time_bkgsub.C` | 299 | NPB timing with background subtraction, saves to analysis note |

### Other Physics Plots
| File | Lines | Description |
|------|-------|-------------|
| `plot_pid.C` | 120 | Particle ID composition in tight-iso region (MC jet) |
| `plot_photonjeteff.C` | 57 | Photon+jet back-to-back efficiency |
| `plot_truth_iso.C` | 205 | Truth isolation studies |
| `plot_xjg.C` | 143 | Photon-jet momentum balance (x_jg) |

### Python Utilities
| File | Lines | Description |
|------|-------|-------------|
| `calc_syst_bdt.py` | 534 | Automated BDT systematic uncertainty pipeline |
| `plot_var_comparison.py` | 142 | Two-variant relative difference comparison |
| `investigate_leakage.py` | 201 | ABCD leakage fraction investigation (uproot + matplotlib) |
| `make_showershape_report.py` | 447 | LaTeX report for shower shape figures |
| `make_selection_report.py` | 413 | LaTeX report for selection/efficiency/purity figures |
| `make_comparison_report.py` | 456 | Side-by-side comparison LaTeX report |

### Shell Scripts
| File | Lines | Description |
|------|-------|-------------|
| `make_all_bdt_selection.sh` | 13 | Loop over all config_bdt*.yaml, run selection plots |
| `make_all_selection.sh` | 13 | Same loop with different config pattern |
| `make_selection_plots.sh` | 35 | Run 11 plot macros in parallel for one suffix |
| `make_all_templatefit.sh` | 22 | Run CONF_plots.C for multiple pT-bin ranges |
| `make_syst.sh` | 24 | Run all syst_*.C macros, then group aggregation + total |

---

## 3. Shared Utilities

### plotcommon.h

**Constants:**
```cpp
const int NptBins = 12;
const float ptRanges[NptBins + 1] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35};
const float pTmin = 10;
const float pTmax = 30;
```

**Frames (created in `init_plot()`):**
- `frame_et_rec` -- TH1F, x: 7-50 GeV, titled `E_T^{gamma,rec}`
- `frame_et_truth` -- TH1F, x: 7-50 GeV, titled `E_T^{gamma,truth}`, y: "Efficiency" (0.2-1.1)
- `frame_isoET` -- TH1F, x: -5 to 15 GeV, titled `E_T^{iso}`
- `frame_iteration` -- TH1F, x: 0-10, for unfolding iteration plots
- `frame_response` -- TH2F, 7-50 x 7-50 GeV, for response matrices
- `lineone` / `linezero` -- TGraph horizontal lines at y=1 and y=0 (dashed)

**Legend strings:**
```cpp
strleg1  = "#bf{#it{sPHENIX}} Internal"
strleg2  = "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV"
strleg2_1 = "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, 16.6 pb^{-1}"
strleg3  = "|#it{#eta^{#gamma}}| < 0.7"
strleg4  = "#it{E}_{T}^{iso, R=0.3} < 4 GeV"
strleg5  = "         = 16.6 pb^{-1}"
strSigMC = "PYTHIA Signal"
strMC    = "PYTHIA8"
strIncMC = "PYTHIA Inclusive Jet"
```

**Helper function:**
```cpp
std::pair<TH1F *, TH1F *> calcDelta(TH1F *h1, TH1F *h2, std::string name)
```
Returns `(h1-h2, (h1-h2)/h2)` -- absolute and relative difference histograms.

### BlairUtils.C

All functions defined (in order of appearance):

| Function | Signature | Purpose |
|----------|-----------|---------|
| `ATLAS_LABEL` | `(x, y, color)` | Draw "ATLAS" label (legacy, unused) |
| `myTGraphErrorsDivide` (TGraphErrors) | `(g1, g2) -> TGraphErrors*` | Divide two TGraphErrors point-by-point with error propagation |
| `myTGraphErrorsDivide` (TGraphAsymmErrors) | `(g1, g2) -> TGraphAsymmErrors*` | Divide two TGraphAsymmErrors with asymmetric error propagation |
| `myMakeBand` | `(g0, g1, g2) -> TGraphAsymmErrors*` | Build asymmetric error band from central (g0) and two variations (g1, g2) |
| `myAddtoBand` | `(g1, g2) -> void` | Add variation g1 to existing band g2 in quadrature |
| `TH1TOTGraph` | `(h1) -> TGraphErrors*` | Convert TH1 to TGraphErrors |
| `myText` | `(x, y, color, text, tsize, isRightAlign)` | Draw text label at NDC position |
| `myBoxText` | `(x, y, boxsize, mcolor, bstyle, lcolor, lstyle, text)` | Draw colored box + text legend entry |
| `myMarkerText` | `(x, y, color, mstyle, text, msize, tsize)` | Draw marker + text legend entry |
| `myMarkerLineText` | `(x, y, msize, mcolor, mstyle, lcolor, lstyle, text, tsize, EX0)` | Draw marker+line + text legend entry |
| `mySmallMarkerLineText` | `(x, y, msize, mcolor, mstyle, lcolor, lstyle, text)` | Small version of marker+line text |
| `myOnlyBoxText` | `(x, y, boxsize, mcolor, lcolor, lstyle, text, tsize, bstyle)` | Box-only legend entry |
| `mySmallBoxText` | `(x, y, boxsize, mcolor, lcolor, lstyle, text)` | Small box legend entry |
| `myBoxTextAlpha` | `(x, y, boxsize, mcolor, falpha, lcolor, lstyle, text)` | Transparent box legend entry |
| `legStyle` | `(TLegend*, margin, textSize, textFont, head)` | Configure legend style (no border, no fill) |

### sPhenixStyle.C (external)

Located at `/cvmfs/sphenix.sdcc.bnl.gov/.../rootmacros/sPhenixStyle.C`. Called via `SetsPhenixStyle()` in `init_plot()`. Configures ROOT's global `gStyle` for sPHENIX publication standards: font (42/Helvetica), pad margins, tick marks, stat box suppression, title suppression, canvas dimensions.

---

## 4. Each Plotting Macro in Detail

### plot_final_selection.C (816 lines) -- PRIMARY FINAL PLOT

**Signature:** `void plot_final_selection(string tune = "bdt_nom")`

**Input files:**
- `efficiencytool/results/Photon_final_{tune}.root` -- data result
- `efficiencytool/results/Photon_final_{tune}_mc.root` -- MC closure
- `plotting/rootFiles/syst_sum.root` -- total systematics
- `NLO/rootFiles/jetPHOX_10.root` -- JETPHOX NLO (mu = E_T)
- `NLO/rootFiles/jetPHOX_05.root` -- JETPHOX NLO (mu = 0.5*E_T)
- `NLO/rootFiles/jetPHOX_20.root` -- JETPHOX NLO (mu = 2*E_T)
- `plotting/sphenix_nlo/photons_newphenix_sc{1,05,2}.dat` -- Werner Vogelsang NLO text files

**Histograms read from Photon_final_{tune}.root:**
- `h_unfold_sub_result` -- final unfolded, efficiency-corrected cross section
- `h_truth_pT_0` -- PYTHIA8 truth spectrum
- `h_common_cluster_0` -- clusters passing preliminary cuts
- `h_tight_iso_cluster_0` -- tight isolated clusters
- `h_data_sub_leak` -- purity-corrected data
- `h_unfold_sub_result_woeff` -- unfolded without efficiency correction

**Histograms read from syst_sum.root:**
- `h_sum_low`, `h_sum_high` -- absolute systematic deviations
- `h_sum_rel_low`, `h_sum_rel_high` -- relative systematic deviations

**Histograms read from NLO files:**
- `h_truth_pT` -- NLO cross section

**Luminosity:** Read from `efficiencytool/config_{tune}.yaml` (key `analysis.lumi`), fallback 49.562 pb-1. Formats legend string dynamically.

**Output figures (all saved to `figures/`):**
1. `final_{tune}.pdf` -- Two-pad: cross-section (log y) with data + JETPHOX + PYTHIA8 + Werner NLO; Theory/Data ratio panel with systematic bands
2. `final_common_cluster_{tune}.pdf` -- Common cluster yield: data vs MC with MC/data ratio
3. `final_tight_iso_cluster_{tune}.pdf` -- Tight-iso cluster yield: data vs MC with ratio
4. `final_all_{tune}.pdf` -- Analysis pipeline stages: common cluster -> tight iso -> purity corrected -> unfolded -> efficiency corrected
5. `final_phenix_{tune}.pdf` -- Data vs PHENIX published cross section (PRD 86 072008)

**PHENIX data:** 18 data points hardcoded (pT 5-26 GeV, E d3sigma/dp3 format), converted via `2*pi*pT` scaling.

**Werner NLO:** 38-point text files, `yield = (direct + fragmentation) * 2*pi*pT`, integrated over bin widths for ratio.

**Constants:** `deta = 1.4` (full eta width). All cross sections divided by deta for d2sigma/deta/dET presentation.

---

### calc_syst_bdt.py (534 lines) -- AUTOMATED SYSTEMATICS

**Usage:** `python calc_syst_bdt.py [--results DIR] [--outdir rootFiles] [--figdir figures] [--nom bdt_nom] [--histogram h_unfold_sub_result]`

**Architecture:** Imports `VARIANTS`, `SYST_TYPES`, `SYST_GROUPS`, `FINAL_SYSTS`, `LUMI_SYST` from `efficiencytool/make_bdt_variations.py`.

**3-tier aggregation:**
1. **Per-type:** For each syst_type in SYST_TYPES, loads variant result files `Photon_final_bdt_{variant_name}.root`, computes delta from nominal. Three modes:
   - `two_sided`: separate up/down variants; mirrors if one direction missing
   - `one_sided`: symmetric (low = high = |delta|)
   - `max`: bin-wise maximum of |delta| over all "max" variants
2. **Per-group:** For each group in SYST_GROUPS (e.g., purity, eff, escale, eres, mbd, unfolding), quadrature sum of member types
3. **Total:** Quadrature sum of all FINAL_SYSTS groups + flat luminosity uncertainty (LUMI_SYST up/down fractions)

**Quadrature combination:** bin-by-bin `sqrt(sum of squares)` for each of 4 histograms: h_al (abs low), h_ah (abs high), h_rl (rel low), h_rh (rel high).

**Output ROOT files (rootFiles/):**
- `syst_bdt_{type}.root` -- per syst_type (h_dev_low/high, h_dev_rel_low/high)
- `syst_bdt_{group}.root` -- per group
- `syst_bdt_total.root` -- total
- `syst_sum.root` -- total in format for plot_final_selection.C (h_sum_low/high, h_sum_rel_low/high)

**Output figures (figures/):**
- `syst_bdt_rel_{type}.pdf` -- per-type ratio plot (up=red, down=blue, auto-Y-scaled)
- `syst_bdt_breakdown.pdf` -- all groups + total on one canvas

**Group colors:**
- purity: kViolet+1
- eff: kAzure+7
- escale: kGreen-2
- eres: kOrange+1
- mbd: kViolet-1
- unfolding: kMagenta+2

---

### plot_efficiency.C (258 lines)

**Signature:** `void plot_efficiency(const std::string suffix = "bdt_nom")`

**Input files:**
- `efficiencytool/results/MC_efficiency_{suffix}.root`
- `efficiencytool/results/Photon_final_{suffix}.root`

**Objects read:**
- `eff_reco_eta_0`, `eff_iso_eta_0`, `eff_id_eta_0`, `eff_all_eta_0` -- TEfficiency objects
- `g_mbd_eff`, `g_mbd_eff_north`, `g_mbd_eff_south` -- TGraphAsymmErrors
- `h_truth_pT_vertexcut_0`, `h_truth_pT_vertexcut_mbd_cut_0`, `h_truth_pT_vertexcut_mbd_only_north_0`, `h_truth_pT_vertexcut_mbd_only_south_0`, `h_truth_pT_vertexcut_mbd_neither_0` -- MBD fraction histograms

**Output figures (`figures/`):**
1. `eff_reco_{suffix}.pdf` -- reconstruction efficiency
2. `eff_iso_{suffix}.pdf` -- isolation efficiency
3. `eff_id_{suffix}.pdf` -- identification efficiency
4. `eff_total_{suffix}.pdf` -- total efficiency
5. `eff_mbd_{suffix}.pdf` -- MBD trigger efficiency (north, south, combined)
6. `eff_mbd_frac_{suffix}.pdf` -- MBD hit configuration fractions
7. `eff_photon_{suffix}.pdf` -- Combined photon efficiency (eps_reco, eps_reco x eps_ID, eps_reco x eps_ID x eps_iso)

---

### plot_sideband_selection.C (235 lines)

**Signature:** `void plot_sideband_selection(const std::string suffix = "nomtest")`

**Input files:**
- `efficiencytool/results/data_histo_{suffix}.root` -- data
- `efficiencytool/results/MC_efficiency_{suffix}.root` -- signal MC
- `efficiencytool/results/MC_efficiency_jet_{suffix}.root` -- background MC

**Histograms read (all with `_0` eta suffix):**
- `h_tight_iso_cluster_0` (A), `h_tight_noniso_cluster_0` (B), `h_nontight_iso_cluster_0` (C), `h_nontight_noniso_cluster_0` (D)

**Output figures (`figures/`):**
1. `et_sbs_{suffix}.pdf` -- ABCD ET spectra (data), log scale
2. `et_sbs_ratio_{suffix}.pdf` -- B/A, C/A, D/A ratios (data + signal MC)
3. `leakage_fraction_et_{suffix}.pdf` -- signal leakage fractions from MC

---

### plot_purity_selection.C (114 lines)

**Signature:** `void plot_purity_selection(const std::string suffix = "bdt_nom")`

**Input files:**
- `efficiencytool/results/Photon_final_{suffix}.root`
- `efficiencytool/results/Photon_final_{suffix}_mc.root`

**Objects read:** `gpurity`, `gpurity_leak`, `grFineConf`, `grFineConf_leak`, `f_purity_fit`, `f_purity_leak_fit`, `g_purity_truth`, `f_purity_mc_fit`

**Output:** `figures/purity_{suffix}.pdf` -- purity with/without leakage correction, fit confidence band

---

### plot_showershapes_variations.C (670 lines)

**Signature:** `static bool plotOneConfig(const std::string &configname, bool use_mixed = false)`

**Config-driven:** Reads `pT_bins`, `bdt_bins` from YAML config.

**Input files (pattern based on config suffix):**
- `efficiencytool/results/data_histoshower_shape_{suffix}.root`
- `efficiencytool/results/MC_efficiencyshower_shape_signal_{suffix}.root`
- `efficiencytool/results/MC_efficiencyshower_shape_jet_inclusive_{suffix}.root`
- For mixed mode: `MC_efficiencyshower_shape_photon10_combined_{suffix}.root`, `MC_efficiencyshower_shape_jet12_combined_inclusive_{suffix}.root`

**Histograms read:** 2D profiles named `h2d_{variable}_{eta}_{pt}_{cut}` for variables: prob, CNN_prob, e ratios (e17_to_e77, e37_to_e77, e32_to_e35, e33_to_e35, e11_to_e33, etc.), tower energies (et1-4), widths (weta, wphi, w32, w52, w72), BDT score, NPB score.

**Output:** Per-config subdirectory `figures/{suffix}/` with one PDF per (prefix, variable, eta, pt, cut) combination. Naming: `{dis|dis_mixed}_{variable}_eta{ieta}_pt{ipt}_{cut}.pdf`

Features ratio panels with chi2/ndf displayed.

---

### plot_unfold_iter.C (356 lines)

**Input files:**
- `efficiencytool/results/Photon_final_bdt_nom.root`
- `efficiencytool/results/MC_response_bdt_nom.root`
- `efficiencytool/results/MC_efficiency_bdt_nom.root`

**Histograms:** `h_pT_truth_response_0`, `h_pT_reco_response_0`, `h_response_full_0`, `h_response_half_0`, etc.

**Output figures (`figures/`):**
- Full/half closure tests (10 iterations)
- `unfold_iter_bdt_nom.pdf` -- iteration convergence (relative change vs iteration)
- Response matrix projections and ratios

---

### CONF_plots.C (228 lines)

**Signature:** `void CONF_plots(string tune = "nom", int lowbin = 6, int highbin = 8)`

**Input files:** `data_histo_{tune}.root`, `MC_efficiency_{tune}.root`

**Histograms:** `h_tight_isoET_0_{ipt}`, `h_nontight_isoET_0_{ipt}` -- isolation ET per pT bin

**Output:** `figures/h1D_iso_{tune}_{lowbin}_{highbin}.pdf` and `h1D_iso_fit_{tune}_{lowbin}_{highbin}.pdf` -- isolation ET template fits with non-tight background shape, exponential+Gaussian models

---

### plot_isoET.C (310 lines)

**Signature:** `void plot_isoET(const std::string suffix, const int etabin = 0, const int rebin = 2)`

**Input:** `data_histo_{suffix}.root`, `MC_efficiency_jet_{suffix}.root`

**Histograms:** `h_tight_isoET_{etabin}_{ipt}`, `h_nontight_isoET_{etabin}_{ipt}`

**Output per pT bin:** `figures/iso_ET_tight_pt{ipt}_{suffix}.pdf`, `figures/iso_ET_nontight_pt{ipt}_{suffix}.pdf`, `figures/iso_ET_comb_pt{ipt}_{suffix}.pdf`

---

### plot_double_interaction.C (630 lines)

**Signature:** `void plot_double_interaction(const std::string &configname, std::vector<int> smear_levels = {})`

**Config-driven:** Reads pT_bins, eta_bins, output filenames from YAML.

**Input files:** `{eff_outfile}_double_interaction_check_{signal,jet}.root`, `{data_outfile}_double_interaction_check.root`

**Histograms:** `h_tight_{iso,noniso}_cluster_npb_{ieta}` (TH2D: ET vs NPB score), shower shape histograms per interaction type.

**Output:** `figures/double_interaction/` -- ABCD yields vs NPB, single vs double shower shapes, purity comparisons

---

### chi2_weight_scan.C (835 lines)

**Signature:** `void chi2_weight_scan(const std::string &configname, int nSteps, double wMin, double wMax, bool useBdtForBest, int ptBinMin, int ptBinMax)`

Scans double-interaction mixing weight, computes chi2/ndf between data and `(1-w)*nominal + w*double` MC for each shower-shape variable.

**Output:** `figures/chi2_weight_scan_{suffix}.pdf`, `figures/chi2_weight_scan_{suffix}.root`

---

### plot_var_comparison.py (142 lines)

**Usage:** `python plot_var_comparison.py VAR1 VAR2 [options]`

Loads `h_unfold_sub_result` from two `Photon_final_{VAR}.root` files, plots relative difference.

**Output:** `figures/compare_{VAR1}_vs_{VAR2}.pdf`

---

### investigate_leakage.py (201 lines)

Uses `uproot` + `matplotlib` to read leakage fractions, signal counts in ABCD regions, purity comparisons.

**Input:** `Photon_final_bdt_nom_mc.root` (h_leak_B/C/D, gpurity, gpurity_leak, g_purity_truth, g_mc_purity_fit_ratio), `MC_efficiency_bdt_nom.root`

**Output:** `figures/leakage_fractions_vs_pT.pdf`, `figures/signal_counts_ABCD.pdf`, `figures/purity_comparison.pdf`, `figures/leakage_sensitivity.pdf`

---

## 5. Systematics Aggregation (calc_syst_bdt.py)

### Three-level hierarchy:

**Level 1 -- Per syst_type:** Each type (e.g., `iso_up`, `tight_bdt_up`, `energyscale26`) is a single config variant that ran through the full pipeline. The script loads `Photon_final_bdt_{variant_name}.root`, computes `delta = variant - nominal` bin-by-bin.

**Level 2 -- Per group (SYST_GROUPS):** Quadrature combination of member types:
- `purity`: iso + ntf + fit variations
- `eff`: tight BDT + fudge variations
- `escale`: energy scale variations
- `eres`: energy resolution variations
- `mbd`: MBD efficiency variations
- `unfolding`: normalization/prior variations

**Level 3 -- Total (FINAL_SYSTS + luminosity):**
Quadrature sum of all groups, then adds flat luminosity uncertainty (from `LUMI_SYST` dict with `up`/`down` fractional values).

### Output files:
- `rootFiles/syst_bdt_{type}.root` -- 4 histograms: h_dev_low, h_dev_high, h_dev_rel_low, h_dev_rel_high
- `rootFiles/syst_bdt_{group}.root` -- same structure for each group
- `rootFiles/syst_bdt_total.root` -- total
- `rootFiles/syst_sum.root` -- h_sum_low/high, h_sum_rel_low/high (consumed by plot_final_selection.C)

### Legacy C++ pipeline (make_syst.sh):
Runs individual `syst_*.C` macros (each writes to `rootFiles/syst_{name}.root`), then `syst_purity.C` and `syst_eff.C` aggregate sub-groups, then `calcSyst.C` combines all groups into `rootFiles/syst_sum.root`. Same quadrature logic but hardcoded in C++.

---

## 6. Final Cross-Section (plot_final_selection.C)

### How it combines everything:

1. Opens `Photon_final_{tune}.root` for the unfolded, efficiency-corrected cross section (`h_unfold_sub_result`)
2. Loads systematic band from `rootFiles/syst_sum.root`
3. Divides cross section by `deta = 1.4` for d2sigma/deta/dET presentation
4. Scales systematics by 1/deta for absolute bands
5. Constructs TGraphAsymmErrors for data systematic band, NLO scale band, and data ratio band
6. JETPHOX NLO: three scale choices (mu = 0.5, 1.0, 2.0 * E_T^gamma) from ROOT files
7. Werner Vogelsang NLO: three scale files from text (direct + fragmentation yields), integrated over bin widths
8. PHENIX comparison: 18 hardcoded data points, converted from E d3sigma/dp3 to d2sigma/deta/dET via `2*pi*pT`
9. Bottom panel: Theory/Data ratio with systematic bands for both data and NLO
10. Reads luminosity from YAML config (`analysis.lumi`) for correct legend text

---

## 7. Report Generation

### make_showershape_report.py

**What it generates:** `showershape_report/showershape_report.tex` master document with per-variant sub-documents.

**Discovery:** Scans `efficiencytool/config_showershape*.yaml` for variants. Extracts suffixes. For non-reference configs, shows unified diff against `config_showershape.yaml`.

**Figure layout:** 3x3 grids of shower-shape variables (e11_to_e33, e32_to_e35, bdt, weta_cogx, wphi_cogx, et1, et2, et3, et4) for each (prefix, cut, pT bin). Copies PDFs from `figures/{suffix}/` to `showershape_report/{suffix}/figures/`.

**Cut labels:** cut0 (no nbkg cut), cut1 (nbkg cut), cut2 (tight photon ID), cut3 (non-tight sideband)

### make_selection_report.py

**What it generates:** `selection_report/selection_report.tex` collecting all selection/efficiency/purity/final plots for each BDT config variant.

**Plot groups in order:**
1. Sideband (Data): et_sbs, et_sbs_ratio, leakage_fraction_et
2. Sideband (MC): et_sbs_sim, et_sbs_ratio_sim, leakage_fraction_et_sim
3. Purity: purity, purity_sim
4. Efficiency: eff_reco, eff_iso, eff_id, eff_total, eff_mbd, eff_photon
5. Final result: final, final_phenix, final_common_cluster, final_tight_iso_cluster, final_all
6. Isolation fit: h1D_iso, h1D_iso_fit for multiple pT ranges
7. IsoET: iso_ET_tight_pt{ipt} per pT bin
8. Catch-all: any remaining figure files matching the suffix

For non-nominal variants, shows unified diff against `config_bdt_nom.yaml`.

### make_comparison_report.py

**What it generates:** `comparison_report/comparison_report.tex` with side-by-side two-column layout comparing pairs of config variants. Same PLOT_GROUPS structure. Each pair gets its own section with config diff.

**Usage:** `python3 make_comparison_report.py --pair bdt_nom bdt_tightbdt50 --pair bdt_nom bdt_npb03`

---

## 8. Output Conventions

### Figure output directory

Primary output: `plotting/figures/` (flat directory for most macros).

Exception subdirectories:
- `figures/{config_suffix}/` -- per-variant shower shape plots (from plot_showershapes_variations.C)
- `figures/double_interaction/` -- double interaction plots

### Naming conventions

Most macros follow the pattern: `{plot_type}_{suffix}.pdf` where suffix is the `var_type` from the analysis config (e.g., `bdt_nom`, `bdt_tightbdt50`).

Examples:
- `final_bdt_nom.pdf`
- `eff_reco_bdt_nom.pdf`
- `purity_bdt_nom.pdf`
- `et_sbs_bdt_nom.pdf`
- `syst_bdt_rel_{type}.pdf`
- `syst_bdt_breakdown.pdf`
- `iso_ET_tight_pt{ipt}_{suffix}.pdf`
- `h1D_iso_{tune}_{lowbin}_{highbin}.pdf`

### Format

All output figures are PDF. A few legacy PNG files exist in the directory but are not generated by current macros.

---

## 9. Config Dependencies

| Macro | Config Required |
|-------|----------------|
| `plot_final_selection.C` | `efficiencytool/config_{tune}.yaml` (reads `analysis.lumi`) |
| `plot_showershapes.C` | `efficiencytool/config_showershape.yaml` (reads `analysis.pT_bins`) |
| `plot_showershapes_variations.C` | Any `config_showershape*.yaml` (reads pT_bins, bdt_bins) |
| `plot_showershapes_selections.C` | `config_showershape*.yaml` |
| `plot_double_interaction.C` | `config_showershape.yaml` (reads output.eff_outfile, output.data_outfile, analysis.pT_bins, analysis.eta_bins) |
| `chi2_weight_scan.C` | `config_showershape*.yaml` |
| `calc_syst_bdt.py` | Imports metadata from `efficiencytool/make_bdt_variations.py` |
| `Closure.C` | `config.yaml` (reads analysis.eta_bins) |

All macros expecting a suffix parameter use it to construct file paths: `efficiencytool/results/Photon_final_{suffix}.root`, `data_histo_{suffix}.root`, `MC_efficiency_{suffix}.root`, etc.

---

## 10. Constants

### Hardcoded Luminosity
- `plot_final.C~` (backup): `datalumi = 16.8468 * 23. / 26.1` pb-1
- `plot_final_backup.C`: `datalumi = 49.562` pb-1
- `plot_final_selection.C`: reads from YAML config, fallback `49.562` pb-1
- `plotcommon.h` legend string: `16.6 pb^{-1}` (legacy string, not used by final plot)

### Style settings (applied by init_plot())
- `SetsPhenixStyle()` -- font 42, standard margins
- `gStyle->SetHatchesLineWidth(4)` -- for hatch patterns in systematic bands

### Legend text
All defined in `plotcommon.h`:
- `strleg1 = "#bf{#it{sPHENIX}} Internal"` -- collaboration tag
- `strleg2 = "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV"` -- collision system
- `strleg3 = "|#it{#eta^{#gamma}}| < 0.7"` -- eta acceptance
- `strleg4 = "#it{E}_{T}^{iso, R=0.3} < 4 GeV"` -- isolation cut

### Cross-section presentation
- `deta = 1.4` -- full pseudorapidity range (2 * 0.7), divides cross sections for d2sigma/deta/dET

---

## 11. Figures Saved to Analysis Note

The following macros save directly to `PPG12-analysis-note/Figures/`:

| Macro | Output Path |
|-------|-------------|
| `plot_showershapes.C` | `../PPG12-analysis-note/Figures/showershapes/` |
| `plot_showershapes_selections.C` | `../PPG12-analysis-note/Figures/showershapes_selections/` |
| `plot_npb_time_bkgsub.C` | `../PPG12-analysis-note/Figures/npb_time_bkgsub/` |
| `plot_vertex_check.C` | `../PPG12-analysis-note/Figures/` (default savePath parameter) |

Additionally, `plot_sideband.C` has `savePath="../PPG12-analysis-note/Figures/analysis/"` hardcoded.

All other macros save to `plotting/figures/` and require manual copying or report generation to reach the analysis note.

---

## 12. Data Files in sphenix_nlo/

NLO theory predictions from Werner Vogelsang, stored as text files:
- `photons_newphenix_sc{1,05,2}.dat` -- direct + fragmentation photon yields at 3 scale choices (mu = E_T, 0.5*E_T, 2*E_T)
- `jets_newphenix_sc{1,05,2}.dat` -- jet NLO predictions
- `pi0_newphenix_sc{1,05,2}.dat` -- pi0 NLO predictions
- `cms_raa_05.dat`, `result_*.dat` -- additional NLO data

Format: 38 rows, 4 columns (pT, direct_yield, frag_yield, other). Read by `plot_final_selection.C`.

---

## 13. Orchestration Scripts

### make_selection_plots.sh
Runs 11 ROOT macros in parallel for a given suffix:
```
plot_sideband_selection.C, plot_sideband_sim_selection.C,
plot_purity_selection.C, plot_purity_sim_selection.C,
plot_isoET.C, plot_efficiency.C, plot_final_selection.C,
CONF_plots.C (4 pT-bin ranges: 0-2, 3-5, 6-8, 9-9)
```

### make_all_bdt_selection.sh
Loops over all `config_bdt*.yaml` files, extracts suffix, calls `make_selection_plots.sh` for each.

### make_syst.sh
Runs all individual syst_*.C in parallel, waits, then runs group aggregators (syst_purity.C, syst_eff.C), waits, then runs calcSyst.C for total.

### make_all_templatefit.sh
Runs CONF_plots.C for 4 pT-bin ranges on specified suffixes.
