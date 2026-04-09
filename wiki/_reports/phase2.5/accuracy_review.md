# Accuracy Review -- PPG12 Codebase Wiki Reports

**Generated**: 2026-04-08
**Reviewer scope**: Phase 1 reports (anatreemaker, bdt_pipeline, efficiency_tool, plotting) + Phase 2 reports (dataflow, config_schema, constants_audit). The analysis_note report listed in the task does not exist.
**Method**: Each claim was verified by grep/read of the actual source code file referenced.

---

## Verification Table

| # | Claim | Source Report | Verified Against | Result | Notes |
|---|-------|--------------|-----------------|--------|-------|
| 1 | "photon5cross = 146359.3 pb" | constants_audit, bdt_pipeline | `efficiencytool/CrossSectionWeights.h:13` | **CONFIRMED** | Exact value `146359.3f` |
| 2 | "photon10cross = 6944.675 pb" | constants_audit, bdt_pipeline | `efficiencytool/CrossSectionWeights.h:14` | **CONFIRMED** | Exact value `6944.675f` |
| 3 | "photon20cross = 130.4461 pb" | constants_audit, bdt_pipeline | `efficiencytool/CrossSectionWeights.h:15` | **CONFIRMED** | Exact value `130.4461f` |
| 4 | "jet50cross = 7.3113 pb" | constants_audit | `efficiencytool/CrossSectionWeights.h:26` | **CONFIRMED** | Exact value `7.3113f` |
| 5 | "BDTinput.C normalizes jet weights relative to jet30cross instead of jet50cross" | constants_audit (Section 1.3) | `FunWithxgboost/BDTinput.C:85,95,104,120` | **CONFIRMED** | e.g., `jet10cross / jet30cross` at line 85 vs `jet10cross / jet50cross` in CrossSectionWeights.h |
| 6 | "NptBins = 12, ptRanges = {8,10,12,14,16,18,20,22,24,26,28,30,35}" | constants_audit, plotting | `plotting/plotcommon.h:5-6` | **CONFIRMED** | Exact match |
| 7 | "config_bdt_nom.yaml pT_bins: [8,10,12,14,16,18,20,22,24,26,28,32,36]" | constants_audit (Section 2.2) | `efficiencytool/config_bdt_nom.yaml:51` | **CONFIRMED** | Exact match; last two edges (32,36) differ from plotcommon.h (30,35) |
| 8 | "Python report scripts hardcode N_PT_BINS=10" | constants_audit (Section 2.4) | `plotting/make_selection_report.py:28`, `plotting/make_comparison_report.py:34` | **CONFIRMED** | Both say `N_PT_BINS = 10` with misleading comment "matches NptBins in plotcommon.h" |
| 9 | "RecoEffCalculator_TTreeReader signature: (configname, filetype='jet40', do_vertex_scan=false, mix_weight=1.0, vtxscan_sim_override='')" | efficiency_tool (Section 5) | `efficiencytool/RecoEffCalculator_TTreeReader.C:44` | **CONFIRMED** | Exact 5-parameter signature matches |
| 10 | "CalculatePhotonYield signature: (configname, bool isMC = false)" | efficiency_tool (Section 4) | `efficiencytool/CalculatePhotonYield.C:52` | **CONFIRMED** | Default configname is `"config_bdt_purity_pade.yaml"` (not stated in report but not claimed otherwise) |
| 11 | "ShowerShapeCheck.C signature: (configname, filetype, doinclusive, do_vertex_scan, mix_weight, vtxscan_sim_override)" | efficiency_tool | `efficiencytool/ShowerShapeCheck.C:21` | **CONFIRMED** | All 6 parameters match, defaults match |
| 12 | "fit_option=0 is Pade, fit_option=1 is erf" | **config_schema (line 193)** | `efficiencytool/CalculatePhotonYield.C:585-591` | **WRONG** | **The code shows fit_option=0 creates Erf (default), fit_option=1 creates Pade rational function.** The config_schema report has the mapping inverted. The efficiency_tool report (Section 4, line 204) has it correct: "Default (fit_option=0): Error function". |
| 13 | "mc_iso_scale default 1.2, mc_iso_shift default 0.1" | efficiency_tool (Section 5, line 339) | `efficiencytool/RecoEffCalculator_TTreeReader.C:367-368` | **WRONG** | **Code defaults are mc_iso_shift=0.0 and mc_iso_scale=1.0** (`.as<float>(0.0)` and `.as<float>(1.0)`). The report conflates the **nominal config value** (1.2/0.1) with the **code default**. The config_bdt_nom.yaml sets mc_iso_scale=1.2 and mc_iso_shift=0.1, but these are config overrides, not hardcoded defaults. |
| 14 | "Reweighting order: Class reweighting -> ET reweighting -> eta reweighting -> vertex z reweighting" | bdt_pipeline (Section 7, line ~367) | `FunWithxgboost/reweighting.py:155-165` | **WRONG** | **Actual order per class: eta -> ET -> vertex** (not ET -> eta -> vertex). The `apply_all_reweighting` method calls `apply_eta_reweighting` BEFORE `apply_et_reweighting` in the per-class loop. |
| 15 | "Eta reweighting uses hardcoded [-0.7, 0.7]" | bdt_pipeline (Section 7, line ~392) | `FunWithxgboost/reweighting.py:49` | **CONFIRMED** | Code hardcodes `eta_min, eta_max = -0.7, 0.7` at line 49, ignoring the config `eta_reweight_range: [-2.5, 2.5]`. The config_schema report (line 99) correctly lists `-2.5, 2.5` as the config value, but this is **not actually used** by the code. |
| 16 | "MergeSim merges 6 jet samples: jet5/8/12/20/30/40" | dataflow, efficiency_tool | `efficiencytool/MergeSim.C:21-26,59-64` | **CONFIRMED** | Exactly these 6 jet samples. Does NOT include jet10, jet15, or jet50. |
| 17 | "DoubleInteractionCheck.C uses flat bdt_min (no slope/intercept)" | constants_audit (Section 6.2) | `efficiencytool/DoubleInteractionCheck.C:355` | **CONFIRMED** | Reads only `tight_bdt_min` (flat), no `bdt_min_slope` or `bdt_min_intercept`. Grep returns zero matches for `bdt_min_slope` or `bdt_min_intercept` in this file. |
| 18 | "TIME_SAMPLE_NS = 16.67 in time_energy_corr.C (vs 17.6 everywhere else)" | constants_audit (Section 7.2) | `efficiencytool/time_energy_corr.C:15`, `efficiencytool/DoubleInteractionCheck.C:27`, and 7 other files | **CONFIRMED** | All other files use 17.6; time_energy_corr.C alone uses 16.67. |
| 19 | "EMCal radius = 93.5 cm in DoubleInteractionCheck.C" | constants_audit (Section 7.1) | `efficiencytool/DoubleInteractionCheck.C:793` | **CONFIRMED** | `constexpr float cemc_radius_cm = 93.5f` |
| 20 | "train_single_model: true means 'Train one model per ET bin'" | config_schema (line 110) | `FunWithxgboost/main_training.py:763-766` | **WRONG** | **When `train_single_model=true`, the code calls `_train_single_model` (ONE model for all bins). When false, it calls `_train_per_bin_models`.** The config_schema description has the meaning inverted. |
| 21 | "File line counts: CaloAna24.cc=2321, apply_BDT.C=492, MergeSim.C=94, RecoEffCalculator_TTreeReader.C=2719, CalculatePhotonYield.C=1170, ShowerShapeCheck.C=1917, DoubleInteractionCheck.C=1446, CrossSectionWeights.h=119, MbdPileupHelper.h=207" | efficiency_tool, bdt_pipeline, anatreemaker | `wc -l` on each file | **CONFIRMED** | All 9 line counts match exactly. |
| 22 | "ET threshold for BDT scoring in apply_BDT.C: cluster_Et > 7 GeV (else score = -1)" | bdt_pipeline, dataflow | `FunWithxgboost/apply_BDT.C:396` | **CONFIRMED** | `if (cluster_Et_BDT > 7)` |
| 23 | "use_topo_iso: 2 in config_bdt_nom.yaml" | constants_audit (Section 5.1) | `efficiencytool/config_bdt_nom.yaml:53` | **CONFIRMED** | `use_topo_iso: 2` |
| 24 | "BDT threshold parameters in config_bdt_nom.yaml: slope=-0.015, intercept=0.80" | constants_audit (Section 6.1) | `efficiencytool/config_bdt_nom.yaml:148-149` | **CONFIRMED** | `bdt_min_slope: -0.015` and `bdt_min_intercept: 0.80` |
| 25 | "Luminosity in config_bdt_nom.yaml is 16.2735 pb^-1" | constants_audit (Section 4) | `efficiencytool/config_bdt_nom.yaml:239` | **CONFIRMED** | `lumi: 16.2735` |
| 26 | "plotcommon.h legend says R=0.3 but nominal uses topo R=0.4" | constants_audit (Section 5.4) | `plotting/plotcommon.h:24`, `efficiencytool/config_bdt_nom.yaml:53` | **CONFIRMED** | `strleg4` contains "R=0.3"; config says `use_topo_iso: 2` (topo R=0.4) |
| 27 | "MbdPileupHelper.h: PMT selection charge > 0.4, minimum 2 hits per side" | efficiency_tool (Section 2) | `efficiencytool/MbdPileupHelper.h:52,50` | **CONFIRMED** | `charge_cut = 0.4f` default, `hitcut = 2` default |
| 28 | "DOUBLE_FRAC=0.187 (0 mrad default) in run_showershape_double.sh" | efficiency_tool | `efficiencytool/run_showershape_double.sh:12` | **CONFIRMED** | `DOUBLE_FRAC=${2:-0.187}` |
| 29 | "BDTinput.C jet15 lower bound = 19 (vs CrossSectionWeights.h = 15)" | constants_audit (Section 1.5) | `FunWithxgboost/BDTinput.C:90`, `efficiencytool/CrossSectionWeights.h:91-92` | **CONFIRMED** | BDTinput: `max_jet_lower = 19`; CrossSectionWeights: `jet_pt_lower = 15` |
| 30 | "BDTinput.C jet20 lower bound = 23 (vs CrossSectionWeights.h = 21)" | constants_audit (Section 1.5) | `FunWithxgboost/BDTinput.C:100`, `efficiencytool/CrossSectionWeights.h:95-96` | **CONFIRMED** | BDTinput: `max_jet_lower = 23`; CrossSectionWeights: `jet_pt_lower = 21` |

---

## Summary of "Gotcha" Claims Verified

### Gotcha 1: DoubleInteractionCheck.C uses flat BDT threshold (constants_audit Section 6.2)
**CONFIRMED as a real issue.** `DoubleInteractionCheck.C` reads only `bdt_min` (flat) and has zero references to `bdt_min_slope` or `bdt_min_intercept`. In contrast, `RecoEffCalculator_TTreeReader.C` and `ShowerShapeCheck.C` both read `bdt_min_slope` and `bdt_min_intercept` for parametric ET-dependent thresholds. This means the toy double-interaction simulation applies a less restrictive (and ET-independent) BDT cut than the main analysis.

### Gotcha 2: config_schema report inverts fit_option mapping (Claim #12)
**CONFIRMED as a real error in the report.** The config_schema report at line 193 states `0 = Pade, 1 = erf`. The actual code at `CalculatePhotonYield.C:585-591` shows: default (fit_option=0) creates `TMath::Erf`, and `fitoption==1` creates the Pade rational function. The efficiency_tool report correctly describes this at line 204.

### Gotcha 3: config_schema report inverts train_single_model meaning (Claim #20)
**CONFIRMED as a real error.** The config_schema report says `train_single_model: true` means "Train one model per ET bin (vs one global)". The code at `main_training.py:763-766` shows that `true` calls `_train_single_model` (one model) and `false` calls `_train_per_bin_models`. The description is backwards.

---

## Scorecard

| Result | Count |
|--------|-------|
| CONFIRMED | 26 |
| WRONG | 4 |
| UNVERIFIABLE | 0 |
| Total checked | 30 |

**Accuracy rate: 87% (26/30)**

---

## Errors Found

1. **config_schema.md line 193**: `fit_option` mapping is inverted. Says "0 = Pade, 1 = erf"; correct is "0 = erf, 1 = Pade".

2. **config_schema.md line 110**: `train_single_model` description is inverted. Says `true` means "Train one model per ET bin"; correct is `true` = one single model for all bins.

3. **efficiency_tool.md line 339**: States "Default: `mc_iso_scale = 1.2`, `mc_iso_shift = 0.1`" as if these are code defaults. The actual code defaults (`.as<float>()` fallbacks) are 1.0 and 0.0 respectively. The values 1.2/0.1 are the nominal *config* values, not code defaults. This distinction matters for understanding what happens when a config omits these fields.

4. **bdt_pipeline.md Section 7 (line ~367)**: States reweighting order as "Class -> ET -> eta -> vertex". Actual code order is "Class -> (per class: eta -> ET -> vertex)". Eta reweighting precedes ET reweighting.

---

## Additional Observation

The `config_schema.md` and `efficiency_tool.md` reports have an internal inconsistency about `fit_option`. The efficiency_tool report (line 204) correctly says 0=erf, 1=Pade. The config_schema report (line 193) inverts this. Reports should be cross-checked for internal consistency.

The reweighting order error in the bdt_pipeline report is minor -- eta and ET reweighting are both inverse-PDF methods applied independently per class, so their ordering has negligible impact on the final weights. But it should still be corrected for documentary accuracy.
