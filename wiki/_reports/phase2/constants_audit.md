# Constants Audit Report -- PPG12 Isolated Photon Analysis

**Generated**: 2026-04-08
**Scope**: All shared constants across the PPG12 pipeline -- cross-section weights, pT bins, eta range, luminosity, isolation parameters, BDT thresholds, geometry, trigger efficiency, truth isolation.

---

## Executive Summary

The codebase has a centralized `CrossSectionWeights.h` header that most efficiency-tool macros now include. However, several inconsistencies remain:

1. **BDTinput.C uses jet30cross as jet reference denominator** instead of jet50cross -- different normalization scheme from the rest of the pipeline (benign for its purpose, but a maintainability risk).
2. **Luminosity is inconsistent across documents**: 16.6 pb^-1 (conf note, plotcommon.h), 16.9 pb^-1 (analysis note), 16.2735 pb^-1 (config_bdt_nom.yaml), 16.8588 pb^-1 (1.5 mrad configs).
3. **plotcommon.h legend says R=0.3** but the nominal reco isolation cone is topo R=0.4 (`use_topo_iso: 2`).
4. **pT bins differ between plotcommon.h and config_bdt_nom.yaml**: plotcommon.h has 12 bins ending at 35, config_bdt_nom has 12 bins ending at 36.
5. **Python report scripts hardcode N_PT_BINS=10**, inconsistent with both plotcommon.h (12) and config_bdt_nom (12).
6. **FunWithTMVA/ has completely obsolete cross-section values** (factor-of-~1 different from current values).
7. **TIME_SAMPLE_NS differs**: 17.6 ns in most files but 16.67 ns in `time_energy_corr.C`.
8. **Showershape configs use different BDT thresholds** (-0.02/0.86) than nominal BDT configs (-0.015/0.80).
9. **Conf note isolation parameters (1.08 + 0.03*ET)** match the older FunWithxgboost/config_nom.yaml, not the current config_bdt_nom.yaml (0.502 + 0.043*ET).
10. **Jet pT boundaries differ between BDTinput.C and CrossSectionWeights.h** for jet15 (19 vs 15) and jet20 (23 vs 21).

---

## 1. Cross-Section Weights

### 1.1 Central Values (CrossSectionWeights.h)

All values from `efficiencytool/CrossSectionWeights.h`:

| Sample    | Cross-section (pb) | Reference denominator |
|-----------|--------------------|-----------------------|
| photon5   | 146359.3           | photon20cross         |
| photon10  | 6944.675           | photon20cross         |
| photon20  | 130.4461           | (reference)           |
| jet5      | 1.3878e+08         | jet50cross            |
| jet8      | 1.15e+07           | jet50cross            |
| jet10     | 3.997e+06          | jet50cross            |
| jet12     | 1.4903e+06         | jet50cross            |
| jet15     | 4.073e+05          | jet50cross            |
| jet20     | 6.2623e+04         | jet50cross            |
| jet30     | 2.5298e+03         | jet50cross            |
| jet40     | 1.3553e+02         | jet50cross            |
| jet50     | 7.3113             | (reference)           |

### 1.2 Files Using CrossSectionWeights.h (Consistent)

All of these `#include "CrossSectionWeights.h"` and use the centralized constants:

- `efficiencytool/RecoEffCalculator_TTreeReader.C` -- uses `GetSampleConfig()` for weights
- `efficiencytool/CalculatePhotonYield.C` -- uses `photon20cross`, `jet50cross`
- `efficiencytool/ShowerShapeCheck.C` -- uses `GetSampleConfig()`
- `efficiencytool/DoubleInteractionCheck.C` -- uses `GetSampleConfig()`
- `efficiencytool/RecoEffCalculator.C` -- inline weight assignments match header values
- `efficiencytool/SaturationStudy.C` -- inline weight assignments match header values
- `efficiencytool/BDTScoreVsET.C` -- inline weight assignments match header values
- `efficiencytool/IsoROC_calculator.C` -- inline weight assignments match header values
- `efficiencytool/TruthJetEff.C` -- includes header
- `efficiencytool/PlotSpectra.C` -- includes header
- `efficiencytool/plot_cluster_time.C` -- includes header
- `efficiencytool/AnalyzeTruthPhotonTowers.C` -- includes header
- `efficiencytool/TruthJetInETWindow.C` -- includes header
- `plotting/plot_SB.C` -- includes `../efficiencytool/CrossSectionWeights.h`
- `FunWithxgboost/BDTinput.C` -- includes `../efficiencytool/CrossSectionWeights.h`

### 1.3 BDTinput.C Jet Normalization Mismatch

`FunWithxgboost/BDTinput.C` includes `CrossSectionWeights.h` and reads the correct absolute cross-section values, but **normalizes jet weights relative to jet30cross instead of jet50cross**:

| Sample | BDTinput.C weight          | CrossSectionWeights.h weight |
|--------|----------------------------|------------------------------|
| jet10  | `jet10cross / jet30cross`  | `jet10cross / jet50cross`    |
| jet15  | `jet15cross / jet30cross`  | `jet15cross / jet50cross`    |
| jet20  | `jet20cross / jet30cross`  | `jet20cross / jet50cross`    |
| jet30  | `1.0` (reference)          | `jet30cross / jet50cross`    |
| jet50  | `jet50cross / jet30cross`  | `1.0` (reference)            |

**Impact**: BDTinput.C uses these weights only for internal per-run luminosity counting, not for training sample weights (the weights are not written to the output txt files). The relative ratios between jet samples are preserved regardless of the reference denominator, so the jet10:jet15:jet20:jet30 ratios are correct. However, the jet30:jet50 ratio is inverted in interpretation. This is a **maintainability concern** rather than a physics bug -- the code works correctly for its purpose (BDT feature extraction does not depend on absolute cross-section weighting), but it violates the convention established by CrossSectionWeights.h.

### 1.4 Legacy FunWithTMVA/ (Obsolete, Different Values)

Files `FunWithTMVA/TMVA_train*.C` define **completely different cross-section constants**:

| Sample    | FunWithTMVA value                    | Computed (pb)  | CrossSectionWeights.h (pb) |
|-----------|--------------------------------------|----------------|----------------------------|
| photon5   | `2.017e+08 * 0.000442571`            | 89,284         | 146,359                    |
| photon10  | `3.690e+07 * 0.000181474`            | 6,696          | 6,945                      |
| photon20  | `1.571e+05 * 0.000673448`            | 105.8          | 130.4                      |
| jet10     | `3.646e-6 * 1e12`                    | 3,646,000      | 3,997,000                  |
| jet15     | `36864930.0 * 0.011059973`           | 407,632        | 407,300                    |
| jet20     | `1392140.9 * 0.042`                  | 58,470         | 62,623                     |
| jet30     | `2.505e-9 * 1e12`                    | 2,505          | 2,530                      |

These are **obsolete** values from an earlier Pythia tune or cross-section calculation. FunWithTMVA/ is a legacy directory not used in the current pipeline. **No physics impact on current results**, but could confuse anyone reusing that code.

### 1.5 Jet pT Window Boundaries

BDTinput.C uses different jet pT window boundaries than CrossSectionWeights.h:

| Sample | BDTinput.C lower | CrossSectionWeights.h lower | BDTinput.C upper (commented) | CrossSectionWeights.h upper |
|--------|------------------|-----------------------------|-----------------------------|-----------------------------|
| jet10  | 10               | 10                          | 19 (commented)              | 15                          |
| jet15  | **19**           | **15**                      | 23 (commented)              | 21                          |
| jet20  | **23**           | **21**                      | 30 (commented)              | 32                          |
| jet30  | 30               | 32                          | 100                         | 42                          |
| jet50  | 50               | 52                          | 70                          | 100                         |

The BDTinput.C boundaries do not match CrossSectionWeights.h. Note that BDTinput.C has the upper bounds commented out, so in practice only the lower bound is enforced. The mismatch in lower bounds (jet15: 19 vs 15; jet20: 23 vs 21) means BDTinput.C excludes a portion of each sample's valid truth jet pT range.

---

## 2. pT Bin Edges

### 2.1 plotcommon.h (Plotting Reference)

```cpp
const int NptBins = 12;
const float ptRanges[NptBins + 1] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35};
```

### 2.2 Nominal BDT Config (config_bdt_nom.yaml)

```yaml
pT_bins: [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
pT_bins_truth: [7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45]
```

### 2.3 Comparison

| Source | N bins | Edges | Discrepancy |
|--------|--------|-------|-------------|
| plotcommon.h | 12 | 8,10,12,14,16,18,20,22,24,26,**28,30,35** | -- |
| config_bdt_nom.yaml | 12 | 8,10,12,14,16,18,20,22,24,26,**28,32,36** | Last 2 edges differ |
| FunWithxgboost/config_nom.yaml | 10 | 8,10,12,14,16,18,20,22,24,26,35 | Fewer bins, different top |
| config_bdt_nom.yaml (truth) | 14 | 7,...,28,32,36,45 | Extended for unfolding |
| Showershape configs | 5 | 10,14,18,22,28,30 | Coarser binning |
| config/config.yaml | 7 | 8,10,12,14,16,20,30,35 | Legacy, different |
| Most config_bdt_*.yaml | 12 | 8,10,12,14,16,18,20,22,24,26,28,32,36 | Match nominal |
| Backup configs | 10 | 8,10,12,14,16,18,20,22,24,26,35 | Older binning |

**The last two bin edges differ between plotcommon.h (28,30,35) and the current nominal config (28,32,36)**. Plotting macros that use `NptBins` and `ptRanges` from plotcommon.h will not match the binning of results produced with config_bdt_nom.yaml.

### 2.4 Python N_PT_BINS Mismatch

| File | Value | plotcommon.h NptBins |
|------|-------|---------------------|
| `plotting/make_selection_report.py` | `N_PT_BINS = 10` | 12 |
| `plotting/make_comparison_report.py` | `N_PT_BINS = 10` | 12 |

Both Python scripts hardcode `N_PT_BINS = 10` with a comment saying "matches NptBins in plotcommon.h", but plotcommon.h actually has `NptBins = 12` and the nominal config also has 12 bins. This means the Python report generators will only process the first 10 of 12 pT bins, **missing the last 2 bins**.

---

## 3. Eta Range

**Consistent across the pipeline.** All sources use `|eta| < 0.7`:

| Source | Value |
|--------|-------|
| All config YAML `eta_bins` | `[-0.7, 0.7]` |
| plotcommon.h `strleg3` | `"|#it{#eta^{#gamma}}| < 0.7"` |
| `CalculatePhotonYield.C` solid_angle | `2 * M_PI * 0.7 * 2` |
| Analysis note | `$|\eta^{\gamma}| < 0.7$` |
| Conf note | `$|\eta^{\gamma}| < 0.7$` |
| Hardcoded in BDTinput.C | `fabs(eta) < 0.7` |
| `investigate_vertex_shift.C` | `eta_min = -0.7, eta_max = 0.7` |

**No discrepancies found.**

---

## 4. Luminosity

Multiple luminosity values appear across the codebase:

| Value (pb^-1) | Where Used | Period |
|---------------|------------|--------|
| 16.2735 | config_bdt_nom.yaml and most config_bdt_*.yaml | 1.5 mrad nominal |
| 16.6 | plotcommon.h (`strleg2_1`, `strleg5`), CLAUDE.md | Rounded for plots/labels |
| 16.8588 | config_showershape_1p5rad.yaml | 1.5 mrad (different calc?) |
| 16.9 | PPG12-analysis-note (introduction.tex, selection.tex, analysis.tex) | Analysis note |
| 32.6574 | config_showershape_0radt2.yaml, config_showershape_0rad.yaml, config_bdt_0rad.yaml | 0 mrad period |
| 49.562 | config_showershape_nom.yaml, config_bdt_all.yaml, config_backup/config_bdt_nom.yaml | All runs combined |
| $16.6~\mathrm{pb}^{-1}$ | ppg12_conf_note/main.tex (lines 60, 260) | Conf note |
| $16.9^{+1.4}_{-1.2}~\mathrm{pb}^{-1}$ | PPG12-analysis-note/selection.tex:135 | Analysis note with uncertainty |
| 15.2036 | plotting/plot_final_backup250327_v1.C (commented-out) | Legacy |
| 49.562 | plotting/plot_final_backup.C | All-run plot |

**Inconsistencies**:

1. The conf note states **16.6 pb^-1** but the analysis note states **16.9 pb^-1** (with asymmetric uncertainties).
2. The actual YAML config for the nominal 1.5 mrad analysis uses **16.2735 pb^-1**, which differs from both documents.
3. `plotcommon.h` hardcodes **16.6 pb^-1** in legend strings.
4. The 1.5 mrad showershape config uses **16.8588 pb^-1**, a fourth distinct value.
5. `CLAUDE.md` states "16.6 pb^-1 (Run 24)" and `yaml-config.md` states "16.8588 pb^-1" for the 1.5 mrad period.

These are not all directly comparable (some include trigger/vertex corrections, some don't), but the lack of a single documented value chain is a source of confusion.

---

## 5. Isolation Parameters

### 5.1 Reco Isolation Cone

The pipeline has evolved from simple cone isolation to topological cluster isolation:

| Parameter | Value | Files |
|-----------|-------|-------|
| `cone_size` | 3 (meaning R=0.3) | All YAML configs |
| `use_topo_iso` | 2 (topo R=0.4, nominal) | Most current configs |
| `use_topo_iso` | 0 (standard cone R=0.3) | config_bdt_none.yaml only |

When `use_topo_iso=2`, the code uses `cluster_iso_topo_04` (R=0.4) regardless of `cone_size`. The `cone_size: 3` field is effectively ignored for the nominal analysis.

### 5.2 Reco Isolation Cut Parameters

Three distinct sets of (intercept, slope) appear:

| Set | intercept | slope | Used in |
|-----|-----------|-------|---------|
| **Nominal BDT** | 0.502095 | 0.0433036 | config_bdt_nom.yaml and all config_bdt_*.yaml |
| **Showershape** | 0.453194 | 0.0360234 | All config_showershape_*.yaml |
| **FunWithxgboost** | 1.08128 | 0.0299107 | FunWithxgboost/config_nom.yaml, config_nom_split.yaml |

The showershape configs using different isolation parameters than the nominal BDT configs is expected (they are separate studies). The FunWithxgboost configs have the oldest values.

### 5.3 Conf Note Isolation Description vs Code

The conf note (ppg12_conf_note/main.tex:103) states:
> "An ET-dependent isoET threshold is chosen to maintain an 80% isolation efficiency, isoET < 1.08 GeV + 0.03 * ET_gamma"

This matches `FunWithxgboost/config_nom.yaml` (reco_iso_max_b=1.08128, reco_iso_max_s=0.0299107), which is the **older** configuration. The current nominal is `config_bdt_nom.yaml` with (0.502095, 0.0433036).

Additionally, the conf note describes the isolation as summing towers with E > 60 MeV within R=0.3, while the actual nominal analysis uses topological cluster isolation with R=0.4 (`use_topo_iso: 2`).

### 5.4 plotcommon.h Legend Mismatch

```cpp
string strleg4 = "#it{E}_{T}^{iso, #kern[-0.2]{#it{R}=0.3}}< 4 GeV";
```

This legend says R=0.3 and a flat 4 GeV cut. The actual nominal analysis uses:
- Topo cluster R=0.4 (not R=0.3)
- Parametric cut: 0.502 + 0.0433*ET (not flat 4 GeV)

The `strleg4` string is a truth-level description (truth isolation is indeed R=0.3, < 4 GeV), but it could be misleading if used on reco-level plots.

---

## 6. BDT Threshold Parameters

### 6.1 Tight BDT Min

Two distinct parameter sets in active configs:

| Set | slope | intercept | Used in |
|-----|-------|-----------|---------|
| **Nominal BDT** | -0.015 | 0.80 | config_bdt_nom.yaml and all config_bdt_*.yaml |
| **Showershape** | -0.02 | 0.86 | All config_showershape_*.yaml |

Both tight and non-tight blocks use the same slope/intercept within a given config. This is by design (the parametric threshold is the same form, the bdt_min/bdt_max flat bounds define the acceptance window).

### 6.2 Legacy Flat BDT Thresholds

Some macros read only a flat `bdt_min` field (no slope/intercept):
- `DoubleInteractionCheck.C`: `tight_bdt_min = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0.0)`
- `ShowerShapeCheck.C` (embed_test version): same
- `NPB_PurityStudy.C`: same

These macros use a **flat threshold** (no ET dependence). When used with showershape configs that have `bdt_min: 0.7` (tight) and `bdt_min_slope: -0.02`, `bdt_min_intercept: 0.86`, only the flat `bdt_min` is read by these macros. The parametric version (slope*ET + intercept) is used only by `RecoEffCalculator_TTreeReader.C`.

This means `DoubleInteractionCheck.C` and `ShowerShapeCheck.C` apply a **less restrictive** BDT threshold than the main analysis code for the same config.

---

## 7. Geometry Constants

### 7.1 EMCal Radius

| File | Value | Consistent? |
|------|-------|-------------|
| `DoubleInteractionCheck.C` | `cemc_radius_cm = 93.5f` | Yes |
| `investigate_vertex_shift.C` | `R_CEMC = 93.5` | Yes |
| Analysis note | `R_EMCal = 93.5 cm` | Yes |
| Reports | `R_EMCal = 93.5 cm` | Yes |

**Consistent.**

### 7.2 TIME_SAMPLE_NS

| File | Value |
|------|-------|
| `BDTinput.C` | 17.6 ns |
| `RecoEffCalculator_TTreeReader.C` | 17.6 ns |
| `ShowerShapeCheck.C` | 17.6 ns |
| `DoubleInteractionCheck.C` | 17.6 ns |
| `plot_cluster_time.C` | 17.6 ns |
| `plot_cluster_jet_time.C` | 17.6 ns |
| `plot_cluster_mbd_time_eta.C` | 17.6 ns |
| `AnalyzeTruthPhotonTowers.C` | 17.6 ns |
| `embed_test/ShowerShapeCheck.C` | 17.6 ns |
| **`time_energy_corr.C`** | **16.67 ns** |

`time_energy_corr.C` uses 16.67 ns instead of 17.6 ns. This ~5% difference would affect any timing calibration derived from that file. The value 16.67 ns = 1/(60 MHz), while 17.6 ns is the actual measured sample period. This appears to be a bug in `time_energy_corr.C`.

---

## 8. Trigger Efficiency / MBD Correction

### 8.1 CalculatePhotonYield.C Constants

```cpp
float mbdcorr = 25.2/42 / 0.57;     // = 1.053
float jetevents = 0.3555 * 1E7;      // effective jet events (trigger fraction)
```

These appear only in `CalculatePhotonYield.C` and are not cross-referenced elsewhere. The `0.57` factor and `25.2/42` factor are not documented in config files -- they are hardcoded physics constants.

### 8.2 MBD Efficiency Scale

From configs:
- `mbd_eff_scale: 0.05` (config_backup/config_mbdeffup.yaml)
- `mbd_eff_scale: -0.05` (config_backup/config_mbdeffdown.yaml)
- Default: 0.0 (most configs)

This is a systematic variation parameter, not a hardcoded constant.

---

## 9. Truth Isolation Cut

**Consistent.** All sources use 4.0 GeV:

| Source | Value |
|--------|-------|
| All YAML configs `truth_iso_max` | 4.0 |
| Analysis note | "isoETtruth < 4 GeV" |
| Conf note | "less than 4 GeV" |
| `investigate_deltaR.C` (hardcoded) | `truth_iso_max = 4.0` |
| plotcommon.h `strleg4` | "< 4 GeV" |

**No discrepancies.**

---

## 10. Summary Table

| Constant | Canonical Value | Discrepant File(s) | Nature of Discrepancy |
|----------|----------------|--------------------|-----------------------|
| **Cross-section weights** | CrossSectionWeights.h values | FunWithTMVA/*.C | Completely different values (obsolete) |
| **Jet weight reference** | jet50cross (denominator) | BDTinput.C | Uses jet30cross as denominator |
| **Jet pT boundaries** | CrossSectionWeights.h | BDTinput.C jet15 (19 vs 15), jet20 (23 vs 21) | Different lower bounds |
| **pT bin edges** | config_bdt_nom: [8..28,32,36] | plotcommon.h: [8..28,30,35] | Last 2 edges differ (30,35 vs 32,36) |
| **N_PT_BINS** | 12 (plotcommon.h, config_bdt_nom) | make_selection_report.py, make_comparison_report.py: 10 | Python reports miss last 2 bins |
| **Luminosity** | 16.2735 pb^-1 (config_bdt_nom) | plotcommon.h: 16.6; conf note: 16.6; analysis note: 16.9; showershape 1p5rad: 16.8588 | 4+ different values |
| **Reco isolation cone** | Topo R=0.4 (use_topo_iso=2) | plotcommon.h strleg4: "R=0.3"; conf note: "R=0.3" | Label/text says 0.3, code uses topo 0.4 |
| **Reco iso parameters** | 0.502/0.0433 (config_bdt_nom) | conf note: "1.08 + 0.03*ET" (matches older config_nom) | Conf note documents old parameters |
| **BDT threshold** | -0.015/0.80 (config_bdt_nom) | showershape configs: -0.02/0.86 | Different study, but same config schema |
| **BDT threshold form** | Parametric (slope*ET+intercept) | DoubleInteractionCheck.C, ShowerShapeCheck.C: flat bdt_min only | Macros ignore slope/intercept |
| **TIME_SAMPLE_NS** | 17.6 ns | time_energy_corr.C: 16.67 ns | ~5% difference |
| **Eta range** | [-0.7, 0.7] | (none) | Consistent everywhere |
| **Truth iso cut** | 4.0 GeV | (none) | Consistent everywhere |
| **EMCal radius** | 93.5 cm | (none) | Consistent everywhere |
| **nsimevents** | 1e7 | (none) | Consistent (CrossSectionWeights.h + CalculatePhotonYield.C) |

### Priority Ranking

Items that could affect physics results if not addressed:

1. **pT bin mismatch** (plotcommon.h vs config_bdt_nom) -- plotting macros using `ptRanges` will bin data incorrectly for results from config_bdt_nom
2. **N_PT_BINS=10 in Python** -- report generation misses last 2 pT bins
3. **Conf note isolation description** -- documents old isolation scheme (R=0.3, 1.08+0.03*ET) instead of current (topo R=0.4, 0.502+0.043*ET)
4. **plotcommon.h legend "R=0.3"** -- misleading on reco-level isolation plots
5. **Luminosity scatter** -- 16.2735 vs 16.6 vs 16.9 creates confusion about which is authoritative
6. **DoubleInteractionCheck.C flat BDT** -- applies less restrictive BDT cut than the parametric threshold in the main analysis
7. **TIME_SAMPLE_NS in time_energy_corr.C** -- if this file is used for any calibration, the 5% timing difference matters
8. **BDTinput.C jet pT boundaries** -- training data may exclude valid events in the 15-19 GeV and 21-23 GeV truth jet pT ranges
