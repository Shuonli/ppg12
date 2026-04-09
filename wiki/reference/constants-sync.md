# Constants Synchronization

This page tracks constants that must match across multiple files and flags current discrepancies.

## Cross-Section Weights

All values centralized in `efficiencytool/CrossSectionWeights.h`:

| Sample | Cross-section (pb) | Reference |
|--------|-------------------|-----------|
| photon5 | 146,359.3 | photon20cross |
| photon10 | 6,944.675 | photon20cross |
| photon20 | 130.4461 | (reference) |
| jet5 | 1.3878e+08 | jet50cross |
| jet8 | 1.15e+07 | jet50cross |
| jet10 | 3.997e+06 | jet50cross |
| jet12 | 1.4903e+06 | jet50cross |
| jet15 | 4.073e+05 | jet50cross |
| jet20 | 6.2623e+04 | jet50cross |
| jet30 | 2.5298e+03 | jet50cross |
| jet40 | 1.3553e+02 | jet50cross |
| jet50 | 7.3113 | (reference) |

**Files including CrossSectionWeights.h:** RecoEffCalculator_TTreeReader.C, CalculatePhotonYield.C, ShowerShapeCheck.C, DoubleInteractionCheck.C, MergeSim.C, BDTinput.C, and others. All use the centralized values.

### Discrepancies

**BDTinput.C jet reference denominator:**
BDTinput.C normalizes jet weights relative to `jet30cross` instead of `jet50cross`. The weights are only used for internal luminosity counting (not written to training files), so this is benign but violates the convention.

**FunWithTMVA/ (obsolete):**
Contains completely different cross-section values from an earlier Pythia tune. Not used in current pipeline.

**Jet pT window boundaries:**

| Sample | BDTinput.C lower | CrossSectionWeights.h lower |
|--------|------------------|-----------------------------|
| jet15 | **19** | **15** |
| jet20 | **23** | **21** |

BDTinput.C excludes valid truth jet pT ranges (15-19 GeV, 21-23 GeV) from training data.

## pT Bin Edges

| Source | N bins | Edges |
|--------|--------|-------|
| `plotcommon.h` | 12 | 8,10,12,14,16,18,20,22,24,26,**28,30,35** |
| `config_bdt_nom.yaml` | 12 | 8,10,12,14,16,18,20,22,24,26,**28,32,36** |
| `config_nom.yaml` (FunWithxgboost) | 10 | 8,10,12,...,24,26,35 |
| Showershape configs | 5 | 10,14,18,22,28,30 |

**The last two edges differ between plotcommon.h (30,35) and config_bdt_nom (32,36).** Plotting macros using `ptRanges` from plotcommon.h will not match results from config_bdt_nom.

**Python N_PT_BINS mismatch:**
`make_selection_report.py` and `make_comparison_report.py` hardcode `N_PT_BINS = 10`, missing the last 2 of 12 pT bins.

## Luminosity

| Value (pb^-1) | Source | Period |
|---------------|--------|--------|
| 16.2735 | config_bdt_nom.yaml | 1.5 mrad nominal |
| 16.6 | plotcommon.h, CLAUDE.md, conf note | Rounded label |
| 16.8588 | config_showershape_1p5rad.yaml | 1.5 mrad (different calc?) |
| 16.9 | PPG12-analysis-note | Analysis note |
| 32.6574 | config_bdt_0rad.yaml | 0 mrad period |
| 49.562 | config_bdt_all.yaml | All runs combined |

Five distinct values appear for the 1.5 mrad period. The actual YAML config uses 16.2735, the analysis note says 16.9 with uncertainties, and plot labels say 16.6.

## Isolation Parameters

| Set | intercept | slope | Used in |
|-----|-----------|-------|---------|
| Nominal BDT | 0.502095 | 0.0433036 | config_bdt_nom.yaml |
| Showershape | 0.453194 | 0.0360234 | config_showershape_*.yaml |
| FunWithxgboost | 1.08128 | 0.0299107 | config_nom.yaml |
| Conf note | 1.08 | 0.03 | ppg12_conf_note/main.tex |

The conf note documents the older FunWithxgboost parameters, not the current nominal.

**plotcommon.h vs code:** Legend says "R=0.3, < 4 GeV" but nominal reco isolation uses topo R=0.4 with parametric cut.

## BDT Thresholds

| Set | slope | intercept | Used by |
|-----|-------|-----------|---------|
| Nominal BDT | -0.015 | 0.80 | RecoEffCalculator (parametric) |
| Showershape | -0.02 | 0.86 | ShowerShapeCheck configs |

**DoubleInteractionCheck.C reads only flat `bdt_min`** (no slope/intercept). It applies a less restrictive, ET-independent threshold compared to the parametric cut in the main analysis.

## Geometry Constants (consistent)

| Constant | Value | Status |
|----------|-------|--------|
| EMCal radius | 93.5 cm | Consistent everywhere |
| Eta range | [-0.7, 0.7] | Consistent everywhere (see [Eta Edge Migration](../concepts/eta-edge-migration.md) for boundary effects) |
| Truth iso cut | 4.0 GeV | Consistent everywhere |
| nsimevents | 1e7 | Consistent everywhere |

## TIME_SAMPLE_NS

| Value | Files |
|-------|-------|
| 17.6 ns | BDTinput.C, RecoEffCalc, ShowerShapeCheck, DoubleInteractionCheck, + 5 others |
| **16.67 ns** | **time_energy_corr.C only** |

The 16.67 ns value (= 1/60 MHz) appears to be a bug. The correct measured value is 17.6 ns.

## Priority Ranking of Discrepancies

1. **pT bin mismatch** (plotcommon.h vs config_bdt_nom) -- plotting macros bin data incorrectly
2. **N_PT_BINS=10 in Python** -- reports miss last 2 bins
3. **Conf note isolation parameters** -- documents old scheme
4. **plotcommon.h "R=0.3" legend** -- misleading on reco plots
5. **Luminosity scatter** -- 5 values for same period
6. **DoubleInteractionCheck.C flat BDT** -- less restrictive than main analysis
7. **BDTinput.C jet pT boundaries** -- training data excludes valid events
