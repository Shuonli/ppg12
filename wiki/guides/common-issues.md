# Common Issues and Gotchas

Known bugs, discrepancies, and pitfalls in the PPG12 codebase.

## Critical Issues

### ~~pT Bin Mismatch Between plotcommon.h and Nominal Config~~ (FIXED 2026-04-08)

**Problem:** `plotcommon.h` had `ptRanges` ending in `30, 35` but config uses `32, 36`.

**Status:** Fixed. `plotcommon.h` now matches config: `{8,10,12,14,16,18,20,22,24,26,28,32,36}`, `pTmax=36`.

> **Note:** Some `config_showershape_*.yaml` files still have `pT_bins_truth` ending in `28, 30, 35` — these are separate configs for the showershape study and may be intentional.

### ~~Python Report Scripts Miss Last 2 pT Bins~~ (FIXED 2026-04-08)

**Problem:** `make_selection_report.py` and `make_comparison_report.py` hardcoded `N_PT_BINS = 10`. The actual count is 12.

**Status:** Fixed. Both files now have `N_PT_BINS = 12`.

### ~~DoubleInteractionCheck.C Uses Flat BDT Threshold~~ (FIXED 2026-04-08)

**Problem:** `DoubleInteractionCheck.C` used flat `tight_bdt_min` instead of parametric `slope * ET + intercept`.

**Status:** Fixed. Both tight and non-tight thresholds are now parametric via a `classifyTightNonTight` lambda called from all 3 classification sites (original, shifted, smeared). Backward compatible — slope defaults to 0.

> **Note:** `ShowerShapeCheck.C` still uses a flat non-tight BDT upper bound (line ~1787). This is a remaining inconsistency.

## Configuration Issues

### Duplicate var_type Causes Silent Overwrite

**Problem:** If two configs have the same `var_type`, their output files will silently overwrite each other.

**Prevention:** Always verify `var_type` is unique across all configs. The `make_bdt_variations.py` script generates unique names automatically.

### bdt_et_bin_edges vs bdt_et_bin_models Off-by-One

**Problem:** `bdt_et_bin_edges` needs N+1 entries for N entries in `bdt_et_bin_models`. Getting this wrong causes undefined behavior.

**Check:** `len(bdt_et_bin_edges) == len(bdt_et_bin_models) + 1`

### mc_iso_scale/mc_iso_shift Defaults vs Config Values

**Problem:** The code defaults for mc_iso_scale and mc_iso_shift are 1.0 and 0.0 respectively. The nominal config sets them to 1.2 and 0.1. If a config omits these fields, the code defaults kick in -- which is different from the nominal analysis.

**Prevention:** Always explicitly set these fields in every analysis config.

### fit_option Mapping

**Corrected:** `fit_option=0` is the error function (erf), `fit_option=1` is the Pade rational function. Earlier documentation had this inverted.

### train_single_model Meaning

**Corrected:** `train_single_model: true` means train ONE model for all ET bins (not per-bin). `false` means train separate models per ET bin.

## Data Processing Issues

### Two-Pass Architecture Order (main cross-section pipeline)

**Problem:** In `RecoEffCalculator_TTreeReader.C`, the vertex scan (Pass 1) MUST complete before the full analysis (Pass 2) begins. Running them out of order produces incorrect vertex reweighting.

**Solution:** Use `oneforall_tree.sh` which enforces the correct order.

**Note:** The showershape DI pipeline no longer uses this two-pass architecture. It runs single-pass with truth-vertex reweighting (`submit_showershape_di.sub` + `TruthVertexReweightLoader.h`); the legacy two-pass reco-vertex script is kept as `run_showershape_double_reco_legacy.sh`.

### MergeSim Expected Samples

**Problem:** MergeSim.C expects exactly 6 jet samples (jet5, jet8, jet12, jet20, jet30, jet40). If you add or remove jet samples from the pipeline, you must also update MergeSim.C.

**Note:** jet10, jet15, and jet50 are used for BDT training but NOT included in MergeSim.

### hadd Must Be Run Manually

**Problem:** After sim condor jobs complete, you must manually run `hadd_combined.sh` to merge per-job outputs into `combined.root`.

### ~~BDTinput.C et4_on Copy-Paste Bug~~ (FIXED 2026-04-08)

**Problem:** Line 194 read `configYaml["analysis"]["et3_on"]` instead of `"et4_on"`, causing et4_on to shadow et3_on's value.

**Status:** Fixed. Now correctly reads `"et4_on"`.

### ~~npb_score_cut Dual YAML Paths~~ (FIXED 2026-04-08)

**Problem:** `DoubleInteractionCheck.C`, `NPB_PurityStudy.C`, and `embed_test/ShowerShapeCheck.C` read `analysis.npb_score_cut` (flat path, doesn't exist in configs, fell back to default 0.5). Main analysis reads `analysis.common.npb_score_cut`.

**Status:** Fixed. All files now read `["analysis"]["common"]["npb_score_cut"]`.

### Division-by-Zero Protection Inconsistency

**Problem:** BDTinput.C uses `safe_div` (checks denominator != 0 and isfinite), while apply_BDT.C checks `numerator > 0`. Functionally equivalent for positive energies, but could theoretically differ for edge cases.

## Constant Discrepancies

### Luminosity Values (5 Different Numbers for Same Period)

The 1.5 mrad period has these luminosity values:
- 16.2735 pb^-1 (config_bdt_nom.yaml)
- 16.6 pb^-1 (plotcommon.h, conf note)
- 16.8588 pb^-1 (config_showershape_1p5rad.yaml)
- 16.9 pb^-1 (analysis note)
- The config value is what actually enters the cross-section calculation.

### Isolation Legend Text

plotcommon.h says "R=0.3, < 4 GeV" but the nominal analysis uses topo R=0.4 with parametric cut. The legend text describes truth-level isolation.

### BDTinput.C Jet pT Boundaries

BDTinput.C uses different lower bounds than CrossSectionWeights.h for jet15 (19 vs 15 GeV) and jet20 (23 vs 21 GeV). This excludes some valid truth jet pT events from training data.

### TIME_SAMPLE_NS

`time_energy_corr.C` uses 16.67 ns (= 1/60 MHz) instead of the correct 17.6 ns used everywhere else. This is likely a bug.

## Build Issues

### YAML-CPP Loading

All C++ macros load yaml-cpp via:
```cpp
gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
```
This is a hardcoded path to a user installation.

### CNN Probability Always -1

`cluster_CNN_prob` is always -1 because the ONNX runtime is disabled in CaloAna24.cc. Do not rely on this branch.

## Workflow Tips

### Checking Condor Job Status

```bash
condor_q  # See running jobs
condor_q -analyze  # Diagnose held jobs
```

### Verifying Pipeline Completeness

```bash
cd efficiencytool/results
# Check all expected output files exist for a given var_type
for f in MC_efficiency_{bdt_nom,jet_bdt_nom} data_histo_bdt_nom MC_response_bdt_nom Photon_final_bdt_nom; do
  ls ${f}.root 2>/dev/null || echo "MISSING: ${f}.root"
done
```

### Debugging ABCD Issues

If purity looks wrong, check:
1. R factor (`h_R` in Photon_final): should be close to 1
2. Leakage fractions (`h_leak_B/C/D`): should be small (< 0.1)
3. Signal region (A) must have more events than sideband regions
4. Non-tight region boundaries should be adjacent to (not overlapping with) tight boundaries
