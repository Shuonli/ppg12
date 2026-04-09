---
description: Verify that constants, filenames, features, and bins are in sync across pipeline codes
---

Check that critical values stay synchronized across the PPG12 analysis pipeline.

Arguments: $ARGUMENTS (optional focus area)
- `xsec` → cross-section weights only
- `features` → BDT feature lists only
- `ptbins` → pT bin definitions only
- `models` → BDT model names only
- `all` or (empty) → check everything

## Cross-section weights (`xsec`)

Grep for `photon5cross`, `photon10cross`, `photon20cross`, `jet10cross`, `jet15cross`, `jet20cross`, `jet30cross`, `jet50cross` across all `.C` files in `efficiencytool/` and `FunWithxgboost/`.

For each constant, compare the value across files. The authoritative values are in `RecoEffCalculator.C` and `MergeSim.C` — all other files must match. Report MATCH or MISMATCH with the exact values and file:line locations.

Known risk areas:
- `jet20cross`, `jet30cross`, `jet50cross` have historically diverged between main pipeline and QA macros (e.g., `ShowerShapeCheck.C`, `IsoROC_calculator.C`)

## BDT feature lists (`features`)

Compare the 25 features in `FunWithxgboost/config.yaml` under `data.features` against the feature vectors hardcoded in `apply_BDT.C` (the `x_list` arrays for each model, lines ~284-391). Each model uses a subset — verify every feature name in each subset exists in the training config.

Also check the NPB score features (lines ~451-477 of `apply_BDT.C`).

## pT bin definitions (`ptbins`)

Compare:
- `plotting/plotcommon.h`: `ptRanges[NptBins+1]` (lines 5-6)
- `efficiencytool/config_bdt_nom.yaml`: `analysis.pT_bins`
- `FunWithxgboost/config.yaml`: `binning.pt_edges`

Note: these serve different purposes (reco analysis vs plotting vs training) and are intentionally different, but flag any unexpected discrepancies.

## BDT model names (`models`)

Compare the `model_names` vector in `apply_BDT.C` (lines ~73-85) against:
- Actual model files in `FunWithxgboost/binned_models/model_*_single_tmva.root`
- `bdt_et_bin_models` entries in `efficiencytool/config_bdt_nom.yaml`

Flag any model referenced in configs that doesn't exist as a file.

## Additional checks (included in `all`)

- **Cluster node name**: Grep `CLUSTERINFO_CEMC` across all `.C` and `.yaml` files. Flag any that hardcode instead of reading from config.
- **Sample names**: Compare sample lists in `MergeSim.C` against actual directory structure at the MC input path.

## Output

Print a sync report with **MATCH/MISMATCH/WARNING** per check category and a final summary.

This is a **read-only** inspection. Do not modify any files.
