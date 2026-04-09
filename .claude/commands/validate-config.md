---
description: Validate an analysis config YAML for common errors (paths, bins, cuts, models)
---

Validate an analysis config YAML for common errors.

Arguments: $ARGUMENTS (config file path or variant shortname, e.g., `bdt_nom`)

## Steps

1. If `$ARGUMENTS` has no path separator, expand to `efficiencytool/config_bdt_$ARGUMENTS.yaml`.

2. Read the config file.

3. Run these validation checks and report **PASS/WARN/FAIL** for each:

### Structural checks
- `var_type` is set and non-empty
- `var_type` is unique: grep other `config_bdt_*.yaml` in the same directory for duplicates
- Required sections exist: `input`, `output`, `analysis`, `analysis.tight`, `analysis.non_tight`, `analysis.common`

### Path checks
- `input.data_file` glob matches at least one file
- `input.photon_jet_file_root_dir` directory exists
- `output.final_outfile` parent directory exists

### pT bin checks
- `analysis.pT_bins` is strictly ascending
- `analysis.pT_bins_truth` is strictly ascending
- `pT_bins_truth` range encompasses `pT_bins` range

### BDT model checks
- `bdt_et_bin_edges` has exactly N+1 entries for N entries in `bdt_et_bin_models`
- Each model in `bdt_et_bin_models` has a corresponding file in `FunWithxgboost/binned_models/`

### Cut consistency checks
- `non_tight.bdt_max` boundary does not overlap `tight.bdt_min` at reference ET=15 GeV
- `reco_iso_max_b + reco_iso_max_s * 10` produces a positive isolation cut at ET=10 GeV
- `reco_noniso_min_shift` > 0

### Physics sanity
- `eta_bins` is `[-0.7, 0.7]` (warn if different)
- `cone_size` is 3 or 4 (warn if other)
- `truth_iso_max` is 4.0 (warn if different)

4. Print a summary: total checks, passes, warnings, failures.

This is a **read-only** validation. Do not modify any files.
