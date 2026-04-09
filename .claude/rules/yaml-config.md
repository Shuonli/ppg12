---
globs: "**/*.yaml"
---

# YAML Config Conventions (PPG12)

- `var_type` must be unique across configs — it becomes the output filename suffix (`Photon_final_{var_type}.root`). Duplicate var_types silently overwrite results.
- `eta_bins` default: `[-0.7, 0.7]` (barrel EMCal coverage)
- `pT_bins` must be in ascending order; `pT_bins_truth` extends beyond reco bins on both sides for unfolding overflow
- Parametric isolation: `reco_iso_max = reco_iso_max_b + reco_iso_max_s * ET` — two fields, not one flat value
- Parametric BDT: `tight_bdt_min_intercept + tight_bdt_min_slope * ET`
- `bdt_et_bin_edges` needs N+1 entries for N models in `bdt_et_bin_models`
- `cone_size: 3` means R = 0.3 (integer encoding)
- All file paths must be absolute GPFS paths (`/sphenix/user/...`)
- `lumi` value must match the `run_min`/`run_max` range — cross-check with `lumi/60cmLumi_fromJoey.list`
- Run periods by crossing angle (split at run 51274):
  - **0 mrad**: `run_min: 47289`, `run_max: 51274`, `lumi: 32.6574` pb⁻¹
  - **1.5 mrad** (nominal): `run_min: 51274`, `run_max: 54000`, `lumi: 16.8588` pb⁻¹
  - **All runs**: `run_min: 47289`, `run_max: 54000`, `lumi: 49.562` pb⁻¹
  - `run_min: -1` / `run_max: -1` means no run range filter (accepts all runs)
- Use `make_bdt_variations.py` to generate systematic variation configs — each variant dict must have `syst_type` and `syst_role` fields (or None for cross-checks)
