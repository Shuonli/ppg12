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
- `lumi` value must match the `run_min`/`run_max` range — cross-check with `lumi/60cmLumi_fromJoey.list` (60cm fiducial) or `lumi/allzLumi_fromJoey.list` (beam-delivered)
- Run periods by crossing angle (split at run 51274), 60cm-fiducial lumi:
  - **0 mrad**: `run_min: 47289`, `run_max: 51274`, `lumi: 32.6574` pb⁻¹
  - **1.5 mrad** (nominal): `run_min: 51274`, `run_max: 54000`, `lumi: 16.2735` pb⁻¹
  - **All runs**: `run_min: 47289`, `run_max: 54000`, `lumi: 48.9309` pb⁻¹ (per-period sum)
  - `run_min: -1` / `run_max: -1` means no run range filter (accepts all runs)
- Beam-delivered (allz) lumi for the all-z fiducial cross-check:
  - **0 mrad allz**: `lumi: 47.2076` pb⁻¹; **1.5 mrad allz**: `lumi: 17.1642` pb⁻¹; **All allz**: `lumi: 64.3718` pb⁻¹
- `lumi_target`: optional field. Default = `lumi` → no scaling. Set to `sum(L_periods)` in merge-feeder configs (`config_bdt_0rad.yaml`, `config_bdt_nom.yaml`, `config_bdt_allz_{0rad,1p5mrad}.yaml`) so per-event `lumi_weight = lumi/lumi_target` pre-scales MC fills and a plain hadd across periods reproduces all-range MC. See `wiki/pipeline/full-run-range.md`.
- `vertex_cut_truth`: optional field. Default = `vertex_cut`. Decouples the truth-vertex cut on the MBD-eff denominator from the analysis-fiducial reco-vertex cut. Set to `9999.0` in allz configs to widen the truth denominator while keeping `vertex_cut: 60` for data and the MBD-eff numerator.
- Filename convention (Plan B, April 2026): the nominal cross-section is the all-range result. For each variant X:
  - `config_bdt_X.yaml` = all-range analysis target (bare name, no period suffix). Reads the lumi-weighted hadd of the two merge-feeder siblings.
  - `config_bdt_X_0rad.yaml` = 0 mrad merge-feeder (run 47289–51274, lumi=32.6574, lumi_target=48.9309, truth-vertex reweight ON with 0mrad file)
  - `config_bdt_X_1p5mrad.yaml` = 1.5 mrad merge-feeder (run 51274–54000, lumi=16.2735, lumi_target=48.9309, truth-vertex reweight ON with 1.5mrad file)
- `oneforall.sh` dispatch:
  - `*_0rad.yaml` / `*_1p5mrad.yaml` → `MergeSim.C` (per-period)
  - bare name WITH companion `_0rad` + `_1p5mrad` on disk → `merge_periods.sh` + `RecoEffCalculator(data)` + `CalculatePhotonYield × 2`
  - bare name WITHOUT companions → standalone `MergeSim.C` (legacy)
- `oneforall_tree_double_dispatch.sh` dispatch (Phase 1, per-sample TTree reader):
  - `*_0rad.yaml` / `*_1p5mrad.yaml` / `*_0mrad.yaml` → run DI TTree pipeline
  - bare name WITH merge-feeder companions → SKIP (no per-sample MC; handled by Phase 2)
  - bare name WITHOUT companions → run DI TTree pipeline (standalone)
- Use `make_bdt_variations.py` to generate variant configs — each variant dict must have `syst_type` / `syst_role` fields (or None for cross-checks). A variant NOT pinning its own `run_min` / `run_max` is auto-expanded into 3 configs (all-range + 2 merge-feeders) using `PER_PERIOD_OVERRIDES`. Period-pinned variants (e.g. `allz_0rad`, ntbdtpair scan) generate a single config.
- `apply_overrides` auto-defaults `lumi_target := lumi` for any variant that doesn't explicitly set it, preventing merge-feeder `lumi_target` inheritance from polluting systematic-variant semantics.
