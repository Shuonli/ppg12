# Configuration Schema Reference

Two distinct YAML schemas exist:
1. **Training config** (`FunWithxgboost/config.yaml`) -- Python BDT training
2. **Analysis config** (`efficiencytool/config_bdt_*.yaml`) -- C++ efficiency/yield pipeline

## Training Config (`config.yaml`)

### binning

| Field | Type | Value | Purpose |
|-------|------|-------|---------|
| `mode` | str | `"pt"` | `"pt"` or `"vertex"` |
| `pt_edges` | list | `[6,10,15,20,25,35]` | ET bin edges (GeV) |
| `pt_labels` | list | `["6_10",...,"25_35"]` | Bin labels |

### data

| Field | Type | Value | Purpose |
|-------|------|-------|---------|
| `features` | list[str] | 25 features | Training feature columns (order matters) |
| `pt_column` | str | `"cluster_Et"` | ET column name |
| `use_single_file_set` | bool | `true` | Single vs per-bin file sets |
| `single_file_set.signal` | list | photon txt files | Signal training data |
| `single_file_set.background` | list | jet txt files | Background training data |

### model.params

| Field | Value | Purpose |
|-------|-------|---------|
| `n_estimators` | 750 | Boosting rounds |
| `max_depth` | 5 | Tree depth |
| `learning_rate` | 0.1 | Step size |
| `subsample` | 0.5 | Row subsampling |
| `colsample_bytree` | 0.6 | Column subsampling |
| `reg_alpha` | 5.0 | L1 regularization |
| `reg_lambda` | 0.3 | L2 regularization |
| `grow_policy` | `"lossguide"` | Leaf-wise growth |

### training

| Field | Value | Purpose |
|-------|-------|---------|
| `global_seed` | 42 | Random seed |
| `train_single_model` | true | **true = ONE model total, false = per-bin models** |
| `train_size` / `val_size` / `test_size` | 0.7 / 0.1 / 0.2 | Data splits |

### reweighting

| Field | Value | Purpose |
|-------|-------|---------|
| `class_reweight` | true | Balance signal/background |
| `et_reweight` | true | Flatten ET (inverse PDF) |
| `eta_reweight` | true | Flatten eta (range hardcoded to [-0.7, 0.7]) |
| `vertex_reweight` | false | Flatten vertex z (disabled) |

**Reweighting order:** class -> per class: (eta -> ET -> vertex)

---

## Analysis Config (`config_bdt_nom.yaml`)

### input

| Field | Nominal | Purpose |
|-------|---------|---------|
| `photon_jet_file_root_dir` | `/sphenix/user/shuhangli/ppg12/FunWithxgboost/` | Sim file root dir |
| `photon_jet_file_branch_dir` | `/bdt_split.root` | Per-sample filename suffix |
| `data_file` | `part_*_with_bdt_split.root` (glob) | Data file pattern |
| `tree` | `"slimtree"` | TTree name |
| `cluster_node_name` | `"CLUSTERINFO_CEMC"` | Branch suffix |
| `bdt_model_name` | `"base_E"` | Default BDT model |
| `bdt_et_bin_edges` | `[8, 15, 35]` | ET bin edges for model selection (N+1) |
| `bdt_et_bin_models` | `["base_v1E", "base_v1E"]` | Model per ET bin (N) |

### output

| Field | Nominal | Purpose |
|-------|---------|---------|
| `eff_outfile` | `results/MC_efficiency` | MC output base path |
| `response_outfile` | `results/MC_response` | Response matrix base path |
| `data_outfile` | `results/data_histo` | Data output base path |
| `final_outfile` | `results/Photon_final` | Final result base path |
| `var_type` | `"bdt_nom"` | **Must be unique per config** |

### analysis (top-level)

| Field | Nominal | Default | Purpose |
|-------|---------|---------|---------|
| `eta_bins` | `[-0.7, 0.7]` | -- | Pseudorapidity window |
| `vertex_cut` | 60.0 | -- | Max |vertexz| (cm) |
| `truth_iso_max` | 4.0 | -- | Truth isolation cut (GeV) |
| `reco_iso_max_b` | 0.502095 | -- | Parametric iso intercept |
| `reco_iso_max_s` | 0.0433036 | -- | Parametric iso slope |
| `reco_noniso_min_shift` | 0.8 | -- | Iso/noniso gap |
| `reco_min_ET` | 5.0 | -- | Min cluster ET |
| `cone_size` | 3 | -- | Cone (3 = R=0.3) |
| `use_topo_iso` | 2 | 0 | 0=tower, 1=topo R=0.3, 2=topo R=0.4 |
| `iso_threshold` | 1 | 0 | 1=threshold-corrected tower iso |
| `mc_iso_scale` | 1.2 | **1.0** | MC iso multiplicative correction |
| `mc_iso_shift` | 0.1 | **0.0** | MC iso additive correction |
| `cluster_escale` | 1.0 | 1.0 | Energy scale factor (systematics) |
| `cluster_eres` | 0.0 | 0.0 | Energy resolution smearing |
| `trigger_used` | 30 | -- | Trigger bit (supports list) |
| `pT_bins` | `[8..28,32,36]` | -- | Reco pT bin edges (12 bins) |
| `pT_bins_truth` | `[7..36,45]` | -- | Truth pT bins (extended) |
| `fit_option` | 1 | **0** | **0=erf (default), 1=Pade** |
| `mc_purity_correction` | 0 | 0 | Apply MC purity correction |
| `lumi` | 16.2735 | 49.562 | Data luminosity (pb⁻¹). Used by CalculatePhotonYield for 1/binwidth/lumi data normalization. |
| `lumi_target` | (= lumi) | (= lumi) | **NEW**: target lumi for the per-event MC scaling `weight *= lumi/lumi_target`. Default = lumi → no scaling (standalone per-period). Set to `sum(L_periods)` in merge-feeder configs so plain hadd reproduces all-range MC. See `wiki/pipeline/full-run-range.md`. |
| `vertex_cut` | 60.0 | -- | Max |z_reco| (cm). Applied to data + MC events and the MBD-eff numerator. |
| `vertex_cut_truth` | (= vertex_cut) | (= vertex_cut) | **NEW**: max |z_truth| (cm) for the MBD-eff denominator only (`RecoEffCalculator_TTreeReader.C:1785`). Decouples the truth-vertex cut from the reco fiducial cut. Set to 9999 for the allz cross-check. |
| `run_min` / `run_max` | 51274 / 54000 | -1 / -1 | Run range filter |
| `vertex_reweight_on` | 1 | 1 | Enable vertex reweighting |
| `n_nt_fail` | 0 | 1 | Cuts to fail for non-tight |
| `mbd_avg_sigma_max` | (not in nom) | 999.0 | Max MBD avg sigma for pileup rejection |
| `mbd_avg_sigma_min` | (not in nom) | 0.0 | Min MBD avg sigma |
| `cluster_mbd_time_min` | -999.0 | -999.0 | Min cluster-MBD timing difference |
| `cluster_mbd_time_max` | 999.0 | 999.0 | Max cluster-MBD timing difference |
| `common_b2bjet_cut` | 0 | 0 | Enable back-to-back jet veto |
| `common_b2bjet_pt_min` | (not in nom) | 7.0 | Min jet pT for b2b veto |

### analysis.tight / analysis.non_tight

Both sections have identical field names. Key fields:

**Current nominal: ntbdtpair `t80_70_h70_70_l60_40` (April 2026).** Tight cut goes 0.80 → 0.70 (ET=10 → 25 GeV); non-tight is bounded above by a flat 0.70 cap and below by a parametric floor 0.60 → 0.40.

| Field | Tight nominal | Non-tight nominal | Purpose |
|-------|---------------|-------------------|---------|
| `bdt_min` | 0.7 | 0.02 | Flat BDT lower bound (legacy fallback when parametric fields absent) |
| `bdt_max` | 1.0 | 0.5 | Flat BDT upper bound (legacy fallback) |
| `bdt_min_slope` | -0.00667 | -0.01333 | Parametric BDT lower-bound slope (tight: cut on tight ABCD; non-tight: floor of non-tight window) |
| `bdt_min_intercept` | 0.8667 | 0.7333 | Parametric BDT lower-bound intercept |
| `bdt_max_slope` | n/a | 0.0 | Parametric upper BDT slope (non-tight: ceiling of non-tight window) |
| `bdt_max_intercept` | n/a | 0.7 | Parametric upper BDT intercept |
| `weta_cogx_min/max` | 0.0 / 1.0 | 0.0 / 1.0 | Width bounds (flat) |
| `weta_cogx_max_b/max_s` | 1.0 / 0.0 | 1.0 / 0.0 | Width bounds (parametric) |

**Code support:** `RecoEffCalculator_TTreeReader.C` and `ShowerShapeCheck.C` both read all six parametric fields (slope+intercept for tight.bdt_min, non_tight.bdt_max, non_tight.bdt_min). Older macros (`RecoEffCalculator.C`, `EtaMigrationStudy.C`, `Cluster_rbr.C`, `NPB_PurityStudy.C`, `DoubleInteractionCheck.C`) still only read the flat `bdt_min` and apply ET-independent cuts.

### analysis.common

| Field | Nominal | Purpose |
|-------|---------|---------|
| `wr_cogx_bound` | 0.0 | Max radial width pre-selection (0 = disabled) |
| `cluster_weta_cogx_bound` | 2.0 | Max weta pre-selection |
| `npb_cut_on` | 1 | Enable NPB score cut |
| `npb_score_cut` | 0.5 | NPB threshold |
| `e11_over_e33_min/max` | 0.0 / 0.98 | e11/e33 pre-selection |

### analysis.unfold

| Field | Nominal | Purpose |
|-------|---------|---------|
| `resultit` | 2 | Bayesian iterations |
| `resultleak` | 1 | Apply leakage correction |
| `reweight` | 1 | Response reweighting (forced to 0 in code) |

## Run Period / Luminosity

|z|<60 fiducial lumi (`60cmLumi_fromJoey.list`):

| Period | run_min | run_max | lumi (pb⁻¹) | lumi_target (pb⁻¹) |
|--------|---------|---------|-------------|--------------------|
| 0 mrad merge-feeder | 47289 | 51274 | 32.6574 | 48.9309 |
| 1.5 mrad merge-feeder (nominal base) | 51274 | 54000 | 16.2735 | 48.9309 |
| All runs (analysis target) | 47289 | 54000 | 48.9309 | (= lumi) |
| No filter | -1 | -1 | depends | — |

Beam-delivered lumi (`allzLumi_fromJoey.list`, used by allz cross-check):

| Period | run_min | run_max | lumi (pb⁻¹) | lumi_target (pb⁻¹) |
|--------|---------|---------|-------------|--------------------|
| 0 mrad allz | 47289 | 51274 | 47.2076 | 64.3718 |
| 1.5 mrad allz | 51274 | 54000 | 17.1642 | 64.3718 |
| All allz (analysis target) | 47289 | 54000 | 64.3718 | (= lumi) |

Per-period merge-feeders use `lumi_target = sum(L_periods)` so per-event
`lumi_weight = lumi/lumi_target` pre-scales MC fills; plain hadd of two
merge-feeder outputs reproduces the all-range expectation. See
`wiki/pipeline/full-run-range.md`.

## Critical Sync Requirements

- `var_type` must be unique per config (prevents silent overwrite)
- `bdt_et_bin_edges` needs N+1 entries for N entries in `bdt_et_bin_models`
- `pT_bins` should match `plotcommon.h:ptRanges` for consistent plotting
- 25 features in training config must match feature vectors in `apply_BDT.C`
- Configs ending in `_all.yaml` are dispatched by `oneforall.sh` to `merge_periods.sh` instead of `MergeSim.C`. The substitution `_all.yaml → _0rad.yaml + _1p5mrad.yaml` (with the bare `config_bdt_all.yaml` mapping to `_nom.yaml`) names the per-period merge-feeder pair.
- Merge-feeder configs MUST set `lumi_target = sum(L_periods)`. `make_bdt_variations.py:apply_overrides` auto-defaults `lumi_target := lumi` for variants WITHOUT explicit override, so systematic variants regenerated from a merge-feeder base do NOT inherit `lumi_target` and stay as standalone per-period systematics.
