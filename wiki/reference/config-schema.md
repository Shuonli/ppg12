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
| `lumi` | 16.2735 | 49.562 | Luminosity (pb^-1) |
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

| Field | Tight nominal | Non-tight nominal | Purpose |
|-------|---------------|-------------------|---------|
| `bdt_min` | 0.7 | 0.02 | Flat BDT lower bound |
| `bdt_max` | 1.0 | 0.5 | Flat BDT upper bound |
| `bdt_min_slope` | -0.015 | -0.015 | Parametric BDT slope |
| `bdt_min_intercept` | 0.80 | 0.80 | Parametric BDT intercept |
| `bdt_max_slope` | (not set) | (not set) | Parametric upper BDT slope |
| `bdt_max_intercept` | (not set) | (not set) | Parametric upper BDT intercept |
| `weta_cogx_min/max` | 0.0 / 1.0 | 0.0 / 1.0 | Width bounds (flat) |
| `weta_cogx_max_b/max_s` | 1.0 / 0.0 | 1.0 / 0.0 | Width bounds (parametric) |

**Config/code mismatch in non_tight:** The config file (config_bdt_nom.yaml) has `bdt_min_slope` and `bdt_min_intercept` under `non_tight` (lines 205-206), but the code (`RecoEffCalculator_TTreeReader.C` line 479) reads `bdt_max_slope` and `bdt_max_intercept` for the non_tight upper boundary. Since these fields do not exist in the config, the code falls back to defaults: `bdt_max_slope=0` (flat) and `bdt_max_intercept=bdt_max` (0.5). This means the non_tight parametric upper BDT boundary is effectively flat at 0.5, not the intended parametric `-0.015*ET + 0.80`.

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

| Period | run_min | run_max | lumi (pb^-1) |
|--------|---------|---------|-------------|
| 0 mrad | 47289 | 51274 | 32.6574 |
| 1.5 mrad (nominal) | 51274 | 54000 | 16.2735 |
| All runs | 47289 | 54000 | 49.562 |
| No filter | -1 | -1 | depends |

## Critical Sync Requirements

- `var_type` must be unique per config (prevents silent overwrite)
- `tight.bdt_min_intercept` should match `non_tight.bdt_max_intercept`
- `bdt_et_bin_edges` needs N+1 entries for N entries in `bdt_et_bin_models`
- `pT_bins` should match `plotcommon.h:ptRanges` for consistent plotting
- 25 features in training config must match feature vectors in `apply_BDT.C`
