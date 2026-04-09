# BDT Pipeline Report -- FunWithxgboost/

## 1. Purpose

The BDT pipeline implements photon identification for the PPG12 isolated photon cross-section analysis in pp collisions at sqrt(s) = 200 GeV at sPHENIX. It trains XGBoost binary classifiers to discriminate direct photon clusters (signal) from QCD jet fragments (background) in the barrel EMCal, using calorimeter shower-shape variables. A secondary BDT (NPB score) rejects non-physics background (NPB) detector artifacts.

The end-to-end flow is:
1. Feature extraction from ROOT slimtrees to space-separated text files (BDTinput.C)
2. XGBoost training with kinematic reweighting (main_training.py + modules)
3. Model export to TMVA ROOT format for C++ inference
4. Score application back to slimtrees (apply_BDT.C)

---

## 2. Key Files

| File | Lines | Description |
|------|-------|-------------|
| `BDTinput.C` | 1040 | Feature extraction from slimtree ROOT files to training text files |
| `main_training.py` | 1295 | Main training pipeline orchestrator (BinnedTrainingPipeline class) |
| `data_loader.py` | 338 | Data loading from txt/csv/ROOT files, signal/background labeling, ET subsampling |
| `model_builder.py` | 77 | sklearn Pipeline construction with XGBoost or HGB classifier |
| `reweighting.py` | 165 | Kinematic reweighting: class, ET, eta, vertex z (inverse-PDF method) |
| `plotting.py` | 1222 | BinnedTrainingPlotter: ROC, AUC bars, score distributions, feature importance, correlations |
| `apply_BDT.C` | 492 | Apply trained TMVA BDT models to slimtree, write score branches |
| `apply_NPB_score.C` | 364 | Standalone NPB score application (legacy, now integrated into apply_BDT.C) |
| `train_npb_score.py` | 1073 | NPBTrainer class: train BDT for NPB rejection |
| `config.yaml` | 200 | Training configuration (features, binning, XGB params, reweighting) |
| `config_nom.yaml` | 200 | Analysis/efficiency tool configuration (cuts, paths, ABCD regions) |
| `config_npb_training.yaml` | 139 | NPB score training configuration |
| `make_split_configs.py` | 107 | Generate per-model-variant training configs for 11 feature-set variants |
| `train_split_model.py` | 8 | Thin wrapper to run BinnedTrainingPipeline on a single config |
| `export_cluster_arrays_csv.C` | 194 | Export raw cluster tower arrays (49 towers) to CSV for CNN/alternative studies |

---

## 3. Feature Extraction (BDTinput.C)

### Function Signature
```cpp
void BDTinput(const std::string &configname = "config_nom.yaml",
              const std::string filetype = "data")
```

### Branches Read from slimtree

All cluster branches are suffixed with `_{clusternodename}` (currently `CLUSTERINFO_CEMC`).

**Event-level:**
- `mbdnorthhit`, `mbdsouthhit` (int) -- MBD hit multiplicity
- `vertexz` (float) -- reconstructed vertex z position
- `energy_scale` (float) -- truth energy scale
- `runnumber` (int) -- run number
- `scaledtrigger[N]` (Bool_t array) -- trigger bits
- `livetrigger[N]` (Bool_t array) -- live trigger bits
- `nparticles` (int) -- number of truth particles (MC only)
- `particle_trkid[nparticles]` (int array) -- truth track IDs
- `particle_photonclass[nparticles]` (int array) -- truth photon classification (pid)
- `njet_truth` (int), `jet_truth_Pt[njet_truth]` (float array) -- truth jets (MC only)

**Per-cluster (all arrays of size ncluster):**
- `cluster_Et`, `cluster_Eta`, `cluster_Phi` -- kinematics
- `cluster_prob` -- chi2 probability
- `cluster_truthtrkID` -- truth-matching track ID (MC)
- `cluster_iso_02`, `cluster_iso_03`, `cluster_iso_04` -- isolation energy in cone R=0.2/0.3/0.4
- `cluster_iso_03_60_emcal`, `cluster_iso_03_60_hcalin`, `cluster_iso_03_60_hcalout` -- threshold isolation (EMCal+HCal components)
- `cluster_et1`, `cluster_et2`, `cluster_et3`, `cluster_et4` -- energy fraction in 1st/2nd/3rd/4th tower ring
- `cluster_weta_cogx`, `cluster_wphi_cogx` -- width in eta/phi (center-of-gravity)
- `cluster_nsaturated` -- number of saturated towers
- `cluster_e11`, `cluster_e22`, `cluster_e13`, `cluster_e15`, `cluster_e17` -- energy in NxM tower windows
- `cluster_e31`, `cluster_e51`, `cluster_e71` -- energy in transposed tower windows
- `cluster_e33`, `cluster_e35`, `cluster_e37`, `cluster_e53`, `cluster_e73`, `cluster_e77` -- energy in larger windows
- `cluster_w32`, `cluster_w52`, `cluster_w72` -- tower width ratios
- `cluster_e32` -- 3x2 tower energy

**NPB-specific (data only, loaded when `npb_cut_on=1`):**
- `mbd_time` (float) -- MBD timing
- `njet` (int), `jet_Pt[njet]` (float), `jet_Phi[njet]` (float) -- reco jets
- `cluster_time_array[ncluster*49]` (float) -- per-tower timing
- `cluster_e_array[ncluster*49]` (float) -- per-tower energy
- `cluster_ownership_array[ncluster*49]` (int) -- tower ownership flags

### Derived Features (computed in BDTinput.C)

All ratios use a `safe_div` lambda that returns 0.0 on div-by-zero or non-finite inputs:

| Feature | Formula |
|---------|---------|
| `e11_over_e33` | `cluster_e11 / cluster_e33` |
| `e32_over_e35` | `cluster_e32 / cluster_e35` |
| `e11_over_e22` | `cluster_e11 / cluster_e22` |
| `e11_over_e13` | `cluster_e11 / cluster_e13` |
| `e11_over_e15` | `cluster_e11 / cluster_e15` |
| `e11_over_e17` | `cluster_e11 / cluster_e17` |
| `e11_over_e31` | `cluster_e11 / cluster_e31` |
| `e11_over_e51` | `cluster_e11 / cluster_e51` |
| `e11_over_e71` | `cluster_e11 / cluster_e71` |
| `e22_over_e33` | `cluster_e22 / cluster_e33` |
| `e22_over_e35` | `cluster_e22 / cluster_e35` |
| `e22_over_e37` | `cluster_e22 / cluster_e37` |
| `e22_over_e53` | `cluster_e22 / cluster_e53` |

### Output Format

Space-separated text files with a header row. 30 columns total:

```
cluster_Et cluster_Eta cluster_Phi vertexz
e11_over_e33 e32_over_e35 e11_over_e22 e11_over_e13
e11_over_e15 e11_over_e17 e11_over_e31
e11_over_e51 e11_over_e71 e22_over_e33
e22_over_e35 e22_over_e37 e22_over_e53
cluster_prob cluster_weta_cogx cluster_wphi_cogx
cluster_et1 cluster_et2 cluster_et3 cluster_et4
cluster_w32 cluster_w52 cluster_w72
recoisoET
is_tight pid
```

**Output filenames:**
- MC: `shapes_split_{filetype}.txt` (e.g., `shapes_split_photon5.txt`)
- Data NPB: `shapes_split_data_npb.txt`
- Legacy (NOSPLIT): `NOSPLIT/shapes_{filetype}.txt`

### Signal/Background Labeling

The `pid` column encodes truth information:
- **Signal (photon MC):** `pid` = `particle_photonclass` value (1 = direct photon, 2 = fragmentation photon)
- **Background (jet MC):** `pid` = `particle_photonclass` value (0 or other for non-photon)
- **Data NPB:** `pid` = -1 (hardcoded for NPB-tagged clusters)

In the Python training code, signal selection is `pid in [1, 2]`; background selection is `pid not in [1, 2]`.

### Selection Cuts Applied in BDTinput.C

- Cluster ET: 6 < ET < 40 GeV (hardcoded)
- Cluster eta: |eta| < 0.7 (hardcoded)
- Vertex: |vertexz| < vertex_cut (from config, default 100 cm)
- MBD: mbdnorthhit >= 1 AND mbdsouthhit >= 1
- Trigger: scaledtrigger[trigger_used] == 1 (data only)
- Skip run numbers: {47698, 51489, 51721, 51725, 53284}
- Skip event 16152886 (hardcoded)
- `common_pass` pre-selection (configurable via `write_outside_common_pass_mc`)
- NPB tagging for data: badtime AND no back-to-back jet AND weta >= npb_weta_min

### Cross-Section Weights

From `CrossSectionWeights.h` (`namespace PPG12`):

| Sample | Cross-section (pb) | Weight (vs reference) |
|--------|-------------------|----------------------|
| photon5 | 146359.3 | photon5cross / photon20cross = 1122.0 |
| photon10 | 6944.675 | photon10cross / photon20cross = 53.24 |
| photon20 | 130.4461 | 1.0 (reference) |
| jet10 | 3.997e+06 | jet10cross / jet30cross = 1580.0 |
| jet15 | 4.073e+05 | jet15cross / jet30cross = 160.9 |
| jet20 | 6.2623e+04 | jet20cross / jet30cross = 24.75 |
| jet30 | 2.5298e+03 | 1.0 (reference) |
| jet50 | 7.3113 | jet50cross / jet30cross = 0.00289 |

Note: BDTinput.C normalizes photon weights relative to photon20, and jet weights relative to jet30. However, these weights are computed in BDTinput.C but are NOT written to the output txt files -- they are only used internally. The cross-section weighting in BDTinput.C applies to per-run luminosity counting, not to training sample weights.

---

## 4. Training Pipeline

### Architecture

```
main_training.py (BinnedTrainingPipeline)
    |
    +-- data_loader.py (DataLoader)
    |       loads txt files, labels signal/background,
    |       applies background subsampling
    |
    +-- reweighting.py (KinematicReweighter)
    |       class balance, ET flattening, eta flattening, vertex z flattening
    |
    +-- model_builder.py (ModelBuilder)
    |       builds sklearn Pipeline: Imputer -> (Scaler) -> XGBClassifier
    |
    +-- plotting.py (BinnedTrainingPlotter)
            ROC curves, AUC bars, score distributions,
            feature importance, correlation analysis
```

### Flow (BinnedTrainingPipeline.run_training)

1. **Load config** (`_load_config`): YAML -> dict
2. **Validate config** (`_validate_config`): checks required keys, binning consistency
3. **Get binning** (`_get_binning_config`): determines bin edges/labels based on mode (pt or vertex)
4. **Load data globally** (`DataLoader.get_all_bin_data`):
   - When `use_single_file_set: true`, loads all signal+background txt files at once
   - Signal: `pid in [1, 2]` -> label=1
   - Background: `pid not in [1, 2]` -> label=0
   - Optional background subsampling below `et_threshold` (flatten or fixed fraction)
5. **Apply reweighting globally** (`KinematicReweighter.apply_all_reweighting`):
   - Class reweighting -> ET reweighting -> eta reweighting -> vertex z reweighting
   - Applied to full dataset BEFORE splitting into bins
6. **Analyze feature correlations** (`analyze_feature_correlations`):
   - Pearson correlation between each training feature and isoET for background samples
7. **Train model** (single or per-bin):
   - **Single model** (`_train_single_model`): one model for all ET bins, evaluated per-bin on test set
   - **Per-bin models** (`_train_per_bin_models`): separate model per ET bin
   - Optional Optuna hyperparameter tuning
   - Train/val/test split: 70/10/20% (configurable), stratified by label

### Post-Training Steps

8. **Evaluate** (`evaluate_models`): AUC on test set, ROC curves, isoET-BDT correlations
9. **Generate plots** (`generate_plots`): reweighting QA, score distributions, AUC bars, ROC, feature importance, 2D correlations
10. **Save results** (`save_results`): joblib models, TMVA ROOT models, metrics CSV, metadata JSON

---

## 5. Config Schema (config.yaml)

### `binning`
| Field | Type | Description |
|-------|------|-------------|
| `mode` | str | `"pt"` or `"vertex"` -- binning dimension |
| `pt_edges` | list[float] | ET bin edges in GeV: `[6, 10, 15, 20, 25, 35]` |
| `pt_labels` | list[str] | Bin labels: `["6_10", "10_15", "15_20", "20_25", "25_35"]` |
| `vertex_edges` | list[float] | Vertex z bin edges in cm: `[-100, -30, 30, 100]` |
| `vertex_labels` | list[str] | Vertex bin labels |
| `edges`, `labels` | list | Backward-compatible aliases for pt_edges/pt_labels |

### `data`
| Field | Type | Description |
|-------|------|-------------|
| `pt_column` | str | Column name for ET: `"cluster_Et"` |
| `iso_column` | str | Column name for isolation: `"recoisoET"` |
| `label_column` | str | Column name for label: `"label"` |
| `features` | list[str] | 25 ordered feature names for training |
| `use_single_file_set` | bool | If true, use `single_file_set`; else `per_bin_file_sets` |
| `single_file_set.signal` | list[str] | Signal txt file paths (photon MC) |
| `single_file_set.background` | list[str] | Background txt file paths (jet MC) |
| `per_bin_file_sets` | dict | Per-ET-bin file sets (legacy mode) |
| `background_subsampling.enabled` | bool | Enable low-ET background subsampling |
| `background_subsampling.et_threshold` | float | ET threshold for subsampling (15.0 GeV) |
| `background_subsampling.flatten_distribution` | bool | Flatten ET dist below threshold |
| `background_subsampling.n_et_bins` | int | Number of bins for flattening (20) |
| `background_subsampling.subsample_fraction` | float | Fixed fraction to keep (0.3) |
| `background_subsampling.random_seed` | int | Seed for reproducibility (42) |

### `model`
| Field | Type | Description |
|-------|------|-------------|
| `classifier` | str | `"xgb"` (XGBoost) or `"hgb"` (HistGradientBoosting) |
| `params` | dict | XGBoost hyperparameters (see section 6) |
| `use_scaler` | bool | If true, add StandardScaler to pipeline (false) |

### `training`
| Field | Type | Description |
|-------|------|-------------|
| `global_seed` | int | Random seed for all splits and models: `42` |
| `stratify` | bool | Stratified train/test split by label (true) |
| `train_single_model` | bool | Train one model for all bins (true) vs per-bin (false) |
| `train_size` | float | Training fraction: `0.7` |
| `val_size` | float | Validation fraction: `0.1` |
| `test_size` | float | Test fraction: `0.2` |
| `enable_hyperparameter_tuning` | bool | Enable Optuna tuning (false) |
| `n_trials` | int | Number of Optuna trials (5) |
| `n_cv_folds` | int | K-fold CV splits (3) |
| `tuning_timeout_seconds` | int/null | Timeout for tuning |
| `tuning_data_fraction` | float | Fraction of data for tuning (1.0) |

### `reweighting`
| Field | Type | Description |
|-------|------|-------------|
| `class_reweight` | bool | Balance signal/background weights (true) |
| `et_reweight` | bool | Flatten ET distribution per class (true) |
| `et_reweight_bins` | int | Number of bins for ET density estimation (20) |
| `et_reweight_max` | float | Maximum ET weight cap (800.0) |
| `eta_reweight` | bool | Flatten eta distribution per class (true) |
| `eta_reweight_bins` | int | Number of bins (20) |
| `eta_reweight_range` | list[float] | Eta range `[-2.5, 2.5]` |
| `vertex_reweight` | bool | Flatten vertex z distribution (false) |
| `vertex_reweight_bins` | int | Number of bins (20) |
| `vertex_reweight_max` | float | Maximum vertex weight cap (50.0) |

### `output`
| Field | Type | Description |
|-------|------|-------------|
| `output_dir` | str | Model output directory: `"binned_models"` |
| `save_models` | bool | Save joblib models (true) |
| `save_metrics` | bool | Save per-bin metrics CSV (true) |
| `save_plots` | bool | Save diagnostic plots (true) |
| `save_tmva` | bool | Save TMVA ROOT format (true) |
| `root_file_prefix` | str | Prefix for ROOT files: `"model_insitu_E"` |

### `plotting`
| Field | Type | Description |
|-------|------|-------------|
| `context` | str | Seaborn context: `"talk"` |
| `density` | bool | Normalize histograms (true) |
| `figure_size` | list[float] | Figure size: `[7.5, 5]` |
| `grid` | bool | Show grid (true) |
| `n_bins` | int | Number of histogram bins (50) |
| `zoom_y_axis` | bool | Zoom 2D correlation y-axis (true) |

### `tuning`
| Field | Type | Description |
|-------|------|-------------|
| `baseline_auc_threshold` | float | Minimum AUC to accept tuned params (0.50) |
| `learning_rate_min/max` | float | Log-uniform range for learning rate |
| `max_depth_min/max` | int | Integer range for tree depth |
| `n_estimators_min/max/step` | int | Range and step for boosting rounds |
| `subsample_min/max` | float | Range for row subsampling |
| `reg_alpha_min/max` | float | Log-uniform range for L1 regularization |
| `reg_lambda_min/max` | float | Log-uniform range for L2 regularization |
| `colsample_bytree_min/max` | float | Range for column subsampling |

---

## 6. Model Architecture

### XGBoost Hyperparameters (config.yaml defaults)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `objective` | `binary:logistic` | Binary classification (set in ModelBuilder) |
| `eval_metric` | `auc` | Evaluation metric (set in ModelBuilder) |
| `n_estimators` | 750 | Number of boosting rounds |
| `max_depth` | 5 | Maximum tree depth |
| `learning_rate` | 0.1 | Step size shrinkage |
| `subsample` | 0.5 | Row subsampling ratio |
| `colsample_bytree` | 0.6 | Column subsampling per tree |
| `colsample_bylevel` | 1.0 | Column subsampling per level |
| `reg_alpha` | 5.0 | L1 regularization |
| `reg_lambda` | 0.3 | L2 regularization |
| `grow_policy` | `lossguide` | Leaf-wise tree growth |
| `tree_method` | `hist` | Histogram-based algorithm |
| `max_bin` | 256 | Maximum histogram bins |
| `n_jobs` | 4 | Parallel threads |
| `random_state` | 42 | Set from global_seed |

### ET Binning

5 bins for training: `[6, 10, 15, 20, 25, 35]` GeV

When `train_single_model: true` (current default), a single model is trained on all ET bins but evaluated per-bin.

### Train/Val/Test Split

- **Train**: 70% -- used for fitting
- **Validation**: 10% -- used for hyperparameter tuning (when enabled) and intermediate evaluation
- **Test**: 20% -- final unbiased evaluation, never seen during training
- Split is stratified by label (signal/background)
- First split separates test set, second split separates train/validation
- Different random seeds for each split: `global_seed` for first, `global_seed + 1` for second

### sklearn Pipeline Structure

```python
Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    # ("scaler", StandardScaler()),  # optional, disabled by default
    ("clf", XGBClassifier(...))
])
```

---

## 7. Reweighting (reweighting.py)

All reweighting is applied globally BEFORE binning/splitting. The `apply_all_reweighting` method chains four steps:

### Class Reweighting (`apply_class_reweighting`)
- Balances signal and background to contribute equally
- `weight_sig = (n_sig + n_bkg) / (2 * n_sig)`
- `weight_bkg = (n_sig + n_bkg) / (2 * n_bkg)`
- Sets initial `df["weight"] = class_weight`

### ET Reweighting (`apply_et_reweighting`, per class)
- Goal: flatten the ET distribution within each class
- Algorithm:
  1. Bin ET values into `et_reweight_bins` (20) bins spanning data range
  2. Compute histogram density estimate, multiply by n_bins to normalize
  3. Fit `UnivariateSpline(s=0.0)` through bin centers
  4. Evaluate spline at each sample's ET
  5. Clip PDF minimum to 1e-3
  6. Weight = 1.0 / PDF
  7. Cap weight at `et_reweight_max` (800.0)
  8. Normalize: weights /= mean(weights)
  9. Multiply into existing `df["weight"]`
- Applied independently for signal (label=1) and background (label=0)

### Eta Reweighting (`apply_eta_reweighting`, per class)
- Same inverse-PDF method as ET reweighting
- Bins: 20 bins in `[-0.7, 0.7]`
- Hardcoded eta range `[-0.7, 0.7]` matching barrel EMCal
- No weight cap (only PDF floor of 1e-3)

### Vertex Z Reweighting (`apply_vertex_reweighting`, per class)
- Disabled by default (`vertex_reweight: false`)
- Same inverse-PDF method
- Column: `vertexz`
- Weight cap: `vertex_reweight_max` (50.0)

### Final weight

The `weight` column in the DataFrame is the product of all enabled reweighting factors:
```
weight = class_weight * et_weight * eta_weight * vertex_weight
```

This weight is passed to XGBoost via `clf__sample_weight` during training.

---

## 8. Model Export

### Formats Saved

1. **Joblib** (`model_{bin_label}.joblib`): Full sklearn Pipeline for Python inference
2. **TMVA ROOT** (`{prefix}_single_tmva.root`): For C++ inference via `TMVA::Experimental::RBDT`
3. **Metrics CSV** (`per_bin_metrics.csv`): AUC, correlations, timing per bin
4. **Metadata JSON** (`training_metadata.json`): Full config, metrics, correlations

### TMVA Export Procedure (`_save_tmva_models`)

```python
booster = xgb_model.get_booster()
booster.feature_names = [f"f{i}" for i in range(len(features))]
ROOT.TMVA.Experimental.SaveXGBoost(xgb_model, "myBDT", tmva_path, num_inputs=len(features))
```

Critical: Feature names are set to `f0, f1, f2, ...` -- the TMVA model uses positional indexing, not named features. The feature ORDER in the training config must exactly match the order used when calling `RBDT.Compute()` in `apply_BDT.C`.

### Model File Naming

- Single model mode: `{output_dir}/{root_file_prefix}_single_tmva.root`
- Per-bin mode: `{output_dir}/{root_file_prefix}_{bin_label}_tmva.root`
- Standard naming: `binned_models/model_{variant_name}_single_tmva.root`
- Split naming: `binned_models/model_{variant_name}_split_single_tmva.root`

### Existing Models on Disk

Standard models (11 variants): `model_{base,base_vr,base_v0,...,base_v3E}_single_tmva.root`
Split models (11 variants): `model_{base,...,base_v3E}_split_single_tmva.root`
Vertex-binned models (3 bins): `model_{-100_-30,-30_30,30_100}_tmva.root`
Other models: `model_single_tmva.root`, `model_insitu_E_single_tmva.root`, `model_vtx_reweight_single_tmva.root`, `model_vtx_noreweight_single_tmva.root`

---

## 9. BDT Application (apply_BDT.C)

### Function Signature
```cpp
void apply_BDT(const std::string &configname = "config_nom.yaml",
               const std::string filetype = "jet12_double",
               const std::string inputfilename = "...")
```

### Input Branches Read

Same cluster node pattern: `{branchname}_{clusternodename}`.

**Core branches (photon-ID BDT inputs):**
- `vertexz` -- vertex z position
- `ncluster_{node}` -- number of clusters
- `cluster_Et_{node}`, `cluster_Eta_{node}`, `cluster_Phi_{node}` -- cluster kinematics
- `cluster_et1_{node}`, `cluster_et2_{node}`, `cluster_et3_{node}`, `cluster_et4_{node}` -- ring energies
- `cluster_weta_cogx_{node}`, `cluster_wphi_cogx_{node}` -- shower widths
- `cluster_e11_{node}`, `cluster_e22_{node}`, `cluster_e33_{node}`, `cluster_e35_{node}`, `cluster_e32_{node}` -- tower energies
- `cluster_e13_{node}`, `cluster_e15_{node}`, `cluster_e17_{node}` -- asymmetric windows
- `cluster_e31_{node}`, `cluster_e51_{node}`, `cluster_e71_{node}` -- transposed windows
- `cluster_e37_{node}`, `cluster_e53_{node}` -- additional windows
- `cluster_prob_{node}` -- chi2 probability
- `cluster_w32_{node}`, `cluster_w52_{node}`, `cluster_w72_{node}` -- width ratios

**NPB score inputs (additional):**
- `cluster_iso_02_{node}`, `cluster_iso_03_{node}`, `cluster_iso_04_{node}` -- isolation
- `cluster_iso_03_60_emcal_{node}`, `cluster_iso_03_60_hcalin_{node}`, `cluster_iso_03_60_hcalout_{node}` -- threshold isolation

### Derived Variables Computed

```cpp
e11_over_e33_BDT = (cluster_e11 > 0) ? cluster_e11 / cluster_e33 : 0;
e32_over_e35_BDT = (cluster_e32 > 0) ? cluster_e32 / cluster_e35 : 0;
```

Note: The divide-by-zero protection in apply_BDT.C checks `numerator > 0`, while BDTinput.C uses `safe_div` which checks `denominator != 0` and `isfinite()`. This is a minor inconsistency but functionally equivalent for positive energies.

### 11 Model Variants and Feature Vectors

All models are loaded from `binned_models/model_{name}[_split]_single_tmva.root` and applied via `TMVA::Experimental::RBDT`.

**ET threshold**: Score is computed only when `cluster_Et > 7 GeV`; otherwise set to -1.

| # | Model Name | Feature Vector (exact order) | N features |
|---|-----------|------------------------------|-----------|
| 0 | `base` | vertexz, Eta, e11/e33, et1, et2, et3, et4 | 7 |
| 1 | `base_vr` | vertexz, Eta, e11/e33, et1, et2, et3, et4 | 7 |
| 2 | `base_v0` | vertexz, Eta, e11/e33, et2, et3, et4 | 6 |
| 3 | `base_v1` | weta, vertexz, Eta, e11/e33, et1, et2, et3, et4 | 8 |
| 4 | `base_v2` | weta, wphi, vertexz, Eta, e11/e33, et1, et2, et3, et4 | 9 |
| 5 | `base_v3` | weta, wphi, vertexz, Eta, e11/e33, et1, et2, et3, et4, e32/e35 | 10 |
| 6 | `base_E` | Et, vertexz, Eta, e11/e33, et1, et2, et3, et4 | 8 |
| 7 | `base_v0E` | Et, vertexz, Eta, e11/e33, et2, et3, et4 | 7 |
| 8 | `base_v1E` | Et, weta, vertexz, Eta, e11/e33, et1, et2, et3, et4 | 9 |
| 9 | `base_v2E` | Et, weta, wphi, vertexz, Eta, e11/e33, et1, et2, et3, et4 | 10 |
| 10 | `base_v3E` | Et, weta, wphi, vertexz, Eta, e11/e33, et1, et2, et3, et4, e32/e35 | 11 |

**Key differences:**
- `base` vs `base_vr`: Identical feature sets (vr = "vertex reweighted" variant trained with vertex reweighting enabled)
- `v0`: Drops `et1` (most correlated with isolation)
- `v1`: Adds `weta_cogx`
- `v2`: Adds `weta_cogx` + `wphi_cogx`
- `v3`: Adds `weta_cogx` + `wphi_cogx` + `e32_over_e35`
- `_E` suffix: Adds `cluster_Et` as first feature (energy-aware models)

### Output Branches Written

For each of the 11 models:
```
cluster_bdt_{clusternodename}_{model_name}[ncluster_{clusternodename}]/F
```
Example: `cluster_bdt_CLUSTERINFO_CEMC_base_v3E[ncluster_CLUSTERINFO_CEMC]/F`

NPB score:
```
cluster_npb_score_{clusternodename}[ncluster_{clusternodename}]/F
```

### Output Files

- **Sim**: `{filetype}/bdt_split.root` (e.g., `photon5/bdt_split.root`)
- **Data**: `{inputfilename_without_ext}_with_bdt_split.root`

The output ROOT file is a full clone of the input slimtree with the BDT score branches appended.

### Model Loading Pattern

```cpp
const bool use_split_bdt = configYaml["analysis"]["use_split_bdt"].as<int>(0) != 0;
const std::string model_suffix = use_split_bdt ? "_split" : "";
// File: "binned_models/model_" + name + model_suffix + "_single_tmva.root"
TMVA::Experimental::RBDT bdt("myBDT", model_file);
// Score: bdt.Compute(feature_vector)[0]
```

---

## 10. NPB Score (train_npb_score.py)

### Purpose

The NPB (Non-Physics Background) score model discriminates real physics clusters from detector timing artifacts. NPB clusters are caused by out-of-time energy deposits that survive reconstruction.

### NPBTrainer Class

```python
class NPBTrainer:
    def __init__(self, config_path)  # Load config
    def load_data(self)              # Load signal + background, apply cuts
    def split_data(self, df)         # Train/val/test split (70/10/20%)
    def compute_weights(self, y, et) # Class balance + optional ET reweighting
    def train_model(...)             # XGBoost training with optional tuning
    def evaluate_model(...)          # Test set evaluation
    def generate_plots(...)          # ROC, scores, importance, confusion matrix
    def save_model(self)             # Save joblib + TMVA ROOT + metadata
    def run(self)                    # Full pipeline
```

### Training Data

- **Signal (label=1)**: All clusters from jet MC (`shapes_jet10.txt` through `shapes_jet50.txt`), filtered by `pid != -1`
- **Background (label=0)**: NPB-tagged clusters from data (`shapes_data_npb.txt`), filtered by `pid == -1`
- Kinematic cuts: 6 <= ET <= 40 GeV, |eta| <= 0.7

### NPB Feature Set (25 features, in order)

```
cluster_Et, cluster_Eta, vertexz,
e11_over_e33, e32_over_e35, e11_over_e22, e11_over_e13,
e11_over_e15, e11_over_e17, e11_over_e31,
e11_over_e51, e11_over_e71, e22_over_e33,
e22_over_e35, e22_over_e37, e22_over_e53,
cluster_weta_cogx, cluster_wphi_cogx,
cluster_et1, cluster_et2, cluster_et3, cluster_et4,
cluster_w32, cluster_w52, cluster_w72
```

This feature order must exactly match `config_npb_training.yaml: features.feature_list` and the `x_npb` vector constructed in `apply_BDT.C` (lines 451-477).

### NPB XGBoost Parameters (config_npb_training.yaml)

| Parameter | Value |
|-----------|-------|
| n_estimators | 750 |
| max_depth | 6 |
| learning_rate | 0.1 |
| subsample | 0.8 |
| colsample_bytree | 0.8 |
| reg_alpha | 1.0 |
| reg_lambda | 1.0 |
| tree_method | hist |
| max_bin | 256 |
| objective | binary:logistic |
| eval_metric | logloss |

### NPB Model Files

- `npb_models/npb_score_tmva.root` -- standard TMVA model
- `npb_models/npb_score_split_tmva.root` -- split variant
- `npb_models/npb_score.joblib` -- Python joblib
- `npb_models/npb_score_metadata.yaml` -- training metadata

### NPB Application in apply_BDT.C

- Applied when model file exists on disk
- Phase-space gate: `6 <= ET <= 40 GeV` AND `|eta| <= 0.7`
- Outside gate: `cluster_npb_score = -1`
- Uses same `TMVA::Experimental::RBDT` pattern as photon-ID models
- Model path selection: `use_split_bdt` config flag switches between split/standard variants

---

## 11. Input/Output Data Contract

### Data Paths

**Training Input (txt files):**
- Standard: `FunWithxgboost/NOSPLIT/shapes_{sample}.txt` (sample = photon5, photon10, photon20, jet10, jet15, jet20, jet30, jet50)
- Split: `FunWithxgboost/shapes_split_{sample}.txt`
- NPB data: `FunWithxgboost/shapes_split_data_npb.txt`

**Sim ROOT files (input to BDTinput.C and apply_BDT.C):**
- Path pattern: `{photon_jet_file_root_dir}/{filetype}/{photon_jet_file_branch_dir}`
- Default: `/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/{filetype}/condorout/combined.root`

**Data ROOT files:**
- `/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/condorout/part_*_with_bdt_split.root`

**BDT Model Output:**
- `FunWithxgboost/binned_models/model_{variant}[_split]_single_tmva.root`
- `FunWithxgboost/npb_models/npb_score[_split]_tmva.root`

**BDT-Scored Trees (output of apply_BDT.C):**
- Sim: `FunWithxgboost/{filetype}/bdt_split.root`
- Data: `{original_path}_with_bdt_split.root`

### Branch Name Patterns

All cluster branches follow: `{branch}_{clusternodename}` where clusternodename = `CLUSTERINFO_CEMC`.

**Branches written by apply_BDT.C:**
- `cluster_bdt_CLUSTERINFO_CEMC_base[ncluster_CLUSTERINFO_CEMC]/F`
- `cluster_bdt_CLUSTERINFO_CEMC_base_vr[ncluster_CLUSTERINFO_CEMC]/F`
- `cluster_bdt_CLUSTERINFO_CEMC_base_v0[ncluster_CLUSTERINFO_CEMC]/F`
- `cluster_bdt_CLUSTERINFO_CEMC_base_v1[ncluster_CLUSTERINFO_CEMC]/F`
- `cluster_bdt_CLUSTERINFO_CEMC_base_v2[ncluster_CLUSTERINFO_CEMC]/F`
- `cluster_bdt_CLUSTERINFO_CEMC_base_v3[ncluster_CLUSTERINFO_CEMC]/F`
- `cluster_bdt_CLUSTERINFO_CEMC_base_E[ncluster_CLUSTERINFO_CEMC]/F`
- `cluster_bdt_CLUSTERINFO_CEMC_base_v0E[ncluster_CLUSTERINFO_CEMC]/F`
- `cluster_bdt_CLUSTERINFO_CEMC_base_v1E[ncluster_CLUSTERINFO_CEMC]/F`
- `cluster_bdt_CLUSTERINFO_CEMC_base_v2E[ncluster_CLUSTERINFO_CEMC]/F`
- `cluster_bdt_CLUSTERINFO_CEMC_base_v3E[ncluster_CLUSTERINFO_CEMC]/F`
- `cluster_npb_score_CLUSTERINFO_CEMC[ncluster_CLUSTERINFO_CEMC]/F`

---

## 12. Constants and Hardcoded Values

### BDTinput.C
- `TIME_SAMPLE_NS = 17.6` -- EMCal time sample period in nanoseconds
- `cluster_Et` cut: `6 < ET < 40 GeV` (hardcoded, lines 582-585)
- `cluster_Eta` cut: `|eta| < 0.7` (hardcoded, line 587)
- Skip runs: `{47698, 51489, 51721, 51725, 53284}` (line 472)
- Skip event: `16152886` (line 505)
- NPB back-to-back jet pT threshold: `5 GeV` (line 630)
- NPB back-to-back dphi threshold: `pi/2` (line 636)
- Cluster array size for NPB timing: `49` towers per cluster (line 602)
- yaml-cpp path: `/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so`
- Lumi file: `/sphenix/user/jocl/projects/LumiList/list_468_RN_18_26_30_2.list`

### apply_BDT.C
- `nclustercontainermx = 4096` -- maximum clusters per event (line 115)
- BDT ET threshold: `cluster_Et > 7 GeV` (line 396) -- below this, score = -1
- NPB ET range: `6 <= ET <= 40 GeV` (lines 64-65)
- NPB eta range: `|eta| <= 0.7` (line 66)
- Model names list: hardcoded 11-element vector (lines 73-85)
- Model file pattern: `"binned_models/model_" + name + model_suffix + "_single_tmva.root"` (line 89)
- NPB model path: `"npb_models/npb_score[_split]_tmva.root"` (lines 40-42)

### Training (Python)
- `global_seed = 42` -- used throughout for reproducibility
- Eta reweighting range: `[-0.7, 0.7]` (hardcoded in reweighting.py line 49)
- PDF floor for inverse-PDF reweighting: `1e-3` (reweighting.py lines 58, 98)
- TMVA feature names: `f0, f1, f2, ...` (positional, not named)

### Cross-Section Weights (CrossSectionWeights.h)
- `photon5cross = 146359.3` pb
- `photon10cross = 6944.675` pb
- `photon20cross = 130.4461` pb
- `jet10cross = 3.997e+06` pb
- `jet15cross = 4.073e+05` pb
- `jet20cross = 6.2623e+04` pb
- `jet30cross = 2.5298e+03` pb
- `jet50cross = 7.3113` pb
- `default_nsimevents = 1e7`

---

## 13. Gotchas and Sync Requirements

### Feature Order Synchronization (CRITICAL)

The most critical sync requirement is that the feature ordering in the training config MUST exactly match the feature vector construction in `apply_BDT.C`. Since TMVA models use positional indexing (`f0, f1, ...`), a mismatch silently produces wrong scores.

**Three-way sync chain:**
1. `make_split_configs.py: MODEL_FEATURES[variant]` -- defines training feature order
2. Training config `data.features` list -- used by `main_training.py` to train the model
3. `apply_BDT.C: x_list[variant_index]` -- constructs feature vector for inference

These three must be identical for each model variant. `make_split_configs.py` was explicitly designed with comments matching `apply_BDT.C x_list exactly`.

### Divide-by-Zero Handling Mismatch

- `BDTinput.C` uses `safe_div()`: checks `isfinite(num) && isfinite(den) && den != 0`
- `apply_BDT.C` uses `(num > 0) ? num / den : 0` -- checks numerator, not denominator
- This is functionally close but not identical. If `num > 0` but `den = 0`, `apply_BDT.C` would produce infinity while `BDTinput.C` would produce 0. In practice, negative cluster energies are rare.

### Model Variant Selection

The `use_split_bdt` flag in config_nom.yaml controls which model set is loaded:
- `use_split_bdt: 1` -> loads `model_{name}_split_single_tmva.root`
- `use_split_bdt: 0` -> loads `model_{name}_single_tmva.root`

Both sets must exist on disk for the configured variant.

### Background Subsampling

When `background_subsampling.enabled: true`, the low-ET background is subsampled before training. This affects class weights and the training distribution. The subsampling seed must be fixed (`random_seed: 42`) for reproducibility.

### config.yaml vs config_nom.yaml

These are DIFFERENT configs for DIFFERENT stages:
- `config.yaml` -- training config (features, XGB params, binning, reweighting)
- `config_nom.yaml` -- analysis/efficiency config (cuts, paths, ABCD regions, isolation thresholds)

`BDTinput.C` reads `config_nom.yaml` to apply cuts during feature extraction.
`main_training.py` reads `config.yaml` for the training pipeline.
`apply_BDT.C` reads `config_nom.yaml` for tree/node names and the `use_split_bdt` flag.

### 25-Feature vs Model Variant Feature Sets

The training config `config.yaml` lists 25 features under `data.features`. However, each model variant uses only a SUBSET of these features (6-11 features). The per-variant feature lists are defined in `make_split_configs.py: MODEL_FEATURES` and override the global list when generating split configs.

### NPB Score Integration

The NPB score model is trained separately from the photon-ID models (`train_npb_score.py` vs `main_training.py`) but applied in the same `apply_BDT.C` pass. The NPB feature list (25 features) is DIFFERENT from the photon-ID feature lists and includes features like `cluster_w32`, `cluster_w52`, `cluster_w72` that some photon-ID variants do not use.

### Luminosity File

BDTinput.C reads a luminosity file (`list_468_RN_18_26_30_2.list`) for per-run luminosity tracking, but this is for counting purposes only -- it does not affect the training data.

### make_split_configs.py Signal/Background Files

The split training uses a different set of MC samples than the NOSPLIT training:
- **Split signal**: `shapes_split_photon{5,10,20}.txt`
- **Split background**: `shapes_split_jet{5,12,20,30,40}.txt` (includes jet5 and jet40, finer binning)
- **NOSPLIT background**: `shapes_jet{10,15,20,30,50}.txt` (different jet samples)

This means the split and standard models are trained on different background samples.
