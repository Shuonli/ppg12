# Stage 2: BDT Training

## Key Files

| File | Purpose |
|------|---------|
| `FunWithxgboost/BDTinput.C` | Feature extraction: slimtree to text files (1040 lines) |
| `FunWithxgboost/main_training.py` | Training orchestrator: BinnedTrainingPipeline (1295 lines) |
| `FunWithxgboost/data_loader.py` | Data loading, signal/background labeling (338 lines) |
| `FunWithxgboost/model_builder.py` | sklearn Pipeline construction (77 lines) |
| `FunWithxgboost/reweighting.py` | Kinematic reweighting: class, ET, eta, vertex (165 lines) |
| `FunWithxgboost/plotting.py` | Diagnostic plots (1222 lines) |
| `FunWithxgboost/config.yaml` | Training configuration |
| `FunWithxgboost/train_npb_score.py` | NPB score training (1073 lines) |
| `FunWithxgboost/config_npb_training.yaml` | NPB training config |
| `FunWithxgboost/make_split_configs.py` | Generate per-model-variant configs |

## Feature Extraction (BDTinput.C)

```cpp
void BDTinput(const std::string &configname = "config_nom.yaml",
              const std::string filetype = "data")
```

Reads slimtree, computes 13 energy ratios (e.g., `e11_over_e33 = cluster_e11 / cluster_e33`), and writes 30-column space-separated text files.

**Output files:** `shapes_split_{filetype}.txt` (e.g., `shapes_split_photon5.txt`)

**30 output columns:**
```
cluster_Et cluster_Eta cluster_Phi vertexz
[13 energy ratios] cluster_prob cluster_weta_cogx cluster_wphi_cogx
cluster_et1-4 cluster_w32 cluster_w52 cluster_w72
recoisoET is_tight pid
```

The last 3 columns (recoisoET, is_tight, pid) are metadata, not training features.

**Signal/background labeling:** `pid` = `particle_photonclass` value. In Python: signal = `pid in [1, 2]`, background = `pid not in [1, 2]`.

**Selection cuts:** 6 < ET < 40 GeV, |eta| < 0.7, |vertexz| < vertex_cut, MBD >= 1 hit per side, trigger fired (data), skip bad runs {47698, 51489, 51721, 51725, 53284}.

## Training Pipeline (main_training.py)

```
BinnedTrainingPipeline.run_training():
  1. Load config (YAML)
  2. Load data (DataLoader: txt files, label sig/bkg)
  3. Apply reweighting (KinematicReweighter, globally before binning)
  4. Analyze feature correlations
  5. Train model (single or per-bin)
  6. Evaluate (AUC, ROC)
  7. Generate plots
  8. Save results (joblib + TMVA ROOT)
```

### 25 Training Features (config.yaml, exact order)

```
cluster_Et, cluster_weta_cogx, cluster_wphi_cogx, vertexz, cluster_Eta,
e11_over_e33, cluster_et1, cluster_et2, cluster_et3, cluster_et4,
e32_over_e35, cluster_w32, cluster_w52, cluster_w72,
e11_over_e22, e11_over_e13, e11_over_e15, e11_over_e17,
e11_over_e31, e11_over_e51, e11_over_e71,
e22_over_e33, e22_over_e35, e22_over_e37, e22_over_e53
```

Feature ORDER matters: TMVA uses positional indexing (`f0, f1, f2, ...`), not named features.

### Reweighting Order

Applied globally before binning. The actual code order (corrected from earlier reports):

1. **Class reweighting** -- balance signal/background: `weight_sig = (n_sig + n_bkg) / (2 * n_sig)`
2. **Per class**, in order:
   - **Eta reweighting** -- inverse PDF in [-0.7, 0.7] (hardcoded range, ignores config `eta_reweight_range`)
   - **ET reweighting** -- inverse PDF via spline, cap at `et_reweight_max` (800.0)
   - **Vertex z reweighting** -- inverse PDF, cap at 50.0 (disabled by default)

### XGBoost Hyperparameters

| Parameter | Value |
|-----------|-------|
| n_estimators | 750 |
| max_depth | 5 |
| learning_rate | 0.1 |
| subsample | 0.5 |
| colsample_bytree | 0.6 |
| reg_alpha | 5.0 |
| reg_lambda | 0.3 |
| grow_policy | lossguide |
| tree_method | hist |
| objective | binary:logistic |
| eval_metric | auc |

### train_single_model

- `true` (default) = train ONE model for all ET bins, evaluate per-bin
- `false` = train separate model per ET bin

### Train/Val/Test Split

70/10/20%, stratified by label. Seeds: `global_seed` (42) for first split, `global_seed + 1` for second.

### ET Bins for Training

5 bins: `[6, 10, 15, 20, 25, 35]` GeV

## Model Export

Models are saved to `binned_models/` in two formats:
- Joblib: `model_{bin_label}.joblib` (Python inference)
- TMVA ROOT: `{prefix}_single_tmva.root` (C++ inference via `TMVA::Experimental::RBDT`)

Standard naming: `binned_models/model_{variant_name}_single_tmva.root`
Split naming: `binned_models/model_{variant_name}_split_single_tmva.root`

11 model variants exist: `base`, `base_vr`, `base_v0`, `base_v1`, `base_v2`, `base_v3`, `base_E`, `base_v0E`, `base_v1E`, `base_v2E`, `base_v3E`

Each uses a different subset of the 25 features -- see [BDT Application](03-bdt-application.md) for exact feature vectors.

## NPB Score Model

**Purpose:** Discriminate real physics clusters from non-physics background (timing artifacts).

- Signal (label=1): jet MC clusters (`pid != -1`)
- Background (label=0): NPB-tagged data clusters (`pid == -1`)
- 25 features (different order from photon-ID BDT)
- Output: `npb_models/npb_score_tmva.root`
- Applied only when 6 < ET < 40 GeV, |eta| < 0.7

## Gotchas

- The 25 features in `config.yaml` must stay in sync with feature vectors in `apply_BDT.C`
- Feature names are replaced with `f0, f1, ...` during TMVA export -- positional order is the contract
- Division protection differs between files: BDTinput.C uses `safe_div` (checks denominator != 0 and isfinite). apply_BDT.C has two patterns: the photon-ID block checks the numerator (`cluster_e11 > 0`), while the NPB score block checks the denominator (`cluster_e33 > 0`). Functionally equivalent for positive energies
- BDTinput.C normalizes jet weights relative to jet30cross (different from CrossSectionWeights.h which uses jet50cross) -- benign for feature extraction but a maintainability concern
