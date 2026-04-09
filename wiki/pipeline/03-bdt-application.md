# Stage 3: BDT Application (apply_BDT.C)

## Key File

`FunWithxgboost/apply_BDT.C` (492 lines)

## Function Signature

```cpp
void apply_BDT(const std::string &configname = "config_nom.yaml",
               const std::string filetype = "jet12_double",
               const std::string inputfilename = "...")
```

## What It Does

1. Opens the input slimtree ROOT file
2. Loads all 11 photon-ID TMVA models + the NPB model
3. Clones the tree and adds new BDT score branches
4. For each event/cluster, constructs the feature vector for each model and calls `RBDT.Compute()`
5. Writes the scored tree to a new ROOT file

## ET Threshold

BDT scores are computed only when `cluster_Et > 7 GeV`; otherwise the score is set to -1.

## 11 Model Variants and Feature Vectors

| # | Name | Features (exact order) | N |
|---|------|----------------------|---|
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

**Key differences between variants:**
- `base` vs `base_vr`: identical features (vr = vertex-reweighted training variant)
- `v0`: drops `et1` (most correlated with isolation)
- `v1/v2/v3`: progressively add weta, wphi, e32/e35
- `_E` suffix: adds `cluster_Et` as first feature (energy-aware)

## Derived Variables

Computed in apply_BDT.C for model input:
```cpp
e11_over_e33_BDT = (cluster_e11 > 0) ? cluster_e11 / cluster_e33 : 0;
e32_over_e35_BDT = (cluster_e32 > 0) ? cluster_e32 / cluster_e35 : 0;
```

## Output Branches

For each of the 11 models:
```
cluster_bdt_{clusternodename}_{model_name}[ncluster_{clusternodename}]/F
```
Example: `cluster_bdt_CLUSTERINFO_CEMC_base_v3E`

NPB score:
```
cluster_npb_score_{clusternodename}[ncluster_{clusternodename}]/F
```

## Output Files

- **Sim:** `{filetype}/bdt_split.root` (e.g., `FunWithxgboost/photon10/bdt_split.root`)
- **Data:** `{inputfilename}_with_bdt_split.root`

The output is a full clone of the input tree with BDT branches appended.

## Model Loading

```cpp
const bool use_split_bdt = configYaml["analysis"]["use_split_bdt"].as<int>(0) != 0;
const std::string model_suffix = use_split_bdt ? "_split" : "";
// File: "binned_models/model_" + name + model_suffix + "_single_tmva.root"
TMVA::Experimental::RBDT bdt("myBDT", model_file);
float score = bdt.Compute(feature_vector)[0];
```

## NPB Score Features (25, in order)

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

## Running

```bash
cd FunWithxgboost
# Apply to a sim sample
root -l -b -q 'apply_BDT.C("config_nom.yaml", "photon5")'

# Apply to data (uses glob pattern from config)
root -l -b -q 'apply_BDT.C("config_nom.yaml", "data")'

# Condor submission for all samples
condor_submit apply_BDT_condor.sub
```

## Gotchas

- Feature order in the C++ code must EXACTLY match the order used during training -- TMVA uses positional indexing
- When adding a new model variant, update both `main_training.py` (or its config) and `apply_BDT.C`
- The `use_split_bdt` flag selects between `model_{name}_single_tmva.root` and `model_{name}_split_single_tmva.root`
