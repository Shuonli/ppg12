# BDT-Based Photon Identification

PPG12 uses a Boosted Decision Tree (BDT) trained with XGBoost to discriminate prompt photon clusters from hadronic background (pi0/eta decays, jet fragments) in the sPHENIX CEMCal. This article covers the XGBoost algorithm, TMVA inference, the PPG12 model configuration, and comparison with photon identification strategies at other experiments.

## XGBoost Training Algorithm

XGBoost (Chen & Guestrin 2016, arXiv:1603.02754) is a gradient-boosted tree ensemble. The prediction for an input **x** is the sum of K additive tree functions:

```
y_hat = sum_{k=1}^{K} f_k(x)
```

### Regularized Objective

The training minimizes a regularized loss:

```
L = sum_i l(y_i, y_hat_i) + sum_k [gamma * T_k + (1/2) * lambda * sum_j w_{kj}^2]
```

where `l` is the per-sample loss (logistic for binary classification), `T_k` is the number of leaves in tree k, `w_{kj}` are leaf weights, `gamma` penalizes tree complexity, and `lambda` is L2 regularization on leaf scores. PPG12 also uses L1 regularization (`reg_alpha`), which encourages sparsity in leaf weights.

### Split Finding

At each boosting step, the loss is Taylor-expanded to second order. The gain from splitting a node is:

```
Gain = (1/2) * [G_L^2/(H_L + lambda) + G_R^2/(H_R + lambda) - (G_L+G_R)^2/(H_L+H_R + lambda)] - gamma
```

where `G` and `H` are sums of first and second derivatives of the loss over instances in the left (L) and right (R) child nodes. A split is only made if `Gain > 0`, meaning `gamma` acts as a minimum-gain threshold.

### Anti-Overfitting Mechanisms

XGBoost provides multiple regularization levers:

| Mechanism | Parameter | PPG12 value | Effect |
|-----------|-----------|-------------|--------|
| L2 on leaf weights | `reg_lambda` | 0.3 | Shrinks leaf scores toward zero |
| L1 on leaf weights | `reg_alpha` | 5.0 | Drives some leaf scores to exactly zero |
| Learning rate | `learning_rate` | 0.1 | Scales each new tree's contribution |
| Row subsampling | `subsample` | 0.5 | Each tree sees 50% of training events |
| Column subsampling | `colsample_bytree` | 0.6 | Each tree uses 60% of features |
| Max tree depth | `max_depth` | 5 | Limits individual tree complexity |
| Leaf-wise growth | `grow_policy: lossguide` | -- | Splits the leaf with highest gain (not level-wise) |

### Histogram-Based Approximate Splitting

PPG12 uses `tree_method: hist`, which buckets continuous feature values into discrete bins (controlled by `max_bin: 256`). This replaces the exact greedy search over all possible split points with an approximate search over histogram boundaries, providing a large speedup with negligible accuracy loss.

## PPG12 Training Configuration

### Training Data

- **Signal:** Photon MC (photon5/10/20), clusters with `particle_photonclass in {1, 2}` (direct photons and fragmentation photons)
- **Background:** Jet MC (jet10/15/20/30/50), clusters with `particle_photonclass not in {1, 2}`
- **Selection:** 6 < ET < 40 GeV, |eta| < 0.7, vertex cuts, trigger requirement
- **Split:** 70% train / 10% validation / 20% test, stratified by label, `global_seed: 42`

### 25 Input Features

```
cluster_Et, cluster_weta_cogx, cluster_wphi_cogx, vertexz, cluster_Eta,
e11_over_e33, cluster_et1, cluster_et2, cluster_et3, cluster_et4,
e32_over_e35, cluster_w32, cluster_w52, cluster_w72,
e11_over_e22, e11_over_e13, e11_over_e15, e11_over_e17,
e11_over_e31, e11_over_e51, e11_over_e71,
e22_over_e33, e22_over_e35, e22_over_e37, e22_over_e53
```

These include shower shape variables (energy ratios, widths, ring fractions), kinematic variables (ET, eta), and vertex position. See [Shower Shape Discriminants](shower-shape-discriminants.md) for the physics of each variable.

**Feature order matters:** TMVA uses positional indexing (`f0, f1, f2, ...`), not named features. The order in `config.yaml` must exactly match the feature vectors in `apply_BDT.C`.

### Reweighting

Training samples are reweighted to equalize signal/background populations and flatten kinematic distributions:

1. **Class reweighting:** balances signal and background counts
2. **Eta reweighting:** inverse PDF in |eta| < 0.7 (per class)
3. **ET reweighting:** inverse PDF via cubic spline (per class), capped at 800.0
4. **Vertex z reweighting:** optional, currently disabled

Reweighting is applied globally before any ET binning.

### ET-Binned Models

PPG12 trains separate BDT models for ET bins [6, 10, 15, 20, 25, 35] GeV. The `bdt_et_bin_edges` and `bdt_et_bin_models` config fields control which model is applied to each cluster based on its reconstructed ET. This allows the BDT to adapt its discrimination strategy across the energy range where the signal-to-background ratio and shower shape distributions change significantly.

When `train_single_model: true` (current default), a single model is trained on all ET bins but evaluated per-bin. This provides a simpler model with slightly less per-bin optimization.

### 11 Model Variants

PPG12 maintains 11 BDT variants with different feature subsets:

| Model | Feature set (beyond vertexz, Eta, e11/e33, et1-4) |
|-------|----------------------------------------------------|
| base | (core set only) |
| base_v1 | + weta |
| base_v2 | + weta, wphi |
| base_v3 | + weta, wphi, e32/e35 |
| base_E | + Et |
| base_v0 | core set minus et1 |
| base_v0E | + Et, minus et1 |
| base_v1E | + Et, weta |
| base_v2E | + Et, weta, wphi |
| base_v3E | + Et, weta, wphi, e32/e35 |
| base_vr | identical to base |

The variants provide systematic cross-checks: dependence on the feature set is a source of BDT systematic uncertainty.

## TMVA Inference

### Model Export

Trained XGBoost models in Python are exported to TMVA-compatible ROOT format:

```
binned_models/model_{variant_name}_single_tmva.root
```

The export maps each XGBoost tree ensemble to the `TMVA::Experimental::RBDT` format within ROOT. Feature names are replaced by positional indices (`f0, f1, ...`).

### C++ Application (apply_BDT.C)

In `apply_BDT.C`, the RBDT model is loaded and evaluated per cluster:

```cpp
TMVA::Experimental::RBDT bdt("bdt_model", "binned_models/model_base_single_tmva.root");
// ... fill feature vector ...
auto score = bdt.Compute({f0, f1, f2, ...});
```

Scores are computed only for clusters with `cluster_Et > 7 GeV`; below this threshold the score is set to -1. Output branches follow the naming convention `cluster_bdt_{nodename}_{modelname}[ncluster]`.

### NPB Score Model

A separate XGBoost model discriminates real physics clusters from non-physics background (beam-induced timing artifacts):

- Signal: jet MC clusters (real physics)
- Background: NPB-tagged data clusters
- 25 features (different order from photon-ID BDT)
- Applied for 6 < ET < 40 GeV, |eta| < 0.7
- Output: `cluster_npb_score_{nodename}[ncluster]`

## ET-Dependent BDT Threshold

The tight/non-tight classification boundary is not a flat BDT score cut. PPG12 uses a parametric threshold:

```
tight_bdt_min = intercept + slope * ET
```

where `intercept` and `slope` are config parameters. This accounts for the ET dependence of the BDT score distribution: at higher ET, photon and jet showers become more similar, and the optimal working point shifts.

## Comparison with Other Experiments

### ATLAS: Rectangular Shower Shape Cuts

ATLAS (arXiv:1012.4389) uses a set of rectangular cuts on individual shower shape variables from a finely segmented liquid-argon calorimeter:

- **Strip layer variables:** lateral width (w_s,tot), energy asymmetry (E_ratio), side energy fraction (F_side), 3-strip width (w_s,3)
- **Middle layer variables:** R_eta (3x7/7x7 energy ratio), R_phi (3x3/3x7), lateral width (w_2), hadronic leakage (R_had)
- Classification: "tight" requires passing all cuts; "non-tight" fails a subset of strip-layer criteria

This approach is transparent and well-understood but does not exploit correlations between variables. ATLAS's fine granularity (strip layer with 0.003 x 0.1 segmentation) provides powerful standalone discrimination without multivariate methods.

### CMS: sigma_eta_eta Template Fit

CMS (arXiv:1108.2044) uses the second moment of the cluster energy distribution in eta:

- **sigma_eta_eta:** computed from the 5x5 crystal matrix around the seed
- Barrel cut: sigma_eta_eta < 0.010; endcap: sigma_eta_eta < 0.028--0.030
- Purity extracted via template fits to isolation distributions or sigma_eta_eta sidebands
- CMS's PbWO4 crystal calorimeter with 0.0174 x 0.0174 granularity provides excellent intrinsic shower shape resolution

CMS's sigma_eta_eta is a single scalar capturing the transverse shower width -- analogous to the weta_cogx feature used in the sPHENIX BDT. The BDT approach used by PPG12 combines the information from many such variables simultaneously.

### ALICE: sigma_long^2 Eigenvalue

ALICE (arXiv:1906.01371) uses the larger eigenvalue of the cluster energy distribution in the eta-phi plane:

```
sigma_long^2 = (sigma_phi_phi + sigma_eta_eta)/2 + sqrt((sigma_phi_phi - sigma_eta_eta)^2/4 + sigma_eta_phi^4)
```

- Narrow cluster selection: 0.1 < sigma_long^2 < sigma_max (ET-dependent upper bound, 0.3--0.4)
- ALICE EMCal is a lead-scintillator sampling calorimeter with cell size 6 x 6 cm^2 (Delta_eta x Delta_phi ~ 0.014 x 0.014), comparable to sPHENIX CEMCal

This is the closest detector analog to sPHENIX. The sigma_long^2 eigenvalue combines eta and phi width information -- similar to using both weta_cogx and wphi_cogx. ALICE's ET-dependent sigma_max threshold parallels PPG12's ET-dependent BDT score threshold.

### D0: Neural Network

D0 (hep-ex/0511054) pioneered multivariate photon identification using a 4-input artificial neural network (JETNET):

1. Number of EM1 cells with E > 400 MeV within R < 0.2
2. Number of EM1 cells with E > 400 MeV within 0.2 < R < 0.4
3. Scalar sum of track pT within 0.05 < R < 0.4
4. Energy-weighted cluster width in the finely segmented EM3 layer

The NN output O_NN is cut at 0.5, yielding 93.7% signal efficiency. Purity is extracted via template fits of the O_NN distribution.

D0's approach is the historical precursor to PPG12's BDT: both use a multivariate classifier trained on MC to separate photons from hadrons, with the classifier score used to define tight/non-tight regions in the ABCD method. The main advances in PPG12 are the much larger feature set (25 vs 4), gradient boosting instead of a shallow neural network, and ET-binned models.

### Summary of Approaches

| Experiment | Detector | Method | Variables | Purity extraction |
|-----------|----------|--------|-----------|-------------------|
| ATLAS | LAr (fine strips) | Rectangular cuts | ~9 shower shape | ABCD with tight/non-tight |
| CMS | PbWO4 crystal | sigma_eta_eta cut | 1 primary | Template fit (isolation) |
| ALICE | Pb-scintillator | sigma_long^2 cut | 1 eigenvalue | ABCD with narrow/wide |
| D0 | LAr/uranium | Neural network (4 inputs) | 4 | Template fit (NN output) |
| sPHENIX (PPG12) | W/SciFi | XGBoost BDT (25 features) | 25 shower shape + kinematic | ABCD with tight/non-tight |

PPG12's BDT approach is the most multivariate of these methods. The W/SciFi calorimeter has coarser granularity than ATLAS or CMS but similar granularity to ALICE EMCal, making the multivariate approach important for extracting maximum discrimination from the available tower information.

## Key References

| Reference | Role |
|-----------|------|
| Chen & Guestrin, arXiv:1603.02754 | XGBoost algorithm |
| Hoecker et al., physics/0703039 | TMVA framework for BDT inference |
| ATLAS, arXiv:1012.4389 | ABCD with rectangular shower shape cuts |
| CMS, arXiv:1108.2044 | sigma_eta_eta template method |
| ALICE, arXiv:1906.01371 | sigma_long^2 ABCD method |
| D0, hep-ex/0511054 | Neural network photon ID |

## See Also

- [Stage 2: BDT Training](../../pipeline/02-bdt-training.md) -- full pipeline details, feature extraction, model export
- [Shower Shape Discriminants](shower-shape-discriminants.md) -- physics of the 25 input features
- [Shower Shape Variables (concept)](../../concepts/shower-shape-variables.md) -- variable definitions and computation
