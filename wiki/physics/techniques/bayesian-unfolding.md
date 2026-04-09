# Bayesian Unfolding

Iterative Bayesian unfolding (D'Agostini method) is the primary technique PPG12 uses to correct the measured photon ET spectrum for detector smearing, bin migration, and acceptance effects. This article covers the algorithm, its implementation in RooUnfold, the SVD alternative, and how PPG12 applies it.

## The Unfolding Problem

A measured (reco-level) spectrum is distorted from the true (particle-level) spectrum by finite detector resolution, acceptance gaps, and reconstruction inefficiency. Formally, the measured distribution **b** is related to the true distribution **x** by a response matrix **A**:

```
b_j = sum_i A_ji * x_i + noise
```

where `A_ji = P(measured in bin j | true in bin i)` encodes the detector smearing. Direct matrix inversion amplifies statistical noise, producing wildly oscillating solutions. Unfolding methods regularize the inversion to obtain a physically meaningful result.

## D'Agostini Iterative Bayesian Method

### Algorithm

The method applies Bayes' theorem to individual cause-effect bin pairs. Given the response matrix `lambda_ji = P(E_j | C_i)`, a prior `P(C_i)`, and efficiency `epsilon_i = sum_j lambda_ji`, the unfolding matrix is:

```
theta_ij = P(C_i | E_j) = lambda_ji * P(C_i) / sum_k [lambda_jk * P(C_k)]
```

The unfolded spectrum is then:

```
x_hat(C_i) = (1 / epsilon_i) * sum_j [theta_ij * n(E_j)]
```

### Iteration as Regularization

The procedure iterates:

1. Start with an initial prior P(C_i) (e.g., MC truth distribution or flat)
2. Compute the unfolding matrix theta_ij from the prior and response matrix
3. Unfold the observed spectrum to obtain x_hat(C_i)
4. Use x_hat(C_i) (normalized) as the new prior
5. Repeat

**The number of iterations controls regularization.** Few iterations bias the result toward the prior but suppress statistical noise. Many iterations converge toward the unregularized (noisy) inverse. Typically 3--5 iterations suffice, and the result is often insensitive to the exact count within this range (D'Agostini 2010, arXiv:1010.0632).

### Convergence

Convergence can be monitored by:

- Chi2 between the re-folded result and observed data (should decrease then plateau)
- Relative change between consecutive iterations (should stabilize)
- Smoothness of the unfolded distribution (increasing oscillation signals over-unfolding)

### Uncertainty Propagation

The original 1995 algorithm used linear error propagation, which underestimates uncertainties for more than one iteration. The improved 2010 version models:

- Response matrix columns via **Dirichlet distributions**
- Observed counts via **Gamma distributions**
- Sharing among cause bins via **Multinomial sampling**
- Full covariance via Monte Carlo integration

This properly handles small-number statistics, finite MC statistics in the response matrix, non-Gaussian posteriors, and full inter-bin correlations.

**Key correction (Adye, RooUnfold):** The error propagation must include the dependence of the updated prior on the data. Without this second-order term, errors are underestimated for more than one iteration. RooUnfold implements this corrected propagation.

## RooUnfold Implementation

[RooUnfold](https://gitlab.cern.ch/RooUnfold/RooUnfold) (Adye 2011, arXiv:1105.1160) provides a common C++ framework within ROOT for multiple unfolding algorithms.

### Software Architecture

```
RooUnfoldResponse  -->  RooUnfold (base class)
  Fill(measured, truth)     |-- RooUnfoldBayes      (iterative Bayesian)
  Miss(truth)               |-- RooUnfoldSvd        (SVD regularization)
                            |-- RooUnfoldTUnfold    (Tikhonov regularization)
                            |-- RooUnfoldInvert     (direct matrix inversion)
                            |-- RooUnfoldBinByBin   (bin-by-bin correction factors)
```

The **response matrix** is an independent object trained from MC, serialized via ROOT I/O, and loaded for unfolding. `Fill(x_measured, x_true)` records migration pairs; `Miss(x_true)` registers events passing truth cuts but failing reconstruction (efficiency loss).

### RooUnfoldBayes Specifics

| Setting | Value |
|---------|-------|
| Initial prior | Training truth distribution (not flat) |
| Default iterations | 4 |
| Smoothing | None (avoids introducing bias) |
| Error propagation | Full covariance with iteration-dependent correction |

### Error Calculation Methods

RooUnfold offers three approaches:

1. **Bin-by-bin errors** -- diagonal only, no correlations (fast but incomplete)
2. **Full covariance from error propagation** -- propagates measurement errors through the unfolding
3. **MC toy covariance** -- generate toy datasets, unfold each, compute covariance from ensemble

Caveat: the covariance matrix can be poorly conditioned. For chi2 calculations involving the inverse, `TDecompSVD` is recommended over direct inversion.

## SVD Unfolding (Alternative)

SVD unfolding (Hocker and Kartvelishvili 1996, hep-ph/9509307) provides an independent regularization approach.

### Core Idea

Factor the response matrix via singular value decomposition: `A = U * S * V^T`. Small singular values correspond to high-frequency (oscillating) components that are amplified by inversion. Regularization suppresses these:

```
z_i^(tau) = d_i * s_i / (s_i^2 + tau)
```

where `tau` is the regularization parameter. For large `s_i` (well-determined components), the filter passes the signal unchanged. For small `s_i` (noise-dominated), the filter suppresses the contribution.

### Key Recommendations from the SVD Paper

- Use the **number-of-events matrix** (actual MC counts), not the probability matrix. This automatically weights well-determined bins more heavily.
- Normalize unknowns by the MC prior: `w_j = x_j / x_j^ini`. If the prior is close to truth, `w` is smooth and requires fewer SVD terms.
- Add a curvature penalty (second-derivative matrix **C**) to suppress oscillating solutions.
- Choose `tau` by plotting `|d_i|` vs component index `i`: the effective rank `k` is where `|d_i|` drops to the noise level (~1), and `tau = s_k^2`.

### Comparison: Bayesian vs SVD

| Aspect | Bayesian (RooUnfoldBayes) | SVD (RooUnfoldSvd) |
|--------|--------------------------|---------------------|
| Regularization control | Number of iterations | Cutoff parameter k |
| Interpretation | Iterative prior update | Frequency-domain filtering |
| Tuning | 3--5 iterations, insensitive | k must be optimized per case |
| Error propagation | Full covariance (corrected) | Built-in, includes MC stats |
| Covariance inverse | tau-dependent | tau-independent (advantage for chi2) |
| Implementation | Iterative (loop) | Single linear algebra pass |

Both methods address the same bias-variance tradeoff. Using both as cross-checks is standard practice.

## How PPG12 Uses Unfolding

### Response Matrix Construction

The response matrix is built in `RecoEffCalculator_TTreeReader.C` from signal MC samples (photon5, photon10, photon20). Only tight + isolated + truth-matched clusters contribute:

```cpp
responses_full[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
```

Events passing truth cuts but failing reconstruction are registered via `Miss()` to encode efficiency loss.

### Binning

- **Reco pT bins:** from config `pT_bins` (nominal: [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36] GeV, 12 bins)
- **Truth pT bins:** from config `pT_bins_truth` (nominal: [7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45] GeV)

Truth bins extend beyond reco bins on both sides to handle overflow and prevent edge bias. This is consistent with RooUnfold's support for non-square response matrices.

### Unfolding Call

In `CalculatePhotonYield.C`, the purity-corrected yield (from the ABCD method with signal leakage corrections) is unfolded:

```cpp
RooUnfoldBayes unfold(response, h_yield, n_iterations);
TH1D *h_unfolded = (TH1D *)unfold.Hreco();
```

The code runs iterations 1--10 and stores each as `h_unfold_sub_{1..10}`. The selected iteration is controlled by config field `analysis.unfold.resultit` (default: 2).

### Efficiency Correction

Efficiency is applied **after** unfolding (not folded into the response matrix):

```cpp
float total_eff = eff_reco * eff_iso * eff_id * vertex_eff;
```

This is the correct methodology -- efficiencies must remain separate from the response to avoid double-counting.

### Closure Tests

Two complementary closure tests validate the procedure:

- **Full closure:** Use entire MC reco distribution as pseudo-data, unfold with the same MC response, compare to truth. Tests for internal bias.
- **Half closure:** Build response from first 50% of MC events, unfold second 50%. Tests generalization (response and pseudo-data are statistically independent).

Good closure means the unfolded/truth ratio is centered at 1.0 across all pT bins, with RMS below 5% and chi2/ndf near 1.0. See `plotting/Closure.C`.

### Systematic Uncertainty

Unfolding systematics come from:

- **Iteration count variation** from the nominal value
- **Response matrix prior dependence** (reweighting studies)
- **SVD cross-check** (running RooUnfoldSvd alongside RooUnfoldBayes)

These are grouped as `unfolding` in the systematic aggregation pipeline.

## Key References

| Reference | Role |
|-----------|------|
| D'Agostini, NIM A 362 (1995) 487 | Original Bayesian unfolding algorithm |
| D'Agostini, arXiv:1010.0632 | Improved algorithm with MC uncertainty propagation |
| Adye, arXiv:1105.1160 | RooUnfold implementation and corrected error propagation |
| Hocker & Kartvelishvili, hep-ph/9509307 | SVD unfolding theory |

## See Also

- [Unfolding (pipeline concept)](../../concepts/unfolding.md) -- config fields, output histograms, closure test details
- [Stage 4: Efficiency and Yield](../../pipeline/04-efficiency-yield.md) -- where unfolding fits in the analysis chain
