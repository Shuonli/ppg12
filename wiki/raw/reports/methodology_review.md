# Methodology Review: PPG12 vs World Standards

Systematic comparison of PPG12's analysis techniques against published isolated photon measurements from ATLAS, CMS, ALICE, and PHENIX. Each section identifies where PPG12 aligns with, departs from, or could improve upon established practice.

---

## 1. Isolation Definition

### PPG12 Implementation

- **Cone radius:** R = 0.3 (truth), R = 0.4 topological clusters (reco, `use_topo_iso: 2`)
- **Threshold:** Parametric, `iso_max = 0.502 + 0.0433 * ET` GeV
- **MC correction:** `recoisoET_sim = recoisoET * 1.2 + 0.1` (scale + shift)
- **Gap:** 0.8 GeV between isolated and non-isolated boundaries
- **Truth isolation:** `iso_ET_truth < 4 GeV`, cone R = 0.3

### Comparison with World Standards

| Aspect | ATLAS | CMS | ALICE (13 TeV) | PHENIX | PPG12 |
|--------|-------|-----|----------------|--------|-------|
| Cone R | 0.4 | 0.4 | 0.4 | 0.5 | 0.3 (truth) / 0.4 (reco) |
| Threshold | 4.8 + 0.0042*ET | ~5 GeV flat | 1.5 GeV (tracks) | 10% of E | 0.50 + 0.043*ET |
| UE subtraction | jet-area (rho*A) | PF + rho | none (kappa_iso) | none | MC scale+shift |
| Core exclusion | 0.125 x 0.175 rect | PF cluster removal | track-only (no core issue) | cone sum | topo-cluster self-subtraction |

### Assessment

**The R = 0.3 truth cone and Catani et al. concern.** PPG12's truth-level isolation uses R = 0.3, while the reco-level uses topological R = 0.4. This is a deliberate choice documented in the wiki: the small truth cone minimizes the measurement's sensitivity to underlying event modeling while the larger reco cone captures hadronic activity more completely. The concern from Catani, Fontannaz, Guillet, and Pilon about log(1/R) corrections at small R applies primarily to R < 0.2, where these logarithms can reach ~10% of the cross section at NLO. At R = 0.3, the log(1/R) correction is modest (log(1/0.3) ~ 1.2, compared to log(1/0.4) ~ 0.9 for ATLAS/CMS), and the JETPHOX NLO calculation used for theory comparison explicitly includes the R = 0.3 cone. **Verdict: The R = 0.3 truth cone is borderline but defensible, especially since JETPHOX handles it correctly. The reco R = 0.4 is standard. This is a reasonable design choice given sPHENIX's lower center-of-mass energy and smaller cluster multiplicity compared to LHC.**

**Parametric threshold.** PPG12's slope (0.043) is an order of magnitude larger than ATLAS's (0.0042). This reflects the very different energy regimes: PPG12 operates at 8-35 GeV where underlying event and pileup are a larger fraction of the photon ET, while ATLAS operates at 125-2500 GeV. At PPG12's nominal 20 GeV, the threshold is ~1.4 GeV; at ATLAS's 200 GeV, it is ~5.6 GeV. Both yield epsilon_h ~ 0.05-0.10 at their respective midpoints, so the fractional tightness is comparable. **Verdict: Well-motivated. The parametric form is now standard (ATLAS adopted it from 8 TeV onward). The steeper slope is appropriate for the lower ET range.**

**MC isolation correction.** The two-parameter correction (scale = 1.2, shift = 0.1 GeV) is functionally equivalent to ATLAS's data-driven peak matching but less transparent. ATLAS performs the correction in the isolation distribution itself before cutting, calibrating the MC isolation peak position to match data. PPG12 applies it as a simple affine transformation on the isolation variable in simulation. **Recommendation: Document the derivation of scale=1.2 and shift=0.1 more explicitly in the analysis note. Consider deriving these from a fit to the isolation peak position/width ratio between data and MC rather than fixed values. The analysis note (analysis.tex) does not currently describe how these correction values were determined.**

**Gap between iso and non-iso regions.** PPG12 uses 0.8 GeV; ATLAS (7 TeV) used 2 GeV. The narrower gap is motivated by the lower ET range and tighter isolation cut -- the isolation resolution at 8-35 GeV is smaller in absolute terms than at LHC energies. **Verdict: Acceptable, but the systematic variation of this gap should be included, and it is (via `noniso_gap_tight/loose` variations).**

---

## 2. ABCD Implementation

### PPG12 Implementation

- **Axis 1 (Isolation):** Parametric ET-dependent threshold (described above)
- **Axis 2 (Identification):** BDT score with parametric threshold `bdt_min = 0.80 - 0.015 * ET`
- **Signal leakage:** Computed from signal MC (photon5/10/20, cross-section weighted)
  - cB, cC, cD per ET bin
  - Propagated via Gaussian smearing of leakage fractions in 20,000 Poisson toy MC
- **R_bg:** Computed two ways (data-driven and MC truth-unmatched), but set to R = 1 in nominal analysis
- **Purity solving:** Root-finding on the implicit ABCD equation via `TF1::GetX()`
- **Purity smoothing:** Pade [1/1] rational `(a + bx)/(1 + cx)` or error function, with 68.3% CL confidence intervals
- **MC purity correction:** Optional closure correction (truth/ABCD ratio from MC)

### Comparison with World Standards

**Signal leakage corrections: PPG12 is fully compliant with the standard.** The implementation in `CalculatePhotonYield.C` (lines 199-218, 384-465) correctly:
1. Computes cB, cC, cD from signal MC per bin
2. Solves the implicit ABCD equation with all three leakage terms
3. Propagates MC statistical uncertainty on leakage via Gaussian sampling
4. Uses Poisson toy MC for data statistical uncertainty (accounting for weighted events via effective N = N^2/sum(w^2))

This is more sophisticated than ATLAS's analytic error propagation but equivalent in intent. The Poisson toy approach with effective event counts is a strength -- it correctly handles the non-Poisson statistics from weighted MC events.

**BDT score as the tight/non-tight axis: A departure from ATLAS, well-motivated but with a caveat.** ATLAS deliberately chooses which shower-shape variables to invert for the non-tight definition, specifically selecting first-layer variables that have minimal correlation with calorimetric isolation. PPG12 uses the BDT score, which combines all shower-shape features. The key question is whether the BDT score is correlated with isolation for background events.

Evidence that PPG12 handles this:
- The BDT is trained without isolation features (confirmed in `config.yaml` -- `recoisoET` is the label, not a feature)
- The isolation variable is not in the 25-feature list
- PPG12 computes R_bg both from data and MC to monitor deviations from the factorization assumption
- 23+ systematic variations are run to cover ABCD parameter space

**However, a concern remains:** Even though isolation is not a BDT input, the BDT score can develop implicit correlations with isolation through intermediate variables. For example, clusters in jets tend to have both larger isolation energy AND broader shower shapes. The BDT might implicitly learn a proxy for isolation via the shower shape variables. This is precisely the correlation that ATLAS avoids by selectively inverting only first-layer variables. **Recommendation: Quantify the BDT-isolation correlation for background events explicitly (e.g., profile plot of mean BDT score vs isolation energy in the non-tight region, or a Pearson correlation coefficient). If the correlation is significant, consider either (a) applying ATLAS's approach of validating R_bg via control region subdivision, or (b) training a decorrelated BDT (adversarial or distance correlation penalty).**

**R_bg = 1 assumption:** PPG12 currently uses R = 1 in the nominal analysis. Both the data-driven and MC-truth R values are computed and monitored, but deviations are propagated only as a systematic. ATLAS (13 TeV, 1908.02746) found R_bg deviations of 15-40% depending on ET and eta, with impact on the cross section of up to 3.6%. **Recommendation: Explicitly report the measured R_bg values per ET bin in the analysis note. If deviations are greater than ~10%, consider using the measured R_bg rather than unity as the nominal value, with the deviation from unity as the systematic.**

**Purity smoothing:** PPG12 smooths per-bin purity values with a functional fit (Pade [1/1] or error function), then uses the fit values for purity correction. ALICE (13 TeV) similarly uses sigmoid smoothing. This is a sound approach for reducing statistical fluctuations in the purity correction, especially at low ET where statistics are limited. The confidence interval propagation via `TVirtualFitter::GetConfidenceIntervals()` at 68.3% CL is standard. **Verdict: Good practice.**

**MC purity correction:** PPG12 implements a closure-based correction where the ABCD purity extracted from MC is compared to truth purity, and the ratio is applied to data. This is a unique feature not used by ATLAS or ALICE. It effectively corrects for violations of the ABCD independence assumption captured by the simulation. **Verdict: This is a valuable innovation, analogous to ALICE's alpha_MC correction for background correlation. Ensure it is clearly documented as an optional correction with its own systematic (on/off toggle).**

---

## 3. Unfolding

### PPG12 Implementation

- **Method:** D'Agostini iterative Bayesian (RooUnfoldBayes)
- **Default iterations:** 2 (`resultit: 2` in config)
- **Scanned iterations:** 1-10 (all stored)
- **Reco bins:** [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36] GeV (12 bins)
- **Truth bins:** [7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45] GeV (14 bins)
- **Response matrix reweighting:** Currently disabled (`reweight = 0`, hardcoded override at line 836)
- **Efficiency correction:** Applied after unfolding (correct methodology)
- **Closure tests:** Full closure and half closure documented in analysis note

### Comparison with World Standards

**Iteration count.** PPG12 uses 2 iterations. The D'Agostini recommendation (arXiv:1010.0632) and general practice is 3-5 iterations. ATLAS typically uses 4 iterations. The analysis note (Section "Unfolding", Figure unfold_iter) shows that relative deviations between successive iterations become smaller than statistical uncertainty after iteration 1, which supports using 2 iterations. **Verdict: 2 iterations is on the low side of standard practice but justified by the analysis-specific convergence study. The relatively good energy resolution of the sPHENIX CEMCal (sigma_E/E ~ 2.8% + 15.5%/sqrt(E)) means less smearing and faster convergence. The 2-10% unfolding correction mentioned in the analysis note confirms that detector smearing is a small effect.**

**Closure tests.** PPG12 performs both full and half closure tests:
- Full closure: reproduces truth within 1%
- Half closure: agrees within 2%

This is exactly the standard practice. ATLAS and CMS both perform analogous tests. **Verdict: Fully compliant.**

**Non-square response matrix.** PPG12 uses 12 reco bins and 14 truth bins, with truth extending beyond reco on both sides (underflow [7,8] GeV and overflow [36,45] GeV). This is correct practice to handle bin migration at the edges, consistent with RooUnfold documentation and ATLAS/CMS usage. **Verdict: Good.**

**Response matrix reweighting.** The code contains a full implementation of response matrix reweighting (data/MC spectrum ratio applied to the response), but this is currently hardcoded off (`reweight = 0` at line 836, overriding the config value). The analysis note describes the reweighting as disabled in the nominal configuration but shows the data/MC ratio fitted by a Pade [2/2] approximant. **Recommendation: The reweight toggle should be controlled by the config, not hardcoded. A systematic variation with reweight on vs off should be included to quantify prior dependence.**

**No SVD cross-check.** PPG12 exclusively uses Bayesian unfolding. The wiki article on unfolding mentions SVD as a cross-check, and it is standard practice at ATLAS and CMS to run both Bayesian and SVD (or TUnfold/matrix inversion) and compare. **Recommendation: Add an SVD cross-check (RooUnfoldSvd) as validation. This requires minimal code changes -- swap `RooUnfoldBayes` for `RooUnfoldSvd` with appropriate regularization parameter k.**

**Error propagation.** PPG12 uses the default RooUnfold error propagation (full covariance with iteration-dependent correction, per Adye). This is the recommended approach. **Verdict: Standard.**

---

## 4. BDT Photon Identification

### PPG12 Implementation

- **Algorithm:** XGBoost, 750 trees, max_depth=5, learning_rate=0.1
- **Features:** 25 (shower shapes + kinematic + vertex position)
- **Training data:** Photon MC (signal) vs jet MC (background), 6-40 GeV, |eta| < 0.7
- **Reweighting:** Class + ET + eta reweighting (vertex z disabled)
- **ET binning:** 5 bins [6, 10, 15, 20, 25, 35] GeV for training; 2 bins [8, 15, 35] for inference
- **Threshold:** Parametric `bdt_min = 0.80 - 0.015 * ET`
- **11 model variants:** Different feature subsets for systematic cross-checks
- **Inference:** TMVA::Experimental::RBDT in ROOT C++

### Comparison with World Standards

**25 features vs other experiments.** PPG12 uses the most extensive feature set of any isolated photon measurement:
- ATLAS: ~9 shower shape variables (rectangular cuts, no multivariate)
- CMS: 1 primary (sigma_eta_eta), with BDT for 13 TeV using a few more
- ALICE: 1 (sigma_long^2 eigenvalue)
- D0: 4 (neural network)
- PPG12: 25 (XGBoost BDT)

**Are all 25 features necessary?** The features fall into groups:
1. **Core discrimination** (6 features): e11/e33, et1-et4, e32/e35 -- these are the sPHENIX analogs of ATLAS R_eta and CMS sigma_eta_eta
2. **Width variables** (5 features): weta_cogx, wphi_cogx, w32, w52, w72 -- transverse shower profile at multiple scales
3. **Extended energy ratios** (10 features): e11/e22, e11/e13, ..., e22/e53 -- multi-scale asymmetry probes
4. **Kinematic** (3 features): cluster_Et, cluster_Eta, vertexz -- enable position-dependent boundaries
5. **Cross ratios** (1 feature): e11_over_e22 -- additional core concentration measure

The 11 model variants with different feature subsets provide evidence for which features matter: the "base" model uses only 7 features (vertexz, Eta, e11/e33, et1-4) and already provides strong discrimination. The extended features provide incremental improvements. **Verdict: 25 features is on the high end and raises a question about overfitting, but the extensive regularization (L1 alpha=5.0, L2 lambda=0.3, subsample=0.5, colsample=0.6) and the 11 variant cross-checks mitigate this concern. The approach is reasonable for sPHENIX's coarser granularity, where no single variable achieves the discrimination power of ATLAS's fine strip layer.**

**ET-dependent threshold.** The parametric threshold `bdt_min = 0.80 - 0.015 * ET` decreases with ET, reflecting the fact that photon and jet showers become more similar at higher ET (pi0 decay photons merge). This is analogous to ALICE's ET-dependent sigma_long^2 upper bound (0.4 at 10 GeV, 0.3 above 16 GeV). ATLAS uses fixed shower-shape cuts without ET dependence for the individual variables, relying instead on the increasing purity at high ET to maintain performance. **Verdict: Well-motivated. The ET-dependent threshold is appropriate for PPG12's lower ET range where the signal/background discrimination changes rapidly with energy.**

**Including ET as a BDT feature.** The "base_E" and all "_E" variants include `cluster_Et` as a BDT input feature. This is a deliberate design choice: it allows the BDT to learn ET-dependent discrimination boundaries internally, rather than requiring separate models per ET bin. However, this can introduce a subtle issue: the BDT score becomes explicitly ET-dependent, which means the tight/non-tight boundary is driven by both the ET-dependent BDT score AND the ET-dependent parametric threshold. **Recommendation: Verify that the ET-dependence of the BDT score does not introduce a spurious correlation with isolation through the ET variable itself. A BDT that strongly depends on ET could effectively become a proxy for the isolation cut (since isolation energy scales with ET).**

**Training data composition.** PPG12 trains on photon MC (pid in {1,2}, direct + fragmentation photons) vs jet MC (pid not in {1,2}). This is standard and correct. The class/ET/eta reweighting ensures balanced training despite the very different cross sections of signal and background. **Verdict: Good practice.**

**TMVA inference.** The export to TMVA::Experimental::RBDT format for C++ evaluation is a pragmatic choice for integration with the ROOT-based analysis chain. Feature ordering must match exactly between training and inference, which is a known fragility. **Recommendation: Add an automated sync check between config.yaml features and apply_BDT.C feature vectors (or reference the existing sync-check skill).**

---

## 5. Potential Improvements

### High Priority

1. **Quantify BDT-isolation correlation for background.** The BDT-based ABCD method's validity rests on the assumption that BDT score and isolation are uncorrelated for background. While the BDT is not trained on isolation, implicit correlations through shower-shape features are possible. Produce profile plots of BDT score vs isolation energy for background-enriched samples (non-tight region, or truth-unmatched MC). If correlation is significant (Pearson |rho| > 0.1), either (a) validate R_bg more rigorously via ATLAS-style control region subdivision, or (b) train a decorrelated BDT.

2. **Report R_bg values per bin.** The nominal R = 1 assumption should be accompanied by the measured values. ATLAS found 15-40% deviations. If PPG12's deviations are similar, using R = 1 introduces a bias that should either be corrected or conservatively bounded.

3. **Add SVD unfolding cross-check.** This is minimal effort (swap RooUnfoldBayes for RooUnfoldSvd) and is standard practice at ATLAS and CMS. The comparison validates that the result is not sensitive to the choice of regularization scheme.

4. **Fix hardcoded reweight=0 override.** Line 836 of CalculatePhotonYield.C overrides the config value for response matrix reweighting. This should be controlled by the config file for reproducibility and systematic variation purposes.

### Medium Priority

5. **Document MC isolation correction derivation.** The analysis note should explain how `mc_iso_scale = 1.2` and `mc_iso_shift = 0.1` were determined. Consider deriving them from a fit to the isolation peak position and width in data vs MC.

6. **Purity fit function stability.** PPG12 offers two fit functions (error function and Pade [1/1]) and uses `fit_option` to select between them. Consider including the fit function choice as a systematic variation (the difference between erf and Pade fits propagated as a purity systematic).

7. **Response matrix reweighting as systematic.** The reweighting infrastructure exists in the code but is disabled. Run the analysis with reweight on and off and include the difference as an unfolding systematic.

8. **Feature importance monitoring.** Track and report the top-5 feature importances from the BDT training. If kinematic variables (ET, eta, vertexz) dominate over shower-shape variables, the BDT may be learning kinematic correlations rather than shower-shape discrimination.

### Low Priority / Future Improvements

9. **Frixione smooth-cone isolation.** For future NLO/NNLO comparisons, consider whether a smooth-cone (Frixione) isolation could be adopted at truth level. This eliminates fragmentation contributions entirely and is increasingly the standard at ATLAS for NNLO comparisons. At PPG12's energies and with JETPHOX, the standard cone with R = 0.3 is adequate, but for future NNLOJET comparisons, smooth cone isolation would reduce theoretical uncertainty.

10. **Adversarial BDT training for decorrelation.** If the BDT-isolation correlation (item 1) is found to be significant, an adversarial training objective (e.g., distance correlation penalty) could be added to the XGBoost training to explicitly minimize the correlation between BDT score and isolation energy for background events.

11. **Template fit cross-check.** ATLAS and CMS both perform template fits (of the isolation distribution or the shower-shape distribution) as an independent cross-check of the ABCD purity. PPG12 could fit the BDT score distribution in the isolated region to signal and background templates extracted from MC (signal) and the non-isolated sideband (background). Agreement with ABCD purity would strengthen confidence in the result.

12. **Covariance matrix for unfolded result.** The current analysis provides bin-by-bin uncertainties on the unfolded spectrum. For future NLO/NNLO chi2 comparisons, the full covariance matrix (from RooUnfold) should be stored and used. This is important because unfolding introduces inter-bin correlations.

---

## Summary Table

| Technique | PPG12 Status | Standard Compliance | Key Concern |
|-----------|-------------|-------------------|-------------|
| Parametric isolation | b + s*ET, R=0.3 truth / R=0.4 reco | Good (matches ATLAS parametric form) | R=0.3 truth borderline for NLO log(R) corrections; handled by JETPHOX |
| ABCD method | BDT-based, full leakage corrections | Good (signal leakage, toy MC errors) | BDT-isolation correlation not explicitly quantified; R_bg = 1 assumed |
| Bayesian unfolding | 2 iterations, full/half closure | Good (convergence justified) | No SVD cross-check; reweighting hardcoded off |
| BDT photon ID | 25 features, ET-dependent threshold | Good (most extensive feature set) | Potential implicit correlation with isolation via shower-shape proxy |
| Purity smoothing | Pade fit with confidence intervals | Good (matches ALICE sigmoid approach) | Fit function choice not varied as systematic |
| MC purity correction | Closure-based correction factor | Innovative (beyond ATLAS/CMS practice) | Ensure clear documentation of on/off systematic |
| Systematic variations | 23+ configs, quadrature aggregation | Good (comprehensive coverage) | Missing SVD and reweighting variations |

---

## Overall Assessment

PPG12's methodology is generally sound and in many areas more sophisticated than comparable measurements at similar energies (ALICE, PHENIX). The BDT-based ABCD method, parametric isolation, Poisson toy MC for uncertainty propagation, and 23+ systematic variation pipeline represent a modern and thorough approach.

The primary areas requiring attention before publication are:
1. Explicit validation of the BDT-isolation decorrelation assumption for background events
2. Reporting of measured R_bg values alongside the R = 1 nominal choice
3. Addition of an SVD unfolding cross-check
4. Removal of the hardcoded reweight override in favor of config-driven control

These are achievable with moderate effort and would bring the analysis fully in line with the most stringent LHC publication standards.
