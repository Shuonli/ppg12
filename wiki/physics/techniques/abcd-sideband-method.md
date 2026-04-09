# The ABCD Sideband Method for Isolated Photon Background Subtraction

## Origin and Principle

The two-dimensional sideband method (commonly called "ABCD") is a data-driven technique for estimating residual background in a signal region defined by two weakly correlated discriminating variables. It was first applied to isolated photon measurements by ATLAS in the 7 TeV analysis (1012.4389, Phys. Rev. D 83, 052005, 2011) and has since become the standard approach at the LHC and now at RHIC/sPHENIX.

### The 2x2 matrix

Two discriminating variables divide the event sample into four regions:

```
                Isolated (pass iso cut)     Non-isolated (fail iso cut)
Tight ID        A = signal region           B = background control
Non-tight ID    C = background control      D = background control
```

- **Axis 1 (horizontal):** Isolation energy -- separates signal (low hadronic activity around photon) from jets (embedded in hadronic debris)
- **Axis 2 (vertical):** Photon identification -- separates photon-like showers (narrow, symmetric) from hadron-like showers (broad, irregular)

### Core assumption

The method assumes that for background events, the two discriminating variables are uncorrelated:

```
P_bkg(tight AND iso) / P_bkg(tight AND noniso) = P_bkg(nontight AND iso) / P_bkg(nontight AND noniso)
```

Equivalently, the background populates the 2D plane factorizably: the probability of passing the tight cut is independent of whether the background candidate is isolated or not. Under this assumption:

```
N_A^bkg / N_B^bkg = N_C^bkg / N_D^bkg
```

and the background in the signal region A can be estimated from the three control regions:

```
N_A^bkg = N_B^bkg * N_C^bkg / N_D^bkg
```

The signal yield is then:

```
N_A^sig = N_A - N_B * N_C / N_D     [naive, ignoring signal leakage]
```

## Signal Leakage Corrections

Real signal photons are not perfectly confined to region A. A small fraction leaks into the control regions due to finite detector resolution:

- **c_B = N_B^sig / N_A^sig:** Signal photons that are tight but fail the isolation cut (tail of isolation distribution)
- **c_C = N_C^sig / N_A^sig:** Signal photons that are isolated but fail the tight ID cut (BDT score or shower shape below threshold)
- **c_D = N_D^sig / N_A^sig:** Signal photons that fail both tight and isolated cuts

These leakage fractions are extracted from signal MC simulation. Typical values from ATLAS (1012.4389):

| Fraction | Typical range | Trend with E_T |
|----------|--------------|----------------|
| c_B | 3-17% | Increases (isolation tail grows) |
| c_C | 2-14% | Decreases (ID efficiency improves) |
| c_D | < 2% | Always small (both cuts fail rarely) |

### The self-consistent equation

Accounting for signal leakage, the background estimate uses corrected control region counts:

```
N_A^sig = N_A - R_bg * (N_B - c_B * N_A^sig)(N_C - c_C * N_A^sig) / (N_D - c_D * N_A^sig)
```

This is implicit in N_A^sig -- the signal yield appears on both sides. It must be solved iteratively or by root-finding. Rearranging:

```
f(N_A^sig) = N_A^sig - N_A + R_bg * (N_B - c_B * N_A^sig)(N_C - c_C * N_A^sig) / (N_D - c_D * N_A^sig) = 0
```

This is a rational equation in N_A^sig with a unique physical root (0 < N_A^sig < N_A). PPG12 solves it numerically via `myfunc()` in `CalculatePhotonYield.C`.

### The R_bg correlation factor

If the background isolation and identification variables are not perfectly uncorrelated, the "double ratio" R_bg corrects the factorization assumption:

```
R_bg = (N_A^bkg * N_D^bkg) / (N_B^bkg * N_C^bkg)
```

R_bg = 1 means perfect decorrelation. Deviations from unity indicate that background candidates passing tight ID are systematically more or less isolated than those failing tight ID.

**ATLAS validation (1908.02746):** The control regions B and D are each subdivided into B'/B'' and D'/D'' using additional E_T^iso boundaries. R_bg is measured in these sub-regions after subtracting signal leakage:

```
R_bg = (N_B'^bkg / N_D'^bkg) / (N_B''^bkg / N_D''^bkg)
```

Deviations of 15-40% from unity are observed depending on ET and eta, and propagated as a systematic uncertainty. For |eta| < 1.37, the effect on the cross section is typically < 0.6%; for 1.52 < |eta| < 1.81, it reaches ~3.6%.

## Purity

The signal purity in region A is:

```
P = N_A^sig / N_A = 1 - N_A^bkg / N_A
```

Equivalently, after solving the ABCD equation:

```
P = N_A^sig / N_A
```

Purity increases with ET (harder photons are better separated from jets) and depends on the tightness of both the ID and isolation cuts:

| Experiment | pT range | Purity range |
|------------|----------|-------------|
| ATLAS (13 TeV) | 125-2500 GeV | > 90% in all bins |
| ALICE (13 TeV) | 7-200 GeV/c | 5% (7 GeV) to 80% (>80 GeV) |
| ALICE (7 TeV) | 10-60 GeV/c | 20% (10 GeV) to 60% (>18 GeV) |
| PPG12 (200 GeV) | 8-35 GeV | Low at 8 GeV, increasing with ET |

The low purity at the bottom of the pT range is the primary experimental challenge for PPG12 and ALICE, where the ABCD correction subtracts a large fraction of region A.

## Evolution of the ABCD Method

### ATLAS 7 TeV (2011): The Pioneer Implementation

**Paper:** 1012.4389

The foundational ATLAS measurement established the template:

- **Axis 1:** Calorimetric isolation E_T^iso (cone R = 0.4, UE-corrected via jet-area method)
  - Isolated: E_T^iso < 3 GeV (flat threshold)
  - Non-isolated: E_T^iso > 5 GeV
- **Axis 2:** Tight shower-shape identification (9 variables from EM calorimeter first and second layers)
  - Non-tight: fails at least one of five first-layer variables (w_{s,tot}, E_ratio, Delta_E, F_side, w_{s,3}), chosen specifically to minimize correlation with isolation

The non-tight definition was carefully designed: only first-layer shower shape variables are inverted, because these are least correlated with calorimetric isolation (which sums over the full cone including hadronic calorimeter). Second-layer variables like R_eta and w_2 were not inverted because they show stronger correlation with isolation for background.

A cross-check was provided by an independent isolation template fit method: the E_T^iso distribution of tight photon candidates is fitted to signal + background templates (signal from Z -> ee electrons, background from non-tight candidates). Results agreed with ABCD within 2%.

### ATLAS 13 TeV (2019): Mature Implementation

**Paper:** 1908.02746

The method was refined with 36 fb^-1:

- **ET-dependent isolation:** E_T^iso < 4.8 + 0.0042*ET GeV (maintained since 8 TeV)
- **Gap between iso and noniso:** E_T^iso,cut + 2 GeV as lower bound for non-isolated region
- **Non-tight:** Four first-layer variables inverted (w_{s,3}, f_side, Delta_E_s, E_ratio)
- **Signal leakage:** From PYTHIA MC; cross-checked with SHERPA
- **R_bg validation:** Systematic subdivision of control regions into B'/B'', D'/D'' with progressive E_T^iso boundary shifts; maximum 15-40% deviations propagated as systematic
- **Signal purity:** > 90% everywhere (benefit of high ET threshold at 125 GeV)
- **NNLO comparison:** First comparison to NNLOJET NNLO QCD, validating the extracted signal yield at the highest precision

### ALICE 7 TeV (2019): Sampling Calorimeter Adaptation

**Paper:** 1906.01371

ALICE adapted the ABCD method for a lead-scintillator sampling EMCal -- a detector more similar to sPHENIX than the crystal calorimeters of ATLAS/CMS:

- **Axis 1:** Combined isolation (EMCal clusters + charged tracks, cone R = 0.4)
  - Isolated: p_T^iso < 2 GeV/c
  - Non-isolated: p_T^iso > 3 GeV/c
- **Axis 2:** Shower shape eigenvalue sigma_long^2 (larger eigenvalue of energy distribution in eta-phi plane)
  - Narrow (photon-like): 0.1 < sigma_long^2 < sigma_max (pT-dependent: 0.4 at 10 GeV, 0.3 above 16 GeV)
  - Wide (background): pT-dependent ranges, e.g. 0.55 < sigma_long^2 < 1.75

A key innovation was the MC correction factor alpha_MC:

```
P = 1 - (N_narrow^iso / N_narrow^{anti-iso}) / (N_wide^iso / N_wide^{anti-iso}) * alpha_MC
```

This correction accounts for the correlation between sigma_long^2 and isolation momentum for background clusters. ALICE found that the double ratio of background isolation fractions in narrow vs wide clusters differs from unity, requiring the MC correction. A linear extrapolation of the double ratio f(sigma_long^2) was used to evaluate the correction in the signal region from the background control region.

### ALICE 13 TeV (2024): Tracks-Only Isolation

**Paper:** 2407.01165

The most recent ALICE measurement refined the method:

- **Axis 1:** Charged-track-only isolation (p_T^iso,ch < 1.5 GeV/c, anti-iso: 2.5-10 GeV/c)
- **Axis 2:** sigma_long^2 with simplified flat threshold (0.1 < sigma_long^2 < 0.3)
- **Purity smoothing:** Fitted measured purity with double-sigmoid functions to reduce statistical fluctuations:
  ```
  P(pT) = a / (1 + exp(-b * (pT - c)))
  ```
  with separate parameters for low-pT (7-30 GeV, a=0.617, b=0.292, c=14.28) and high-pT (30-200 GeV, a=0.852, b=0.034, c=2.4) ranges

## PPG12's BDT-Based Variant

PPG12 replaces the traditional shower-shape variable with a multivariate BDT score, while keeping the same 2x2 ABCD structure.

### Axis definitions

```
              Isolated                        Non-isolated
              (isoET < b + s*ET)              (isoET > b + s*ET + gap)

Tight         A = signal                      B = bkg template
(BDT > intercept + slope*ET)

Non-tight     C = bkg control                 D = bkg normalization
(BDT < intercept + slope*ET)
```

- **Tight/non-tight boundary:** Parametric BDT threshold
  ```
  bdt_min = tight_bdt_min_intercept + tight_bdt_min_slope * ET
  ```
  Nominal: intercept = 0.80, slope = -0.015 (decreasing with ET because the BDT becomes more discriminating at higher ET)

- **Isolated/non-isolated boundary:** Parametric isolation threshold
  ```
  iso_max = reco_iso_max_b + reco_iso_max_s * ET
  ```
  Nominal: b = 0.502 GeV, s = 0.0433

### Comparison: ATLAS shower shape vs PPG12 BDT

| Aspect | ATLAS shower shape | PPG12 BDT |
|--------|-------------------|-----------|
| Discriminant | Individual shower shape variables with rectangular cuts | XGBoost classifier combining 7-11 shower shape features |
| Non-tight definition | Fail specific first-layer variables | BDT score below parametric threshold |
| ET dependence | Fixed cuts (no ET dependence in individual variables) | ET-dependent threshold (intercept + slope*ET) |
| Correlation with isolation | Minimized by choice of inverted variables | Minimized by training without isolation features |
| Interpretability | Transparent (each variable has physical meaning) | Opaque (multivariate combination) |
| Discrimination power | Limited by rectangular cuts in multi-dimensional space | Exploits correlations between features |

The BDT approach is more powerful (higher signal efficiency at fixed background rejection) but requires more careful validation of the decorrelation assumption: the BDT score can develop subtle correlations with isolation that are not present in individual shower-shape variables.

### Signal leakage in PPG12

PPG12 computes leakage fractions from merged signal MC (photon5 + photon10 + photon20 weighted by cross section):

```
cB = h_tight_noniso_cluster_signal / h_tight_iso_cluster_signal
cC = h_nontight_iso_cluster_signal / h_tight_iso_cluster_signal
cD = h_nontight_noniso_cluster_signal / h_tight_iso_cluster_signal
```

These are ET-dependent (computed per pT bin) and used in the self-consistent ABCD equation.

### R factor in PPG12

PPG12 computes R_bg two ways:

1. **Data-driven:** Background yields estimated as (data - MC signal leakage) per region:
   ```
   R = (A_bkg * D_bkg) / (B_bkg * C_bkg)
   ```

2. **MC truth:** Using truth-unmatched histograms:
   ```
   R_notmatch = (A_notmatch * D_notmatch) / (B_notmatch * C_notmatch)
   ```

Both are monitored for consistency. PPG12 currently uses R = 1 (uncorrelated assumption) in the default analysis, with deviations propagated as a systematic uncertainty.

### Purity extraction: Toy MC approach

PPG12 uses a Poisson toy MC to propagate statistical uncertainties on the purity:

1. Compute effective events per bin: N_eff = N^2 / sum(w^2) (accounts for weighted MC)
2. Draw 20,000 Poisson realizations of each region count
3. Solve the ABCD equation for each realization
4. Fit the resulting purity distribution with a Gaussian to extract central value and uncertainty
5. Repeat with Gaussian-smeared leakage fractions (c_B, c_C, c_D) to include MC statistical uncertainty on the leakage

The bin-by-bin purity values are then smoothed with a functional fit:
- Error function: P(ET) = p0 * Erf((ET - p1) / p2)
- Or Pade rational: P(ET) = (p0 + p1*ET) / (1 + p2*ET)

Confidence intervals at 68.3% CL are computed via `TVirtualFitter::GetConfidenceIntervals()`. The `fittingerror` config parameter shifts the purity by +/- 1 confidence interval for systematic studies.

### MC purity correction

When `mc_purity_correction = 1`, PPG12 runs an additional closure test:
1. Apply the ABCD method to MC as if it were data
2. Compare the extracted purity to the truth purity (from truth matching)
3. Compute a correction factor: mc_corr = truth_purity / extracted_purity
4. Apply: purity_corrected = extracted_purity * mc_corr

This corrects for violations of the ABCD independence assumption that are captured by the simulation.

## Systematic Uncertainty Sources

The ABCD method introduces several systematic uncertainties, common across experiments:

| Source | ATLAS treatment | PPG12 treatment |
|--------|----------------|-----------------|
| Non-tight definition | Vary which variables are inverted | Vary BDT threshold (intercept +/- delta) |
| Isolation boundary | Vary E_T^iso cut by +/- 1 GeV | Vary intercept b and slope s |
| Iso/noniso gap | Vary lower noniso boundary | Vary reco_noniso_min_shift |
| Signal leakage (c_B, c_C, c_D) | Pythia vs Sherpa MC comparison | Gaussian smearing in toy MC |
| R_bg correlation | Subdivide B, D into sub-regions | Data-driven vs MC truth comparison |
| Purity fitting | Not applicable (per-bin) | Vary fit function (erf vs Pade) |
| MC purity correction | Not applicable | On/off toggle |

PPG12 generates 23+ systematic variation configs via `make_bdt_variations.py`, each varying one or more ABCD parameters. These are grouped into "purity" and "efficiency" systematic categories and aggregated in quadrature by `calc_syst_bdt.py`.

## Mathematical Summary

### Full ABCD equation with all corrections

```
N_A^sig = N_A - R_bg * [(N_B - c_B * N_A^sig) * (N_C - c_C * N_A^sig)] / (N_D - c_D * N_A^sig)
```

### Purity

```
P = N_A^sig / N_A
```

### Cross section (PPG12)

```
d sigma / d ET = (1 / L) * (1 / Delta_ET) * (N_A * P) / epsilon * U
```

where L is integrated luminosity, epsilon is the total efficiency (trigger * reconstruction * ID * isolation), and U is the unfolding correction (Bayesian, not bin-by-bin).

### Equivalent formulas used across experiments

| Experiment | Cross section formula |
|------------|---------------------|
| ATLAS | d sigma / d ET = N_yield * U / (L * Delta_ET * epsilon_trig * epsilon_reco * epsilon_ID) |
| ALICE | d^2 sigma / (dpT deta) = (1/epsilon_trig) * (1/L) * P / (epsilon_gamma^iso * kappa_iso) * dN / dpT |
| CMS | d^2 sigma / (dET deta) = U(N_gamma) / (Delta_ET * Delta_eta * epsilon * SF * L) |
| PPG12 | d sigma / d ET = N_yield * P * U / (L * Delta_ET * epsilon) |

## Key Code Locations in PPG12

| File | Role |
|------|------|
| `efficiencytool/CalculatePhotonYield.C` | ABCD equation solver (`myfunc()`), leakage computation, purity fitting, toy MC |
| `efficiencytool/RecoEffCalculator_TTreeReader.C` | Fills ABCD region histograms from data and MC |
| `efficiencytool/ShowerShapeCheck.C` | Shower shape distributions per ABCD region |
| `FunWithxgboost/make_bdt_variations.py` | Generates systematic variation configs for ABCD parameters |
| `plotting/calc_syst_bdt.py` | Aggregates systematic variations into purity/efficiency groups |

## Key References

| Reference | Contribution |
|-----------|-------------|
| ATLAS, 1012.4389 (2011) | First ABCD implementation for isolated photons; defined non-tight criteria, leakage corrections, template fit cross-check |
| ATLAS, 1908.02746 (2019) | Mature implementation with ET-dependent isolation, R_bg validation via control region subdivision, NNLO comparison |
| ALICE, 1906.01371 (2019) | Adapted ABCD for sampling calorimeter (sigma_long^2); introduced alpha_MC correction for background correlation |
| ALICE, 2407.01165 (2024) | Tracks-only isolation, sigmoid purity smoothing, lowest-xT isolated photon measurement |
| CMS, 1108.2044 (2011) | Template fit alternative (1D), sigma_eta_eta shower shape, PF isolation |

## See Also

- [ABCD Method](../../concepts/abcd-method.md) -- PPG12-specific implementation: histogram naming, config fields, code walkthrough
- [Isolation Cuts](../../concepts/isolation-cuts.md) -- Parametric isolation details, branch names, MC corrections
- [Isolation Theory](../theory/isolation-theory.md) -- Frixione vs standard cone, fragmentation suppression, IR safety
- [Isolation Methods](isolation-methods.md) -- Experimental isolation implementations compared across experiments
