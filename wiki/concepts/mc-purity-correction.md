# MC Purity Correction

## Overview

The MC purity correction factor compensates for bias in the [ABCD method](abcd-method.md) purity estimate. It is defined as the ratio of truth purity (from generator-level matching) to the Pade fit of the leakage-corrected ABCD purity:

```
correction = truth_purity / f_purity_leak_fit(pT)
```

The correction is stored as `g_mc_purity_fit_ratio` and applied to data purity when `mc_purity_correction: 1` is set in the config.

A notable feature is that the correction crosses 1.0 at pT ~ 13--14 GeV:
- **Below the crossing:** ABCD overestimates background, underestimates purity, so correction > 1
- **Above the crossing:** ABCD underestimates background, overestimates purity, so correction < 1

**Technical detail:** the stored ratio uses truth purity *data points* divided by a *smooth Pade fit* to the ABCD purity. The correction therefore inherits truth purity statistical fluctuations. When applied to data via `TGraphErrors::Eval()`, linear interpolation between these noisy points is used.

## Correction Factor Table

> **Note:** The table below shows the current ET-dependent non-tight configuration (2026-04-08). For the previous fixed non-tight (0.50) values, see `reports/mc_purity_correction_etdep_nontight.md` Table 4.

| pT (GeV) | Truth Purity | ABCD Purity (leak) | Ratio |
|-----------|-------------|---------------------|-------|
| 9 | 0.576 | 0.351 | 1.702 |
| 11 | 0.609 | 0.453 | 1.258 |
| 13 | 0.632 | 0.603 | 1.080 |
| 15 | 0.637 | 0.659 | 0.966 |
| 17 | 0.647 | 0.714 | 0.904 |
| 19 | 0.714 | 0.741 | 0.938 |
| 21 | 0.739 | 0.820 | 0.927 |
| 23 | 0.713 | 0.857 | 0.862 |
| 25 | 0.732 | 0.753 | 0.859 |
| 27 | 0.756 | 0.869 | 0.865 |
| 30 | 0.761 | 0.894 | 0.845 |
| 34 | 0.777 | 0.741 | 0.837 |

The ABCD purity points are non-monotonic at high pT, reflecting ABCD limitations and poor MC statistics. The truth purity rises from ~0.58 to ~0.78 across 9--34 GeV.

## Root Cause: ABCD Independence Violation

The ABCD method assumes that the BDT tight/non-tight classification and the isolation criterion are uncorrelated for background. This assumption is violated.

The true background R-factor:

```
R_true = (A_bkg * D_bkg) / (B_bkg * C_bkg)
```

should equal 1.0 if the two axes are independent. In practice, R_true varies strongly with pT:

> **Note:** Table updated to current ET-dependent non-tight config (2026-04-08). R_true = h_R from `Photon_final_bdt_nom_mc.root`, computed as (A_bkg × D_bkg)/(B_bkg × C_bkg) after truth signal subtraction from jet MC.

| pT bin (GeV) | R_true | Deviation from 1.0 |
|--------------|--------|-------------------|
| [8, 10] | 0.675 | -33% |
| [10, 12] | 0.729 | -27% |
| [12, 14] | 0.904 | -10% |
| [14, 16] | 1.058 | +6% |
| [16, 18] | 1.261 | +26% |
| [18, 20] | 0.991 | -1% |
| [20, 22] | 1.273 | +27% |
| [22, 24] | 1.779 | +78% |
| [24, 26] | 1.151 | +15% |
| [26, 28] | 1.617 | +62% |
| [28, 32] | 1.854 | +85% |
| [32, 36] | 0.827 | -17% |

R_true crosses 1.0 between the [12,14] and [14,16] GeV bins. The mechanism:

- **Below pT = 13 GeV:** R_true < 1 means ABCD overestimates background, so purity is too low, correction > 1
- **Above pT = 13 GeV:** R_true > 1 means ABCD underestimates background, so purity is too high, correction < 1

Fluctuations at high pT (e.g., R = 0.76 at [24, 26] GeV followed by R = 2.09 at [26, 28] GeV) reflect limited MC statistics.

**Corrected interpretation (2026-04-08):** The R crossover is a genuine background physics effect. After removing signal contamination from the jet MC (see "Signal Contamination in Jet MC" section), the pure background R (h_R) ranges 0.67-1.85, crossing 1.0 at pT ~ 14 GeV. At low pT, tight and non-tight background have nearly identical isolation properties; at high pT, photon-like background jets become marginally more isolated. Earlier diagnostics using the signal-contaminated isolation fraction ratio (0.77-1.77) overstated the intrinsic background correlation because the jet MC tight region contains 58-78% truth-matched prompt photons.

## ET-Dependent BDT Threshold and Non-Tight Boundary

The tight BDT threshold is ET-dependent: `bdt_min = 0.80 - 0.015 * ET`.

### Current configuration (ET-dependent non-tight upper bound)

As of 2026-04-08, the non-tight upper bound also uses `0.80 - 0.015 * ET`, matching the tight lower bound exactly. The BDT axis is a clean partition: non-tight = [0.02, threshold), tight = [threshold, 1.0], with no dead zone at any ET.

### Previous configuration (fixed non-tight at 0.50)

Previously the non-tight BDT window was fixed at [0.02, 0.50], creating a "dead zone" between the two regions:

| ET (GeV) | Tight threshold | Dead zone (old) | Width (old) | Dead zone (new) |
|----------|----------------|-----------------|-------------|-----------------|
| 8 | 0.68 | [0.50, 0.68] | 0.18 | 0 |
| 12 | 0.62 | [0.50, 0.62] | 0.12 | 0 |
| 20 | 0.50 | vanishes | 0 | 0 |
| 25 | 0.425 | overlap | --- | 0 |

Eliminating the dead zone was expected to improve ABCD independence. Instead, it **increased the BDT-isolation correlation** (see isolation fraction ratio below). The dead zone had accidentally excluded intermediate-BDT jets whose isolation properties fell between tight and non-tight populations, masking the intrinsic correlation at low pT.

This is **not** caused by the BDT model boundary at ET = 15 GeV: both ET bins use the same model (`base_v1E`).

## Pade Fit Quality Issues

The ABCD purity is fitted with a Pade function `f(x) = (p0 + p1*x) / (1 + p2*x)` using 11 of 12 pT bins (the last bin at pT = 34 GeV lies outside the fit range [8, 32] GeV):

- **Leak fit:** chi2/NDF = 18.9/8 = 2.36, with a pull of 4.67 at pT = 13 GeV (poor fit quality)
- **Truth fit:** chi2/NDF = 6.1/8 = 0.76 (good fit quality)
- The Pade form forces a monotonic rise, but the ABCD points are non-monotonic (e.g., 0.56 at pT = 25, 0.92 at pT = 27 GeV)

Because the stored correction uses truth *points* divided by the smooth leak *fit*, the fit quality primarily affects the denominator. The poor leak fit means the denominator may be systematically biased at individual pT bins.

## Statistical Limitations at High pT

The effective number of events in region D (the denominator of `B * C / D`) drops rapidly at high pT:

| pT (GeV) | D_eff |
|-----------|-------|
| < 14 | > 1000 |
| 23 | 107 |
| 25 | 71 |
| 27 | 19 |
| 34 | 7 |

Above pT = 25 GeV, the ABCD estimate is statistically unreliable. This explains the wild fluctuations in the ABCD purity at high pT and correspondingly large uncertainties on the correction factor.

## Leakage Fractions

Signal leakage fractions (from signal MC, current ET-dependent non-tight config):

- **cB** = 5.7--8.8% (signal leaking to tight + non-isolated)
- **cC** = 2.2--6.7% (signal leaking to non-tight + isolated; increased at low pT due to wider non-tight window)
- **cD** < 0.8% (signal leaking to non-tight + non-isolated)

The leakage correction moves the ABCD purity in the right direction but cannot compensate for the fundamental independence violation.

## Isolation Fraction Ratio: Tight vs Non-Tight Background

The isolation fraction f_iso = N_isolated / N_total for background jets in each BDT region. ABCD independence requires f_iso(tight) / f_iso(nontight) = 1.0.

Results from jet-only MC (current ET-dependent non-tight config, 2026-04-08):

> **Important caveat:** The values below are from the **total** jet MC, which includes signal contamination (58-78% of region A counts are truth-matched prompt photons; see "Signal Contamination in Jet MC" below). The R_check column (AD/BC from jet MC totals) is the signal-contaminated R_raw, NOT the pure background R. The pure background R (= h_R, after truth subtraction) ranges 0.67-1.85 and crosses 1.0 at pT ~ 14 GeV.

| pT (GeV) | f_iso(tight) | f_iso(nontight) | Ratio | R_check (AD/BC, signal-contaminated) |
|-----------|-------------|-----------------|-------|--------------------------------------|
| [8, 10] | 0.526 | 0.431 | 1.22 | 1.47 |
| [10, 12] | 0.539 | 0.403 | 1.34 | 1.73 |
| [12, 14] | 0.560 | 0.357 | 1.57 | 2.29 |
| [14, 16] | 0.542 | 0.311 | 1.74 | 2.62 |
| [16, 18] | 0.539 | 0.281 | 1.92 | 3.00 |
| [20, 22] | 0.545 | 0.222 | 2.46 | 4.21 |
| [26, 28] | 0.603 | 0.240 | 2.52 | 4.82 |
| [28, 32] | 0.625 | 0.213 | 2.94 | 6.16 |

Previous (fixed non-tight at 0.50): ratio ranged 0.77--1.77, straddling 1.0.
Current (ET-dependent non-tight): ratio ranges **1.22--2.94, all above 1.0**.

The tight-BDT background is systematically more isolated (mean isoET 2.3--2.9 GeV) than non-tight (2.8--6.0 GeV). See `plotting/figures/isoET_tight_vs_nontight_jet_bdt_nom.pdf` for per-pT-bin isoET overlay plots. For the corrected background-only comparison (after truth signal subtraction), see `plotting/figures/bkg_isoET_tight_vs_nontight_jet_bdt_nom.pdf`.

## Signal Contamination in Jet MC

The jet MC samples (jet10, jet15, jet20, jet30, jet50) contain truth-matched prompt photons that pass all analysis cuts. These signal photons dominate region A:

| pT bin (GeV) | Signal fraction in A |
|:---:|:---:|
| [8, 10] | 57.6% |
| [10, 12] | 60.9% |
| [12, 14] | 63.2% |
| [14, 16] | 63.7% |
| [16, 18] | 64.7% |
| [18, 20] | 71.4% |
| [20, 22] | 73.9% |
| [22, 24] | 71.3% |
| [24, 26] | 73.2% |
| [26, 28] | 75.6% |
| [28, 32] | 76.1% |
| [32, 36] | 77.7% |

Signal fractions in other regions are lower (B: 4-11%, C: 2-21%, D: <1%).

This contamination means that **any diagnostic computed from raw jet MC totals** (isolation fractions, R from AD/BC, isoET overlays) is dominated by signal in the tight-isolated region. Two distinct R factors exist:

- **R_raw** (signal-contaminated, from jet MC totals): 1.47-6.16, always > 1, monotonically increasing with pT
- **R_bkg** (= h_R, pure background after truth subtraction): 0.675-1.854, crosses 1.0 at pT ~ 14 GeV

The earlier conclusion that the "pure background R" was always > 1 and that the sign flip in h_R was caused by signal leakage corrections was incorrect -- those R values were computed from signal-contaminated jet MC totals.

The pure background isoET distributions (after truth subtraction) show that tight and non-tight background have nearly identical isolation at low pT, with only a subtle divergence at high pT. See `plotting/figures/R_decomposition_bdt_nom.pdf` for the full R factor decomposition.

## Recommendations

1. The MC purity correction addresses a **real bias** (ABCD independence violation), but the correction factor has large statistical uncertainties from truth purity, especially at high pT.
2. The sign change at pT ~ 14 GeV is a **genuine background physics effect** (R_bkg crossing 1.0), not a fit artifact or signal contamination artifact. It persists regardless of the non-tight boundary definition (confirmed with both fixed and ET-dependent configs).
3. The intrinsic BDT-isolation correlation in pure background is subtle: tight and non-tight background have nearly identical isolation at low pT, with a slight divergence at high pT. The large apparent correlation seen in raw jet MC diagnostics is dominated by signal contamination (58-78% of region A).
4. **Jet MC diagnostics must account for signal contamination**: always use truth-subtracted quantities (h_R, background-only isoET) for physics interpretation of background behavior.
5. **Consider:** using `h_R` from MC to correct the ABCD R-factor directly, rather than correcting the final purity post-hoc.
6. High-pT bins (> 25 GeV) have very poor MC statistics -- the correction factor there is unreliable.
7. ~30 systematic variation configs still use the old fixed non-tight bound (missing `bdt_max_slope`/`bdt_max_intercept`). Regenerate if ET-dependent bound is intended for all variations.

## Key Files

| File | Purpose |
|------|---------|
| `reports/mc_purity_correction_investigation.tex` | Full investigation report (original, fixed non-tight) |
| `reports/mc_purity_correction_etdep_nontight.md` | Re-evaluation with ET-dependent non-tight (2026-04-08) |
| `efficiencytool/CalculatePhotonYield.C` | ABCD equation solver, purity fitting, MC correction (`g_mc_purity_fit_ratio`) |
| `efficiencytool/RecoEffCalculator_TTreeReader.C` | Tight/non-tight classification, ABCD region filling |
| `plotting/figures/isoET_tight_vs_nontight_jet_bdt_nom.pdf` | isoET shape: tight vs non-tight background (per pT bin) |
| `plotting/figures/isoET_tight_jetonly_vs_inclusive_bdt_nom.pdf` | Tight isoET: jet-only vs inclusive MC |
| `plotting/figures/purity_fit_diagnostics.pdf` | Pade fit curves, residuals, ratio crossing analysis |
| `plotting/figures/leakage_fractions_vs_pT.pdf` | Signal leakage fractions vs pT |
| `plotting/figures/purity_comparison.pdf` | Truth purity vs leakage-corrected ABCD purity |
| `plotting/figures/R_decomposition_bdt_nom.pdf` | R factor decomposition: R_raw vs R_bkg vs signal fraction |
| `plotting/figures/bkg_isoET_tight_vs_nontight_jet_bdt_nom.pdf` | Pure background isoET overlay (truth-subtracted) |
| `plotting/plot_R_decomposition.C` | R decomposition plotting macro |
| `plotting/plot_bkg_isoET_tight_vs_nontight.C` | Background-only isoET plotting macro |

## See Also

- [ABCD Method](abcd-method.md) -- Background subtraction, leakage corrections, R-factor, purity fitting
- [Systematic Variations](systematic-variations.md) -- Purity-related systematic variations
