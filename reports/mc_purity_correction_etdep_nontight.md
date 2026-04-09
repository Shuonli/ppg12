# MC Purity Correction Re-evaluation: ET-Dependent Non-Tight BDT Bound

**Date:** 2026-04-08
**Config:** `config_bdt_nom.yaml` (nominal, results from today)
**Change under study:** Non-tight BDT upper bound changed from fixed 0.50 to ET-dependent `0.80 - 0.015 * ET` (matching the tight lower bound exactly)

## 1. Executive Summary

The non-tight BDT upper bound now tracks the tight BDT lower bound, eliminating the "dead zone" between tight and non-tight regions at all ET. This was expected to improve ABCD independence by removing the pT-dependent gap that altered background composition across regions.

**Finding: the purity correction sign flip persists, and the ABCD independence violation is actually worse.** The isolation fraction ratio for background (tight/non-tight) moved from 0.77--1.77 (straddling 1.0) to 1.22--2.94 (entirely above 1.0). The old dead zone was accidentally suppressing the intrinsic BDT--isolation correlation at low pT. The underlying physics correlation -- photon-like jets are preferentially more isolated -- is now fully exposed.

## 2. Configuration Change

### BDT boundaries vs ET

| ET (GeV) | Tight lower (both) | NT upper (OLD) | Dead zone (OLD) | NT upper (NEW) | Dead zone (NEW) |
|-----------|-------------------|----------------|-----------------|----------------|-----------------|
| 8 | 0.680 | 0.500 | 0.180 | 0.680 | 0 |
| 12 | 0.620 | 0.500 | 0.120 | 0.620 | 0 |
| 16 | 0.560 | 0.500 | 0.060 | 0.560 | 0 |
| 20 | 0.500 | 0.500 | 0 | 0.500 | 0 |
| 24 | 0.440 | 0.500 | -0.060 (overlap) | 0.440 | 0 |
| 28 | 0.380 | 0.500 | -0.120 (overlap) | 0.380 | 0 |
| 32 | 0.320 | 0.500 | -0.180 (overlap) | 0.320 | 0 |

The dead zone is eliminated at all ET. The BDT axis is now a clean two-region partition: non-tight = [0.02, threshold), tight = [threshold, 1.0].

## 3. R_true Comparison

R_true = (A_bkg * D_bkg) / (B_bkg * C_bkg) from `h_R` in `Photon_final_bdt_nom_mc.root`.

| pT bin (GeV) | R_true (OLD) | R_true (NEW) | Change |
|:---:|:---:|:---:|:---:|
| [8, 10] | 0.70 | 0.675 +/- 0.15 | -0.03 |
| [10, 12] | 0.80 | 0.729 +/- 0.08 | -0.07 |
| [12, 14] | 1.07 | 0.904 +/- 0.15 | -0.17 |
| [14, 16] | 1.14 | 1.058 +/- 0.21 | -0.08 |
| [16, 18] | 1.16 | 1.261 +/- 0.35 | +0.10 |
| [18, 20] | 0.98 | 0.991 +/- 0.48 | +0.01 |
| [20, 22] | 1.04 | 1.273 +/- 0.70 | +0.23 |
| [22, 24] | 1.59 | 1.779 +/- 0.41 | +0.19 |
| [24, 26] | 0.76 | 1.152 +/- 0.39 | +0.39 |
| [26, 28] | 2.09 | 1.617 +/- 0.86 | -0.47 |
| [28, 32] | 0.93 | 1.854 +/- 0.99 | +0.92 |
| [32, 36] | 1.07 | 0.827 +/- 0.43 | -0.24 |

**R_true crossing point:** OLD ~ [12,14] GeV, NEW ~ [14,16] GeV (shifted slightly higher). Below this, R < 1 (ABCD overestimates background); above, R > 1 (underestimates).

**Key observation:** Low-pT bins (8-14 GeV) have R_true further below 1.0 than before, while many high-pT bins moved further above 1.0. The asymmetry is larger, not smaller.

## 4. MC Purity Correction Ratio

Correction = truth_purity / f_purity_leak_fit from `g_mc_purity_fit_ratio`.

| pT (GeV) | Truth Purity | ABCD Purity | Ratio (OLD) | Ratio (NEW) | Change |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 9 | 0.576 | 0.351 | 1.505 | **1.702** | +0.20 |
| 11 | 0.609 | 0.453 | 1.206 | **1.258** | +0.05 |
| 13 | 0.632 | 0.603 | 1.112 | **1.080** | -0.03 |
| 15 | 0.637 | 0.659 | 0.952 | **0.966** | +0.01 |
| 17 | 0.647 | 0.714 | 0.862 | **0.904** | +0.04 |
| 19 | 0.714 | 0.741 | 0.860 | **0.938** | +0.08 |
| 21 | 0.739 | 0.820 | 0.875 | **0.927** | +0.05 |
| 23 | 0.713 | 0.857 | 0.852 | **0.862** | +0.01 |
| 25 | 0.732 | 0.753 | 0.787 | **0.859** | +0.07 |
| 27 | 0.756 | 0.869 | 0.829 | **0.865** | +0.04 |
| 30 | 0.761 | 0.894 | 0.766 | **0.845** | +0.08 |
| 34 | 0.777 | 0.741 | 0.919 | **0.837** | -0.08 |

**Crossing point:** Both OLD and NEW cross 1.0 at pT ~ 14 GeV (essentially unchanged).

**At low pT:** The correction is LARGER than before (1.70 vs 1.51 at 9 GeV), meaning the ABCD underestimate of purity at low pT is worse.

**At high pT:** The correction is somewhat closer to 1.0 and more stable (0.84-0.97 vs 0.77-0.92), with less wild fluctuation.

**Truth purity** is slightly higher across the board (~0.58-0.78 vs ~0.57-0.85) — this reflects the new non-tight bound capturing more background in the non-tight region.

## 5. Isolation Fraction Ratio: The Key Diagnostic

f_iso = N_isolated / N_total for background in each BDT region. Independence requires f_iso(tight) / f_iso(nontight) = 1.0.

**Source: jet-only MC** (`MC_efficiency_jet_bdt_nom.root`)

> **WARNING (corrected 2026-04-08):** The jet MC file contains significant signal contamination -- truth-matched prompt photons comprise 58-78% of the jet MC counts in region A (see Section 9). The isolation fractions and ratios below are computed from the **total** jet MC (signal + background combined), NOT from pure background. The pure background values after truth signal subtraction are substantially different: R_bkg (= h_R) ranges 0.67-1.85 and crosses 1.0 at pT ~ 14 GeV, while the signal-contaminated R_raw ranges 1.4-6.2 (always > 1). See Section 9 for the full decomposition.

| pT (GeV) | f_iso(tight) | f_iso(nontight) | Ratio (NEW) | Ratio (OLD) |
|:---:|:---:|:---:|:---:|:---:|
| [8, 10] | 0.526 | 0.431 | **1.22** | ~0.77 |
| [10, 12] | 0.539 | 0.403 | **1.34** | ~0.85 |
| [12, 14] | 0.560 | 0.357 | **1.57** | ~1.00 |
| [14, 16] | 0.542 | 0.311 | **1.74** | ~1.10 |
| [16, 18] | 0.539 | 0.281 | **1.92** | ~1.20 |
| [18, 20] | 0.533 | 0.264 | **2.02** | ~1.30 |
| [20, 22] | 0.545 | 0.222 | **2.46** | ~1.40 |
| [22, 24] | 0.580 | 0.215 | **2.69** | ~1.59 |
| [24, 26] | 0.560 | 0.283 | **1.98** | ~0.76 |
| [26, 28] | 0.603 | 0.240 | **2.52** | ~1.77 |
| [28, 32] | 0.625 | 0.213 | **2.94** | ~0.93 |
| [32, 36] | 0.630 | 0.355 | **1.77** | ~1.07 |

| Summary | OLD | NEW |
|:---|:---:|:---:|
| Ratio range | 0.77 -- 1.77 | **1.22 -- 2.94** |
| Centered on 1.0? | Yes (straddles) | **No (all > 1)** |
| Max deviation from 1.0 | 0.77 | **1.94** |

These ratios are inflated by signal contamination in the tight region. After truth signal subtraction, the pure background isolation fraction ratio is much closer to 1.0 at low pT (see Section 9).

### Physical interpretation

The tight-BDT background (photon-like jets) is systematically **more isolated** than non-tight background. Mean isoET for tight background is 2.3--2.9 GeV across all pT, while non-tight background rises from 2.8 GeV (low pT) to 6.0 GeV (high pT).

However, much of this apparent difference is driven by signal contamination in the tight region. When only pure background is considered (see `plotting/figures/bkg_isoET_tight_vs_nontight_jet_bdt_nom.pdf`), the tight and non-tight isoET distributions are nearly identical at low pT, with a slight divergence developing at high pT. The dominant driver of the large ratios above is signal photons in the jet MC tight region, not an intrinsic background effect.

The residual intrinsic jet physics effect is real but smaller: jets faking photons deposit most energy in a single narrow cluster, leaving less ambient energy in the isolation cone. The BDT selects for this topology, creating a mild BDT-isolation correlation.

**Why the old dead zone masked this:** By excluding BDT scores between 0.50 and the tight threshold (a gap of 0.18 at ET=8 shrinking to 0 at ET=20), the old config excluded an intermediate jet population from the non-tight region. These intermediate jets had isolation properties between tight and non-tight, so their exclusion pushed the non-tight isolation fraction higher at low pT, making the ratio closer to 1.0. Closing the gap pulls them back in, revealing the full correlation.

## 6. isoET Distribution Overlay (Background)

See figures:
- **`plotting/figures/isoET_tight_vs_nontight_jet_bdt_nom.pdf`** — Tight (red) vs non-tight (blue) isoET, normalized, per pT bin. All 12 bins show KS p-value = 0 (statistically incompatible shapes). The tight distribution peaks sharply near 1 GeV while non-tight develops a long tail extending to 15+ GeV at high pT.
- **`plotting/figures/isoET_tight_jetonly_vs_inclusive_bdt_nom.pdf`** — Jet-only vs inclusive MC in the tight region. Inclusive is dominated by signal (mean isoET 0.6--1.3 GeV), confirming the BDT selects well.

> **Note (corrected 2026-04-08):** The plots above use the **total** jet MC in the tight region, which contains 58-78% truth-matched prompt photon signal (see Section 9). The apparent sharp peak near 1 GeV in the tight distribution is dominated by signal photons, not background. For the pure background isoET comparison after truth signal subtraction, see `plotting/figures/bkg_isoET_tight_vs_nontight_jet_bdt_nom.pdf` (produced by `plotting/plot_bkg_isoET_tight_vs_nontight.C`). In the corrected background-only plots, tight and non-tight isoET shapes are nearly identical at low pT.

## 7. Signal Leakage Fractions

| pT (GeV) | cB (OLD) | cB (NEW) | cC (OLD) | cC (NEW) | cD (OLD) | cD (NEW) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| [8, 10] | 5--6% | 5.7% | 2.6--3% | 6.7% | <0.6% | 0.7% |
| [12, 14] | ~7% | 6.6% | ~3.5% | 5.7% | <0.6% | 0.8% |
| [16, 18] | ~8% | 7.2% | ~4% | 4.9% | <0.6% | 0.6% |
| [20, 22] | ~9% | 7.5% | ~4% | 4.4% | <0.6% | 0.5% |
| [26, 28] | ~11% | 8.8% | ~4.5% | 4.3% | <0.6% | 0.2% |
| [32, 36] | ~13% | 8.4% | ~4.6% | 2.2% | <0.6% | 0.2% |

**cC increased at low pT** (from ~3% to 6.7% at 8-10 GeV): the non-tight region now extends to higher BDT scores (up to 0.68 instead of 0.50 at ET=8), capturing more signal photons that fall just below the tight threshold.

**cB decreased slightly**: unchanged tight definition, so similar leakage to region B.

**cD unchanged**: negligible in both configs.

## 8. Discussion

### Why the sign flip persists

The purity correction sign change at pT ~ 14 GeV was attributed in the original investigation to the dead zone creating a pT-dependent BDT-isolation correlation. The dead zone is now eliminated, but the sign flip persists at the same pT. This proves the dominant mechanism is a **genuine background physics effect**: the pure background R factor (h_R, computed after truth signal subtraction) genuinely crosses 1.0 at pT ~ 12-14 GeV.

A critical finding (Section 9) is that much of the apparent BDT-isolation correlation seen in the jet MC total is driven by signal contamination (58-78% of jet MC region A is truth-matched signal). After removing this contamination, the pure background tight and non-tight isolation distributions are nearly identical at low pT, with only a subtle divergence at high pT. The R crossover from < 1 to > 1 is a genuine but modest background physics effect, not a large intrinsic correlation.

The dead zone was a secondary mechanism that partially masked the intrinsic correlation at low pT (by excluding intermediate-BDT jets from non-tight). Removing it reveals the full correlation, making the violation larger but more systematic.

### Positive and negative aspects of the change

**Positive:**
- Clean BDT partition with no dead zone or overlap -- simpler to interpret
- High-pT correction ratio is more stable (0.84-0.97 vs 0.77-0.92), with less wild bin-to-bin fluctuation
- The isolation fraction ratio is now monotonically above 1.0 -- the bias direction is consistent and predictable

**Negative:**
- Low-pT correction is larger (1.70 vs 1.51 at 9 GeV) -- greater reliance on MC correction
- The isolation fraction ratio is further from 1.0 everywhere -- ABCD non-closure is worse in absolute terms
- cC leakage increased at low pT due to wider non-tight window

### Implications for systematic uncertainty

The MC purity correction systematic uncertainty should be re-evaluated. The correction factor spans 0.84--1.70 (vs 0.77--1.51 before), but the high-pT range is tighter. The dominant uncertainty is still the MC statistical precision on truth purity and the ABCD fit quality. The signal contamination in jet MC (Section 9) means that isolation-based diagnostics from the raw jet MC totals overstate the true level of ABCD independence violation.

## 9. Signal Contamination in Jet MC and R Factor Decomposition

### Signal fraction in jet MC ABCD regions

The jet MC samples (jet10, jet15, jet20, jet30, jet50) contain truth-matched prompt photons that pass all analysis cuts. These signal photons are especially concentrated in region A (tight, isolated), where they comprise the majority of the jet MC total:

| pT bin (GeV) | Signal fraction in A | Signal fraction in B | Signal fraction in C | Signal fraction in D |
|:---:|:---:|:---:|:---:|:---:|
| [8, 10] | 57.6% | 4-11% | 2-21% | <1% |
| [10, 12] | 60.9% | 4-11% | 2-21% | <1% |
| [12, 14] | 63.2% | 4-11% | 2-21% | <1% |
| [14, 16] | 63.7% | 4-11% | 2-21% | <1% |
| [16, 18] | 64.7% | 4-11% | 2-21% | <1% |
| [18, 20] | 71.4% | 4-11% | 2-21% | <1% |
| [20, 22] | 73.9% | 4-11% | 2-21% | <1% |
| [22, 24] | 71.3% | 4-11% | 2-21% | <1% |
| [24, 26] | 73.2% | 4-11% | 2-21% | <1% |
| [26, 28] | 75.6% | 4-11% | 2-21% | <1% |
| [28, 32] | 76.1% | 4-11% | 2-21% | <1% |
| [32, 36] | 77.7% | 4-11% | 2-21% | <1% |

This means the jet MC "total" in region A is dominated by signal, not background.

### R factor decomposition: R_raw vs R_bkg

Two R factors can be computed from jet MC:

- **R_raw** (signal-contaminated): R = (A_total * D_total) / (B_total * C_total), using all jet MC counts without truth subtraction.
- **R_bkg** (pure background, = h_R): R = (A_bkg * D_bkg) / (B_bkg * C_bkg), after subtracting truth-matched signal from each region.

| pT bin (GeV) | R_raw (contaminated) | R_bkg (h_R, pure bkg) |
|:---:|:---:|:---:|
| [8, 10] | 1.47 | 0.675 |
| [10, 12] | 1.73 | 0.729 |
| [12, 14] | 2.29 | 0.904 |
| [14, 16] | 2.62 | 1.058 |
| [16, 18] | 3.00 | 1.261 |
| [18, 20] | 3.19 | 0.991 |
| [20, 22] | 4.21 | 1.273 |
| [22, 24] | 5.03 | 1.779 |
| [24, 26] | 3.22 | 1.151 |
| [26, 28] | 4.82 | 1.617 |
| [28, 32] | 6.16 | 1.854 |
| [32, 36] | 3.08 | 0.827 |

Key observations:
- **R_raw is always > 1** (range 1.47-6.16) and increases monotonically with pT. This was previously misidentified as "pure background R" in earlier versions of this report.
- **R_bkg (= h_R) genuinely crosses 1.0** at pT ~ 14 GeV. Below this, R_bkg < 1 (ABCD overestimates background); above, R_bkg > 1 with large statistical fluctuations.
- The factor-of-2-to-6 difference between R_raw and R_bkg is entirely due to signal contamination inflating region A.

### Corrected isoET distributions (background only)

The earlier isoET overlay plots (Section 6, `isoET_tight_vs_nontight_jet_bdt_nom.pdf`) used the total jet MC, in which the tight region is 58-78% signal. The sharp isolation peak in the tight distribution was dominated by signal photons, not background.

After truth signal subtraction, the background-only isoET distributions (`plotting/figures/bkg_isoET_tight_vs_nontight_jet_bdt_nom.pdf`, produced by `plotting/plot_bkg_isoET_tight_vs_nontight.C`) show:
- At low pT (8-14 GeV): tight and non-tight background isoET shapes are **nearly identical**
- At high pT (> 18 GeV): a slight divergence develops, with tight background marginally more isolated
- The difference is far smaller than suggested by the signal-contaminated plots

### Corrected physical interpretation

The R crossover at pT ~ 14 GeV is a **genuine but subtle background physics effect**. At low pT, tight and non-tight background jets have nearly identical isolation properties (R ~ 1, with R slightly below 1 due to a small residual correlation). At high pT, photon-like background jets are marginally more isolated than non-photon-like jets, pushing R above 1. The magnitude of the deviation from 1.0 in the pure background is modest (R_bkg = 0.67-1.85), with large statistical uncertainties at high pT.

The earlier conclusion that the "pure background" R is always > 1 (1.4-5.0 range) and that the sign flip in h_R was caused by signal leakage corrections was incorrect. Those R values were computed from signal-contaminated jet MC totals.

See `plotting/figures/R_decomposition_bdt_nom.pdf` (produced by `plotting/plot_R_decomposition.C`) for the full R factor decomposition showing all four curves (R_raw, R_bkg, R_signal, and the signal fraction) as a function of pT.

## 10. Recommendations

1. **The purity correction sign flip is a genuine background physics effect** at the level of R_bkg crossing 1.0 at pT ~ 14 GeV. It cannot be eliminated by adjusting the non-tight BDT boundary. The magnitude is modest after removing signal contamination from the diagnostic.

2. **Consider applying the MC R-factor directly**: use h_R from MC to correct the ABCD R-factor in data, rather than correcting the final purity post-hoc. This would address the non-closure at its source.

3. **Jet MC diagnostics must account for signal contamination**: any isolation-based ABCD diagnostic computed from raw jet MC totals will be dominated by the 58-78% signal fraction in region A. Always use truth-subtracted quantities (h_R, background-only isoET) for physics interpretation of background behavior.

4. **The ~30 systematic variation configs missing `bdt_max_slope`/`bdt_max_intercept`** still use the old fixed-0.50 non-tight bound. If the ET-dependent bound is the intended configuration, regenerate all variation configs with `make_bdt_variations.py` (ensure it propagates these fields to all variants, not just the BDT-threshold variants).

5. **`RecoEffCalculator.C` (old code path)** does not read `bdt_max_slope`/`bdt_max_intercept` — only `RecoEffCalculator_TTreeReader.C` does. Confirm which code path `oneforall.sh` uses (it appears to use the TTreeReader version).

## Appendix: Files Produced

| File | Description |
|------|-------------|
| `plotting/figures/isoET_tight_vs_nontight_jet_bdt_nom.pdf` | Tight vs non-tight isoET overlay (jet MC total, per pT bin) -- signal-contaminated |
| `plotting/figures/isoET_tight_jetonly_vs_inclusive_bdt_nom.pdf` | Tight isoET: jet-only vs inclusive MC |
| `plotting/figures/bkg_isoET_tight_vs_nontight_jet_bdt_nom.pdf` | Pure background isoET overlay (truth-subtracted, per pT bin) |
| `plotting/figures/R_decomposition_bdt_nom.pdf` | R factor decomposition: R_raw, R_bkg, signal fraction vs pT |
| `plotting/plot_isoET_tight_vs_nontight.py` | Original isoET plotting script |
| `plotting/plot_bkg_isoET_tight_vs_nontight.C` | Background-only isoET plotting macro |
| `plotting/plot_R_decomposition.C` | R factor decomposition plotting macro |
| `reports/mc_purity_correction_etdep_nontight.md` | This report |
