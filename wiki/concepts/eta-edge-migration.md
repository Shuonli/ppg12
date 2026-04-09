# Eta Edge Migration and Fiducial Boundary Effects

## Summary

Truth photons outside the fiducial |eta| < 0.7 can have reconstructed clusters that migrate inside the acceptance, and vice versa. In the unfolding framework, inward-migrating clusters should be treated as fakes but are currently entering the response matrix as signal in `RecoEffCalculator_TTreeReader.C`.

## The Problem

Both the truth-level fiducial cut and the reco-level cluster cut use `eta_bins: [-0.7, 0.7]` from the config. In `RecoEffCalculator_TTreeReader.C`:

1. Truth photons are registered in `photon_reco` map at line 1526 **before** the eta check at line 1538
2. If the truth photon fails the eta check (`etabin == -1`), `continue` skips truth histograms but the photon remains in the map
3. When a reco cluster matches this truth photon, it passes the `photon_reco.find()` check at line 2514
4. The matched pair enters the response matrix at line 2588 using the cluster's reco eta bin
5. The truth pT from outside fiducial contaminates the response matrix

## Study: EtaMigrationStudy.C

A dedicated macro (`efficiencytool/EtaMigrationStudy.C`) quantifies this effect:

- Loops over photon MC samples (photon5, photon10, photon20)
- Registers truth photons with **no eta filter** (full CEMC range)
- Classifies each truth-reco matched pair into 4 categories:
  - **Good**: both truth and reco inside |eta| < 0.7
  - **Fake (inward)**: truth outside, reco inside -- should be unfolding fakes
  - **Lost (outward)**: truth inside, reco outside -- missed signal
  - **Both outside**: neither in fiducial

### Single Interaction Results

| Metric | Value |
|--------|-------|
| Overall inward migration rate | 1.2% |
| Signal region (tight+iso) contamination | 1.24% |
| High-pT contamination (32-36 GeV) | 2.57% |
| Outward loss rate | 0.9% |
| Region C (nontight+iso) contamination | up to 52% at high pT |

The contamination rises with pT because higher-pT clusters at the boundary have more energy and are more likely to reconstruct inside fiducial.

### Double Interaction Results

| Metric | Single | Double | Ratio |
|--------|--------|--------|-------|
| Overall raw inward rate | 1.2% | 13.7% | 11.4x |
| Signal region contamination | 1.24% | 0.34% | 0.27x |
| Region C contamination | up to 52% | -- | -- |

The vertex shift in double interaction amplifies raw migration by ~11x, but the tight+iso selection is **more effective** at rejecting double-interaction fakes (0.34% vs 1.24%). The shifted vertex distorts shower shapes enough that the BDT rejects them.

### Signal Leakage Comparison

The signal leakage plot (`eta_migration_signal_leakage.pdf`) shows the fraction of inward-migrating fakes that pass tight+iso (leak into region A) vs those absorbed into region C (nontight+iso), as a function of ET, comparing single and double interaction:

- **Single interaction**: ~5-10% of inward fakes leak into region A (tight+iso). The fakes come from truth photons barely outside eta=0.7, so their shower shapes are similar to fiducial photons.
- **Double interaction**: Dramatically lower region A leakage. The vertex shift adds additional shower shape distortion on top of the edge effect, so the BDT rejects them more aggressively into region C.
- **Region C** absorbs the bulk of edge fakes in both cases, with double interaction having an even larger C fraction than single.

This demonstrates that the tight BDT selection discriminates between "near-boundary" fakes (single, similar shower shapes) and "vertex-shifted" fakes (double, compounded distortion).

**Caveat**: Double interaction tight+iso statistics are marginal (77 unweighted events total).

## Physics Mechanism

The eta migration arises from:

1. **Cluster position resolution**: The reconstructed cluster centroid differs from the true photon direction due to shower fluctuations and finite tower granularity
2. **Vertex resolution**: The reconstructed vertex differs from truth, shifting the calculated cluster eta via `eta = asinh(z/R_EMCal)` where `R_EMCal = 93.5 cm`
3. **Double interaction vertex shift**: Pileup averages two collision vertices, amplifying the eta displacement by `delta_eta ~ -dz / (R_EMCal * cosh(eta))`

The effect is asymmetric: inward migration (fake rate ~1.2-6.5% per bin) exceeds outward migration (loss rate ~0.9-2%) at high pT because the CEMC extends to |eta| ~ 1.1, providing a larger source population outside the 0.7 boundary.

## Impact on Analysis

### Response Matrix
The 1.24% signal region contamination directly enters the response matrix as spurious (reco ET, truth pT) entries where truth pT is from outside the fiducial region. This biases the unfolding by associating incorrect truth pT values with reconstructed clusters.

### ABCD Background Estimation
Region C (nontight+iso) contamination reaching 52% at high pT could affect the R-factor used in ABCD background subtraction. Edge-migrated clusters preferentially fail the tight BDT selection.

## Remediation Options

1. **Response matrix fix**: Check truth eta when filling the response matrix; count outside-fiducial matches as `RooUnfoldResponse::Fake()` instead of `Fill()`
2. **Tighter reco cut**: Apply reco |eta| < 0.6 with truth |eta| < 0.7, creating a buffer zone
3. **Systematic uncertainty**: Evaluate impact on unfolded cross-section with/without the fix

## Key Files

| File | Purpose |
|------|---------|
| `efficiencytool/EtaMigrationStudy.C` | Migration study macro |
| `efficiencytool/run_eta_migration.sh` | Orchestration script |
| `plotting/plot_eta_migration.C` | Plotting macro (9 diagnostic plots, incl. signal leakage comparison) |
| `efficiencytool/reports/eta_migration_study.tex` | Full LaTeX report |
| `plotting/figures/eta_migration_*.pdf` | Generated plots |

## Output ROOT Files

- Single: `results/eta_migration_signal_{var_type}.root`
- Double: `results/eta_migration_double_signal_{var_type}.root`
- Per-sample: `results/eta_migration_{photon5,photon10,photon20}_{var_type}.root`

## See Also

- [Shower Shape Variables](shower-shape-variables.md) -- BDT features that suppress edge fakes
- [Double-Interaction Efficiency](double-interaction-efficiency.md) -- Vertex shift effects on cluster kinematics
- [Systematic Variations](systematic-variations.md) -- Framework for adding edge migration as a systematic
- [Constants Sync](../reference/constants-sync.md) -- eta_bins must be consistent across all macros
