# ABCD Background Subtraction

## Region Definitions

Two orthogonal discriminators divide clusters into four regions:

```
            Isolated (iso)       Non-isolated (noniso)
Tight       A = signal region    B = bkg template
Non-tight   C = bkg control      D = bkg normalization
```

- **Tight vs Non-tight:** BDT score threshold (parametric: `bdt_min = intercept + slope * ET`)
- **Isolated vs Non-isolated:** Isolation ET within parametric cone cut (`iso_max = b + s * ET`)

## Histogram Naming

In ROOT files produced by `RecoEffCalculator_TTreeReader.C`:

| Region | Data | Signal MC | Background MC |
|--------|------|-----------|---------------|
| A (tight, iso) | `h_tight_iso_cluster_0` | `h_tight_iso_cluster_signal_0` | `h_tight_iso_cluster_notmatch_0` |
| B (tight, noniso) | `h_tight_noniso_cluster_0` | `h_tight_noniso_cluster_signal_0` | `h_tight_noniso_cluster_notmatch_0` |
| C (nontight, iso) | `h_nontight_iso_cluster_0` | `h_nontight_iso_cluster_signal_0` | `h_nontight_iso_cluster_notmatch_0` |
| D (nontight, noniso) | `h_nontight_noniso_cluster_0` | `h_nontight_noniso_cluster_signal_0` | `h_nontight_noniso_cluster_notmatch_0` |

The `_0` suffix is the eta bin index.

## Signal Leakage

True signal photons can leak into non-signal regions. Leakage fractions (computed from merged signal MC):

```
cB = B_sig / A_sig    (h_tight_noniso_cluster_signal / h_tight_iso_cluster_signal)
cC = C_sig / A_sig    (h_nontight_iso_cluster_signal / h_tight_iso_cluster_signal)
cD = D_sig / A_sig    (h_nontight_noniso_cluster_signal / h_tight_iso_cluster_signal)
```

## Self-Consistent ABCD Equation

With leakage corrections, the signal yield in region A satisfies:

```
NsigA = NA - R * (NB - cB * NsigA) * (NC - cC * NsigA) / (ND - cD * NsigA)
```

where R = 1 (the double ratio is embedded in the formula).

This is solved numerically by finding the root of `f(NsigA)`:
```
f(NsigA) = NsigA - (NA - R * (NB - cB*NsigA) * (NC - cC*NsigA) / (ND - cD*NsigA))
```

Implemented as `myfunc()` in `CalculatePhotonYield.C`.

## R Factor

Two calculations:

1. **Data-driven R:** Background yields estimated as `data - sim_signal_leak` per region.
   ```
   R = (A_bkg / B_bkg) * (D_bkg / C_bkg)
   ```

2. **MC truth R:** Uses `_notmatch` histograms directly.
   ```
   R_notmatch = (A_notmatch / B_notmatch) * (D_notmatch / C_notmatch)
   ```

## Purity Extraction (Toy MC)

For each pT bin, 20,000 Poisson toy experiments:

1. Compute effective events: `A_eff = A^2 / sum(w^2)`
2. Draw Poisson: `NA_eff = Poisson(A_eff)`, rescale: `NA = A * NA_eff / A_eff`
3. Solve ABCD equation (without leakage) -> purity = NsigA / NA
4. Repeat with Gaussian-smeared leakage (cB, cC, cD) -> purity with leakage correction
5. Fit each purity distribution with Gaussian -> central value and uncertainty

## Purity Fitting

The bin-by-bin purity values are smoothed:

- `fit_option=0` (default): Error function `[0]*Erf((x-[1])/[2])` over [8, 32] GeV
- `fit_option=1`: Pade rational function `([0]+[1]*x)/(1+[2]*x)`

Confidence intervals at 68.3% CL via `TVirtualFitter::GetConfidenceIntervals()`.

The `fittingerror` config parameter shifts purity by +/- 1 confidence interval for systematic studies.

## MC Purity Correction

When `mc_purity_correction = 1`:
1. Run MC closure (isMC=true) to get truth purity
2. Fit MC truth purity with same functional form
3. Compute ratio: `mc_corr = truth_purity / fitted_purity`
4. Apply: `purity_corrected = fitted_purity * mc_corr`

## Parametric Thresholds

Both the tight BDT and isolation cuts are ET-dependent (not flat):

**BDT threshold:**
```
bdt_min = bdt_min_intercept + bdt_min_slope * ET
```
Nominal: intercept = 0.80, slope = -0.015

**Isolation cut:**
```
iso_max = reco_iso_max_b + reco_iso_max_s * ET
```
Nominal: b = 0.502095, s = 0.0433036

## Key Config Fields

| Field | Purpose | Nominal |
|-------|---------|---------|
| `tight.bdt_min_intercept` | Tight BDT threshold intercept | 0.80 |
| `tight.bdt_min_slope` | Tight BDT threshold slope | -0.015 |
| `non_tight.bdt_max_intercept` | Non-tight upper BDT (should match tight min) | 0.80 |
| `non_tight.bdt_max_slope` | Non-tight upper BDT slope | -0.015 |
| `analysis.reco_iso_max_b` | Isolation intercept | 0.502095 |
| `analysis.reco_iso_max_s` | Isolation slope | 0.0433036 |
| `analysis.reco_noniso_min_shift` | Gap between iso and noniso boundaries | 0.8 |
| `analysis.n_nt_fail` | Number of tight cuts to fail for non-tight | 0 |

## Code Locations

- ABCD equation: `efficiencytool/CalculatePhotonYield.C` (`myfunc()`)
- Leakage computation: `CalculatePhotonYield.C`
- Cut application: `RecoEffCalculator_TTreeReader.C` (tight/nontight classification)
- Purity fitting: `CalculatePhotonYield.C`
