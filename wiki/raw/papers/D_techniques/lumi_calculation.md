---
title: "sPHENIX Luminosity Calculation"
authors: "sPHENIX Collaboration (internal note)"
category: D_techniques
priority: 3
tags: [luminosity, pileup, MBD, Poisson, scaler, trigger, prescale, vernier-scan, double-interaction]
ppg12_relevance: "Defines the luminosity calculation method used in PPG12 cross-section normalization, including pileup correction and MBD z-vertex selection effects"
source_pdf: "wiki/raw/papers/lumi_calculation-5.pdf"
---

# sPHENIX Luminosity Calculation

## Summary

This internal sPHENIX note documents the full procedure for calculating integrated luminosity from pp collision data. It covers three main topics: (1) extracting raw luminosity from DST scalers and trigger prescales, (2) correcting for pileup (multiple collisions per bunch crossing) using a Poisson model with GEANT4-derived MBD firing probabilities, and (3) handling the more complex case where an MBD coincidence and z-vertex selection are required.

## 1. Basic Luminosity from Scalers

### Input quantities per segment per trigger

Nine values are extracted from each DST segment for each trigger:
1. Raw scalers at beginning/end of segment
2. Live scalers at beginning/end of segment
3. Scaled scalers at beginning/end of segment
4. Number of triggers fired in all produced events
5. Start/end BCO (beam crossing offset, single number per segment)

### Average prescale

Summing over the entire run (~60 minutes):

```
prescale = total_live_scalers / total_scaled_scalers
```

Note: sPHENIX DAQ "scaledown" has the relation `prescale = scaledown + 1` (scaledown = 0 means every live event is recorded).

### Sampled MBD collisions

```
N_mbd_sampled = N_live_MBD(N&S>=1) / prescale
```

If the physics trigger (e.g., Jet10) does **not** require an MBD N&S>=1 coincidence, divide by the fraction of the 42 mb pp inelastic cross-section that fires MBD N&S>=1.

## 2. Pileup Correction

### Poisson model

All collisions assumed uncorrelated → Poisson distribution:

```
P(n, mu) = exp(-mu) * mu^n / n!
```

Mean collisions per bunch crossing:

```
mu = R_true / f
f = 111 bunches × 78.196 kHz = 8.86 MHz
```

### MBD firing probabilities (from PYTHIA 8 + GEANT4)

**Single collision (Table 1):**

| Outcome | Probability |
|---------|-------------|
| Fire both N&S | 0.562 |
| Fire N only | 0.175 |
| Fire S only | 0.175 |
| Fire neither | 0.088 |

These probabilities are z-vertex dependent (see Figure 2 in the note). The MBD efficiency for firing both sides drops off sharply for |z_vtx| > ~100 cm.

**Two collisions (Table 2):** 4×4 matrix of joint probabilities for (collision 1 outcome) × (collision 2 outcome). Key: any combination that fires at least one tube on each side triggers MBD N&S>=1 (starred entries in the table).

### Mapping R_raw_MBD to R_true

```
R_raw_MBD = R_true * sum_n [ P(MBD|n) * P(n, mu) ]
```

where P(MBD|n) is the probability of firing MBD N&S>=1 given n collisions. A lookup table from R_raw_MBD → R_true is built by evaluating many R_true values.

### Correction factor

```
F_corr = sum_{n=1}^{inf} [ n * P(n, mu) ] / sum_{n=1}^{inf} [ P(n, mu) ]
```

This upweights multi-interaction bunch crossings by the number of interactions.

### Corrected luminosity

```
L_int = F_corr * L_int(uncorrected)
```

## 3. Final Luminosity Formula

For a trigger that does **not** require MBD coincidence (e.g., Jet10 alone):

```
L_int = N_live_mbd * F_corr / (P_Jet10 * sigma_MBD)
```

where:
- `N_live_mbd` = live MBD N&S>=1 triggers in the good run list segments
- `P_Jet10` = average prescale of the Jet10 trigger
- `sigma_MBD` = 25.2 mb (from sPHENIX Vernier Scan)

For minimum bias coincidence only:

```
L_int = N_scaled_mbd * F_corr / sigma_MBD
```

The difference between "MB + Jet10" and "Jet10 alone" luminosity arises only when their prescales differ.

### Available luminosity lists

Published at: `github.com/sPHENIX-Collaboration/analysis/blob/master/LuminosityCounterGoodRuns/run/allzlumi.list`

These are calculated **without** a z-vertex selection (no vertex cut in denominator).

## 4. Run-to-Run QA

### Njet / Lint stability check

Figures 3-6 in the note show Njet/Lint vs run number as a QA metric:
- Jets defined: dijet + timing cuts (Δφ > 3π/4, |t_lead + 2ns| < 6ns, |Δt| < 3ns), uncalibrated pT > 15 GeV
- 0th-order polynomial fit across all runs

### Key observations

- **Pileup correction is significant only during 0 mrad running** (runs 51400-53864)
- Remaining non-statistical deviations likely from jet selection / energy calibration
- **0 mrad vs 1.5 mrad mean values differ** (Figure 6): ~10% effect in data, ~15% in simulation (from truth z-vtx within 10 cm vs 50 cm comparison)
- Runs 51400-52200 have notably higher Njet/Lint (mean 7660 pb vs 6550 pb overall)
- Separate 0 mrad / 1.5 mrad calculations at low pT may reveal a systematic uncertainty

## 5. MBD Z-Vertex Selection Case

When the analysis requires an MBD coincidence **and** a z-vertex cut (as in photon analyses), the correction becomes more complex because the numerator (particle counts) also needs correction.

### Algorithm

For each segment with collision rate R:

1. Calculate Poisson probabilities: Prob1, Prob2, Prob3, ... for 1, 2, 3, ... collisions
2. Using PYTHIA 8 events (γ for signal, mb for overlay):
   - **Pass1**: probability that a single γ event passes MBD coincidence + z-vtx cut
   - **Pass2**: probability that γ + 1 mb overlay passes MBD coincidence + z-vtx cut
   - **Pass3**: probability that γ + 2 mb overlays passes
3. Correction factor for the numerator:

```
f_corr = [1 / (1 - Prob0)] × [Prob1 × Pass1 + Prob2 × Pass2 + Prob3 × Pass3 + ...]
```

where `1 - Prob0 = Prob1 + Prob2 + ...`

### Practical considerations

- For R = 2 MHz (mu = 0.225): Prob1 = 0.179, Prob2 = 0.020, Prob3 = 0.0015, Prob4 = 8.2e-5 → first three terms sufficient at any pp luminosity
- The correction is **probe-dependent** (depends on pT of the photon/jet)
- The correction is **segment-dependent** (z-vtx distribution changes store-to-store and within a store)
- Two MBD collision events can conspire to:
  - Satisfy MBD coincidence when neither alone would
  - Produce a z-vertex near zero even when neither collision was in the fiducial range

## Key Numbers for PPG12

| Parameter | Value |
|-----------|-------|
| MBD cross-section (σ_MBD) | 25.2 mb |
| pp inelastic cross-section | 42 mb |
| Bunch crossing frequency | 8.86 MHz (111 bunches × 78.196 kHz) |
| Single-collision MBD N&S efficiency | 0.562 |
| Single-collision MBD-neither probability | 0.088 |
| Vernier scan source | sPHENIX Vernier Scan (reference TBD) |

## PPG12 Connection

- `LumiCalculator.C` implements this procedure: parses per-run scaler data from `60cmLumi_fromJoey.list`, sums by trigger bit (Bit10/18/22/30) with corrected/uncorrected columns, returns Bit30Corr as primary luminosity
- `CalculatePhotonYield.C` reads `analysis.lumi` from config and uses it to convert yields to differential cross-sections: `dσ/dpT = N / (Δ pT × L_int)`
- The pileup correction (F_corr) is already applied in the luminosity list — the "Corr" columns include it
- The Njet/Lint QA approach (Figures 3-6) could be adapted for photon yield stability checks
- Section 3 (MBD z-vtx case) is relevant for understanding PPG12's vertex-dependent efficiency corrections
- `calc_pileup_range.C` uses σ_MB = 25.2 × 10^9 pb to compute L_int = N_mbdc / σ_MB per run
