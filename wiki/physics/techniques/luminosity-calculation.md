# Luminosity Calculation

## Overview

Converting raw photon yields into a differential cross-section requires dividing by integrated luminosity. In sPHENIX pp collisions, luminosity is derived from MBD (Minimum Bias Detector) scaler counts, corrected for trigger prescales and pileup. This article synthesizes the sPHENIX luminosity calculation procedure and its implications for PPG12.

Source: sPHENIX internal note (`wiki/raw/papers/lumi_calculation-5.pdf`), full summary at `wiki/raw/papers/D_techniques/lumi_calculation.md`.

## Luminosity Formula

For the Jet10 trigger (no MBD coincidence required):

```
L_int = N_live_mbd × F_corr / (P_Jet10 × σ_MBD)
```

| Symbol | Meaning | Value |
|--------|---------|-------|
| N_live_mbd | Live MBD N&S>=1 scaler count | Per-run from DSTs |
| F_corr | Pileup correction factor | ~1.0 at 1.5 mrad, ~1.05-1.10 at 0 mrad |
| P_Jet10 | Average prescale of Jet10 trigger | = live_scalers / scaled_scalers |
| σ_MBD | MBD cross-section (Vernier scan) | 25.2 mb |

The prescale accounts for the DAQ scaledown: `prescale = scaledown + 1`.

## Pileup Correction

Multiple pp collisions per bunch crossing follow Poisson statistics with mean:

```
mu = R_true / f,    f = 111 bunches × 78.196 kHz = 8.86 MHz
```

The correction factor upweights multi-interaction crossings:

```
F_corr = <n | n>=1> = sum(n × P(n,mu)) / sum(P(n,mu))   for n >= 1
```

To compute F_corr, one needs the mapping from the observed MBD rate to the true collision rate. This uses GEANT4-derived MBD firing probabilities:

| Single collision outcome | Probability |
|--------------------------|-------------|
| Fire both N&S | 0.562 |
| Fire N only | 0.175 |
| Fire S only | 0.175 |
| Fire neither | 0.088 |

For two collisions, a 4×4 probability matrix is evaluated — any combination firing at least one PMT on each side triggers MBD N&S>=1. These probabilities are z-vertex dependent, dropping sharply for |z_vtx| > 100 cm.

**Impact by crossing angle:**
- **0 mrad**: Higher rates → mu ≈ 0.2-0.3, pileup correction significant (~5-10%)
- **1.5 mrad**: Lower rates → mu ≈ 0.05-0.1, correction nearly negligible

## Z-Vertex Selection Complication

When the analysis applies an MBD coincidence + z-vertex cut (as PPG12 does), the **numerator** (photon counts) also needs correction. Two pileup collisions can conspire to:
- Satisfy the MBD coincidence when neither alone would
- Shift the reconstructed vertex into the fiducial range when neither collision vertex was inside it

The numerator correction factor is:

```
f_corr = [1/(1-Prob0)] × [Prob1 × Pass1 + Prob2 × Pass2 + Prob3 × Pass3 + ...]
```

where Pass_k is the probability that a signal event + (k-1) overlay mb events passes the MBD + z-vtx selection. This correction is:
- **Probe-dependent** (depends on photon pT)
- **Segment-dependent** (z-vtx distribution changes within a store)
- Practically, only the first ~3 terms matter at any pp luminosity (Prob4 < 10^-4 for mu < 0.3)

## Run-to-Run QA

The note establishes Njet/Lint vs run number as a QA metric (Figures 3-6). Key findings:
- After pileup correction, Njet/Lint is approximately flat across the run
- 0 mrad and 1.5 mrad periods show ~10% difference in mean Njet/Lint, consistent with z-vertex distribution effects
- Runs 51400-52200 show elevated Njet/Lint (mean ~7660 vs ~6550 overall)

This technique could be adapted for photon yield stability (Nphoton/Lint vs run).

## PPG12 Implementation

| Component | How luminosity enters |
|-----------|-----------------------|
| `lumi/60cmLumi_fromJoey.list` | Per-run luminosity (8 trigger columns), pileup correction already applied in "Corr" columns |
| `LumiCalculator.C` | Parses lumi list, sums over run range, returns Bit30Corr total |
| Config `analysis.lumi` | Stores pre-computed total for the run range (e.g., 16.86 pb^-1 for 1.5 mrad) |
| `CalculatePhotonYield.C` | Reads `analysis.lumi`, scales yields: `dσ/dpT = N / (ΔpT × L_int)` |
| `calc_pileup_range.C` | Per-run L_int = N_mbdc / σ_MB for pileup studies |

### Luminosity values used in PPG12

| Period | Run range | L_int (pb^-1) | Source |
|--------|-----------|---------------|--------|
| 0 mrad | 47289-51274 | 32.66 | Bit30Corr sum |
| 1.5 mrad | 51274-54000 | 16.86 | Bit30Corr sum |
| All | 47289-54000 | 49.56 | Bit30Corr sum |

### Luminosity systematic uncertainty

From `make_bdt_variations.py`:
- Central: 25.2 mb (σ_MBD)
- Down: 23.5 mb → -6.7%
- Up: 27.5 mb → +9.1%

This asymmetric uncertainty propagates as a flat fractional uncertainty on the cross-section.

## Comparison with Other Experiments

- **PHENIX**: Used BBC (Beam-Beam Counter) cross-section from Vernier scan, similar Poisson pileup correction. PPG12's MBD is the successor detector with improved z-vertex resolution.
- **ATLAS/CMS**: Use van der Meer scans (equivalent to Vernier scans) with LUCID/BCM/HF luminosity detectors. Pileup corrections are more complex at LHC (mu ~ 20-50 vs < 0.3 at RHIC).
- **STAR**: Similar RHIC conditions, uses ZDC and VPD for luminosity measurement.
