# RHIC Accelerator

## Overview

The Relativistic Heavy Ion Collider (RHIC) at Brookhaven National Laboratory provides the pp collisions at sqrt(s) = 200 GeV analyzed in PPG12. RHIC consists of two independent superconducting storage rings (Blue and Yellow) sharing a common tunnel of 3.83 km circumference, with six interaction regions. The machine is capable of accelerating heavy ions up to gold (Au) and polarized protons, with maximum center-of-mass energies of sqrt(s_NN) = 200 GeV for Au+Au and sqrt(s) = 510 GeV for pp.

Reference: NIM A 499 (2003) 235 (Harrison, Ludlam, Ozaki)

## Run 24 Parameters for PPG12

PPG12 uses pp collision data from the RHIC 2024 Run. The run is divided into two periods distinguished by beam crossing angle:

### Crossing Angle Periods

| Parameter | 0 mrad period | 1.5 mrad period |
|-----------|---------------|------------------|
| Run range | 47289--51274 | 51274--54000 |
| Crossing angle | 0 mrad (head-on) | 1.5 mrad |
| Integrated luminosity | 32.66 pb^-1 | 16.86 pb^-1 |
| Double-interaction fraction | 18.7% | 7.2% |
| Pileup rate | Higher | Lower |

The crossing angle change at run 51274 significantly reduced the pileup rate. The 1.5 mrad crossing angle reduces the luminous region overlap between bunches, lowering the probability of two independent pp collisions in the same bunch crossing from 18.7% to 7.2%.

### Combined Run Statistics

| Parameter | Value |
|-----------|-------|
| Total luminosity (all runs) | ~49.6 pb^-1 |
| 1.5 mrad luminosity (nominal PPG12 dataset) | 16.86 pb^-1 |
| Center-of-mass energy | sqrt(s) = 200 GeV |
| Species | proton + proton (unpolarized for PPG12) |

PPG12 primarily uses the 1.5 mrad period (16.86 pb^-1) as the nominal dataset due to its lower pileup contamination. The 0 mrad period can be included for cross-checks with appropriate pileup corrections, bringing the total to ~50 pb^-1.

## Luminosity Determination

The integrated luminosity is determined by counting minimum-bias events recorded by the MBD trigger, corrected for prescale factors and vertex reconstruction criteria:

```
L_int = N_MB / sigma_MBD
```

where:
- N_MB is the prescale-corrected number of MBD-triggered events
- sigma_MBD = 25.2 +2.3/-1.7 mb is the MBD trigger cross-section in pp at sqrt(s) = 200 GeV

The MBD trigger cross-section was measured via Vernier scans during Run 24 (documented in the PPG09 IAN, sPHENIX-Invenio). The ~9% uncertainty on sigma_MBD is one of the dominant systematic uncertainties in the PPG12 cross-section.

## Bunch Structure and Pileup

RHIC operates with 111 filled bunches per ring (out of 120 possible slots), with a bunch-crossing rate of approximately 9.4 MHz. The sPHENIX calorimeter electronics digitize at 56.4 MHz (6x the crossing rate), providing multiple samples per crossing for waveform analysis.

### Double Interactions

At the luminosities achieved in Run 24, there is a significant probability of two independent pp collisions occurring in the same bunch crossing. This "pileup" or "double interaction" produces:

- A reconstructed vertex that is the average of the two collision vertices
- Distorted cluster kinematics (ET, eta) from the shifted vertex
- Additional energy in the isolation cone from the second collision
- Modified shower-shape BDT scores from the vertex-displaced reconstruction

The double-interaction fractions (18.7% at 0 mrad, 7.2% at 1.5 mrad) are derived from RHIC beam parameters and validated by MBD timing measurements.

### Pileup Rejection

Several handles exist to identify and reject pileup events:

1. **INTT timing:** The silicon strip tracker provides 60 ns timing resolution, enabling single bunch-crossing readout to suppress out-of-time pileup.

2. **MBD timing spread:** The timing spread (avg sigma_t) across the 128 MBD PMTs is broader for pileup events because the two collision vertices produce different time-of-flight distributions. The `mbd_avg_sigma_max` config cut rejects events with excessive timing spread.

3. **Vertex distribution:** Pileup events tend to have vertices closer to z = 0 (average of two uniformly distributed vertices) than single-interaction events.

## MBD Trigger Detector

The Minimum Bias Detector (MBD, formerly the PHENIX Beam-Beam Counter or BBC) is the primary trigger and luminosity detector:

| Parameter | Value |
|-----------|-------|
| Location | +/-250 cm from IP along beam axis |
| Coverage | 3.51 < \|eta\| < 4.61 |
| Channels | 128 PMTs (64 per side), 3 concentric rings, full 2pi |
| Radiator | Quartz (Cherenkov) |
| Timing resolution | ~50 ps per PMT |
| Trigger condition | >= 2 PMTs fired per side (north AND south coincidence) |
| Trigger efficiency | 92--94% for inelastic pp (Vernier scan measurement) |
| Vertex resolution | ~2 cm (from north-south time difference) |

The MBD was originally designed and built for PHENIX (NIM A 499 (2003) 549), then inherited by sPHENIX with new readout and trigger electronics. It sits in the same location but now has updated coverage (3.51 < |eta| < 4.61 vs. the original PHENIX 3.0 < |eta| < 3.9).

### MBD Pileup Metric

The MBD timing information provides a pileup metric computed in `MbdPileupHelper.h`:

- Select PMTs with |time| < 25 ns, charge > 0.4 MIP, and >= 2 hits per side
- Compute the timing RMS for north and south arms separately
- **avgsigma** = average of north and south timing RMS (primary metric)
- Additional metrics: prodsigma, maxsigma, proddelta, avgdelta, maxdelta

Events with large avgsigma are flagged as likely pileup and rejected by the `mbd_avg_sigma_max` cut.

## Accelerator Complex

The RHIC injection chain:

1. **OPPIS/LINAC:** Polarized proton source and linear accelerator (for pp)
2. **Booster:** First circular accelerator stage
3. **AGS (Alternating Gradient Synchrotron):** Accelerates to injection energy for RHIC
4. **RHIC rings:** Two independent superconducting rings (Blue and Yellow)
   - 1,740 superconducting magnets operating at 4.6 K
   - Six interaction regions (sPHENIX at IR-8)

For heavy-ion running, the injection chain uses the Tandem Van de Graaff or EBIS as ion sources instead of OPPIS.

## Relevance to PPG12

The RHIC accelerator parameters directly affect PPG12 through:

1. **Luminosity:** Determines the statistical reach of the measurement. The cross-section is normalized by L_int, so the ~9% sigma_MBD uncertainty propagates directly.

2. **Pileup rate:** The crossing-angle-dependent double-interaction fraction (18.7% vs 7.2%) motivates the nominal choice of the 1.5 mrad period and the systematic studies of pileup effects on cluster kinematics and isolation.

3. **Bunch structure:** The 9.4 MHz crossing rate and INTT 60 ns timing determine the pileup suppression capability.

4. **Vertex distribution:** The longitudinal vertex distribution (driven by bunch length) affects the acceptance correction and the vertex-dependent kinematics recalculation. The |z_vtx| < 10 cm cut ensures tracking acceptance.

5. **Run-to-run stability:** Energy scale drifts and detector conditions vary across the ~7000-run dataset, requiring the run-range filtering (`run_min`/`run_max` in config) and the MBD T0 corrections (`mbd_t0_correction_file`).

## References

- RHIC overview: NIM A 499 (2003) 235 (Harrison, Ludlam, Ozaki)
- MBD (BBC) design: NIM A 499 (2003) 549 (Allen et al.)
- PPG09 luminosity IAN: sPHENIX-Invenio
- BaBar solenoid: IEEE Trans. Appl. Supercond. 9 (1999) 847

## Cross-Links

- [sPHENIX Overview](sphenix-overview.md) -- detector built at RHIC IR-8
- [Simulation Tools](simulation-tools.md) -- MC production uses RHIC beam conditions
- [Tree Making](../../pipeline/01-tree-making.md) -- run-range and vertex cuts from RHIC conditions
- [MC Samples](../../reference/mc-samples.md) -- MC produced at Run 28 conditions
