# sPHENIX Detector Overview

## From PHENIX to sPHENIX

sPHENIX is a major upgrade to the PHENIX experiment at the Relativistic Heavy Ion Collider (RHIC), designed to exploit the full jet physics capabilities of RHIC's Au+Au and pp collision program. While PHENIX had limited calorimetric coverage (two arms, each spanning ~pi/2 in azimuth and |eta| < 0.35), sPHENIX provides hermetic calorimetry and tracking over |eta| < 1.1 and full 2pi in azimuth. This transformation from a specialized single-arm spectrometer to a barrel detector enables jet reconstruction, gamma-jet correlations, and the isolated photon cross-section measurement of PPG12.

The upgrade proposal (arXiv:1501.06197) identified three pillars of the physics program:

1. **Jet quenching** -- inclusive jet R_AA, dijets, gamma-jet, fragmentation functions out to 70 GeV in Au+Au
2. **Upsilon spectroscopy** -- separate 1S, 2S, 3S states via e+e- invariant mass (~100 MeV resolution)
3. **Heavy-flavor tagged jets** -- b-jet R_AA via displaced vertex reconstruction

PPG12 contributes to pillar (1) by measuring the pp baseline cross-section for isolated direct photons at sqrt(s) = 200 GeV, which serves as a calibrated probe of hard-scattering processes unmodified by the QGP.

## Detector Layout

sPHENIX is a compact, cylindrical spectrometer built around the recycled BaBar superconducting solenoid. Moving radially outward from the beam pipe:

| Layer | Subsystem | Radial Extent | Technology |
|-------|-----------|---------------|------------|
| 1 | Beam pipe | ~2 cm | Beryllium (1.8% X0) |
| 2 | MVTX | 2.5--4.0 cm | MAPS silicon (3 layers, 27 um pitch) |
| 3 | INTT | 7.2--10.3 cm | Silicon strip (2 barrel layers, 78 um pitch) |
| 4 | TPC | up to 80 cm | Gas (Ar-CF4-Isobutane), quad-GEM readout |
| 5 | TPOT | on EMCal | Micromegas (8 modules, TPC distortion correction) |
| 6 | EMCal | ~90--115 cm | W/SciFi (20 X0, 0.025 x 0.025 towers) |
| 7 | IHCal | 116--137 cm | Al plates + scintillating tiles (~1 lambda_I) |
| 8 | Solenoid | 140--173 cm | BaBar superconducting coil (1.4 T operational) |
| 9 | OHCal | 182--269 cm | Steel plates + scintillating tiles (~4 lambda_I) |

### Forward and Global Detectors

| Detector | Coverage | Function |
|----------|----------|----------|
| MBD | 3.51 < \|eta\| < 4.61 | Minimum-bias trigger, vertex, centrality |
| sEPD | 2.0 < \|eta\| < 4.9 | Event plane determination |
| ZDC | ~18 m from IP, < 2 mrad | Neutral energy, beam-background rejection |

## Magnetic Field

The superconducting solenoid is the former BaBar magnet from SLAC's PEP-II facility. Key parameters:

- **Design field:** 1.5 T
- **Operational field:** 1.4 T (as used in Run 24)
- **Operating current:** 4596 A
- **Stored energy:** ~27 MJ
- **Cryostat inner radius:** 140 cm; outer radius: 173 cm
- **Field uniformity:** < 5% variation over tracking volume

The cryostat wall represents ~1.4 radiation lengths of material between the inner and outer hadronic calorimeters. The OHCal steel absorber doubles as the magnetic flux return.

## Acceptance

The barrel acceptance is uniform in azimuth and covers:

- **Calorimeters:** |eta| < 1.1, full 2pi
- **Tracking:** |eta| < 1.1 for |z_vtx| < 10 cm
- **PPG12 fiducial region:** |eta| < 0.7 (tighter cut to ensure full shower containment in EMCal)

## Tracking System

The tracking system provides charged-particle momentum measurement and vertex reconstruction:

- **MVTX:** Three MAPS layers (copy of ALICE ITS inner layers) at r = 2.5, 3.2, 4.0 cm. Primary vertex and displaced vertex reconstruction. 27 um pixel pitch.
- **INTT:** Two silicon strip barrel layers at r ~ 7--10 cm. Provides 60 ns timing resolution for single bunch-crossing readout, critical for pileup suppression in pp running.
- **TPC:** 2.1 m long, 80 cm outer radius gas volume with quad-GEM readout. ~4 m^3 active volume. Space-point resolution < 200 um. Provides the primary momentum measurement.
- **TPOT:** Eight Micromegas modules at the TPC outer radius, used for space-charge distortion corrections.

For PPG12, the tracking system is relevant primarily through the INTT pileup timing and the vertex reconstruction that feeds into cluster kinematics (ET = E / cosh(eta), where eta is recalculated from the vertex-corrected tower positions).

## Triggering and DAQ

### Minimum Bias Trigger

The MBD provides the Level-1 minimum-bias trigger by requiring coincidence of hits in both the north and south arms (>=2 PMTs per side). The MBD trigger cross-section in pp at sqrt(s) = 200 GeV is 25.2 +2.3/-1.7 mb, measured via Vernier scans.

### Calorimeter Trigger

The calorimeter electronics provide a high-energy cluster trigger for photon and jet analyses. The calorimeter trigger rejection exceeds 100 for E_gamma > 10 GeV in pp.

### Electronics

All calorimeter subsystems and the MBD share common readout electronics:

- **Digitization:** 14-bit precision at 56.4 MHz (6x RHIC bunch-crossing rate)
- **Waveform recording:** 12-sample waveforms above 2-sigma pedestal noise threshold; 2 summary samples otherwise
- **Zero suppression:** Hardware-level at 2-sigma pedestal from start of Run 2024
- **Streaming readout:** Tracking detectors use continuous streaming mode, enabling unbiased pp data collection

## What Makes sPHENIX Capable of PPG12

Several design features make the isolated photon measurement possible:

1. **EMCal granularity (0.025 x 0.025):** Provides the shower-shape discrimination needed to separate single photons from pi0 decay pairs. The tower size (~2.5 cm) is close to the Moliere radius (~2.3 cm), so photon showers are concentrated in a small 3x3 cluster while pi0 pairs at high pT produce resolvable or partially merged two-cluster signatures.

2. **Hermetic calorimetry (EMCal + IHCal + OHCal):** Enables measurement of the total energy in the isolation cone (R = 0.3). The ~6 lambda_I total depth contains ~97% of hadronic energy, so the isolation sum captures both electromagnetic and hadronic contributions.

3. **EMCal depth (20 X0):** Ensures full electromagnetic shower containment across the PPG12 pT range (8--35 GeV). Back-leakage into the IHCal is small (~0.83 lambda_I hadronic interaction length).

4. **Energy resolution (sigma/E = 2.8% + 15.5%/sqrt(E)):** Sufficient for measuring the steeply falling photon pT spectrum with bin-by-bin bin-migration corrections via Bayesian unfolding.

5. **Full azimuthal coverage:** The 2pi acceptance provides uniform isolation cone coverage without edge effects (except at the eta boundaries, mitigated by the |eta| < 0.7 fiducial cut).

6. **INTT timing for pileup rejection:** The 60 ns single-crossing timing suppresses out-of-time pileup events that would distort cluster kinematics and contaminate the isolation cone.

7. **Streaming readout in pp:** Enables accumulation of 16.6 pb^-1 (1.5 mrad period) or ~50 pb^-1 (full Run 24) without prescale-limited statistics.

## References

- Upgrade proposal: arXiv:1501.06197
- BaBar solenoid: IEEE Trans. Appl. Supercond. 9 (1999) 847
- sPHENIX TDR: BNL internal document (2019)
- First results: arXiv:2410.18031, arXiv:2504.02240, arXiv:2504.02242
- Predictions: arXiv:2305.15491

## Cross-Links

- [EMCal Performance](emcal-performance.md) -- detailed EMCal specifications and beam test results
- [Calorimeter System](calorimeter-system.md) -- full EMCal + HCal system and isolation measurement
- [RHIC Accelerator](rhic-accelerator.md) -- beam parameters and luminosity
- [Tree Making](../../pipeline/01-tree-making.md) -- how detector data flows into the analysis
- [MC Samples](../../reference/mc-samples.md) -- simulated detector response
