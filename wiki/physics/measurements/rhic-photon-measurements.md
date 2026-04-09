---
title: "RHIC Direct and Isolated Photon Measurements"
scope: "History of direct/isolated photon measurements at RHIC, from PHENIX's first pp result through thermal photon discovery, the direct photon puzzle, and how PPG12 represents the next generation with sPHENIX"
source_papers: [hep-ex/0502006, nucl-ex/0503003, 0804.4168, 1006.1347, 1105.4126, 1205.5533, 1208.1234, 1405.3940, 2202.08158, 2203.17187]
related_articles:
  - world-data-landscape.md
  - ../../concepts/abcd-method.md
  - ../../concepts/isolation-cuts.md
  - ../../concepts/shower-shape-variables.md
  - ../../pipeline/overview.md
---

# RHIC Direct and Isolated Photon Measurements

The Relativistic Heavy Ion Collider has produced a comprehensive program of direct photon measurements spanning two decades (2001--2024), using the PHENIX detector across pp, d+Au, and Au+Au collision systems at sqrt(s_NN) = 200 GeV, and pp at 510 GeV. These measurements established the NLO pQCD baseline, discovered thermal photon radiation from the QGP, revealed the direct photon flow puzzle, and motivated the next-generation sPHENIX measurement (PPG12). This article traces that evolution chronologically and thematically.

## Timeline of RHIC Direct Photon Measurements

| Year | Collaboration | System | sqrt(s) | Technique | pT range | Key result | Reference |
|------|--------------|--------|---------|-----------|----------|------------|-----------|
| 2005 | PHENIX | Au+Au | 200 GeV | Calorimetric R_gamma | 4--13 GeV/c | R_AA ~ 1 (Ncoll scaling) | [nucl-ex/0503003](../../raw/papers/A_rhic/nucl-ex_0503003.md) |
| 2007 | PHENIX | pp | 200 GeV | Calorimetric R_gamma | 5.5--7 GeV/c | First pp signal, 3 data points | [hep-ex/0502006](../../raw/papers/A_rhic/hep-ex_0502006.md) |
| 2010 | PHENIX | Au+Au, pp | 200 GeV | Internal conversion | 1--5 GeV/c | Thermal photon discovery, T = 221 MeV | [0804.4168](../../raw/papers/A_rhic/0804.4168.md) |
| 2010 | PHENIX | pp | 200 GeV | Isolated + correlations | 5--15 GeV/c | kT measurement, qg Compton dominance | [1006.1347](../../raw/papers/A_rhic/1006.1347.md) |
| 2012 | PHENIX | pp | 200 GeV | Calorimetric pi0-tagging | 5.25--25 GeV/c | Precision cross section, NLO agreement | [1205.5533](../../raw/papers/A_rhic/1205.5533.md) |
| 2012 | PHENIX | Au+Au | 200 GeV | Calorimetric R_gamma | 1--13 GeV/c | Direct photon v2 puzzle | [1105.4126](../../raw/papers/A_rhic/1105.4126.md) |
| 2013 | PHENIX | d+Au, pp | 200 GeV | Three methods combined | 1--16 GeV/c | R_dA ~ 1 (no CNM effects) | [1208.1234](../../raw/papers/A_rhic/1208.1234.md) |
| 2015 | PHENIX | Au+Au | 200 GeV | External conversion (HBD) | 0.4--5 GeV/c | Npart^1.38 scaling, T ~ 240 MeV | [1405.3940](../../raw/papers/A_rhic/1405.3940.md) |
| 2023 | PHENIX | pp | 510 GeV | Calorimetric pi0-tagging | 6--30 GeV/c | ALL resolves gluon spin sign | [2202.08158](../../raw/papers/A_rhic/2202.08158.md) |
| 2024 | PHENIX | Au+Au | 200 GeV | External conversion (VTX) | 0.8--10 GeV/c | pT-dependent slope, dNch scaling | [2203.17187](../../raw/papers/A_rhic/2203.17187.md) |
| 2025+ | sPHENIX | pp | 200 GeV | BDT + ABCD isolation | 8--35 GeV/c | Isolated photon cross section | PPG12 |

## Phase I: Establishing the pQCD Baseline (2001--2007)

### First Au+Au measurement (Run 2, 2001)

The first RHIC direct photon result came from Au+Au collisions, not pp. Using Run 2 data (30 million minimum-bias events), PHENIX measured the direct photon spectrum up to pT ~ 13 GeV/c in nine centrality bins ([nucl-ex/0503003](../../raw/papers/A_rhic/nucl-ex_0503003.md)). The key finding was that direct photon R_AA is consistent with unity at all centralities for pT > 6 GeV/c, while pi0 R_AA drops to ~0.2--0.3 in central collisions. This simultaneously confirmed:

1. High-pT direct photons scale with Ncoll (binary collision scaling), as expected for an initial-state process.
2. The strong pi0 suppression is a final-state effect (jet quenching in the QGP), not an initial-state modification.
3. NLO pQCD (Gordon-Vogelsang, CTEQ6, GRV FF) correctly predicts the direct photon production rate.

The analysis used the double-ratio technique R_gamma = (gamma/pi0)_meas / (gamma/pi0)_bkg, where the background ratio is computed from a hadronic cocktail based on the measured pi0 spectrum with eta/pi0 = 0.45 +/- 0.05 at high pT. This method, with its 10--14% systematic uncertainty on the gamma/pi0 ratio, became the standard PHENIX approach for a decade.

### First pp measurement (Run 2)

The pp baseline was measured with far less statistics: 16.7 million minimum-bias events (L = 0.77 nb^{-1}) plus EMCal-triggered events (L = 40.3 nb^{-1}) ([hep-ex/0502006](../../raw/papers/A_rhic/hep-ex_0502006.md)). The direct photon signal was observed in only a narrow pT window (5.5--7 GeV/c), with three data points at pT = 5.75, 6.25, 6.75 GeV/c. Statistical uncertainties were 42--75% and systematic uncertainties 37--56%. Despite these limitations, the cross section values were consistent with NLO pQCD, validating the theoretical framework.

The methodology established key PHENIX conventions: the BBC trigger cross section (21.8 mb, ~50% of inelastic pp), the bias correction factor f = 0.75, and the energy scale precision of 1.5%. No isolation cut was applied -- the measurement was fully inclusive.

## Phase II: Thermal Photons and the Direct Photon Puzzle (2008--2012)

### Discovery of thermal photon radiation

The landmark measurement of 2010 used a completely different technique: internal conversions of virtual photons to e+e- pairs ([0804.4168](../../raw/papers/A_rhic/0804.4168.md)). For low-mass pairs (m_ee < 0.3 GeV/c^2) at 1 < pT < 5 GeV/c, the yield above the hadronic cocktail was fitted with a two-component model separating the direct photon fraction r from hadronic decay contributions. The key advantage is a 5x improvement in signal-to-background above the pi0 Dalitz cutoff at m_ee ~ 135 MeV, enabling detection of a 10% direct photon signal as a 50% excess in pairs.

In central Au+Au (0--20%), the excess above Ncoll-scaled pp is exponential with inverse slope T = 221 +/- 19(stat) +/- 19(syst) MeV. Hydrodynamical models require initial temperatures T_init = 300--600 MeV (at tau_0 = 0.6--0.15 fm/c) to reproduce the yield -- well above the lattice QCD deconfinement temperature of ~170 MeV. This was the first direct evidence for thermal radiation from a QGP formed at RHIC.

### Precision pp cross section and isolation (Run 5+6)

With 533 million photon-triggered events from Runs 5 and 6 (L ~ 14.5 pb^{-1}), PHENIX published the definitive pp measurement ([1205.5533](../../raw/papers/A_rhic/1205.5533.md)), extending the pT range to 25 GeV/c with 18 data points (from 5.25 GeV/c). The measurement used pi0-tagging statistical subtraction: N_dir = N_incl - (1+A)(1+R) N_pi0, where A = 0.235 accounts for non-pi0 hadronic decays and R corrects for untagged pi0 photons.

A key result was the isolated-to-inclusive ratio: above 10 GeV/c, more than 90% of direct photons pass the isolation cut (cone R = 0.5, E_cone < 10% E_gamma). This confirmed that the fragmentation photon component is small at PHENIX kinematics, where quark-gluon Compton scattering dominates (~80% at pT = 10 GeV/c). The NLO pQCD prediction of this ratio agreed with data.

The companion correlation measurement ([1006.1347](../../raw/papers/A_rhic/1006.1347.md)) used the same dataset to study isolated direct photon-triggered azimuthal correlations. The near-side yield for direct photon triggers was consistent with zero (< 15% of pi0 near-side yield), confirming negligible fragmentation photon contamination in the isolated sample. The away-side x_E distribution slope (b = 8.2 +/- 0.3) matched the quark fragmentation function prediction, validating the qg Compton picture.

### The direct photon v2 puzzle

The measurement of direct photon elliptic flow v2 in Au+Au ([1105.4126](../../raw/papers/A_rhic/1105.4126.md)) created one of the enduring puzzles of RHIC heavy-ion physics. At low pT (1--4 GeV/c), the direct photon v2 is large (~0.1--0.15), comparable in magnitude to pi0 v2. This is paradoxical:

- The **large yield** (exponential excess over Ncoll-scaled pp) implies early emission when the medium temperature is highest.
- The **large v2** requires late emission when collective flow has developed and the spatial anisotropy has been converted to momentum anisotropy.

No theoretical model has simultaneously explained both observations. Models with early thermalization (tau_0 = 0.2--0.4 fm/c) reproduce the yield shape but severely underestimate v2.

Above pT ~ 5 GeV/c, direct photon v2 drops to zero, consistent with dominance of initial hard scattering. This confirms that PPG12's operating range (8--35 GeV/c) is firmly in the prompt-dominated regime where thermal contamination is negligible and NLO pQCD applies.

## Phase III: Cold Nuclear Matter and Extended Kinematics (2013--2015)

### d+Au baseline

The d+Au measurement ([1208.1234](../../raw/papers/A_rhic/1208.1234.md)) combined three independent analysis methods -- virtual photon (1--5 GeV/c), pi0-tagging (3--16 GeV/c), and statistical subtraction (3--6 GeV/c) -- to cover 1--16 GeV/c. The nuclear modification factor R_dA is consistent with unity across the full pT range, ruling out large cold nuclear matter effects (Cronin enhancement, shadowing, isospin, initial-state energy loss) on direct photon production at 200 GeV.

This result is important for PPG12 because it validates that NLO pQCD predictions for pp are directly applicable without nuclear modification corrections, and that any excess observed in Au+Au is a genuine medium effect.

The pp data from this analysis provided an improved empirical fit function: E d^3sigma/dp^3 = a * pT^{-(b + c*ln(xT))} * (1 - xT^2)^n, with effective power-law index b + c*ln(xT) = 4.6--5.5 for 0.01 < xT < 0.1. NLO pQCD shows a preference for scale mu = 0.5 pT.

### Low-pT thermal photons with external conversions

PHENIX pushed the measurement to pT = 0.4 GeV/c using photon conversions in the HBD readout plane ([1405.3940](../../raw/papers/A_rhic/1405.3940.md)). The thermal photon spectrum has an inverse slope of ~240 MeV/c (centrality-independent within uncertainties), and the integrated yield scales as Npart^{1.38 +/- 0.03}, faster than Ncoll ~ Npart^{4/3} scaling but slower than Npart^2 (surface emission).

## Phase IV: Precision Era and New Physics Channels (2023--2024)

### Gluon spin at 510 GeV

The 510 GeV measurement ([2202.08158](../../raw/papers/A_rhic/2202.08158.md)) combined the cross section with the first measurement of the double-helicity asymmetry A_LL for direct photons. With 108 pb^{-1} of polarized pp data, the isolated photon A_LL in the range 0.02 < x_gluon < 0.08 discriminates between positive and negative gluon polarization scenarios at > 2.8 sigma, resolving a fundamental ambiguity that pi0 and jet data could not.

For the cross section, a critical finding was that NLO pQCD underestimates the **inclusive** direct photon cross section by up to 3x at pT < 12 GeV/c (attributed to multiparton interactions), but agrees well with the **isolated** cross section. This strongly motivates the use of isolation as a primary analysis tool -- the approach adopted by PPG12.

The pileup correction procedure (logarithmic extrapolation to zero collision rate, r_pileup ~ 0.8 for inclusive, ~0.9 for isolated) is relevant to PPG12's double-interaction pileup studies at RHIC.

### Definitive Au+Au thermal measurement with VTX conversions

The final PHENIX Au+Au measurement ([2203.17187](../../raw/papers/A_rhic/2203.17187.md)) used the silicon vertex tracker (VTX) as a passive photon converter, achieving 10x the statistics of previous measurements (12.5 billion events from the 2014 run). The continuous coverage from 0.8 to 10 GeV/c and fine centrality binning (nine 10% bins) revealed a new feature: the nonprompt (thermal) photon inverse slope increases with pT, from Teff ~ 0.26 GeV/c at low pT to ~0.38 GeV/c above 2 GeV/c. This suggests that higher-pT nonprompt photons probe earlier, hotter stages of the collision.

The integrated nonprompt yield scales as (dNch/deta)^{1.1}, and the data are consistent with Ncoll-scaled pp for pT > 5 GeV/c, further confirming the prompt-dominated regime above the thermal contribution.

## Methodological Evolution at RHIC

### From statistical subtraction to isolation

PHENIX used statistical subtraction (double-ratio R_gamma with pi0 tagging) for all measurements through 2023. This method requires:

- Precise knowledge of the pi0 spectrum (5--8% systematic on yield)
- Hadronic decay cocktail (eta/pi0, omega/pi0 ratios; mT-scaling assumption)
- Large systematic amplification at low pT (W = N_incl/N_dir >> 1)

The isolation cut was introduced at PHENIX in [1006.1347](../../raw/papers/A_rhic/1006.1347.md) and [1205.5533](../../raw/papers/A_rhic/1205.5533.md) (cone R = 0.5, E_cone < 10% E_gamma), but was applied on top of the statistical subtraction rather than as the primary background rejection tool.

PPG12 fundamentally changes the approach by using isolation as one axis of the ABCD method and BDT shower-shape discrimination as the other. This eliminates the need for pi0 spectral knowledge and the large systematic amplification inherent in R_gamma. See [ABCD Method](../../concepts/abcd-method.md) and [Pipeline Overview](../../pipeline/overview.md).

### Detector evolution: PHENIX to sPHENIX

| Feature | PHENIX EMCal | sPHENIX CEMCal |
|---------|-------------|----------------|
| Technology | PbSc + PbGl | W-SciFi (projective) |
| |eta| coverage | < 0.35 | < 1.1 (barrel) |
| Azimuthal coverage | 2 x pi/2 arms | Full 2pi |
| Tower size | 5.5 x 5.5 cm (PbSc) | 2.4 x 2.4 cm |
| Energy resolution | 8.1%/sqrt(E) + 5.7% | ~15%/sqrt(E) |
| pi0 merging onset | 12--17 GeV/c | ~8 GeV/c |
| Acceptance for R=0.3 cone | Marginal | Full |

The sPHENIX CEMCal's full azimuthal coverage and wider eta acceptance provide roughly 20x the geometric acceptance of PHENIX, while the finer tower granularity (2.4 cm vs 5.5 cm) extends the two-photon resolution for pi0 discrimination. However, the sampling calorimeter resolution (~15%/sqrt(E)) is worse than PHENIX PbSc, making shower-shape BDT discrimination more important. The full 2pi coverage ensures that the isolation cone (R = 0.3) is always fully contained within the detector acceptance for |eta| < 0.7.

## How PPG12 Builds on This Legacy

PPG12's isolated photon cross section measurement at sqrt(s) = 200 GeV with sPHENIX Run 24 data (16.6 pb^{-1}) represents the convergence of two decades of RHIC photon physics:

1. **pT range**: 8--35 GeV/c (12 bins), extending the PHENIX measurement ([1205.5533](../../raw/papers/A_rhic/1205.5533.md)) from 25 to 35 GeV/c and replacing statistical subtraction with the ABCD method.

2. **Background subtraction**: The ABCD method with BDT-based tight/non-tight classification and parametric isolation (reco_iso_max = b + s * ET, cone R = 0.3) replaces the pi0-tagging double-ratio. Signal leakage corrections c_B, c_C, c_D from MC account for ABCD independence violations.

3. **Photon identification**: A 25-feature XGBoost BDT trained on GEANT simulation, applied as ET-dependent thresholds (tight_bdt_min = intercept + slope * ET), replacing PHENIX's simple shower-shape and ToF cuts.

4. **Theory comparison**: JETPHOX NLO pQCD with modern PDFs (CTEQ-TEA) and scale variations, building on the Vogelsang/Gordon NLO framework used throughout the PHENIX era.

5. **Systematic framework**: 23+ systematic variations generated programmatically, run on HTCondor, and aggregated into groups (purity, efficiency, energy scale, energy resolution) via quadrature -- a formalization of the uncertainty assessment developed across all PHENIX measurements.

6. **pp baseline for future sPHENIX Au+Au**: Just as the PHENIX pp measurements ([hep-ex/0502006](../../raw/papers/A_rhic/hep-ex_0502006.md), [1205.5533](../../raw/papers/A_rhic/1205.5533.md)) provided the reference for Au+Au R_AA and thermal photon extraction, PPG12 will serve as the prompt photon baseline for sPHENIX Au+Au data.

## Open Questions That PPG12 Addresses

- **NLO pQCD at moderate xT**: Does the JETPHOX NLO prediction describe the isolated photon cross section at sqrt(s) = 200 GeV in the xT = 0.08--0.35 range, where PHENIX's inclusive measurement showed agreement but isolation-based measurements have not been performed?

- **Low-ET excess**: CDF observed a steeper-than-NLO slope at ET < 50 GeV; PHENIX at 510 GeV saw a factor ~3 excess for inclusive photons at pT < 12 GeV. Does PPG12 at 8--35 GeV/c see any tension between isolated photon data and NLO at the low end of its pT range?

- **Isolation efficiency at RHIC**: With the smaller underlying event at sqrt(s) = 200 GeV compared to LHC energies, the isolation cut efficiency may differ significantly. PPG12's parametric cut (ET-dependent threshold) is designed to handle this, but the isolated-to-inclusive ratio at RHIC energies provides a new data point for understanding fragmentation photon contributions.

- **sPHENIX detector performance**: Can the sPHENIX CEMCal's shower-shape BDT achieve sufficient photon/pi0 discrimination despite the sampling calorimeter resolution? PPG12 is the first physics measurement to validate this capability.

## See Also

- [World Data Landscape](world-data-landscape.md) -- Global context of isolated/direct photon measurements
- [ABCD Method](../../concepts/abcd-method.md) -- PPG12's background subtraction
- [Isolation Cuts](../../concepts/isolation-cuts.md) -- PPG12's parametric isolation definition
- [Shower Shape Variables](../../concepts/shower-shape-variables.md) -- BDT input features
- [Pipeline Overview](../../pipeline/overview.md) -- PPG12 analysis pipeline end-to-end

## Source Papers

- [PHENIX Au+Au 200 GeV calorimetric (2005)](../../raw/papers/A_rhic/nucl-ex_0503003.md) -- First RHIC direct photon, Ncoll scaling
- [PHENIX pp 200 GeV first measurement (2007)](../../raw/papers/A_rhic/hep-ex_0502006.md) -- 3-point pp baseline
- [PHENIX Au+Au thermal photon discovery (2010)](../../raw/papers/A_rhic/0804.4168.md) -- Internal conversion, T = 221 MeV
- [PHENIX pp 200 GeV isolated correlations (2010)](../../raw/papers/A_rhic/1006.1347.md) -- kT, Compton dominance
- [PHENIX Au+Au v2 puzzle (2012)](../../raw/papers/A_rhic/1105.4126.md) -- Large direct photon flow
- [PHENIX pp 200 GeV precision cross section (2012)](../../raw/papers/A_rhic/1205.5533.md) -- Definitive PHENIX pp measurement
- [PHENIX d+Au 200 GeV (2013)](../../raw/papers/A_rhic/1208.1234.md) -- Cold nuclear matter baseline, R_dA ~ 1
- [PHENIX Au+Au low-pT conversions (2015)](../../raw/papers/A_rhic/1405.3940.md) -- Npart^1.38 scaling
- [PHENIX pp 510 GeV ALL (2023)](../../raw/papers/A_rhic/2202.08158.md) -- Gluon spin, NLO for isolated vs inclusive
- [PHENIX Au+Au VTX conversions (2024)](../../raw/papers/A_rhic/2203.17187.md) -- pT-dependent slope, definitive Au+Au result
