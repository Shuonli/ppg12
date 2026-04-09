# Nuclear Modification of Photon Production

Prompt photons are produced in the initial hard scattering and do not interact strongly with the QCD medium. Their nuclear modification factor R_AA is therefore expected to be unity at high pT, in sharp contrast to the factor-of-five suppression observed for hadrons. This makes photons a calibrated probe of the initial state: any deviation from binary (N_coll) scaling signals either initial-state nuclear effects (shadowing, anti-shadowing, isospin, energy loss) or final-state thermal radiation. This article reviews R_AA and R_pA measurements across experiments and energies, from the foundational PHENIX Au+Au result through the LHC program, and explains how the PPG12 pp baseline enables the RHIC nuclear modification program.

---

## 1. Theoretical Framework

### 1.1 Binary scaling and the nuclear modification factor

In the absence of nuclear effects, hard-process yields in nucleus-nucleus collisions scale with the number of binary nucleon-nucleon collisions N_coll (or equivalently with the nuclear thickness function T_AA = N_coll / sigma_pp^inel):

R_AA = (1 / T_AA) * (d^2 N_AA / dpT dy) / (d^2 sigma_pp / dpT dy)

For prompt photons at high pT, R_AA = 1 is the null hypothesis. Deviations arise from:

- **Initial-state effects** (cold nuclear matter, CNM):
  - Nuclear parton distribution functions (nPDFs): shadowing (suppression at low x), anti-shadowing (enhancement at moderate x), EMC effect (suppression at large x)
  - Isospin: neutron content of nuclei modifies the quark flavour composition (u/d ratio)
  - Initial-state parton energy loss: parton radiates energy before the hard scattering
  - Cronin effect: multiple soft scatterings broaden parton k_T, enhancing yield at moderate pT

- **Final-state effects** (hot nuclear matter):
  - Jet quenching: does NOT suppress photons directly (they are colour-neutral), but suppresses the fragmentation photon component
  - Thermal radiation: QGP and hadronic matter emit thermal photons, producing an excess above N_coll-scaled pQCD at low pT (< 3-4 GeV/c)

### 1.2 The photon as an initial-state calibrator

The fact that R_AA^photon ~ 1 while R_AA^pi0 ~ 0.2 in central Au+Au was the key early evidence that hadron suppression is a **final-state effect** (jet quenching in the QGP), not an initial-state effect (modified parton distributions). If the initial-state parton flux were suppressed, both photons and hadrons would be equally affected.

### 1.3 R_pA and cold nuclear matter isolation

Proton-nucleus (p+A) or deuteron-nucleus (d+A) collisions produce no QGP but include all initial-state effects. Measuring R_dA or R_pPb for photons isolates CNM contributions:

R_pA = (1 / (A * sigma_pp)) * d sigma_pA / dpT

If R_pA = 1, CNM effects are negligible. Any R_AA deviation from unity in A+A that is absent in p+A must be a final-state (medium) effect.

---

## 2. RHIC: The Foundational Measurements

### 2.1 N_coll scaling established in Au+Au (PRL 94, 232301, 2005)

[arXiv: nucl-ex/0503003]

PHENIX measured direct photon production in Au+Au at sqrt(s_NN) = 200 GeV up to pT ~ 13 GeV/c in nine centrality bins. This was the first direct test of binary scaling for photons at RHIC.

**Key results:**

- **R_AA consistent with unity** for pT > 6 GeV/c at all centralities, from peripheral (60-80%) to most central (0-10%).
- **pi0 R_AA ~ 0.2-0.3** in central collisions, demonstrating a factor-of-five suppression.
- NLO pQCD (Gordon & Vogelsang, CTEQ6 PDFs) describes the N_coll-scaled yields well.
- No evidence for thermal photon excess above pQCD at pT > 4 GeV/c, in contrast to earlier SPS Pb+Pb results.

**Method:** Statistical subtraction via the double ratio R_gamma = (gamma/pi0)_measured / (gamma/pi0)_background. The photon is not identified event-by-event; rather, the excess of photons above the decay background is extracted statistically. Energy-scale uncertainties partially cancel in this ratio.

This measurement established the paradigm: **photons see through the QGP; hadrons do not**. The NLO pQCD baseline validated here is the ancestor of the JETPHOX predictions that PPG12 compares against.

### 2.2 Thermal photon discovery via internal conversions (PRL 104, 132301, 2010)

[arXiv: 0804.4168]

PHENIX exploited the relation between real photon production and low-mass e+e- pair production (virtual photon internal conversions) to access the direct photon signal at low pT (1-5 GeV/c), where the decay background overwhelms the calorimetric measurement:

- **Enhanced yield in central Au+Au** for pT < 2.5 GeV/c above N_coll-scaled pp.
- The excess is exponential with **inverse slope T = 221 +/- 19(stat) +/- 19(syst) MeV** in 0-20% centrality.
- Hydrodynamical models require **T_init = 300-600 MeV** (at tau_0 = 0.6-0.15 fm/c) to reproduce the data, well above the QCD phase transition temperature T_c ~ 170 MeV.

This was the first direct evidence for thermal QGP radiation at RHIC. The internal conversion technique (fitting the invariant mass spectrum of e+e- pairs for a direct photon component with S -> 1 form factor) achieved a factor-of-five improvement in signal-to-background over calorimetric methods.

### 2.3 Centrality dependence of thermal photons (PRC 91, 064904, 2015)

[arXiv: 1405.3940]

Using photon conversions in the HBD readout plane, PHENIX extended the measurement down to pT = 0.4 GeV/c with fine centrality bins:

- **Inverse slope T_eff ~ 240 MeV/c is centrality-independent** (0-20%: 239 +/- 25, 20-40%: 260 +/- 33, 40-60%: 225 +/- 28, 60-92%: 238 +/- 50 MeV/c for the 0.6-2.0 GeV/c range).
- **Integrated yield scales as N_part^alpha with alpha = 1.38 +/- 0.03(stat) +/- 0.07(syst)**, intermediate between N_part (soft) and N_coll (hard) scaling.
- The centrality-independent shape but centrality-dependent yield constrains the space-time evolution of thermal radiation.

### 2.4 Nonprompt direct photons with 10x statistics (PRC 109, 044912, 2024)

[arXiv: 2203.17187]

The definitive PHENIX measurement used 12.5 billion Au+Au events from 2014 with VTX external conversions:

- Continuous pT coverage from 0.8 to 10 GeV/c in nine 10%-wide centrality bins.
- **New finding: pT-dependent inverse slope.** T_eff increases from ~0.26 GeV/c at low pT (0.8-1.9 GeV/c) to ~0.38 GeV/c at higher pT (2.0-4.0 GeV/c), suggesting that higher-pT nonprompt photons probe earlier, hotter stages of the collision.
- Integrated yield scales as (dN_ch/deta)^alpha with alpha ~ 1.1.
- R_gamma agrees with earlier PHENIX results from HBD conversions, internal conversions, and EMCal, validating consistency across independent techniques.

### 2.5 The direct photon v2 puzzle (PRL 109, 122302, 2012)

[arXiv: 1105.4126]

PHENIX measured the elliptic flow v2 of direct photons in Au+Au:

- **Large v2 ~ 0.1-0.15 at low pT (1-3 GeV/c)**, comparable in magnitude to pi0 flow.
- **v2 consistent with zero above 5 GeV/c**, as expected for initial hard scattering.

The simultaneously large thermal photon yield (requiring early emission when T is high and flow is small) and large v2 (requiring late emission when collective flow is fully developed) creates the **"direct photon puzzle"** -- one of the outstanding problems in heavy-ion physics. Models with early thermalization (tau_0 = 0.2-0.4 fm/c) can reproduce the yield but severely underestimate the v2 magnitude.

This result demonstrates that R_AA and v2 provide complementary constraints on the space-time dynamics of the QGP. The prompt photon component (pT > 5-6 GeV/c, where v2 = 0) is cleanly separated from the thermal component, validating that PPG12's operating range (8-35 GeV/c) is firmly in the prompt-dominated regime.

### 2.6 Cold nuclear matter: d+Au baseline (PRC 87, 054907, 2013)

[arXiv: 1208.1234]

PHENIX measured direct photons in d+Au at sqrt(s_NN) = 200 GeV using three independent methods (virtual photon, pi0-tagging, statistical subtraction) over 1-16 GeV/c:

- **R_dA consistent with unity** across the full pT range.
- Standard CNM calculations (Cronin enhancement, isospin, EKS98 shadowing, initial-state energy loss) are all consistent with data but cannot be individually resolved.
- **Comparison with Au+Au R_AA**: R_AA >> 1 (reaching > 7) below 2 GeV/c in Au+Au, while R_dA ~ 1. This proves unambiguously that the Au+Au excess is a **final-state medium effect** (thermal radiation), not an initial-state nuclear modification.
- Power-law index at high pT: n = 7.17 +/- 0.76 in d+Au, consistent with pp (7.08 +/- 0.09) and Au+Au (7.18 +/- 0.14).

The d+Au result validates the assumption in PPG12 that NLO pQCD predictions for pp are directly applicable at sqrt(s) = 200 GeV without large initial-state nuclear corrections.

---

## 3. LHC: Photon R_AA at Higher Energy

### 3.1 ATLAS: Isolated photons in Pb+Pb at 2.76 TeV (PRC 93, 034914, 2016)

[arXiv: 1506.08552]

The first comprehensive LHC isolated photon measurement in heavy-ion collisions:

- T_AA-scaled yields measured in four centrality classes (0-10%, 10-20%, 20-40%, 40-80%) for 22 < pT < 280 GeV in two pseudorapidity intervals (|eta| < 1.37 and 1.52 < |eta| < 2.37).
- **Yields agree with JETPHOX pp NLO predictions** within uncertainties in all centrality and eta bins -- photon R_AA is consistent with unity.
- Forward-to-central yield ratio R_FC shows slight preference for calculations including isospin effects (neutron content of Pb), though all three scenarios (pp, Pb+Pb isospin, EPS09 nPDF) are consistent with data.
- No evidence for modification of photon production by the hot nuclear medium.

**ABCD method in heavy-ion collisions:** This paper is the definitive demonstration that the double-sideband (tight/non-tight x isolated/non-isolated) method works in heavy-ion environments with large underlying-event fluctuations. The signal leakage factor c_B increases from < 1% in peripheral events to 8-11% in 0-10% central Pb+Pb due to UE contamination of the isolation cone -- directly relevant to PPG12's pileup/double-interaction studies.

**Isolation in heavy ions:** ATLAS uses R = 0.3 (smaller than the R = 0.4 used in pp) with E_T^iso < 6 GeV after UE subtraction. The iterative UE subtraction, accounting for elliptic flow modulation, handles the ~60 GeV of UE energy that enters the isolation cone in central events. PPG12 also uses R = 0.3.

### 3.2 ALICE: Isolated photons in pp and Pb-Pb at 5.02 TeV (EPJC 85, 553, 2025)

[arXiv: 2409.12641]

The most recent comprehensive isolated photon measurement:

- Cross section measured in pp and five Pb-Pb centrality classes (0-10% through 70-90%) at |eta| < 0.67 with both R = 0.2 and R = 0.4 isolation cones.
- **pT reach extends to 10-14 GeV** (lower than ATLAS/CMS at same energy), using charged-particle isolation (pT_iso^ch < 1.5 GeV/c) which is more robust in high-multiplicity than calorimetric isolation.
- **R_AA consistent with unity for 0-70% centrality**, confirming photons are unaffected by the QGP.
- **Peripheral 70-90% centrality: R_AA tends below unity (~0.8-0.9)**, consistent with the HG-PYTHIA centrality bias model (predicted R_AA = 0.82). This is NOT physics suppression but a centrality determination artifact: hard processes bias the centrality classification in peripheral collisions, causing the T_AA denominator to be overestimated. The same effect is observed in CMS Z boson measurements.
- **Cone-radius ratio** (R = 0.4 / R = 0.2) ~ 0.9, consistent with JETPHOX NLO, with no centrality dependence -- first measurement of this ratio at low pT at the LHC.
- Purity is challenging at low pT: 20-30% in Pb-Pb at 10-12 GeV/c (due to pi0 contamination), rising to 50-60% at 18-40 GeV/c. The low-pT challenges are directly relevant to PPG12's lowest bins (8-10 GeV).

### 3.3 ATLAS: Prompt photons in p+Pb at 8.16 TeV (PLB 796, 230, 2019)

[arXiv: 1903.02209]

The first prompt photon measurement in p+Pb:

- Differential cross section measured from 20 GeV to 550 GeV in three nucleon-nucleon centre-of-mass pseudorapidity regions.
- **R_pPb consistent with unity** in all pseudorapidity regions within uncertainties.
- At forward rapidity (probing gluon x ~ 10^-2): nPDF predictions show slight shadowing suppression at low ET and anti-shadowing enhancement at high ET, both within data uncertainties.
- At backward rapidity (probing quark x ~ 0.2): isospin effect causes small suppression; nPDF effects straddle the anti-shadowing to EMC crossover.
- **Initial-state energy loss model with lambda_q = 1.5 fm disfavoured**; data prefer no or limited parton energy loss before the hard scattering.
- JETPHOX NLO systematically underestimates the cross section by 10-20% at low ET (a known NLO-vs-NNLO discrepancy), but this largely cancels in R_pPb.

---

## 4. The Experimental Picture: R_AA and R_pA Summary

| Experiment | System | sqrt(s_NN) | pT range | R_AA or R_pA | Key finding |
|---|---|---|---|---|---|
| PHENIX (2005) | Au+Au | 200 GeV | 6-13 GeV | R_AA ~ 1 | N_coll scaling; pi0 suppressed |
| PHENIX (2010) | Au+Au | 200 GeV | 1-5 GeV | R_AA >> 1 (low pT) | Thermal photon T = 221 MeV |
| PHENIX (2015) | Au+Au | 200 GeV | 0.4-5 GeV | R_AA >> 1 (low pT) | N_part^1.38 scaling, T ~ 240 MeV |
| PHENIX (2024) | Au+Au | 200 GeV | 0.8-10 GeV | R_AA >> 1 (low pT) | pT-dependent T_eff |
| PHENIX (2012) | Au+Au | 200 GeV | 1-13 GeV | v2 ~ 0.1-0.15 (low pT) | Photon flow puzzle |
| PHENIX (2013) | d+Au | 200 GeV | 1-16 GeV | R_dA ~ 1 | CNM negligible |
| ATLAS (2016) | Pb+Pb | 2.76 TeV | 22-280 GeV | R_AA ~ 1 | T_AA scaling confirmed |
| ALICE (2025) | Pb-Pb | 5.02 TeV | 10-140 GeV | R_AA ~ 1 (0-70%) | Peripheral bias in 70-90% |
| ATLAS (2019) | p+Pb | 8.16 TeV | 20-550 GeV | R_pPb ~ 1 | Energy loss disfavoured |

### The consistent picture

1. **High-pT direct/isolated photons scale with N_coll** across all centralities and all collision energies probed so far. There is no evidence that the QGP suppresses prompt photon production.

2. **Low-pT thermal photon excess** (pT < 3-4 GeV/c) is observed exclusively in nucleus-nucleus collisions (not in d+Au or p+Pb), confirming it as a final-state medium effect. The excess is exponential with T_eff ~ 220-240 MeV at RHIC, with a newly discovered pT-dependent slope suggesting sensitivity to the early hot stages.

3. **Cold nuclear matter effects are small** for photons at both RHIC (R_dA ~ 1) and LHC (R_pPb ~ 1) energies. Isospin, shadowing, and initial-state energy loss are all within current experimental uncertainties.

4. **The direct photon v2 puzzle remains unresolved**: large yield (early emission) + large v2 (late emission) cannot be simultaneously explained by standard hydrodynamic models. Proposed solutions include pre-equilibrium photon production, magnetic field effects, and viscous corrections, but none fully resolve the tension.

5. **Centrality determination biases** in peripheral heavy-ion collisions (observed by ALICE in 70-90% Pb-Pb) produce artificial R_AA suppression unrelated to physics. This is important context for interpreting all peripheral-centrality nuclear modification measurements.

---

## 5. Relevance of PPG12

### 5.1 The pp reference at sqrt(s) = 200 GeV

Every nuclear modification measurement at RHIC requires a pp reference spectrum. The PHENIX pp direct photon cross section was measured via statistical subtraction with 20-40% total uncertainties in the 5-15 GeV/c range. PPG12 provides the first sPHENIX isolated photon cross section at sqrt(s) = 200 GeV in the 8-35 GeV range using modern techniques (shower-shape BDT, parametric isolation, ABCD background subtraction). This directly enables:

- **R_AA for sPHENIX Au+Au**: The PPG12 cross section is the denominator in R_AA = (yield_AuAu / T_AA) / sigma_pp. The ATLAS result (1506.08552) demonstrated that R_AA precision is dominated by the pp reference uncertainty.
- **Consistency check with PHENIX**: PPG12 measures the isolated photon cross section, while PHENIX measured the inclusive direct photon cross section. The two are related by the isolation efficiency for fragmentation photons, providing a non-trivial cross-check.

### 5.2 Establishing the prompt photon regime

The PHENIX v2 measurement (1105.4126) showed that direct photon v2 = 0 above 5 GeV/c, confirming that prompt hard-scattered photons dominate. PPG12 operates at 8-35 GeV, firmly within this regime. By validating the isolated photon cross section against NLO JETPHOX at these pT values, PPG12 confirms that the prompt photon production mechanism is well-understood at RHIC energy, establishing the clean baseline needed for thermal photon and jet-quenching measurements.

### 5.3 Constraining NLO pQCD at RHIC energy

JETPHOX NLO has been validated at LHC energies (1801.00112, ATLAS pp 13 TeV) and found to underpredict by 10-20% at low ET in p+Pb (1903.02209). At sqrt(s) = 200 GeV, where the fragmentation component is relatively larger (35-50% at pT ~ 5-10 GeV vs 15-35% at the LHC), the NLO prediction is more sensitive to fragmentation function uncertainties. PPG12 provides the first rigorous test of JETPHOX at RHIC energy with an isolated photon observable.

### 5.4 Cold nuclear matter baseline at RHIC

The PHENIX d+Au measurement (1208.1234) confirmed R_dA ~ 1 with ~20-30% uncertainties. PPG12's precise pp cross section, combined with future sPHENIX d+Au or p+Au data, will tighten the CNM constraints at RHIC kinematics (x ~ 0.04-0.2 for 8-35 GeV photons at sqrt(s) = 200 GeV). The ATLAS p+Pb result (1903.02209) constrains nPDFs at higher energy; RHIC provides complementary coverage at higher x values.

### 5.5 Thermal photon programme context

Although PPG12 operates well above the thermal photon regime (pT >> 4 GeV/c), understanding the prompt cross section precisely is essential for extracting the nonprompt component in Au+Au. The nonprompt yield is obtained by subtracting N_coll-scaled pp from the Au+Au measurement:

gamma_nonprompt = gamma_AuAu - T_AA * sigma_pp^prompt

Any uncertainty in sigma_pp propagates directly into the thermal photon extraction. The PHENIX 2024 result (2203.17187) with pT-dependent inverse slope will benefit from a more precise pp baseline at the high-pT end of the thermal region.

---

## 6. Open Questions

1. **What produces the large thermal photon v2?** Pre-equilibrium radiation, late hadronic sources, magnetic field effects, and viscous corrections have all been proposed. Higher-statistics measurements at both RHIC and LHC are needed to resolve this.

2. **Is the pT-dependent inverse slope a signature of early-time dynamics?** The PHENIX 2024 observation (2203.17187) that T_eff increases from ~260 to ~380 MeV with pT needs confirmation and theoretical interpretation. Measurements at different sqrt(s) can distinguish between scenarios.

3. **How large are nPDF effects at RHIC kinematics?** Current RHIC data (R_dA ~ 1) lack the precision to discriminate between shadowing models. sPHENIX with improved pp and p+A baselines can constrain nPDFs in the intermediate-x regime.

4. **Can the centrality bias in peripheral collisions be controlled at RHIC?** The ALICE observation (2409.12641) of R_AA < 1 in 70-90% Pb-Pb due to centrality determination artifacts is relevant for any sPHENIX peripheral measurement.

5. **Is there a smooth transition from thermal to prompt photons?** The PPG12 pT range (8-35 GeV) and the thermal photon range (0.4-6 GeV) overlap only marginally. Understanding the 4-10 GeV transition region, where both prompt and thermal contribute, requires combining both measurements with controlled uncertainties.

---

## Cross-References

- [ABCD Method](../../concepts/abcd-method.md) -- Background subtraction used at ATLAS, ALICE, and PPG12
- [Isolation Cuts](../../concepts/isolation-cuts.md) -- Isolation definitions across experiments
- [Shower Shape Variables](../../concepts/shower-shape-variables.md) -- Photon ID from PHENIX shower shapes to sPHENIX BDT
- [Photon-Jet Correlations](photon-jet-correlations.md) -- Complementary gamma-jet programme
- [Pipeline Overview](../../pipeline/overview.md) -- PPG12 analysis chain

## Bibliography

### RHIC
- nucl-ex/0503003: PHENIX direct photon Au+Au, PRL 94, 232301 (2005)
- 0804.4168: PHENIX thermal photon excess, PRL 104, 132301 (2010)
- 1105.4126: PHENIX direct photon v2, PRL 109, 122302 (2012)
- 1208.1234: PHENIX direct photon d+Au, PRC 87, 054907 (2013)
- 1405.3940: PHENIX low-pT thermal photon centrality, PRC 91, 064904 (2015)
- 2203.17187: PHENIX nonprompt photon Au+Au, PRC 109, 044912 (2024)

### LHC
- 1506.08552: ATLAS isolated photon Pb+Pb 2.76 TeV, PRC 93, 034914 (2016)
- 2409.12641: ALICE isolated photon pp+PbPb 5.02 TeV, EPJC 85, 553 (2025)
- 1903.02209: ATLAS prompt photon p+Pb 8.16 TeV, PLB 796, 230 (2019)
