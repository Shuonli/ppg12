# Photon-Jet Correlations

Photon-tagged jet measurements exploit the fact that a prompt photon participates in the hard scattering but does not interact with the QCD medium. Its transverse momentum therefore provides a calibrated tag of the initial parton kinematics -- a "standard candle" that dijet measurements cannot supply, because both jets lose energy in the QGP. This article traces the experimental program from PHENIX photon-hadron correlations through STAR photon+jet to the LHC gamma-jet measurements, and explains why the PPG12 pp isolated photon cross section at sqrt(s) = 200 GeV is the essential baseline for the RHIC jet-quenching program.

---

## 1. Physics Motivation

At leading order, the dominant process for direct photon production at RHIC and LHC is the QCD Compton process qg -> q gamma. The photon and the outgoing quark are produced back-to-back in azimuth. In vacuum (pp), the jet momentum approximately equals the photon momentum, modulo NLO radiation and intrinsic k_T. In a quark-gluon plasma, the away-side parton loses energy, and the measurable consequences include:

- **Momentum imbalance** x_Jgamma = pT_jet / pT_gamma shifts downward.
- **Jet yield suppression** R_Jgamma (fraction of photons with an associated jet above threshold) decreases.
- **Fragmentation function modification** D(z_T) with z_T = pT_hadron / pT_gamma is suppressed at high z_T (hard fragments lost) and enhanced at low z_T (energy redistributed to soft particles).
- **Intra-jet broadening** detected by comparing jet yields at different resolution parameters R.

Because the photon is unmodified by the medium, the ratio of Au+Au to pp per-trigger yields (I_AA) or the nuclear modification factor R_AA for the recoil jet directly quantifies the medium-induced modification.

### Key observables

| Observable | Definition | What it measures |
|---|---|---|
| I_AA | Y_AuAu / Y_pp (per-trigger yield ratio) | Jet yield suppression |
| x_Jgamma | pT_jet / pT_gamma | Momentum imbalance |
| R_Jgamma | Fraction of photons with jet above threshold | Jet survival probability |
| D(z_T) | (1/N_trig) dN/dz_T | Fragmentation function |
| xi | ln(1/z_T) | Fragmentation variable (emphasises soft particles) |
| R_{0.2/0.5} | Yield(R=0.2) / Yield(R=0.5) | Angular redistribution of jet energy |
| S_loss | Fractional energy loss | Absolute energy loss |

---

## 2. PHENIX: Photon-Hadron Correlations at RHIC

PHENIX pioneered the gamma-jet program at RHIC using two-particle correlations between direct photon triggers and associated charged hadrons. Without full jet reconstruction capability, PHENIX used z_T = pT_hadron / pT_gamma as a proxy for the fragmentation function.

### 2.1 First observation of away-side suppression (PRC 80, 024908, 2009)

[arXiv: 0903.3399]

Using 5 < pT_gamma < 15 GeV/c triggers and 1-5 GeV/c associated hadrons at sqrt(s_NN) = 200 GeV, PHENIX observed strong suppression of back-to-back correlations in central Au+Au:

- **I_AA = 0.32 +/- 0.12(stat) +/- 0.09(syst)** for 3 < pT_h < 5 GeV/c hadrons opposite 5-15 GeV/c photons in 0-20% centrality.
- In pp, per-trigger yields scale with z_T, confirming that the correlation measures the away-side fragmentation function.
- Suppression decreases from central to peripheral collisions.

**Photon identification**: EMCal shower-shape cuts provided >98% photon candidate purity above 5 GeV. Direct photon signals were extracted via statistical subtraction using the double ratio R_gamma = (gamma/pi0)_measured / (gamma/pi0)_background, since isolation cuts cannot be applied in the high-multiplicity Au+Au environment. Decay photon-hadron correlations were estimated from measured pi0-hadron correlations via MC pair-by-pair weighting.

### 2.2 Medium modification of the fragmentation function (PRL 111, 032301, 2013)

[arXiv: 1212.3323]

With improved statistics (3.9 billion MB events from 2007 + 2.9 billion from 2010), PHENIX provided the first evidence for medium modification of the jet fragmentation function at RHIC:

- **Suppression at high z_T (low xi)**: I_AA < 1 for xi < 0.8, consistent with parton energy loss depleting hard fragments.
- **Enhancement at low z_T (high xi)**: I_AA > 1 for xi > 0.8, with chi2/dof = 17.6/4 for I_AA = 1 (probability < 0.1%).
- The enhancement is predominantly at **large angles**: the ratio I_AA(full)/I_AA(narrow) = 1.9 +/- 0.3(stat) +/- 0.3(syst) for xi > 0.8, indicating energy redistribution to wide-angle soft particles.

Critically, this paper used a **photon isolation cut in pp** (cone R = 0.3, total transverse energy < 10% of photon energy) -- the closest RHIC precedent for the parametric isolation strategy used in PPG12. In Au+Au, statistical subtraction via R_gamma remained necessary.

### 2.3 Definitive photon-hadron measurement with d+Au (PRC 102, 054902, 2020)

[arXiv: 2005.14270]

The final PHENIX photon-hadron result combined 11.2 billion Au+Au events with the first d+Au gamma-hadron measurement:

- **d+Au: I_dA consistent with unity** across all xi bins (chi2/ndf = 4.0-10.0/5 for the three photon pT bins), confirming negligible cold nuclear matter effects on jet fragmentation.
- **Differential photon pT dependence**: The xi value where I_AA exceeds unity shifts to higher xi (lower z_T) with increasing photon pT. This means the enhancement appears at approximately **fixed low hadron pT (~2 GeV/c)**, not at fixed z_T -- suggesting the excess is a medium response effect rather than a simple redistribution along the fragmentation chain.
- **Angular dependence**: The jet core (|Delta_phi - pi| < pi/6) is suppressed similarly for all photon pT bins, while the broad-angle enhancement grows preferentially for lower-pT triggers.
- **CoLBT-hydro model** (jet transport + hydrodynamic medium response) reproduces the data trends, supporting the interpretation that lost jet energy excites the medium and reappears as soft hadrons.

---

## 3. STAR: Photon+Jet with Reconstructed Jets at RHIC

STAR extended the RHIC gamma-jet program to full jet reconstruction using anti-k_T charged-particle jets in the Time Projection Chamber.

### 3.1 Photon-hadron correlations (PLB 760, 689, 2016)

[arXiv: 1604.01117]

Using the BEMC+BSMD transverse shower profile (TSP) for gamma/pi0 discrimination, STAR measured direct photon and pi0-triggered away-side hadron yields at sqrt(s_NN) = 200 GeV with 23 pb^-1 in pp and 2.8 nb^-1 in Au+Au:

- **I_AA suppression by factor 3-5** at high z_T for both gamma_dir and pi0 triggers in 0-12% central Au+Au.
- **Similar suppression for both trigger types** as a function of z_T, meaning the expected differences from colour-factor effects (quark vs gluon jets on the away side) and surface bias are not resolved at RHIC statistics.
- **Trigger-energy independence**: I_AA for gamma_dir triggers is flat across 8-20 GeV/c trigger pT at fixed z_T, indicating that parton energy loss is not strongly sensitive to initial parton energy in this range.
- **Energy redistribution at low pT**: Lost energy reappears at fixed low hadron pT (< 2 GeV/c) rather than at fixed z_T, consistent with the PHENIX finding.

The TSP method for gamma/pi0 discrimination (E_cluster / sum_i(e_i * r_i^1.5) in BSMD) is the most direct predecessor to the sPHENIX EMCal shower-shape BDT approach used in PPG12.

### 3.2 Direct photon+jet PRL (PRL 134, 012301, 2025)

[arXiv: 2309.00156]

The first gamma_dir+jet measurement with reconstructed jets at RHIC:

- Anti-k_T charged jets at R = 0.2 and R = 0.5, with 15 < ET_trig < 20 GeV photon triggers and recoil jets measured down to pT_jet ~ 10 GeV/c.
- **Medium-induced intra-jet broadening**: R_{0.2/0.5} drops from 0.50 +/- 0.06(syst) in pp to 0.26 +/- 0.09(syst) in 0-15% Au+Au at 10 < pT_jet < 15 GeV/c, demonstrating that quenched jets are significantly broader.
- **Larger suppression at smaller R**: I_AA is markedly lower for R = 0.2 than for R = 0.5, quantifying the angular scale of jet energy redistribution.
- Five jet-quenching models compared (Jet-fluid, LBT, CoLBT-hydro, SCET, Hybrid) -- none fully consistent with all measured distributions, indicating that the interplay of energy loss, medium response, and angular redistribution is not yet theoretically understood.

### 3.3 Companion paper with full details (PRC 111, 014911, 2025)

[arXiv: 2309.00145]

The comprehensive companion paper reports:

- Semi-inclusive recoil-jet distributions for multiple trigger kinematics (gamma_dir 9-20 GeV, pi0 9-15 GeV).
- Quark fraction of recoil partons is significantly larger for gamma_dir triggers (from PYTHIA-6), confirming the colour-charge selection power of the photon tag.
- pp baseline well described by PYTHIA-6 STAR tune.
- Detailed mixed-event procedure for Au+Au background subtraction and closure tests with embedded PYTHIA jets.

---

## 4. LHC: Gamma-Jet at Multi-TeV Energies

The LHC gamma-jet program operates at higher pT (> 30 GeV jets, > 50-60 GeV photons) and in a hotter, longer-lived QGP.

### 4.1 CMS: The "golden channel" established (PLB 718, 773, 2013)

[arXiv: 1205.0206]

The first gamma-jet measurement in heavy-ion collisions:

- Isolated photons (pT > 60 GeV, R = 0.4 isolation cone, UE-subtracted E_T^iso < 1 GeV) paired with anti-k_T R = 0.3 jets (pT > 30 GeV).
- **Mean x_Jgamma drops from 0.86 to 0.73** in 0-10% central PbPb at sqrt(s_NN) = 2.76 TeV.
- **R_Jgamma drops from ~0.69 to 0.49**: an additional 20% of recoil jets fall below the 30 GeV threshold due to energy loss.
- **No significant azimuthal broadening**: jets remain back-to-back with the photon even in central collisions (sigma(Delta_phi) unchanged).

This paper launched the paradigm that isolated photons are the ideal calibrated probe of parton energy loss.

### 4.2 CMS: Photon-tagged jet fragmentation (PRL 121, 242301, 2018)

[arXiv: 1801.04895]

First measurement of fragmentation functions for photon-tagged jets in heavy-ion collisions at sqrt(s_NN) = 5.02 TeV:

- Fragmentation functions measured in two variables: xi_jet (using reconstructed jet pT) and xi_gamma_T (using photon pT as proxy for initial parton pT).
- **xi_gamma_T reveals the full medium modification** because the photon pT is unmodified, while xi_jet dilutes the effect due to jet energy loss.
- **Excess of soft particles and depletion of hard particles** in 0-10% central PbPb relative to pp, with crossover around pT_track ~ 3 GeV/c.
- The chi2 p-value for PbPb = pp is 10^-20 for xi_gamma_T in 0-10% centrality.
- At LHC energies, photon+jet is dominated by quark jets, making this a measurement specifically of quark energy loss and shower modification.

### 4.3 ATLAS: Colour-charge-dependent jet quenching (PLB 846, 138154, 2023)

[arXiv: 2303.10090]

ATLAS compared photon-tagged jet suppression to inclusive jet suppression at sqrt(s_NN) = 5.02 TeV:

- **Photon-tagged jets (75-80% quark-initiated at 50 GeV) show R_AA = 0.60-0.75** in 0-10% central PbPb, significantly higher than inclusive jet R_AA (~0.5).
- The difference directly demonstrates sensitivity to the colour charge of the initiating parton: quarks lose less energy than gluons (C_F/C_A = 4/9 in the Casimir scaling limit).
- **Fractional energy loss S_loss** is 0.07-0.10 for photon-tagged jets vs 0.10-0.14 for inclusive jets at low pT.
- At pT > 200 GeV, where quark fractions converge, S_loss values for both samples converge.
- Five theoretical models tested; all qualitatively reproduce the ordering but none fully describes all observables.

### 4.4 ATLAS: Gamma+jet pp cross section at 13 TeV (PLB 780, 578, 2018)

[arXiv: 1801.00112]

The closest LHC analogue to PPG12's methodology:

- Isolated photon+jet cross section measured as functions of photon ET (> 125 GeV), jet pT (> 100 GeV), Delta_phi, invariant mass, and scattering angle.
- **ABCD background subtraction** with signal leakage fractions f_B, f_C, f_D -- the same mathematical structure as PPG12's c_B, c_C, c_D corrections.
- **ET-dependent parametric isolation**: E_T^iso < 4.2e-3 * ET + 10 GeV, directly analogous to PPG12's reco_iso_max = b + s * ET.
- **NLO JETPHOX agrees within ~10%** across six orders of magnitude in cross section.
- Confirms dominance of t-channel quark exchange (qg -> q gamma) via the cos(theta*) distribution.

---

## 5. Synthesis: The Experimental Landscape

The table below summarises the key gamma-jet measurements, ordered by the level of jet reconstruction:

| Experiment | System | sqrt(s_NN) | Method | Photon pT | Key observable |
|---|---|---|---|---|---|
| PHENIX (2009) | pp, AuAu | 200 GeV | gamma-hadron | 5-15 GeV | I_AA = 0.32 |
| PHENIX (2013) | pp, AuAu | 200 GeV | gamma-hadron | 5-9 GeV | D(xi), enhancement at high xi |
| PHENIX (2020) | pp, AuAu, dAu | 200 GeV | gamma-hadron | 5-12 GeV | I_dA = 1, pT-dependent D(xi) |
| STAR (2016) | pp, AuAu | 200 GeV | gamma-hadron | 8-20 GeV | I_AA factor 3-5, trigger-pT independent |
| STAR (2025) | pp, AuAu | 200 GeV | gamma+jet (anti-kT) | 9-20 GeV | R_{0.2/0.5} broadening |
| CMS (2013) | pp, PbPb | 2.76 TeV | gamma+jet | > 60 GeV | x_Jgamma shift, R_Jgamma drop |
| CMS (2018) | pp, PbPb | 5.02 TeV | gamma+jet FF | > 60 GeV | D(xi_gamma_T) modification |
| ATLAS (2023) | pp, PbPb | 5.02 TeV | gamma+jet R_AA | > 50 GeV | Colour-charge dependent S_loss |

### Common findings across experiments and energies

1. **Recoil jets opposite photons are strongly suppressed** in central heavy-ion collisions (factor 3-5 at RHIC, factor ~2 at LHC at accessible pT).
2. **Soft particle excess at wide angles**: Lost energy reappears as soft hadrons at fixed low pT (< 2-3 GeV/c), not at fixed z_T, across both RHIC and LHC.
3. **No significant azimuthal broadening** of the jet axis relative to the photon direction, even in central collisions (CMS 2.76 TeV).
4. **Intra-jet broadening observed at RHIC** via R_{0.2/0.5} (STAR 2025), consistent with energy redistribution from the jet core to periphery.
5. **Cold nuclear matter effects are negligible**: d+Au I_dA ~ 1 (PHENIX 2020), confirming that the modification is a final-state QGP effect.
6. **Theoretical models remain challenged**: No single model simultaneously reproduces the yield suppression, angular redistribution, and fragmentation function modification at both RHIC and LHC energies.

---

## 6. Relevance of PPG12 to the Gamma-Jet Program

PPG12 measures the isolated photon cross section in pp at sqrt(s) = 200 GeV using sPHENIX. This measurement is the essential foundation for the RHIC jet-quenching program for several interconnected reasons:

### 6.1 The pp baseline for R_AA and I_AA

Every gamma-jet nuclear modification measurement -- I_AA at RHIC, x_Jgamma and R_Jgamma at the LHC -- requires a precise pp reference. At RHIC, the pp photon spectrum has been measured by PHENIX via statistical subtraction (R_gamma technique), which yields inclusive direct photon cross sections with 20-40% uncertainties at pT = 5-15 GeV. PPG12 replaces this with an isolated photon cross section measured via shower-shape BDT + parametric isolation + ABCD background subtraction, providing event-by-event photon identification rather than statistical extraction. The isolated cross section is both more precise and more directly comparable to NLO pQCD predictions (JETPHOX), since the isolation criterion suppresses the fragmentation component.

### 6.2 Photon pT calibration for fragmentation functions

The CMS measurement (1801.04895) demonstrated that xi_gamma_T = ln(pT_gamma^2 / pT_track dot pT_gamma) reveals the full magnitude of medium modification because the photon pT is unmodified. This calibration power depends on knowing the photon spectrum precisely in pp. At RHIC, where sPHENIX extends the photon pT reach to 35 GeV (roughly double PHENIX), PPG12 establishes this spectrum.

### 6.3 Complementary quark/gluon composition at RHIC

The ATLAS paper (2303.10090) showed that photon-tagged jets are 75-80% quark-initiated at 50 GeV at the LHC, while inclusive jets are only 30-40% quark-initiated. At RHIC sqrt(s) = 200 GeV, the quark/gluon mixture differs due to different parton distribution function sampling (higher x values probed at lower sqrt(s)). The sPHENIX gamma-jet program, with PPG12 as the photon baseline, will test whether the colour-charge dependence of energy loss (C_F vs C_A) observed at the LHC persists in a different kinematic regime.

### 6.4 Testing theoretical models at different QGP conditions

The QGP created at RHIC is cooler and shorter-lived than at the LHC. The STAR 2025 result showed that no jet-quenching model simultaneously reproduces all observables at sqrt(s) = 200 GeV. Extending the photon pT reach to 35 GeV with sPHENIX provides a wider lever arm for discriminating models. PPG12 validates the photon reconstruction, isolation efficiency, and purity extraction that are prerequisites for the photon-tagged jet analysis.

### 6.5 Methodological continuity

PPG12's analysis methods -- ABCD background subtraction, parametric ET-dependent isolation (R = 0.3), shower-shape BDT -- are the same techniques used at the LHC (ATLAS 1801.00112, 1506.08552; ALICE 2409.12641). By demonstrating these methods at RHIC, PPG12 establishes the bridge for sPHENIX to perform the full gamma-jet program with compatible analysis strategies.

---

## 7. Open Questions and the sPHENIX Outlook

The gamma-jet program has identified several open questions that sPHENIX, with PPG12 as the pp photon baseline, is positioned to address:

1. **Is I_AA truly trigger-pT independent at RHIC?** STAR (1604.01117) found flat I_AA vs trigger pT for 8-20 GeV; sPHENIX extends to 35 GeV photons, testing whether this persists at higher parton energies where radiative vs collisional energy loss contributions differ.

2. **What is the angular scale of jet energy redistribution?** The STAR R_{0.2/0.5} measurement (2309.00156) provides the first RHIC constraint. sPHENIX, with full electromagnetic and hadronic calorimetry plus tracking, can measure calorimetric jets at multiple R values with better angular resolution.

3. **Does the modified fragmentation function "universality" break down?** At RHIC, gamma_dir and pi0 triggers show similar I_AA (1604.01117). With higher statistics and jet reconstruction, sPHENIX can test whether the quark/gluon flavour dependence seen at the LHC (2303.10090) emerges at RHIC.

4. **What is the role of medium response (wake)?** PHENIX (2005.14270) and CMS (1801.04895) observe soft-particle excess consistent with medium excitation. Full jet reconstruction at sPHENIX, calibrated by PPG12 photons, can quantify the wake contribution to jet observables.

5. **Can JETPHOX NLO predictions be validated at sqrt(s) = 200 GeV?** ATLAS (1801.00112) validated JETPHOX at 13 TeV; ATLAS p+Pb (1903.02209) noted a 10-20% NLO deficit at low ET. PPG12 provides the first rigorous test at RHIC energy, where the fragmentation component is relatively larger.

---

## Cross-References

- [ABCD Method](../../concepts/abcd-method.md) -- Background subtraction framework shared with LHC experiments
- [Isolation Cuts](../../concepts/isolation-cuts.md) -- Parametric isolation in PPG12 vs fixed cuts at LHC
- [Shower Shape Variables](../../concepts/shower-shape-variables.md) -- BDT-based photon ID evolved from PHENIX/STAR shower-shape cuts
- [Nuclear Modification](nuclear-modification.md) -- Photon R_AA and R_pA measurements
- [Pipeline Overview](../../pipeline/overview.md) -- PPG12 analysis chain

## Bibliography

- 0903.3399: PHENIX photon-hadron, PRC 80, 024908 (2009)
- 1212.3323: PHENIX fragmentation modification, PRL 111, 032301 (2013)
- 2005.14270: PHENIX definitive gamma-h, PRC 102, 054902 (2020)
- 1604.01117: STAR gamma-hadron, PLB 760, 689 (2016)
- 2309.00156: STAR gamma+jet PRL, PRL 134, 012301 (2025)
- 2309.00145: STAR gamma+jet PRC, PRC 111, 014911 (2025)
- 1205.0206: CMS first gamma-jet quenching, PLB 718, 773 (2013)
- 1801.04895: CMS gamma-jet fragmentation, PRL 121, 242301 (2018)
- 2303.10090: ATLAS photon-tagged jet R_AA, PLB 846, 138154 (2023)
- 1801.00112: ATLAS gamma+jet pp 13 TeV, PLB 780, 578 (2018)
