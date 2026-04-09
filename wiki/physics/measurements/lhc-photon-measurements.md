# LHC Isolated Photon Measurements

The LHC isolated photon program spans three experiments (ATLAS, CMS, ALICE), four collision energies (2.76--13 TeV), and three collision systems (pp, p+Pb, Pb+Pb). Over 2011--2025, these measurements established the methodology that PPG12 adapts -- ET-dependent parametric isolation, the ABCD signal-leakage equation, BDT-based photon identification, and the progression from NLO JETPHOX comparisons to NNLO NNLOJET. This article traces the experimental evolution with emphasis on techniques and quantitative results most relevant to PPG12.

## ATLAS: The ABCD Method and Parametric Isolation

### First 7 TeV measurements (2011--2014)

ATLAS published two measurements at 7 TeV with dramatically different datasets:

**35 pb^-1 (1108.0253):** The foundational paper. Measured the cross section for 45 < ET < 400 GeV in four eta bins (|eta| < 0.6, 0.6--1.37, 1.52--1.81, 1.81--2.37). Used a **fixed** isolation cut of ET^iso < 3 GeV in a cone R = 0.4, with a non-isolated sideband at ET^iso > 5 GeV. Shower-shape identification split candidates into tight and non-tight categories using lateral and longitudinal calorimeter variables. The 2D sideband (ABCD) method -- tight/non-tight vs isolated/non-isolated -- was established here as the standard ATLAS approach. Signal leakage into control regions reached 17% (tight-iso leaking to non-iso) and 6% elsewhere. Purity ranged from 91% at 45 GeV to near 100% at 200 GeV. NLO JETPHOX (CTEQ 6.6) agreed over four orders of magnitude, with scale uncertainty of ~10% and PDF uncertainty of ~5%.

**4.6 fb^-1 (1311.1440):** Extended the reach to 100--1000 GeV, adopting a **looser** isolation of ET^iso < 7 GeV to retain efficiency at high ET. This paper introduced the explicit ABCD signal-leakage equation that PPG12 directly implements:

```
N_A^sig = N_A - R_BKG * (N_B - c_B * N_A^sig)(N_C - c_C * N_A^sig) / (N_D - c_D * N_A^sig)
```

where c_K = N_K^sig / N_A^sig are signal leakage fractions from MC. The non-tight definition inverted four first-layer EM shower-shape variables (w_s3, F_side, Delta_E_s, E_ratio) -- chosen specifically because they are least correlated with isolation. Fragmentation photons contributed ~20% at low ET, diminishing to negligible above 500 GeV. This observation is important for PPG12, where JETPHOX comparisons must account for fragmentation at lower sqrt(s).

### The 8 TeV breakthrough: ET-dependent isolation (2016)

The 8 TeV measurement (1605.03495, 20.2 fb^-1) introduced the single most important methodological advance for PPG12: **ET-dependent parametric isolation**:

```
ET^iso < 4.8 GeV + 4.2 x 10^-3 x ET^gamma
```

with a non-isolated boundary at ET^iso > 7.8 + 0.0042 x ET (a 3 GeV gap to limit signal leakage). This replaced the fixed thresholds of earlier measurements and optimized the trade-off between signal retention and background rejection across the full ET range of 25--1500 GeV. PPG12 adopts the identical functional form (reco_iso_max = b + s x ET) with sPHENIX-optimized parameters.

The measurement spanned ten orders of magnitude in cross section. Purity rose from ~60% at 25 GeV to ~100% at 300 GeV -- the low-ET behavior foreshadows the challenge PPG12 faces at 8--35 GeV. The R_bkg correlation parameter was validated using subdivision of the B and D control regions into B', B'', D', D'' with shifted ET^iso boundaries; a 10% deviation from unity was observed and propagated as a systematic.

A comparison with the PeTeR resummation code (NLO + NNNLL threshold logarithms, approximately NNLO) was performed for the first time, alongside MCFM. The data showed a ~20% excess over JETPHOX NLO at low ET, decreasing at high ET. This persistent low-ET NLO deficit is a recurring theme across all experiments and energies.

### 13 TeV: first measurement and NNLO era (2017--2019)

**3.2 fb^-1 (1701.06882):** The first 13 TeV result covered ET > 125 GeV. Cross sections were ~2x larger than 8 TeV at low ET and ~10x at high ET, as expected. The same parametric isolation (4.8 + 0.0042 x ET) was retained. The R_bg validation procedure matured: subdivision into B'/B''/D'/D'' with progressive ET^iso shifts of +/-1 GeV, revealing deviations from unity of 15--40% depending on ET and eta. These were propagated as systematic uncertainties. Pythia required 10--14% upward normalization and Sherpa ~30%, illustrating the limitations of LO generators for absolute rates.

**36.1 fb^-1 (1908.02746):** The definitive ATLAS paper. ET range 125--2500 GeV spanning seven orders of magnitude. Three theory comparisons were presented:

| Program | Order | Scale uncertainty | Description |
|---------|-------|-------------------|-------------|
| JETPHOX 1.3.1_2 | NLO | 10--15% | Direct + fragmentation, fixed-cone isolation |
| Sherpa 2.2.2 (ME+PS@NLO) | NLO+PS | Improved normalization | gamma+jet at NLO, gamma+2j/3j at LO, parton shower |
| NNLOJET | NNLO | **0.6--5%** | Direct only, hybrid-cone isolation |

The NNLO comparison was transformative: scale uncertainties dropped by a factor of 2--20 compared to NLO, from 10--15% down to 0.6--5%. NNLOJET provided an excellent description of the data. JETPHOX with mu = ET tended to underestimate data by 5--15%; using mu = ET/2 improved the normalization. The 76-parameter photon energy scale model -- propagated separately to preserve correlations -- represented the state of the art in systematic treatment.

For PPG12, the key lesson is that NLO predictions carry a 10--15% scale uncertainty that will dominate the theory comparison, and that NNLO corrections generally improve both the central value and the uncertainty band.

### Nuclear collisions: Pb+Pb and p+Pb

**Pb+Pb at 2.76 TeV (1506.08552):** Measured for 22 < pT < 280 GeV in four centrality classes. T_AA-scaled yields agreed with JETPHOX pp NLO, confirming binary collision scaling -- photons are unmodified by the QGP. Used R = 0.3 isolation (smaller than pp to reduce UE sensitivity) with ET^iso < 6 GeV. Signal leakage c_B increased from <1% (peripheral) to 8--11% (0--10% central) due to UE contamination of the isolation cone. This directly informs PPG12's pileup studies where additional collision vertices similarly contaminate isolation.

**p+Pb at 8.16 TeV (1903.02209):** Measured R_pPb from 20--550 GeV in three eta* regions. R_pPb was consistent with unity, indicating modest nuclear effects. Data disfavoured large initial-state energy loss (lambda_q = 1.5 fm model rejected). Purity ranged from ~45% at 20 GeV to 99% at 300 GeV. The Bjorken-x range probed (3 x 10^-3 to 0.4) overlaps with RHIC kinematics at sqrt(s) = 200 GeV.

## CMS: Template Fitting and BDT Innovation

### 7 TeV: two complementary methods (2011)

The first CMS measurement (1108.2044, 36 pb^-1) demonstrated two independent purity extraction approaches:

1. **Photon conversion method:** Fitted the ET/pT ratio (calorimeter energy / tracker momentum) for converted photons, where the signal peaks at ET/pT ~ 1 and background extends above 1. Limited to ET < 200 GeV by conversion statistics.

2. **Isolation template method:** Fitted the combined isolation variable (ISO = Iso_TRK + Iso_ECAL + Iso_HCAL) with a signal parametrization (exponential x Gaussian convolution) and background parametrization (threshold function). Signal yields quoted for ISO < 5 GeV.

Both methods yielded consistent results, combined via a weighted average. CMS used sigma_eta_eta as the primary shower-shape variable (< 0.010 in barrel for tight selection). NLO JETPHOX (CT10) agreed with data across 25--400 GeV. Purity rose from 50--60% at 25 GeV to 95--100% above 100 GeV. The isolation cone was R = 0.4 with a generator-level threshold of 5 GeV.

### 13 TeV: BDT-based photon identification (2019)

The CMS 13 TeV measurement (1807.00782, 2.26 fb^-1) introduced **BDT-based photon identification** to cross-section measurements -- the most directly analogous methodology to PPG12:

- TMVA BDT trained separately per rapidity and ET bin
- Input variables: eta, phi, energy, 3x3/total energy ratio, E22/E55, sigma_eta_eta, covariance sigma_eta_phi, energy-weighted spreads sigma_eta and sigma_phi, plus preshower variables in endcaps and median energy density rho
- Signal template from PYTHIA 8, validated with Z -> mu+mu-gamma and Z -> e+e- data
- Background template from data sideband with relaxed hadron isolation (7--13 GeV barrel, 6--12 GeV endcaps)
- Two-template binned likelihood fit to BDT output distribution

This mirrors PPG12's approach of using a trained BDT (XGBoost, 25 features) as the shower-shape discriminator, with signal templates from simulation and background templates from control regions. The sideband approach validates that isolation-relaxed data provides an accurate background shape. CMS measured for ET > 190 GeV, well above PPG12's range, but the methodology is scale-invariant.

### PbPb + pp at 5.02 TeV: sigma_eta_eta template fitting (2020)

CMS measured isolated photons in both PbPb and pp at 5.02 TeV (2003.12797), finding R_AA consistent with unity across all ET and centrality bins. The purity extraction used a two-component template fit of the sigma_eta_eta distribution: signal from MC, background from a nonisolated data sideband (1 < I < 5 GeV). Example purity: 64.1% for 40 < ET < 50 GeV in 10--30% centrality PbPb. Pileup contributed 0--11% systematic uncertainty in pp.

## ALICE: Sampling Calorimetry and Low-pT Reach

ALICE provides the closest detector and methodology analog to PPG12/sPHENIX: the ALICE EMCal is a **lead-scintillator sampling calorimeter** (like sPHENIX CEMCal) with cell granularity Delta_eta x Delta_phi ~ 0.014 x 0.014 (vs ~0.024 for sPHENIX). It shares analogous challenges: APD neutron anomalous signals, T-Card electronics cross-talk, and limited intrinsic resolution compared to crystal calorimeters.

### 7 TeV: first ALICE measurement (2019)

The ALICE 7 TeV measurement (1906.01371, 473 nb^-1) extended the isolated photon pT reach down to 10 GeV/c at |eta| < 0.27 -- below the ATLAS/CMS lower limit of ~25 GeV. The shower-shape discriminant was **sigma_long^2**, the larger eigenvalue of the energy distribution ellipse in the eta-phi plane:

```
sigma_long^2 = (sigma_phi_phi + sigma_eta_eta)/2
    + sqrt((sigma_phi_phi - sigma_eta_eta)^2/4 + sigma_eta_phi^4)
```

Signal (single photons) peaks at sigma_long^2 ~ 0.25; merged pi0 decays create a distinct pT-dependent band at higher values. The narrow cluster selection (0.1 < sigma_long^2 < sigma_max, with sigma_max pT-dependent: 0.4 for 10--14 GeV, 0.35 for 14--16 GeV, 0.3 above 16 GeV) is the direct analog of PPG12's ET-dependent BDT threshold.

Purity extraction used the **ABCD method** with (sigma_long^2, isolation pT) as the two axes -- the same 2D classification PPG12 uses with (BDT score, isolation ET). An MC correction factor alpha_MC accounted for the correlation between sigma_long^2 and isolation for background clusters. Purity was ~20% at 10 GeV/c, saturating at ~50--60% above 18 GeV/c due to merged pi0 contamination. This low-pT purity challenge closely parallels PPG12.

Isolation used both charged tracks and neutral EMCal clusters in R = 0.4, with piso_T < 2 GeV/c. The restricted eta acceptance (|eta| < 0.27) was necessary to contain the full isolation cone within both EMCal and tracker. Energy resolution was sigma_E/E = 1.7% + 11.3%/sqrt(E) + 4.8%/E. The dominant systematic was isolation probability at 20% (low pT), reflecting the same signal-leakage physics as PPG12's purity systematic.

### 13 TeV: lowest xT measurement (2024)

The ALICE 13 TeV measurement (2407.01165, 10.79 pb^-1) achieved the lowest xT (= 2pT/sqrt(s)) of any isolated photon measurement: xT down to 1.08 x 10^-3, with pT coverage 7--200 GeV/c at |eta| < 0.67. Key methodological improvements over 7 TeV:

- **Tracks-only isolation:** Switched from tracks+clusters to charged-track isolation (piso_ch < 1.5 GeV/c in R = 0.4), dramatically expanding eta acceptance from |eta| < 0.27 to |eta| < 0.67 at the cost of missing neutral isolation energy.
- **Simplified sigma_long^2 cut:** 0.1 < sigma_long^2 < 0.3 (flat), replacing the pT-dependent selection.
- **Sigmoid purity fit:** Purity fitted with modified sigmoid functions to smooth statistical fluctuations:
  - Low pT (7--30 GeV/c): parameters a1 = 0.617, b1 = 0.292, c1 = 14.28
  - High pT (30--200 GeV/c): parameters a2 = 0.852, b2 = 0.034, c2 = 2.4
- **Multi-trigger combination:** MB, L1-gamma-low (~5.5 GeV), L1-gamma-high (~12 GeV), with rejection factors RF ~ 471 and 12.65 respectively.

Purity was ~5% at 7 GeV/c, ~50% at 18 GeV/c, rising to ~80% above 80 GeV/c. Total efficiency was 30--45%. NLO JETPHOX (NNPDF4.0) agreed with data. The xT scaling compilation (ISR to LHC) with n ~ 4.5 provides a universal reference for PPG12 at sqrt(s) = 200 GeV.

### 5.02 TeV: pp and Pb-Pb with cone-radius ratio (2025)

ALICE measured at 5.02 TeV (2409.12641) in both pp and five Pb-Pb centrality classes (0--10% through 70--90%), extending to pT = 10--14 GeV/c. Two isolation cone radii (R = 0.2 and R = 0.4) were used simultaneously for the first time, yielding a cone-radius cross-section ratio f_R(0.4/0.2) ~ 0.9, consistent with JETPHOX NLO. R_AA was consistent with unity for 0--70% centrality, but in 70--90% peripheral events R_AA tended below unity (~0.8--0.9), explained by the HG-PYTHIA centrality bias model rather than physics suppression.

### p+Pb: shadowing hints (2025)

The most recent ALICE result (2502.18054) measured in pp at 5.02 and 8 TeV and p-Pb at 5.02 and 8.16 TeV. For pT > 20 GeV/c, R_pA was consistent with unity at both energies. Below 20 GeV/c, a 2.3 sigma hint of suppression appeared at 8.16 TeV, with R_pA < 1 at 1.8 sigma significance, consistent with nuclear shadowing from nNNPDF3.0 and nCTEQ15HQ. This probes Bjorken-x down to ~2.9 x 10^-3. The tight ALICE isolation (piso_ch < 1.5 GeV/c) was found to reduce fragmentation contribution compared to the looser ATLAS parametric isolation, improving NLO description.

## Cross-Experiment Comparison of Methods

| Feature | ATLAS | CMS | ALICE | PPG12 (sPHENIX) |
|---------|-------|-----|-------|-----------------|
| Calorimeter | LAr/Pb, 3 layers | PbWO4 crystals | Pb-scintillator | Pb-scintillator |
| Granularity (eta x phi) | 0.003 x 0.098 (strip) | 0.017 x 0.017 | 0.014 x 0.014 | 0.024 x 0.024 |
| Shower-shape variable | 9 discrete variables | sigma_eta_eta (BDT at 13 TeV) | sigma_long^2 | BDT (25 features) |
| Purity extraction | ABCD with signal leakage | Template fit / BDT fit | ABCD with alpha_MC | ABCD with c_B, c_C, c_D |
| Isolation cone | R = 0.4 (0.3 in Pb+Pb) | R = 0.4 (0.3 for charged) | R = 0.4 (0.2 also at 5 TeV) | R = 0.3 |
| Isolation type | Calorimetric, ET-dependent | Calorimetric + track-based | Charged tracks only (13 TeV) | Calorimetric, ET-dependent |
| Isolation cut | 4.8 + 0.0042 x ET | fixed (varies by analysis) | 1.5 GeV/c (charged) | b + s x ET (parametric) |
| NLO program | JETPHOX, Sherpa, NNLOJET | JETPHOX | JETPHOX | JETPHOX |
| Lowest pT | 25 GeV (8 TeV) | 25 GeV (7 TeV) | 7 GeV/c (13 TeV) | 8 GeV |

## Evolution of Theory Comparisons

The LHC photon program drove a clear progression in QCD precision:

1. **NLO only (2011--2016):** JETPHOX with CTEQ6.6/CT10/CT14 PDFs. Scale uncertainties 10--15%. Data consistently 10--20% above NLO at low ET, converging at high ET. The choice of mu = ET vs mu = ET/2 shifted predictions by ~15%.

2. **NLO + resummation (2016):** PeTeR (NLO + NNNLL threshold logarithms) provided approximately NNLO normalization, improving agreement with data normalization.

3. **NNLO (2019):** NNLOJET provided the first full NNLO calculation for direct photon production. Scale uncertainties dropped to 0.6--5%, a factor 2--20 improvement over NLO. The central value moved closer to data. This calculation uses hybrid-cone isolation (Frixione at small Delta R = 0.1, fixed cone at Delta R = 0.4) and includes only the direct component.

4. **Modern PDF era (2020s):** NNPDF4.0, nNNPDF3.0, CT18 replace earlier sets. PDF uncertainties reduced to 1--6% at LHC kinematics. Data constrains the gluon PDF.

For PPG12 at sqrt(s) = 200 GeV, the relevant comparison remains NLO JETPHOX, where 10--15% scale uncertainty dominates. The LHC NNLO results demonstrate the magnitude of missing higher-order corrections and motivate eventual NNLO comparisons at RHIC kinematics.

## Key Lessons for PPG12

1. **Parametric isolation works.** The ATLAS E_T^iso < 4.8 + 0.0042 x ET functional form, retained across 8 TeV and both 13 TeV papers, is the direct ancestor of PPG12's reco_iso_max = b + s x ET. The gap between isolated and non-isolated boundaries (3 GeV at ATLAS, 0.8 GeV shift at PPG12) limits signal leakage.

2. **Low-pT purity is the dominant challenge.** ATLAS purity of 60% at 25 GeV, ALICE purity of 5--20% at 7--10 GeV/c, and CMS purity of 64% at 40--50 GeV in PbPb all demonstrate that the purity systematic dominates at low ET -- the regime where PPG12 operates.

3. **R_bg validation requires control regions.** The ATLAS B'/B''/D'/D'' procedure with ET^iso boundary shifts (revealing 15--40% deviations from R_bg = 1) provides a concrete recipe for PPG12 to validate its ABCD correlation assumption.

4. **Energy scale dominates at high ET.** Across all experiments, the photon energy scale is the dominant systematic at high ET (5--16% on the cross section from a 1--2% absolute scale uncertainty), while purity dominates at low ET.

5. **NLO undershoots data at low ET.** The persistent ~20% deficit of NLO JETPHOX at low ET across Tevatron and LHC energies implies PPG12 should expect a similar or larger discrepancy at sqrt(s) = 200 GeV.

6. **Sampling calorimeters can do this measurement.** ALICE demonstrates that lead-scintillator EMCal with ~0.014 granularity achieves useful photon identification via sigma_long^2 and ABCD purity extraction down to 7 GeV/c, validating the sPHENIX CEMCal approach.

## Source Papers

| arXiv | Experiment | System | Energy | Key Contribution |
|-------|------------|--------|--------|-----------------|
| 1108.0253 | ATLAS | pp | 7 TeV | Established ABCD method, R=0.4 fixed isolation |
| 1108.2044 | CMS | pp | 7 TeV | Conversion + isolation template fit methods |
| 1311.1440 | ATLAS | pp | 7 TeV | Explicit signal-leakage equation, 100--1000 GeV |
| 1506.08552 | ATLAS | Pb+Pb | 2.76 TeV | T_AA scaling, ABCD in heavy ions, R=0.3 |
| 1605.03495 | ATLAS | pp | 8 TeV | ET-dependent isolation (4.8+0.0042xET), PeTeR |
| 1701.06882 | ATLAS | pp | 13 TeV | First 13 TeV, R_bg validation procedure |
| 1807.00782 | CMS | pp | 13 TeV | BDT-based photon ID for cross section |
| 1903.02209 | ATLAS | p+Pb | 8.16 TeV | R_pPb, nPDF constraints, energy loss limits |
| 1906.01371 | ALICE | pp | 7 TeV | sigma_long^2 ABCD, pT down to 10 GeV |
| 1908.02746 | ATLAS | pp | 13 TeV | NNLO (NNLOJET) comparison, precision benchmark |
| 2003.12797 | CMS | pp+PbPb | 5.02 TeV | R_AA ~ 1, sigma_eta_eta template fit |
| 2407.01165 | ALICE | pp | 13 TeV | Lowest xT, tracks-only isolation, pT from 7 GeV |
| 2409.12641 | ALICE | pp+PbPb | 5.02 TeV | Cone radius ratio, centrality bias, low pT |
| 2502.18054 | ALICE | pp+pPb | 5/8 TeV | R_pA shadowing hint at 2.3 sigma |

## Cross-Links

- [Isolation Cuts](../../concepts/isolation-cuts.md) -- PPG12 parametric isolation implementation
- [ABCD Background Subtraction](../../concepts/abcd-method.md) -- PPG12 signal-leakage equation
- [Shower Shape Variables](../../concepts/shower-shape-variables.md) -- sPHENIX BDT features
- [Prompt Photon Production](../theory/prompt-photon-production.md) -- NLO/NNLO theory framework
- [Isolation Theory](../theory/isolation-theory.md) -- Frixione vs cone isolation
- [Tevatron Photon Measurements](tevatron-photon-measurements.md) -- Predecessor measurements
