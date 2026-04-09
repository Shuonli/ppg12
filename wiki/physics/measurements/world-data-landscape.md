---
title: "World Data Landscape: Isolated and Direct Photon Measurements"
scope: "Overview of all isolated/direct photon measurements worldwide, spanning sqrt(s) from 200 GeV to 13 TeV, across RHIC, Tevatron, and LHC experiments, with emphasis on NLO pQCD agreement and where PPG12 fits in the global picture"
source_papers: [hep-ex/0502006, 1205.5533, 1006.1347, 2202.08158, hep-ex/0511054, 0910.3623, CDF_PRD48_1993, 1108.0253, 1108.2044, 1311.1440, 1605.03495, 1701.06882, 1908.02746, 1807.00782, 1906.01371, 2407.01165, 2003.12797]
related_articles:
  - ../techniques/isolation-methods.md
  - ../../concepts/isolation-cuts.md
  - ../../concepts/abcd-method.md
  - ../../pipeline/overview.md
  - rhic-photon-measurements.md
---

# World Data Landscape: Isolated and Direct Photon Measurements

Isolated and direct photon production in hadronic collisions is one of the cleanest probes of perturbative QCD. The photon couples directly to the hard scattering, avoiding fragmentation ambiguities that affect hadron measurements. Over four decades, measurements spanning center-of-mass energies from sqrt(s) ~ 20 GeV (fixed target) to 13 TeV (LHC) have established that NLO pQCD describes isolated photon cross sections to better than 10-20% across this entire kinematic range. PPG12, the sPHENIX isolated photon measurement at sqrt(s) = 200 GeV, sits at the low-energy frontier of this global dataset, probing the highest xT values accessible with modern isolation techniques.

## Energy Landscape and Kinematic Coverage

The world dataset of isolated/direct photon measurements spans three orders of magnitude in sqrt(s):

| Facility  | Experiment | sqrt(s)          | pT range (GeV/c)  | |eta| coverage | Key reference |
|-----------|-----------|------------------|--------------------|--------------------|---------------|
| RHIC      | PHENIX    | 200 GeV (pp)     | 5.5--25            | < 0.25             | [1205.5533](../../raw/papers/A_rhic/1205.5533.md) |
| RHIC      | PHENIX    | 510 GeV (pp)     | 6--30              | < 0.25             | [2202.08158](../../raw/papers/A_rhic/2202.08158.md) |
| RHIC      | sPHENIX   | 200 GeV (pp)     | 8--35              | < 0.7              | PPG12 (this analysis) |
| Tevatron  | D0        | 1.96 TeV (ppbar) | 23--300            | < 0.9              | [hep-ex/0511054](../../raw/papers/A_tevatron/hep-ex_0511054.md) |
| Tevatron  | CDF       | 1.96 TeV (ppbar) | 30--400            | < 1.0              | [0910.3623](../../raw/papers/A_tevatron/0910.3623.md) |
| LHC       | ALICE     | 7 TeV (pp)       | 10--60             | < 0.27             | [1906.01371](../../raw/papers/A_lhc/1906.01371.md) |
| LHC       | CMS       | 7 TeV (pp)       | 25--400            | < 2.5              | [1108.2044](../../raw/papers/A_lhc/1108.2044.md) |
| LHC       | ATLAS     | 7 TeV (pp)       | 45--400            | < 2.37             | [1108.0253](../../raw/papers/A_lhc/1108.0253.md) |
| LHC       | CMS       | 5.02 TeV (pp)    | 25--200            | < 1.44             | [2003.12797](../../raw/papers/A_lhc/2003.12797.md) |
| LHC       | ATLAS     | 13 TeV (pp)      | 125--2500          | < 2.37             | [1908.02746](../../raw/papers/A_lhc/1908.02746.md) |
| LHC       | CMS       | 13 TeV (pp)      | 190--1000          | < 2.5              | [1807.00782](../../raw/papers/A_lhc/1807.00782.md) |
| LHC       | ALICE     | 13 TeV (pp)      | 7--200             | < 0.67             | [2407.01165](../../raw/papers/A_lhc/2407.01165.md) |

The xT = 2pT/sqrt(s) variable unifies these measurements onto a common scaling curve. At xT ~ 0.05--0.15 (PPG12's range), the cross section probes valence quarks and moderate-x gluons, where PDFs are well-constrained. At xT ~ 10^{-3} (ALICE 13 TeV at low pT), the measurement probes the low-x gluon and becomes sensitive to PDF uncertainties. The xT scaling analysis by PHENIX ([1205.5533](../../raw/papers/A_rhic/1205.5533.md)) with effective power n_eff ~ 4.5 demonstrates universality of the pQCD subprocess picture from ISR to LHC energies.

## NLO pQCD Agreement Across Energies

A central result across the world dataset is that NLO pQCD, as implemented in JETPHOX, describes isolated photon cross sections to within 10--20% across the full kinematic range. The level of agreement improves with increasing pT relative to the measurement's lower bound.

### RHIC (200--510 GeV)

At sqrt(s) = 200 GeV, the PHENIX inclusive direct photon measurement ([1205.5533](../../raw/papers/A_rhic/1205.5533.md)) found NLO pQCD (Vogelsang, CTEQ6M, BFGII, scales mu = pT/2, pT, 2pT) in good agreement with the data across 5.25--25 GeV/c. The power-law fit to the data (n = 7.08 +/- 0.09 for 8 < pT < 25 GeV/c) is consistent with the pQCD prediction. The ratio of isolated to inclusive direct photons exceeds 90% above 10 GeV/c, confirming that the fragmentation component is small.

At sqrt(s) = 510 GeV ([2202.08158](../../raw/papers/A_rhic/2202.08158.md)), an important observation emerged: NLO pQCD underestimates the **inclusive** (non-isolated) direct photon cross section by up to a factor of ~3 at pT < 12 GeV/c, attributed to multiparton interactions and parton shower contributions not captured at fixed order. However, NLO agrees well with the **isolated** photon cross section across the full pT range. This underscores that isolation removes the soft/multiparton contributions that fixed-order calculations cannot describe, making isolated photon measurements the cleanest comparison to NLO theory.

### Tevatron (1.96 TeV)

Both D0 ([hep-ex/0511054](../../raw/papers/A_tevatron/hep-ex_0511054.md)) and CDF ([0910.3623](../../raw/papers/A_tevatron/0910.3623.md)) measured isolated photon cross sections at sqrt(s) = 1.96 TeV with JETPHOX NLO comparisons. D0 found agreement across 23--300 GeV/c. CDF found agreement for ET > 50 GeV (chi2 p-value = 21%), but observed a persistent discrepancy at ET < 50 GeV where data has a steeper slope than NLO. This low-ET excess is a recurring feature across experiments and energy scales, suggesting missing higher-order or non-perturbative contributions at the lower edge of the pT spectrum.

The Tevatron measurements also established the non-perturbative hadronization correction C_had ~ 0.91 +/- 0.03 (CDF), quantifying the effect of underlying event energy deposited in the isolation cone. This correction, approximately ET-independent, reduces the NLO-predicted cross section by ~9%.

### LHC (5.02--13 TeV)

The LHC measurements represent the highest precision and widest kinematic reach. ATLAS at 13 TeV ([1908.02746](../../raw/papers/A_lhc/1908.02746.md)) with 36 fb^{-1} spans 125--2500 GeV across four eta regions, covering seven orders of magnitude in cross section. The key innovations at this precision level include:

- **NNLO comparison**: NNLOJET provides the first NNLO QCD calculation for direct photon production. Scale uncertainties are reduced by a factor of 2--20 compared to NLO (0.6--5% vs 10--15%), with excellent data description.
- **Scale choice sensitivity**: NLO JETPHOX with mu = ET tends to underestimate data by 5--15%; mu = ET/2 improves normalization. The Sherpa ME+PS@NLO generator, which includes parton shower matching, provides improved normalization over fixed-order JETPHOX.
- **PDF discrimination**: At the percent-level precision of NNLO, the data can discriminate between PDF sets (MMHT2014, CT14, NNPDF3.0, HERAPDF2.0, ABMP16).

ALICE extends the LHC pT reach to the lowest values: 10 GeV/c at 7 TeV ([1906.01371](../../raw/papers/A_lhc/1906.01371.md)) and 7 GeV/c at 13 TeV ([2407.01165](../../raw/papers/A_lhc/2407.01165.md)). The 13 TeV measurement achieves the lowest xT (1.08 x 10^{-3}) of any isolated photon measurement, probing the low-x gluon PDF. Both ALICE measurements find NLO JETPHOX agreement within uncertainties.

CMS contributes measurements at both 7 TeV ([1108.2044](../../raw/papers/A_lhc/1108.2044.md)) and 13 TeV ([1807.00782](../../raw/papers/A_lhc/1807.00782.md)), with the latter being the first CMS measurement to employ BDT-based photon identification -- directly analogous to PPG12's approach.

## Evolution of Background Subtraction Techniques

The techniques for separating prompt photons from the dominant hadronic decay background have evolved significantly across experiments.

### Statistical subtraction (PHENIX, early)

The earliest RHIC measurements ([hep-ex/0502006](../../raw/papers/A_rhic/hep-ex_0502006.md)) used the double-ratio R_gamma = (gamma_cand/pi0)_meas / (gamma_cand/pi0)_bkg to extract the direct photon signal statistically. This requires precise knowledge of the pi0 spectrum and the gamma/pi0 ratio from hadronic decays (eta/pi0 = 0.48, mT-scaling for other mesons). The method is elegant but suffers from large systematic amplification at low pT where the signal/background ratio is small (amplification factor W = N_incl/N_dir).

### Isolation template fitting (D0, CDF, CMS)

D0 ([hep-ex/0511054](../../raw/papers/A_tevatron/hep-ex_0511054.md)) introduced neural network-based photon discrimination with 4 shower-shape variables, extracting purity via template fitting of the NN output. CDF ([0910.3623](../../raw/papers/A_tevatron/0910.3623.md)) developed isolation-distribution template fitting, exploiting the sharp peak at low isolation for signal versus the flat distribution for backgrounds. CMS ([1108.2044](../../raw/papers/A_lhc/1108.2044.md)) used two complementary methods (conversion-based ET/pT fitting and isolation template fitting).

### ABCD / 2D sideband (ATLAS, ALICE, PPG12)

The most refined approach, used by ATLAS in all papers from 7 TeV onward ([1108.0253](../../raw/papers/A_lhc/1108.0253.md) through [1908.02746](../../raw/papers/A_lhc/1908.02746.md)) and by ALICE ([1906.01371](../../raw/papers/A_lhc/1906.01371.md), [2407.01165](../../raw/papers/A_lhc/2407.01165.md)), defines four regions in the (shower-shape, isolation) plane:

```
            Isolated     Non-isolated
Tight       A (signal)   B
Non-tight   C            D
```

The background in A is estimated as R_bg * B * C / D, with signal leakage corrections. ATLAS validates R_bg using control region subdivisions, finding deviations of 15--40% from unity that are propagated as systematics. ALICE applies an MC correction factor alpha_MC to account for correlations between shower shape and isolation for background clusters.

PPG12 adopts the same ABCD framework with BDT score (tight/non-tight) and parametric isolation (isolated/non-isolated), including signal leakage corrections c_B, c_C, c_D. See [ABCD Method](../../concepts/abcd-method.md) for PPG12-specific details.

### BDT-based identification (CMS 13 TeV, PPG12)

CMS at 13 TeV ([1807.00782](../../raw/papers/A_lhc/1807.00782.md)) was the first experiment to use a BDT for photon identification in a cross section measurement, training separate BDTs per rapidity and ET bin using shower-shape variables. PPG12 extends this approach with a 25-feature XGBoost BDT trained on GEANT MC, applied as the tight/non-tight axis in the ABCD method.

## Isolation Definitions Across Experiments

The isolation cone definition varies across experiments but follows a common principle: sum the energy/momentum of all particles (or tracks) in a cone around the photon candidate, excluding the photon itself.

| Experiment | Cone R | Isolation criterion | Type |
|------------|--------|---------------------|------|
| PHENIX (200 GeV) | 0.5 | E_cone < 10% E_gamma | Fractional |
| PHENIX (510 GeV) | 0.5 | E_cone < 10% E_gamma | Fractional |
| D0 | 0.4 | (E_tot(0.4) - E_EM(0.2))/E_EM(0.2) < 10% | Fractional |
| CDF | 0.4 | E_T^iso < 2 GeV | Absolute |
| ATLAS (7 TeV early) | 0.4 | E_T^iso < 3 GeV | Absolute |
| ATLAS (8--13 TeV) | 0.4 | E_T^iso < 4.8 + 0.0042 ET | Parametric |
| CMS (7 TeV) | 0.4 | sum ET < 5 GeV (gen-level) | Absolute |
| CMS (5.02 TeV) | 0.4 | I < 1 GeV (after UE subtraction) | Absolute |
| ALICE (7 TeV) | 0.4 | p_T^iso < 2 GeV (tracks + clusters) | Absolute |
| ALICE (13 TeV) | 0.4 | p_T^iso,ch < 1.5 GeV (tracks only) | Absolute |
| PPG12 (sPHENIX) | 0.3 | reco_iso_max = b + s * ET | Parametric |

The evolution toward ET-dependent (parametric) isolation cuts, pioneered by ATLAS at 8 TeV and adopted by PPG12, provides more robust performance across a wide ET range by scaling the allowed isolation energy with the photon energy. This naturally accommodates the ET dependence of underlying event fluctuations and photon shower leakage. See [Isolation Cuts](../../concepts/isolation-cuts.md) for PPG12-specific definitions.

## Systematic Uncertainty Hierarchy

Across all experiments, the dominant systematic uncertainties follow a consistent hierarchy:

1. **Energy scale** (5--15%): The steeply falling cross section amplifies small absolute energy calibration uncertainties into large cross-section uncertainties. ATLAS uses 76 independent energy scale components; PPG12 uses a few energy scale variations.
2. **Purity / signal extraction** (3--30%): Dominant at low pT where signal/background is small. At ALICE, the isolation probability systematic reaches 28% at 7 GeV/c. At CDF, the photon fraction systematic is 13% at low ET. PPG12 faces similar challenges in the 8--12 GeV/c range.
3. **Luminosity** (1.5--10%): A global normalization uncertainty. PPG12 uses 16.6 pb^{-1} with the sPHENIX luminosity monitoring system.
4. **Photon identification efficiency** (1--7%): Efficiency corrections and their data/MC scale factors.
5. **Theory scale uncertainty** (0.6--15%): The dominant theoretical systematic. Ranges from 10--15% at NLO to 0.6--5% at NNLO.

## Where PPG12 Fits

PPG12 occupies a unique position in this landscape:

- **Lowest sqrt(s) with modern isolation techniques**: At 200 GeV, PPG12 probes xT = 0.08--0.35, the highest xT values measured with the ABCD method. Previous measurements at this energy (PHENIX) used statistical subtraction rather than isolation-based techniques.
- **Bridge between RHIC and Tevatron/LHC**: PPG12 tests NLO pQCD at an energy scale where fixed-target experiments show tensions but modern collider measurements have not yet reached with isolation techniques.
- **Detector analog to ALICE**: The sPHENIX CEMCal (lead-scintillator sampling, cell size ~2.4 cm, |eta| < 0.7) is more similar to the ALICE EMCal (lead-scintillator, cell size ~6 cm, |eta| < 0.7) than to the ATLAS/CMS crystal calorimeters. PPG12's shower-shape challenges -- cross-talk, merged clusters, moderate resolution -- are shared with ALICE.
- **BDT + ABCD synergy**: PPG12 combines the BDT-based identification pioneered by CMS 13 TeV with the ABCD background subtraction developed by ATLAS, integrating both state-of-the-art techniques.
- **pp baseline for heavy-ion program**: Unlike LHC isolated photon measurements in PbPb (CMS [2003.12797](../../raw/papers/A_lhc/2003.12797.md) found R_AA ~ 1), the RHIC program has shown R_AA ~ 1 for high-pT direct photons but a large thermal excess at low pT. PPG12 establishes the precise pp baseline that future sPHENIX Au+Au measurements will reference.

## xT Scaling and Universality

When the cross section is expressed as a function of xT and scaled by sqrt(s)^{n_eff}, measurements across all energies collapse onto a universal curve. PHENIX determined n_eff = 4.5 from a global fit spanning ISR (sqrt(s) = 19.4--63 GeV), RHIC (200 GeV), Tevatron (1.8--1.96 TeV), and LHC (7 TeV) data ([1205.5533](../../raw/papers/A_rhic/1205.5533.md)). ALICE extended this scaling to 13 TeV ([2407.01165](../../raw/papers/A_lhc/2407.01165.md)). The consistency of this scaling confirms that the underlying pQCD subprocess (dominated by qg Compton scattering at RHIC energies, with ~80% contribution at pT = 10 GeV/c) is correctly described across the entire energy range.

PPG12's measurement at xT = 0.08--0.35 will add high-precision sPHENIX data points to this scaling curve, testing universality at the highest accessible xT values with a modern detector.

## See Also

- [RHIC Photon Measurements](rhic-photon-measurements.md) -- Detailed history of RHIC direct/isolated photon measurements
- [ABCD Method](../../concepts/abcd-method.md) -- PPG12's background subtraction implementation
- [Isolation Cuts](../../concepts/isolation-cuts.md) -- PPG12's parametric isolation definition
- [Pipeline Overview](../../pipeline/overview.md) -- PPG12 analysis pipeline

## Source Papers

### RHIC
- [PHENIX pp 200 GeV first measurement (2007)](../../raw/papers/A_rhic/hep-ex_0502006.md)
- [PHENIX pp 200 GeV inclusive/isolated (2012)](../../raw/papers/A_rhic/1205.5533.md)
- [PHENIX pp 200 GeV isolated photon correlations (2010)](../../raw/papers/A_rhic/1006.1347.md)
- [PHENIX pp 510 GeV cross section and ALL (2023)](../../raw/papers/A_rhic/2202.08158.md)

### Tevatron
- [D0 ppbar 1.96 TeV (2006)](../../raw/papers/A_tevatron/hep-ex_0511054.md)
- [CDF ppbar 1.96 TeV (2009)](../../raw/papers/A_tevatron/0910.3623.md)
- [CDF ppbar 1.8 TeV (1993)](../../raw/papers/A_tevatron/CDF_PRD48_1993.md)

### LHC
- [ATLAS pp 7 TeV 35 pb^-1 (2011)](../../raw/papers/A_lhc/1108.0253.md)
- [CMS pp 7 TeV 36 pb^-1 (2011)](../../raw/papers/A_lhc/1108.2044.md)
- [ALICE pp 7 TeV (2019)](../../raw/papers/A_lhc/1906.01371.md)
- [ATLAS pp 13 TeV 36 fb^-1 (2019)](../../raw/papers/A_lhc/1908.02746.md)
- [CMS pp 13 TeV with BDT (2019)](../../raw/papers/A_lhc/1807.00782.md)
- [ALICE pp 13 TeV (2024)](../../raw/papers/A_lhc/2407.01165.md)
- [CMS PbPb+pp 5.02 TeV (2020)](../../raw/papers/A_lhc/2003.12797.md)
