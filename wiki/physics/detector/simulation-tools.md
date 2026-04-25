# Simulation Tools

## Overview

PPG12's Monte Carlo simulation chain has two stages: event generation with PYTHIA 8.2 (producing particle-level final states) and detector simulation with GEANT4 (propagating particles through the sPHENIX detector material). The combined output is processed through the same reconstruction software as real data, producing the slimtree ROOT files used for BDT training, efficiency calculation, and purity estimation.

## PYTHIA 8.2

PYTHIA 8.2 is the primary event generator for all PPG12 MC production. It provides a coherent set of physics models for the complete evolution from a few-body hard QCD process to a complex multiparticle final state.

Reference: arXiv:1410.3012 (Sjostrand et al., Comput. Phys. Commun. 191 (2015) 159)

### Hard Processes for Photon Production

PPG12 signal MC uses PYTHIA's electroweak hard processes for prompt photon production:

- **QCD Compton scattering:** qg -> q gamma (dominant at RHIC energies)
- **Quark-antiquark annihilation:** q qbar -> g gamma

These produce "direct" photons that emerge directly from the hard scattering, carrying the full momentum of one side of the 2->2 process. PYTHIA also produces "fragmentation" photons through the parton shower and hadronization, but these are suppressed by the isolation cut in PPG12.

### Background Processes

PPG12 background MC uses QCD 2->2 hard processes:

- gg -> gg, q qbar -> gg, qg -> qg, qq -> qq, etc.

These jet events produce background photon candidates through pi0 -> gamma gamma and eta -> gamma gamma decays in the fragmentation products. At high pT, the pi0 decay photons can merge in the EMCal to form a single cluster that mimics a direct photon shower shape.

### Parton Distribution Functions

PYTHIA uses leading-order PDFs for the partonic initial state. The Detroit tune uses **NNPDF3.1 LO** with alpha_s(M_Z) = 0.130, updated from the Monash default of NNPDF2.3 LO. PYTHIA's two-PDF option allows separate PDFs for the hard process and for showers/MPI.

### Parton Showers: ISR and FSR

PYTHIA's parton showers use dipole-style pT-ordered evolution:

- **ISR (Initial-State Radiation):** Backward evolution of the incoming partons, producing additional radiation before the hard scatter. ISR photons and gluons contribute to the underlying event activity in the isolation cone.
- **FSR (Final-State Radiation):** Forward evolution of outgoing partons, producing the shower that develops into jets. FSR determines the spatial distribution of energy within jets, affecting shower shapes when jet fragments hit the EMCal.

The splitting kernels include standard LO DGLAP: q -> qg, g -> gg, g -> q qbar (QCD) and f -> f gamma, gamma -> f fbar (QED). The effective alpha_s(M_Z) is ~0.139 (larger than the MSbar value due to CMW scheme translation).

### Multiparton Interactions (MPI)

In PYTHIA, multiple parton-parton interactions occur within a single pp collision because protons are composite objects. MPI produces the underlying event (UE) -- soft particles that fill the isolation cone and affect the photon isolation measurement.

The MPI cross section is regulated by the replacement:

```
alpha_s^2/pT^4 -> alpha_s^2(pT^2 + pT0^2) / (pT^2 + pT0^2)^2
```

where pT0 is the key tunable parameter that controls the amount of soft MPI activity. The Detroit tune sets pT0Ref = 1.40 GeV at sqrt(s) = 200 GeV.

MPI, ISR, and FSR are combined into one interleaved pT-ordered evolution sequence, where the next step is determined by whichever mechanism produces the highest pT emission.

### String Fragmentation

PYTHIA uses the Lund string fragmentation model exclusively for hadronization. Colour flux tubes between quarks are modeled as relativistic strings with tension kappa ~ 1 GeV/fm. String breaking produces new quark-antiquark pairs via quantum tunneling, with flavour suppression u:d:s:c ~ 1:1:0.3:10^-11.

Gluons appear as kinks on the string, carrying energy and momentum. The colour reconnection mechanism (Detroit tune: range = 5.4) rearranges colour flow between different MPI systems to minimize total string length, reducing particle multiplicity.

## Detroit Tune

The Detroit tune (arXiv:2110.09447, Phys. Rev. D 105 (2022) 016011) is a PYTHIA 8 underlying event tune developed by STAR specifically for RHIC energies. It supersedes the default Monash tune for sqrt(s) = 200 GeV applications. All PPG12 PYTHIA MC uses the Detroit tune.

### Tuned Parameters

| Parameter | Monash (default) | Detroit | Physical meaning |
|-----------|-----------------|---------|-----------------|
| pT0Ref | 2.28 GeV (at 7 TeV) | 1.40 GeV (at 200 GeV) | MPI regularization scale |
| ecmPow | 0.215 | 0.135 | Energy scaling of pT0 |
| coreRadius | 0.4 | 0.56 | Proton core radius (double Gaussian) |
| coreFraction | 0.5 | 0.78 | Fraction in proton core |
| CR:range | 1.8 | 5.4 | Colour reconnection range |

### Key Differences from Monash

1. **Reference energy changed from 7 TeV to 200 GeV.** The Monash tune was optimized for LHC, then extrapolated to RHIC with a power-law formula. The Detroit tune directly targets 200 GeV, eliminating extrapolation uncertainty at RHIC.

2. **Double Gaussian proton shape.** Monash uses a default profile (bProfile=3); Detroit uses a double Gaussian (bProfile=2) that better describes the impact-parameter-dependent MPI rate at low sqrt(s).

3. **Increased colour reconnection.** Range = 5.4 (3x Monash) allows more extensive rearrangement of colour flow, reducing string lengths and lowering particle multiplicity.

4. **pT0 energy scaling is flatter.** ecmPow = 0.135 (vs 0.215) means pT0 varies more slowly with energy, preventing overshoot when extrapolating to higher energies.

### Performance at sqrt(s) = 200 GeV

- Charged pion pT spectra: well described across full measured pT range
- Underlying event multiplicities: significant improvement over Monash in all regions (toward, transverse, away)
- Jet substructure (SoftDrop groomed observables): improved agreement
- **Forward rapidity limitation:** Both Monash and Detroit fail to describe pion data at eta > 2.9 (not relevant for PPG12 at |eta| < 0.7)

### Eigentunes

Ten eigentune variations (5 principal directions, +/- each) are provided for systematic uncertainty assessment:

| Parameter | Range across eigentunes |
|-----------|------------------------|
| pT0Ref | 1.37--1.44 GeV |
| ecmPow | 0.119--0.150 |
| coreRadius | 0.41--0.77 |
| coreFraction | 0.60--0.90 |
| CR:range | 3.61--7.50 |

PPG12 currently does not use eigentune variations for generator-level systematics, but they are available for future studies.

### Impact on PPG12

The Detroit tune directly affects several PPG12 observables:

- **Isolation cone energy:** MPI activity in the R = 0.3 cone determines the UE contribution that the parametric isolation cut (iso_max = b + s * ET) must absorb.
- **Jet shower shapes:** Colour reconnection and MPI control how background jet fragments are distributed in the EMCal, determining the BDT training distributions for all 25 shower-shape features.
- **Signal/background yields:** Cross-section weights depend on PYTHIA's total cross-sections at 200 GeV.
- **Fragmentation photon rate:** The rate of high-pT photons from jet fragmentation (pi0 -> gamma gamma) affects the irreducible background in the isolation signal region.

## GEANT4

GEANT4 (NIM A 506 (2003) 250, Agostinelli et al.) provides the full detector simulation for all sPHENIX MC. It propagates PYTHIA-generated particles through the detailed sPHENIX detector geometry, modeling all physics interactions with detector materials.

### Physics Processes

Key electromagnetic processes for EMCal shower simulation:
- Photoelectric effect
- Compton scattering
- Pair production (gamma -> e+e-)
- Bremsstrahlung (e -> e gamma)
- Ionization and multiple scattering

Hadronic processes for HCal and isolation cone:
- Nuclear interactions at multiple energy scales (QGSP_BERT_HP physics list)
- Birks' law for scintillation quenching: kB = 0.0794 mm/MeV

### Digitization Chain

The simulation includes the full digitization chain from energy deposition to detector signals:

1. Energy deposition in active volumes (scintillating fibers/tiles)
2. Birks' law quenching correction
3. Visible energy to photoelectrons (fiber attenuation, light guide efficiency)
4. SiPM pixel saturation (Poisson model)
5. ADC conversion

Timing window for energy summation: 0--60 ns.

### Additional Smearing

An additional 11.9% energy smearing is applied to EMCal clusters in simulation to match the pi0 mass peak width in data. This accounts for data-MC resolution differences from effects not fully modeled in the simulation (cable delays, electronics noise, imperfect calibration).

## MC Sample Production

### Generation Workflow

1. **PYTHIA 8** generates hard-scattering events with Detroit tune parameters at sqrt(s) = 200 GeV
2. **GEANT4** propagates particles through the sPHENIX detector geometry
3. **Reconstruction** produces calibrated clusters from energy deposits
4. **CaloAna24.cc** writes slimtree ROOT files with cluster-level variables

The Fun4All framework (sPHENIX software framework) orchestrates the PYTHIA -> GEANT4 -> reconstruction chain. MC samples are produced at "Run 28 conditions" with CDB global tag MDC2.

### Signal Samples

| Sample | pT hat range | Physics | Purpose |
|--------|-------------|---------|---------|
| photon5 | 0--14 GeV | Prompt photon (Compton + annihilation) | Low-pT signal |
| photon10 | 14--30 GeV | Same | Mid-pT signal (primary) |
| photon20 | 30+ GeV | Same | High-pT signal |

The pT hat ranges define the PYTHIA-level hard-process pT window. The three samples are combined with cross-section weights to form a smooth signal spectrum across the full PPG12 range (8--35 GeV).

### Background Samples

| Sample | Jet pT range | Purpose |
|--------|-------------|---------|
| jet5 | 7--9 GeV | Low-pT background |
| jet8 | 9--14 GeV | |
| jet12 | 14--21 GeV | |
| jet20 | 21--32 GeV | |
| jet30 | 32--42 GeV | |
| jet40/jet50 | 42+ GeV | High-pT background |

Background samples use QCD 2->2 processes. Clusters from these events that pass photon-ID cuts are predominantly from pi0/eta -> gamma gamma decays where the decay photons merge in the EMCal.

### Double-Interaction Samples

Full GEANT4 double-interaction samples overlay two independent pp collisions. After the Apr 2026 reprocess, 8 SI/DI pairs are available:

| DI sample | Partner SI | Purpose |
|-----------|-----------|---------|
| photon5_double | photon5 | Signal with pileup (0--14 GeV) |
| photon10_double | photon10 | Signal with pileup (14--30 GeV) |
| photon20_double | photon20 | Signal with pileup (30+ GeV) |
| jet8_double | jet8 | Background with pileup (9--14 GeV jet) |
| jet12_double | jet12 | Background with pileup (14--21 GeV jet) |
| jet20_double | jet20 | Background with pileup (21--32 GeV jet) |
| jet30_double | jet30 | Background with pileup (32--42 GeV jet) |
| jet40_double | jet40 | Background with pileup (42--100 GeV jet) |

Mixing fractions (all pairs): 22.4% (0 mrad), 7.9% (1.5 mrad), cluster-weighted.

These are blended with single-interaction samples in the single-pass truth-vertex-reweight showershape pipeline (`submit_showershape_di.sub`) to study pileup effects on shower shapes and purity. See [Double-Interaction Efficiency](../../concepts/double-interaction-efficiency.md) and [MC Samples](../../reference/mc-samples.md#double-interaction-samples).

### Events per Sample

Default: 10^7 events per sample (from `CrossSectionWeights.h`), sufficient for the BDT training and efficiency calculations across all pT bins.

## Photon Production Processes in Detail

At sqrt(s) = 200 GeV, prompt photon production proceeds through:

### Direct Photons (hard process)

1. **QCD Compton:** qg -> q gamma -- dominant at RHIC because of the gluon-rich proton. The outgoing photon is balanced by a quark jet.
2. **Annihilation:** q qbar -> g gamma -- subdominant, requires valence quark-antiquark pair. Becomes relatively more important at forward rapidity.

At mid-rapidity (|eta| < 0.7), the Compton process dominates for pT > 8 GeV, probing the gluon PDF at x ~ 2*pT/sqrt(s) ~ 0.08--0.35.

### Fragmentation Photons

Photons radiated from final-state partons during the QCD shower (q -> q gamma in FSR). These carry a fraction z of the parent parton's momentum and are typically accompanied by nearby hadronic activity, making them suppressible by the isolation cut. However, at small cone sizes (R = 0.3) and high z, some fragmentation photons pass the isolation requirement and constitute an irreducible background.

### Decay Photons (background)

The dominant background: neutral mesons (pi0, eta, omega) produced in jet fragmentation decay to photon pairs. At pT > 20 GeV, the opening angle of pi0 -> gamma gamma becomes smaller than the EMCal tower size, causing the two photons to merge into a single cluster. This merged cluster is the primary target of the BDT shower-shape discrimination.

## Modified Power-Law Spectrum Fitting

The combined MC photon pT spectrum is fit with the Spousta-Cole modified power law (Eur. Phys. J. C 76 (2016) 50):

```
d sigma / d pT = A * pT^(-n) * (1 + pT/p0)^(-m)
```

This captures the curvature of the spectrum in log-log space arising from alpha_s running and PDF evolution. A pure power law (pT^(-n)) is insufficient across the PPG12 range of 8--35 GeV. The fit provides smooth interpolation for ratio construction and bin-width corrections.

## References

- PYTHIA 8.2: arXiv:1410.3012 (Comput. Phys. Commun. 191 (2015) 159)
- Detroit tune: arXiv:2110.09447 (Phys. Rev. D 105 (2022) 016011)
- GEANT4: NIM A 506 (2003) 250 (Agostinelli et al.)
- Modified power law: Eur. Phys. J. C 76 (2016) 50 (Spousta, Cole)
- PHENIX clustering heritage: IEEE TNS 47 (2000) 1982 (David et al.)

## Cross-Links

- [sPHENIX Overview](sphenix-overview.md) -- detector geometry for GEANT4 simulation
- [EMCal Performance](emcal-performance.md) -- simulated vs measured resolution
- [Tree Making](../../pipeline/01-tree-making.md) -- CaloAna24.cc processes MC and data identically
- [BDT Training](../../pipeline/02-bdt-training.md) -- training uses PYTHIA+GEANT4 MC
- [MC Samples](../../reference/mc-samples.md) -- sample inventory, cross-sections, file locations
