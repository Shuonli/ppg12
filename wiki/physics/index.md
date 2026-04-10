# PPG12 Physics Knowledge Base

Domain knowledge compiled from ~70 research papers for the sPHENIX PPG12 isolated photon cross-section analysis.

## Measurements

- [World Data Landscape](measurements/world-data-landscape.md) -- Overview of all isolated/direct photon measurements worldwide, spanning sqrt(s) from 200 GeV to 13 TeV, and where PPG12 fits in the global picture
- [RHIC Photon Measurements](measurements/rhic-photon-measurements.md) -- History of direct/isolated photon measurements at RHIC, from PHENIX's first pp result through thermal photon discovery and the direct photon puzzle
- [LHC Photon Measurements](measurements/lhc-photon-measurements.md) -- ATLAS, CMS, and ALICE isolated photon program at 2.76--13 TeV, establishing the ABCD and parametric isolation methods PPG12 adapts
- [Tevatron Photon Measurements](measurements/tevatron-photon-measurements.md) -- D0 and CDF precision measurements at sqrt(s) = 1.8--1.96 TeV, introducing neural-network photon ID, template fitting, and the low-ET NLO discrepancy
- [Photon-Jet Correlations](measurements/photon-jet-correlations.md) -- Photon-tagged jet measurements from PHENIX through STAR and the LHC, and why the PPG12 pp cross section is the essential baseline for jet quenching
- [Nuclear Modification](measurements/nuclear-modification.md) -- R_AA and R_pA measurements of prompt photons across experiments and energies, from PHENIX Au+Au to the LHC heavy-ion program

## Theory

- [Prompt Photon Production](theory/prompt-photon-production.md) -- Cross-section factorization into direct and fragmentation components, LO subprocesses (Compton, annihilation), and gluon PDF sensitivity at RHIC kinematics
- [Isolation Theory](theory/isolation-theory.md) -- Theoretical foundations of photon isolation: fixed-cone vs Frixione smooth-cone, infrared safety, and residual fragmentation contribution
- [NLO Calculations](theory/nlo-calculations.md) -- JETPHOX and Vogelsang NLO codes used for PPG12 theory comparison, scale choices, PDF and FF inputs, and uncertainty estimation
- [Fragmentation Functions](theory/fragmentation-functions.md) -- BFG set II photon fragmentation functions: anomalous vs non-perturbative components, VMD modeling, and role in NLO predictions

## Techniques

- [Isolation Methods](techniques/isolation-methods.md) -- Comparison of experimental isolation implementations across ATLAS, CMS, ALICE, PHENIX, and sPHENIX, including cone definitions, UE subtraction, and ET-dependent thresholds
- [Shower Shape Discriminants](techniques/shower-shape-discriminants.md) -- Physics of photon vs hadron shower shapes in the sPHENIX W/SciFi EMCal, variable definitions, and how they enter the PPG12 BDT
- [ABCD Sideband Method](techniques/abcd-sideband-method.md) -- Two-dimensional sideband background subtraction with signal leakage corrections, from ATLAS origin to PPG12 implementation
- [BDT Photon Identification](techniques/bdt-photon-id.md) -- XGBoost algorithm, TMVA inference, PPG12 model configuration, and comparison with photon ID strategies at other experiments
- [Bayesian Unfolding](techniques/bayesian-unfolding.md) -- D'Agostini iterative Bayesian unfolding algorithm, RooUnfold implementation, SVD alternative, and PPG12 application
- [Luminosity Calculation](techniques/luminosity-calculation.md) -- MBD scaler-based luminosity with Poisson pileup correction, Vernier scan cross-section, z-vertex selection effects

## Detector

- [sPHENIX Overview](detector/sphenix-overview.md) -- Detector layout, subsystems, and physics program; the upgrade from PHENIX to full barrel calorimetry and tracking
- [EMCal Performance](detector/emcal-performance.md) -- W/SciFi SPACAL design, energy resolution, tower segmentation, calibration, and performance relevant to photon measurement
- [Calorimeter System](detector/calorimeter-system.md) -- Full EMCal + IHCal + OHCal system: specifications, isolation energy sum, and hadronic containment
- [RHIC Accelerator](detector/rhic-accelerator.md) -- RHIC Run 24 parameters, crossing angle periods (0 mrad vs 1.5 mrad), luminosity, and pileup rates
- [Simulation Tools](detector/simulation-tools.md) -- PYTHIA 8.2 event generation and GEANT4 detector simulation chain for PPG12 MC production

## Raw Paper Summaries

Individual paper summaries are in [raw/papers/](../raw/papers/) organized by category.
See the [manifest](../raw/manifest.md) for the complete paper inventory.
