# EMCal Performance

## Design Overview

The sPHENIX electromagnetic calorimeter (EMCal) is a tungsten powder / scintillating fiber (W/SciFi) sampling calorimeter of the SPACAL (Spaghetti Calorimeter) type. It is the primary detector subsystem for PPG12, providing the photon energy measurement, cluster position, and shower-shape variables that drive the BDT photon identification.

### Construction

Each EMCal block contains scintillating fibers embedded in a matrix of tungsten powder and epoxy:

| Parameter | Value |
|-----------|-------|
| Fibers per block | 2668 (Saint Gobain BCF-12, 0.47 mm diameter, single cladding) |
| Fiber emission peak | 435 nm; decay time 3.2 ns |
| Fiber attenuation length | >= 1.6 m |
| Absorber | Tungsten powder (THP Technon 100 mesh, >= 99% W purity) + EPO-TEK 301 epoxy |
| Block density | 9.2--9.8 g/cm^3 (approximately 9.5 g/cm^3) |
| Sampling fraction | ~2.1% (electromagnetic showers) |

### Geometry

| Parameter | Value |
|-----------|-------|
| Tower segmentation | Delta_eta x Delta_phi = 0.025 x 0.025 |
| Tower volume | ~2.5 x 2.5 x 14 cm^3 |
| Total channels | ~25,000 (256 in phi x 96 in eta) |
| Coverage | \|eta\| < 1.1, full 2pi azimuth |
| Depth | ~20 radiation lengths (X0 ~ 7 mm, total ~14 cm) |
| Hadronic interaction length | ~0.83 lambda_I |
| Radiation length (X0) | ~7 mm |
| Moliere radius (R_M) | ~2.3 cm |
| Projectivity | 2D projective (tapered in both eta and phi) |

The 2D projective design means each tower points approximately toward the interaction point. Fibers are tilted slightly off the direct IP line to prevent channeling effects where particles could travel along fiber axes. The tower cross-section of approximately (1.1 R_M)^2 means most of a single photon's electromagnetic shower energy is contained within a 3x3 tower cluster.

### Light Collection and Readout

- **Light guides:** UV-transmitting acrylic, trapezoidal shape, epoxied to front end of blocks
- **Back reflectors:** Aluminum, epoxied to block rear
- **SiPMs:** Hamamatsu S12572-015P, 4 per tower (2x2 array), 3x3 mm^2 active area, 40,000 pixels (15 um pitch)
- **SiPM operation:** 4 V above breakdown, gain ~2.3 x 10^5, PDE = 25%
- **Light yield:** ~500 photoelectrons/GeV
- **Light guide efficiency:** 71% relative to direct PMT coupling
- **Position-dependent non-uniformity:** ~30% variation across light guide input cross-section, producing ~20% position-dependent energy response within a tower
- **Effective scintillation attenuation length in blocks:** L_eff = 125 +/- 11 cm (much longer than block length, so longitudinal response is approximately uniform)
- **Temperature dependence:** -3.68 +/- 0.29 %/degC uncorrected; 0.33 +/- 0.30 %/degC after offline compensation (effectively negligible)

## Beam Test Results

Two beam test campaigns at the Fermilab Test Beam Facility (T-1044) established the EMCal performance:

### 2016 Beam Test: 1D Projective Prototype at eta = 0

Tested with the original 1D projective geometry (tapered in phi only), at normal incidence.

**Electron energy resolution (full tower, 2.5 x 2.5 cm^2 area):**

```
sigma(E)/E = 2.8% oplus 15.5%/sqrt(E [GeV])
```

This is the "operational" resolution that includes the effects of position-dependent response across the full tower face. At center-of-tower (1.0 x 0.5 cm^2 hodoscope selection), the intrinsic resolution is better: 1.6% + 12.7%/sqrt(E) for the UIUC blocks.

The fit form accounts for beam momentum spread: sigma(E)/E = sqrt((delta_p/p)^2 + a^2 + b^2/E), with delta_p/p ~ 2% unfolded.

Reference: arXiv:1704.01461

### 2018 Beam Test: 2D Projective Prototype at eta = 1

Tested the final sPHENIX design with 2D projectivity (tapered in both eta and phi), centered at eta = 1.

**Electron energy resolution (tower-averaged, full tower area):**

| Correction method | Constant term a | Stochastic term b |
|-------------------|-----------------|-------------------|
| Hodoscope correction | 3.5 +/- 0.1% | 13.3 +/- 0.2% sqrt(GeV) |
| Cluster correction | 3.0 +/- 0.1% | 15.4 +/- 0.3% sqrt(GeV) |

At the center of the tower (1.0 x 0.5 cm^2 cut), the resolution improves to ~2.3--2.7% + 12.3--13.4%/sqrt(E).

**Key improvements over 2016:**
- 2D projectivity provides better shower containment at forward rapidity
- Stochastic term improved by ~2.5% (from 15.5% to 13.3%)
- Constant term increased slightly (from 2.8% to 3.5%)
- Linearity improved by ~1%

**Nonlinearity:** A quadratic nonlinearity term c ~ -10^-3 GeV^-1 from SiPM saturation affects the energy scale at the highest PPG12 pT bins (above ~25 GeV), corrected during calibration.

Reference: arXiv:2003.13685

### Summary of Resolution Measurements

The beam test resolutions bracket the in-situ performance:

```
Full tower, eta=0 (2016):  sigma/E = 2.8% + 15.5%/sqrt(E)
Full tower, eta=1 (2018):  sigma/E = 3.0-3.5% + 13.3-15.4%/sqrt(E)
```

The first measurement (eta = 0, full tower) is the most commonly cited resolution for the sPHENIX EMCal.

## In-Situ Calibration

### Pi0 Mass Calibration

The primary energy scale calibration uses the pi0 -> gamma gamma decay:

1. **Relative phi-symmetry calibration:** Tower-by-tower energy calibration for uniform phi response
2. **Cluster pair reconstruction:** Photon candidates form diphoton invariant mass M_gamma_gamma
3. **Iterative eta-dependent calibration:** The pi0 mass peak position is shifted to the expected value (134.98 MeV) based on the eta of the most energetic tower in the higher-energy cluster
4. **Iteration until convergence:** The procedure is repeated until calibration constants are stable across all eta rings

This calibration procedure was validated with the first Au+Au collision data (arXiv:2410.18031) and establishes the EMCal energy scale used for PPG12 photon cluster reconstruction.

### Additional Cluster Energy Smearing

In simulation, an additional 11.9% energy smearing is applied to EMCal clusters to match the pi0 mass peak width observed in data (arXiv:2504.02242). This data-MC resolution mismatch is the basis for the energy resolution systematic variations (`eres_up`/`eres_down`) in PPG12.

## Impact on PPG12

### Shower Shape Variables

The tower granularity (0.025 x 0.025) and Moliere radius (~2.3 cm) define the shower-shape variables used in the BDT:

- **Energy ratios** (e11/e33, e32/e35, e11/e22, e11/e13, etc.): Ratios of energy in NxM tower grids centered on the cluster. Single photon showers are more compact (higher e11/e33) than pi0 or hadronic clusters.
- **Tower energy fractions** (et1, et2, et3, et4): Fraction of cluster energy in the 1st through 4th highest-energy towers. Photons concentrate energy in fewer towers.
- **Cluster widths** (weta_cogx, wphi_cogx): Second moments of the tower energy distribution in eta and phi, weighted by log(E_tower/E_cluster). Narrower for photons.

### Energy Scale Systematic

The beam-test-established resolution and the pi0 calibration procedure set the photon energy scale uncertainty. PPG12 evaluates this through `escale_up`/`escale_down` config variations that shift cluster energies by the calibration uncertainty.

### Energy Resolution Systematic

The gap between beam-test resolution and in-situ performance (quantified by the 11.9% additional smearing) drives the `eres_up`/`eres_down` variations that alter the unfolding response matrix.

### Position-Dependent Response

The ~20% energy variation across a tower from light guide non-uniformity is corrected during cluster reconstruction, but residual non-uniformity contributes to the constant term of the resolution and enters the energy scale systematic.

## Technical Specifications Summary

| Specification | Value | PPG12 Impact |
|---------------|-------|-------------|
| Energy resolution (operational) | 2.8% + 15.5%/sqrt(E) | Unfolding response matrix, pT bin migration |
| Tower size | 0.025 x 0.025 | Shower shape variable granularity |
| Depth | 20 X0 | Full EM shower containment up to 35 GeV |
| Moliere radius | ~2.3 cm | Cluster size, e11/e33 discrimination |
| Hadronic interaction length | 0.83 lambda_I | Small hadronic leakage to IHCal |
| Light yield | 500 p.e./GeV | Statistical contribution to resolution |
| Temperature stability | 0.33%/degC (corrected) | Energy scale stability |
| Additional smearing (data-MC) | 11.9% | Energy resolution systematic |

## References

- 2016 beam test: arXiv:1704.01461 (IEEE TNS 65 (2018) 2901)
- 2018 beam test: arXiv:2003.13685 (IEEE TNS 68 (2021) 173)
- First calorimeter results: arXiv:2504.02242
- Upgrade proposal: arXiv:1501.06197
- sPHENIX TDR: BNL internal document (2019)

## Cross-Links

- [sPHENIX Overview](sphenix-overview.md) -- full detector layout
- [Calorimeter System](calorimeter-system.md) -- combined EMCal + HCal performance
- [Shower Shape Variables](../../concepts/shower-shape-variables.md) -- variable definitions and BDT usage
- [Tree Making](../../pipeline/01-tree-making.md) -- cluster reconstruction from calibrated towers
- [MC Samples](../../reference/mc-samples.md) -- GEANT4 simulation of EMCal response
