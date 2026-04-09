# Calorimeter System

## Overview

The sPHENIX calorimeter system provides hermetic coverage at |eta| < 1.1 with three concentric subsystems: the electromagnetic calorimeter (EMCal), the inner hadronic calorimeter (IHCal), and the outer hadronic calorimeter (OHCal). Together they total approximately 6 nuclear interaction lengths, containing ~97% of hadronic energy below 50 GeV. For PPG12, this system serves two distinct roles: the EMCal provides the photon energy and shower-shape measurement, while the full EMCal + IHCal + OHCal system provides the calorimetric isolation energy sum.

## Subsystem Specifications

### EMCal

See [EMCal Performance](emcal-performance.md) for full details. Summary:

| Parameter | Value |
|-----------|-------|
| Technology | W/SciFi SPACAL |
| Depth | 20 X0, 0.83 lambda_I |
| Tower size | 0.025 x 0.025 |
| Energy resolution (electrons) | 2.8% + 15.5%/sqrt(E) |
| Coverage | \|eta\| < 1.1, full 2pi |

### Inner HCal (IHCal)

The IHCal sits inside the solenoid bore, immediately behind the EMCal:

| Parameter | Value |
|-----------|-------|
| Technology | Aluminum absorber plates + scintillating tiles |
| Inner radius | 116 cm |
| Outer radius | 137 cm |
| Tower size | 0.1 x 0.1 (Delta_eta x Delta_phi) |
| Absorber plates | 21, ASTM A36 steel, 1.02--1.47 cm thick |
| Tilt angle | 32 deg relative to radius |
| Scintillator | 0.7 cm thick extruded polystyrene tiles, embedded 1 mm WLS fiber (Kuraray Y11) |
| Gap width | 0.85 cm |
| Sampling fraction | 0.060--0.078 (radius dependent) |
| Depth | ~1 lambda_I |

The tile orientation is tilted at 32 degrees from the radial direction to reduce channeling, where particles could pass through only scintillator or only absorber. This tilting ensures each particle traverses multiple alternating layers regardless of incidence angle.

### Outer HCal (OHCal)

The OHCal sits outside the solenoid cryostat and doubles as the magnetic flux return:

| Parameter | Value |
|-----------|-------|
| Technology | Steel absorber plates + scintillating tiles |
| Inner radius | 182 cm |
| Outer radius | 269 cm |
| Tower size | 0.1 x 0.1 (Delta_eta x Delta_phi) |
| Absorber plates | 21, ASTM A36 steel, 2.62--4.25 cm thick |
| Tilt angle | 12 deg relative to radius |
| Scintillator | 0.7 cm thick tiles with serpentine WLS fiber |
| Light yield | ~400 p.e./GeV (~12 p.e./MIP) |
| Tile attenuation length | ~2--2.5 m |
| Sampling fraction | 0.028--0.037 (radius dependent) |
| Depth | ~4 lambda_I |
| Punch-through probability | 1% for particles > 2 GeV at normal incidence |

The OHCal tiles are tilted at 12 degrees (opposite direction from IHCal), with plates staggered by half-fin thickness to minimize gaps. The steel plates simultaneously absorb hadronic showers and return the magnetic flux from the 1.4 T solenoid.

### Readout

Both hadronic calorimeters use wavelength-shifting (WLS) fibers embedded in scintillating tiles, read out by SiPMs. The scintillation light from each tile is collected by the WLS fiber and transmitted to the SiPM at the tile edge. The HCal towers (0.1 x 0.1) are coarser than EMCal towers (0.025 x 0.025) by a factor of 4 in each direction.

**HCal calibration** uses cosmic ray muon MIP energy deposition: single-tower triggers select vertically traversing muons, with offline cuts ensuring azimuthal-direction selection. The measured MIP energy is then corrected to the full absorber+scintillator response using single-hadron GEANT4 simulations validated by beam test data. Temperature-dependent corrections account for SiPM gain variation during collision running.

## Combined System Performance

### Hadronic Energy Resolution

Measured at the 2016 Fermilab beam test (arXiv:1704.01461):

| Configuration | Constant term | Stochastic term |
|---------------|---------------|-----------------|
| Combined EMCal+HCal (all showers) | 13.5% | 64.9%/sqrt(E) |
| HCal-only (pass through EMCal) | 14.5% | 74.9%/sqrt(E) |
| OHCal-only (pass through EMCal+IHCal) | 17.1% | 75.5%/sqrt(E) |
| HCal standalone (hadrons) | 11.8% | 81.1%/sqrt(E) |
| Design target (single hadron) | -- | < 100%/sqrt(E) |

The combined resolution (13.5% + 64.9%/sqrt(E)) comfortably meets the design specification of better than 100%/sqrt(E) for single hadrons. The improvement in the combined system over HCal-only comes from measuring the electromagnetic component of hadronic showers in the EMCal.

### Total System Depth

| Subsystem | Interaction lengths (lambda_I) |
|-----------|-------------------------------|
| EMCal | ~0.83 |
| IHCal | ~1.0 |
| Solenoid cryostat | ~0.1 (aluminum equivalent) |
| OHCal | ~4.0 |
| **Total** | **~5.9** |

At ~6 lambda_I total depth, approximately 97% of hadronic energy below 50 GeV is contained, and the punch-through probability for the full system is 0.6% for protons.

### Energy Sharing Between Subsystems

Measured in-situ with Au+Au collision data (arXiv:2504.02242):

- **EMCal sees ~70% of total ET** -- absorbs all electromagnetic energy and a large fraction of hadronic energy
- **HCal sees ~20% of total ET** -- captures hadronic shower tails leaking beyond the EMCal
- **Full calorimeter MC correction factor ~0.9** -- some energy is lost in inactive material (cryostat, support structures)

The EMCal-only and HCal-only dET/deta measurements are independently consistent with previous PHENIX and STAR results, validating the separate calibration chains.

### Projective Tower Geometry

All calorimeter towers are projective in eta from z = 0 (the nominal detector center). Tower positions (eta, phi) are defined relative to z_vtx = 0. For events with |z_vtx| up to 10 cm, the tower ET is recalculated as:

```
ET_tower = E_tower * sin(theta_tower)
```

where theta is computed from the vertex-corrected tower position. This projective assumption introduces vertex-dependent biases in cluster kinematics that are studied in the double-interaction analysis (`DoubleInteractionCheck.C`).

## Jet Reconstruction

The hermetic calorimetry enables full jet reconstruction in sPHENIX. Design performance from GEANT4 simulation:

- **Jet energy resolution (pp):** Substantially better than 120%/sqrt(E) for R = 0.2 and R = 0.4
- **Jet energy scale uncertainty:** < 3%
- **No significant quark vs gluon jet resolution difference** despite different fragmentation functions
- **In central Au+Au:** Jet resolution dominated by underlying event fluctuations (~3.5 GeV RMS in R = 0.2 cone, ~7 GeV in R = 0.4), not intrinsic detector resolution

For PPG12, jet reconstruction is not directly used, but the jet energy scale and resolution are relevant as cross-checks for the calorimetric isolation cone energy.

## Isolation Energy Measurement

PPG12 measures the isolation energy in a cone of R = 0.3 around each photon candidate cluster, using all three calorimeter subsystems. The isolation energy sum E_T^iso is the total transverse energy in calorimeter towers within the cone, excluding the photon cluster itself.

### Tower Contributions to Isolation Cone

| Subsystem | Tower size | Towers in R = 0.3 cone | Energy contribution |
|-----------|-----------|----------------------|-------------------|
| EMCal | 0.025 x 0.025 | ~113 full towers | EM energy + some hadronic |
| IHCal | 0.1 x 0.1 | ~7 full towers | Hadronic shower tails |
| OHCal | 0.1 x 0.1 | ~7 full towers | Hadronic shower remnants |

### Isolation Performance Considerations

1. **EMCal dominates the isolation sum.** Since photon showers are fully contained in the EMCal, nearby pi0 decays and hadronic activity deposit most of their electromagnetic energy in EMCal towers within the isolation cone.

2. **HCal contribution matters for hadronic isolation.** Charged hadrons (pions, kaons, protons) from the underlying event or nearby jet fragments deposit energy across all three subsystems. Including IHCal and OHCal towers in the cone captures this hadronic component.

3. **Underlying event energy.** In pp collisions at 200 GeV, the mean UE energy per EMCal tower is small (~53 MeV extrapolated from Au+Au; much less in pp). The parametric isolation cut (iso_max = b + s * ET) absorbs the UE contribution through the offset parameter b.

4. **Zero suppression.** Hardware zero suppression at 2-sigma pedestal noise removes noise towers from the isolation sum but introduces a systematic (1--6% for EMCal, smaller for HCal) from discarding low-energy real energy deposits.

5. **Noise contamination.** Residual noise towers above the zero-suppression threshold can bias the isolation energy upward. This is mitigated by the offline peak-minus-pedestal algorithm with subsystem-dependent thresholds.

## Systematic Uncertainties from Calorimetry

PPG12 calorimetry systematics (from arXiv:2504.02242, Table 2):

| Source | EMCal | HCal | Combined |
|--------|-------|------|----------|
| Calibration | 2.6% | 2.7% | 2.1% |
| Hadronic response | 4.1% | 6.6% | 4.7% |
| Modeling | 1.4--2.1% | 2.5--3.8% | 1.6--2.2% |
| Zero suppression | 1.0--5.8% | 0.2--0.3% | 0.8--4.4% |

These feed into PPG12 through the energy scale (`escale_up`/`escale_down`) and energy resolution (`eres_up`/`eres_down`) systematic variations.

## References

- 2016 beam test: arXiv:1704.01461
- 2018 beam test: arXiv:2003.13685
- First calorimeter results (dET/deta): arXiv:2504.02242
- First sPHENIX highlights: arXiv:2410.18031
- Upgrade proposal: arXiv:1501.06197

## Cross-Links

- [EMCal Performance](emcal-performance.md) -- detailed EMCal specifications
- [sPHENIX Overview](sphenix-overview.md) -- full detector layout
- [ABCD Method](../../concepts/abcd-method.md) -- how isolation defines signal region
- [Isolation Cuts](../../concepts/isolation-cuts.md) -- parametric isolation, cone definition
- [Tree Making](../../pipeline/01-tree-making.md) -- tower and cluster readout
- [MC Samples](../../reference/mc-samples.md) -- GEANT4 calorimeter simulation
