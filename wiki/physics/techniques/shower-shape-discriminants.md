# Shower Shape Discriminants

Electromagnetic shower shape variables are the primary observables used to discriminate prompt photons from hadronic background (pi0, eta, and other meson decays) in calorimeters. This article explains the physics origin of discrimination power for each variable class, how the sPHENIX CEMCal W/SciFi granularity affects them, and how they enter the PPG12 BDT.

## Physics of Photon vs Hadron Showers

A prompt photon deposits its energy in a single, narrow electromagnetic shower centered on the impact point. The dominant background -- pi0 -> gamma gamma -- produces two overlapping showers whose separation depends on the pi0 boost:

```
angular separation ~ 2 / (E_pi0 / m_pi0)
```

At sPHENIX energies (8--35 GeV), the two photons from a pi0 are separated by 1--3 CEMCal tower widths, producing a broader, more elongated, or double-peaked energy deposit compared to a single photon. As energy increases, the two showers merge and become harder to distinguish -- this is why photon purity decreases at high ET and why ET-dependent cuts are necessary.

Other backgrounds (eta -> gamma gamma, jet fragments with leading pi0) produce similar but typically even broader shower profiles due to larger opening angles or multiple overlapping particles.

## The sPHENIX CEMCal

The sPHENIX Central Electromagnetic Calorimeter (CEMCal) is a tungsten-powder/scintillating-fiber (W/SciFi) sampling calorimeter:

- **Tower size:** approximately 0.024 x 0.024 in Delta_eta x Delta_phi
- **Barrel coverage:** |eta| < 1.1
- **Energy resolution:** sigma_E/E ~ 2.8% + 15.5%/sqrt(E) (full-tower operational, E in GeV)
- **Moliere radius:** approximately 2.3 cm (W/SciFi), spanning roughly one tower width

The granularity is coarser than ATLAS liquid-argon (0.003 x 0.1 strips, 0.025 x 0.025 middle) or CMS PbWO4 crystals (0.0174 x 0.0174), but comparable to ALICE EMCal (0.014 x 0.014 lead-scintillator). The W/SciFi technology provides compact, dense showers with good energy resolution but limited transverse segmentation. This makes multivariate methods essential: no single tower-level variable provides enough discrimination on its own.

### 7x7 Tower Matrix

All shower shape variables in PPG12 are computed from a 7x7 tower energy matrix centered on the cluster center-of-gravity (CoG), as implemented in `anatreemaker/source/CaloAna24.cc`. Towers with energy below 0.07 GeV are zeroed. This matrix captures the full transverse shower profile over approximately 3 Moliere radii in each direction.

## Variable Classes

### Tower Energy Ratios

Energy ratios measure how concentrated the shower is in the core relative to the periphery. For a NxM tower window centered on the CoG, `eNM` is the total energy within that window.

#### Core Concentration: e11/e33

```
e11_over_e33 = E(1x1 seed tower) / E(3x3 around CoG)
```

A single photon deposits most energy in the seed tower, giving e11/e33 ~ 0.4--0.6. A pi0 with separated decay photons spreads energy across neighboring towers, giving lower e11/e33. This is the single most important shower shape feature across all experiments -- it appears as R_eta in ATLAS, is related to sigma_eta_eta in CMS, and drives the sigma_long^2 eigenvalue in ALICE.

#### Strip Concentration: e32/e35

```
e32_over_e35 = E(3 eta x 2 phi strip) / E(3 eta x 5 phi)
```

Measures how well the shower is contained in a narrow phi strip. Sensitive to the phi asymmetry of pi0 decays (which can produce an elongated profile in either eta or phi depending on the decay plane orientation).

#### Asymmetry Ratios: e11/e13, e11/e31, e11/e15, e11/e51, etc.

```
e11_over_e13 = E(1x1) / E(1 eta x 3 phi)
e11_over_e31 = E(1x1) / E(3 eta x 1 phi)
```

These probe asymmetries between the eta and phi directions. A single photon shower is symmetric (e11/e13 ~ e11/e31), while a pi0 with one decay photon displaced along eta produces e11/e31 < e11/e13. The full set (e11/e13, e11/e15, e11/e17, e11/e31, e11/e51, e11/e71) provides multi-scale probes of asymmetry.

#### 2x2 Ratios: e22/e33, e22/e35, e22/e37, e22/e53

```
e22_over_e33 = E(2x2 quadrant) / E(3x3)
```

The 2x2 window captures the CoG quadrant. These ratios are sensitive to the shower centroid position within the seed tower and to substructure at the one-tower scale.

### Width Variables

#### weta_cogx and wphi_cogx

Energy-weighted second moment of the tower energy distribution in eta (phi), excluding the seed tower:

```
weta_cogx = sum_{i != seed} [E_i * (eta_i - eta_CoG)^2] / sum_{i != seed} E_i
```

Excluding the seed tower makes these variables more sensitive to the shower tails. A single photon produces a narrow distribution (small weta_cogx), while a pi0 with separated showers produces broader tails (larger weta_cogx).

These are the sPHENIX analogs of:
- ATLAS R_eta, R_phi (energy ratios encoding the same width information)
- CMS sigma_eta_eta (second moment of the 5x5 crystal matrix)
- ALICE sigma_long^2 (eigenvalue combining both eta and phi widths)

The "cogx" (center-of-gravity, excluding seed) variants are the primary BDT features because they isolate the tail shape without being dominated by the seed tower energy.

#### Strip Widths: w32, w52, w72

Second moment in eta for a 2-column-wide strip in phi centered on the cluster:

```
w32: |delta_eta| <= 1 (3 eta towers, 2 adjacent phi)
w52: |delta_eta| <= 2 (5 eta towers, 2 adjacent phi)
w72: |delta_eta| <= 3 (7 eta towers, 2 adjacent phi)
```

These provide multi-scale measurements of the shower width in the eta direction, restricted to the core phi region. The progression from w32 to w72 probes the shower at increasing distances from the core.

### Ring Energy Fractions (et1--et4)

```
et1 = E(4 adjacent towers) / E_total     (innermost ring)
et2 = E(2nd ring) / E_total
et3 = E(3rd ring) / E_total
et4 = E(outermost ring) / E_total
```

The rings are defined around the seed tower in the 7x7 matrix. Together they form a radial energy profile:

- Single photon: steeply falling (et1 >> et2 >> et3 >> et4)
- pi0 with separated photons: flatter profile (more energy in outer rings)
- Broad jet fragment: even flatter

The ring fractions are complementary to the energy ratios: ratios compare rectangular windows, while rings probe circular symmetry. The combination captures both symmetric and asymmetric shower broadening.

### Kinematic and Position Variables

The BDT also uses non-shower-shape features:

| Variable | Role |
|----------|------|
| `cluster_Et` | Shower shape distributions are strongly ET-dependent (pi0 merging) |
| `cluster_Eta` | Shower profile varies with incidence angle and material budget |
| `vertexz` | Vertex position affects the cluster eta/ET reconstruction and thus the shower shape |

Including these allows the BDT to learn ET-dependent and position-dependent discrimination boundaries, rather than requiring separate flat cuts for each kinematic region.

## How Granularity Affects Discrimination

The key physical scale is the ratio of the pi0 decay photon separation to the tower size. At sPHENIX:

- **Well-separated regime (low ET):** Two photons hit different towers; shower shape variables can resolve the double structure. Discrimination is strong.
- **Partially merged regime (medium ET, ~10--20 GeV):** Photons are ~1--2 towers apart; the shower is broader than a single photon but not fully resolved as two peaks. The BDT exploits subtle broadening.
- **Fully merged regime (high ET, >25 GeV):** Photons hit the same tower or adjacent towers; the shower approaches a single-photon shape. Discrimination is weakest. This is the regime where ATLAS's fine strip layer provides a significant advantage.

Compared to ATLAS (0.003 strip pitch, effectively 8x finer than sPHENIX) and CMS (0.017, ~1.4x finer), sPHENIX has less ability to resolve substructure at high ET. The multivariate BDT partially compensates by combining information from many variables that are individually weak discriminants.

Compared to ALICE EMCal (0.014, ~1.7x finer), the sPHENIX granularity is coarser per tower but the overall approach is similar: both rely on energy distribution moments and both use ABCD methods for purity extraction.

## Feature Importance in the BDT

While exact feature importances depend on the ET bin and model variant, the general hierarchy is:

1. **Core energy ratios** (e11/e33, e11/e22): most discriminating at all energies
2. **Ring fractions** (et1, et2): strong discrimination from radial profile
3. **Width variables** (weta_cogx, wphi_cogx): important for asymmetric showers
4. **Extended ratios** (e11/e31, e11/e51, etc.): provide multi-scale information
5. **Kinematic variables** (ET, eta, vertexz): enable position-dependent boundaries

The BDT learns to combine these variables in ways that are not accessible to rectangular cuts. For example, it can impose a tighter e11/e33 cut for showers that also have large weta_cogx, capturing the correlation between core concentration and width.

## 11 Variables in ShowerShapeCheck.C

The shower shape monitoring code tracks distributions for 11 variables across ABCD regions, interaction types (single/double), and kinematic bins:

```
weta_cogx, wphi_cogx, wr_cogx, et1, et2, et3, et4,
e11/e33, e32/e35, bdt_score, npb_score
```

This set covers the primary width, ratio, and ring fraction observables plus the two multivariate scores.

## Key References

| Reference | Relevant shower shape content |
|-----------|-------------------------------|
| ATLAS, arXiv:1012.4389 | Appendix A: full list of shower shape variables, tight cut definitions |
| CMS, arXiv:1108.2044 | sigma_eta_eta definition, shower shape preselection |
| ALICE, arXiv:1906.01371 | sigma_long^2 eigenvalue, narrow/wide cluster classification |
| D0, hep-ex/0511054 | EM3 cluster width, cell counting in EM1 |

## See Also

- [Shower Shape Variables (concept)](../../concepts/shower-shape-variables.md) -- full variable definitions, branch names, computation code
- [BDT-Based Photon Identification](bdt-photon-id.md) -- how these variables enter the BDT training and inference
- [Stage 2: BDT Training](../../pipeline/02-bdt-training.md) -- feature extraction pipeline, 25-feature list, model variants
