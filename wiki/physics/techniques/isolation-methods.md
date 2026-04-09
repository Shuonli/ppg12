# Experimental Isolation Methods Across Experiments

## Overview

Every isolated photon measurement faces the same problem: define an isolation criterion that (a) efficiently rejects hadronic background, (b) retains high signal efficiency, (c) is stable against pileup and underlying event (UE), and (d) matches a theoretically calculable definition for NLO/NNLO comparison. Different experiments make different trade-offs depending on detector granularity, tracking coverage, and collision environment.

This article compares isolation implementations across ATLAS, CMS, ALICE, PHENIX, and sPHENIX/PPG12, focusing on how each handles the cone definition, energy threshold, UE subtraction, and residual fragmentation contamination.

## ATLAS: Topological Cluster Isolation with ET-Dependent Threshold

**References:** 1012.4389 (7 TeV, 880 nb^-1), 1908.02746 (13 TeV, 36 fb^-1)

### Implementation

- **Cone radius:** R = 0.4 in the eta-phi plane
- **Energy sum:** Topological clusters of calorimeter cells (both EM and hadronic calorimeters)
- **Core exclusion:** A rectangular region Delta_eta x Delta_phi = 0.125 x 0.175 around the photon barycentre is excluded from the sum to remove the photon's own energy deposition
- **Photon leakage correction:** A pT-dependent correction accounts for photon shower energy spilling outside the excluded core but inside the isolation cone
- **UE/pileup correction:** Event-by-event jet-area method (rho * A_cone), with typical corrections of 3.5 GeV at |eta| < 1.37, decreasing to 1.3 GeV at larger |eta|
- **Data-driven MC correction:** E_T^iso peak position in MC is shifted to match data

### Threshold

The ATLAS parametric isolation cut (used consistently across the 8 TeV and both 13 TeV measurements):

```
E_T^iso < 4.8 + 0.0042 * E_T^gamma   [GeV]
```

This ET-dependent form maintains roughly constant signal efficiency (~95%) across the full ET range (25-2500 GeV at 13 TeV). The slope of 0.0042 accounts for the growth of the UE pedestal with harder interactions.

At the early 7 TeV measurement (1012.4389), a flat threshold was used:
- Isolated: E_T^iso < 3 GeV
- Non-isolated control: E_T^iso > 5 GeV

### Particle-level (unfolded) isolation

For the cross-section definition, ATLAS specifies a particle-level isolation: the sum of transverse energy of all stable particles (lifetime > 10 ps, excluding muons and neutrinos) in a cone R = 0.4, with the same parametric threshold. The unfolding corrects detector-level to this definition.

### Fragmentation survival

With R = 0.4 and epsilon_h ~ 0.1-0.15, approximately 10-20% of the isolated photon cross section at midrapidity comes from the fragmentation mechanism at NLO (from Ichou & d'Enterria, 1005.4529). The residual fragmentation photons that survive are those at very high z (> 0.85-0.9), carrying most of the parent parton momentum.

## CMS: Particle-Flow Isolation

**Reference:** 1108.2044 (7 TeV, 36 pb^-1), 1807.00782 (13 TeV)

### Implementation

- **Cone radius:** R = 0.4
- **Energy sum:** Particle-flow (PF) reconstruction combining tracks, ECAL, HCAL, and muon system. Three isolation components computed separately:
  - Iso_TRK: sum of pT of charged hadrons from primary vertex
  - Iso_ECAL: sum of ET of ECAL PF candidates (excluding the photon cluster)
  - Iso_HCAL: sum of ET of HCAL PF candidates
- **Total isolation:** ISO = Iso_TRK + Iso_ECAL + Iso_HCAL
- **Signal region:** ISO < 5 GeV (7 TeV); variable by analysis

CMS's particle-flow approach provides excellent charged/neutral separation and pileup mitigation (charged hadrons from pileup vertices are explicitly removed). This is not available at sPHENIX, which lacks a magnetic tracker of comparable capability for this purpose.

### Purity extraction

CMS uses template fitting of the ISO distribution:
- Signal template: from Z -> e+e- electrons or MC photons (peaks sharply at ISO ~ 0)
- Background template: from shower-shape sideband (sigma_eta_eta > 0.011 for barrel)
- Unbinned extended maximum likelihood fit

This is a one-dimensional approach, unlike the 2D ABCD method used by ATLAS, ALICE, and PPG12. CMS also developed a BDT-based photon ID for the 13 TeV measurement, combining shower shape and isolation into a single discriminant.

### Generator-level isolation

CMS defines generator-level isolation as: sum of ET of all stable particles in R = 0.4 around the photon, threshold < 5 GeV. This is close to the ATLAS definition but with a flat rather than ET-dependent threshold.

## ALICE: Charged-Track Isolation

**References:** 1906.01371 (7 TeV, 473 nb^-1), 2407.01165 (13 TeV, 10.79 pb^-1)

### Evolution from 7 TeV to 13 TeV

ALICE's isolation strategy changed significantly between the two measurements:

**7 TeV (1906.01371):**
- **Cone radius:** R = 0.4
- **Energy sum:** Neutral EMCal clusters + charged TPC/ITS tracks in the cone
- **Threshold:** p_T^iso < 2 GeV/c
- **Eta acceptance:** Restricted to |eta| < 0.27 (full isolation cone must fit within both EMCal and tracker)

**13 TeV (2407.01165):**
- **Cone radius:** R = 0.4
- **Energy sum:** Charged tracks only (pT > 0.15 GeV/c, |eta| < 0.9)
- **Threshold:** p_T^iso,ch < 1.5 GeV/c
- **Eta acceptance:** Extended to |eta| < 0.67 (only tracker containment required)

The switch to tracks-only isolation was motivated by the dramatic improvement in eta acceptance (from 0.27 to 0.67). The lower threshold (1.5 vs 2 GeV/c) compensates for the missing neutral component. For clusters near the eta boundary (0.5 < |eta| < 0.67), the isolation momentum is scaled up by the fraction of cone area outside the tracker acceptance.

### UE treatment

ALICE does not subtract the underlying event from the isolation momentum at the detector level. Instead, the effect is corrected at the cross-section level through the kappa_iso factor:

```
kappa_iso = (fraction of generated isolated photons that are also isolated after UE)
```

This factor divides the measured cross section, effectively correcting for UE-induced isolation failures. Typical values are kappa_iso ~ 0.85-0.90.

### Fragmentation and purity

ALICE's low-pT reach (down to 7 GeV/c at 13 TeV) means the measurement operates in a regime of heavy contamination: only 5% purity at 7 GeV/c, rising to 50% at 18 GeV/c and ~80% above 80 GeV/c. The ABCD method is essential for extracting the signal from this large background. The isolation criterion is less effective at suppressing fragmentation at low pT because the fragmentation contribution peaks at low z (soft photon), and the absolute isolation threshold becomes less discriminating relative to the UE level.

## PHENIX: Fixed Cone with Proportional Threshold

**Reference:** 1205.5533 (pp 200 GeV, 8 pb^-1)

### Implementation

- **Cone radius:** R = 0.5 (larger than all other experiments discussed here)
- **Energy sum:** EMCal clusters (E > 0.15 GeV) + charged tracks (0.2 < p < 15 GeV/c) within the cone
- **Threshold:** Proportional -- cone energy < 10% of photon energy (epsilon_h = 0.10)
- **Detector coverage:** Two arms, each covering Delta_phi = pi/2 and |eta| < 0.35; fiducial |eta| < 0.25

### Background subtraction

PHENIX does not use the ABCD method. Instead, it employs statistical pi0-tagging:

```
N_dir = N_incl - (1 + A)(1 + R) N_pi0
```

where N_pi0 counts photons paired with another cluster to form a pi0 invariant mass, R corrects for untagged pi0 partners, and A accounts for non-pi0 decay photons (eta, omega, eta'). The isolation cut is applied to the same formula:

```
N_iso_dir = N_iso_incl - (n_iso_pi0 + N_iso_pi0 * R) - A_iso * (1+R) * N_iso_pi0
```

### Isolated/inclusive ratio

PHENIX finds that the ratio of isolated to inclusive direct photons exceeds 90% for pT > 10 GeV/c, consistent with NLO pQCD. The isolation cut rejects approximately 60% of hadronic decay background while retaining nearly all direct photons.

### Fragmentation survival

With R = 0.5 and epsilon_h = 0.10 (z_c = 0.91), fragmentation photons at z < 0.91 are rejected. The larger cone radius compared to PPG12 provides stronger fragmentation suppression but loses more signal at the margins.

## sPHENIX / PPG12: Parametric Cone Isolation

### Implementation

- **Cone radius:** R = 0.3 (reco); R = 0.3 (truth)
- **Energy sum:** Three options, selected by config:
  - Tower-level: sum of EMCal + IHCal + OHCal tower ET within R, with cluster self-subtraction (towers above 0.06 GeV threshold)
  - Topological cluster R = 0.3
  - **Topological cluster R = 0.4** (current nominal, `use_topo_iso: 2`)
- **Threshold:** ET-dependent parametric cut:
  ```
  reco_iso_max = b + s * ET
  ```
  Nominal: b = 0.502 GeV, s = 0.0433
- **MC correction:** Simulation isolation is adjusted before cut application:
  ```
  recoisoET_sim = recoisoET * mc_iso_scale + mc_iso_shift
  ```
  Nominal: scale = 1.2, shift = 0.1 GeV
- **Region boundaries:** A gap of 0.8 GeV separates the isolated and non-isolated regions (config: `reco_noniso_min_shift = 0.8`) to prevent sharp boundary effects in the ABCD method

### Truth-level isolation

The cross section is defined with truth isolation: iso_ET_truth < 4 GeV in a cone R = 0.3 around the truth photon. This is a loose cut compared to the reco-level parametric threshold, ensuring that the unfolding correction for isolation efficiency is moderate.

### Comparison table

| Experiment | R | Threshold | Energy sum | UE subtraction | Eta coverage |
|------------|---|-----------|------------|----------------|-------------|
| ATLAS | 0.4 | 4.8 + 0.0042*ET GeV | Topo-clusters (EM+HAD) | Jet-area method | |eta| < 2.37 |
| CMS | 0.4 | ~5 GeV (flat) | Particle-flow | PV charged + rho correction | |eta| < 2.5 |
| ALICE (13 TeV) | 0.4 | 1.5 GeV (tracks only) | Charged tracks | None (kappa_iso correction) | |eta| < 0.67 |
| PHENIX | 0.5 | 10% of E_gamma | EMCal + tracks | None | |eta| < 0.25 |
| **PPG12** | **0.3** | **0.50 + 0.043*ET GeV** | **Topo-clusters (EM+HAD)** | **MC scale+shift** | **|eta| < 0.7** |

### What fraction of fragmentation survives?

The fraction of the fragmentation cross section surviving isolation depends on the effective z_c = 1/(1 + epsilon_h). At PPG12's nominal cut (ET = 20 GeV: E_T^max ~ 1.37 GeV, epsilon_h ~ 0.069):

- **z_c ~ 0.94** -- only fragmentation at z > 0.94 contributes
- At NLO (Ichou & d'Enterria, JETPHOX with R = 0.3), fragmentation is roughly 10-15% of the total isolated cross section at midrapidity
- This is slightly lower than ATLAS/CMS (z_c ~ 0.87-0.91 with epsilon_h ~ 0.10-0.15) due to PPG12's tighter effective threshold
- PHENIX's proportional threshold gives z_c = 0.91, but the larger R = 0.5 cone captures more hadronic activity, providing additional suppression

The residual fragmentation is fully included in the JETPHOX NLO calculation used for theory comparison.

### PPG12 design choices: rationale

1. **R = 0.3 rather than R = 0.4:** sPHENIX's EMCal granularity (Delta_eta x Delta_phi ~ 0.024 x 0.024) and the desire to minimize contamination from neighboring clusters at relatively low pT (8-35 GeV) favor a smaller cone. The trade-off is borderline NLO reliability at R = 0.3 (Catani et al. warn about log(R) corrections), but the actual reco isolation uses topo R = 0.4 while truth uses R = 0.3.

2. **Parametric rather than flat threshold:** At low ET (8-10 GeV), UE and pileup contribute roughly 0.5-1.0 GeV to the cone, while at high ET (30-35 GeV) the UE pedestal grows. A flat threshold would either cut too aggressively at low ET (losing signal) or too loosely at high ET (admitting background). The linear parametrization b + s*ET balances efficiency and purity across the full range.

3. **MC iso correction (scale + shift):** Simulation does not perfectly reproduce the isolation distribution in data (different UE, noise, tower thresholds). The multiplicative scale (1.2) and additive shift (0.1 GeV) are calibrated to match the isolation peak and tail in data, ensuring that the ABCD region populations are correctly modeled. This is analogous to ATLAS's data-driven peak matching but implemented as a simple two-parameter correction.

4. **Gap between iso and noniso regions:** The 0.8 GeV gap between the isolated and non-isolated boundaries (iso: E_T^iso < reco_iso_max; noniso: E_T^iso > reco_iso_max + 0.8 GeV) prevents events near the boundary from migrating between ABCD regions due to resolution effects. ATLAS uses a similar gap (e.g., iso < 3 GeV, noniso > 5 GeV in the early 7 TeV analysis).

## Systematic Variations of Isolation in PPG12

PPG12 generates systematic variations of the isolation parameters via `make_bdt_variations.py`:

- **iso_tight** / **iso_loose**: intercept b varied by +/- ~20%, slope s varied correspondingly
- **noniso_gap_tight** / **noniso_gap_loose**: non-isolated region boundary shift varied
- **mc_iso_scale** / **mc_iso_shift** variations: MC correction parameters varied

These are propagated through the full ABCD + unfolding + cross-section pipeline and aggregated in quadrature by `calc_syst_bdt.py` under the "purity" systematic group.

## See Also

- [Isolation Theory](../theory/isolation-theory.md) -- Frixione smooth cone, IR safety, fragmentation elimination
- [Isolation Cuts](../../concepts/isolation-cuts.md) -- PPG12-specific implementation details, config fields, branch names
- [ABCD Sideband Method](abcd-sideband-method.md) -- How isolation serves as one axis of the ABCD background subtraction
