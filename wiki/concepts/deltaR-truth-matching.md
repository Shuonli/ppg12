# Delta-R Truth Matching Study

## Overview

The PPG12 efficiency calculation previously used two-step truth matching for
associating reconstructed EMCal clusters with generator-level photons:

1. **GEANT track ID**: `cluster_truthtrkID` matched to `particle_trkid` (ground truth from simulation)
2. **deltaR < 0.1**: geometric cut between cluster position (at reco vertex) and truth particle direction (at true vertex)

A study of single-interaction MC showed the deltaR cut rejects approximately
3.5% of well-reconstructed photons. The rejection is driven by vertex
resolution -- the mismatch between reconstructed and true vertex positions --
not by genuine reconstruction failure or misidentification.

**Decision**: The deltaR cut was removed. Track ID matching alone is now used
for truth association in the production efficiency calculation.

## Truth Matching in the Pipeline

**Where it is used**: The primary consumer is `RecoEffCalculator_TTreeReader.C`,
which computes reconstruction efficiency and builds the response matrix.
Approximately 15 additional diagnostic macros read `eff_dR` from the analysis
config for optional cross-checks.

**Track ID matching**: The GEANT simulation propagates a unique track ID through
the detector response. The slimtree stores `cluster_truthtrkID` (per cluster)
and `particle_trkid` (per truth particle). Matching these gives an unambiguous
cluster-to-truth association with no geometric approximation.

**deltaR computation**: The geometric distance was computed as

    deltaR = sqrt(delta_eta^2 + delta_phi^2)

between the cluster centroid (evaluated at the reconstructed vertex) and the
truth particle four-vector (evaluated at the true vertex). This introduces a
systematic bias whenever reco and true vertices differ.

**Config field**: `eff_dR: 0.1` was set in all analysis YAML configs under
`analysis.eff_dR`. The field is retained for use by diagnostic macros but is no
longer applied in the production efficiency calculation.

**Note**: `BDTinput.C` loads `eff_dR` from the config but never applies it --
training labels use track ID matching only. This was already the case before the
study.

## Vertex Mismatch Mechanism

The cluster eta is computed from the cluster z-position on the EMCal surface and
the reconstructed vertex. When the reconstructed vertex differs from the true
vertex by dz, the cluster eta shifts:

    delta_eta ~ -dz / (R_CEMC * cosh(eta))

where `R_CEMC = 93.5 cm` is the EMCal cylinder radius.

For the deltaR < 0.1 cut to fail from vertex mismatch alone (at delta_phi ~ 0),
the critical vertex displacement is:

    |dz|_crit = 0.1 * R_CEMC * cosh(eta) = 9.35 cm   (at eta = 0)

**Single-interaction MC vertex resolution**: The cross-section-weighted dz RMS
is approximately 2.9 cm, well below the critical threshold. However, 2.12% of
photons have |dz| > 9.35 cm, landing in the tail where the deltaR cut fails.

**Step-function behavior**: The failure probability is effectively a step
function in |dz|: 0% below |dz| ~ 9 cm, rising sharply to 100% above
|dz| ~ 15 cm. This means the cut does not degrade gracefully -- it either passes
or fails based almost entirely on vertex displacement.

## Key Results (Single-Interaction MC)

### Failure rates per signal sample

| Sample   | No trkID match | Fail dR < 0.1 (of trkID-matched) |
|----------|----------------|-----------------------------------|
| photon5  | 3.42%          | 2.70%                             |
| photon10 | 2.37%          | 2.98%                             |
| photon20 | 1.72%          | 4.01%                             |

The deltaR failure rate increases with pT because higher-pT photons are more
forward-boosted, where cosh(eta) is larger and the critical dz threshold is
smaller.

### Energy quality of failing clusters

73% of clusters that pass track ID matching but fail deltaR < 0.1 have
ET/pT(truth) in [0.8, 1.2]. These are well-reconstructed photons with good
energy response -- they fail only because the vertex displacement pushes the
geometric distance above threshold.

### Single vs double interaction comparison

| Condition         | Fail dR < 0.1 | dz RMS  |
|-------------------|----------------|---------|
| photon10 single   | 3.06%          | 9.8 cm  |
| photon10 double   | 65.0%          | 46.2 cm |

In double-interaction MC, the vertex displacement is much larger (the
reconstructed vertex is the average of two collision vertices), causing the
deltaR cut to reject the majority of truth-matched clusters. This confirmed
that the cut tests vertex resolution, not calorimeter performance.

## Decision: Remove the deltaR Cut

Five arguments supported removing the cut:

1. **Track ID is ground truth** -- the GEANT track ID provides an unambiguous
   truth association. The deltaR cut adds no additional truth content beyond
   what trkID already provides.

2. **Cut tests vertex resolution, not calorimetry** -- the failure mechanism is
   entirely driven by dz, not by shower reconstruction quality.

3. **73% of rejected clusters are good photons** -- well-reconstructed clusters
   with ET/pT in [0.8, 1.2] are discarded, biasing the efficiency downward.

4. **Response matrix captures vertex effects** -- kinematic shifts from vertex
   mismatch are already encoded in the unfolding response matrix. Removing
   clusters before they enter the matrix discards information.

5. **Removes ~3% vertex-dependent efficiency bias** -- the bias depends on the
   dz distribution, which differs between data and MC, introducing a systematic
   that is difficult to control.

The production code (`RecoEffCalculator_TTreeReader.C`) now uses track ID
matching only. The `eff_dR` config field is retained so that diagnostic macros
can still apply the cut for cross-checks.

## Key Files

| File | Purpose |
|------|---------|
| `efficiencytool/study_deltaR_single.C` | Main study macro (all 3 signal MC samples) |
| `plotting/plot_deltaR_single_study.C` | Plotting macro (8 diagnostic figures) |
| `efficiencytool/run_deltaR_study.sh` | Runner script |
| `efficiencytool/results/deltaR_single_study.root` | Study output histograms |
| `efficiencytool/reports/deltaR_matching_study.tex` | Full LaTeX report |
| `plotting/figures/deltaR_study/` | Study figures (8 PDFs) |
| `efficiencytool/RecoEffCalculator_TTreeReader.C` | Production efficiency code (deltaR cut removed) |

## See Also

- [Double-Interaction Efficiency](double-interaction-efficiency.md) -- The companion study for double-interaction MC
- [Config Schema](../reference/config-schema.md) -- `eff_dR` field documentation
- [Efficiency and Yield](../pipeline/04-efficiency-yield.md) -- Where truth matching is used in the pipeline
