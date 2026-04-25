# Efficiency Delta-R Comparison

## Overview

This study compares the factorized photon efficiency chain (reco, isolation, photon ID, total) with and without the `deltaR < 0.1` truth-matching cut, across three MC sample types: single-interaction, double-interaction, and physics-mixed (0 mrad: 77.6% single + 22.4% double, cluster-weighted). All samples use the `photon10` generator (truth pT 14--30 GeV).

**Key finding**: Removing the deltaR cut improves total efficiency for all interaction types. The reconstruction stage dominates the gain. Isolation efficiency is unaffected. Photon ID (BDT) efficiency degrades in double-interaction MC due to vertex-shifted kinematics, but the net effect is positive because the reconstruction gain outweighs the BDT loss.

Quantitatively, total efficiency improves by +2.1% (single), +11.5% (mixed), and nearly doubles for pure double-interaction MC.

## Efficiency Tables

### Absolute Efficiencies (noDR Overlay: Single vs Double vs Mixed)

| Stage | Single | Double | Mixed 0 mrad |
|-------|--------|--------|--------------|
| Reco  | 93.5%  | 73.1%  | 88.4%        |
| Iso   | 74.2%  | 71.0%  | 73.5%        |
| ID    | 82.7%  | 69.4%  | 80.0%        |
| All   | 57.4%  | 36.0%  | 52.0%        |

With the deltaR cut removed, the remaining differences between single and double MC reflect genuine pileup physics effects rather than truth-matching artifacts.

### withDR vs noDR Comparison

| Stage | Single withDR | Single noDR | Double withDR | Double noDR | Mixed withDR | Mixed noDR |
|-------|--------------|-------------|---------------|-------------|-------------|------------|
| Reco  | 91.4%        | 93.5%       | 30.7%         | 73.1%       | 76.2%       | 88.4%      |
| Iso   | 74.2%        | 74.2%       | 71.1%         | 71.0%       | 73.9%       | 73.5%      |
| ID    | 82.9%        | 82.7%       | 82.6%         | 69.4%       | 82.9%       | 80.0%      |
| All   | 56.2%        | 57.4%       | 18.0%         | 36.0%       | 46.6%       | 52.0%      |

### Ratio (noDR / withDR)

| Stage | Single | Double | Mixed 0 mrad |
|-------|--------|--------|--------------|
| Reco  | 1.023  | 2.386  | 1.161        |
| Iso   | 1.000  | 0.998  | 0.995        |
| ID    | 0.998  | 0.840  | 0.966        |
| All   | 1.021  | 1.999  | 1.115        |

## Per-Stage Discussion

### Reconstruction

The reconstruction stage captures the full impact of deltaR removal. Single MC improves by 2.3% (91.4% to 93.5%), reflecting the small population of vertex-displaced clusters that fail the geometric cut. Double MC improves by a factor of ~2.4 (30.7% to 73.1%), recovering the 65% of clusters that failed deltaR due to vertex displacement.

The residual gap between double (73.1%) and single (93.5%) noDR efficiency reflects two genuine effects: (1) cluster loss from merging with pileup energy deposits (no-trkID-match rate: 7.9% double vs 2.4% single, ~5.5 pp gap), and (2) common cut failures from shifted kinematics (~15 pp gap).

### Isolation

Isolation efficiency is unaffected by the deltaR removal (ratios within 0.5% of unity for all samples). Comparing single to double MC, the ~3 pp absolute degradation (74.2% to 71.0%) reflects extra hadronic energy deposited in the isolation cone by the second interaction. This is a genuine, small pileup effect that does not interact with the matching requirement.

### Photon ID (BDT)

This is the most physically interesting stage. Single MC shows negligible change (-0.2 pp, ratio 0.998). Double MC shows a significant 13 pp reduction (82.6% to 69.4%, ratio 0.840). The mechanism is the ET-dependent BDT threshold: `BDT_min(ET) = intercept + slope * ET`. In double-interaction events, the vertex shift modifies the reconstructed ET = E / cosh(eta_reco), causing the cluster to evaluate the threshold at a shifted working point. The shifted shower-shape variables (particularly w_eta, which is sensitive to eta) further degrade the BDT score.

The mixed 0 mrad sample shows a 3.4% decrease (82.9% to 80.0%), consistent with the 22.4% cluster-weighted double-interaction admixture.

### Total Efficiency

The total efficiency improves for all interaction types because the reconstruction gain dominates the BDT degradation:
- Single: 56.2% to 57.4% (+2.1% relative)
- Double: 18.0% to 36.0% (nearly doubling)
- Mixed 0 mrad: 46.6% to 52.0% (+11.5% relative)

For mixed MC, the improvement is essential for pileup efficiency studies: it removes the dominant truth-matching artifact that previously inflated the apparent pileup degradation from ~17% (with deltaR) to ~9% (without deltaR, genuine physics effect only).

## New Clusters BDT Pass Rate

Clusters recovered by removing the deltaR cut ("new clusters": pass trkID, fail deltaR < 0.1) have a systematically lower BDT pass rate due to vertex-shifted kinematics:

| Sample | BDT pass rate (new clusters) | BDT pass rate (all noDR) |
|--------|------------------------------|--------------------------|
| Single | 73.9%                        | 82.7%                    |
| Double | 59.8%                        | 69.4%                    |
| Mixed  | 61.7%                        | 80.0%                    |

The new clusters show a 10--13 pp lower BDT pass rate than the full population across all sample types. In single MC, only ~2.3% of matched clusters are "new," so their lower pass rate has negligible impact on the integrated efficiency. In double MC, ~58% of matched clusters are "new," driving the 13 pp overall BDT efficiency degradation. The mixed sample's new-cluster BDT pass rate is dominated by the double-interaction contribution.

## Key Files

| File | Purpose |
|------|---------|
| `efficiencytool/compare_efficiency_deltaR.C` | Main analysis macro: computes per-stage efficiency for 3 sample sets x 2 matching modes |
| `plotting/plot_efficiency_comparison_deltaR.C` | Plotting macro: 6 canvases (per-stage, ratios, BDT focus, summary table) |
| `efficiencytool/results/efficiency_comparison_deltaR.root` | Output ROOT file with all TEfficiency objects |
| `efficiencytool/reports/efficiency_deltaR_comparison.tex` | Full LaTeX report |
| `plotting/figures/efficiency_comparison/` | Output figures directory |

## See Also

- [Delta-R Truth Matching](deltaR-truth-matching.md) -- The companion study demonstrating the deltaR cut is a vertex-resolution artifact
- [Double-Interaction Efficiency](double-interaction-efficiency.md) -- Per-stage efficiency investigation that first identified the matching artifact
- [Systematic Variations](systematic-variations.md) -- Pileup treated as a systematic uncertainty
