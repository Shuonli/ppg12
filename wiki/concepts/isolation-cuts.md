# Isolation Cuts

## Overview

Isolation energy measures the activity around a photon candidate within a cone. Genuine prompt photons are isolated (low cone energy), while jet fragments have nearby hadronic activity.

## Reco Isolation: Three Methods

### 1. Standard Cone Isolation (tower-level)

Computed by `CaloAna24.cc::calculateET()`. Sums tower energies within deltaR of the cluster, separately for EMCal, IHCal, OHCal. The cluster's own ET is subtracted from the EMCal sum.

Branches: `cluster_iso_03_{emcal,hcalin,hcalout}_{node}`

With tower energy threshold: `cluster_iso_03_60_{emcal,hcalin,hcalout}_{node}` (towers > 0.06 GeV)

Total reco isolation (config `iso_threshold=1`):
```
recoisoET = cluster_iso_03_60_emcal + cluster_iso_03_60_hcalin + cluster_iso_03_60_hcalout
```

### 2. RawCluster Isolation (module-level)

Computed by the `ClusterIso` Fun4All module at fixed cone sizes.

Branches: `cluster_iso_02`, `cluster_iso_03`, `cluster_iso_04`

### 3. Topological Cluster Isolation

Uses topological clusters instead of individual towers. More robust against noise.

Branches: `cluster_iso_topo_03`, `cluster_iso_topo_04`

Selected via config:
- `use_topo_iso: 0` -- standard cone R=0.3 (legacy)
- `use_topo_iso: 1` -- topo R=0.3
- `use_topo_iso: 2` -- topo R=0.4 (**current nominal**)

## Parametric Isolation Cut

The isolation threshold is ET-dependent (not flat):

```
reco_iso_max = reco_iso_max_b + reco_iso_max_s * ET
```

| Config set | Intercept (b) | Slope (s) | Used in |
|------------|---------------|-----------|---------|
| Nominal BDT | 0.502095 | 0.0433036 | config_bdt_nom.yaml |
| Showershape | 0.453194 | 0.0360234 | config_showershape_*.yaml |
| FunWithxgboost | 1.08128 | 0.0299107 | config_nom.yaml (older) |

## Region Boundaries

- **Isolated:** `reco_iso_min < isoET < reco_iso_max`
  - Default: `reco_iso_min = -20.0`
- **Non-isolated:** `reco_iso_max + reco_noniso_min_shift < isoET < reco_noniso_max`
  - Default: shift = 0.8, max = 20.0
  - The gap prevents the iso/noniso boundary from being sharp

## MC Isolation Correction

Simulation isolation is adjusted to match data:
```cpp
recoisoET = recoisoET * mc_iso_scale;   // code default: 1.0, nominal config: 1.2
recoisoET += mc_iso_shift;              // code default: 0.0, nominal config: 0.1
```

## Truth Isolation

Computed in `CaloAna24.cc` by summing truth particle energies within a cone around the photon. No detector effects.

Branches: `particle_truth_iso_02`, `particle_truth_iso_03`, `particle_truth_iso_04`

**Truth isolation cut:** `iso_ET_truth < 4 GeV` (R=0.3)

This cut is applied when selecting truth photons for efficiency calculation. The cone size depends on `cone_size` config (3 = R=0.3 for truth).

## Isolation Variables Available in Slimtree

| Branch pattern | R | Threshold | Type |
|---------------|---|-----------|------|
| `cluster_iso_{02,03,04}` | 0.2/0.3/0.4 | none | RawCluster module |
| `cluster_iso_03_{emcal,hcalin,hcalout}` | 0.3 | none | Tower-level per calo |
| `cluster_iso_03_60_{emcal,hcalin,hcalout}` | 0.3 | 0.06 GeV | Tower-level with threshold |
| `cluster_iso_03_70_{emcal,hcalin,hcalout}` | 0.3 | 0.07 GeV | Tower-level with threshold |
| `cluster_iso_03_120_{emcal,hcalin,hcalout}` | 0.3 | 0.12 GeV | Tower-level with threshold |
| `cluster_iso_04_{emcal,hcalin,hcalout}` | 0.4 | none | Tower-level per calo |
| `cluster_iso_{005,01,02}_70_emcal` | 0.05/0.1/0.2 | 0.07 GeV | Small-cone EMCal |
| `cluster_iso_03_sub1_{emcal,hcalin,hcalout}` | 0.3 | none | Background-subtracted |
| `cluster_iso_04_sub1_{emcal,hcalin,hcalout}` | 0.4 | none | Background-subtracted |
| `cluster_iso_topo_{03,04}` | 0.3/0.4 | n/a | Topological cluster |
| `cluster_iso_topo_soft_{03,04}` | 0.3/0.4 | n/a | Soft topo cluster |

## Config Fields

| Field | Default | Purpose |
|-------|---------|---------|
| `analysis.cone_size` | 3 | Cone size integer (3 = R=0.3) |
| `analysis.use_topo_iso` | 2 | 0=tower, 1=topo R=0.3, 2=topo R=0.4 |
| `analysis.iso_threshold` | 1 | 1=use threshold-corrected tower iso |
| `analysis.iso_hcalonly` | 0 | 1=HCal-only isolation |
| `analysis.iso_emcalinnerr` | 0.0 | Inner EMCal ring radius (0=full cone) |
| `analysis.reco_iso_max_b` | 0.502095 | Parametric iso intercept |
| `analysis.reco_iso_max_s` | 0.0433036 | Parametric iso slope |
| `analysis.reco_iso_min` | -20.0 | Min reco iso ET |
| `analysis.reco_noniso_min_shift` | 0.8 | Gap between iso and noniso |
| `analysis.reco_noniso_max` | 20.0 | Max noniso ET |
| `analysis.mc_iso_scale` | 1.0 (code) / 1.2 (config) | MC iso multiplicative correction |
| `analysis.mc_iso_shift` | 0.0 (code) / 0.1 (config) | MC iso additive correction |
| `analysis.truth_iso_max` | 4.0 | Truth isolation cut (GeV) |

## Known Discrepancy

The plotcommon.h legend string says `R=0.3` and `< 4 GeV`, but the nominal analysis uses topo R=0.4 with a parametric cut. The legend text describes truth-level isolation, which could be misleading on reco-level plots. See [Constants Sync](../reference/constants-sync.md).
