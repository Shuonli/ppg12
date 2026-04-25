# Converted vs Unconverted Photon Fractions (PPG12)

Sample: single-particle photon10 MC, `photon10_with_bdt_vtx_noreweight_single.root`, tree `slimtree`. Processed 9,998,541 events.

`particle_converted` convention from CaloAna24: 0 = clean photon, 1 = converted (e+/e- daughters), 2 = bad (non-ee secondary with >40% photon momentum).

Cluster match = at least one reco cluster in `CLUSTERINFO_CEMC` with `cluster_truthtrkID == particle_trkid` and `cluster_Et > 5` GeV.

## Inclusive (|eta|<0.7, 14 <= pT < 30 GeV)

Total truth photons in this selection: **855,380**

| Category | N | Fraction |
|---|---:|---:|
| unconverted | 707,047 | 82.659% +- 0.041% |
| converted (e+e-) | 148,129 | 17.317% +- 0.041% |
| bad secondary | 204 | 0.024% +- 0.002% |

### Cluster-match efficiency

| Category | Matched | Total | Efficiency |
|---|---:|---:|---:|
| unconverted | 481,943 | 708,879 | 67.987% +- 0.055% |
| converted (e+e-) | 98,434 | 148,490 | 66.290% +- 0.123% |
| bad secondary | 80 | 205 | 39.024% +- 3.407% |

### Per-pT-bin fractions (|eta|<0.7)

| pT [GeV] | N | clean | converted | bad |
|---|---:|---:|---:|---:|
|    8-  10 |   611720 | 82.615% | 17.355% | 0.030% |
|   10-  12 |  2430660 | 82.569% | 17.406% | 0.025% |
|   12-  14 |  1022054 | 82.630% | 17.344% | 0.025% |
|   14-  16 |   446977 | 82.638% | 17.334% | 0.027% |
|   16-  18 |   206018 | 82.663% | 17.315% | 0.022% |
|   18-  20 |    99821 | 82.749% | 17.233% | 0.018% |
|   20-  22 |    50681 | 82.530% | 17.456% | 0.014% |
|   22-  24 |    26160 | 82.898% | 17.091% | 0.011% |
|   24-  26 |    13977 | 82.829% | 17.135% | 0.036% |
|   26-  28 |     7564 | 82.113% | 17.861% | 0.026% |
|   28-  32 |     6376 | 83.156% | 16.813% | 0.031% |
|   32-  36 |     2015 | 82.184% | 17.816% | 0.000% |

### Per-|eta|-bin fractions (14 <= pT < 30 GeV, |eta|<0.7)

| |eta| | N | clean | converted | bad |
|---|---:|---:|---:|---:|
| 0.00-0.10 |   125675 | 84.612% | 15.369% | 0.019% |
| 0.10-0.20 |   125600 | 84.251% | 15.729% | 0.020% |
| 0.20-0.30 |   125254 | 83.741% | 16.237% | 0.022% |
| 0.30-0.40 |   123178 | 82.858% | 17.117% | 0.025% |
| 0.40-0.50 |   121160 | 82.102% | 17.871% | 0.026% |
| 0.50-0.60 |   118921 | 80.865% | 19.107% | 0.028% |
| 0.60-0.70 |   115592 | 79.849% | 20.124% | 0.027% |
| 0.70-0.80 |        0 | 0.000% | 0.000% | 0.000% |
| 0.80-0.90 |        0 | 0.000% | 0.000% | 0.000% |
| 0.90-1.00 |        0 | 0.000% | 0.000% | 0.000% |
| 1.00-1.10 |        0 | 0.000% | 0.000% | 0.000% |
| 1.10-1.20 |        0 | 0.000% | 0.000% | 0.000% |
| 1.20-1.30 |        0 | 0.000% | 0.000% | 0.000% |
| 1.30-1.40 |        0 | 0.000% | 0.000% | 0.000% |
| 1.40-1.50 |        0 | 0.000% | 0.000% | 0.000% |

## Inclusive (|eta|<1.5, 14 <= pT < 30 GeV)

Total truth photons in this selection: **1,510,683**

| Category | N | Fraction |
|---|---:|---:|
| unconverted | 1,193,005 | 78.971% +- 0.033% |
| converted (e+e-) | 317,252 | 21.001% +- 0.033% |
| bad secondary | 426 | 0.028% +- 0.001% |

### Cluster-match efficiency

| Category | Matched | Total | Efficiency |
|---|---:|---:|---:|
| unconverted | 676,650 | 1,195,400 | 56.604% +- 0.045% |
| converted (e+e-) | 149,153 | 317,780 | 46.936% +- 0.089% |
| bad secondary | 110 | 427 | 25.761% +- 2.116% |

### Per-pT-bin fractions (|eta|<1.5)

| pT [GeV] | N | clean | converted | bad |
|---|---:|---:|---:|---:|
|    8-  10 |  1199042 | 78.351% | 21.617% | 0.032% |
|   10-  12 |  4873769 | 78.152% | 21.819% | 0.029% |
|   12-  14 |  1968269 | 78.395% | 21.577% | 0.028% |
|   14-  16 |   822918 | 78.683% | 21.285% | 0.031% |
|   16-  18 |   362443 | 79.017% | 20.959% | 0.024% |
|   18-  20 |   167696 | 79.337% | 20.642% | 0.021% |
|   20-  22 |    80901 | 79.665% | 20.306% | 0.028% |
|   22-  24 |    39891 | 80.133% | 19.847% | 0.020% |
|   24-  26 |    20455 | 80.308% | 19.653% | 0.039% |
|   26-  28 |    10706 | 80.282% | 19.681% | 0.037% |
|   28-  32 |     8597 | 81.738% | 18.239% | 0.023% |
|   32-  36 |     2503 | 81.143% | 18.857% | 0.000% |

### Per-|eta|-bin fractions (14 <= pT < 30 GeV, |eta|<1.5)

| |eta| | N | clean | converted | bad |
|---|---:|---:|---:|---:|
| 0.00-0.10 |   125675 | 84.612% | 15.369% | 0.019% |
| 0.10-0.20 |   125600 | 84.251% | 15.729% | 0.020% |
| 0.20-0.30 |   125254 | 83.741% | 16.237% | 0.022% |
| 0.30-0.40 |   123178 | 82.858% | 17.117% | 0.025% |
| 0.40-0.50 |   121160 | 82.102% | 17.871% | 0.026% |
| 0.50-0.60 |   118921 | 80.865% | 19.107% | 0.028% |
| 0.60-0.70 |   115592 | 79.849% | 20.124% | 0.027% |
| 0.70-0.80 |   110487 | 78.795% | 21.177% | 0.028% |
| 0.80-0.90 |   104868 | 76.980% | 22.989% | 0.031% |
| 0.90-1.00 |    97374 | 75.501% | 24.465% | 0.034% |
| 1.00-1.10 |    89106 | 73.951% | 26.020% | 0.029% |
| 1.10-1.20 |    79427 | 72.713% | 27.256% | 0.030% |
| 1.20-1.30 |    69065 | 71.169% | 28.792% | 0.039% |
| 1.30-1.40 |    57990 | 69.229% | 30.717% | 0.053% |
| 1.40-1.50 |    46986 | 67.482% | 32.482% | 0.036% |

## Physics discussion

Inclusive converted fraction at |eta|<0.7 (14 <= pT < 30 GeV) is **17.317%**, and at |eta|<1.5 it is **21.001%**.

The EMCal front face is preceded by the TPC, TPC outer field cage / support, INTT, and services. In the central barrel the integrated material is usually quoted at the ~10-20% X0 level, rising at larger |eta|. A photon conversion probability of ~(7/9) * X/X0 predicts ~7-15% converted fraction in the central region, increasing toward |eta|~1.5 due to extra path length and edge material.

The observed central-barrel value is in the ballpark of this expectation; if anything it sits a little higher, suggesting the integrated material in front of the EMCal is closer to the upper end of the ~10-20% X0 range. The sharp rise at |eta| > 1 in the per-|eta| table is consistent with endcap-style material growth.

### Cluster-match efficiency per category (14 <= pT < 30 GeV)

- clean: |eta|<0.7 = **67.987%**, |eta|<1.5 = 56.604%
- converted: |eta|<0.7 = **66.290%**, |eta|<1.5 = 46.936%
- bad: |eta|<0.7 = **39.024%**, |eta|<1.5 = 25.761%

Converted photons have a match efficiency that is **-1.697%** (absolute) / **-2.50%** (relative) lower than clean photons in the fiducial region at |eta|<0.7. Most converted photons still produce a single merged cluster (the two e+e- daughters typically deposit within the same tower cluster at these ETs, so the deficit is modest). At |eta|<1.5 the converted match efficiency is substantially lower than clean (~10 percentage points), reflecting that conversions at forward rapidity are more likely to spread over multiple towers or escape the fiducial cluster. The 'bad' category -- where a substantial fraction of the photon momentum is carried by a non-ee secondary -- shows the largest deficit, as expected.

## Concerns / caveats

- The matching uses the **primary truth track ID** of the photon. A conversion can produce e+e- daughters with different track IDs whose cluster may be labeled with the *daughter* trkid instead of the parent photon; such events will be counted as 'not matched' even when a physical cluster is present. The converted-category match efficiency should be read as a **lower bound**.
- The 'bad' fraction is small and dominated by statistics in several bins; bin-to-bin noise is statistical, not physics.
- The cluster Et cut (5 GeV) is below the analysis threshold; tighter Et cuts would widen the clean/converted gap, especially at low pT.
- Single-particle MC: no pileup, no underlying event, no run-range or vertex reweighting.
- The per-pT-bin tables include **all** truth photons that appear in the event record, not just the primary generator photon. Truth bins below ~14 GeV are dominated by low-pT photons produced in the shower / rescattering of the primary (pi0 decays, bremsstrahlung photons, etc.). The inclusive numbers at 14 <= pT < 30 GeV are the ones relevant to isolated-photon physics.
- The absolute match-efficiency values (~65-70% for clean photons at |eta|<0.7) look low at first glance for a truth photon in a single-particle sample. This is because `cluster_truthtrkID` labels the *dominant contributing truth particle* of the cluster, which after showering can be a secondary e+/e- rather than the original photon -- so many physical matches are missed by the strict track-ID equality rule. The *relative* comparison between categories is more meaningful than the absolute number.
