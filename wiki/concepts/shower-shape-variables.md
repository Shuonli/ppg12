# Shower Shape Variables

All shower shape variables are computed in `anatreemaker/source/CaloAna24.cc` from a 7x7 tower energy matrix centered on the cluster center-of-gravity (CoG).

## 7x7 Tower Matrix

- Centered on CoG position (not lead tower): `maxieta = floor(avg_eta)`, `maxiphi = floor(avg_phi)`
- Towers below `m_shower_shape_min_tower_E = 0.07 GeV` are zeroed
- Only `isGood()` towers contribute
- Stored in `cluster_e_array_{node}[n*49]` (row-major, eta outer)

## Ring Energy Fractions (et1-et4)

| Variable | Branch | Definition |
|----------|--------|------------|
| et1 | `cluster_et1` | Energy fraction in innermost ring (4 towers adjacent to center) |
| et2 | `cluster_et2` | 2nd ring |
| et3 | `cluster_et3` | 3rd ring |
| et4 | `cluster_et4` | Outermost ring |

From `showershape[0-3]`. Each is normalized by total cluster energy.

## Width Variables

### weta, wphi (integer distance)

Energy-weighted second moment using integer distance from center tower (index 3), considering only owned towers:
```
weta = sum(E[i][j] * (i-3)^2) / sum(E[i][j])
```

### weta_cog, wphi_cog (floating-point CoG distance)

Same but using floating-point distance from the CoG position:
```
cog_eta = 3 + (avg_eta - floor(avg_eta) - 0.5)
weta_cog = sum(E * (i - cog_eta)^2) / sum(E)
```

### weta_cogx, wphi_cogx (CoG, excluding seed)

Same as cog variants but **excluding the seed tower** (i=3, j=3). These are the primary BDT features because they are more sensitive to shower shape differences between photons and hadrons.

Branches: `cluster_weta_cogx_{node}`, `cluster_wphi_cogx_{node}`

### wr_cogx

Radial width: `wr = sqrt(weta_cogx^2 + wphi_cogx^2)`

Not stored as a branch -- computed on the fly where needed.

## Strip Width Variables (w32, w52, w72)

Second moment in eta for a 2-column-wide strip in phi centered on the cluster:

| Variable | Strip | Definition |
|----------|-------|------------|
| w32 | 3x2 | `|deta| <= 1`, 2 adjacent phi columns |
| w52 | 5x2 | `|deta| <= 2`, same phi |
| w72 | 7x2 | `|deta| <= 3`, same phi |

The phi column selection uses `signphi = +1 if fractional phi > 0.5, else -1` to pick the adjacent column closest to the CoG.

## Energy Window Sums (eNM)

Energy sums in NxM tower windows centered on the CoG. N = towers in eta, M = towers in phi.

| Branch | Window | Description |
|--------|--------|-------------|
| `cluster_e11` | 1x1 | Seed tower energy |
| `cluster_e22` | 2x2 | Sum of 4 quadrants (e1+e2+e3+e4) |
| `cluster_e13` | 1x3 | 1 eta x 3 phi |
| `cluster_e15` | 1x5 | |
| `cluster_e17` | 1x7 | |
| `cluster_e31` | 3x1 | 3 eta x 1 phi |
| `cluster_e51` | 5x1 | |
| `cluster_e71` | 7x1 | |
| `cluster_e33` | 3x3 | |
| `cluster_e35` | 3x5 | |
| `cluster_e37` | 3x7 | |
| `cluster_e53` | 5x3 | |
| `cluster_e55` | 5x5 | |
| `cluster_e73` | 7x3 | |
| `cluster_e77` | 7x7 | Full matrix sum |
| `cluster_e32` | 3x2 | Used in strip width |

Computed via: `|di| <= (N-1)/2 && |dj| <= (M-1)/2`

## Energy Ratios (derived, computed in BDTinput.C and apply_BDT.C)

These are NOT raw branches -- they are computed at analysis time:

| Ratio | Formula | Physics |
|-------|---------|---------|
| e11_over_e33 | e11 / e33 | Core concentration (3x3) |
| e32_over_e35 | e32 / e35 | Strip concentration |
| e11_over_e22 | e11 / e22 | 1x1 vs 2x2 |
| e11_over_e13 | e11 / e13 | Eta/phi asymmetry |
| e11_over_e15 | e11 / e15 | |
| e11_over_e17 | e11 / e17 | |
| e11_over_e31 | e11 / e31 | Transposed asymmetry |
| e11_over_e51 | e11 / e51 | |
| e11_over_e71 | e11 / e71 | |
| e22_over_e33 | e22 / e33 | 2x2 vs 3x3 |
| e22_over_e35 | e22 / e35 | |
| e22_over_e37 | e22 / e37 | |
| e22_over_e53 | e22 / e53 | |

Division protection: BDTinput.C uses `safe_div` (denominator != 0 && isfinite). apply_BDT.C has two patterns: the photon-ID block checks the numerator (`cluster_e11 > 0`), while the NPB score block checks the denominator (`cluster_e33 > 0`, `cluster_e35 > 0`, etc.).

## Other Cluster Variables

| Branch | Type | Description |
|--------|------|-------------|
| `cluster_prob` | Float_t | Chi2 probability from clustering |
| `cluster_CNN_prob` | Float_t | CNN classification (always -1, disabled) |
| `cluster_nsaturated` | Int_t | Number of saturated towers |
| `cluster_detamax` | Int_t | Max delta-eta extent |
| `cluster_dphimax` | Int_t | Max delta-phi extent |
| `cluster_ietacent` | Float_t | Continuous eta index of CoG |
| `cluster_iphicent` | Float_t | Continuous phi index of CoG |

## 11 Shower Shape Variables in ShowerShapeCheck.C

ShowerShapeCheck fills distributions for: `weta_cogx`, `wphi_cogx`, `wr_cogx`, `et1`, `et2`, `et3`, `et4`, `e11/e33`, `e32/e35`, `bdt_score`, `npb_score`.

## Usage in BDT

The 25 BDT training features include: `cluster_Et`, `cluster_Eta`, `vertexz`, `weta_cogx`, `wphi_cogx`, `e11_over_e33`, `et1-4`, `e32_over_e35`, `w32`, `w52`, `w72`, and 8 additional energy ratios. Different model variants use subsets -- see [BDT Application](../pipeline/03-bdt-application.md).
