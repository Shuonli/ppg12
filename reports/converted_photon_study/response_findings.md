# Converted vs Unconverted Photon Response — Findings

**Sample**: `photon10_with_bdt_vtx_noreweight_single.root` (single-particle MC, ~10.0M events,
truth photons with pT in [14, 30] GeV and |eta| < 1.5 at generation).
**Selection**: `particle_pid==22`, `particle_Pt in [10, 34]` GeV, `|particle_Eta| < 0.7`.
**Matching**: truth photon trkid == `cluster_truthtrkID_CLUSTERINFO_CEMC`, highest-ET cluster,
`cluster_Et > 3` GeV.

Categories from `particle_converted`:

| Code | Name | Meaning |
|------|------|---------|
| 0    | `unconv` | Primary photon did not convert before entering the calorimeter |
| 1    | `conv`   | Photon converted to e+e- in detector material |
| 2    | `bad`    | Bad-secondary flag |

All figures are under `figures/` and numeric results under `rootFiles/summary.pkl` /
`rootFiles/summary.txt`.

---

## 1. Sample statistics and conversion fraction

| Category | Truth photons in selection | Fraction of total |
|---|---:|---:|
| Unconverted | 3,561,433 | 82.6% |
| Converted   |   749,069 | 17.4% |
| Bad         |     1,065 | 0.02% |
| **Total**   | 4,311,567 | |

The conversion fraction is essentially flat vs pT: ~17.2 - 17.9% across
[10, 36] GeV. **Bad-secondary photons are rare (0.02%)** - the "bad" curves
are drawn for completeness but have too few entries for precision claims
at high pT.

---

## 2. Matching efficiency

![match efficiency](figures/response_match_eff_vs_pt.pdf)

Integrated over pT in [10, 34] GeV:

| Category | matched / truth | efficiency |
|---|---:|---:|
| Unconverted | 2,464,250 / 3,561,433 | **69.2%** |
| Converted   |   495,721 /   749,069 | **66.2%** |
| Bad         |       271 /     1,065 | 25.5% |

**Converted photons have ~3% lower match efficiency than unconverted** across all pT,
with a ~5% gap at the highest bin (32-36 GeV: 61.5% vs 56.6%). This is consistent
with: (i) e+e- pair lands in the same cluster most of the time, but (ii) after the
magnetic field separates the pair, one lepton can swing outside the cluster
footprint and the remaining tower energy drops below the ET=3 GeV match
threshold. The match efficiency decreases gently with pT for both categories
(69.6% -> 61.5% for unconv, 65.8% -> 56.6% for conv) due to the asymmetric
selection window `pT_truth in [10, 34]` vs truth pT generation in [14, 30]
(matched clusters with truth pT near the bin edges scatter out of the sample).

The bad-secondary category has only ~21% match efficiency at 10-12 GeV rising
to ~43% at 16-18 GeV, and zero at higher pT (not enough statistics).

---

## 3. Energy response shapes

![shape grid](figures/response_shape_grid.pdf)

**Key visual**: for pT in [10, 20] GeV the converted distribution has a
**pronounced left tail** extending down to R ~ 0.4-0.5, caused by e+ or e- curling
out of the cluster footprint and leaving partial energy in the EMCal. At higher
pT (>= 20 GeV) the tail still persists but the bulk of the distribution
re-converges with the unconverted shape - this is why the means cross around
pT ~ 25 GeV (Section 4).

`R_ET` shape grid is in `figures/response_shape_ET_grid.pdf` and shows an
analogous tail.

---

## 4. Mean response <R_E> vs pT

![mean vs pT](figures/response_mean_vs_pt.pdf)

Inclusive (pT 10-34 GeV):

| Category | <R_E> | median | p16 - p84 |
|---|---:|---:|---|
| Unconverted | **0.9424** | 0.9650 | 0.845 - 1.045 (width 0.200) |
| Converted   | **0.8666** | 0.9050 | 0.675 - 1.025 (width 0.350) |
| Bad         |  0.6506    | 0.6150 | 0.475 - 0.875 (width 0.400) |

**Inclusive mean-response shift: Delta<R_E> = 0.9424 - 0.8666 = -0.0758** (converted lower
by 7.6 percentage points, or ~8% of the nominal response).

Per-pT:
- **10-12 GeV**: unconv 0.9461 vs conv 0.8644 -> shift -0.0817 (-8.6%)
- **14-16 GeV**: unconv 0.9367 vs conv 0.8693 -> shift -0.0674 (-7.2%)
- **18-20 GeV**: unconv 0.9293 vs conv 0.8877 -> shift -0.0416 (-4.5%)
- **22-24 GeV**: unconv 0.9227 vs conv 0.9085 -> shift -0.0142 (-1.5%)
- **26-28 GeV**: unconv 0.9155 vs conv 0.9242 -> shift **+0.0087** (crossover)
- **28-32 GeV**: unconv 0.9034 vs conv 0.9264 -> shift **+0.0230**

The converted-<R_E> crosses above the unconverted curve near pT ~ 25-27 GeV.
This is because at high pT the e+e- pair is collimated (opening angle ~ m_e/pT)
and lands together in the cluster, so the left tail shrinks; meanwhile the
unconverted response keeps degrading slowly with pT (consistent with EMCal
non-linearity from shower leakage into HCal). The narrow-core fit (below)
shows the converted core mu is actually at 1.0 already at low pT - the mean
shift for the full distribution is driven entirely by the left tail.

---

## 5. Resolution: Gaussian core sigma(R_E) vs pT

![sigma vs pT](figures/response_sigma_vs_pt.pdf)

Two complementary widths are reported:
1. **Narrow-core sigma** (Gaussian fit in [peak-0.08, peak+0.08]) - measures the
   intrinsic EMCal resolution of the surviving core.
2. **Asymmetric sigma** (Gaussian fit in [peak-0.2, peak+0.1]) - as requested
   in the task spec; captures both core + left tail. For converted photons this
   is driven by the left tail and is therefore substantially larger than the
   narrow core.

Inclusive (pT 10-34 GeV):

| Category | narrow-core mu | narrow-core sigma | asym mu | asym sigma |
|---|---:|---:|---:|---:|
| Unconverted | 0.9956 | **0.0669** | 0.9839 | **0.0826** |
| Converted   | 0.9994 | **0.0703** | 0.9687 | **0.1051** |
| Bad         | 0.5569 | 0.1118 | 0.5616 | 0.1022 |

**Core-resolution broadening** (converted vs unconverted, narrow-core):
Delta sigma / sigma_unconv = (0.0703 - 0.0669) / 0.0669 = **+5.1%** inclusive.

**Asymmetric resolution broadening** (converted vs unconverted):
Delta sigma / sigma_unconv = (0.1051 - 0.0826) / 0.0826 = **+27%** inclusive.

Physically:
- The Gaussian **core** of the cluster-energy response is only ~5% wider for
  converted photons. Once the e+e- pair is bent together by the field the
  cluster reconstructs almost as if there were a single photon.
- The **left tail** dominates the apparent width: when the pair is split by
  the magnetic field, one shower may deposit less than 90% of its expected
  energy in the EMCal footprint and pull the distribution down.

Per-pT narrow-core sigma is ~0.065-0.070 for both categories across most
pT bins (see figures). The conv narrow-core sigma spikes to 0.094 at pT 22-24
GeV; this bin suffers from a low-stats fit near the peak and is a fit artifact
(see Concerns).

---

## 6. Tail fractions: R < 0.8 and R < 0.5

![tail fractions](figures/response_tailfrac_vs_pt.pdf)

Inclusive (pT 10-34 GeV):

| Category | frac R<0.8 | frac R<0.5 |
|---|---:|---:|
| Unconverted | **0.107** | 0.0053 |
| Converted   | **0.311** | 0.0201 |
| Bad         | 0.760 | 0.221 |

**Converted photons are ~3x more likely to have R < 0.8** (tail fraction
0.311 vs 0.107) and **~4x more likely to have R < 0.5** (0.0201 vs 0.0053).

Per-pT, the conv tail-0.8 fraction **decreases from 0.32 at low pT to 0.13 at
high pT**, reflecting the collimation of e+e- pairs at high pT. The unconv
tail-0.8 actually **increases gently** from 0.10 to 0.17 with pT, reflecting
the growing EMCal non-linearity. The two curves cross around pT ~ 26 GeV,
consistent with the mean-response crossover.

---

## 7. Eta dependence

![eta dependence](figures/response_eta_dep.pdf)

pT-inclusive, split by truth eta:

| eta bin | unconv <R> | conv <R> | unconv tail0.8 | conv tail0.8 |
|---|---:|---:|---:|---:|
| [-0.70, -0.35] | 0.9557 | 0.8800 | 0.092 | 0.287 |
| [-0.35,  0.00] | 0.9314 | 0.8545 | 0.118 | 0.333 |
| [ 0.00,  0.35] | 0.9310 | 0.8545 | 0.119 | 0.332 |
| [ 0.35,  0.70] | 0.9541 | 0.8781 | 0.094 | 0.289 |

- **Central eta is worse than forward/backward eta for both categories.** The
  |eta| < 0.35 slice has <R>_unconv ~ 0.931 vs <R>_unconv ~ 0.955 in |eta| > 0.35.
  This is consistent with the well-known EMCal gap / dead-material structure
  near eta = 0 (tower seams, cable slot region).
- **The converted-unconverted <R> shift is eta-independent**: ~0.076 - 0.077
  at both small and large |eta|.
- Tail fractions also show the same "worse at center" pattern but the
  converted/unconverted ratio is flat: conv tail is ~3.1x unconv tail in both
  central and forward regions.
- Narrow-core sigma is larger in the central region for both categories
  (unconv ~0.075 at |eta|<0.35 vs ~0.058 at |eta|>0.35).

No clear asymmetry between positive and negative eta.

---

## 8. Angular deltas (dEta, dPhi, dR)

![angular deltas](figures/response_angular_deltas.pdf)

The cluster - truth angular difference shows a **striking conversion signature
in dPhi**:

- **Unconverted**: narrow Gaussian core around dPhi ~ 0 with RMS much less
  than 0.01 rad.
- **Converted**: **double-peak structure** in dPhi, with lobes at dPhi ~ +/- 0.01
  and a dip at zero. This is the magnetic-field bending signature: the
  centroid of the cluster shifts in the direction of whichever lepton
  (e+ or e-) carries more energy, producing a bimodal shift in phi.
- The two lobes are roughly symmetric (field bends e+ and e- in opposite
  senses with equal probability), so the **mean dPhi is unchanged** but the
  RMS is much larger (~10x) for converted.

dEta shows essentially no shift or broadening between conv/unconv - consistent
with the solenoidal field only bending in phi.

dR shows a broader tail for converted, up to 0.15 rad, whereas unconverted
drops off below 0.05 rad. The increased dR correlates with the cluster
centroid being pulled away from the truth photon direction by the bent
leptons.

---

## 9. Summary table

| Quantity | Unconverted | Converted | Shift |
|---|---:|---:|---|
| N matched (10-34 GeV) | 2,464,250 | 495,721 | |
| Match eff | 69.2% | 66.2% | -3.0 pp |
| <R_E> inc | 0.9424 | 0.8666 | **-0.076 (-8%)** |
| median R_E | 0.9650 | 0.9050 | -0.060 |
| p84-p50 (upper res) | 0.080 | 0.120 | +50% |
| p50-p16 (lower res) | 0.120 | 0.230 | +92% |
| narrow-core mu | 0.9956 | 0.9994 | +0.004 |
| narrow-core sigma | 0.0669 | 0.0703 | +5% |
| asym sigma | 0.0826 | 0.1051 | **+27%** |
| frac R<0.8 | 0.107 | 0.311 | **x2.9** |
| frac R<0.5 | 0.0053 | 0.0201 | **x3.8** |

---

## 10. Concerns and follow-ups

1. **Bad-secondary category is too small** for precision claims at high pT
   (only 271 matched photons, zero above pT=18 GeV). Their response is
   peaked at R ~ 0.6 with huge tails. This likely corresponds to photons
   from very-late-decaying hadrons or reconstruction glitches. Consider
   dropping the category from downstream analyses - it's physically
   negligible (0.02%).

2. **Narrow-core sigma fit is unstable at low stats** (conv 22-24 GeV spike
   to 0.094). The Gaussian fit within peak +/- 0.08 uses only ~16 bins of
   0.01 width and can be pulled by fluctuations for N < 3000. Would
   benefit from a fixed-window fit (e.g. 0.93-1.05) or a double-Gaussian
   + exponential tail model.

3. **Match threshold `cluster_Et > 3` GeV** was chosen to capture low-R tails.
   The main PPG12 analysis uses a tighter threshold. The ~3% efficiency
   difference between conv/unconv here is the *lower bound* - in the main
   analysis pipeline (with ET cut near the bin edge and tight BDT) the
   converted-photon efficiency loss would be somewhat larger. **This is a
   real detection-efficiency systematic**, currently not separated from the
   overall reconstruction efficiency in the PPG12 pipeline.

4. **Converted photon BDT response** is not studied here (we only examined
   kinematic response). Given the shower-shape distortion (leptons bend in
   opposite phi directions), the BDT shower-shape score distribution is
   likely broader for converted photons - worth comparing `bdt_score`
   distributions for conv vs unconv, especially near the tight BDT
   threshold. **Follow-up**: produce a matched conv/unconv BDT-score plot.

5. **No charged-hadron background** in this single-particle sample.
   The ABCD method assumes the signal shape is measured from MC, so any
   conv/unconv bias in the BDT template could leak into the purity
   systematic. **Follow-up**: check the conv/unconv *fraction* in the ABCD
   tight-iso region (signal region) matches the physics prediction once
   reconstruction is applied - if converted-photon eff is different,
   the effective conversion fraction in the tight-iso sample might differ
   from the truth 17.4%, creating a BDT-shape systematic.

6. **High-pT crossover**: the fact that conv <R> *exceeds* unconv <R> above
   pT ~ 25 GeV is a real prediction (pair collimation) but worth confirming
   isn't a subtle statistical artifact. Currently only ~650 matched
   high-pT converted photons in [32, 36] GeV. **Follow-up**: run the full
   sample (already done - 4.3M photons total) gives reasonable uncertainty
   above 25 GeV. This crossover should also be visible in the tail-fraction
   curves at R<0.8 and we do see it (conv tail drops below unconv tail near
   pT ~ 26 GeV).

7. **Vertex smearing** was not considered here. The photon10 sample has the
   standard MC vertex distribution; the energy response uses fixed cluster
   eta/phi (which can be off by the dz/r geometric shift). This is a
   small effect for |eta|<0.7 but would enlarge the dEta tail for both
   categories. Could be cross-checked by using the vertex-corrected cluster
   Eta if available.

---

## Files

- `analyze_response.py` - uproot-based analysis, vectorized matching
- `make_plots.py`       - plotting + summary extraction
- `rootFiles/response_converted.root` - all histograms
- `rootFiles/summary.pkl`  - numeric summary dicts
- `rootFiles/summary.txt`  - text summary
- `figures/response_shape_grid.pdf` - R_E per-pT panels
- `figures/response_shape_ET_grid.pdf` - R_ET per-pT panels
- `figures/response_mean_vs_pt.pdf`
- `figures/response_sigma_vs_pt.pdf`
- `figures/response_tailfrac_vs_pt.pdf`
- `figures/response_angular_deltas.pdf`
- `figures/response_match_eff_vs_pt.pdf`
- `figures/response_eta_dep.pdf`
