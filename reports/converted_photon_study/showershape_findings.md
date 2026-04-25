# Converted vs Unconverted Photon Shower-Shape Study

**Sample**: `photon10_with_bdt_vtx_noreweight_single.root` (9,998,541 events,
full file processed).

**Input file:** `/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_with_bdt_vtx_noreweight_single.root`

**Script:** `reports/converted_photon_study/analyze_showershape.py`
(full 10M-event run, ~2m20s wall on a single core; **no subsampling applied**).

**Tree:** `slimtree`.
**Cluster nodes:**
- shower-shape + CNN prob -> `CLUSTERINFO_CEMC`
- BDT score -> `CLUSTERINFO_CEMC_NO_SPLIT` (only node where BDT was written).

---

## 1. Selection and category definitions

Truth-photon selection (`reports/converted_photon_study/analyze_showershape.py:59-63,390-398`):

```
particle_pid == 22
|particle_Eta| < 0.7
particle_truth_iso_03 < 4 GeV        (fiducial truth isolation)
8 <= particle_Pt <= 36 GeV
```

Reco match (`reports/converted_photon_study/analyze_showershape.py:248-333`, `match_clusters_vec`):

```
cluster_truthtrkID == particle_trkid
Pick highest-ET cluster match
Require cluster_Et > 5 GeV
```

pT binning follows the PPG12 analysis pT bins
(`ptRanges = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]`, 12 bins)
using the **matched cluster ET** as the x-axis
(`reports/converted_photon_study/analyze_showershape.py:42-43`).

Conversion categories follow `CaloAna24.cc:953-959`:

```
particle_converted == 0 -> unconverted   (no "bad" or converted flag)
                   == 1 -> converted     (truth e+e- from material)
                   == 2 -> "bad" photon  (secondary non-e+e- from material
                                          carrying >40% of the photon momentum)
```

The sample is not a single-particle gun — it is a Pythia photon-jet MC with
`photonclass == 1` (direct), `2` (fragmentation), or `3` (decay) making up the
photons. The conversion/bad flag is independent of photonclass. We do *not*
filter on `particle_photonclass` because the converted/not-converted split is
a pre-reco detector-material effect applicable to any truth photon.

## 2. Category counts

After the selection and matching described above, the number of matched
truth photons in each category, broken down by matched cluster ET:

| pT bin [GeV] | Unconverted | Converted | Bad photon |
|--------------|------------:|----------:|-----------:|
|  8  -  10    |     723,749 |   148,305 |         68 |
| 10  -  12    |   1,027,045 |   157,409 |         26 |
| 12  -  14    |     464,818 |    71,947 |         19 |
| 14  -  16    |     197,106 |    31,035 |          1 |
| 16  -  18    |      87,700 |    14,872 |          4 |
| 18  -  20    |      41,277 |     7,116 |          2 |
| 20  -  22    |      19,913 |     3,715 |          0 |
| 22  -  24    |      10,072 |     1,913 |          0 |
| 24  -  26    |       5,179 |     1,015 |          0 |
| 26  -  28    |       2,643 |       572 |          0 |
| 28  -  32    |       2,018 |       454 |          0 |
| 32  -  36    |         373 |       106 |          0 |
| **sum 8-36** | **2,581,893** | **438,459** | **120** |

Matched photons whose matched cluster ET falls outside [8, 36] GeV are excluded
from this breakdown (cluster ET used for binning, and some truth photons with
`particle_Pt >= 8 GeV` get reconstructed with cluster ET < 8 GeV due to energy
losses and resolution). Adding those, **the total matched in the selection is:**
- unconverted: 2,783,265
- converted:   555,613
- bad:         281

**Fractions** (pass-selection + matched, before pT-binning cut): unconv=83.34%,
conv=16.65%, bad=0.01%. The 16.6% conversion fraction is higher than the naive
`1-exp(-7/9 * 0.15) ~ 11%` estimate for 0.15 X0 of tracker material; the
extra ~5% likely comes from additional passive material (TPC gas, MVTX/INTT
cables, support structure) and from photons that convert late within the
EMCal tower support (still flagged as conversion upstream of the active
volume). Worth cross-checking against the simulation X0 map.

**Bad photons are statistically negligible (281 total).** I keep them in the
ROOT output but suppress them from the 1D shower-shape panels and
inclusive plots (threshold of 200 entries per panel; see
`MIN_CAT_ENTRIES`, `analyze_showershape.py:543`). Mean-vs-pT plots still show
bad at low pT to illustrate the separation magnitude.

## 3. Inclusive (8 <= pT < 36 GeV) shower-shape shifts

Means from the 8-36 GeV inclusive histograms (category weighted by cluster ET):

| Variable             | unconv mean | conv mean | Delta (conv - unconv) | Delta/unconv |
|----------------------|------------:|----------:|----------------------:|-------------:|
| weta_cogx            | 0.191       | 0.225     | +0.034                | +17.8%       |
| wphi_cogx            | 0.128       | 0.244     | +0.116                | **+91.2%**   |
| ET1 (leading tower)  | 0.939       | 0.911     | -0.028                | -3.0%        |
| ET2                  | 0.542       | 0.498     | -0.043                | -8.0%        |
| ET3                  | 0.661       | 0.586     | -0.075                | -11.3%       |
| ET4                  | 0.042       | 0.050     | +0.008                | +17.8%       |
| e11/e33              | 0.649       | 0.597     | -0.052                | -8.0%        |
| e32/e35              | 0.980       | 0.946     | -0.034                | -3.4%        |
| et2/et1              | 0.572       | 0.539     | -0.034                | -5.9%        |
| et3/et1              | 0.705       | 0.639     | -0.067                | -9.5%        |
| et4/et1              | 0.045       | 0.056     | +0.011                | +24.1%       |
| CNN prob             | 0.638       | 0.464     | **-0.174**            | **-27.3%**   |
| cluster_prob (chi2)  | 0.565       | 0.415     | -0.149                | -26.5%       |
| merged_prob          | 0.804       | 0.673     | -0.131                | -16.3%       |
| **BDT score**        | **0.736**   | **0.616** | **-0.121**            | **-16.4%**   |

### Key observations

- **wphi_cogx is the single biggest discriminator**: converted photons widen
  in phi by ~91% relative, because the e+e- pair opens in the bending plane.
  This also explains why the PPG12 BDT's parametric tight cut
  (`tight_wphi_cogx_max = tight_wphi_cogx_max_b + tight_wphi_cogx_max_s*ET`,
  see `FunWithxgboost/BDTinput.C:754-755`) rejects converted photons
  preferentially.

- **weta_cogx widens by ~18%** inclusively (approx +0.03 absolute). The effect
  is largest in the lowest pT bin (20% at 9 GeV) and shrinks at high pT where
  the e+e- opening angle is small (table below).

- **ET1 / ET1-ring fractions shift downward** for converted photons
  (e.g. et3 drops 11% and et3/et1 by 9.5%), consistent with energy leaking
  out of the central 3x1 ring into wider rings when the shower is composed of
  two close photons.

- **CNN photon probability** is pushed strongly toward 0: mean drops from 0.638
  to 0.464 (-27.3% rel.). Looking at the inclusive 1D (fig
  `showershape_CNN_prob_inclusive.pdf`), converted photons have a ~1.6x
  higher peak at prob~0 and a ~1.7x lower peak at prob~1. This is by far the
  most discriminating single variable for converted/unconverted.

- **BDT score (NO_SPLIT node)** drops by ~16% in mean (0.736 -> 0.616). The
  inclusive 1D (fig `showershape_bdt_score_nosplit_inclusive.pdf`) shows the
  converted-photon distribution has a long low-score tail (density ~0.6 at
  BDT=0.3) that is essentially absent in unconverted (density ~0.1 at
  BDT=0.3). **Converted photons fail the tight BDT selection at a much higher
  rate than unconverted photons**, which is the efficiency driver to quantify
  in the main analysis.

## 4. pT dependence of key shifts

Selected pT-bin centers, unconv vs conv mean (signed difference, relative %):

| pT centre [GeV] | weta_cogx conv-unconv (rel%) | wphi_cogx (rel%) | e11/e33 (rel%) | BDT (rel%)   | CNN (rel%)    |
|-----------------|-----------------------------:|-----------------:|----------------:|--------------:|--------------:|
|  9              | +0.041 (+20.9%)              | +0.173 (+134%)   | -0.050 (-8.1%) | -0.135 (-18.4%) | -0.221 (-33.5%) |
| 11              | +0.031 (+16.2%)              | +0.087  (+69%)   | -0.047 (-7.2%) | -0.109 (-14.7%) | -0.150 (-22.9%) |
| 13              | +0.029 (+15.2%)              | +0.083  (+65%)   | -0.049 (-7.3%) | -0.106 (-14.4%) | -0.146 (-23.2%) |
| 15              | +0.028 (+15.0%)              | +0.088  (+68%)   | -0.052 (-7.9%) | -0.116 (-15.9%) | -0.150 (-24.8%) |
| 17              | +0.029 (+14.9%)              | +0.095  (+72%)   | -0.063 (-9.4%) | -0.131 (-18.1%) | -0.162 (-28.4%) |
| 19              | +0.029 (+15.3%)              | +0.103  (+77%)   | -0.069 (-10.5%)| -0.151 (-21.1%) | -0.178 (-34.2%) |
| 21              | +0.030 (+15.8%)              | +0.105  (+78%)   | -0.076 (-11.5%)| -0.162 (-22.7%) | -0.167 (-38.1%) |
| 23              | +0.030 (+15.2%)              | +0.111  (+81%)   | -0.081 (-12.3%)| -0.178 (-25.3%) | -0.140 (-43.7%) |
| 25              | +0.026 (+13.2%)              | +0.103  (+76%)   | -0.088 (-13.2%)| -0.176 (-25.0%) | -0.108 (-52.3%) |
| 27              | +0.031 (+15.9%)              | +0.100  (+72%)   | -0.084 (-12.7%)| -0.171 (-24.7%) | -0.047 (-40.9%) |
| 30              | +0.021 (+10.5%)              | +0.097  (+68%)   | -0.069 (-10.4%)| -0.166 (-24.2%) | -0.010 (-23.6%) |
| 34              | +0.005 (+2.1%)               | +0.071  (+48%)   | -0.058 (-9.2%) | -0.168 (-25.8%) | +0.007 (+63%) |

- **weta_cogx**: steady ~15% shift from 11-30 GeV, tails off at 34 GeV
  (poor stats there, 106 conv / 373 unconv).
- **wphi_cogx**: biggest at low pT (134% at 9 GeV) due to larger e+e-
  opening angle, decreases but remains ~70% over 11-30 GeV.
- **BDT shift grows** with pT (from -14% at 11 GeV to -26% at 23-34 GeV).
  Converted-photon inefficiency worsens at high pT.
- **CNN shift changes sign** at 34 GeV but with only 106 converted
  entries there the statistic is noisy; above 28 GeV the CNN is
  close-to-zero for everyone so relative shifts are unstable.

## 5. Discriminating variables ranking (inclusive)

Ranked by `|Delta mean| / sigma_unconv` (separation power in "sigmas",
bigger = more useful for rejecting converted photons), at inclusive 8-36 GeV:

| Variable           | |shift| | sigma_unconv | `|shift|/sigma_u` |
|--------------------|--------:|-------------:|------------------:|
| e32/e35            | 0.0337  | 0.0145       | **2.33**          |
| wphi_cogx          | 0.1163  | 0.0540       | **2.15**          |
| BDT score          | 0.1205  | 0.1621       | 0.74              |
| ET1                | 0.0283  | 0.0512       | 0.55              |
| merged_prob        | 0.1307  | 0.2385       | 0.55              |
| cluster_prob       | 0.1493  | 0.3134       | 0.48              |
| CNN prob           | 0.1738  | 0.4302       | 0.40              |
| e11/e33            | 0.0522  | 0.1688       | 0.31              |
| ET3                | 0.0750  | 0.2442       | 0.31              |
| weta_cogx          | 0.0341  | 0.1277       | 0.27              |
| et3/et1            | 0.0667  | 0.2592       | 0.26              |
| et4/et1            | 0.0109  | 0.0597       | 0.18              |
| ET2                | 0.0434  | 0.2836       | 0.15              |
| ET4                | 0.0075  | 0.0554       | 0.14              |
| et2/et1            | 0.0337  | 0.2925       | 0.12              |

The most-used BDT features (weta_cogx, wphi_cogx, et1, e11/e33, et4/et1,
see `FunWithxgboost/config.yaml:63-87`) all shift in the direction that
makes converted photons look more background-like, i.e. the BDT training is
implicitly rejecting converted photons even though no truth-level conversion
information is passed to the BDT.

**The top two discriminators (e32/e35 and wphi_cogx) both test for shower
opening in the bending plane** — the e+e- pair produced by conversions is
bent by the ~1.4 T sPHENIX magnetic field in the phi direction, splitting
the shower across 3x5 vs 3x2 towers. e32/e35 tightens that into a very
clean ~2.3-sigma separator.

## 6. Surprises and concerns

1. **16.6% conversion fraction is high.** This is the fraction of truth
   photons flagged as e+e- converted (`particle_converted==1`). Given the
   sPHENIX material budget quoted at ~0.15 X0 in front of the EMCal, the
   naive conversion probability is `1 - exp(-(7/9) * 0.15) ~ 12%`. The
   measured 16.6% is a bit higher and may include secondary productions.
   Worth cross-checking against detector simulation expectations
   documented elsewhere.

2. **BDT score does not visibly separate converted/unconverted at high score**
   (both peak near 0.85) but has a large low-score tail for converted. This
   means the BDT efficiency is NOT a simple convolution — the tight BDT cut
   at `tight_bdt_min = intercept + slope * ET` (see
   `efficiencytool/config_bdt_nom.yaml`) will preferentially reject converted
   photons, and the efficiency loss increases with pT (mean BDT shift grows
   from -15% at 11 GeV to -26% at 34 GeV).

3. **"Bad" photons are essentially absent** (281 total in 10M events, all
   below 20 GeV). Their shower-shape distributions are consistent with being
   mislabeled secondaries: very low BDT scores (~0.2-0.4), broad weta_cogx.
   The efficiency impact on the cross-section measurement is negligible.

4. **CNN prob saturates to 0 above ~25 GeV for both categories.** Mean drops
   from 0.64 (9 GeV) to 0.04 (30 GeV) for unconverted. The CNN was trained
   on a different pT range and its output saturates here — interesting
   sanity check, but it is not used in the analysis selection so it's
   informational only. See `showershape_CNN_prob_mean.pdf`.

5. **wphi_cogx is strongly bimodal for converted photons** in each pT bin
   (see `showershape_wphi_cogx.pdf`) — there are two populations, one
   matching unconverted (unopened pairs) and one at high wphi
   (wide-opening pairs). This bimodality is visible also in `wphi_cogx_mean`
   plot where the mean is high but the RMS (std) is 5x larger for conv
   (0.263) than unconv (0.054).

## 7. Files produced

```
reports/converted_photon_study/
  analyze_showershape.py            <-- main script
  showershape_findings.md           <-- this file
  rootFiles/
    showershape_converted.root      <-- all histograms + mean/rms/nent per (var, cat)
  figures/
    showershape_<VAR>.pdf           <-- 4x3 multi-panel (one per pT bin)
    showershape_<VAR>_inclusive.pdf <-- single canvas 8-36 GeV
    showershape_<VAR>_mean.pdf      <-- mean vs pT, one line per category
```

`VAR` ranges over: `weta_cogx`, `wphi_cogx`, `et1`-`et4`, `CNN_prob`,
`cluster_prob`, `merged_prob`, `e11_over_e33`, `e32_over_e35`,
`et2_over_et1`, `et3_over_et1`, `et4_over_et1`, `bdt_score_nosplit`.

## 8. Next steps

- Compute the **BDT tight-selection efficiency** separately for converted /
  unconverted using the parametric tight cut
  (`tight_bdt_min_intercept + tight_bdt_min_slope * ET` from
  `efficiencytool/config_bdt_nom.yaml`); this is the actual efficiency
  systematic that enters the photon cross-section measurement.
- Investigate whether re-training the BDT with `converted` as a target-flag
  feature could equalize efficiency (at the cost of introducing truth-level
  information into the classifier).
- Repeat this study on the 1-5 GeV (`photon5`) and 20-40 GeV (`photon20`)
  samples to confirm pT dependence of the conversion efficiency loss.
- Check whether the data/MC show consistent e11/e33 and wphi_cogx
  bimodality shape vs MC prediction (data does not carry truth conv flag,
  so shape integrals vs the unconv-only MC expectation).
