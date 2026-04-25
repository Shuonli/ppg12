# Converted vs Unconverted Photon Isolation Study

**Sample**: `photon10_with_bdt_vtx_noreweight_single.root` (single photon with truth pT in
[14, 30] GeV, one fired photon per event), tree `slimtree`.
**Subsample**: 2,000,000 events out of ~10M total (runtime-based subsampling; the
converted statistics in the top few pT bins are tail-limited but adequate).
**Analysis config cut values** (from `efficiencytool/config_bdt_nom.yaml`):
parametric reco iso cut `iso_max = 0.502095 + 0.0433036 * E_T` [GeV], applied to
`cluster_iso_03_CLUSTERINFO_CEMC`.

## Selection summary

- Truth selection: `particle_pid == 22`, `|particle_Eta| < 0.7`,
  `particle_truth_iso_03 < 4 GeV`.
- Cluster match: highest-`cluster_Et` cluster whose
  `cluster_truthtrkID_CLUSTERINFO_CEMC` matches the truth `particle_trkid`.
- Additional: `cluster_Et > 5 GeV`.

The input sample file does **not** contain `cluster_iso_topo_03/04` branches,
so we use the calo-tower iso (`cluster_iso_03`) for the parametric cut instead
of the production `use_topo_iso: 2` setting. Conclusions on converted vs
unconverted shape differences are expected to carry over to topo iso because
the differences come from the EMCal tower pattern.

## Counts per category (inclusive 8 <= cluster_Et < 36 GeV)

| Category    | Definition                                              | N (2M events) | Fraction |
|-------------|---------------------------------------------------------|--------------:|---------:|
| unconverted | `particle_converted == 0`                               | 517,472 | 85.50% |
| converted   | e+e- pair conversion (`particle_converted == 1`)         |  87,717 | 14.49% |
| bad         | non-ee secondary with >40% momentum (`particle_converted == 2`) |      18 |  0.003% |
| **total**   |                                                         | 605,207 | 100.00% |

The "bad" category is vanishingly rare in photon10, so quantitative
claims below focus on unconverted vs converted only.

## Parametric iso pass fraction per pT bin (cluster_iso_03)

Cut: `iso_03 < 0.502095 + 0.0433036 * cluster_Et`.

| cluster_Et [GeV] | unconv | converted | delta (conv - unconv) |
|-----------------:|-------:|----------:|----------------------:|
| 8 - 10           | 0.562 |   0.334   | **-0.228** |
|10 - 12           | 0.673 |   0.561   | -0.112     |
|12 - 14           | 0.691 |   0.605   | -0.086     |
|14 - 16           | 0.695 |   0.631   | -0.064     |
|16 - 18           | 0.696 |   0.643   | -0.053     |
|18 - 20           | 0.707 |   0.659   | -0.048     |
|20 - 22           | 0.718 |   0.665   | -0.053     |
|22 - 24           | 0.730 |   0.686   | -0.044     |
|24 - 26           | 0.708 |   0.727   | +0.019     |
|26 - 28           | 0.731 |   0.746   | +0.016     |
|28 - 32           | 0.715 |   0.660   | -0.054     |
|32 - 36           | 0.772 |   0.880   | +0.108 (N_conv = 25)  |
| **8 - 36**       | **0.649** | **0.503** | **-0.146** |

The high-pT bins (>= 24 GeV) are limited by converted-photon statistics
(a few hundred events and fewer); their fluctuations around zero are
consistent with no remaining shape difference there.

## Mean isolation per pT bin (cluster_iso_03, GeV)

| cluster_Et [GeV] | unconv | converted | delta |
|-----------------:|-------:|----------:|------:|
| 8 - 10  | 0.881 | 1.847 | **+0.966** |
|10 - 12  | 0.739 | 1.147 | +0.409 |
|12 - 14  | 0.788 | 1.124 | +0.335 |
|14 - 16  | 0.869 | 1.159 | +0.290 |
|16 - 18  | 0.952 | 1.176 | +0.224 |
|18 - 20  | 1.014 | 1.213 | +0.199 |
|20 - 22  | 1.081 | 1.263 | +0.182 |
|22 - 24  | 1.163 | 1.280 | +0.117 |
|24 - 26  | 1.251 | 1.220 | -0.030 |
|26 - 28  | 1.305 | 1.274 | -0.031 |
|28 - 32  | 1.432 | 1.581 | +0.149 |
|32 - 36  | 1.316 | 1.261 | -0.055 |
| **8 - 36** | **0.816** | **1.386** | **+0.570** |

## Mean isolation per pT bin (cluster_iso_04, GeV)

Same shape as R=0.3 but larger absolute values, with similar inclusive
shift (+0.574 GeV).

| cluster_Et [GeV] | unconv | converted | delta |
|-----------------:|-------:|----------:|------:|
| 8 - 10  | 1.029 | 1.998 | **+0.969** |
|10 - 12  | 0.875 | 1.290 | +0.415 |
|12 - 14  | 0.928 | 1.267 | +0.339 |
|14 - 16  | 1.009 | 1.292 | +0.284 |
|16 - 18  | 1.090 | 1.313 | +0.223 |
|18 - 20  | 1.167 | 1.364 | +0.197 |
|20 - 22  | 1.218 | 1.434 | +0.216 |
|22 - 24  | 1.304 | 1.447 | +0.144 |
|24 - 26  | 1.402 | 1.356 | -0.046 |
|26 - 28  | 1.454 | 1.456 | +0.002 |
|28 - 32  | 1.575 | 1.706 | +0.131 |
|32 - 36  | 1.507 | 1.330 | -0.177 |
| **8 - 36** | **0.957** | **1.531** | **+0.574** |

## Which calorimeter drives the shift? (EMCal, inclusive 8-36)

Inclusive means (GeV) of the iso sub-components. Note that the
`_emcal/_hcalin/_hcalout` sub-branches from `CaloAna24.cc` are defined
differently from the total `cluster_iso_03`:

- `cluster_iso_03`: from `RawCluster::get_et_iso(2, false, true)`, i.e.
  the cone integral with the cluster's own towers removed via the
  internal tower threshold scheme (typical mean ~0.8-1.5 GeV).
- `cluster_iso_03_emcal = emcalET(R=0.3) - cluster_ET`: the full EMCal
  cone minus the bare cluster ET. This over-subtracts slightly
  differently, leaving the "halo" of the shower still inside the cone
  (typical mean ~4.7 GeV).
- `cluster_iso_03_hcalin`, `cluster_iso_03_hcalout`: full HCal cones
  with no own-cluster subtraction.

The sub-components should therefore be read as a **diagnostic of where
the converted/unconverted difference sits**, not as a literal
decomposition of the total. The shift pattern is what matters.

### Inclusive shifts (converted - unconverted), GeV

| Iso variable               | unconv | converted | delta    |
|----------------------------|-------:|----------:|---------:|
| cluster_iso_03 (total)     | 0.816 | 1.386     | **+0.570** |
| cluster_iso_03_emcal       | 4.683 | 5.151     | **+0.468** |
| cluster_iso_03_hcalin      | 0.269 | 0.213     | -0.056   |
| cluster_iso_03_hcalout     | 0.287 | 0.253     | -0.033   |
| cluster_iso_04 (total)     | 0.957 | 1.531     | **+0.574** |
| cluster_iso_04_emcal       | 8.166 | 8.484     | **+0.318** |
| cluster_iso_04_hcalin      | 0.297 | 0.241     | -0.056   |
| cluster_iso_04_hcalout     | 0.434 | 0.398     | -0.036   |

### Conclusion — the shift is driven by the EMCal

- Converted photons show ~0.47 GeV (R=0.3) / ~0.32 GeV (R=0.4) more
  EMCal energy in the iso cone than unconverted. This is what drives
  the total iso shift.
- **Both HCal sub-components actually go the other way** (converted
  deposit slightly less than unconverted in the HCal cones,
  delta ~ -0.05 GeV), consistent with the converted e+e- pair showering
  earlier and more compactly in the EMCal, leaving less energy to leak
  to the HCal.
- Physics picture: the e+e- pair from a conversion before the EMCal is
  bent in opposite directions by the solenoid B-field. The cluster
  algorithm reconstructs the leading arm but the other arm's
  electromagnetic shower lands a few towers away — still inside the
  R=0.3 cone for low-pT photons (where the opening after bending is
  largest relative to the cone size). That extra EMCal energy enters
  the iso cone, inflating `iso_03_emcal` and consequently the total
  `iso_03`.
- The effect is largest at low pT (Et_cluster 8-10 GeV: delta_iso03 =
  +0.97 GeV, pass-fraction drop 0.56 -> 0.33) and fades at Et_cluster
  >= 24 GeV, because at high pT the post-B-field opening of the e+e-
  pair is smaller relative to a fixed R=0.3 cone, so more of the
  conversion energy ends up inside the cluster rather than outside.

## Concerns and next steps

- **Ongoing: this file lacks `cluster_iso_topo_03/04`.** The nominal
  analysis uses topo-cluster isolation at R=0.4 (`use_topo_iso: 2`),
  which will give smaller absolute iso values and possibly a different
  scale for the conversion shift (topo subtracts a calibrated event UE
  background). The per-category ordering (converted worse than
  unconverted at low pT) should survive, but the magnitude numbers here
  apply to the calo-tower iso and should be re-measured on the
  topo-iso branches once they are available in photon10_with_bdt
  production. The `apply_BDT.C` or slim-tree builder may need to be
  re-run to populate those branches for this sample.
- The `particle_converted == 2` ("bad") category has only 18 entries in
  2M events, i.e. too few to conclude anything quantitative. Mean iso
  is higher (~2.0 GeV) and pass fraction lower (~0.33), but the
  statistical uncertainty is large. Larger samples (jet MC or the full
  photon5/10/20 set combined) would be needed to quantify this.
- **Implication for the analysis purity / efficiency**: the parametric
  iso cut rejects ~14.6 percentage points more converted than
  unconverted photons inclusively, rising to ~23 pp in the lowest pT
  bin. Since the nominal photon signal shape in the ABCD A region
  already contains conversions at their natural fraction (~14.5% here),
  the reco iso cut is implicitly acting as a partial anti-conversion
  cut at low pT. Two follow-ups are worth considering:
  - **Conversion-aware efficiency systematic**: verify that the signal
    cross-section weights are insensitive to the conversion fraction
    modeling in GEANT by comparing with a data-driven conversion
    tagger (e.g. cluster multiplicity or MVTX hit match, if available).
  - **BDT behavior**: the BDT training likely learns the converted
    shower shape (wider in phi, different wetacogx, etc.). The fact
    that low-pT converted photons also fail the iso cut more means
    low-pT signal efficiency has a double penalty (shower shape + iso).
    Worth checking whether the BDT scale factor varies vs converted
    fraction across pT bins.
- **Cross-check the matching**: I used highest-Et cluster matched on
  `cluster_truthtrkID`. For converted photons, the largest Et arm of
  the e+e- pair is the one that got the largest tower deposits. This
  is the standard PPG12 matching and matches `RecoEffCalculator_TTreeReader.C`
  behavior.

## Output artifacts

- `reports/converted_photon_study/analyze_isolation.py` - analysis code
- `reports/converted_photon_study/rootFiles/isolation_converted.root` -
  12 pT bins + inclusive, 3 categories, 8 iso variables = 312 TH1D
- `reports/converted_photon_study/summary.json` - full numeric summary
  (all means, all pass fractions, all pT bins)
- `reports/converted_photon_study/figures/isolation_cluster_iso_03.pdf` -
  1D normalized overlays per pT bin (13 panels)
- `reports/converted_photon_study/figures/isolation_cluster_iso_03_emcal.pdf`
  (+ `_hcalin`, `_hcalout`) - sub-component overlays per pT bin
- `reports/converted_photon_study/figures/isolation_cluster_iso_04*.pdf` -
  same for R=0.4
- `reports/converted_photon_study/figures/isolation_cluster_iso_03_cumulative.pdf`,
  `..._iso_04_cumulative.pdf` - CDFs with the parametric iso cut drawn
- `reports/converted_photon_study/figures/isolation_mean_vs_pt_R003.pdf`,
  `..._R004.pdf` - mean iso vs cluster_Et, 4 sub-panels (total, EMCal,
  HCalIn, HCalOut), 3 categories overlaid
- `reports/converted_photon_study/figures/isolation_pass_fraction_vs_pt.pdf` -
  parametric iso pass fraction vs cluster_Et, 3 categories overlaid

Regeneration: `python3 analyze_isolation.py --max-events 2000000` (~4 min
wall-time on a cvmfs interactive node).
