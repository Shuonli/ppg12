# Double-Interaction Efficiency

## Overview

Full GEANT double-interaction MC shows a dramatic drop in overall photon reconstruction efficiency: ~94% in single-interaction MC falls to ~31% in double-interaction MC. This study demonstrates that the drop is **predominantly a truth-matching artifact** caused by vertex displacement, not a genuine physics inefficiency. Signal clusters remain intact in the calorimeter; the matching failure arises because the reconstructed vertex shifts the cluster eta beyond the `deltaR < 0.1` matching threshold.

Key findings:
- BDT shower-shape identification is pileup-robust (identical 81--92% efficiency)
- True isolation degradation is ~4% absolute
- Single-interaction MC is appropriate for efficiency corrections

Cluster-weighted pileup fractions are higher than event-level fractions because double-interaction events produce more clusters on average: 22.4% of clusters at 0 mrad (vs 11.1% event-level) and 7.9% at 1.5 mrad (vs 3.9% event-level). See [Pileup Fractions Reference](../reference/pileup-fractions.md) for the full calculation.

## Available Full-GEANT DI MC Samples

Eight SI/DI sample pairs are available for blending. DI sample cross-sections mirror their SI partners (see `efficiencytool/CrossSectionWeights.h:54-134`). Six of the eight `_double`/`_nom` aliases reuse the SI partner's kinematic window verbatim; `photon10_double` and `jet12_double` carry wider windows (10--100 GeV) as legacy pre-expansion behavior.

| Pair | SI sample | DI sample | SI truth pT window |
|------|-----------|-----------|--------------------|
| photon5 | `photon5` | `photon5_double` | 0--14 GeV |
| photon10 | `photon10` | `photon10_double` | 14--30 GeV |
| photon20 | `photon20` | `photon20_double` | 30+ GeV |
| jet8 | `jet8` | `jet8_double` | 9--14 GeV truth jet pT |
| jet12 | `jet12` | `jet12_double` | 14--21 GeV |
| jet20 | `jet20` | `jet20_double` | 21--32 GeV |
| jet30 | `jet30` | `jet30_double` | 32--42 GeV |
| jet40 | `jet40` | `jet40_double` | 42--100 GeV |

SI filetypes in DI blending contexts (job lists, `ShowerShapeCheck.C`) are passed as `{sample}_nom` (e.g. `photon10_nom`), which `GetSampleConfig` treats identically to `{sample}`.

Paths: `/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/{sample}/condorout/combined.root`.

**Skipped samples:**
- `jet5` — no DI partner produced (`jet5_double` does not exist).
- `jet50_double` — no SI partner in the main pipeline sample set (`jet50` is only used for BDT training, not in `FunWithxgboost/mc_samples.list` or `MergeSim.C`).

## Showershape DI Pipeline (Single-Pass Truth-Vertex Reweight)

The showershape DI pipeline has been migrated to single-pass truth-vertex reweighting and is packaged for condor via `efficiencytool/submit_showershape_di.sub` with job lists `showershape_di_jobs_{0rad,1p5rad}.list` (8 SI/DI pairs + 1 data = 17 jobs per crossing-angle config).

Per-event weight applied in `ShowerShapeCheck.C:1211-1216`:
- Single MC: `w_truth = h_w_iterative.Interpolate(z_truth)`
- Double MC: `w_truth = h_w_iterative.Interpolate(z_hard) * h_w_iterative.Interpolate(z_mb)` (factorized over the two truth vertices)
- A multiplicative `mix_weight` (cluster-weighted DI fraction or its SI complement; 1.0 for data) is also applied.

`h_w_iterative` is derived once per crossing-angle period by the iterative data-driven fit in `efficiencytool/truth_vertex_reweight/fit_truth_vertex_reweight.py`, driven by `truth_vertex_reweight/config.yaml:79-93` (uses SI and DI `photon10`/`jet12` for the fit) and written to `truth_vertex_reweight/output/{0mrad,1p5mrad}/reweight.root`. Runtime loading is handled by `efficiencytool/TruthVertexReweightLoader.h`. The reweight is sample-invariant; the same `h_w_iterative` is applied to every MC sample.

The previous two-pass reco-vertex approach (fill `h_vertexz`, hadd, then rescan with override) is retired. The executable script has been renamed `run_showershape_double.sh` → `run_showershape_double_reco_legacy.sh` and is kept only for reproducing historical results.

The nominal cross-section pipeline (`RecoEffCalculator_TTreeReader.C`) also supports `truth_vertex_reweight_on` / `truth_vertex_reweight_file` and is being extended to the same 8 SI/DI pairs for DI-aware efficiency estimation.

## Truth Photon Content

Before examining reconstruction, the truth-level photon content was verified to be identical between single and double-interaction samples (100k events each from `photon10` / `photon10_double`).

| Quantity | Single | Double |
|----------|--------|--------|
| Truth photons per event | 0.493 | 0.495 |
| Fraction with >= 2 photons | 0.22% | 0.21% |
| Mean truth isolation ET | 0.33 GeV | 0.39 GeV |

The small increase in truth isolation energy (0.33 to 0.39 GeV) is far below the 4 GeV isolation cut and has negligible impact. The photon multiplicity and pT spectra are statistically identical, confirming that the truth-level photon content is not affected by pileup.

## Vertex Shift Mechanism

In double-interaction events, the reconstructed vertex is the average of the two collision vertices:

```
z_vtx_reco = (z_vtx_1 + z_vtx_2) / 2
```

This introduces a shift `dz = z_vtx_reco - z_vtx_true` that propagates to the reconstructed cluster pseudorapidity. For a cluster at the EMCal barrel radius R_EMCal = 93.5 cm:

```
d_eta ~ -dz / (R_EMCal * cosh(eta))
```

The critical vertex shift that causes the `deltaR < 0.1` matching to fail is:

```
|dz|_crit = 0.1 * R_EMCal * cosh(eta) = 9.35 cm   (at eta = 0)
```

When `|dz|` exceeds this threshold, the cluster eta shifts by more than 0.1 and truth-reco matching fails even though the cluster is physically present and well-reconstructed.

| Property | Single | Double |
|----------|--------|--------|
| Vertex z RMS | 7.6 cm | 34.5 cm |
| Fraction with \|dz\| > 9.35 cm | --- | 65% |

The double-interaction vertex RMS of 34.5 cm is far larger than the single-interaction value of 7.6 cm, and 65% of double-interaction events exceed the critical threshold. This directly explains the bulk of the efficiency loss.

## Matching Failure Decomposition

Unmatched truth photons in double-interaction MC are categorized by failure reason:

| Category | Fraction | Description |
|----------|----------|-------------|
| deltaR failure | 92% | Cluster with correct `trkID` exists, but `deltaR > 0.1` due to vertex shift. Signal cluster is present and intact. |
| Missing cluster | 8% | No cluster with matching `trkID` found. Cluster was not reconstructed or was merged. |
| Extra photon | ~0% | Additional truth photons from the second interaction. Negligible. |

The matching failure rate as a function of vertex shift:

| \|dz\| range | Failure rate |
|-------------|-------------|
| < 5 cm | ~10% |
| 5--15 cm | ~41% (spans the sharp transition at 9.35 cm) |
| 15--30 cm | ~100% |
| > 30 cm | 100% |

The sharp transition around the critical threshold and the dominance of deltaR failures confirm that the efficiency drop is overwhelmingly a matching artifact.

## Per-Stage Efficiency Breakdown

Each stage of the photon selection was examined separately to identify which cuts are genuinely affected by pileup. Ranges are quoted for the `photon10` truth pT range (14--30 GeV, bins 3+).

| Efficiency stage | Single | Double | Comment |
|------------------|--------|--------|---------|
| Reconstruction (matching) | 93--94% | 27--32% | Dominated by deltaR matching artifact |
| BDT shower-shape | 81--92% | 81--92% | Identical -- BDT is pileup-robust |
| Isolation | 61--74% | 58--71% | ~4% absolute degradation from extra energy in cone |
| All combined | 52--57% | 15--19% | Drop dominated by reconstruction stage |

Key observations:

- **Reconstruction efficiency** shows the dramatic drop (94% to 31%), driven entirely by the deltaR matching artifact.
- **BDT shower-shape efficiency** is identical between single and double-interaction MC, demonstrating that the BDT is pileup-robust. Clusters that pass the matching cut have shower shapes unaffected by the second interaction.
- **Isolation efficiency** shows a modest ~4% absolute degradation in double-interaction events, from additional energy deposited in the isolation cone by the second collision.
- **Combined efficiency** drops from 52--57% to 15--19%, but this is dominated by the reconstruction matching artifact rather than genuine physics inefficiency.

Note: the first two pT bins fall below the `photon10` truth threshold and have significantly lower efficiencies due to poor statistics and edge effects (e.g., reco ~88%, BDT ~45--68% for single-interaction MC).

## TEfficiency Bug Fix

During this investigation, a bug was identified in ROOT's `TEfficiency` handling of weighted events.

**Broken pattern:**
```cpp
eff->SetWeight(weight);
eff->Fill(passed, x);
```

This stores unweighted event counts internally because `SetWeight` followed by `Fill` does not correctly propagate weights into the numerator and denominator histograms.

**Correct pattern:**
```cpp
eff->SetUseWeightedEvents();
eff->FillWeighted(passed, weight, x);
```

`SetUseWeightedEvents()` configures `TEfficiency` to use `Sumw2()`-enabled internal histograms, and `FillWeighted` correctly accumulates weighted counts. This fix is essential for the blended single+double interaction efficiency study where mixing weights (e.g., 0.224 and 0.776 for the 0 mrad configuration) must be correctly propagated.

## Conclusions

1. **The efficiency drop is a matching artifact, not a physics inefficiency.** Reconstruction efficiency falls from ~94% to ~31% because vertex displacement shifts cluster eta beyond `deltaR < 0.1`. 92% of unmatched photons have their signal cluster present and intact.

2. **Signal clusters are well-reconstructed.** BDT shower-shape efficiency is identical (81--92%) for matched clusters in both samples, demonstrating pileup-robust calorimeter response.

3. **True pileup effects are small.** Isolation degrades by ~4% absolute; missing clusters account for 8% of the unmatched fraction (~5% of all double-interaction truth photons).

4. **Single-interaction MC is appropriate for efficiency corrections.** The efficiency drop is an artifact of truth-matching, not detector response. Real pileup effects are small and treated as systematic uncertainties after weighting by cluster-weighted double-interaction fractions (22.4% at 0 mrad, 7.9% at 1.5 mrad).

## Key Files

| File | Purpose |
|------|---------|
| `efficiencytool/reports/double_interaction_efficiency_report.tex` | Full investigation report with all figures |
| `efficiencytool/DoubleInteractionCheck.C` | Toy vertex-shift simulation (~1446 lines) |
| `efficiencytool/ShowerShapeCheck.C` | Shower shape analysis with `mix_weight` + truth-vertex reweight |
| `efficiencytool/TruthVertexReweightLoader.h` | Loads `h_w_iterative` and applies per-vertex weight (including double-vertex factorization) |
| `efficiencytool/truth_vertex_reweight/` | Iterative data-driven fit producing per-period `reweight.root` |
| `efficiencytool/submit_showershape_di.sub` + `showershape_di_jobs_{0rad,1p5rad}.list` | Condor packaging for the 17-job single-pass DI showershape pipeline |
| `efficiencytool/run_showershape_di_job.sh` | Per-job executable (calls `ShowerShapeCheck.C` with `do_vertex_scan=false`) |
| `efficiencytool/run_showershape_double_reco_legacy.sh` | Retired two-pass reco-vertex blending (kept for history) |
| `efficiencytool/run_double_interaction.sh` | Toy simulation orchestration |
| `efficiencytool/CrossSectionWeights.h` | MC cross-section weights |
| `plotting/plot_double_interaction.C` | Plotting macro for toy sim results |

## See Also

- [ABCD Method](abcd-method.md) -- Background subtraction regions affected by pileup
- [Shower Shape Variables](shower-shape-variables.md) -- BDT input variables shown to be pileup-robust
- [Systematic Variations](systematic-variations.md) -- Pileup treated as a systematic uncertainty
- [Eta Edge Migration](eta-edge-migration.md) -- Fiducial boundary fakes amplified by pileup vertex shift
