# Unfolding

## Purpose

Bayesian unfolding corrects the measured (reco-level) photon ET spectrum for detector smearing effects, mapping it to the truth-level pT spectrum. This accounts for finite energy resolution, bin migration, and acceptance effects.

## Method: RooUnfoldBayes

The analysis uses iterative Bayesian unfolding via `RooUnfoldBayes` from the RooUnfold package.

### Response Matrix

Built in `RecoEffCalculator_TTreeReader.C` for tight + isolated + signal-matched clusters only:

```cpp
responses_full[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
```

Stored as `RooUnfoldResponse` objects:
- `response_matrix_full_0` -- full-sample response (used for data unfolding)
- `response_matrix_half_0` -- first-half-of-events response (used for closure tests)

Also stored as TH2D: `h_response_full_0`, `h_response_half_0`

Truth and reco projections: `h_pT_truth_response_0`, `h_pT_reco_response_0`

### Binning

- **Reco pT bins:** from config `pT_bins` (nominal: `[8,10,12,14,16,18,20,22,24,26,28,32,36]`)
- **Truth pT bins:** from config `pT_bins_truth` (nominal: `[7,8,10,12,14,16,18,20,22,24,26,28,32,36,45]`)

Truth bins extend beyond reco bins on both sides for unfolding overflow handling.

### Iterations

RooUnfoldBayes is run for iterations 1-10, each stored separately:
- `h_unfold_sub_{1..10}` -- without leakage correction
- `h_unfold_sub_leak_{1..10}` -- with leakage correction

The selected iteration is `analysis.unfold.resultit` (default: 2).

### Response Reweighting

Optional reweighting of the response matrix to match data prior vs MC prior. Controlled by `analysis.unfold.reweight` (currently forced to 0 in code).

## Closure Tests

### Full Closure

Run `CalculatePhotonYield.C` with `isMC=true`: uses merged jet MC as "data" and compares unfolded result to truth. Output: `Photon_final_{var_type}_mc.root`.

### Half Closure

Response built from first half of events (`response_matrix_half_0`), unfolded spectrum from second half. Tests whether the response matrix correctly describes the data it wasn't trained on.

Implemented in `plotting/Closure.C` and `plot_unfold_iter.C`.

### Convergence

`plot_unfold_iter.C` shows the relative change between consecutive iterations. Typically converges by iteration 2-3.

## Leakage Correction in Unfolding

When `resultleak = 1` (default), the purity-corrected yield used for unfolding already accounts for signal leakage via the self-consistent ABCD equation. When `resultleak = 0`, the simpler (no-leakage) ABCD solution is used.

## Efficiency Correction

After unfolding, the spectrum is divided by total efficiency per truth bin:

```cpp
float total_eff = eff_reco * eff_iso * eff_id * vertex_eff;
```

The vertex efficiency includes MBD efficiency: `h_truth_pT_vertexcut_mbd_cut / h_truth_pT_vertexcut + mbd_eff_scale`.

## Config Fields

| Field | Default | Purpose |
|-------|---------|---------|
| `unfold.resultit` | 2 | Number of Bayesian iterations to use |
| `unfold.resultleak` | 1 | 0=no leakage, 1=with leakage correction |
| `unfold.reweight` | 1 (config) / 0 (code) | Response reweighting (currently disabled) |
| `pT_bins` | [8..36] | Reco bin edges |
| `pT_bins_truth` | [7..45] | Truth bin edges (wider) |

## Key Output Histograms

| Histogram | File | Description |
|-----------|------|-------------|
| `h_unfold_sub_result` | Photon_final_*.root | Final cross-section (THE primary result) |
| `h_unfold_sub_result_woeff` | Photon_final_*.root | Without efficiency correction |
| `h_unfold_sub_{1..10}` | Photon_final_*.root | Each iteration |
| `h_response_full_0` | MC_response_*.root | Response matrix (TH2D) |
| `response_matrix_full_0` | MC_response_*.root | Response matrix (RooUnfoldResponse) |
| `h_pT_truth_response_0` | MC_response_*.root | Truth projection |
| `h_pT_reco_response_0` | MC_response_*.root | Reco projection |

## Plotting

- `plot_unfold_iter.C` -- iteration convergence, closure tests
- `plot_response.C` -- response matrix visualization
- `plot_reweight.C` -- unfolding prior comparison
- `Closure.C` -- full/half closure test driver

## Systematic Uncertainty

Unfolding systematics come from:
- Number of iterations (varied from nominal)
- Response matrix reweighting (prior dependence)
- Grouped as `unfolding` in the systematic aggregation

See [Systematic Variations](systematic-variations.md).

---

## Closure Test Details

### What Closure Tests Validate

Closure tests confirm that the unfolding procedure is unbiased and correctly recovers the truth spectrum. Two complementary tests are implemented in `Closure.C`:

- **Full closure:** Use the entire MC reconstructed distribution as pseudo-data, unfold with the response matrix built from the same MC, and compare the unfolded result to the MC truth. This tests whether the machinery itself introduces bias. A perfect result is unfolded/truth = 1.0 in every bin.
- **Half closure:** Split the MC into two independent halves by event index (first 50% and second 50%). Build the response matrix from the first half. Use the second half reco distribution as pseudo-data. Unfold and compare to the second half truth. This is a more realistic test because the response matrix and the pseudo-data are statistically independent.

### Running Closure Tests

The closure test macro lives in `plotting/Closure.C`. It requires a merged response matrix file produced by the efficiency pipeline (Stage 4).

Prerequisite -- generate and merge response matrices:

```bash
cd efficiencytool
root -l -b -q 'RecoEffCalculator_TTreeReader.C("config_bdt_none.yaml", "photon5")'
root -l -b -q 'RecoEffCalculator_TTreeReader.C("config_bdt_none.yaml", "photon10")'
root -l -b -q 'RecoEffCalculator_TTreeReader.C("config_bdt_none.yaml", "photon20")'
root -l -b -q 'MergeSim.C("config_bdt_none.yaml")'
```

Run the closure test:

```bash
cd plotting
root -l -b -q 'Closure.C("config_bdt_none.yaml")'
```

Note: the current `Closure.C` only accepts a single `configname` parameter and has a hardcoded response file path that may need manual editing. See "Closure.C Current State and Planned Extensions" below for details.

### Closure Test Output

The current `Closure.C` produces a single PNG plot per eta bin:

- `figure/h_unfold_reco_full_1.png` -- full closure with ratio panel (1 iteration)

The plot has two panels: the top panel shows unfolded vs truth distributions (with reco for reference), the bottom panel shows the ratio unfolded/truth (target: 1.0).

The planned extensions (see below) would add chi2-vs-iteration graphs, half closure plots, and a comprehensive ROOT output file.

### Interpreting Closure Results

Good closure:
- Ratio plots centered at 1.0 across all pT bins
- RMS of ratio < 5% for most bins
- Mean ratio stable after 2-4 iterations
- Chi2/ndf close to 1.0

Poor closure indicators:
- Systematic deviation from 1.0 (bias in the unfolding)
- Large bin-to-bin fluctuations (insufficient MC statistics)
- Ratio trending with pT (bin migration not well modeled)
- Strong iteration dependence (over-unfolding or under-unfolding at different iteration counts)

Optimal iteration selection: choose the iteration where bias is minimized and chi2/ndf is closest to 1.0. Typically 2-4 iterations for Bayesian unfolding. Too few leaves regularization bias; too many amplifies statistical noise.

### Systematic Closure Studies

Once `Closure.C` is extended with multi-iteration and config-driven support, the closure test can be repeated for each systematic variation config to check robustness.

Recommended checks when running closure studies:
- Iteration stability: compare iterations 1-10, find the plateau
- Bin migration: inspect the 2D response matrix for off-diagonal strength
- Statistical tests: chi2/ndf should be near 1.0 for every config
- Reweighting impact: compare `reweight: 0` vs `reweight: 1` in the config

## Code Review Findings

A review of the unfolding code (`RecoEffCalculator_TTreeReader.C` and `CalculatePhotonYield.C`) confirmed the following.

### Verified Correct

- **Response matrix fill order.** `responses_full[etabin]->Fill(cluster_Et, particle_Pt, weight * response_reweight)` uses (reco, truth, weight), matching RooUnfold conventions.
- **Efficiency corrections.** Applied after unfolding in `CalculatePhotonYield.C`. The unfolded distribution is divided by `total_eff = eff_reco * eff_iso * eff_id * vertex_eff`. This is the correct methodology -- efficiencies must not be folded into the response matrix.
- **Purity correction.** The ABCD self-consistent equation with signal leakage (cB, cC, cD) is solved numerically via `myfunc()`. The leakage-corrected yield feeds into unfolding.
- **Response reweighting.** When enabled, adjusts MC response bins by the ratio of data-to-MC reco spectra. This is a valid technique for prior dependence studies.

### Historical Bug: Second Half Histogram Weighting (Fixed)

A code review found that the second-half marginal histograms (`h_pT_truth_secondhalf_response`, `h_pT_reco_secondhalf_response`) were incorrectly filled with `weight * response_reweight`, while the first-half marginals used plain `weight`. This has been fixed in the current code -- both halves now use plain `weight`.

The bug affected half-closure test validation but did not affect the main physics result, because the full-sample response matrix and the data unfolding path do not use these histograms. The half-sample response matrix itself was always correct.

## Closure.C Current State and Planned Extensions

The current `Closure.C` has a partial implementation: full closure for a single iteration with a stub comment for half closure. The function signature is:

```cpp
void Closure(const std::string &configname = "config.yaml")
```

Current limitations:
- Hardcoded response file path (does not read from config)
- Only runs 1 iteration (`niterations = 1`)
- Half closure test is a stub (`//half closure test` comment only)
- No chi2 calculation or diagnostic output ROOT file

A code review (documented in `CLOSURE_EXTENSION_SUMMARY.md`) proposed extending the macro with config-driven paths, multi-iteration loops, half closure implementation, chi2 diagnostics, and error handling. These extensions have not yet been applied to the code.

Note: `run_closure_test.sh` in `efficiencytool/` was written for the planned extended signature and will not work with the current `Closure.C`.

## Troubleshooting

### Cannot open response file
- Ensure `RecoEffCalculator_TTreeReader.C` has been run for all signal MC samples (photon5, photon10, photon20)
- Ensure `MergeSim.C` has been run to produce the merged response file
- Check that `var_type` in the config matches your analysis suffix
- Verify the `response_outfile` path in the config

### Cannot find response_matrix_full_X or response_matrix_half_X
- The response file exists but does not contain the expected objects
- Open the file in ROOT and run `.ls` to inspect contents
- Expected objects: `response_matrix_full_0`, `response_matrix_half_0`, truth/reco projections

### Large deviations in half closure
- May indicate insufficient MC statistics -- consider using higher-statistics samples
- Check for bugs in event selection logic between the two halves

### Chi2/ndf much larger than 1
- Errors may be underestimated in the unfolding
- Indicates poor closure -- look for systematic biases in the ratio plots
- May need a different iteration count or a different unfolding method
