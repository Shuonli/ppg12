# Unfolding Closure Tests

## Overview

This document describes the closure test implementation for validating the unfolding procedure in the photon isolation analysis.

## Issues Found in Current Code

### 1. **Bug in Second Half Response Matrix Weighting** (Line 2190-2191 in RecoEffCalculator_TTreeReader.C)

**Issue**: The second half histograms are filled with `response_reweight` applied, while the first half is not.

```cpp
// First half (lines 2183-2186) - CORRECT
if (ientry < (nentries / 2))
{
    h_pT_truth_half_response[etabin]->Fill(particle_Pt[iparticle], weight);
    h_pT_reco_half_response[etabin]->Fill(cluster_Et[icluster], weight);
    responses_half[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
    h_response_half_list[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
}
// Second half (lines 2188-2192) - BUG!
else
{
    h_pT_truth_secondhalf_response[etabin]->Fill(particle_Pt[iparticle], weight * response_reweight); // WRONG
    h_pT_reco_secondhalf_response[etabin]->Fill(cluster_Et[icluster], weight * response_reweight);   // WRONG
}
```

**Fix**: Remove `response_reweight` from lines 2190-2191:

```cpp
else
{
    h_pT_truth_secondhalf_response[etabin]->Fill(particle_Pt[iparticle], weight);           // FIXED
    h_pT_reco_secondhalf_response[etabin]->Fill(cluster_Et[icluster], weight);             // FIXED
}
```

**Impact**: This bug causes the second half marginal distributions to have incorrect normalization, which could affect half-closure test results. However, if the second half histograms are only used for pseudo-data (not for response matrix initialization), the impact is minimal since the pseudo-data is taken directly from the filled histograms.

### 2. **Response Matrix Filling Order** - CORRECT ✓

The response matrix is filled correctly:
```cpp
responses_full[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
```

Order: (reco, truth, weight) - This matches RooUnfold conventions.

### 3. **Efficiency Corrections** - CORRECT ✓

Efficiency corrections are applied AFTER unfolding (CalculatePhotonYield.C lines 902-916), which is the correct methodology. The unfolded distribution should be corrected for:
- Reconstruction efficiency
- Isolation efficiency
- ID efficiency (tight cuts)
- Vertex efficiency

### 4. **Reweighting Procedure** - VALID ✓

The optional reweighting (when config `reweight: 1`) adjusts the MC response matrix to match the data distribution. This is a valid technique to account for MC/data discrepancies.

## Closure Tests

### What are Closure Tests?

Closure tests validate that the unfolding procedure is unbiased and working correctly.

### Full Closure Test

**Procedure:**
1. Use MC reconstructed distribution as "pseudo-data"
2. Unfold using response matrix built from the same MC
3. Compare unfolded result to MC truth distribution
4. **Expected result**: Perfect closure (unfolded = truth) if unfolding is unbiased

**Interpretation:**
- Deviation from 1.0 indicates bias in the unfolding
- Should test multiple iterations to find optimal iteration number
- Chi2/ndf close to 1.0 indicates good closure

### Half Closure Test

**Procedure:**
1. Split MC sample into two independent halves (first 50% and second 50% of events)
2. Build response matrix from first half
3. Use second half reconstructed distribution as "pseudo-data"
4. Unfold second half with first half's response matrix
5. Compare unfolded result to second half truth

**Expected result:** Good closure even with statistically independent samples

**Interpretation:**
- Tests robustness against statistical fluctuations
- Deviation from 1.0 indicates:
  - Insufficient statistics
  - Bias from using different samples
  - Iteration-dependent bias
- More realistic test than full closure

## Running the Closure Tests

### Basic Usage

```bash
# Run with default settings (10 iterations, both full and half closure)
root -l -b -q 'Closure.C("config_bdt_none.yaml")'

# Run with custom number of iterations
root -l -b -q 'Closure.C("config_bdt_none.yaml", 5)'

# Run and save all iteration histograms to ROOT file
root -l -b -q 'Closure.C("config_bdt_none.yaml", 10, true)'

# Using the shell script wrapper
./run_closure_test.sh config_bdt_none.yaml 10 false
```

### Function Signature

```cpp
void Closure(const std::string &configname = "config_bdt_none.yaml",
             int niterations_max = 10,
             bool save_all_iterations = false)
```

**Parameters:**
- `configname`: Path to YAML configuration file
- `niterations_max`: Maximum number of Bayesian iterations to test (default: 10)
- `save_all_iterations`: Save all iteration histograms to output ROOT file (default: false)

### Requirements

The closure test macro requires:
1. Response matrix file (from `RecoEffCalculator_TTreeReader.C`)
2. Configuration file specifying eta bins and file paths

### Output

The macro produces:
1. **ROOT file**: `ClosureTest_[configname].root` containing:
   - Chi2 vs iteration plots for full and half closure
   - Comparison plots
   - Optional: All iteration histograms (if `save_all_iterations = true`)

2. **PNG plots** in `figure/` directory:
   - **Per iteration plots** (iterations 1, 2, 4):
     - `FullClosure_eta[N]_iter[I].png` - Full closure with ratio panel
     - `HalfClosure_eta[N]_iter[I].png` - Half closure with ratio panel

   - **Chi2 summary plots**:
     - `FullClosure_chi2_eta[N].png` - Chi2/ndf vs iteration for full closure
     - `HalfClosure_chi2_eta[N].png` - Chi2/ndf vs iteration for half closure
     - `ClosureComparison_eta[N].png` - Overlay of full and half closure chi2

Each iteration plot contains 2 panels:
- Top: Unfolded vs Truth distributions (with reco for reference)
- Bottom: Ratio unfolded/truth (target: ~1.0)

### Interpreting Results

**Good Closure:**
- Ratio plots centered at 1.0 across all pT bins
- RMS of ratio < 5% for most bins
- Mean ratio vs iteration stable after 2-4 iterations
- Chi2/ndf ≈ 1.0

**Poor Closure:**
- Systematic deviation from 1.0 (bias)
- Large fluctuations (insufficient statistics)
- Ratio trends with pT (bin migration issues)
- Strong iteration dependence (over-unfolding or under-unfolding)

**Optimal Iteration Selection:**
- Choose iteration where bias is minimized
- Typically 2-4 iterations for Bayesian unfolding
- Too few iterations: bias from regularization
- Too many iterations: amplified statistical fluctuations

## Systematic Studies

### Testing Different Configurations

The closure test can be run with different configurations to estimate systematic uncertainties:

```bash
# Nominal
root -l -b -q 'ClosureTest.C("config_nom.yaml", 0)'

# Energy scale variations
root -l -b -q 'ClosureTest.C("config_scale0.yaml", 0)'
root -l -b -q 'ClosureTest.C("config_scale10.yaml", 0)'

# Energy resolution variations
root -l -b -q 'ClosureTest.C("config_eres.yaml", 0)'

# Different BDT cuts
root -l -b -q 'ClosureTest.C("config_bdt_test.yaml", 0)'
```

### Recommended Checks

1. **Iteration stability**: Compare iterations 1-10 to find plateau region
2. **Bin migration**: Check 2D response matrix for off-diagonal elements
3. **Statistical tests**: Ensure chi2/ndf ≈ 1.0 for good closure
4. **Systematic variations**: Run with different configs to assess systematic uncertainty
5. **Reweighting impact**: Compare with/without reweighting (`reweight: 0` vs `reweight: 1`)

## Example Workflow

```bash
# 1. Generate response matrices (if not already done)
root -l -b -q 'RecoEffCalculator_TTreeReader.C("config_bdt_none.yaml", "photon5")'
root -l -b -q 'RecoEffCalculator_TTreeReader.C("config_bdt_none.yaml", "photon10")'
root -l -b -q 'RecoEffCalculator_TTreeReader.C("config_bdt_none.yaml", "photon20")'

# 2. Merge response matrices from different MC samples
root -l -b -q 'MergeSim.C("config_bdt_none.yaml")'

# 3. Run closure tests with 10 iterations
root -l -b -q 'Closure.C("config_bdt_none.yaml", 10, false)'

# Or using the shell script
./run_closure_test.sh config_bdt_none.yaml 10 false

# 4. Inspect output
root -l ClosureTest_config_bdt_none.root
# Look at plots in figure/ directory
ls -lh figure/*Closure*.png
```

## Troubleshooting

### "Cannot open response file"
- Ensure you've run `RecoEffCalculator_TTreeReader.C` first
- Check that `var_type` in config matches your analysis (e.g., "photon", "jet")
- Verify the `response_outfile` path in config is correct
- Example: config should have `MC_response` which becomes `MC_response_photon.root`

### "Cannot find response_matrix_full_X"
- The response file exists but doesn't contain expected histograms
- Verify that the response file contains the expected histograms:
  ```bash
  root -l MC_response_photon.root
  .ls
  ```
- Should see: `response_matrix_full_0`, `response_matrix_half_0`, etc.

### Large deviations in half-closure test
- May indicate insufficient statistics in MC
- Consider combining more MC samples or using higher statistics sample
- Check for bugs in event selection between first and second halves

### Chi2/ndf >> 1
- Under-estimated errors in unfolding
- Indicates poor closure - check for systematic biases
- May need different iteration number or different unfolding method

### All ratios systematically above/below 1.0
- Check normalization of pseudo-data and truth histograms
- Verify efficiency corrections are applied consistently
- Check for bugs in response matrix filling

## References

- RooUnfold documentation: http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html
- Bayesian unfolding: G. D'Agostini, NIM A 362 (1995) 487
- Closure tests: https://arxiv.org/abs/1910.14215

## Contact

For questions about the closure test implementation, check the sPHENIX photon working group wiki or contact the analysis team.
