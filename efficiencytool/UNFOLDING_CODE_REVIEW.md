# Unfolding Code Review and Validation

## Executive Summary

I reviewed the unfolding implementation in `efficiencytool/` and found:
- ✅ **Overall methodology is CORRECT** - Bayesian unfolding with RooUnfold
- ✅ **Response matrix construction is CORRECT** - Proper (reco, truth, weight) ordering
- ✅ **Efficiency corrections are CORRECT** - Applied after unfolding
- ✅ **Reweighting procedure is VALID** - Legitimate technique for MC/data matching
- ⚠️ **ONE BUG FOUND** - Incorrect weighting in second half histograms (see details below)

## Code Review Details

### 1. Response Matrix Construction ✅

**File**: `RecoEffCalculator_TTreeReader.C`

**Lines 958, 982**: Response matrix initialization
```cpp
responses_full.push_back(new RooUnfoldResponse(
    (const TH1 *)h_pT_reco_response[ieta],      // Measured histogram
    (const TH1 *)h_pT_truth_response[ieta],     // Truth histogram
    h_response_full,                             // 2D response matrix
    Form("response_matrix_full_%d", ieta), "", false));
```
**Status**: ✅ CORRECT - Proper argument order for RooUnfoldResponse

**Line 2179**: Response matrix filling
```cpp
responses_full[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
```
**Status**: ✅ CORRECT - Fill(reco, truth, weight) matches RooUnfold conventions

### 2. Unfolding Procedure ✅

**File**: `CalculatePhotonYield.C`

**Lines 825-891**: Bayesian unfolding implementation
```cpp
for (int i = 0; i < niterations_total; i++)
{
    RooUnfoldBayes *unfold = new RooUnfoldBayes(response_reweighted, h_data_sub, i + 1);
    TH1D *h_unfold_sub = (TH1D *)unfold->Hunfold()->Clone(Form("h_unfold_sub_%d", i + 1));
    // ... scale by luminosity
    h_unfold_sub_list.push_back(h_unfold_sub);
}
```
**Status**: ✅ CORRECT - Standard Bayesian iterative unfolding

**Configuration**:
- Default iterations: 10 (result selected from iteration 2)
- Supports both standard and leakage-corrected unfolding
- Multiple iterations saved for systematic studies

### 3. Efficiency Corrections ✅

**Lines 902-918**: Efficiency correction after unfolding
```cpp
float total_eff = eff_reco_val * eff_iso_val * eff_id_val * vertex_eff_val;
h_unfold_sub->SetBinContent(ibin, h_unfold_sub->GetBinContent(ibin) / total_eff);
```
**Status**: ✅ CORRECT - Efficiencies applied AFTER unfolding (proper methodology)

**Efficiency components**:
- Reconstruction efficiency (matching cluster to truth)
- Isolation efficiency (isolation cut efficiency)
- ID efficiency (tight photon selection)
- Vertex efficiency (MBD trigger vertex cut)

### 4. Purity Correction ✅

**Lines 304-524**: Monte Carlo purity extraction
```cpp
// 2x2 sideband method
// A = tight + iso (signal)
// B = tight + noniso (leakage)
// C = nontight + iso (contamination)
// D = nontight + noniso (background)
// Purity = A / (A - R * B * C / D) where R is correlation factor
```
**Status**: ✅ CORRECT - Standard ABCD sideband method with MC sampling

### 5. Response Matrix Reweighting ✅

**Lines 745-821**: Optional reweighting to match data
```cpp
if (reweight == 1)
{
    float reweight_val = data_sub_val / reco_val * normcount;
    // Apply to each bin of response matrix
    h_response_reweighted->SetBinContent(ibin, jbin, response_val * reweight_val);
}
```
**Status**: ✅ VALID - Legitimate technique to account for MC/data shape differences

**Purpose**: Adjust MC response matrix when MC doesn't perfectly model data distribution

### 6. 🐛 BUG FOUND: Second Half Histogram Weighting ⚠️

**File**: `RecoEffCalculator_TTreeReader.C`
**Lines**: 2181-2192

**Issue**: Inconsistent weighting between first and second halves

**Current Code** (INCORRECT):
```cpp
if (ientry < (nentries / 2))
{
    // First half - CORRECT
    h_pT_truth_half_response[etabin]->Fill(particle_Pt[iparticle], weight);
    h_pT_reco_half_response[etabin]->Fill(cluster_Et[icluster], weight);
    responses_half[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
}
else
{
    // Second half - BUG HERE!
    h_pT_truth_secondhalf_response[etabin]->Fill(particle_Pt[iparticle], weight * response_reweight); // ❌
    h_pT_reco_secondhalf_response[etabin]->Fill(cluster_Et[icluster], weight * response_reweight);   // ❌
}
```

**Fixed Code**:
```cpp
else
{
    // Second half - FIXED
    h_pT_truth_secondhalf_response[etabin]->Fill(particle_Pt[iparticle], weight);           // ✅
    h_pT_reco_secondhalf_response[etabin]->Fill(cluster_Et[icluster], weight);             // ✅
}
```

**Impact**:
- The second half marginal distributions have incorrect normalization
- This affects half-closure test validation but NOT the main analysis
- The 2D response matrix for the half-sample is filled correctly (line 2185)
- If second half is used as pseudo-data, impact is minimal since filled histograms are used directly

**Severity**: LOW (affects validation only, not physics results)

**Recommendation**: Fix for consistency and correct closure tests

## Validation: Closure Tests

I've created comprehensive closure test code to validate the unfolding:

### Created Files:

1. **`ClosureTest.C`** - Main closure test macro
   - Tests both full and half closure
   - Produces diagnostic plots and chi2 statistics
   - Saves results for all iterations

2. **`CLOSURE_TEST_README.md`** - Complete documentation
   - How to run closure tests
   - How to interpret results
   - Troubleshooting guide

3. **`run_closure_test.sh`** - Simple runner script
   - Convenient execution wrapper
   - Error checking

### Running Closure Tests:

```bash
cd /gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool

# Method 1: Using the shell script
./run_closure_test.sh config_nom.yaml 0

# Method 2: Direct ROOT execution
root -l -b -q 'ClosureTest.C("config_nom.yaml", 0)'
```

### Test Types:
- `0` = Both full and half closure (default)
- `1` = Full closure only
- `2` = Half closure only

### Output:
- `ClosureTest_config_nom.root` - All histograms and ratios
- `FullClosure_eta[N].pdf` - Full closure plots per eta bin
- `HalfClosure_eta[N].pdf` - Half closure plots per eta bin

### Expected Results for Good Closure:

✅ **Full Closure Test**:
- Unfolded / Truth ratio centered at 1.0 (± 1-2%)
- Chi2/ndf ≈ 1.0
- Minimal iteration dependence after 2-4 iterations

✅ **Half Closure Test**:
- Unfolded / Truth ratio centered at 1.0 (± 3-5%)
- Tests statistical independence
- More realistic than full closure

## Recommendations

### Immediate Actions:

1. **Fix the bug** in `RecoEffCalculator_TTreeReader.C` lines 2190-2191
   - Remove `* response_reweight` from second half histogram fills
   - This ensures consistent weighting for closure tests

2. **Run closure tests** on current configuration
   ```bash
   ./run_closure_test.sh config_nom.yaml 0
   ```

3. **Validate iteration selection**
   - Check closure vs iteration plots
   - Confirm iteration 2 is optimal (current default)

### Long-term Validation:

4. **Systematic studies**
   - Run closure tests for all systematic variations (config_*.yaml)
   - Document closure quality for each variation

5. **Statistical tests**
   - Check chi2/ndf across all configurations
   - Ensure closure is robust against systematic variations

6. **Bin migration studies**
   - Plot 2D response matrices to check for off-diagonal elements
   - Ensure bin edges are appropriate for detector resolution

## Code Quality Assessment

### Strengths:
- ✅ Proper use of RooUnfold framework
- ✅ Correct Bayesian iterative unfolding
- ✅ Comprehensive efficiency corrections
- ✅ Sophisticated purity extraction with MC sampling
- ✅ Flexible configuration system with YAML
- ✅ Multiple systematic variations supported
- ✅ Half-sample infrastructure already exists

### Areas for Improvement:
- ⚠️ Fix weighting bug in second half histograms
- 📝 Add more inline documentation
- 📝 Document iteration selection procedure
- ✅ Closure tests now implemented (this work)

## Conclusion

The unfolding code is **fundamentally correct** and uses proper methodology. The one bug found is minor and only affects validation (not physics results). With the closure test framework now implemented, you can:

1. Validate the unfolding is unbiased
2. Select optimal iteration number
3. Estimate systematic uncertainties
4. Document closure quality for publications

The implementation is production-ready after fixing the minor weighting bug.

---

**Reviewer**: Claude (Anthropic)
**Date**: 2025-12-26
**Files Reviewed**:
- `RecoEffCalculator_TTreeReader.C` (2304 lines)
- `CalculatePhotonYield.C` (1056 lines)
- `RecoEffCalculator.C` (2229 lines)
- Configuration files (config_*.yaml)

**Files Created**:
- `ClosureTest.C` - Closure test implementation
- `CLOSURE_TEST_README.md` - Documentation
- `UNFOLDING_CODE_REVIEW.md` - This review
- `run_closure_test.sh` - Runner script
