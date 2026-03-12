# Closure.C Extension Summary

## What Was Extended

The existing `Closure.C` file has been comprehensively extended to provide full closure test validation for the unfolding procedure.

## Changes Made

### 1. **Extended Function Signature**
```cpp
// OLD
void Closure(const std::string &configname = "config.yaml")

// NEW
void Closure(const std::string &configname = "config_bdt_none.yaml",
             int niterations_max = 10,
             bool save_all_iterations = false)
```

**New Parameters:**
- `niterations_max`: Test multiple iterations (1-10) instead of just 1
- `save_all_iterations`: Optionally save all iteration histograms to ROOT file

### 2. **Configuration File Handling**
```cpp
// OLD - Hardcoded path
std::string infilename = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_response_photon10.root";

// NEW - Config-driven
std::string var_type = configYaml["output"]["var_type"].as<std::string>();
std::string infilename = configYaml["output"]["response_outfile"].as<std::string>() + "_" + var_type + ".root";
```

### 3. **Full Closure Test - Enhanced**

**Added:**
- ✅ Chi2 calculation for each iteration
- ✅ Chi2/ndf printed to console
- ✅ All 10 iterations tested (previously only 1)
- ✅ Plots for iterations 1, 2, and 4 (key milestones)
- ✅ Chi2 vs iteration plot
- ✅ Improved ratio range (0.8-1.2) for better visibility

**Kept:**
- ✅ Original plotting style with `draw_1D_multiple_plot_ratio`
- ✅ sPHENIX style formatting

### 4. **Half Closure Test - Completed**

**Previously:** Only a comment `//half closure test` (line 108)

**Now:** Fully implemented with:
- ✅ Unfold second half using first half response matrix
- ✅ Chi2 calculation against second half truth
- ✅ All 10 iterations tested
- ✅ Plots for iterations 1, 2, and 4
- ✅ Chi2 vs iteration plot
- ✅ Same plotting style as full closure

### 5. **New Diagnostic Plots**

**Per Eta Bin:**
1. **Full Closure Plots** (`FullClosure_eta[N]_iter[I].png`)
   - Unfolded vs truth comparison with ratio panel
   - Iterations: 1, 2, 4

2. **Half Closure Plots** (`HalfClosure_eta[N]_iter[I].png`)
   - Second half unfolded vs second half truth
   - Uses first half response matrix
   - Iterations: 1, 2, 4

3. **Chi2 Evolution Plots**:
   - `FullClosure_chi2_eta[N].png` - Chi2/ndf vs iteration (full closure)
   - `HalfClosure_chi2_eta[N].png` - Chi2/ndf vs iteration (half closure)
   - `ClosureComparison_eta[N].png` - Overlay of both tests

### 6. **Output Organization**

**ROOT File** (`ClosureTest_[config].root`):
- Chi2 vs iteration graphs
- Comparison plots
- Optional: All iteration histograms

**Plot Directory** (`figure/`):
- Auto-created if doesn't exist
- All PNG plots organized by test type

### 7. **Enhanced Diagnostics**

**Console Output:**
```
========================================
UNFOLDING CLOSURE TEST
========================================
Configuration: config_bdt_none.yaml
Response file: MC_response_photon.root
Max iterations: 10

Processing eta bin 0: [-0.7, 0.7]
--- FULL CLOSURE TEST ---
  Iteration 1: chi2/ndf = 12.34/8 = 1.54
  Iteration 2: chi2/ndf = 8.76/8 = 1.10
  ...
--- HALF CLOSURE TEST ---
  Iteration 1: chi2/ndf = 15.23/8 = 1.90
  Iteration 2: chi2/ndf = 10.45/8 = 1.31
  ...
```

### 8. **Error Handling**

**Added:**
- Check if response file exists
- Check if response matrices loaded successfully
- Check if histograms loaded successfully
- Informative error messages with file paths

## Key Features Preserved

✅ **Original plotting style** - Uses `draw_1D_multiple_plot_ratio` from `draw.C`
✅ **ATLAS/sPHENIX style** - Maintains SetAtlasStyle() formatting
✅ **Same color scheme** - Uses existing color and marker style vectors
✅ **Backward compatible** - Can still call with just `Closure("config.yaml")`

## Usage Examples

### Basic (10 iterations, both tests)
```bash
root -l -b -q 'Closure.C("config_bdt_none.yaml")'
```

### Custom iterations
```bash
root -l -b -q 'Closure.C("config_bdt_none.yaml", 5)'
```

### Save all iterations
```bash
root -l -b -q 'Closure.C("config_bdt_none.yaml", 10, true)'
```

### Using shell script
```bash
./run_closure_test.sh config_bdt_none.yaml 10 false
```

## Output Files

**Before Extension:**
- `figure/h_unfold_reco_full_1.png` (single iteration only)

**After Extension:**
- `ClosureTest_config_bdt_none.root` (comprehensive results)
- `figure/FullClosure_eta0_iter1.png`
- `figure/FullClosure_eta0_iter2.png`
- `figure/FullClosure_eta0_iter4.png`
- `figure/HalfClosure_eta0_iter1.png`
- `figure/HalfClosure_eta0_iter2.png`
- `figure/HalfClosure_eta0_iter4.png`
- `figure/FullClosure_chi2_eta0.png`
- `figure/HalfClosure_chi2_eta0.png`
- `figure/ClosureComparison_eta0.png`

## Benefits

1. **Iteration Selection** - Can now systematically choose optimal iteration based on chi2
2. **Statistical Independence** - Half closure test validates against independent sample
3. **Complete Validation** - Both full and half closure for comprehensive testing
4. **Publication Ready** - Professional plots with chi2 statistics
5. **Systematic Studies** - Easy to run for multiple configs
6. **Config Driven** - No hardcoded paths, flexible for different analyses

## Files Modified

1. ✅ **Closure.C** - Extended from ~110 lines to ~353 lines
2. ✅ **run_closure_test.sh** - Updated for new parameters
3. ✅ **CLOSURE_TEST_README.md** - Updated documentation

## Files Deleted

- ❌ **ClosureTest.C** - Standalone version no longer needed (functionality merged into Closure.C)

## Next Steps

1. Run closure test on nominal configuration:
   ```bash
   ./run_closure_test.sh config_bdt_none.yaml
   ```

2. Review chi2 plots to select optimal iteration

3. Run for systematic variations:
   ```bash
   for config in config_*.yaml; do
       root -l -b -q "Closure.C(\"$config\", 10, false)"
   done
   ```

4. Document iteration selection in analysis note

---

**Date**: 2025-12-26
**Extension by**: Claude (Anthropic)
**Original File**: Closure.C (partial implementation)
**Extended File**: Closure.C (complete closure test framework)
