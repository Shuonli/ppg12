# Tier 2 Fixes Applied

Date: 2026-04-08

## Fix 1: Update plotcommon.h pT bins to match config

**File:** `plotting/plotcommon.h`

**Problem:** The last two pT bin edges in `ptRanges` were `30, 35`, but `config_bdt_nom.yaml` uses `pT_bins: [8,10,12,14,16,18,20,22,24,26,28,32,36]`. Additionally, `pTmax` was set to 30 instead of 36.

**Changes:**
- Line 6: `ptRanges` last two entries changed from `30, 35` to `32, 36`
- Line 9: `pTmax` changed from `30` to `36`

**Before:**
```cpp
const float ptRanges[NptBins + 1] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35};
const float pTmax = 30;
```

**After:**
```cpp
const float ptRanges[NptBins + 1] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36};
const float pTmax = 36;
```

`NptBins = 12` and `pTmin = 10` remain unchanged.

---

## Fix 2: Add parametric BDT threshold to DoubleInteractionCheck.C

**File:** `efficiencytool/DoubleInteractionCheck.C`

**Problem:** The tight BDT cut used a flat `tight_bdt_min` threshold, while `RecoEffCalculator_TTreeReader.C` uses an ET-dependent parametric form `intercept + slope * ET`. Similarly, the non-tight BDT upper bound was flat instead of parametric.

**Changes (config reading, lines 356-357):**
Added after `tight_bdt_min`:
```cpp
float tight_bdt_min_slope = configYaml["analysis"]["tight"]["bdt_min_slope"].as<float>(0);
float tight_bdt_min_intercept = configYaml["analysis"]["tight"]["bdt_min_intercept"].as<float>(tight_bdt_min);
```

**Changes (config reading, lines 378-379):**
Added after `non_tight_bdt_min`:
```cpp
float non_tight_bdt_max_slope = configYaml["analysis"]["non_tight"]["bdt_max_slope"].as<float>(0);
float non_tight_bdt_max_intercept = configYaml["analysis"]["non_tight"]["bdt_max_intercept"].as<float>(non_tight_bdt_max);
```

**Changes (classification lambda, lines 921-922):**
Replaced flat tight BDT check:
```cpp
// Before:
bool is_bdt_tight = (bdt_score > tight_bdt_min) && (bdt_score < tight_bdt_max);

// After:
float tight_bdt_min_et = tight_bdt_min_slope * clusterET + tight_bdt_min_intercept;
bool is_bdt_tight = (bdt_score > tight_bdt_min_et) && (bdt_score < tight_bdt_max);
```

**Changes (non-tight region, line 935):**
Replaced flat non-tight BDT upper bound:
```cpp
// Before:
bdt_score > non_tight_bdt_min && bdt_score < non_tight_bdt_max)

// After:
bdt_score > non_tight_bdt_min && bdt_score < non_tight_bdt_max_slope * clusterET + non_tight_bdt_max_intercept)
```

**Scope:** The `classifyTightNonTight` lambda is called from 3 locations (lines 1210, 1304, 1352 -- original, double-shifted, and smeared kinematics). All three call sites pass the appropriate `clusterET` to the lambda, so all are now correctly using the parametric threshold without any additional changes.

**Backward compatibility:** Default values for slope (0) and intercept (falls back to the flat value) ensure that configs without the new fields behave identically to before.
