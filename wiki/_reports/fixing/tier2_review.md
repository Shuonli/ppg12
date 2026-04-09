# Tier 2 Fix Review

Date: 2026-04-08
Reviewer: Claude Opus 4.6 (automated)

---

## Fix 1: plotcommon.h pT bins

### Verification Checklist

1. **ptRanges updated** -- CONFIRMED
   - `plotting/plotcommon.h` line 6: `const float ptRanges[NptBins + 1] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36};`
   - Last three entries are now `28, 32, 36` (previously `28, 30, 35`).

2. **NptBins unchanged** -- CONFIRMED
   - `plotting/plotcommon.h` line 5: `const int NptBins = 12;` (unchanged, correct: 13 edges - 1 = 12 bins).

3. **pTmax updated** -- CONFIRMED
   - `plotting/plotcommon.h` line 9: `const float pTmax = 36;` (previously 30).

4. **config_bdt_nom.yaml pT_bins match** -- CONFIRMED
   - `efficiencytool/config_bdt_nom.yaml` line 51: `pT_bins: [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]`
   - Exact match with `ptRanges` in plotcommon.h.

5. **No plotting/*.C files hardcode old values** -- CONFIRMED
   - Grep for `30, 35` and `30., 35.` in `plotting/*.C` returns zero matches.
   - All plotting macros reference `ptRanges[]` and `pTmax` from `plotcommon.h`, so they inherit the fix automatically.

### Observations (not blocking)

- Several `config_showershape_*.yaml` files still use `pT_bins_truth: [7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35]` (old values). These are shower-shape study configs with different (coarser) `pT_bins: [10, 14, 18, 22, 28, 30]`, so they may be intentionally different from the main analysis binning. However, the `28, 30, 35` truth bins here are inconsistent with the updated `28, 32, 36` in the main pipeline. This does not affect Fix 1 correctness but is worth noting for future cleanup.

- The `config_backup/config_bdt_nom.yaml` still has old truth bins `[7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35]`. This is expected for a backup.

### Rating: **CORRECT**

---

## Fix 2: Parametric BDT in DoubleInteractionCheck.C

### 2.1 Config Reads

**Tight BDT (lines 354-357):**
```cpp
float tight_bdt_max = configYaml["analysis"]["tight"]["bdt_max"].as<float>(1.0);
float tight_bdt_min = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0.0);
float tight_bdt_min_slope = configYaml["analysis"]["tight"]["bdt_min_slope"].as<float>(0);
float tight_bdt_min_intercept = configYaml["analysis"]["tight"]["bdt_min_intercept"].as<float>(tight_bdt_min);
```

**Non-tight BDT (lines 376-379):**
```cpp
float non_tight_bdt_max = configYaml["analysis"]["non_tight"]["bdt_max"].as<float>(1.0);
float non_tight_bdt_min = configYaml["analysis"]["non_tight"]["bdt_min"].as<float>(0.0);
float non_tight_bdt_max_slope = configYaml["analysis"]["non_tight"]["bdt_max_slope"].as<float>(0);
float non_tight_bdt_max_intercept = configYaml["analysis"]["non_tight"]["bdt_max_intercept"].as<float>(non_tight_bdt_max);
```

**Comparison with RecoEffCalculator_TTreeReader.C (lines 427-430, 477-480):**
- Tight: identical YAML paths, identical default values (slope=0, intercept falls back to `tight_bdt_min`). MATCH.
- Non-tight: identical YAML paths, identical default values (slope=0, intercept falls back to `non_tight_bdt_max`). MATCH.

**Backward compatibility:** When slope=0 and intercept defaults to the flat value, the parametric formula degenerates to the original flat threshold. CONFIRMED.

### 2.2 All BDT Threshold Application Sites

The fix centralizes all tight/non-tight classification into a single lambda `classifyTightNonTight` (defined at line 905). This lambda is called from exactly 3 locations:

| Call site | Line | Context | ET variable passed | Correct? |
|-----------|------|---------|--------------------|----------|
| Original kinematics | 1215 | Single-interaction, original vertex | `clusterET` (original) | YES |
| Double-interaction shifted | 1309-1310 | Toy double vertex shift | `clusterET_shifted` (recalculated) | YES |
| Smeared double-interaction | 1357-1358 | Double vertex + Gaussian smearing | `clusterET_sm` (recalculated) | YES |

Inside the lambda (line 921-922):
```cpp
float tight_bdt_min_et = tight_bdt_min_slope * clusterET + tight_bdt_min_intercept;
bool is_bdt_tight = (bdt_score > tight_bdt_min_et) && (bdt_score < tight_bdt_max);
```

The `clusterET` parameter in the lambda receives the correct ET for each context. CONFIRMED.

### 2.3 Non-tight Boundary

The non-tight BDT upper bound (line 935) now uses the parametric form:
```cpp
bdt_score < non_tight_bdt_max_slope * clusterET + non_tight_bdt_max_intercept
```

This matches RecoEffCalculator_TTreeReader.C line 2263:
```cpp
bdt_score < non_tight_bdt_max_slope * cluster_Et[icluster] + non_tight_bdt_max_intercept
```

CONFIRMED: both tight and non-tight BDT boundaries are now parametric in DoubleInteractionCheck.C.

### 2.4 Comparison Direction

- **Tight region**: `bdt_score > tight_bdt_min_et` (line 922) -- correct, tight = above threshold.
- **Non-tight region**: `bdt_score < non_tight_bdt_max_slope * clusterET + non_tight_bdt_max_intercept` (line 935) -- correct, non-tight = below the tight threshold.

Matches RecoEffCalculator_TTreeReader.C exactly (lines 2228-2231 for tight, 2262-2263 for non-tight).

### 2.5 Logical Consistency (No Overlap or Gap)

From `config_bdt_nom.yaml`:
- **Tight lower bound**: `intercept + slope * ET = 0.80 + (-0.015) * ET`
- **Non-tight upper bound**: `intercept + slope * ET = 0.80 + (-0.015) * ET`

These are the same function, so the boundary between tight and non-tight is seamless:
- Tight: `bdt_score > threshold(ET)` (strict greater-than)
- Non-tight: `bdt_score < threshold(ET)` (strict less-than)

Clusters with `bdt_score == threshold(ET)` fall into neither category, which is standard floating-point boundary behavior and consistent with RecoEffCalculator_TTreeReader.C.

No overlap. No gap (except the measure-zero boundary itself). CONFIRMED.

### 2.6 Full Comparison with RecoEffCalculator_TTreeReader.C

| Aspect | RecoEffCalculator_TTreeReader.C | DoubleInteractionCheck.C | Match? |
|--------|--------------------------------|--------------------------|--------|
| Tight slope config path | `analysis.tight.bdt_min_slope` | `analysis.tight.bdt_min_slope` | YES |
| Tight intercept config path | `analysis.tight.bdt_min_intercept` | `analysis.tight.bdt_min_intercept` | YES |
| Tight formula | `tight_bdt_min_slope * cluster_Et + tight_bdt_min_intercept` | `tight_bdt_min_slope * clusterET + tight_bdt_min_intercept` | YES |
| Non-tight slope config path | `analysis.non_tight.bdt_max_slope` | `analysis.non_tight.bdt_max_slope` | YES |
| Non-tight intercept config path | `analysis.non_tight.bdt_max_intercept` | `analysis.non_tight.bdt_max_intercept` | YES |
| Non-tight formula | `non_tight_bdt_max_slope * cluster_Et + non_tight_bdt_max_intercept` | `non_tight_bdt_max_slope * clusterET + non_tight_bdt_max_intercept` | YES |
| Default slope | 0 | 0 | YES |
| Default intercept (tight) | `tight_bdt_min` | `tight_bdt_min` | YES |
| Default intercept (non-tight) | `non_tight_bdt_max` | `non_tight_bdt_max` | YES |

All aspects match. CONFIRMED.

### Observation (not blocking, separate from this fix)

`ShowerShapeCheck.C` still uses a flat `non_tight_bdt_max` (line 1787) without slope/intercept for the non-tight upper bound. It does have the parametric tight lower bound (lines 440-441, 1752-1755). This is a pre-existing inconsistency in `ShowerShapeCheck.C` and is outside the scope of this tier 2 fix, but should be addressed separately.

### Rating: **CORRECT**

---

## Summary

| Fix | Rating | Notes |
|-----|--------|-------|
| Fix 1: plotcommon.h pT bins | **CORRECT** | ptRanges, pTmax, and config_bdt_nom.yaml all consistent. No hardcoded old values in plotting macros. |
| Fix 2: Parametric BDT in DoubleInteractionCheck.C | **CORRECT** | All config reads, formula applications, and comparison directions match RecoEffCalculator_TTreeReader.C exactly. All 3 call sites pass the correct ET variable. Backward-compatible defaults. Non-tight boundary also made parametric. |

### Items for Future Attention (not blocking)

1. `ShowerShapeCheck.C` non-tight BDT upper bound is still flat (line 1787) -- should be made parametric to match the other two analysis macros.
2. Several `config_showershape_*.yaml` files retain old `pT_bins_truth: [..., 28, 30, 35]` values. If these should track the main analysis binning, they need updating.
