# Tier 1 Bug Fix Review

Date: 2026-04-08
Reviewer: code-reviewer agent

---

## Fix 1: BDTinput.C et4_on copy-paste bug

**Rating: CORRECT**

Verified at `FunWithxgboost/BDTinput.C` line 194:
```cpp
int et4_on = configYaml["analysis"]["et4_on"].as<int>(1);
```

The git diff confirms exactly one line changed (`et3_on` -> `et4_on`). Grepped the entire file for `et3_on` and `et4_on` -- the only occurrences are:
- Line 193: `int et3_on = configYaml["analysis"]["et3_on"].as<int>(1);` (reads et3, correct)
- Line 194: `int et4_on = configYaml["analysis"]["et4_on"].as<int>(1);` (reads et4, now correct)
- Line 830: `nfail += et3_on;` (uses et3_on variable, correct)
- Line 842: `nfail += et4_on;` (uses et4_on variable, correct)

No other `et3_on` instances that should be `et4_on`. Fix is minimal and correct.

---

## Fix 2: N_PT_BINS in Python report scripts

**Rating: CORRECT**

Verified current state:
- `plotting/make_selection_report.py` line 28: `N_PT_BINS = 12  # noqa: matches NptBins in plotcommon.h`
- `plotting/make_comparison_report.py` line 34: `N_PT_BINS = 12`

Grepped all `*.py` files in the repo for `N_PT_BINS` -- only these two files define or use this constant. Both are set to 12, matching `NptBins` in `plotcommon.h`.

Note: These files are untracked by git (new files), so there is no git diff to confirm the before/after. The fix report states the change was from 10 to 12, and the current state is 12, which is correct.

---

## Fix 3: npb_score_cut YAML path standardization

**Rating: CORRECT**

Verified all three files now use the nested path:

1. `efficiencytool/DoubleInteractionCheck.C` line 65:
   ```cpp
   float npb_score_cut = configYaml["analysis"]["common"]["npb_score_cut"].as<float>(0.5);
   ```
   Git diff confirms change from `["analysis"]["npb_score_cut"]` to `["analysis"]["common"]["npb_score_cut"]`. Default fallback `0.5` preserved.

2. `efficiencytool/NPB_PurityStudy.C` line 227:
   ```cpp
   const float npb_score_cut = configYaml["analysis"]["common"]["npb_score_cut"].as<float>(0.5);
   ```
   Git diff confirms the same path fix. Default fallback `0.5` preserved.

3. `efficiencytool/embed_test/ShowerShapeCheck.C` line 57:
   ```cpp
   float npb_score_cut = configYaml["analysis"]["common"]["npb_score_cut"].as<float>(0.5);
   ```
   This file is untracked (no git baseline), but its current state uses the correct nested path with `0.5` fallback.

**Repo-wide grep for flat path `["analysis"]["npb_score_cut"]` (without "common"):**
Searched all `*.C` and `*.h` files in `efficiencytool/`, `FunWithxgboost/`, and `plotting/` -- zero matches found. The only remaining hits for the flat path are in the wiki report markdown files (documentation, not code). All code files consistently use `["analysis"]["common"]["npb_score_cut"]`.

---

## Additional Changes (beyond Tier 1 scope)

The fixer agent applied changes beyond the three requested Tier 1 fixes. These are NOT bugs in the fixes, but are additional refactoring that should be reviewed separately:

1. **CrossSectionWeights.h** (+92 lines): Added `SampleConfig` struct and `GetSampleConfig()` lookup function to centralize the per-sample if/else chains.

2. **DoubleInteractionCheck.C** (-89 net lines): Beyond the npb_score_cut fix, replaced the ~90-line if/else sample config block with a call to `GetSampleConfig()`. Semantic equivalence verified: for `filetype == "data"`, `issim = false` so `max_photon_*` / `max_jet_*` / `cluster_ET_upper` are never used; default values in `SampleConfig` (`weight=1.0`, `isbackground=false`) match the original defaults.

3. **RecoEffCalculator_TTreeReader.C** (-191/+92 lines): Multiple changes:
   - Added `mix_weight` and `vtxscan_sim_override` function parameters
   - Added `_nom` alias handling for input file paths
   - Replaced `TEfficiency::SetWeight()` + `Fill()` with `SetUseWeightedEvents()` + `FillWeighted()`
   - Replaced sample if/else chain with `GetSampleConfig()`
   - Added R=0.4 jet branch support via config flag

4. **ShowerShapeCheck.C** (-146/+12 lines): Replaced sample if/else chain with `GetSampleConfig()`.

5. **plotting/plot_final_selection.C** (+27/-6 lines): Changed hardcoded luminosity to read from YAML config at runtime with fallback.

6. **CLAUDE.md** (+4 lines): Added wiki reference.

These additional changes appear to be beneficial refactoring (DRY principle for cross-section weights, proper weighted efficiency handling), but they go beyond the Tier 1 bug fix scope and were not requested.

---

## Summary

| Fix | File(s) | Rating |
|-----|---------|--------|
| 1. et4_on copy-paste | BDTinput.C | CORRECT |
| 2. N_PT_BINS = 12 | make_selection_report.py, make_comparison_report.py | CORRECT |
| 3. npb_score_cut path | DoubleInteractionCheck.C, NPB_PurityStudy.C, embed_test/ShowerShapeCheck.C | CORRECT |

All three Tier 1 fixes are verified correct. The fixer agent also performed additional refactoring across 5 other files/locations that should be reviewed as a separate change set.
