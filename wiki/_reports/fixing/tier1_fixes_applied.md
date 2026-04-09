# Tier 1 Bug Fixes Applied

Date: 2026-04-08

## Fix 1: BDTinput.C copy-paste bug (line 194)

**File:** `FunWithxgboost/BDTinput.C`

The `et4_on` variable was reading from the wrong YAML key (`et3_on` instead of `et4_on`),
meaning the et4 feature toggle was always shadowing the et3 setting.

```
- int et4_on = configYaml["analysis"]["et3_on"].as<int>(1);
+ int et4_on = configYaml["analysis"]["et4_on"].as<int>(1);
```

## Fix 2: N_PT_BINS = 10 -> 12 in Python report scripts

**Files:** `plotting/make_selection_report.py`, `plotting/make_comparison_report.py`

`NptBins` in `plotcommon.h` is 12 (pT range 8-35 GeV, 12 bins). Both Python scripts
had `N_PT_BINS = 10`, causing them to skip the last two pT bins when generating reports.

```
- N_PT_BINS = 10
+ N_PT_BINS = 12
```

## Fix 3: npb_score_cut YAML path standardization

**Files:**
- `efficiencytool/DoubleInteractionCheck.C` (line 65)
- `efficiencytool/NPB_PurityStudy.C` (line 227)
- `efficiencytool/embed_test/ShowerShapeCheck.C` (line 57)

These files read `configYaml["analysis"]["npb_score_cut"]` (flat path), but the
canonical location used by `RecoEffCalculator_TTreeReader.C` and `ShowerShapeCheck.C`
is `configYaml["analysis"]["common"]["npb_score_cut"]` (nested under `common`).
Without this fix, these tools would silently fall back to the default value (0.5)
instead of reading the config, breaking systematic variations that modify this cut.

```
- configYaml["analysis"]["npb_score_cut"]
+ configYaml["analysis"]["common"]["npb_score_cut"]
```

## Verification

All six edits verified by reading back the modified lines after applying.
