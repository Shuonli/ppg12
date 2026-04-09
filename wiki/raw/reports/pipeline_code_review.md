# Pipeline Code Review Report

Date: 2026-04-08

This report documents a systematic review of the PPG12 analysis pipeline for cross-section weight consistency, isolation cut implementation, BDT threshold application, config-code sync, and pT bin consistency.

---

## 1. Cross-Section Weights Sync

### Status: Mostly Centralized

The shared header `efficiencytool/CrossSectionWeights.h` provides a `PPG12::GetSampleConfig()` function that is used by the three main pipeline macros:

- `efficiencytool/RecoEffCalculator_TTreeReader.C:96` -- uses `PPG12::GetSampleConfig(filetype)`
- `efficiencytool/ShowerShapeCheck.C:124` -- uses `PPG12::GetSampleConfig(filetype)`
- `efficiencytool/DoubleInteractionCheck.C:98` -- uses `PPG12::GetSampleConfig(filetype)`
- `efficiencytool/CalculatePhotonYield.C:5-6` -- includes header, uses `photon20cross` and `jet50cross` directly via `using namespace PPG12`

`MergeSim.C` does not include `CrossSectionWeights.h` at all, which is fine because it only merges already-weighted ROOT files (no cross-section logic needed).

### Remaining Inline Weight Issues

Several secondary macros still use inline if/else chains instead of `GetSampleConfig()`:

| File | Line | Issue |
|------|------|-------|
| `SaturationStudy.C` | 61-67 | Inline weights; `jet15.upper=20` vs header's 21; `jet20.upper=30` vs header's 32 |
| `BDTScoreVsET.C` | 143-149 | Inline weights; same jet15/jet20 boundary discrepancy |
| `embed_test/ShowerShapeCheck.C` | 111-165 | Inline weights (legacy copy, includes header but doesn't use `GetSampleConfig`) |
| `showershape_vertex_check.C` | 67-117 | Inline weights |
| `plot_cluster_time.C` | 119-149 | Inline weights |
| `PlotSpectra.C` | 47-59 | Inline weights |
| `TruthJetEff.C` | 58-76 | Inline weights; normalizes jets to `jet30cross` (line 76) instead of `jet50cross` |

### BDTinput.C Jet Reference Denominator

`FunWithxgboost/BDTinput.C:85,95,104,120` normalizes all jet sample weights relative to `jet30cross` instead of `jet50cross`. This only affects internal luminosity counting for training data assembly (not written to training files), so it is benign but violates the project convention.

### Jet pT Window Discrepancies (inline files only)

| Sample | `CrossSectionWeights.h` | `SaturationStudy.C` / `BDTScoreVsET.C` |
|--------|------------------------|-----------------------------------------|
| jet15 | lower=15, upper=**21** | lower=15, upper=**20** |
| jet20 | lower=**21**, upper=**32** | lower=**20**, upper=**30** |

Files using `GetSampleConfig()` are unaffected.

---

## 2. Isolation Cut Consistency

### Parametric Formula: CONSISTENT

All macros that apply the isolation cut use the same formula:

```cpp
float recoiso_max = recoiso_max_b + recoiso_max_s * clusterET;
```

Verified at:
- `RecoEffCalculator_TTreeReader.C:2139`
- `ShowerShapeCheck.C:1402` (via `recoisoET` processing at line 1396-1400 for MC shift/scale, then implicit via 2D histogram filling -- ShowerShapeCheck does NOT apply a boolean iso/noniso cut itself; it fills shower shape variables vs continuous `recoisoET` on the Y axis)
- `DoubleInteractionCheck.C:1192`
- `Cluster_rbr.C:421`
- `Cluster_rbr_TTreeReader.C:306`

All read `recoiso_max_b` and `recoiso_max_s` from config. No hardcoded isolation thresholds found in the main analysis chain.

### MC Isolation Correction: CONSISTENT

All three macros apply `recoisoET = recoisoET * mc_iso_scale + mc_iso_shift` for simulation, using config values.

---

## 3. BDT Threshold Consistency

### Tight BDT Min (parametric): CONSISTENT

All macros apply `tight_bdt_min_slope * ET + tight_bdt_min_intercept`:
- `RecoEffCalculator_TTreeReader.C:2228`
- `ShowerShapeCheck.C:1752`
- `DoubleInteractionCheck.C:921`
- `Cluster_rbr.C:493`

### Non-Tight BDT Max: **BUG IN ShowerShapeCheck.C**

**`ShowerShapeCheck.C:1787`** uses flat `non_tight_bdt_max` (value 0.5 from config):
```cpp
bdt_score < non_tight_bdt_max)
```

**`RecoEffCalculator_TTreeReader.C:2263`** uses parametric form:
```cpp
bdt_score < non_tight_bdt_max_slope * cluster_Et[icluster] + non_tight_bdt_max_intercept)
```

**`DoubleInteractionCheck.C:935`** uses parametric form (correct):
```cpp
bdt_score < non_tight_bdt_max_slope * clusterET + non_tight_bdt_max_intercept)
```

**`Cluster_rbr.C:528`** uses parametric form (correct).

With config values `bdt_max_slope=-0.015`, `bdt_max_intercept=0.80`:
- At ET=10: parametric gives 0.65, flat gives 0.50 (ShowerShapeCheck too restrictive by 0.15)
- At ET=20: parametric gives 0.50, flat gives 0.50 (agrees)
- At ET=30: parametric gives 0.35, flat gives 0.50 (ShowerShapeCheck too permissive by 0.15)

**Root cause:** `ShowerShapeCheck.C` reads `non_tight_bdt_max` (line 488) but never reads `non_tight_bdt_max_slope` or `non_tight_bdt_max_intercept`. The parametric config fields are silently ignored.

**Impact:** Shower shape distributions in the non-tight ABCD regions will have a different cluster population than the main analysis, affecting shower shape comparison plots and potentially purity studies.

---

## 4. Common Cuts Consistency

### `wr_cogx_bound`: **INCONSISTENT ACROSS MACROS**

- **`RecoEffCalculator_TTreeReader.C:2157`** applies `wr_cogx > common_wr_cogx_bound` in common cuts
- **`ShowerShapeCheck.C:1691`** has `wr_cogx` cut **commented out**: `// (!(wr_cogx < common_wr_cogx_bound && ...))`; instead applies only `cluster_weta_cogx < common_cluster_weta_cogx_bound`
- **`DoubleInteractionCheck.C:1206-1209`** does not include `wr_cogx` in common cuts at all

With `wr_cogx_bound: 0.0` in nominal config, this difference is benign (the RecoEffCalculator cut `wr_cogx > 0.0` is always true for physical values). However, if this parameter were changed to a nonzero value, the macros would diverge.

### `pscut` (extra common e32/et1 cut): Only in ShowerShapeCheck.C

`ShowerShapeCheck.C:1696-1701` defines an additional `pscut` gating on `common_e32_over_e35` and `common_et1` ranges. This does not appear in `RecoEffCalculator_TTreeReader.C` or `DoubleInteractionCheck.C`. Checking if it affects ABCD fills would require reading further into the ShowerShapeCheck fill logic.

---

## 5. pT Bin Consistency

### Current State: plotcommon.h and config_bdt_nom.yaml NOW MATCH

| Source | Edges |
|--------|-------|
| `plotcommon.h:6` | `{8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36}` |
| `config_bdt_nom.yaml:51` | `[8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]` |

These are consistent (12 bins, edges 8-36 GeV).

### CLAUDE.md is STALE

`CLAUDE.md:47` documents: `[8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35]`

This does not match the current code. The last two edges are `30, 35` in CLAUDE.md vs `32, 36` in both plotcommon.h and config. CLAUDE.md should be updated.

### Wiki constants-sync.md is STALE

`wiki/reference/constants-sync.md:47-48` claims plotcommon.h has `28, 30, 35` -- this is incorrect. Both plotcommon.h and config now have `28, 32, 36`.

### Python N_PT_BINS: FIXED

Both `make_selection_report.py:28` and `make_comparison_report.py:34` now use `N_PT_BINS = 12`, matching the current bin count.

---

## 6. Config-Code Sync

### Luminosity Values: 4+ DISTINCT VALUES

| Value (pb^-1) | Source |
|---------------|--------|
| 16.2735 | `config_bdt_nom.yaml:239` (actual YAML config used in analysis) |
| 16.6 | `plotcommon.h:21,25` (hardcoded legend strings) |
| 16.9 | `PPG12-analysis-note/selection.tex:135`, `introduction.tex:24`, `analysis.tex:264` |
| 49.562 | `CalculatePhotonYield.C:74` (fallback default, for all-runs period) |

The final cross-section is computed using the config value (16.2735), which is correct. Plot labels showing 16.6 are cosmetically inconsistent. The analysis note's 16.9 includes a different calculation methodology (per the note text).

### plotcommon.h Isolation Legend: MISLEADING

`plotcommon.h:24` legend string says:
```
E_T^{iso, R=0.3} < 4 GeV
```

But nominal analysis uses:
- Topo cluster iso R=0.4 (`use_topo_iso: 2` = `cluster_iso_topo_04`)
- Parametric cut: `reco_iso_max = 0.502 + 0.0433 * ET` (not flat 4 GeV)

The analysis note correctly states R=0.4 (`analysis.tex:6`).

### CalculatePhotonYield.C Hardcoded Jet Events

`CalculatePhotonYield.C:63` has `jetevents = 0.3555 * 1E7`. This magic number (35.55% of 10M) represents the fraction of jet50 events passing some selection, but is not derived from any config. It only affects the MC closure test luminosity (`isMC=true` path), not the data cross-section.

---

## 7. Summary of Issues by Severity

### HIGH (affects physics results if relevant config values change)

1. **ShowerShapeCheck.C non-tight BDT: flat vs parametric** (line 1787)
   - Uses `non_tight_bdt_max` instead of `slope * ET + intercept`
   - Currently differs from RecoEffCalculator at low and high ET

### MEDIUM (cosmetic/documentation, or latent bugs)

2. **plotcommon.h isolation legend** (line 24): says R=0.3 < 4 GeV, should be R=0.4, parametric
3. **CLAUDE.md pT bins** (line 47): says `28, 30, 35`, should be `28, 32, 36`
4. **wiki/reference/constants-sync.md pT bins** (lines 47-48): stale description
5. **ShowerShapeCheck.C `wr_cogx` common cut** (line 1691): commented out, while active in RecoEffCalculator (line 2157). Benign only because `wr_cogx_bound=0.0`.
6. **Luminosity scatter**: 16.2735 / 16.6 / 16.9 across config / plots / analysis note

### LOW (secondary macros, benign)

7. **6+ macros with inline cross-section weights** instead of `GetSampleConfig()` -- risk of drift
8. **SaturationStudy.C, BDTScoreVsET.C jet pT boundaries** differ from `CrossSectionWeights.h`
9. **BDTinput.C jet reference denominator** (`jet30cross` vs `jet50cross`)
10. **CalculatePhotonYield.C hardcoded `jetevents = 0.3555 * 1E7`** -- only affects MC closure

---

## 8. Recommendations

1. **Fix ShowerShapeCheck.C non-tight BDT**: Add `non_tight_bdt_max_slope` and `non_tight_bdt_max_intercept` reads (mirroring RecoEffCalculator_TTreeReader.C:479-480) and use `slope * ET + intercept` at line 1787.
2. **Update plotcommon.h strleg4**: Change to `R=0.4` and remove the flat "< 4 GeV" (or make it parametric-aware).
3. **Update CLAUDE.md pT bins**: Replace `[8,...,28,30,35]` with `[8,...,28,32,36]`.
4. **Migrate remaining inline-weight macros** to use `PPG12::GetSampleConfig()`.
5. **Reconcile luminosity values** or document that different values serve different purposes.
