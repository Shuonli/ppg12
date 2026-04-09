# Wiki Critic Review -- Phase 3

Reviewed: 2026-04-08
Reviewer: Claude Opus 4.6 (wiki critic agent)
Scope: All 19 wiki articles across index, pipeline, concepts, reference, and guides

---

## 1. Factual Spot-Checks (21 claims verified)

### CONFIRMED

1. **Claim (pipeline/01-tree-making.md):** CaloAna24.cc is 2321 lines, CaloAna24.h is 385 lines.
   **Code:** `wc -l` confirms exactly 2321 and 385 lines.
   **Verdict:** CONFIRMED

2. **Claim (pipeline/03-bdt-application.md):** apply_BDT.C is 492 lines, signature has default `filetype = "jet12_double"`.
   **Code:** `wc -l` = 492; signature matches exactly.
   **Verdict:** CONFIRMED

3. **Claim (pipeline/03-bdt-application.md):** ET threshold for BDT scoring is `cluster_Et > 7 GeV`.
   **Code:** Line 396: `if (cluster_Et_BDT > 7)`.
   **Verdict:** CONFIRMED

4. **Claim (pipeline/03-bdt-application.md):** 11 model variants with exact feature vectors listed.
   **Code:** Lines 284-391 of apply_BDT.C match every listed feature vector exactly.
   **Verdict:** CONFIRMED

5. **Claim (pipeline/03-bdt-application.md):** NPB score feature vector is 25 features in the listed order.
   **Code:** Lines 451-477 of apply_BDT.C match the wiki's listing exactly.
   **Verdict:** CONFIRMED

6. **Claim (pipeline/04-efficiency-yield.md):** RecoEffCalculator_TTreeReader.C is 2719 lines; signature matches.
   **Code:** `wc -l` = 2719; signature at line 44 matches.
   **Verdict:** CONFIRMED

7. **Claim (pipeline/04-efficiency-yield.md):** MergeSim.C uses exactly 6 jet samples: jet5, jet8, jet12, jet20, jet30, jet40.
   **Code:** Lines 21-26 and 59-64 of MergeSim.C confirm exactly these 6 samples.
   **Verdict:** CONFIRMED

8. **Claim (pipeline/04-efficiency-yield.md):** CalculatePhotonYield constants `mbdcorr = 25.2/42/0.57`, `nsimevents = 1E7`, `jetevents = 0.3555*1E7`.
   **Code:** Lines 54, 57, 63 of CalculatePhotonYield.C match exactly.
   **Verdict:** CONFIRMED

9. **Claim (reference/constants-sync.md):** Cross-section values in CrossSectionWeights.h -- photon5cross=146359.3, jet50cross=7.3113, etc.
   **Code:** Lines 13-27 of CrossSectionWeights.h match every listed value.
   **Verdict:** CONFIRMED

10. **Claim (pipeline/05-plotting-systematics.md):** plotcommon.h has `NptBins=12`, `ptRanges = {8,10,...,28,30,35}`.
    **Code:** Lines 5-6 of plotcommon.h match exactly.
    **Verdict:** CONFIRMED

11. **Claim (pipeline/05-plotting-systematics.md):** `N_PT_BINS = 10` hardcoded in Python report scripts.
    **Code:** make_selection_report.py line 28: `N_PT_BINS = 10`; make_comparison_report.py line 34: `N_PT_BINS = 10`.
    **Verdict:** CONFIRMED

12. **Claim (pipeline/01-tree-making.md):** RawClusterBuilderTemplate threshold is 0.070 GeV.
    **Code:** Fun4All_runDST.C line 219: `ClusterBuilder->set_threshold_energy(0.070)`.
    **Verdict:** CONFIRMED

13. **Claim (reference/config-schema.md):** Nominal config values: `reco_iso_max_b=0.502095`, `reco_iso_max_s=0.0433036`, `use_topo_iso=2`, `trigger_used=30`, `lumi=16.2735`.
    **Code:** config_bdt_nom.yaml lines 33-34, 53, 50, 238-239 -- all match.
    **Verdict:** CONFIRMED

14. **Claim (concepts/abcd-method.md):** fit_option=0 is Erf, fit_option=1 is Pade.
    **Code:** CalculatePhotonYield.C line 585: default is Erf (`TMath::Erf`); line 587: `if(fitoption==1)` switches to Pade.
    **Verdict:** CONFIRMED

15. **Claim (reference/constants-sync.md):** BDTinput.C normalizes jet weights to jet30cross (not jet50cross).
    **Code:** BDTinput.C lines 85, 95, 104, 120 -- all divide by `jet30cross`.
    **Verdict:** CONFIRMED

16. **Claim (pipeline/04-efficiency-yield.md):** mc_iso_scale code default is 1.0, mc_iso_shift code default is 0.0.
    **Code:** RecoEffCalculator_TTreeReader.C lines 367-368: `.as<float>(0.0)` and `.as<float>(1.0)`.
    **Verdict:** CONFIRMED

17. **Claim (reference/config-schema.md):** n_nt_fail nominal config value is 0, code default is 1.
    **Code:** config_bdt_nom.yaml line 57: `n_nt_fail: 0`; RecoEffCalculator line 282: `.as<int>(1)`.
    **Verdict:** CONFIRMED

### WRONG or OUTDATED

18. **Claim (concepts/systematic-variations.md):** Variant names include `iso_up`, `iso_down`, `ntf_up`, `ntf_down`, `fit_up`, `fit_down`, `fudge_up`, `fudge_down`, `mbd_eff_up`, `mbd_eff_down`, `unfold_nor`.
    **Code:** NONE of these variant names exist in make_bdt_variations.py. The actual names are completely different: `noniso04`, `noniso10`, `ntbdtmin05`, `ntbdtmin10`, `purity_pade`, `mciso_noscaleshift`, etc.
    **Verdict:** WRONG -- variant names are fabricated, not from the actual code

19. **Claim (concepts/systematic-variations.md):** Energy scale variants are `energyscale24` (0.024) and `energyscale26` (0.026); energy resolution variants are `energyresolution02` (0.02), `energyresolution04` (0.04), `energyresolution06` (0.06).
    **Code:** Actual variants are `energyscale26up` (1.026), `energyscale26down` (0.974), `energyresolution5` (0.05), `energyresolution7` (0.07), `energyresolution8` (0.08). The values, names, and semantics are all different from what the wiki states.
    **Verdict:** WRONG -- variant names, values, and count are all incorrect

20. **Claim (concepts/systematic-variations.md):** SYST_GROUPS contains groups `purity: [iso, ntf, fit]` and `eff: [tight_bdt, fudge]`.
    **Code:** Actual SYST_GROUPS (lines 167-174):
    - `purity: [noniso, nt_bdt, purity_fit, mc_purity_correction]`
    - `eff: [tight_bdt, npb_cut]`
    Neither the member names nor the group compositions match.
    **Verdict:** WRONG -- group compositions are fabricated

21. **Claim (concepts/systematic-variations.md and guides/adding-a-systematic.md):** VARIANTS is a dict with keys like `"my_new_var_up": {"overrides": {...}}`.
    **Code:** VARIANTS is actually a LIST of dicts, each with `name` key and flat override keys (not nested under `"overrides"`). The wiki presents the entirely wrong data structure.
    **Verdict:** WRONG -- VARIANTS is a list (not dict) and the override structure differs

---

## 2. Cross-Reference Consistency

### 2a. ABCD article vs efficiency-yield article

The ABCD article (concepts/abcd-method.md) and the efficiency-yield article (pipeline/04-efficiency-yield.md) are consistent on:
- Region naming (tight/nontight, iso/noniso)
- Histogram names (`h_tight_iso_cluster_0`, etc.)
- Signal leakage formula and code location
- Purity fitting options (fit_option 0=erf, 1=pade)
- Cross-reference link from 04-efficiency-yield.md to abcd-method.md resolves correctly

**Status:** Consistent

### 2b. Config-schema vs pipeline articles

The config-schema article (reference/config-schema.md) is generally consistent with pipeline articles on field names and types. However:

- **Disagreement on non_tight BDT fields:** The config-schema article correctly states the code reads `bdt_max_slope` and `bdt_max_intercept` under non_tight. However, the ACTUAL config file (config_bdt_nom.yaml) has `bdt_min_slope` and `bdt_min_intercept` under non_tight (lines 205-206). This means the code falls back to defaults (0 slope, bdt_max as intercept). The config-schema article should flag this config file bug/discrepancy more prominently.

- **Nominal fit_option:** config-schema says `fit_option: 1` is nominal with code default 0 -- this is correct and consistent with the actual config.

**Status:** Mostly consistent, but the non_tight BDT field naming mismatch between config file and code is under-documented

### 2c. Constants-sync vs mc-samples

The cross-section values in constants-sync.md match mc-samples.md exactly. Both reference CrossSectionWeights.h as the source of truth. The jet pT window boundaries in mc-samples.md match CrossSectionWeights.h.

**Status:** Consistent

### 2d. Systematic-variations vs make_bdt_variations.py

This is the worst cross-reference problem. The systematic-variations article contains fabricated variant names and group compositions that do not match the actual code at all. See findings 18-21 above.

**Status:** Severely inconsistent

### 2e. Markdown cross-reference links

All cross-reference links resolve to existing files:
- `concepts/abcd-method.md` -> exists
- `concepts/systematic-variations.md` -> exists
- `reference/constants-sync.md` -> exists
- `pipeline/03-bdt-application.md` -> exists
- `guides/adding-a-systematic.md` -> exists

**Status:** All links valid

---

## 3. Usability Review

### 3a. File Paths

All file paths mentioned are absolute GPFS paths or relative within the repo. The paths are accurate and a fresh agent could locate files. Specific good examples:
- Data flow article lists concrete file patterns
- Running the full pipeline guide uses absolute paths throughout

**Minor issue:** The wiki never mentions the `/gpfs/mnt/gpfs02` prefix that maps to `/sphenix/user`. A new user might not know these are equivalent GPFS mount points.

### 3b. Branch Names

Branch names are documented accurately with correct case sensitivity:
- `cluster_bdt_{node}_{model_name}` -- correct
- `cluster_npb_score_{node}` -- correct
- All shower shape branches documented correctly
- The `cluster_e_array_{node}[n*49]` notation is correct

### 3c. Entry Points

Good entry points provided:
- Each pipeline stage article starts with a key files table
- Function signatures are provided
- The running-full-pipeline guide gives copy-paste commands

**Gaps:**
- The tree-making article lacks a quick-check for verifying existing slimtrees are healthy
- No article explains how to verify BDT model files exist before running apply_BDT.C

### 3d. Self-Containedness

Most articles are self-contained enough for a focused task. The ABCD article is particularly strong -- a fresh session could understand the full method from that article alone.

Weaker on self-containedness:
- The systematic-variations article would mislead a fresh session due to the fabricated variant names
- The adding-a-systematic guide shows the wrong data structure for VARIANTS

---

## 4. Missing Content

### 4a. Important Files Not in File Inventory

The following files exist in the codebase but are not in the file inventory:

**efficiencytool/ (12 missing files):**
- `time_energy_corr.C` -- mentioned in constants-sync (16.67 ns bug) but not in inventory
- `calc_pileup_range.C` -- pileup rate calculation
- `investigate_truth_photons.C` -- truth photon debugging
- `investigate_vertex_shift.C` -- vertex shift studies
- `investigate_deltaR.C` -- deltaR matching studies
- `BDTScoreVsET.C` -- BDT score visualization
- `showershape_vertex_check.C` -- vertex-dependent shower shape
- `CheckPhotonPtRange.C` -- truth photon pT range validation
- `PlotTruthJetWindows.C` -- truth jet window visualization
- `TruthJetInETWindow.C` -- jet truth in ET windows
- `MakeSatPlots.C` -- saturation plotting
- `MakeEffPlots.C` -- efficiency plotting

**plotting/ (15+ missing files):**
- `plot_reweight.C` -- unfolding reweighting (referenced in concepts/unfolding.md but not in inventory)
- `plot_mc_purity_correction.C` -- MC purity correction visualization
- `CONF_plots.C` -- conference/publication plots
- `plot_showershapes_selections.C` -- shower shape with selections
- `plot_npb_time_bkgsub.C` -- NPB timing background subtraction
- `plot_background_recoisoET_overlay.C` -- background isolation overlay
- `plot_cluster_rbr_QA.C` -- run-by-run cluster QA
- `plot_mbd_sigma_efficiency.C` -- MBD sigma efficiency
- `plot_bdt_variations.C` -- BDT variation comparison
- `plot_vertex_check.C` -- vertex distributions
- `Closure.C` (in plotting/) -- distinct from the one listed under efficiencytool
- `plot_SB.C` -- signal/background plots
- `plot_cluster_timing.C`, `plot_cluster_mbd_time.C`, `plot_tower_timing.C` -- timing analysis

### 4b. Config Fields Not in Schema

The config-schema article is missing several fields that exist in the actual configs and code:
- `analysis.mbd_avg_sigma_max` / `mbd_avg_sigma_min` -- MBD pileup metric cut
- `analysis.cluster_mbd_time_min` / `cluster_mbd_time_max` -- cluster timing window
- `analysis.common_b2bjet_cut` -- back-to-back jet cut
- `analysis.common.wr_cogx_bound` -- radial width bound (only `cluster_weta_cogx_bound` is listed)
- `analysis.cluster_escale` and `analysis.cluster_eres` -- mentioned in passing but not formally documented in the tight/non_tight/common field tables
- `input.use_split_bdt_npb` -- NPB model variant selection (apply_BDT.C reads this)

### 4c. Pipeline Steps Not in Overview

- The **double-interaction** pipeline (run_showershape_double.sh, DoubleInteractionCheck.C) is not mentioned in the pipeline overview or flow diagram, despite being a significant sub-pipeline with its own two-pass architecture
- The **shower shape analysis** pipeline (ShowerShapeCheck.C, run_showershape_double.sh, plot_showershapes_variations.C) is not represented as a pipeline stage
- The **NPB purity study** (NPB_PurityStudy.C) is not documented in any pipeline or concepts article

---

## 5. Suggested Fixes

### Priority 1: CRITICAL -- Systematic Variations Article Contains Fabricated Data

**Article:** `concepts/systematic-variations.md`
**Sections:** "Key Variation Types", "Systematic Groups (SYST_GROUPS)"
**Fix:** Replace the entire "Key Variation Types" section with the actual variant definitions from `make_bdt_variations.py`. Specifically:

1. Replace the variant name tables with actual names from the code:
   - `noniso04`/`noniso10` (not `iso_up`/`iso_down`)
   - `ntbdtmin05`/`ntbdtmin10` (not `ntf_up`/`ntf_down`)
   - `tightbdt50`/`tightbdt70` (not `tightbdt_up`/`tightbdt_down`)
   - `energyscale26up`/`energyscale26down` (not `energyscale24`/`energyscale26`)
   - `energyresolution5` (not `energyresolution02/04/06`)
   - `purity_pade` (not `fit_up`/`fit_down`)
   - No `fudge_up`/`fudge_down` or `mbd_eff_up`/`mbd_eff_down` or `unfold_nor` variants exist

2. Replace the SYST_GROUPS table with actual values:
   ```
   purity: [noniso, nt_bdt, purity_fit, mc_purity_correction]
   eff:    [tight_bdt, npb_cut]
   escale: [escale]
   eres:   [eres]
   mbd:    [mbd]        # placeholder, no variants yet
   unfolding: [reweight]
   ```

3. Replace the syst_type table with actual types from code (line 148-164 of make_bdt_variations.py).

### Priority 2: CRITICAL -- Adding-a-Systematic Guide Shows Wrong Data Structure

**Article:** `guides/adding-a-systematic.md`
**Section:** "Step 1: Define the Variation", subsection 1a
**Fix:** Replace the dict-based VARIANTS example with the actual list-of-dicts format:

```python
# WRONG (current wiki):
VARIANTS = {
    "my_new_var_up": {
        "overrides": {"my_config_key": new_value_up},
        "syst_type": "my_type",
        "syst_role": "up",
    },
}

# CORRECT (actual code):
VARIANTS = [
    dict(name="my_new_var_up", my_config_key=new_value_up,
         syst_type="my_type", syst_role="up"),
]
```

Also update section 1b: the OVERRIDE_MAP uses tuple values `(["analysis", "section"], "leaf_key")`, not simple string paths like `"analysis.my_config_key"`.

### Priority 3: HIGH -- Config-Schema Non-Tight BDT Field Discrepancy

**Article:** `reference/config-schema.md`
**Section:** "analysis.tight / analysis.non_tight"
**Fix:** Add a warning note that the ACTUAL config file (config_bdt_nom.yaml) has `bdt_min_slope` and `bdt_min_intercept` under `non_tight` instead of the expected `bdt_max_slope` and `bdt_max_intercept` that the code reads. This means the parametric non_tight upper BDT boundary falls back to code defaults (slope=0, intercept=bdt_max). This is likely a config bug, not a code bug.

### Priority 4: HIGH -- Division Protection Claim is Oversimplified

**Article:** `pipeline/02-bdt-training.md` (Gotchas) and `concepts/shower-shape-variables.md`
**Fix:** The claim "apply_BDT.C uses `numerator > 0` check" is only true for the photon-ID BDT block (line 271: `cluster_e11 > 0`). The NPB score block (lines 436-448) checks the DENOMINATOR (`cluster_e33 > 0`, `cluster_e35 > 0`, etc.). Change the wiki text to:

> apply_BDT.C has two division protection patterns: the photon-ID block checks the numerator (`cluster_e11 > 0`), while the NPB score block checks the denominator (`cluster_e33 > 0`). BDTinput.C uses `safe_div` (checks denominator != 0 and isfinite).

### Priority 5: MEDIUM -- Missing Config Fields in Schema

**Article:** `reference/config-schema.md`
**Fix:** Add the following fields to the schema:
- `analysis.mbd_avg_sigma_max` -- MBD average sigma timing cut (used by pileup rejection variants)
- `analysis.cluster_mbd_time_min` / `cluster_mbd_time_max` -- cluster timing window (used by timing systematics)
- `analysis.common_b2bjet_cut` -- back-to-back jet cut toggle
- `analysis.common.wr_cogx_bound` -- radial width preselection bound
- `input.use_split_bdt_npb` -- NPB TMVA model variant flag

### Priority 6: MEDIUM -- File Inventory Incomplete

**Article:** `reference/file-inventory.md`
**Fix:** Add the following files to the inventory:

Under `efficiencytool/` study macros:
- `time_energy_corr.C` -- Tower timing vs energy correction (contains 16.67 ns bug)
- `calc_pileup_range.C` -- Pileup rate calculation for crossing angle periods
- `BDTScoreVsET.C` -- BDT score vs cluster ET visualization
- `showershape_vertex_check.C` -- Vertex-dependent shower shape study
- `investigate_truth_photons.C` -- Truth photon debugging
- `investigate_vertex_shift.C` -- Vertex shift studies
- `investigate_deltaR.C` -- DeltaR matching studies

Under `plotting/` analysis plots:
- `plot_reweight.C` -- Unfolding prior comparison (already referenced in concepts/unfolding.md)
- `plot_mc_purity_correction.C` -- MC purity correction visualization
- `CONF_plots.C` -- Conference/publication-ready plots
- `plot_bdt_variations.C` -- BDT variation comparison
- `plot_mbd_sigma_efficiency.C` -- MBD sigma efficiency

### Priority 7: MEDIUM -- Pipeline Overview Missing Sub-Pipelines

**Article:** `pipeline/overview.md`
**Fix:** Add a section "Auxiliary Pipelines" covering:
1. **Double-interaction blending** (`run_showershape_double.sh`) -- two-pass pipeline for pileup studies
2. **Shower shape analysis** (`ShowerShapeCheck.C`) -- distributions in ABCD regions
3. **Toy double-interaction simulation** (`DoubleInteractionCheck.C` + `run_double_interaction.sh`)

These are significant analysis workflows (~3000+ lines of code combined) not mentioned in the pipeline overview.

### Priority 8: LOW -- plotcommon.h Line Count

**Article:** `pipeline/05-plotting-systematics.md`
**Claim:** plotcommon.h is 93 lines.
**Code:** `wc -l` = 92 lines.
**Fix:** Change "93 lines" to "92 lines".

### Priority 9: LOW -- DoubleInteractionCheck Line Count in File Inventory

**Article:** `reference/file-inventory.md`
**Claim:** DoubleInteractionCheck.C is 1446 lines.
**Code:** `wc -l` = 1446 lines. CONFIRMED.

### Priority 10: LOW -- BlairUtils.C Line Count

**Article:** `pipeline/05-plotting-systematics.md`
**Claim:** BlairUtils.C is 544 lines.
**Code:** `wc -l` = 543 lines.
**Fix:** Change "544 lines" to "543 lines".

---

## Summary

| Category | Count | Critical |
|----------|-------|----------|
| Claims checked | 21 | -- |
| CONFIRMED | 17 | -- |
| WRONG | 4 | Yes (all in systematic-variations and adding-a-systematic) |
| Cross-ref issues | 2 | 1 critical (systematic-variations vs code) |
| Missing files in inventory | 27+ | No |
| Missing config fields | 5 | No |
| Missing pipeline steps | 3 | No |
| Total suggested fixes | 10 | 2 critical, 1 high |

**Overall Assessment:** The wiki is highly accurate on factual details about source code structure, branch names, function signatures, constants, and line counts. The pipeline articles (01 through 05) and the concepts articles on ABCD, isolation, shower shapes, and unfolding are excellent -- well-sourced and verified against code.

The critical weakness is the systematic-variations article and the adding-a-systematic guide, which contain fabricated variant names, wrong data structures, and incorrect group compositions that do not match the actual `make_bdt_variations.py` code at all. A fresh session relying on these articles to add a systematic variation would produce non-functional code. These two articles require immediate rewriting from the actual source code.
