# Phase 4: Independent Final Review of PPG12 Wiki

**Reviewer:** Claude Opus 4.6 (fresh session, no prior context)
**Date:** 2026-04-08
**Method:** Attempted to answer 5 test questions using only the wiki, then verified every claim against source code.

---

## Test Question Results

### Q1: "What branches does apply_BDT.C write?"

**Wiki answer (03-bdt-application.md, lines 59-68):** 11 branches named `cluster_bdt_{clusternodename}_{model_name}[ncluster/F]` plus one `cluster_npb_score_{clusternodename}[ncluster/F]`. Lists all 11 model names and their exact feature vectors.

**Code verification:** Confirmed at `apply_BDT.C` lines 234-248. Branch format matches exactly. Model name list at lines 73-85 matches the wiki's 11 names in order: base, base_vr, base_v0, base_v1, base_v2, base_v3, base_E, base_v0E, base_v1E, base_v2E, base_v3E.

**Verdict: CORRECT.** Excellent coverage -- gives branch naming pattern, all model names, feature vectors, and ET threshold (7 GeV).

---

### Q2: "What pT bins does the analysis use?"

**Wiki answer:** Multiple places give this information:
- `plotcommon.h` section (05-plotting-systematics.md, line 27): `{8,10,12,14,16,18,20,22,24,26,28,30,35}` (12 bins)
- `config-schema.md` line 105: config_bdt_nom uses `[8..28,32,36]`
- `constants-sync.md` lines 45-52: explicit comparison table
- `common-issues.md` lines 7-13: flagged as critical issue

**Code verification:**
- `plotcommon.h` line 6: `{8,10,12,14,16,18,20,22,24,26,28,30,35}` -- matches wiki
- `config_bdt_nom.yaml` line 51: `[8,10,12,14,16,18,20,22,24,26,28,32,36]` -- matches wiki

**Verdict: CORRECT and WELL-DOCUMENTED.** The mismatch is flagged in 3 separate places (constants-sync, common-issues, plotting-systematics). A new Claude session would immediately understand this is a known issue.

---

### Q3: "How do I add a new systematic variation?"

**Wiki answer (guides/adding-a-systematic.md):** Step-by-step guide with 5 sub-steps: (1a) add to VARIANTS list, (1b) add to OVERRIDE_MAP, (1c) add to SYST_TYPES, (1d) add to SYST_GROUPS, (1e) add to FINAL_SYSTS. Then generate, run, aggregate, verify.

**Code verification against `make_bdt_variations.py`:**
- VARIANTS list structure: wiki says "list of dicts" with flat override keys -- CORRECT (line 34 `VARIANTS = [...]`)
- Override keys go directly in dict, not nested -- CORRECT
- Uses ruamel.yaml -- CORRECT (line 28 `from ruamel.yaml import YAML`)
- OVERRIDE_MAP structure: wiki example matches code (lines 185-223)
- validate_variants() runs at import time -- CORRECT (line 266)
- Checklist at end is comprehensive

**Issue found:** The guide's Step 1d says to add the type to `SYST_GROUPS`. However, there is a subtle gap between `SYST_TYPES` and `SYST_GROUPS` in the actual code. See SYST_GROUPS discrepancy below.

**Verdict: MOSTLY CORRECT.** The step-by-step guide would work. The example code is syntactically correct. However, a user following only this guide might miss that `SYST_GROUPS` membership is what actually matters for aggregation, not the `group` field in `SYST_TYPES`.

---

### Q4: "What cross-section weight does photon10 use?"

**Wiki answer (constants-sync.md, line 13 and mc-samples.md, line 8):** photon10cross = 6,944.675 pb. Weight relative to photon20 = 53.24.

**Code verification (`CrossSectionWeights.h` line 14):** `constexpr float photon10cross = 6944.675f` -- matches exactly.

Weight calculation (line 60): `photon10cross / photon20cross = 6944.675 / 130.4461 = 53.24` -- matches.

**Verdict: CORRECT.** Value is easily findable in both the constants-sync and mc-samples articles.

---

### Q5: "What's the ABCD method?"

**Wiki answer (concepts/abcd-method.md):** 134 lines covering: region definitions, histogram naming, signal leakage fractions (cB, cC, cD), self-consistent equation with R factor, toy MC purity extraction, purity fitting, MC purity correction, parametric thresholds, config fields, code locations.

**Code verification:**
- `myfunc()` in CalculatePhotonYield.C lines 8-37: matches wiki's equation exactly
- R=1 confirmed at lines 433 and 486: `f->SetParameter(7, 1)`
- Leakage fraction definitions match
- Histogram naming `h_tight_iso_cluster_0` etc. confirmed in the code
- Purity fit: erf default, Pade for fit_option=1 -- confirmed lines 585-589
- 20,000 toy experiments -- confirmed in the code

**Verdict: CORRECT and THOROUGH.** The explanation is detailed enough to understand the code without reading it first. The self-consistent equation, leakage formulas, and histogram naming are all accurate.

---

## Structural Review

### Link Resolution

All 19 links in `index.md` resolve to existing files. All inter-article cross-references (e.g., `../concepts/abcd-method.md`, `../reference/constants-sync.md`) resolve correctly. No broken links found.

### Index Completeness

Index lists 19 articles (excluding itself). File system contains exactly 19 article files (excluding `_reports/`). Every article that exists is indexed, and every index entry exists.

### Duplicate Information

Several facts appear in multiple articles, which is expected for a wiki. However, some duplications are internally inconsistent:

**1. SYST_GROUPS "eff" membership (CRITICAL DISCREPANCY)**

- `systematic-variations.md` table at lines 65-66 says: `eff | tight_bdt, npb_cut` (correct)
- BUT the SYST_TYPES table at lines 47-50 lists `vtx_reweight`, `bdt_model`, `b2bjet`, `timing` all with group = `eff`

This creates a misleading impression that all 6 types contribute to the eff group. In reality, `SYST_GROUPS["eff"]` only contains `["tight_bdt", "npb_cut"]`. The other 4 types have `group="eff"` in their SYST_TYPES metadata, but that field is purely decorative -- `calc_syst_bdt.py` never reads it.

Furthermore, the actual variant entries for `b2bjet` and `timingcut_2/5` have `syst_type=None, syst_role=None` in the code, making them cross-checks, not systematics. The wiki's SYST_TYPES table correctly shows them as having a syst_type (because they ARE defined in SYST_TYPES), but the wiki's all-variants table correctly shows b2bjet and timing under "Cross-Check Variants" -- so the information is correct but split across tables in a confusing way.

**2. Luminosity values -- consistently documented as inconsistent**

Five values (16.2735, 16.6, 16.8588, 16.9, and others) appear across constants-sync.md, common-issues.md, config-schema.md, and plotting-systematics.md. All articles acknowledge the discrepancy, which is good.

**3. fit_option mapping**

- abcd-method.md line 100: `fit_option=0 (default): Error function` -- CORRECT
- common-issues.md line 54: "Corrected: fit_option=0 is the error function (erf)" -- CORRECT
- config-schema.md line 107: `fit_option` nominal = 1, default = 0, "0=erf (default), 1=Pade" -- CORRECT

All three agree. Good.

**4. purity_pade variant naming**

The variant named `purity_pade` sets `fit_option=0`, which is the **erf** function, not Pade. This is confusing naming in the codebase, not in the wiki. The wiki correctly documents `fit_option=0` but does not explicitly flag that the variant name is misleading. The nominal config uses `fit_option=1` (Pade), so `purity_pade` actually switches FROM Pade TO erf. A new user following the wiki might not notice this inversion.

---

## Spot-Check Results (5 specific facts)

### Spot-Check 1: Division protection in apply_BDT.C

**Wiki claim (03-bdt-application.md, lines 53-55):** "e11_over_e33_BDT = (cluster_e11 > 0) ? ... : 0"

**Code (`apply_BDT.C` line 271):** `e11_over_e33_BDT = (cluster_e11[icluster] > 0) ? cluster_e11[icluster] / cluster_e33[icluster] : 0;`

**Verdict: CORRECT.** The division protection checks the numerator (e11 > 0), not the denominator, exactly as documented.

### Spot-Check 2: N_PT_BINS hardcoded in Python reports

**Wiki claim (common-issues.md lines 16-19):** `make_selection_report.py` and `make_comparison_report.py` hardcode `N_PT_BINS = 10`.

**Code:**
- `plotting/make_selection_report.py` line 28: `N_PT_BINS = 10  # noqa: matches NptBins in plotcommon.h`
- `plotting/make_comparison_report.py` line 34: `N_PT_BINS = 10`

**Verdict: CORRECT.** Both files hardcode 10 (not 12). The comment in make_selection_report.py is itself wrong ("matches NptBins" -- NptBins is 12, not 10).

### Spot-Check 3: MergeSim jet sample count

**Wiki claim (04-efficiency-yield.md, line 91):** MergeSim uses exactly 6 jet samples: jet5, jet8, jet12, jet20, jet30, jet40.

**Code verification:**

```
grep -c "jet" efficiencytool/MergeSim.C
```

I verified `MergeSim.C` (94 lines) merges jet5+jet8+jet12+jet20+jet30+jet40 = 6 samples.

**Verdict: CORRECT.**

### Spot-Check 4: n_nt_fail code default vs config value

**Wiki claim (config-schema.md, line 112):** `n_nt_fail` nominal = 0, code default = 1.

**Code:**
- `RecoEffCalculator_TTreeReader.C` line 282: `.as<int>(1)` (default = 1)
- `config_bdt_nom.yaml` line 57: `n_nt_fail: 0`

**Verdict: CORRECT.** The code default (1) differs from the config value (0). This means if the field is omitted from a config, behavior changes silently.

### Spot-Check 5: solid_angle and deta in cross-section calculation

**Wiki claim (04-efficiency-yield.md, lines 109):** `solid_angle = 2 * M_PI * 0.7 * 2`

**Wiki claim (05-plotting-systematics.md, line 54):** "Divides by deta = 1.4"

**Code:**
- `CalculatePhotonYield.C` line 55: `float solid_angle = 2 * M_PI * 0.7 * 2;`
- `plot_final_selection.C` line 40: `float deta = 1.4;`

**Verdict: CORRECT.**

---

## Errors and Issues Found

### Severity: HIGH

**1. SYST_GROUPS "eff" membership is misleadingly documented**

The systematic-variations.md SYST_TYPES table (lines 47-50) lists `vtx_reweight`, `bdt_model`, `b2bjet`, and `timing` with group = "eff". This suggests all 6 types contribute to the "eff" group's quadrature sum. In reality, `SYST_GROUPS["eff"]` only contains `["tight_bdt", "npb_cut"]`.

The SYST_GROUPS table on line 65 is correct (`eff | tight_bdt, npb_cut`), but because the SYST_TYPES table appears first and is more detailed, readers will likely believe the first table.

**Recommendation:** Add a note to the SYST_TYPES table clarifying that the "group" column shows the SYST_TYPES metadata field, which is informational only. The actual group membership is determined by SYST_GROUPS. Alternatively, remove `vtx_reweight, bdt_model, b2bjet, timing` from the SYST_TYPES table entirely since they are not in any SYST_GROUPS.

**2. SYST_TYPES "group" field mismatch: "unfold" vs "unfolding"**

The `reweight` type in SYST_TYPES has `group: "unfold"`, but the actual SYST_GROUPS key is `"unfolding"`. Since calc_syst_bdt.py never reads the `group` field from SYST_TYPES, this is a code-level inconsistency that doesn't affect behavior. However, the wiki documents both the SYST_TYPES group column (showing "unfold" implicitly via the "unfolding" group link) and the SYST_GROUPS key ("unfolding") without noting this mismatch.

This is minor because it's a code issue, not a wiki issue, but the wiki should note the discrepancy.

### Severity: MEDIUM

**3. purity_pade variant name is misleading (not flagged in wiki)**

The variant `purity_pade` sets `fit_option=0`, which is the **erf** function. The nominal config uses `fit_option=1` (Pade). So `purity_pade` actually switches **away from** Pade. The wiki documents the fit_option value correctly but does not call out the naming inversion. A note like "Despite its name, purity_pade switches from the nominal Pade to erf" would prevent confusion.

**4. b2bjet and timing variants: wiki lists as systematics, code has syst_type=None**

The wiki's systematic-variations.md lists these under SYST_TYPES with actual syst_types (`b2bjet`, `timing`). The code does define them in SYST_TYPES. However, the VARIANTS entries for `b2bjet`, `timingcut_2`, and `timingcut_5` all have `syst_type=None, syst_role=None`. They are cross-checks that HAPPEN to have types defined in SYST_TYPES, but no variant actually activates those types. The wiki's all-variants table correctly places them under "Cross-Check Variants" but the SYST_TYPES table gives the impression they are active systematics.

### Severity: LOW

**5. Line count for DoubleInteractionCheck.C**

Overview.md says "1446 lines" and file-inventory.md says "1446 lines". Verified: the file has exactly 1446 lines. However, overview.md also says "1535 lines" via a reference to the CLAUDE.md documentation. This is just a CLAUDE.md discrepancy, not a wiki issue -- the wiki's own number (1446) is correct.

*(Upon re-reading: the 1535 number does not appear in the wiki. The wiki consistently says 1446. This is correct.)*

---

## Overall Assessment

### Accuracy: 9/10

Nearly every fact verified against source code was correct. The only substantive accuracy issue is the SYST_GROUPS membership confusion (high severity but localized to one table). Cross-section values, branch names, feature vectors, file paths, config defaults, code defaults, and formula descriptions are all accurate.

### Completeness: 9/10

The wiki covers the full pipeline comprehensively. Every major file has documentation. The ABCD method, unfolding, isolation, shower shapes, and systematic variations are explained with enough depth for code modification. Minor gaps:
- No documentation of the `common_b2bjet_pt_min` default (7.0 appears in config-schema but not in the b2bjet systematic variant description)
- The `fittingerror` config parameter is mentioned in abcd-method.md (line 89) but not documented in config-schema.md
- NPB model features are documented in 03-bdt-application.md but not cross-referenced from 02-bdt-training.md's NPB section

### Navigability: 9.5/10

Excellent organization. The three-tier structure (Pipeline / Concepts / Reference + Guides) makes it easy to find information. Cross-references between articles work well. The index.md provides a clear table of contents. The constants-sync.md page is particularly valuable as a single place to check for known discrepancies.

### Actionability: 9/10

The wiki provides specific file paths, line numbers (approximate but useful), exact config field names, function signatures, and code patterns. The adding-a-systematic guide is genuinely step-by-step with a working example. The running-full-pipeline guide has copy-paste commands. The common-issues page flags real pitfalls.

The one actionability gap is that the SYST_GROUPS discrepancy could lead someone to believe adding a type to SYST_TYPES with `group="eff"` is sufficient to include it in the eff group's quadrature sum, when they also need to add it to the SYST_GROUPS list.

---

## Summary

This is a high-quality wiki. Across 20 articles and ~2,800 lines of documentation, I found:
- **1 high-severity issue** (SYST_GROUPS membership confusion in systematic-variations.md)
- **2 medium-severity issues** (purity_pade naming not flagged, b2bjet/timing appear as active systematics)
- **1 low-severity issue** (SYST_TYPES "group" field mismatch unfold/unfolding)
- **0 broken links**
- **0 missing articles**
- **0 factual errors** in cross-section values, branch names, feature vectors, config schema, or code patterns

A new Claude session arriving at this codebase would be able to navigate and modify the code correctly by following this wiki. The most dangerous trap would be the SYST_GROUPS membership issue, which could lead to a systematic type being defined but silently excluded from the final uncertainty.
