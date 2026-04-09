# Impact Analysis of Reported Issues

Each issue is assessed against the actual codebase as of 2026-04-08.

---

## Issue 1: pT Bin Mismatch (plotcommon.h vs config_bdt_nom.yaml)

**Severity: HIGH**
**Verdict: BUG -- but partially mitigated by how the plotting code actually works**

**Evidence:**
- `plotting/plotcommon.h` line 6: `ptRanges = {8,10,12,14,16,18,20,22,24,26,28,30,35}` (NptBins=12)
- `efficiencytool/config_bdt_nom.yaml` line 51: `pT_bins: [8,10,12,14,16,18,20,22,24,26,28,32,36]` (also 12 bins)
- The last two bin edges differ: plotcommon.h has (30, 35) while config_bdt_nom has (32, 36).

**How the plotting macros actually use bins:**
- `plot_final_selection.C` (the cross-section plot): reads histograms directly from `Photon_final_{var_type}.root` via `Get("h_unfold_sub_result")`. These histograms were produced by `CalculatePhotonYield.C`, which reads pT bins from the YAML config. The histogram binning is baked into the ROOT file. `plot_final_selection.C` does NOT re-bin using `ptRanges` from plotcommon.h. It only uses plotcommon.h for axis style, frames, and legend strings (via `init_plot()`). Therefore, **the actual cross-section data points are NOT affected**.
- `plot_truth_iso.C` DOES use `ptRanges[ipt]` in legend labels (lines 104, 144, 186), so the pT range labels on those plots would say "28-30 GeV" and "30-35 GeV" when the actual data covers "28-32 GeV" and "32-36 GeV". This is misleading but does not alter the data.
- Several other plotting macros (`plot_background_recoisoET_overlay.C`, `plot_double_interaction.C`, `plot_showershapes_variations.C`) read pT bins from the YAML config directly -- they are NOT affected.

**Impact:** The published cross-section values are computed correctly (they use YAML-driven bins). The mismatch causes wrong pT-range labels in some per-bin plots (like truth-iso plots) and potentially wrong bin-width normalization if any macro divides by `ptRanges[i+1] - ptRanges[i]`. However, the main cross-section macro (`plot_final_selection.C`) uses the histogram's native binning. Still HIGH because the mismatch creates a latent trap: any new plot using `ptRanges` will silently mislabel or mis-normalize the last two bins.

---

## Issue 2: Multiple Luminosity Values

**Severity: MEDIUM**
**Verdict: INTENTIONAL (multiple run periods) + cosmetic inconsistency in labels**

**Evidence -- five distinct values found:**
1. **16.2735 pb^-1** -- `config_bdt_nom.yaml` line 239 and all `config_bdt_*.yaml` systematic variations. This is the value that enters `CalculatePhotonYield.C` for the actual cross-section calculation.
2. **16.6 pb^-1** -- `plotcommon.h` lines 21, 25 (legend strings `strleg2_1`, `strleg5`); `ppg12_conf_note/main.tex` lines 36, 60, 260. This is a rounded display value used in plot labels and the conference note.
3. **16.8588 pb^-1** -- `config_showershape_1p5rad.yaml` line 246. Used in the showershape study for the 1.5 mrad period (possibly a different luminosity calculation method or run selection).
4. **16.9 pb^-1** -- `PPG12-analysis-note/selection.tex` line 135, `introduction.tex` line 24, `analysis.tex` line 264. Analysis note quotes this with asymmetric uncertainties (+1.4/-1.2).
5. **49.562 pb^-1** -- used for all-runs combined analysis (0 mrad + 1.5 mrad).

**Which is "correct"?**
- The actual cross-section uses 16.2735 from `config_bdt_nom.yaml`. This is the precise value for the 1.5 mrad period (runs 51274-54000).
- 16.6 in plotcommon.h and conf note is a rounded label. The label says "16.6 pb^-1" while the calculation uses 16.2735.
- 16.9 in the analysis note appears to be a different calculation methodology or a slightly different run selection.
- 16.8588 in the showershape config could reflect a different luminosity method (it is used only for the showershape study, not the main cross-section).

**Impact:** The cross-section calculation uses the config value (16.2735) consistently across all systematic variations. The scatter in label values (16.6, 16.9) is cosmetic but could confuse reviewers who compare the analysis note text against the actual computation. The showershape study uses 16.8588 which is a ~3.5% difference from the nominal -- this does not feed into the published cross-section but affects showershape study normalization.

---

## Issue 3: BDTinput.C Line 194 Copy-Paste Bug (et4_on reads et3_on)

**Severity: LOW**
**Verdict: BUG -- but zero physics impact given current configuration**

**Evidence:**
- `FunWithxgboost/BDTinput.C` line 194: `int et4_on = configYaml["analysis"]["et3_on"].as<int>(1);`
- This reads `et3_on` instead of `et4_on`. The default is 1 (enabled).
- In `config_bdt_nom.yaml` lines 70-71: `et3_on: 0` and `et4_on: 0`. So `et4_on` would get the value of `et3_on` = 0, which happens to equal the intended value of `et4_on` = 0.

**How et4_on is used:** In BDTinput.C, `et4_on` controls whether failing the et4 tight cut increments the non-tight fail counter (line 842: `nfail += et4_on`). In the nominal config, both `et3_on` and `et4_on` are 0, so the bug is invisible.

**When it would matter:** If someone set `et3_on: 1` and `et4_on: 0` (enable et3 cuts but disable et4 cuts), BDTinput.C would incorrectly enable et4 in the non-tight fail counting. However, BDTinput.C is used only for feature extraction to training text files, and the individual `_on` flags control non-tight classification (not BDT training inputs). The main analysis code `RecoEffCalculator_TTreeReader.C` line 280 reads `et4_on` correctly: `configYaml["analysis"]["et4_on"].as<int>(1)`.

**Impact:** Zero effect on current results. It is a latent bug that only activates if `et3_on != et4_on` in a BDTinput.C run, and even then only affects non-tight classification in the training-data selection (not the final analysis).

---

## Issue 4: npb_score_cut Dual YAML Paths

**Severity: MEDIUM**
**Verdict: BUG -- but currently masked by matching default values**

**Evidence:**
- `RecoEffCalculator_TTreeReader.C` line 494: reads `analysis.common.npb_score_cut` (default 0.5)
- `ShowerShapeCheck.C` line 78: reads `analysis.common.npb_score_cut` (default 0.5)
- `DoubleInteractionCheck.C` line 65: reads `analysis.npb_score_cut` (default 0.5)
- `NPB_PurityStudy.C` line 227: reads `analysis.npb_score_cut` (default 0.5)
- `embed_test/ShowerShapeCheck.C` line 57: reads `analysis.npb_score_cut` (default 0.5)

In `config_bdt_nom.yaml`, the field is defined at `analysis.common.npb_score_cut: 0.5` (line 213). There is NO `analysis.npb_score_cut` field at the top level.

**What happens:**
- RecoEffCalculator and ShowerShapeCheck find the value at `analysis.common.npb_score_cut` = 0.5. Correct.
- DoubleInteractionCheck looks for `analysis.npb_score_cut`, does not find it, falls back to default 0.5. Also 0.5, so it happens to match.
- If someone varied the NPB score cut by changing `analysis.common.npb_score_cut` to e.g. 0.3 (like in `config_bdt_npb03.yaml` line 213), DoubleInteractionCheck would NOT pick up the change -- it would still use its default 0.5.

**Impact:** For the nominal analysis (npb_score_cut=0.5 everywhere), there is no effect. For the `npb03` systematic variation, the main analysis correctly uses 0.3, but DoubleInteractionCheck (if run with that config) would silently use 0.5 instead. Since DoubleInteractionCheck is a study macro (not part of the cross-section pipeline), the published result is unaffected, but the double-interaction study results would be inconsistent with the claimed NPB cut variation.

---

## Issue 5: Non-tight BDT Config/Code Mismatch

**Severity: LOW**
**Verdict: BUG in study macros, correctly handled in main analysis**

**Evidence:**
The config `config_bdt_nom.yaml` defines parametric non-tight BDT upper bound (lines 203-206):
```yaml
non_tight:
  bdt_max: 0.5
  bdt_max_slope: -0.015
  bdt_max_intercept: 0.80
  bdt_min: 0.02
```

How each macro reads these:
- **RecoEffCalculator_TTreeReader.C** (main analysis, lines 477-480): Reads `non_tight_bdt_max`, `non_tight_bdt_min`, `non_tight_bdt_max_slope`, `non_tight_bdt_max_intercept` and applies parametric cut at line 2263: `bdt_score < non_tight_bdt_max_slope * cluster_Et[icluster] + non_tight_bdt_max_intercept`. **Correct.**
- **ShowerShapeCheck.C** (lines 488-489): Reads only `non_tight_bdt_max` and `non_tight_bdt_min` (no slope/intercept). Uses flat `bdt_score < non_tight_bdt_max` at line 1787. **Missing parametric cut.**
- **DoubleInteractionCheck.C** (lines 374-375): Same as ShowerShapeCheck -- flat `non_tight_bdt_max` only, at line 930. **Missing parametric cut.**

**Impact on main analysis:** None. RecoEffCalculator correctly applies the parametric non-tight BDT boundary. The published cross-section, efficiency, and ABCD background subtraction all use the parametric version.

**Impact on study macros:** ShowerShapeCheck and DoubleInteractionCheck use `non_tight_bdt_max = 0.5` as a flat cut. With the parametric formula, the actual upper bound is `0.80 - 0.015 * ET`, which gives 0.68 at ET=8 GeV and 0.35 at ET=30 GeV. The flat cut at 0.5 is more permissive at low ET and tighter at high ET. This means the non-tight sample in study macros has a somewhat different composition than in the main analysis.

---

## Issue 6: Missing MBDxsec Bib Entry

**Severity: MEDIUM**
**Verdict: BUG (missing reference)**

**Evidence:**
- `PPG12-analysis-note/systematics.tex` line 126: `\cite{MBDxsec}` references the MBD trigger cross-section of 25.2 mb.
- `PPG12-analysis-note/cite.bib`: grep for "MBDxsec" returns zero matches.

This means the LaTeX compilation produces a "[?]" in the text where the MBDxsec citation should appear. The BibTeX run would emit a warning about an undefined citation key.

**Impact:** This is not purely cosmetic. The MBD cross-section (25.2 mb) directly enters the luminosity calculation, which in turn scales the entire cross-section result. The luminosity systematic uncertainty (+9.1%/-6.7%) is one of the dominant systematics. Reviewers will specifically want to verify this reference. A missing citation makes it impossible for readers to check the source, which is a significant issue for an analysis note under review. However, it does not affect the actual computed numbers.

---

## Issue 7: Flat vs Parametric BDT in DoubleInteractionCheck/ShowerShapeCheck

**Severity: LOW**
**Verdict: UNCLEAR -- likely intentional simplification, but creates inconsistency**

**Evidence:**
- **RecoEffCalculator_TTreeReader.C** (main analysis, lines 428-430): Reads `tight_bdt_min_slope` and `tight_bdt_min_intercept`, applies parametric tight BDT cut at line 2228: `tight_bdt_min_et = tight_bdt_min_slope * cluster_Et[icluster] + tight_bdt_min_intercept`.
- **ShowerShapeCheck.C** (lines 439-441): Also reads `tight_bdt_min_slope` and `tight_bdt_min_intercept`, applies parametric tight BDT at line 1752: `tight_bdt_min_et = tight_bdt_min_slope * cluster_Et[icluster] + tight_bdt_min_intercept`. **Correct for tight.**
- **DoubleInteractionCheck.C** (line 355): Reads only flat `tight_bdt_min`. Uses flat threshold at line 917: `bdt_score > tight_bdt_min`. **Missing parametric cut.**

So the situation is:
- ShowerShapeCheck: tight BDT = parametric (correct), non-tight BDT = flat (missing parametric)
- DoubleInteractionCheck: tight BDT = flat (missing parametric), non-tight BDT = flat (missing parametric)

**The flat tight_bdt_min value in config_bdt_nom.yaml is 0.7.** The parametric formula gives `0.80 - 0.015 * ET`, which equals 0.68 at ET=8 and 0.35 at ET=30. So the flat 0.7 is tighter than parametric at all ET > 6.67 GeV. In DoubleInteractionCheck, the flat 0.7 cut rejects more clusters at high ET than the actual analysis does.

**Impact:** These are study macros, not the main cross-section pipeline. The purpose of DoubleInteractionCheck is to assess purity degradation from pileup -- qualitative conclusions about "purity drops under double interaction" would remain valid. However, the quantitative purity numbers would differ from the main analysis, making direct comparisons imprecise. ShowerShapeCheck has the tight BDT correct but the non-tight wrong. Whether this is intentional simplification or oversight is unclear without user input.

---

## Issue 8: N_PT_BINS=10 vs NptBins=12

**Severity: MEDIUM**
**Verdict: BUG**

**Evidence:**
- `plotting/make_selection_report.py` line 28: `N_PT_BINS = 10  # noqa: matches NptBins in plotcommon.h`
- `plotting/make_comparison_report.py` line 34: `N_PT_BINS = 10`
- `plotting/plotcommon.h` line 5: `NptBins = 12`
- `config_bdt_nom.yaml` has 12 pT bins.

The comment "matches NptBins in plotcommon.h" is itself wrong -- NptBins is 12, not 10.

**What the Python scripts do:** They iterate `for ipt in range(N_PT_BINS)` when collecting per-pT-bin isolation ET plots (make_selection_report.py line 145, make_comparison_report.py lines 238, 271). With N_PT_BINS=10, they only collect iso_ET plots for pT bins 0-9, missing bins 10 and 11 (the last two pT bins: 28-32 and 32-36 GeV).

**Are the last 2 bins populated?** The config has `pT_bins: [8,10,12,...,28,32,36]` and the analysis processes data up to 36 GeV. The trigger (photon-30) provides events in this range. These bins ARE populated in the data.

**Impact:** The automated LaTeX reports are missing real physics content. The isolation ET distributions in the 28-32 and 32-36 GeV bins are not shown in the selection and comparison reports. These are exactly the highest-pT bins where statistics are lowest and the analysis is most sensitive to systematic effects. Reviewers examining the reports would not see these bins, which could hide problems in the high-ET region. However, the underlying data analysis is unaffected -- only the report generation is incomplete.

---

## Summary Table

| Issue | Severity | Verdict | Affects Published Cross-Section? |
|-------|----------|---------|----------------------------------|
| 1. pT bin mismatch | HIGH | BUG | No (main plot uses ROOT file bins), but wrong labels on some per-bin plots |
| 2. Luminosity scatter | MEDIUM | INTENTIONAL + cosmetic | No (config value 16.2735 used consistently in calculation) |
| 3. BDTinput.C et4_on copy-paste | LOW | BUG | No (et3_on == et4_on == 0 in current config) |
| 4. npb_score_cut dual paths | MEDIUM | BUG | No (defaults match; only affects study macros with non-default NPB cuts) |
| 5. Non-tight BDT parametric | LOW | BUG (study macros only) | No (main analysis uses parametric correctly) |
| 6. Missing MBDxsec bib | MEDIUM | BUG | No (numbers correct; reference missing for reviewers) |
| 7. Flat tight BDT in DI study | LOW | UNCLEAR | No (study macro, not cross-section pipeline) |
| 8. N_PT_BINS=10 | MEDIUM | BUG | No (reports miss last 2 bins; data analysis unaffected) |

**Bottom line:** None of these issues affect the published cross-section values. The main analysis pipeline (RecoEffCalculator + CalculatePhotonYield + plot_final_selection) is internally consistent. The issues cluster in three categories: (a) study/diagnostic macros using simplified cuts, (b) cosmetic labeling inconsistencies, and (c) report generation missing high-pT bins. The most actionable items are #1 (synchronize plotcommon.h), #6 (add the bib entry), and #8 (fix N_PT_BINS to 12).
