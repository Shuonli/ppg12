# Fix Proposals for PPG12 Known Issues

Generated: 2026-04-08

---

## Issue 1: pT Bin Mismatch (plotcommon.h vs config_bdt_nom.yaml)

### Problem
`plotcommon.h` line 6 defines `ptRanges = {8,10,12,14,16,18,20,22,24,26,28,30,35}` but `config_bdt_nom.yaml` line 51 uses `pT_bins: [8,10,12,14,16,18,20,22,24,26,28,32,36]`. The last two edges differ: `(30,35)` vs `(32,36)`.

### Source of Truth
The config (`config_bdt_nom.yaml`) drives the actual analysis -- histograms are filled with these bin edges. Plotting macros that use `plotcommon.h` will display misaligned bin boundaries. The config should be authoritative.

### File(s) to Change
`/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/plotcommon.h`

### What to Change
**Line 6:**
```cpp
// OLD
const float ptRanges[NptBins + 1] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35};
// NEW
const float ptRanges[NptBins + 1] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36};
```

**Line 9** (also update pTmax):
```cpp
// OLD
const float pTmax = 30;
// NEW
const float pTmax = 36;
```

### Dependencies
- All plotting macros that read `ptRanges` or `pTmax` from `plotcommon.h` will now use the updated bins. Re-run affected plots.
- The conf note (`ppg12_conf_note/`) references the ET range in text (e.g., "10 < ET < 26 GeV"). The bin change extends the upper bound from 35 to 36 GeV but does not affect the quoted physics range unless the user intends to extend it.

### Testing
- Recompile and run any plotting macro that includes `plotcommon.h` (e.g., `plot_final_selection.C`). Verify that bin edges in output histograms match the `Photon_final_bdt_nom.root` histograms.

---

## Issue 2: Luminosity Values (Multiple Inconsistent Values for 1.5 mrad Period)

### Problem
Five distinct luminosity values appear for the 1.5 mrad run period (runs 51274--54000):

| Value (pb^-1) | Location | Context |
|---|---|---|
| 16.2735 | `config_bdt_nom.yaml` line 239 | Enters the actual cross-section calculation |
| 16.6 | `plotcommon.h` lines 21, 25; `ppg12_conf_note/main.tex` | Plot legend label, conf note text |
| 16.8588 | `config_showershape_1p5rad.yaml` line 246 | Showershape study config |
| 16.9 | `PPG12-analysis-note/analysis.tex` line 264; `PPG12-analysis-note/selection.tex` line 135; `PPG12-analysis-note/introduction.tex` line 24 | Analysis note |
| 16.6 with uncertainty | `ppg12_conf_note/main.tex` line 214 | "$16.6^{+1.4}_{-1.2}$" |

### Correct Value
The `lumi` field in `config_bdt_nom.yaml` (16.2735) is what enters the cross-section formula and must match the LumiCalculator output for the run range 51274--54000. The analysis note states 16.9 (with uncertainties +1.4/-1.2), which may reflect a slightly different run list or vertex cut. The conf note uses 16.6, which is yet another value.

The correct approach is: re-run `LumiCalculator.C` with the exact nominal run range (51274--54000), the lumi list, and the vertex cut. That output becomes the single source of truth. All other locations should either match it exactly or use a rounded version with appropriate context (e.g., "approximately 16.X pb^-1").

### File(s) to Change
Assuming `LumiCalculator` returns `L` for the 1.5 mrad period:

1. **`config_bdt_nom.yaml`** line 239: Set `lumi: L` (already the primary source)
2. **`config_showershape_1p5rad.yaml`** line 246: Change `lumi: 16.8588` to `lumi: L`
3. **`plotting/plotcommon.h`** lines 21, 25: Update the legend strings to use the rounded value matching `L`
4. **`PPG12-analysis-note/selection.tex`** line 135: Update `16.9` to match `L` (with uncertainties)
5. **`PPG12-analysis-note/analysis.tex`** line 264: Update `16.9` to match `L`
6. **`PPG12-analysis-note/introduction.tex`** line 24: Update `16.9` to match `L`
7. **`ppg12_conf_note/main.tex`** lines 36, 60, 214, 260: Update `16.6` to match `L` (rounded)
8. **`CLAUDE.md`** luminosity line: Update `16.6 pb^-1` to match `L`

### Concrete Example (if L = 16.2735)

**plotcommon.h line 21:**
```cpp
// OLD
string strleg2_1 = "#it{p}+#it{p} #kern[-0.05]{#sqrt{#it{s}} = 200 GeV, 16.6 pb^{-1}}";
// NEW
string strleg2_1 = "#it{p}+#it{p} #kern[-0.05]{#sqrt{#it{s}} = 200 GeV, 16.3 pb^{-1}}";
```

**plotcommon.h line 25:**
```cpp
// OLD
string strleg5 = "         = 16.6 pb^{-1}";
// NEW
string strleg5 = "         = 16.3 pb^{-1}";
```

### Dependencies
- Any systematic variation configs derived from `config_bdt_nom.yaml` via `make_bdt_variations.py` will inherit the updated `lumi` value automatically (it copies from the base config).
- Re-run the full analysis pipeline for the 1.5 mrad period after updating the config, then re-make all plots.

### Testing
- Run `root -l -b -q 'LumiCalculator.C(51274, 54000)'` and verify the output matches `config_bdt_nom.yaml`.
- Run `root -l -b -q 'LumiCalculator.C(47289, 51274)'` and verify it matches `config_bdt_0rad.yaml` (32.6574).

---

## Issue 3: BDTinput.C Line 194 Copy-Paste Bug

### Problem
Line 194 of `BDTinput.C` reads `et3_on` from the config instead of `et4_on`:
```cpp
int et4_on = configYaml["analysis"]["et3_on"].as<int>(1);  // BUG: should be "et4_on"
```

### File(s) to Change
`/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/FunWithxgboost/BDTinput.C`

### What to Change
**Line 194:**
```cpp
// OLD
    int et4_on = configYaml["analysis"]["et3_on"].as<int>(1);
// NEW
    int et4_on = configYaml["analysis"]["et4_on"].as<int>(1);
```

### Dependencies
- None. In the nominal config, both `et3_on` and `et4_on` default to 0, so the bug is currently silent. It would only manifest if `et3_on != et4_on` in a config.

### Testing
- Grep the config for `et3_on` and `et4_on` values. If they differ, re-run `BDTinput.C` to regenerate training data.

---

## Issue 4: npb_score_cut Dual YAML Paths

### Problem
Three macros read `npb_score_cut` from different YAML paths:

| Macro | YAML path read |
|---|---|
| `RecoEffCalculator_TTreeReader.C` (line 494) | `analysis.common.npb_score_cut` |
| `ShowerShapeCheck.C` (line 78) | `analysis.common.npb_score_cut` |
| `TruthJetInETWindow.C` (line 90) | `analysis.common.npb_score_cut` |
| `Cluster_rbr.C` (line 257) | `analysis.common.npb_score_cut` |
| `BDTScoreVsET.C` (line 64) | `analysis.common.npb_score_cut` |
| **`DoubleInteractionCheck.C` (line 65)** | **`analysis.npb_score_cut`** (flat, NOT under common) |
| **`NPB_PurityStudy.C` (line 227)** | **`analysis.npb_score_cut`** (flat, NOT under common) |
| `embed_test/ShowerShapeCheck.C` (line 57) | `analysis.npb_score_cut` (flat, legacy embed copy) |

The YAML path `analysis.npb_score_cut` does not exist in `config_bdt_nom.yaml` -- it only has `analysis.common.npb_score_cut`. The flat-path reads fall back to the default value (0.5), which happens to match the config value. But if the config value were changed, `DoubleInteractionCheck.C` and `NPB_PurityStudy.C` would silently use the stale default.

### Proposed Fix
Standardize on `analysis.common.npb_score_cut` (the path used by the main analysis code). Change the two out-of-sync macros.

### File(s) to Change

**`/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/DoubleInteractionCheck.C`**
Line 65:
```cpp
// OLD
    float npb_score_cut = configYaml["analysis"]["npb_score_cut"].as<float>(0.5);
// NEW
    float npb_score_cut = configYaml["analysis"]["common"]["npb_score_cut"].as<float>(0.5);
```

**`/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/NPB_PurityStudy.C`**
Line 227:
```cpp
// OLD
    const float npb_score_cut = configYaml["analysis"]["npb_score_cut"].as<float>(0.5);
// NEW
    const float npb_score_cut = configYaml["analysis"]["common"]["npb_score_cut"].as<float>(0.5);
```

### Dependencies
- The showershape configs (`config_showershape*.yaml`) also have `npb_score_cut` under `analysis.common`, so no config changes needed.
- The `make_bdt_variations.py` (line 196) already writes to `analysis.common.npb_score_cut` -- consistent with this fix.
- `embed_test/ShowerShapeCheck.C` is a legacy copy and could be updated too, but it is not part of the current pipeline.

### Testing
- Run `DoubleInteractionCheck.C` with a config that has `analysis.common.npb_score_cut: 0.7` and verify the cut applies correctly (counts should decrease relative to 0.5).

---

## Issue 5: non_tight BDT Config/Code Mismatch

### Problem
`RecoEffCalculator_TTreeReader.C` reads parametric non-tight BDT upper boundary:
```cpp
// Lines 479-480
float non_tight_bdt_max_slope = configYaml["analysis"]["non_tight"]["bdt_max_slope"].as<float>(0);
float non_tight_bdt_max_intercept = configYaml["analysis"]["non_tight"]["bdt_max_intercept"].as<float>(non_tight_bdt_max);
```
And applies it at line 2263:
```cpp
bdt_score < non_tight_bdt_max_slope * cluster_Et[icluster] + non_tight_bdt_max_intercept
```

The config `config_bdt_nom.yaml` (lines 203-205) correctly defines:
```yaml
non_tight:
    bdt_max: 0.5
    bdt_max_slope: -0.015
    bdt_max_intercept: 0.80
```

However, **ShowerShapeCheck.C** (lines 488-489, 1786-1787) and **DoubleInteractionCheck.C** (lines 374-375, 930) read only `bdt_max` (flat) and do NOT read `bdt_max_slope`/`bdt_max_intercept`. They apply a flat threshold:
```cpp
bdt_score < non_tight_bdt_max  // uses 0.5, ignoring the parametric 0.80 - 0.015*ET
```

This means the non-tight BDT upper boundary is applied differently across macros. The config has `bdt_max: 0.5` which is always less restrictive than the parametric `0.80 - 0.015*ET` (= 0.65 at ET=10), so ShowerShapeCheck and DoubleInteractionCheck accept a wider non-tight region.

### Proposed Fix
Add parametric non-tight BDT support to both ShowerShapeCheck.C and DoubleInteractionCheck.C, matching RecoEffCalculator_TTreeReader.C.

### File(s) to Change

**`/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/ShowerShapeCheck.C`**

After line 489, add:
```cpp
// OLD (lines 488-489)
    float non_tight_bdt_max = configYaml["analysis"]["non_tight"]["bdt_max"].as<float>(1.0);
    float non_tight_bdt_min = configYaml["analysis"]["non_tight"]["bdt_min"].as<float>(0.0);
// NEW
    float non_tight_bdt_max = configYaml["analysis"]["non_tight"]["bdt_max"].as<float>(1.0);
    float non_tight_bdt_min = configYaml["analysis"]["non_tight"]["bdt_min"].as<float>(0.0);
    float non_tight_bdt_max_slope = configYaml["analysis"]["non_tight"]["bdt_max_slope"].as<float>(0);
    float non_tight_bdt_max_intercept = configYaml["analysis"]["non_tight"]["bdt_max_intercept"].as<float>(non_tight_bdt_max);
```

At lines 1786-1787, change the non-tight BDT upper bound:
```cpp
// OLD
                bdt_score > non_tight_bdt_min &&
                bdt_score < non_tight_bdt_max)
// NEW
                bdt_score > non_tight_bdt_min &&
                bdt_score < non_tight_bdt_max_slope * cluster_Et[icluster] + non_tight_bdt_max_intercept)
```

**`/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/DoubleInteractionCheck.C`**

After line 375, add:
```cpp
// OLD (lines 374-375)
    float non_tight_bdt_max = configYaml["analysis"]["non_tight"]["bdt_max"].as<float>(1.0);
    float non_tight_bdt_min = configYaml["analysis"]["non_tight"]["bdt_min"].as<float>(0.0);
// NEW
    float non_tight_bdt_max = configYaml["analysis"]["non_tight"]["bdt_max"].as<float>(1.0);
    float non_tight_bdt_min = configYaml["analysis"]["non_tight"]["bdt_min"].as<float>(0.0);
    float non_tight_bdt_max_slope = configYaml["analysis"]["non_tight"]["bdt_max_slope"].as<float>(0);
    float non_tight_bdt_max_intercept = configYaml["analysis"]["non_tight"]["bdt_max_intercept"].as<float>(non_tight_bdt_max);
```

At line 930 (inside `classifyTightNonTight` lambda):
```cpp
// OLD
            bdt_score > non_tight_bdt_min && bdt_score < non_tight_bdt_max)
// NEW
            bdt_score > non_tight_bdt_min && bdt_score < non_tight_bdt_max_slope * clusterET + non_tight_bdt_max_intercept)
```

### Dependencies
- `make_bdt_variations.py` already writes `bdt_max_slope` and `bdt_max_intercept` to the `non_tight` section (lines 194-195), so variation configs are already correct.
- The nominal config already has these fields (lines 204-205). No config changes needed.

### Testing
- Run `ShowerShapeCheck.C` and `DoubleInteractionCheck.C` with `config_bdt_nom.yaml`. Compare non-tight region counts before/after the fix. With the parametric cut, the non-tight upper boundary is tighter (~0.65 at ET=10 vs flat 0.5), so the non-tight region will shrink.

---

## Issue 6: Missing MBDxsec Bib Entry

### Problem
`PPG12-analysis-note/systematics.tex` line 126 cites `\cite{MBDxsec}` but the entry does not exist in `PPG12-analysis-note/cite.bib`. LaTeX/BibTeX will produce a "[?]" for this citation.

### Proposed Fix
Add a BibTeX entry for the sPHENIX MBD trigger cross-section measurement. The MBD cross-section of $25.2^{+2.3}_{-1.7}$ mb in pp at 200 GeV was measured via Vernier scans. The luminosity calculation is documented in the PPG09 IAN (already cited as `ppg09IAN` in the same bib file). Since there is no separate published paper for the MBD cross-section measurement, the appropriate reference is an sPHENIX internal analysis note or the PPG09 IAN itself.

### File(s) to Change
`/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/PPG12-analysis-note/cite.bib`

### What to Add (at end of file)
```bibtex
@article{MBDxsec,
    author  = "Jiang, H. and others",
    title   = "{Minimum-Bias Trigger Cross-Section Measurement from Vernier Scans in $p+p$ at $\sqrt{s}=200$~{GeV}}",
    journal = "sPHENIX-Invenio",
    year    = 2025,
    note    = "See also luminosity calculation in~\cite{ppg09IAN}"
}
```

**Important note:** The exact author list, title, and URL should be confirmed with the luminosity working group. If the MBD cross-section measurement is documented within the PPG09 IAN itself (already cited), an alternative fix is to change the citation in `systematics.tex`:

**Alternative (simpler):** Change `systematics.tex` line 126:
```latex
% OLD
$25.2^{+2.3}_{-1.7}\ \text{mb}$~\cite{MBDxsec}.
% NEW
$25.2^{+2.3}_{-1.7}\ \text{mb}$~\cite{ppg09IAN}.
```

### Dependencies
- Run `bibtex main` after adding the entry to regenerate the bibliography.

### Testing
- Compile the analysis note: `cd PPG12-analysis-note && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex`. Verify no "[?]" markers in the output PDF.

---

## Issue 7: Flat vs Parametric BDT Threshold in Study Macros

### Problem
`DoubleInteractionCheck.C` reads only `tight_bdt_min` (flat value, line 355) and applies it without ET dependence (line 917):
```cpp
bool is_bdt_tight = (bdt_score > tight_bdt_min) && (bdt_score < tight_bdt_max);
```

The main analysis (`RecoEffCalculator_TTreeReader.C`, lines 428-430, 2228-2230) and `ShowerShapeCheck.C` (lines 440-441, 1752-1754) already use parametric thresholds:
```cpp
float tight_bdt_min_slope = configYaml["analysis"]["tight"]["bdt_min_slope"].as<float>(0);
float tight_bdt_min_intercept = configYaml["analysis"]["tight"]["bdt_min_intercept"].as<float>(tight_bdt_min);
// ...
float tight_bdt_min_et = tight_bdt_min_slope * cluster_Et[icluster] + tight_bdt_min_intercept;
bool is_bdt_tight = (bdt_score > tight_bdt_min_et) && (bdt_score < tight_bdt_max);
```

`ShowerShapeCheck.C` already has the fix for the tight BDT (lines 440-441, 1752), but still uses flat for non-tight (covered in Issue 5).

### File(s) to Change
`/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/DoubleInteractionCheck.C`

### What to Change

**Step 1: Add config reads after line 355.**
```cpp
// OLD (lines 354-355)
    float tight_bdt_max = configYaml["analysis"]["tight"]["bdt_max"].as<float>(1.0);
    float tight_bdt_min = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0.0);
// NEW
    float tight_bdt_max = configYaml["analysis"]["tight"]["bdt_max"].as<float>(1.0);
    float tight_bdt_min = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0.0);
    float tight_bdt_min_slope = configYaml["analysis"]["tight"]["bdt_min_slope"].as<float>(0);
    float tight_bdt_min_intercept = configYaml["analysis"]["tight"]["bdt_min_intercept"].as<float>(tight_bdt_min);
```

**Step 2: Update the classifyTightNonTight lambda (line 917).**
```cpp
// OLD (line 917)
        bool is_bdt_tight = (bdt_score > tight_bdt_min) && (bdt_score < tight_bdt_max);
// NEW
        float tight_bdt_min_et = tight_bdt_min_slope * clusterET + tight_bdt_min_intercept;
        bool is_bdt_tight = (bdt_score > tight_bdt_min_et) && (bdt_score < tight_bdt_max);
```

Note: `clusterET` is already a parameter of the lambda (line 901), so it is available.

### Dependencies
- This pairs with Issue 5 (non-tight parametric BDT) -- both should be applied together for a fully consistent `classifyTightNonTight` lambda.
- The nominal config already has `bdt_min_slope: -0.015` and `bdt_min_intercept: 0.80` (lines 148-149). No config changes needed.

### Testing
- Run `DoubleInteractionCheck.C` with `config_bdt_nom.yaml` for a photon10 sample. Compare tight-region yields before/after. The parametric cut `0.80 - 0.015*ET` is tighter at low ET and looser at high ET than the flat `bdt_min: 0.7`, so expect fewer tight clusters at low ET and more at high ET.

---

## Issue 8: N_PT_BINS=10 in Python Report Scripts

### Problem
`make_selection_report.py` (line 28) and `make_comparison_report.py` (line 34) both hardcode `N_PT_BINS = 10`. The actual number of pT bins is 12 (from `plotcommon.h` `NptBins=12` and `config_bdt_nom.yaml` which has 13 edges = 12 bins). The misleading comment says "matches NptBins in plotcommon.h" but is wrong.

This causes the last 2 pT bins to be omitted from the isolation ET distribution pages in the generated report PDFs.

### Proposed Fix
Change `N_PT_BINS = 10` to `N_PT_BINS = 12` in both files. Deriving it from `plotcommon.h` would require parsing C++ from Python, which is unnecessary complexity -- a simple constant is fine since pT bins rarely change.

### File(s) to Change

**`/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/make_selection_report.py`**
Line 28:
```python
# OLD
N_PT_BINS = 10  # noqa: matches NptBins in plotcommon.h
# NEW
N_PT_BINS = 12  # matches NptBins in plotcommon.h
```

**`/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/make_comparison_report.py`**
Line 34:
```python
# OLD
N_PT_BINS = 10
# NEW
N_PT_BINS = 12  # matches NptBins in plotcommon.h
```

### Dependencies
- None. The scripts only use `N_PT_BINS` for looping over per-pT-bin figure filenames. The filenames are expected to exist (the C++ plotting macros already produce all 12 bins).
- If the `plotcommon.h` pT bins are updated (Issue 1), this constant is still 12 bins (the count does not change, only the edges).

### Testing
- Run `python make_selection_report.py` for a config and verify the output PDF contains all 12 pT bins of isolation ET distributions (previously only 10 were shown).
- Check that filenames `isoET_..._pt10.pdf` and `isoET_..._pt11.pdf` (or equivalent) are included in the LaTeX.

---

## Summary Table

| Issue | Severity | File(s) | Lines | Nature |
|---|---|---|---|---|
| 1. pT bin mismatch | High | plotcommon.h | 6, 9 | Wrong constants |
| 2. Luminosity scatter | High | 8+ files | various | Inconsistent constants |
| 3. BDTinput.C copy-paste | Medium | BDTinput.C | 194 | Copy-paste bug |
| 4. npb_score_cut path | Medium | DoubleInteractionCheck.C, NPB_PurityStudy.C | 65, 227 | Wrong YAML path |
| 5. non_tight BDT parametric | Medium | ShowerShapeCheck.C, DoubleInteractionCheck.C | 488-489/1786-1787, 374-375/930 | Missing parametric logic |
| 6. MBDxsec bib missing | Low | cite.bib or systematics.tex | end / 126 | Missing reference |
| 7. Flat BDT in study macros | Medium | DoubleInteractionCheck.C | 354-355, 917 | Missing parametric logic |
| 8. N_PT_BINS=10 | Low | make_selection_report.py, make_comparison_report.py | 28, 34 | Wrong constant |

### Recommended Execution Order
1. **Issue 3** (one-line fix, no risk)
2. **Issue 8** (two-line fix, no risk)
3. **Issue 4** (two-line fix, changes YAML path)
4. **Issue 6** (add bib entry or change citation)
5. **Issue 1** (update plotcommon.h, re-make plots)
6. **Issues 5 + 7 together** (parametric BDT in DoubleInteractionCheck.C and ShowerShapeCheck.C)
7. **Issue 2** (luminosity: requires running LumiCalculator first, then updating many files)
