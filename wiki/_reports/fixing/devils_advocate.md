# Devil's Advocate Review of Reported Issues

Date: 2026-04-08

This review challenges the assumptions in `common-issues.md` and `constants-sync.md`
that the listed discrepancies are bugs requiring fixes.  For each issue the
counter-argument, the risk of a premature fix, and a verdict are provided.

---

## 1. pT Bin Mismatch (plotcommon.h vs config_bdt_nom.yaml)

**Claimed bug:** `plotcommon.h` has `{..., 28, 30, 35}` while `config_bdt_nom.yaml`
has `[..., 28, 32, 36]`.  Last two edges differ.

### Counter-argument

The pT bins that matter for the cross-section result are the ones in
`config_bdt_nom.yaml`.  `RecoEffCalculator_TTreeReader.C` reads bins from the
YAML config (line 324: `pT_bins = configYaml["analysis"]["pT_bins"]`), and
`CalculatePhotonYield.C` does not reference `plotcommon.h` at all -- it reads
histograms from ROOT files whose binning was set by the config.  The data flow
is: config -> RecoEff -> ROOT histograms -> CalculatePhotonYield.  At no point
does the cross-section chain touch `plotcommon.h`.

The plotting macros that include `plotcommon.h` (e.g. `plot_final_selection.C`,
`plot_isoET.C`) use `ptRanges` only for frame axis ranges and legend formatting,
not for rebinning the actual result histograms.  `plot_final_selection.C` reads
the data from `Photon_final_bdt_nom.root`, which already has the config binning
baked in.  The frame range is set to `[10, 32]` via `lowerx`/`upperx` hardcoded
on line 43-44, so the last two bins (28-32, 32-36) are only partially visible
regardless.

Where `ptRanges` is actually used for loop bounds (e.g. `plot_isoET.C` iterates
`for ipt = 0; ipt < NptBins`), it iterates over histogram bin indices.  The
histograms themselves have the correct (config) binning.  The iteration count
NptBins=12 matches the config's 12 bins.  The only mismatch would be in axis
labeling if a plot tried to print bin edges from the `ptRanges` array, but the
plotting code does not do this -- it uses `GetBinCenter()/GetBinLowEdge()` from
the histograms.

**Where it could bite:** If someone adds a new plotting macro that creates fresh
histograms using the `ptRanges` array as bin edges and then tries to divide them
by the result histograms, the bin edges would not match and ROOT would complain
or silently give wrong ratios.  This is a future-proofing concern, not a current
bug.

### Risk of fixing

Changing `plotcommon.h` to `{..., 28, 32, 36}` would affect all 54+ plotting
macros that include it.  Some of those macros may use `ptRanges` for purposes
where the original edges are correct (e.g. overlaying older results, PHENIX
comparisons, or shower-shape studies that use different binning).  A blind
search-and-replace risks breaking working plots.

### Verdict: ASK USER

The actual cross-section is unaffected.  Before touching `plotcommon.h`, verify:
(a) no plotting macro creates new histograms from `ptRanges` for arithmetic with
result histograms, and (b) the shower-shape configs that use 5 coarser bins are
also unaffected.  The safest approach may be to leave `plotcommon.h` as-is and
ensure plotting macros that need exact bin edges read them from the ROOT file.

---

## 2. Luminosity Values (5 Different Numbers)

**Claimed bug:** Five values appear for the 1.5 mrad period: 16.2735, 16.6,
16.8588, 16.9, and the analysis note says "16.9 pb^-1".

### Counter-argument

These are not five versions of the same number.  They serve different purposes
and may each be correct in context:

- **16.2735 pb^-1** (`config_bdt_nom.yaml`, `lumi` field): This is the value
  that enters the actual cross-section calculation via `CalculatePhotonYield.C`
  line 74.  It was computed by `LumiCalculator.C` summing Bit30Corr entries from
  the Joey lumi list for runs 51274-54000.  This is presumably the most precise
  calculation with the specific run range and vertex cut used in the analysis.

- **16.8588 pb^-1** (`config_showershape_1p5rad.yaml`): This config is for
  shower-shape studies, which may use a slightly different run range or vertex
  cut.  Shower-shape studies do not produce cross-section results -- they produce
  shape comparisons.  The lumi value here is used only for display labels, not
  for physics normalization.

- **16.6 pb^-1** (`plotcommon.h`, `strleg2_1`): This is a display string for
  plot legends.  It is a rounded value used in presentations.  It does not enter
  any calculation.

- **16.9 pb^-1** (analysis note, `analysis.tex` line 264 and `selection.tex`
  line 135): The note says "16.9 pb^-1" with asymmetric uncertainties
  (+1.4/-1.2).  Given those uncertainties are ~8%, the difference between 16.3
  and 16.9 (3.7%) is well within the luminosity systematic.  The note value may
  include a different MBD cross-section normalization or a different vertex cut
  than the config value.

The critical question is: **does the code use the right number?**  Yes.
`CalculatePhotonYield.C` reads from the config YAML, which has 16.2735.  The
analysis note's 16.9 is for human readers and carries stated uncertainties.

**Where it matters:** If the analysis note says 16.9 but the code uses 16.2735,
the paper result would be ~3.7% higher than what the note implies. This is
within the lumi systematic, but it would be misleading documentation.

### Risk of fixing

"Fixing" the lumi by making all values identical would be incorrect if they
legitimately correspond to different calculations.  Blindly setting 16.9
everywhere would change the cross-section result.  Blindly setting 16.2735
everywhere would put a too-precise number in a note where rounded + uncertainty
is standard.

### Verdict: ASK USER

The user should confirm which luminosity calculation is the intended final value
for the cross-section, then:
(a) Ensure `config_bdt_nom.yaml` has that value.
(b) Update the analysis note to match (with appropriate rounding and
uncertainties).
(c) Leave `plotcommon.h` with a rounded label.
(d) Leave `config_showershape_1p5rad.yaml` alone if it uses a genuinely
different run selection.

---

## 3. BDTinput.C et3_on/et4_on Bug

**Claimed bug:** Line 194 reads `et3_on` instead of `et4_on`:
```cpp
int et3_on = configYaml["analysis"]["et3_on"].as<int>(1);
int et4_on = configYaml["analysis"]["et3_on"].as<int>(1);  // should be et4_on
```

### Counter-argument

This bug is real but **has zero practical impact** for two independent reasons:

1. **Both et3_on and et4_on are 0 in the nominal config.**
   `config_bdt_nom.yaml` lines 70-71: `et3_on: 0`, `et4_on: 0`.  Since
   `et4_on` reads the `et3_on` YAML field (which is 0), and the correct value
   would also be 0, the result is identical.

2. **BDTinput.C is a training data extractor, not part of the analysis chain.**
   The `et3_on`/`et4_on` variables in BDTinput.C control the `nfail` count for
   the non-tight classification of training data.  Since the BDT nominal config
   uses `bdt_on: 1` and all shower-shape `_on` flags are 0 (the BDT subsumes
   them), this nfail logic is effectively dormant.  The training data selection
   is dominated by the BDT score cut, not individual shower-shape variables.

3. **BDTinput.C may be largely obsolete.**  It is step 2 in the pipeline
   (slimtree -> training text files), and the training has already been done.
   Unless re-training is planned with individual shower-shape cuts enabled, this
   code path is dead.

Even if `et4_on` were set to 1 while `et3_on` were 0, the impact would only be
on which clusters are classified as non-tight in the BDT training data.  This
would slightly change the training signal/background boundary, not the final
cross-section.

### Risk of fixing

Minimal risk -- the fix is a one-character change.  But it could create a false
sense of progress, and if someone re-runs BDTinput.C after "fixing" it, the
training data would be (trivially) different, potentially requiring re-training
and re-validation for no physics reason.

### Verdict: FIX (low priority)

Fix it as a code-quality improvement, but verify first that no downstream
process depends on the exact training data files currently produced by
BDTinput.C.  If the training data is frozen, this fix is cosmetic.

---

## 4. npb_score_cut Config Path Inconsistency

**Claimed bug:** `RecoEffCalculator_TTreeReader.C` reads
`analysis.common.npb_score_cut`, while `NPB_PurityStudy.C` and
`DoubleInteractionCheck.C` read `analysis.npb_score_cut` (top-level).

### Counter-argument

I verified that **no config file** has `npb_score_cut` at the top-level
`analysis` path.  All configs place it under `analysis.common.npb_score_cut`.
Therefore:

- `NPB_PurityStudy.C` line 227: `configYaml["analysis"]["npb_score_cut"].as<float>(0.5)` -- this YAML key is MISSING in every config, so the default value **0.5** is always used.
- `DoubleInteractionCheck.C` line 65: same situation, always falls back to **0.5**.
- `RecoEffCalculator_TTreeReader.C` line 494: reads `analysis.common.npb_score_cut`, which is explicitly **0.5** in every config.

**All three macros end up using the same value: 0.5.**  The inconsistent config
path is cosmetically wrong but produces identical physics results.

However, this is a latent bug.  If someone changes `npb_score_cut` in the
`common` section to a different value (say 0.6), `RecoEffCalculator_TTreeReader.C`
would use 0.6, but `NPB_PurityStudy.C` and `DoubleInteractionCheck.C` would
still use the default 0.5.  This would create a silent inconsistency.

The `make_bdt_variations.py` OVERRIDE_MAP routes `npb_score_cut` to
`analysis.common.npb_score_cut` (line 196), which means systematic variations
of the NPB cut would not propagate to the study macros.

### Risk of fixing

Changing the study macros to read from `analysis.common.npb_score_cut` is safe
and trivially correct.  The only risk is if someone intentionally wanted the
study macros to use a different cut than the main analysis -- but there is no
evidence of this, and the hardcoded `npb_score_abcd_cut = 0.5f` on line 228
of NPB_PurityStudy.C suggests the author simply expected 0.5 everywhere.

### Verdict: FIX

This is a real latent bug with zero current impact.  Fix `NPB_PurityStudy.C`
and `DoubleInteractionCheck.C` to read from `analysis.common.npb_score_cut`.
Verify that the default fallback remains 0.5 in case old configs lack the field.

---

## 5. non_tight BDT Config Mismatch (DoubleInteractionCheck.C)

**Claimed bug:** `DoubleInteractionCheck.C` reads only flat `bdt_max` and
`bdt_min` for non_tight (no `_slope`/`_intercept`), while
`RecoEffCalculator_TTreeReader.C` uses parametric
`bdt_max_slope * ET + bdt_max_intercept`.

### Counter-argument

In `RecoEffCalculator_TTreeReader.C`, the non_tight BDT upper bound is:
```cpp
float non_tight_bdt_max_slope = configYaml[...]["bdt_max_slope"].as<float>(0);
float non_tight_bdt_max_intercept = configYaml[...]["bdt_max_intercept"].as<float>(non_tight_bdt_max);
```
The default for `bdt_max_slope` is **0**, and the default for
`bdt_max_intercept` is `non_tight_bdt_max` (the flat value).

In the nominal config (`config_bdt_nom.yaml` lines 203-206):
```yaml
non_tight:
    bdt_max: 0.5
    bdt_max_slope: -0.015
    bdt_max_intercept: 0.80
    bdt_min: 0.02
```

So the nominal config explicitly sets parametric non_tight BDT.  At ET=20 GeV:
- `RecoEffCalculator`: `bdt_max = -0.015 * 20 + 0.80 = 0.50`
- `DoubleInteractionCheck`: `bdt_max = 0.5` (flat)

These happen to agree at ET=20.  But at ET=10: RecoEff gives 0.65, while
DoubleInteractionCheck gives 0.5.  At ET=30: RecoEff gives 0.35, while
DoubleInteractionCheck gives 0.5.

**This IS a real discrepancy** for the double-interaction study, especially at
high and low ET.  However, the question is whether this matters for the study's
conclusions.

The toy double-interaction study compares single-interaction vs double-interaction
ABCD yields at the **same** threshold.  Since the flat cut is applied
identically to both single and double interaction events, the **relative**
comparison (which is the point of the study) is still internally consistent.
The absolute purity numbers differ from the main analysis, but the study is
asking "how much does pileup shift the purity?" -- and the answer to that
question does not depend strongly on the exact BDT threshold, because pileup
degrades shower shapes regardless of where you draw the line.

### Risk of fixing

Adding parametric BDT to `DoubleInteractionCheck.C` requires:
(a) Adding slope/intercept config reading (4 new lines).
(b) Changing two threshold calculations (tight and non_tight).
(c) Re-running the double-interaction study.
(d) Re-generating all shower-shape comparison plots.

This is a nontrivial validation effort.  If the study is already complete and
the conclusions are incorporated into the analysis note, re-running it for a
second-order correction to an already second-order study (pileup effects) may
not be worth the effort.

### Verdict: ASK USER

If the double-interaction study results feed quantitatively into the final
systematic uncertainty, fix it and re-run.  If the study is qualitative (showing
that pileup effects are small), the flat threshold is adequate and the internal
self-consistency of the comparison is more important than matching the main
analysis exactly.

---

## 6. Missing MBDxsec Bibliography Entry

**Claimed bug:** `\cite{MBDxsec}` in `systematics.tex` line 126 has no matching
entry in `cite.bib`.

### Counter-argument

This will produce a LaTeX warning ("Citation `MBDxsec' undefined") and render
as "[?]" in the compiled PDF.  It is a real omission.

Checking the context: the citation appears in `\subsection{Luminosity}`
(`systematics.tex` line 123-128), which states the MBD cross-section is
"25.2^{+2.3}_{-1.7} mb".  This is a final-result number that requires a
reference.  The section is active (not commented out), and the luminosity
systematic is one of the dominant uncertainties.

There is no indication this section will be removed before publication.  The
analysis note structure shows it is integral to the systematics discussion.

### Risk of fixing

Adding a bib entry for the MBD cross-section measurement is zero-risk.  The
only question is which paper to cite.  This is likely the sPHENIX MBD
cross-section measurement, which may be an internal note or preliminary result.
The user needs to provide the correct reference.

### Verdict: FIX

Add the `MBDxsec` entry to `cite.bib`.  The user must supply the correct
bibliographic reference (author, title, journal/note number, year).

---

## 7. Flat BDT Threshold in DoubleInteractionCheck.C (Tight Cut)

**Claimed bug:** `DoubleInteractionCheck.C` uses flat `tight_bdt_min` instead of
parametric `tight_bdt_min_slope * ET + tight_bdt_min_intercept`.

### Counter-argument

This is essentially the same issue as #5 but for the tight cut.  The nominal
config has `bdt_min: 0.7`, `bdt_min_slope: -0.015`, `bdt_min_intercept: 0.80`.

At ET=20 GeV: parametric gives `0.80 - 0.015*20 = 0.50`, while the flat value
is `0.7`.  So the flat cut is actually **tighter** than the parametric one at
ET=20.  At ET=10: parametric = 0.65, flat = 0.7 (flat is tighter).  At ET=30:
parametric = 0.35, flat = 0.7 (flat is much tighter).

The flat threshold being tighter than the parametric one means the
double-interaction study operates in a **more restrictive regime** than the main
analysis.  This is conservative for a pileup study -- if pileup effects are
small under a tighter cut, they are even smaller under a looser (parametric) cut.

Moreover, the study's core question is: "does the BDT classification change when
a cluster's kinematics are shifted by pileup?"  The answer depends on how far
the BDT score moves relative to the threshold, which is primarily a property of
the BDT model and the kinematic shift, not the threshold position.  Using a flat
threshold simplifies the interpretation: a cluster that was tight under a flat
cut but becomes non-tight after pileup shift has genuinely degraded, without
confounding ET-dependent threshold effects.

### Risk of fixing

Same as #5.  Additionally, switching to parametric in the study would make the
threshold looser at high ET, which could change the study's conclusions about
pileup effects in the high-ET region.  If those conclusions are already quoted
in the analysis note, changing them could create inconsistencies.

### Verdict: LEAVE (or combine with #5 if that is fixed)

The flat threshold is arguably better for this specific study because it provides
a cleaner separation of "pileup effect on BDT score" from "ET-dependent
threshold effect".  If #5 is fixed for the non-tight cut, fix this too for
consistency, but do not fix this in isolation.

---

## 8. N_PT_BINS = 10 in Python Report Scripts

**Claimed bug:** `make_selection_report.py` and `make_comparison_report.py`
hardcode `N_PT_BINS = 10`, missing the last 2 of 12 pT bins.

### Counter-argument

The nominal config has 12 bins: `[8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28,
32, 36]`.  The last two bins cover 28-32 and 32-36 GeV.

These high-pT bins typically have very low statistics in pp at 200 GeV.  The
PHENIX measurement (overlaid in `plot_final_selection.C`) extends only to ~25
GeV.  The `pTmax` variable in `plotcommon.h` is set to 30 GeV (line 9), and
`plot_final_selection.C` sets `upperx = 32` (line 44).  This suggests the
analysis intentionally de-emphasizes bins above 30 GeV.

The report scripts generate LaTeX documents with figures arranged in a grid.
Including 12 pT bins instead of 10 means 2 extra pages of mostly-empty
isolation ET plots.  The developer may have intentionally stopped at 10 bins
because:
(a) The last 2 bins have too few events for meaningful isolation ET shapes.
(b) The purity fit fails or is unreliable above 28 GeV.
(c) The report is used for internal review, not publication, and the high-pT
    bins are inspected separately.

However, the comment in the code says `# matches NptBins in plotcommon.h` but
NptBins is 12, not 10.  This makes it more likely that 10 was a mistake from
an earlier iteration when there were fewer bins, and the comment is stale.

The plot_isoET.C macro that generates the actual figures iterates over all
NptBins=12, so the figures for bins 10 and 11 exist on disk.  The report
scripts simply do not include them.

### Risk of fixing

Changing to 12 could cause the report to include figures that:
(a) Have too few entries to be meaningful (empty or noisy iso ET distributions).
(b) Show purity fits that failed or have huge error bars.
(c) Make the report longer without adding useful information.

None of these would break anything, but the report quality might decrease.

### Verdict: FIX with verification

Change `N_PT_BINS = 10` to `N_PT_BINS = 12`, but first check that the figures
for bins 10 and 11 actually exist and are non-empty.  If the underlying ROOT
histograms have zero or near-zero entries in those bins, the report scripts
should either skip empty bins gracefully or the user should confirm they want
them included.

---

## Summary Table

| # | Issue | Current Impact | Verdict |
|---|-------|---------------|---------|
| 1 | pT bin mismatch | None (plotcommon.h not used in cross-section chain) | ASK USER |
| 2 | Luminosity scatter | None if config value is correct; note discrepancy is documentation | ASK USER |
| 3 | BDTinput.C et3/et4 | Zero (both flags are 0 in config; training is frozen) | FIX (low priority) |
| 4 | npb_score_cut path | Zero (default 0.5 matches config value 0.5) | FIX |
| 5 | non_tight BDT flat | Study is internally consistent; conservative | ASK USER |
| 6 | Missing MBDxsec bib | Renders as [?] in PDF | FIX |
| 7 | Tight BDT flat in study | Tighter than parametric; conservative for pileup study | LEAVE |
| 8 | N_PT_BINS = 10 | Reports miss 2 high-pT bins | FIX (with verification) |

### Issues with zero current physics impact: 1, 3, 4, 7
### Issues that are documentation/cosmetic: 2, 6, 8
### Issues requiring user judgment: 1, 2, 5
