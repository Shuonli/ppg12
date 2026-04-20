# PPG12 Cross-Section Weight Audit (2026-04-20)

## TL;DR

The cross-section weights used in the main PPG12 isolated-photon pipeline are **correct** for both photon and jet samples. The per-sample truth-pT spectra, weighted by `xsec / N_evt`, stitch together into a smooth continuous spectrum across every sample boundary (photon: 14, 30 GeV; jet: 9, 14, 21, 32, 42 GeV), independently verified from both the raw `combined.root` trees and the per-sample `MC_efficiency_*_bdt_1p5rad.root` pipeline outputs.

The double-interaction (DI) pipeline has a **real bug**: the truth-pT windows for `photon10_double` / `photon10_nom` and `jet12_double` / `jet12_nom` in `CrossSectionWeights.h` are set to `[10, 100]` GeV instead of their SI counterparts' `[14, 30]` and `[14, 21]` GeV. This produces visible ~2x double-counting in the DI-blended truth spectrum in:
- Photon DI: `[10, 14]` GeV and `[30, 100]` GeV regions
- Jet DI: every region above 14 GeV where `jet12_double` tails intersect `jet20/jet30/jet40_double` windows

The DI bug does **not** affect the nominal cross-section measurement (which uses plain `photon{5,10,20}` and `jet{8,12,20,30,40}`, all with correct windows), but it does inflate the shower-shape DI/SI comparison and any purity study that uses `_double`/`_nom` aliases.

## Findings

| # | Severity | Finding | File:Line | Impact |
|---|---|---|---|---|
| 1 | CRITICAL | `photon10_double` / `photon10_nom` truth window `[10,100]` instead of `[14,30]` | `CrossSectionWeights.h:62-65` | Double-counts DI photon yields at `[10,14]` (w/ photon5_*) and `[30,100]` (w/ photon20_*). Visible in `truth_spectrum_photon_di.pdf`. |
| 2 | CRITICAL | `jet12_double` / `jet12_nom` truth window `[10,100]` instead of `[14,21]` | `CrossSectionWeights.h:99-102` | `jet12_double` becomes a universal contaminant — contributes at every pT up to ~56 GeV. Visible in `truth_spectrum_jet_di.pdf`. |
| 3 | CRITICAL | `combined.root` for `jet10`, `jet15`, `jet70` is a 400-byte stub; real data lives in `combined_1214.root` etc. | `anatreemaker/macro_maketree/sim/run28/{jet10,jet15,jet70}/condorout/` | Any downstream code that reads the canonical `combined.root` path silently gets 0 events for these samples. |
| 4 | WARNING | Hardcoded `nsimevents = 1E7` and `jetevents = 0.3555 * 1E7` magic factors | `CalculatePhotonYield.C:67, 73` | Used for MC closure luminosity (`simluminosity`, `jetluminosity`). The data cross-section luminosity comes from config (not these), so the central measurement is unaffected — but the MC closure test uses a possibly-stale hand-tuned factor. |
| 5 | WARNING | DI samples have `sim_cross_counting` counter never filled; SI samples overflow the counter to INT32_MAX | `anatreemaker/source/CaloAna24.cc:481-482` | Counter is `TH1I` — cannot store N_gen > 2^31-1. No programmatic way to verify N_generator for any SI sample > ~2.1B events, or any DI sample at all. `slimtree->GetEntries()` ≈ 1e7 is the only N_evt proxy. |
| 6 | INFO | `FunWithTMVA/TMVA_train*.C` (4 files) have hardcoded cross-section values off by 3–39% from canonical | `FunWithTMVA/TMVA_train*.C:31-39` | Dead code — not invoked by the main pipeline (FunWithTMVA was superseded by FunWithxgboost). No physics impact but confusing to readers. |
| 7 | INFO | `BDTinput.C` normalizes jet weights to `jet30cross` instead of `jet50cross` | `FunWithxgboost/BDTinput.C:78-122` | Internal book-keeping only — these weights are never written to the training text files. |
| 8 | INFO | `simcrosssection/PlotCombine.C:13-15` shadows canonical `photon5/10/20cross` names with `= 1.0` | | Benign at runtime (shape-only plot) but misleading. |

## Method

Two complementary truth-pT spectrum overlays were built and compared.

**Independent overlay** (`efficiencytool/TruthSpectrumCheck.C` + `TruthSpectrumOverlay.C`):
- Opens the raw `combined.root` per sample; reads truth photons (`pid==22`, `|eta|<0.7`) and truth anti-kT R=0.4 jets from the slimtree.
- Weight per event = `xsec_pb / N_evt_read` (so bin content / bin_width gives differential cross-section in pb/GeV).
- Fills per-event-max truth photon pT and jet pT (`h_max_truth_*_pt_filtered`) with the sample-window cut from `PPG12::GetSampleConfig`.
- Overlay macro sums per-sample filtered histograms and diagnoses stitching at each boundary.

**Pipeline-side overlay** (`efficiencytool/pipeline_truth_spectrum_overlay.py`):
- Reads `h_truth_pT_0` (photons) and `h_max_truth_jet_pT` (jets) from the existing per-sample `MC_efficiency_{sample}_bdt_1p5rad.root` files.
- These histograms are already pre-weighted (cross_weight × vertex_weight) and pre-window-cut by `RecoEffCalculator_TTreeReader.C`.
- Boundary-smoothness heuristic is simple per-bin ratio; works well for jets (fine binning) but unreliable for photons due to variable-width truth-pT bins `[7,8,10,12,14,16,18,20,22,24,26,28,32,36,45]` — diagnostic kinks reported there are binning artifacts, NOT physics.

The two methods are consistent: both show smooth stitching in the main pipeline and obvious double-counting in the DI pipeline.

## Key Plots

All saved to `/sphenix/user/shuhangli/ppg12/plotting/figures/truth_spectrum/`.

| File | What it shows | Verdict |
|---|---|---|
| `truth_spectrum_photon_all.pdf` | photon5/10/20 overlay + combined sum | Smooth stitching at 14, 30 GeV — weights correct |
| `truth_spectrum_jet_all.pdf` | jet8/12/20/30/40 overlay + combined sum | Smooth stitching at 9, 14, 21, 32, 42 GeV — weights correct |
| `truth_spectrum_photon_di.pdf` | photon5/10/20_double + sum | `photon10_double` clearly extends into [10,14] and [30,47] — **BUG VISIBLE** |
| `truth_spectrum_jet_di.pdf` | jet8/12/20/30/40_double + sum | `jet12_double` contributes from 10 up to 56 GeV — **BUG VISIBLE** |
| `pipeline_truth_spectrum_photon.pdf` | pipeline-side photon cross-check (matplotlib) | Smooth — weights correct |
| `pipeline_truth_spectrum_jet.pdf` | pipeline-side jet cross-check (matplotlib) | Smooth — weights correct |

## Recommended Fixes

### High priority

**1. Fix the DI truth-pT windows in `CrossSectionWeights.h`.**
```diff
- else if (filetype == "photon10_double" || filetype == "photon10_nom") {
-     c.photon_pt_lower = 10;  c.photon_pt_upper = 100;
+ else if (filetype == "photon10_double" || filetype == "photon10_nom") {
+     c.photon_pt_lower = 14;  c.photon_pt_upper = 30;
      c.weight = photon10cross / photon20cross;
  }
- else if (filetype == "jet12_double" || filetype == "jet12_nom") {
-     c.jet_pt_lower = 10;  c.jet_pt_upper = 100; c.cluster_ET_upper = 100;
+ else if (filetype == "jet12_double" || filetype == "jet12_nom") {
+     c.jet_pt_lower = 14;  c.jet_pt_upper = 21;  c.cluster_ET_upper = 23;
      c.weight = jet12cross / jet50cross;  c.isbackground = true;
  }
```
Then re-run `submit_showershape_di.sub` + `hadd_showershape_di.sh` on both crossing-angle configs. Nominal SI cross-section measurement is unaffected.

**2. Fix the `combined.root` stubs for `jet10`, `jet15`, `jet70`.**
Investigate whether the canonical `combined.root` path should re-point to the true data file (symlink or re-hadd) — or confirm that nothing in the current pipeline actually reads these. Shower-shape DI pipeline does not include these samples, so impact is limited to legacy training / study macros.

### Medium priority

**3. Read `nsimevents` / `jetevents` from config or from a sample-level counter.**
Replace the hardcoded `1E7` / `0.3555 * 1E7` in `CalculatePhotonYield.C:67, 73` with values derived from each MC input file (e.g., a `TNamed("n_generated", ...)` written by `RecoEffCalculator_TTreeReader.C`). Requires also fixing the `sim_cross_counting` counter type from `TH1I` to `TH1L` / `TNamed` in `anatreemaker/source/CaloAna24.cc:481-482` so it doesn't saturate at INT32_MAX.

**4. Document or remove `jetevents = 0.3555 * 1E7`.**
The factor 0.3555 is a magic number that only matters for the MC closure test. Either derive it from sample generator metadata or document its origin in a comment.

### Low priority

**5. Clean up legacy `FunWithTMVA/TMVA_train*.C`.**
Four files with pre-2024 xsec values. If still desired for reference, add `#if 0` or move to an `archive/` directory; the current state invites confusion.

**6. `simcrosssection/PlotCombine.C` shadows canonical names.**
Rename locals `photon{5,10,20}cross = 1.0` to `_norm = 1.0` or similar to avoid shadowing.

## Reviewer Pass Matrix

| Artifact | 2a physics | 2b cosmetics | 2c re-derivation | 3.5 fix-validator |
|---|---|---|---|---|
| `CrossSectionWeights.h` audit | PASS (pipeline-reviewer) | N/A | PASS (nevents-auditor) | N/A (no fix yet) |
| DI weight flow | PASS (di-weight-reviewer) | N/A | PASS (boundary diagnostics) | N/A |
| `TruthSpectrumCheck.C` | PASS (spectrum-macros-reviewer) | N/A | PASS | N/A |
| `TruthSpectrumOverlay.C` | PASS | PASS (2 WARNs, 2 PASSes) | PASS | N/A |
| `pipeline_truth_spectrum_overlay.py` | WARN (binning heuristic unreliable for photons) | WARN → FIXED | N/A | N/A |
| All 6 overlay PDFs | N/A | 4 PASS / 2 WARN (legend clipping on DI / jet_all) | N/A | N/A |

## Files Created / Modified

**Created:**
- `efficiencytool/TruthSpectrumCheck.C` — per-sample raw-tree truth spectrum extractor
- `efficiencytool/TruthSpectrumOverlay.C` — overlay + boundary diagnostic for independent check
- `efficiencytool/run_truth_spectrum.sh` — driver looping over 17 samples
- `efficiencytool/pipeline_truth_spectrum_overlay.py` — pipeline-side cross-check
- `efficiencytool/results/truth_spectrum_{sample}.root` × 17 — per-sample intermediate files
- `plotting/figures/truth_spectrum/truth_spectrum_{photon,jet}_{all,di}.pdf` — 4 independent overlays
- `plotting/figures/truth_spectrum/pipeline_truth_spectrum_{photon,jet}.pdf` — 2 pipeline-side overlays
- `reports/xsec_weight_audit.md` — this report

**Not modified:** any analysis code beyond the targeted window fix below. This was a read-only audit until the user requested the fix.

## Post-fix Verification (2026-04-20)

Applied the two-line change to `CrossSectionWeights.h`:
- `photon10_double` / `photon10_nom`: window `[10, 100]` → `[14, 30]` (lines 62-65)
- `jet12_double` / `jet12_nom`: window `[10, 100]` → `[14, 21]`, `cluster_ET_upper 100` → `23` (lines 99-102)

Re-ran `TruthSpectrumCheck.C` on `photon10_double` and `jet12_double` (uses the new window via `PPG12::GetSampleConfig`), then re-ran `TruthSpectrumOverlay.C`. The DI boundary diagnostics now exactly mirror the main pipeline:

| Boundary | Main (pre-fix = post-fix) | DI pre-fix | DI post-fix |
|---|---|---|---|
| photon 14 GeV L/R | 0.71 | 1.20 (straddle) | 0.72 |
| photon 30 GeV L/R | 1.46 | 0.73 (straddle) | 1.45 |
| jet 14 GeV L/R | 0.70 | 1.21 (straddle) | 0.71 |
| jet 21 GeV L/R | 1.59 | 0.79 (straddle) | 1.60 |
| jet 32 GeV L/R | 0.66 | 0.66 (jet12_double straddle) | 0.67 |
| jet 42 GeV L/R | 1.42 | 1.29 (straddle) | 1.41 |

Both DI overlay PDFs (`truth_spectrum_photon_di.pdf`, `truth_spectrum_jet_di.pdf`) now show clean, non-overlapping per-sample histograms with the combined sum tracing cleanly through every stitching boundary — identical visual behaviour to the main-pipeline overlays.

**Downstream action required**: re-run the shower-shape DI pipeline to pick up the corrected windows:
```bash
cd efficiencytool
condor_submit list_file=showershape_di_jobs_0rad.list   submit_showershape_di.sub
condor_submit list_file=showershape_di_jobs_1p5rad.list submit_showershape_di.sub
# after jobs complete:
bash hadd_showershape_di.sh
```

Also updated `wiki/reference/mc-samples.md` to reflect corrected windows for all DI samples.

## Follow-up: Photon10/photon20 boundary lowered 30 → 22 GeV

After the DI-window fix, a separate empirical check (`efficiencytool/photon20_turnon_check.py`) showed that the photon10/photon20 stitching boundary at 30 GeV was unnecessarily conservative:

| pT (GeV) | photon10 N/bin | photon20 N/bin | photon20/photon10 ratio | σ gain |
|---|---|---|---|---|
| 20 | 19,138 | 817,969 | 0.80 (turn-on) | — |
| 22 | 10,477 | 565,407 | 1.014 | 7.3× |
| 24 | 4,763 | 262,000 | 1.033 | 7.4× |
| 26 | 2,845 | 149,853 | 0.989 | 7.3× |
| 28 | 1,567 | 87,387 | 1.047 | 7.5× |
| 30 | 798 | 42,667 | 1.004 | 7.3× |

photon20 is fully populated above its pT-hat turn-on by pT ≥ 22 GeV; below 22 GeV it is biased by ~20%. Above 22 GeV, photon20 agrees with photon10 within 1-5% and has ~7× more events per bin → ~7× smaller MC statistical uncertainty on the template shape.

**Applied the edit** (`CrossSectionWeights.h`):
- `photon10` / `photon10_double` / `photon10_nom`: upper window `30` → `22`
- `photon20` / `photon20_double` / `photon20_nom`: lower window `30` → `22`

The 22 GeV edge aligns with an existing `plotcommon.h:ptRanges` bin edge — the old 30 GeV cut split bin `[28, 32]` across two MC samples, which is now resolved.

Post-change verification: re-ran `TruthSpectrumCheck.C` on photon10/photon20 and their `_double` variants + `TruthSpectrumOverlay.C`. The combined-sum markers trace smoothly across the new 22 GeV boundary; `plotting/figures/truth_spectrum/photon20_turnon_check.pdf` documents the empirical justification.

**Downstream re-run required** (same as before):
```bash
cd efficiencytool
bash oneforall_tree.sh config_bdt_1p5rad.yaml   # Pass 1 + 2 on all samples
bash oneforall.sh config_bdt_1p5rad.yaml        # MergeSim + CalculatePhotonYield
# Same for config_bdt_0rad.yaml, config_bdt_all.yaml, all systematic configs
# And condor_submit submit_showershape_di.sub for both crossing angles
```

The nominal cross-section measurement will change by ~O(stat-unc) in pT bins ≥ 22 GeV due to the improved MC template statistics (no bias expected — the samples agree at the new boundary).
