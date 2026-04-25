# Physics Review: Converted-Photon Energy-Response and Resolution Study

**Date:** 2026-04-16  
**Reviewer:** Wave2-C  
**Deliverable:** Wave1-C converted-photon response analysis  
**Status:** PASS (12/13 checks pass; 1 minor WARN properly flagged)

---

## Overview

This review covers the energy-response (R_E = cluster_E/particle_E) and resolution (σ) characterization of converted vs unconverted photons in sPHENIX pp at 200 GeV, using single-particle MC (photon10). The analysis matches reconstructed CEMC clusters to truth photons via truthtrkID, stratifies by conversion category, and measures response shapes, means, widths, tail fractions, angular displacements, and matching efficiency as functions of pT and η.

**Key Physics Finding:** Converted photons exhibit a pronounced left-tail response distribution (mean response -7.6% vs unconverted) caused by magnetic-field bending of e+/e- pairs. At high pT (>25 GeV), the pair collimation compresses the tail and the converted response asymptotically exceeds the unconverted (non-linear EMCal degradation effect).

---

## Detailed Check Results

### 1. Response Definitions Correct ✓ PASS

**Requirement:** R_E = cluster_E / particle_E (not divided by particle_Pt by accident). R_ET = cluster_Et / particle_Pt.

**Check:**
- `analyze_response.py` line 292–293:
  ```python
  R_E = matched_E / photon_E
  R_ET = matched_ET / photon_pt
  ```
- Verified against reference: `particle_E` and `particle_Pt` from CaloAna24.h are energy and transverse momentum respectively.
- Summary values match: unconv <R_E> = 0.9424, conv = 0.8666 (from summary.txt).

**Result:** ✓ PASS — Definitions are physics-correct and calculations are accurate.

---

### 2. Matching Criteria Correct ✓ PASS

**Requirement:** Cluster matched to truth photon via `cluster_truthtrkID_CLUSTERINFO_CEMC == particle_trkid`, highest cluster_Et, with cluster_Et > 3 GeV threshold.

**Check:**
- `analyze_response.py` line 228–276:
  - Line 243: Cartesian match via `pairs["p"] == pairs["c"]` (trkid equality).
  - Line 252: `argmax(c_ET_masked, axis=2)` selects highest-ET cluster per photon.
  - Line 271: `match_valid = matched_ET > MATCH_ET_MIN` with `MATCH_ET_MIN = 3.0` (line 50).
- Cluster suffix: `_CLUSTERINFO_CEMC` (line 52, 235–239).

**Result:** ✓ PASS — Matching logic is correct and threshold is applied consistently.

---

### 3. Sample Selection Correct ✓ PASS

**Requirement:** `particle_pid==22`, `|particle_Eta|<0.7`, `particle_Pt ∈ [10, 34]` GeV.

**Check:**
- `analyze_response.py` line 207–220:
  - Line 216: `pid == PID_PHOTON` with `PID_PHOTON = 22` (line 49).
  - Line 219: `np.abs(eta) < ETA_ABS_MAX` with `ETA_ABS_MAX = 0.7` (line 48).
  - Line 217–218: `pt >= PT_LO and pt < PT_HI` with `PT_LO = 10.0, PT_HI = 34.0` (line 47).
- Sample statistics (findings.md line 25–30): 3.56M unconverted, 0.75M converted, 0.001M bad-secondary in selection.

**Result:** ✓ PASS — Selection window matches specification exactly.

---

### 4. pT Binning Correct ✓ PASS

**Requirement:** Canonical edges [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36] GeV covering [10, 34] GeV.

**Check:**
- `analyze_response.py` line 39:
  ```python
  PT_BINS = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
  ```
- Confirmed in summary.txt and per-pT results in response_findings.md (Section 4).
- 11 bins in the [10, 34] GeV range (analysis window).

**Result:** ✓ PASS — Binning is canonical and analysis range is covered.

---

### 5. Core σ Extraction: Gaussian Fit Stability ⚠ WARN

**Requirement:** Gaussian fit restricted to [peak − 0.2, peak + 0.1] window should be stable; low-stat bins may show artifacts.

**Check:**
- `make_plots.py` line 165: `gauss_fit_core(counts, centers, window=(-0.2, 0.1))` for asymmetric window.
- `make_plots.py` line 166: `gauss_fit_core(counts, centers, window=(-0.08, 0.08))` for narrow-core window.
- **Concern identified in findings.md Section 5 (Concerns #2):** Narrow-core sigma for converted photons **spikes to 0.094 at pT = 22–24 GeV** (vs. ~0.070 at adjacent bins; see response_sigma_vs_pt.pdf left panel).
  - Root cause: The ±0.08 window contains only ~16 bins of 0.01 width. For ~3000 counts in this bin, Gaussian fit can be pulled by statistical fluctuations.
  - **Finding properly flagged:** "would benefit from a fixed-window fit (e.g. 0.93–1.05) or a double-Gaussian + exponential tail model."
- **Impact:** The narrow-core sigma is elevated but the asymmetric sigma (task spec window) is stable. Quantitative claims in downstream use should cite the asymmetric window, not the narrow-core.

**Result:** ⚠ WARN (acknowledged) — Narrow-core fit shows one low-stat artifact at pT=22–24 GeV; asymmetric window (task spec) is stable. Authors properly identified concern. **Recommendation:** For pT = 22–24 GeV, cite asym sigma (0.1051) rather than narrow-core (0.0703±spike) if high precision needed in that bin.

---

### 6. Asymmetric-σ Definition Matches Task Spec ✓ PASS

**Requirement:** Asymmetric-sigma fit window [peak − 0.2, peak + 0.1] per task spec.

**Check:**
- `make_plots.py` line 165: `window=(-0.2, 0.1)` is correctly hardcoded.
- Numerical result in summary.txt: unconverted asym sigma = 0.0826, converted = 0.1051.
- Findings.md Section 5 reports: "asym sigma converted = 0.1051 vs unconverted 0.0826" → **match verified**.
- Broadening calculation: (0.1051 − 0.0826) / 0.0826 = **27.2%** vs. reported "27%" → **match confirmed**.

**Result:** ✓ PASS — Definition, window, and reported values all correct.

---

### 7. Tail-Fraction Arithmetic ✓ PASS

**Requirement:** R < 0.8 tail fraction = 31% converted / 11% unconverted (cross-check from summary).

**Check:**
- Summary.txt: unconv R<0.8 = 0.1067, conv = 0.3106.
- Ratio: 0.3106 / 0.1067 = **2.91x** vs. findings claim "~3x more likely (0.311 vs 0.107)" → **match confirmed**.
- R < 0.5 tail: unconv = 0.0053, conv = 0.0201 → ratio = **3.79x** vs. findings "~4x more likely" → **match confirmed**.
- Per-pT crossover (findings.md Section 6): converted tail drops below unconverted near pT ~ 26 GeV, consistent with high-pT pair collimation.

**Result:** ✓ PASS — Tail fractions are arithmetically correct and physically consistent with pT trends.

---

### 8. Match-Efficiency Definition ✓ PASS

**Requirement:** (# truth photons with cluster_Et > 3 GeV match) / (# truth photons passing selection) per category, per pT bin. Denominator must be post-selection, not pre-selection.

**Check:**
- `analyze_response.py` line 116–144:
  - `fill_truth()` counts all truth photons already in selection (line 116–121).
  - `fill_matched()` counts matched truth photons (line 123–180).
  - Both use the *same* selection (`pid==22, pT∈[10,34], |η|<0.7`) applied *upstream* at line 207–220.
- Summary.txt match efficiency (10–34 GeV integrated):
  - Unconverted: 2,464,250 / 3,561,433 = **69.2%** ✓
  - Converted: 495,721 / 749,069 = **66.2%** ✓
  - Matches findings.md Section 2: "Unconverted 69.2%, Converted 66.2%".
- Per-pT trend (findings.md Section 2): efficiency decreases ~8 pp from 69.6% → 61.5% (unconv) and 65.8% → 56.6% (conv) across [10, 36] GeV, consistent with asymmetric selection window (truth pT ∈ [14, 30], analysis pT ∈ [10, 34]).

**Result:** ✓ PASS — Efficiency denominator is correctly post-selection. ~3% gap (conv lower) is explained by B-field pair separation and Et threshold loss.

---

### 9. Angular Deltas: 2π Wrap Handling ✓ PASS

**Requirement:** dEta, dPhi defined with proper 2π wrap (e.g., `np.arctan2(sin(dphi), cos(dphi))` or `((dphi + π) % (2π)) − π`).

**Check:**
- `analyze_response.py` line 72–74:
  ```python
  def dphi_wrap(dphi: np.ndarray) -> np.ndarray:
      """Wrap phi difference to [-pi, pi]."""
      return (dphi + np.pi) % (2.0 * np.pi) - np.pi
  ```
- This is the **correct modular-arithmetic formula** for 2π wrapping, equivalent to `arctan2`.
- Applied at line 295: `dPhi = dphi_wrap(matched_phi - photon_phi)`.

**Result:** ✓ PASS — Angular delta wrapping is correctly implemented.

---

### 10. Bimodal dPhi for Converted Photons ✓ PASS

**Requirement:** Physically expected from B-field bending of e+ and e−. PDF should show two peaks; report locations.

**Check:**
- **Visual inspection of `figures/response_angular_deltas.pdf` (center panel, dPhi):**
  - **Unconverted** (black solid): Narrow Gaussian core around dPhi ≈ 0 with RMS < 0.01 rad.
  - **Converted** (red dashed): **Clear double-peak structure** with lobes at dPhi ≈ **±0.01 rad** and a dip at zero.
- **Physics:** Magnetic field bends e+ in one phi direction and e− in opposite direction. The cluster centroid shifts toward whichever lepton carries more energy, creating bimodal distribution with mean still ≈ 0 (symmetric w.r.t. field direction).
- Findings.md Section 8 states: "lobes at dPhi ~ +/- 0.01 and a dip at zero" → **matches PDF observation exactly**.
- dEta: No bimodality (solenoidal field bends only in phi). dR: Broader tail for converted (to ~0.15 rad).

**Result:** ✓ PASS — Bimodal dPhi signature is physically expected and visually confirmed in PDF at correct location (±0.01 rad).

---

### 11. Eta Dependence: EMCal Crack/Gap Structure ✓ PASS

**Requirement:** <R> is 0.02–0.03 lower at |η| < 0.35 than |η| > 0.35 for both categories (known EMCal dead material near η=0).

**Check:**
- Findings.md Section 7 (response_eta_dep.pdf):
  ```
  |eta| < 0.35 (central):  unconv <R> = 0.931, conv = 0.8545
  |eta| > 0.35 (forward):  unconv <R> = 0.955, conv = 0.8781
  Difference: 0.024 (unconv), 0.0236 (conv)
  ```
- Shift magnitude: Δ<R> ≈ **0.02–0.03** per category → **matches requirement**.
- **Physical interpretation:** The |η| < 0.35 region includes the central η=0 crack/seam (known sPHENIX EMCal dead-material region). Response degradation is **eta-independent** in structure (both categories affected similarly).
- Tail fractions also show "worse at center" pattern (conv tail 3.1x unconv in both central and forward).

**Result:** ✓ PASS — Eta dependence is consistent with known sPHENIX EMCal geometry (η=0 dead-material structure) and is properly interpreted.

---

### 12. Bad-Secondary Category: Quantitation and Claims ✓ PASS

**Requirement:** Only 271 matched → quantitative claims there are noise. Should be called out in findings.

**Check:**
- Summary.txt: bad category N matched = 271, N truth = 1,065 → eff = 25.5%.
- **Findings.md Section 1 (Concerns #1):** Explicitly states:
  - "Bad-secondary photons are rare (0.02%)" of total sample.
  - "zero above pT=18 GeV (not enough statistics)."
  - **Recommendation:** "Consider dropping the category from downstream analyses - it's physically negligible (0.02%)."
- Per-pT: Only [10–12], [12–14], [14–16], [16–18] bins have entries; [18+] GeV bins have zero matches.
- The figures appropriately show bad-secondary (blue triangles) only where statistics permit (N > 500 threshold for plotting; see make_plots.py line 225–226).

**Result:** ✓ PASS — Bad-secondary category is properly quantified, flagged as noise, and recommendations for handling are clearly stated. No inflated claims made.

---

### 13. PDF Rendering Quality ✓ PASS

**Requirement:** Spot-check 1–2 PDFs to ensure proper rendering.

**Check:**
- **`response_shape_grid.pdf`:** 3×4 grid of per-pT R_E distributions (pT=[10–36] GeV). Unconverted (black solid) and Converted (red dashed) overlaid per bin. **Renders correctly.** Clear visual of left-tail growth in converted (low pT) → collapse at high pT.
- **`response_angular_deltas.pdf`:** Three panels (dEta, dPhi, dR) with log-y scale. Bimodal dPhi for converted clearly visible. **Renders correctly.**
- **`response_sigma_vs_pt.pdf`:** Two panels (narrow-core and asymmetric sigma vs pT). Spike at pT=22–24 GeV in narrow-core visible. **Renders correctly.**
- All axis labels, legends, and simulation watermarks ("sPHENIX Simulation, photon10, |η|<0.7") present and correct.

**Result:** ✓ PASS — PDFs render correctly with proper physics information and visual clarity.

---

## Summary Table

| Check | Item | Status | File:Line |
|-------|------|--------|-----------|
| 1 | Response definitions (R_E, R_ET) | PASS | analyze_response.py:292–293 |
| 2 | Matching (trkid, highest-ET, Et>3 GeV) | PASS | analyze_response.py:228–276 |
| 3 | Sample selection (pid, eta, pT) | PASS | analyze_response.py:207–220 |
| 4 | pT binning [8,10,...,36] | PASS | analyze_response.py:39 |
| 5 | Core σ extraction (fit stability) | WARN | make_plots.py:60–100; findings.md:150–151 |
| 6 | Asymmetric σ definition | PASS | make_plots.py:165; findings.md:121 |
| 7 | Tail-fraction arithmetic | PASS | findings.md:162–168; summary.txt |
| 8 | Match-efficiency definition | PASS | analyze_response.py:116–144 |
| 9 | Angular deltas (dEta, dPhi wrap) | PASS | analyze_response.py:72–74,294–295 |
| 10 | Bimodal dPhi for converted | PASS | figures/response_angular_deltas.pdf |
| 11 | Eta dependence (crack structure) | PASS | findings.md:178–203; figures/response_eta_dep.pdf |
| 12 | Bad-secondary quantitation | PASS | findings.md:252–259; summary.txt |
| 13 | PDF rendering | PASS | figures/response_*.pdf |

---

## Overall Physics Assessment

### Verdict: **PASS** ✓

The Wave1-C converted-photon response analysis is **physics-sound and ready for downstream use**. The analysis correctly measures the energy response and resolution of converted vs. unconverted photons, with proper handling of matching, selection, binning, and angular definitions. All numerical results are consistent between scripts, findings, and summary outputs.

### Strengths

1. **Comprehensive sample:** 4.3M truth photons (3.56M unconv, 0.75M conv, 0.001M bad) provides excellent statistical power.
2. **Proper physics setup:** Magnetic-field pair bending causes left-tail response ⟹ -7.6% mean shift and bimodal dPhi signature. Both effects observed and correctly interpreted.
3. **Rigorous concern audit:** Authors identify 7 potential issues (narrow-core fit stability, low-stat bad category, Et threshold systematic, BDT shape bias, ABCD conv/unconv fraction mismatch, high-pT crossover validation, vertex smearing) and recommend follow-ups. Shows critical analysis.
4. **Well-documented:** Findings markdown is detailed, per-pT and inclusive numbers match summary, and PDFs are clear.

### Cautions (Not Blockers)

1. **Narrow-core fit artifact at pT=22–24 GeV:** The ±0.08 window is marginal for Gaussian fits at low stats. **Mitigation:** Use asymmetric window (task spec) for pT=22–24 bin if high precision required. Authors acknowledge this.
2. **Et > 3 GeV threshold:** The ~3% efficiency gap (conv vs unconv) is a *lower bound* on the true reconstruction-efficiency systematic in the main analysis pipeline. PPG12 should investigate separately whether converted-photon efficiency loss is larger at tighter Et cuts or with BDT selection.
3. **BDT shape not studied here:** Converted photons have distorted showers (leptons bend opposite phi). BDT scores for conv vs unconv likely differ near tight threshold. Recommend follow-up: produce matched conv/unconv BDT-score distributions.

### Recommendations for Downstream Use

1. **Cross-section analysis:** Apply this response study separately for conv/unconv categories in the isolated-photon yield extraction. The 27% wider asymmetric-σ for converted photons affects unfolding stability.
2. **ABCD purity method:** Check that the effective conversion fraction in the tight-iso region matches the truth 17.4% once reconstruction efficiency is applied. If converted-photon efficiency is 3% lower (as shown here), the tight-iso sample may be biased toward unconverted.
3. **BDT training:** Consider training separate BDT classifiers for conv vs unconv if shower-shape features strongly correlate with conversion status.

---

## Files Reviewed

- `/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/converted_photon_study/analyze_response.py`
- `/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/converted_photon_study/make_plots.py`
- `/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/converted_photon_study/response_findings.md`
- `/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/converted_photon_study/rootFiles/response_converted.root` (histograms)
- `/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/converted_photon_study/rootFiles/summary.txt` (numeric summary)
- `/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/converted_photon_study/figures/response_*.pdf` (6 PDFs)

---

**Review completed:** 2026-04-16  
**Reviewer:** Wave2-C  
**Status:** APPROVED FOR DOWNSTREAM USE with noted cautions.
