---
reviewer: critic
date: 2026-04-16
target: /gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/converted_photon_study.tex (820 lines)
verdict: APPROVE (revision round 1 landed; F-1 fixed at lines 631-633)
---

# Cross-verification of numbers, paths, and claims

## Summary

The report is numerically accurate: every inclusive and per-pT number I
cross-checked against the four findings.md files and the corresponding
ROOT/summary.pkl outputs matches to within rounding (<0.1% in all cases).
All 16 `\includegraphics` paths exist under
`converted_photon_study/figures/`.  Physics interpretation is internally
consistent across the four sections and matches the three reviewer MDs.

There is ONE genuine FAIL: a factually wrong feature name and line-range
citation pointing into `FunWithxgboost/config.yaml`. This propagates from
`showershape_findings.md` Sec. 5, so the finding doc is also wrong, but
since the .tex repeats it verbatim the .tex must be fixed.

## Numbered claim-by-claim verification

### PASS (verified against findings and/or ROOT histograms)

1. **Abstract** (lines 31--58):
   - 17.3% converted — findings: 17.317% ✓ (fractions)
   - $\Delta\langle R_E\rangle \simeq -0.076$ ($-8\%$) — summary.pkl gives
     `<R>unconv=0.9424`, `<R>conv=0.8666`, $\Delta=-0.0758$, $-8.04\%$ ✓
   - Iso pass fraction: $-14.6$ pp inclusive / $-22.8$ pp at 8--10 GeV —
     isolation: 0.649 $\to$ 0.503 ($\Delta-0.146$); 0.562 $\to$ 0.334 ($\Delta-0.228$) ✓
   - BDT mean shift $-16\%$, $w_\phi$ $+91\%$ — showershape ROOT inclusive
     means: BDT 0.7362 $\to$ 0.6157 ($-16.37\%$), $w_\phi$ 0.1275 $\to$ 0.2438 ($+91.22\%$) ✓

2. **Sec. 3 Fractions** (lines 162--240):
   - 855,380 truth photons, 17.317% converted, 0.024% bad — `fractions_findings.md` ✓
   - Table counts 707,047 / 148,129 / 204 ✓
   - Match efficiencies 67.987 / 66.290 / 39.024% ✓
   - Per-|eta| 15.369% at |eta|<0.1 and 20.124% at 0.6-0.7 ✓
   - "$\sim 32\%$ at $|\eta|\simeq 1.4$-$1.5$" — findings: 30.717% at 1.3-1.4, 32.482% at 1.4-1.5 ✓
   - "$\sim 2$~pp lower match efficiency" (caption, line 219) — actual 1.697 pp
     (rounds to 2, OK)
   - "relative $\sim 2.5\%$ deficit" (line 238) — actual -2.50% ✓

3. **Sec. 4 Response** (lines 244--374):
   - Table~\ref{tab:response_inclusive} (lines 298--322):
     $\langle R_E\rangle$ 0.9424 / 0.8666 ($-0.0758$, $-8\%$) ✓ (ROOT verified)
     median 0.9650 / 0.9050 ($-0.060$) ✓
     $p_{84}-p_{50}$: 0.080 / 0.120 ($+50\%$) ✓
     $p_{50}-p_{16}$: 0.120 / 0.230 ($+92\%$) ✓
     narrow-core $\mu$: 0.9956 / 0.9994 ✓ (summary.pkl fit_mu_core)
     narrow-core $\sigma$: 0.0669 / 0.0703 ($+5\%$) ✓
     asymmetric $\sigma$: 0.0826 / 0.1051 ($+27\%$) ✓
     tail $R<0.8$: 0.107 / 0.311 ($\times 2.9$) ✓ (ROOT: 0.1067 / 0.3106)
     tail $R<0.5$: 0.0053 / 0.0201 ($\times 3.8$) ✓
   - Per-$p_T$ shifts $-0.082$ / $-0.042$ / $-0.014$ / $+0.023$ — findings
     Sec. 4: $-0.0817$ / $-0.0416$ / $-0.0142$ / $+0.0230$ ✓
   - Crossover $\sim$26-27 GeV, $+0.023$ at 28-32 GeV ✓
   - $\Delta\phi$ bimodal lobes at $\pm 0.01$ rad ✓ (matches findings Sec. 8)

4. **Sec. 5 Isolation** (lines 377--481):
   - Eq. 1 parametric cut 0.502095 + 0.0433036·ET — matches
     `efficiencytool/config_bdt_nom.yaml:33-34` ✓
   - "nominal uses topo R=0.4 (use_topo_iso: 2)" — matches config line 53 ✓
   - Inclusive pass fractions 0.649 / 0.503 / $-0.146$ ✓ (findings + ROOT)
   - 8-10 GeV: 0.562 / 0.334 / $-0.228$ ✓
   - Table~\ref{tab:iso_decomposition} all rows ✓:
     iso_03 (0.816, 1.386, +0.570), iso_03_emcal (4.683, 5.151, +0.468),
     iso_03_hcalin (0.269, 0.213, $-0.056$), iso_03_hcalout (0.287, 0.253, $-0.033$),
     and matching R=0.4 rows
   - "$+0.57$ GeV shift driven by EMCal ($+0.47$ GeV)" rounded correctly ✓

5. **Sec. 6 Shower Shape** (lines 484--637):
   - Table~\ref{tab:showershape_inclusive} (lines 559--583):
     All 10 rows cross-checked against findings Sec. 3 table and/or ROOT
     mean histograms: wetacog 0.191/0.225, wphicog 0.128/0.244, e11/e33
     0.649/0.597, e32/e35 0.980/0.946, ET1 0.939/0.911, ET3 0.661/0.586,
     CNN 0.638/0.464, cluster_prob 0.565/0.415, merged_prob 0.804/0.673,
     BDT 0.736/0.616 — all ✓
   - Separation powers: $e_{32}/e_{35}$ = 2.33 σ, $w_\phi$ = 2.15 σ (top two) ✓
   - BDT pT table (lines 598--629): all 12 rows match findings Sec. 4 verbatim ✓
     ($-0.135/-18.4\%$ at 9 GeV, ..., $-0.168/-25.8\%$ at 34 GeV)
   - "$\mathrm{BDT} > 0.80 - 0.015 \cdot E_T$" cited against config lines
     154-156 — matches `config_bdt_nom.yaml` (bdt_min_intercept=0.80,
     bdt_min_slope=$-0.015$; note the config key is `bdt_min_*` not
     `tight_bdt_min_*`; the text correctly uses `bdt_min_intercept`) ✓

6. **Bad photons** statistical caveat raised in Sec. 1 (line 108-110) and
   Sec. 3 table (204 entries) ✓ — matches reviewer MDs `_review_response.md`
   and `_review_showershape.md`.

7. **Narrow-core WARN** acknowledged — the text on line 325-326 states
   "narrow-core $\sigma$ broadens by only $\sim 5\%$" (inclusive) and
   explicitly uses asymmetric $\sigma$ as the broadening figure of merit
   ($+27\%$). This addresses the pT 22-24 GeV narrow-core-fit WARN from
   `_review_response.md` by preferring the asymmetric number downstream.

8. **Topo-iso caveat** — Sec. 2 paragraph "Caveat on isolation"
   (lines 143--152) and Sec. 5 intro (lines 386--395) explicitly flag that
   the input file lacks topo-iso branches, calo-tower iso is used as a
   proxy, and the nominal cut is at R=0.4 topo. Caveat fully surfaced.

9. **Match-efficiency lower-bound caveat** — Sec. 3 (lines 232--239)
   explicitly states that the match efficiency is a lower bound because
   `cluster_truthtrkID` can be labelled with a post-shower daughter $e^\pm$
   instead of the parent photon. Matches `fractions_findings.md` Concerns.

10. **Figure paths** — all 16 `\includegraphics{converted_photon_study/
    figures/*.pdf}` references point to files that exist on disk
    (verified with `ls` for each path).

### FAIL

**F-1 (line 631-632): Incorrect config feature list and line range.**
The text reads:

  "The most-used BDT features (\wetacog, \wphicog, $\mathrm{ET1}$,
   $e_{11}/e_{33}$, $e_{t4}/e_{t1}$, see \texttt{FunWithxgboost/config.yaml}
   lines 63--87) ..."

Cross-checked against `FunWithxgboost/config.yaml`:

- Feature list lives at lines **48--73** (25 features under
  `data.features`), NOT lines 63--87.
- `e_{t4}/e_{t1}` (`et4_over_et1`) is **not a training feature at all** —
  grep confirms no occurrence in `config.yaml`.
- The features lines 63--87 actually enumerate are
  `e11_over_e22` through `e22_over_e53` (ratio features), which are
  different variables from the ones cited.

This factually wrong citation originates in `showershape_findings.md`
Sec. 5 (same line-range, same `et4/et1`); the writer propagated it
verbatim. Because the report cites the config for reproducibility, this
needs to be corrected. Suggested fix:

  "The most-used BDT features (\wetacog, \wphicog, $\mathrm{ET1}$,
   $e_{11}/e_{33}$, $e_{32}/e_{35}$, see \texttt{FunWithxgboost/config.yaml}
   lines 48--73) ..."

($e_{32}/e_{35}$ is an actual training feature at config line 59 and was
independently identified in the separation-power ranking as a leading
discriminator; substituting it for the non-existent $e_{t4}/e_{t1}$
makes the sentence both literally correct and physically consistent.)

### WARN (minor, writer may choose to keep or adjust)

**W-1 (line 91): line-range hint for CaloAna24.cc**
Text: "G4 shower particle list of each truth photon
(\texttt{CaloAna24.cc} lines $\sim 807$--843)"
I could not verify this line range from the `anatreemaker/source/CaloAna24.cc`
file in this review pass (not required by the task spec; the tilde prefix
is honest about its approximate nature). Writer may wish to double-check,
but this is not a hard FAIL — the `\sim` prefix already signals
approximation.

**W-2 (line 118-119): truth-pT window statement**
Text says "truth $p_T \in [14,30]$~GeV at generation".
`response_findings.md` Sec. 2 says "truth $p_T \in [10, 34]$" is the
reported analysis window (generator actually in $[14, 30]$, but selection
expands to $[10, 34]$). The report's formulation is consistent with the
generator-level statement and Sec. 4 consistently states the inclusive
window is "$10 \le p_T \le 34$" (line 300), so no contradiction — this is
correct. No change needed.

## Reviewer-MD cross-check

- `_review_showershape.md`: PASS, all 10 checks — matches ✓
- `_review_isolation.md`: PASS, no FAILs — matches ✓
- `_review_response.md`: PASS with one WARN on narrow-core σ at
  pT=22-24 GeV — report uses asymmetric σ as primary broadening figure
  (27%) and merely mentions narrow-core (5%) with proper framing. WARN
  adequately addressed. ✓

## Overall verdict

**APPROVE** (revision round 1)

Round-1 update: writer replaced `et4/et1` with `e32/e35` and corrected
the line range to 48--73; re-reading .tex lines 631-633 confirms the fix
is in. `e32/e35` is a real training feature at config.yaml line 59 and
also the top-ranked separation variable (2.33σ), so the sentence is now
both literally correct and physically stronger than before.

All FAIL/WARN items resolved. W-1 (CaloAna24 line range) was advisory
only and not blocking. Report is cleared for compilation (task #3).
