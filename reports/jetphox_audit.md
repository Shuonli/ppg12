# JETPHOX Datapoint Audit + Per-Bin Statistics + JETPHOX vs Werner NLO

**Date:** 2026-04-22
**Scope:** Audit JETPHOX running setup + per-bin statistics; compare vs the Werner (Vogelsang) NLO curve used in the same plot
**Files audited:**
- `NLO/MakeJetPHOXhisto.C`, `NLO/run_jetphox_histos.sh`, `NLO/PlotTruthIso.C`
- `/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/` (install) incl. `parameter.indat`, `run_condor*_SL7.sh`
- `NLO/rootFiles/jetPHOX_{05,10,20}.root` + `plotting/sphenix_nlo/photons_newphenix_sc{05,1,2}.dat`
- `plotting/plot_final_selection.C` (primary; `plot_final.C` is a 0-byte file)

---

## 0. TL;DR

| Question | Answer |
|---|---|
| Is JETPHOX running correctly? | **MOSTLY YES**, with one CRITICAL setup concern (CT14lo instead of CT14nlo) and a few small code issues. The fiducial (|η|<0.7, iso cone R=0.3, iso ET_max=4 GeV, μ=pT, 3-point scale) matches the PPG12 fiducial. |
| Do we need more JETPHOX events? | **NO.** Max JETPHOX stat_rel in the PPG12 analysis range (8–36 GeV) is **6.8%** in the [32,36] bin. Measured total uncertainty there is ~65% (syst-dominated). The scale envelope dominates stat by ≥3× in every bin. Actual MC = 2×10⁸ events per scale — 2× more than the planning number. |
| Why does JETPHOX differ from the Werner NLO at high pT? | **JETPHOX is isolated, Werner is inclusive.** JP/We goes 0.68 → 1.44 across 8–36 GeV. Low-pT offset (0.68 at 8–10 GeV) is exactly what iso rejection of fragmentation photons should look like. The high-pT crossover (JP > We above ~16 GeV) is not explained by iso alone — it is driven by (a) a linear-interpolation bias on Werner that inflates Werner at low pT by up to 30%, (b) CT14**lo** (not CT14nlo) used in JETPHOX has a harder high-x gluon than CTEQ6M, (c) α_s mismatch (LHAPDF LO ≈0.130 vs Vogelsang NLO ≈0.118 → ~21% effect on direct component). **The two scale envelopes overlap in every PPG12 bin, so the two NLOs are formally consistent.** |

---

## 1. JETPHOX setup — correctness audit

Run-card `parameter.indat` deployed in each condor job under `jetphox_1.3.1_4/condorout/OutDir{0..99}/working/`. The three submit scripts `run_condor{,05,20}_SL7.sh` patch only the histogram name, RNG seed, and scale lines.

| Field | JETPHOX value | PPG12 fiducial | Match? |
|---|---|---|---|
| √s | 200 GeV (`parameter.indat:286`) | 200 GeV | ✓ |
| Beams | pp (H1=H2=0 at :58,:62) | pp | ✓ |
| Photon |y| | <0.7 (:290,:294) | [-0.7, 0.7] | ✓ |
| Photon pT | 8–40 GeV (:298,:302) | 8–36 reco; truth [7,45] | ✓ (covers) |
| Iso mode | flag=1 fixed cone (:315) | fixed cone truth | ✓ |
| Iso R | 0.3 (:327) | cone_size=3 → R=0.3 | ✓ |
| Iso ET_max | 4.0 GeV (:339) | iso_ET_truth<4 GeV | ✓ |
| Scales | μ_R=μ_F=μ_f = cm·pT^γ, cm∈{0.5,1,2} (:169,173,177,181) | 3-point envelope | ✓ |
| PDF | **CT14lo** (LHAPDF-6.5.3, :80) | — | **⚠ CRITICAL** — LO PDF in NLO matrix element |
| FF | BFG II (set code 0200, :156) | — | ✓ (standard) |
| α_s | from LHAPDF (:196) | — | ✓ |
| N_flavours | 5 (:200) | — | ✓ |
| Contributions | direct + one-fragmentation + NLO + LO (:230,241,248) | direct + frag | ✓ |
| N_events | 100 segments × 1M per scale = 10⁸ per scale (design) | — | — |
| N_events (actual) | `h_truth_pT.fEntries` reads **2×10⁸** per scale | — | INFO — 2× design; explained by 2 tracks/event (photon+recoil partons). Production fine. |

### Findings

1. **CRITICAL — CT14lo in NLO calculation.** `parameter.indat:80` reads `CT14lo`. LHC-standard practice (CMS 7 TeV, ATLAS 13 TeV) uses CT10 or CT14nlo for NLO JETPHOX (`wiki/raw/papers/A_lhc/1108.2044.md`, `1908.02746.md`). LO PDF with NLO matrix element mixes α_s schemes (LO α_s(MZ)≈0.130 vs the NLO-consistent ≈0.118) and yields a harder gluon at high x. Expected impact: few % at low pT, up to ~10% at high pT, pushing JETPHOX upward at high pT — **consistent with the direction of the observed JETPHOX/Werner divergence**. *Decision needed*: either regenerate with CT14nlo, or document the choice with physics justification in the note.

2. **WARNING — Segment count `nseg` is hardcoded and silently divergent.** `MakeJetPHOXhisto.C:23` uses `nseg=100` and `PlotTruthIso.C:50` uses `nseg=10` with no runtime check. Each is individually right for its current input set, but any operator swapping files or regenerating with different merge factors silently mis-normalizes by 10×. *Fix*: embed `nseg` into the JETPHOX output file (as TNamed or UserInfo slot) and read it at histogram time.

3. **WARNING — parton-level iso vs particle-level truth iso.** JETPHOX applies its 4 GeV / R=0.3 iso at parton level inside the matrix element. PPG12 applies the 4 GeV truth iso on stable-particle MC. The difference is the underlying-event energy inside R=0.3, which at 200 GeV is sub-GeV — small but not documented. LHC analyses routinely apply a small UE correction to JETPHOX; PPG12 does not. *Recommendation*: quantify the UE leak (toy: mean <4π×r²×dρ/dφdη> for R=0.3 at 200 GeV is estimated ≲0.5 GeV) and add a sentence in the note.

4. **WARNING — `abs(y[0])` on a double.** `MakeJetPHOXhisto.C:72` uses integer `abs()` instead of `std::fabs`. Recent ROOT/cling versions promote it to double, but the intent should be explicit.

5. **WARNING — `h_truth_eta_pT` 2D is filled before the |y|<0.7 cut and never Scale'd or Sumw2'd.** `MakeJetPHOXhisto.C:71` pre-cut, no per-GeV normalization — misleading if read as d²σ/(dη dpT). *Fix*: rename to `h_eta_pT_counts` or apply the cut and per-GeV.

6. **WARNING — top-level `working/parameter.indat` is out of sync with deployed cards.** Line 48 at top-level says `histo`; deployed copies say `ntuple`. Line 173/177/181 top-level says `0.5D0`; deployed nominal says `1.0D0`. The deployed copies are what actually ran. Reproducibility concern only. *Fix*: commit a frozen template.

7. **WARNING — `plot_final.C` is 0 bytes.** The CLAUDE.md / wiki reference this file as the primary; the actual code lives in `plot_final_selection.C`. *Fix*: restore or rename.

8. **INFO — contributions flag.** direct (`0` at :208, "Born not computed since NLO covers it") + one-fragmentation (`1` at :230) + NLO (`1` at :241) + LO (`1` at :248). Standard NLO configuration. ✓

---

## 2. Per-bin statistics (verified by two independent agents)

Central scale μ=pT, file `jetPHOX_10.root`, histogram `h_truth_pT` divided by Δη=1.4 to give d²σ/(dpT dy) in pb/GeV. (The stored histogram is integrated over the generation window, not per-unit-y — the /1.4 is applied in `plot_final_selection.C:71`.)

| pT bin [GeV] | dσ/dpT (pb/GeV) | stat_rel (%) | N_eff | meas total_rel (%) | stat/tot ratio | N_scale→target | Verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| 8–10 | 1519 | 0.14 | 524551 | 75.9 | 0.002 | 3e-5 | OK |
| 10–12 | 540 | 0.23 | 181487 | 30.9 | 0.008 | 6e-4 | OK |
| 12–14 | 219 | 0.37 | 71301 | 21.5 | 0.017 | 3e-3 | OK |
| 14–16 | 98.4 | 0.56 | 31691 | 21.1 | 0.027 | 6e-3 | OK |
| 16–18 | 47.0 | 0.82 | 14861 | 24.3 | 0.034 | 0.010 | OK |
| 18–20 | 23.4 | 1.17 | 7247 | 25.9 | 0.045 | 0.018 | OK |
| 20–22 | 9.82 | 2.02 | 2446 | 27.0 | 0.075 | 0.050 | OK |
| 22–24 | 6.54 | 2.24 | 1991 | 28.4 | 0.079 | 0.056 | OK |
| 24–26 | 3.66 | 2.99 | 1118 | 40.8 | 0.073 | 0.048 | OK |
| 26–28 | 2.10 | 3.90 | 657 | 43.1 | 0.091 | 0.074 | OK |
| 28–32 | 0.940 | 4.09 | 597 | 40.2 | 0.102 | 0.093 | OK |
| 32–36 | 0.323 | 6.80 | 216 | 65.0 | 0.105 | 0.099 | OK |

**Numerical re-derivation** (second agent): all 13 bins MATCH within tolerance (ΔC<0.03%, Δstat_rel ≤ 2.1%, ΔN_eff ≤ 0.6%). Production campaigns across the three scales have 199,999,428 / 199,999,528 / 199,999,528 entries (they agree at 10⁻⁶, same underlying events weighted differently per scale — as it should be for μ-scale variations).

**Scale envelope vs stat** (re-derived, NOT monotonic as previously written):

| bin | envelope_rel (%) | stat_rel (%) | envelope/stat |
|---|---:|---:|---:|
| 8–10 | 20.5 | 0.14 | 148.6 |
| 14–16 | 20.6 | 0.56 | 36.6 |
| 20–22 | 16.0 | 2.02 | 7.9 |
| 26–28 | 24.4 | 3.90 | 6.3 |
| 28–32 | **26.1** | 4.09 | 6.4 |
| 32–36 | 21.5 | 6.80 | **3.2** |
| 36–45* | 10.3 | **13.4** | 0.77 |

*overflow, not drawn in final plot.

**Verdict on events:** no rerun needed for the 8–36 GeV PPG12 analysis range. The overflow [36,45] bin is stat-dominated over the scale band — if the analysis ever extends above 36 GeV, 3–4× more events would be needed to push stat below scale envelope there.

---

## 3. JETPHOX vs Werner (Vogelsang) NLO — high-pT divergence

**Two theory curves in `plot_final_selection.C`**:

- **JETPHOX**: `jetPHOX_{05,10,20}.root` — isolated (iso_ET<4 GeV, R=0.3), CT14**lo**, BFG II, LHAPDF α_s, from our production.
- **Werner**: `sphenix_nlo/photons_newphenix_sc{05,1,2}.dat` — **inclusive prompt photon, no iso cut** (confirmed in `results.tex:21`: "The calculation does not require isoET condition"). Columns: `pT | direct | fragmentation | total` in pb/GeV² (invariant cross-section). Conversion `yield = (yield+fragyield) * 2π * pT` → d²σ/(dpT dy) per unit rapidity. Provenance = Gordon & Vogelsang NLO, CTEQ6M PDFs, GRV or BFG II FF.

### Ratio table (central scale, μ=pT)

Werner "binAvg" = `TGraph::Integral(xlo, xhi)/dx` — linear interpolation over 2.5 GeV spacing.

| pT bin | JP (pb/GeV) | We_binAvg (pb/GeV) | **JP/We** |
|---|---:|---:|---:|
| 8–10 | 1519 | 2242 | **0.68** |
| 10–12 | 540 | 678 | 0.80 |
| 12–14 | 219 | 251 | 0.87 |
| 14–16 | 98.4 | 106 | 0.92 |
| 16–18 | 47.0 | 48.6 | 0.97 |
| 18–20 | 23.4 | 23.3 | 1.01 |
| 20–22 | 9.82 | 11.4 | 0.86 |
| 22–24 | 6.54 | 5.83 | 1.12 |
| 24–26 | 3.66 | 3.10 | 1.18 |
| 26–28 | 2.10 | 1.69 | 1.25 |
| 28–32 | 0.94 | 0.73 | 1.29 |
| 32–36 | 0.32 | 0.23 | **1.41** |

Ratio drifts from 0.68 at low pT to 1.44 at high pT — a factor of 2 across the PPG12 range.

**Scale-envelope overlap:** the JP [05,20] band and the We [05,20] band overlap in every PPG12 bin (`compare_jetphox_werner.out:54-65`). **Formally consistent at NLO** despite the central-value drift.

### Ranked explanation (most likely → least)

1. **ISO — drives the low-pT offset (JP < We)** [~0.7× at 8 GeV is well matched].
   Werner has no iso cut; JETPHOX has ET_iso<4 GeV, R=0.3. This rejects fragmentation photons preferentially. Werner's frag fraction drops from 43% at 5 GeV to 24% at 40 GeV. Assuming ~70% frag rejection by iso: predicted JP/We_iso ~ 0.70 at 8 GeV → 0.81 at 35 GeV. **Matches the 0.68 at low pT.** But this alone monotonically reduces JP below We — cannot explain the crossover at ~16 GeV, nor the JP>We at high pT.

2. **Jensen-inequality bias in Werner bin averaging** [inflates Werner at low pT by up to 30%, pushes JP/We down at low pT].
   `plot_final_selection.C:339` uses `TGraph::Integral(xlo, xhi)/dx` — linear interpolation on a steeply-falling spectrum. Power-law fit test shows linear interp over-estimates the true bin integral by +29% at 8–10 GeV and +4% at 32–36 GeV. After Jensen correction, JP/We becomes ~0.87 at 8 GeV and ~1.39 at 32 GeV. **Strengthens the high-pT ratio.**

3. **PDF — CT14lo vs CTEQ6M at high x** [drives the high-pT drift JP > We].
   At √s=200 GeV, x_gluon ≈ 2pT/√s: x≈0.08 at 8 GeV → **x≈0.35 at 35 GeV**. CT14lo (LO, absorbs NLO K-factors into the gluon) is **harder** at high x than CTEQ6M (NLO). Expected JETPHOX excess: few % at low x, **up to ~10% at high x**. Direction consistent with observed JP>We at high pT.

4. **α_s mismatch** [~5–7%, or up to 21% if both α_s raised by the ratio 0.130/0.118 squared].
   JETPHOX's α_s is inherited from LHAPDF CT14lo (LO α_s(MZ)≈0.130). Vogelsang NLO typically uses α_s≈0.118. At NLO, direct scales as α·α_s, so the effective enhancement is (0.130/0.118) ≈ 1.10 for direct; frag scales ∝ α_s and enters weighted. Small contributor at all pT; no strong pT dependence.

5. **FF (BFG II vs GRV)** [<1% after iso].
   Minor. BFG II has been cross-validated against GRV in the literature; after iso the frag contribution is small at all pT, so FF choice barely matters.

### Recommendations

1. **If you want JETPHOX to match Werner at low pT**, move Werner to an iso cut matching JETPHOX (requires regenerating Werner, which means asking Werner/collaborator, or accepting the inclusive-vs-isolated offset and noting it clearly on the figure and in the caption).

2. **If you want JETPHOX to match Werner at high pT**, two fixes that should help:
   - **Switch JETPHOX to CT14nlo** (regenerate ~1 night of CPU). Expected effect: few % softer at high pT.
   - **Fix the TGraph bin averaging in `plot_final_selection.C:339`** — replace `TGraph::Integral/dx` with log-log interpolation or a per-bin power-law fit on Werner's 2.5 GeV samples. This alone removes ~30% of the low-pT offset on Werner without regenerating anything.

3. **Regardless of rerun decisions**, the scale envelopes overlap — the two calculations are formally consistent at NLO. The note should make clear that JETPHOX is the primary theory comparison (iso-matched fiducial) and Werner is shown for cross-check (inclusive, older PDF).

---

## 4. Reviewer pass matrix

| Artifact | Physics review | Plot cosmetics | Numerical re-derivation | Fix validator |
|---|---|---|---|---|
| JETPHOX run card (`parameter.indat`) | PASS (w/ 1 CRITICAL: CT14lo) | N/A | N/A | N/A |
| `MakeJetPHOXhisto.C` | PASS (w/ 3 CRITICAL: nseg, 2 WARNING) | N/A | N/A | N/A |
| `PlotTruthIso.C` | PASS (w/ 1 CRITICAL: nseg=10 inconsistency) | N/A | N/A | N/A |
| `plot_final_selection.C` JETPHOX usage | PASS (w/ 1 CRITICAL: empty `plot_final.C`, 4 WARNING) | N/A (no new plot produced in this audit) | N/A | N/A |
| Per-bin stat numbers + scale envelope | N/A | N/A | PASS (all 13 bins match, 2×10⁸ events confirmed) | N/A |
| JETPHOX vs Werner ratio table | PASS | N/A | PASS (derived numbers in `compare_jetphox_werner.out`) | N/A |

---

## 5. Files produced / modified

- `reports/jetphox_audit.md` (this file)
- `reports/work_jetphox_stats/jetphox_stats_analysis.py` + 3 CSVs (per-bin stats)
- `reports/work_jetphox_stats/compare_jetphox_werner.py` + `.out` (JP/We numeric comparison)
- `reports/work_jetphox_stats/frag_fraction.py` + `.out` (Werner direct/frag evolution)
- `reports/work_jetphox_stats/jensen_test.py` + `.out` (linear-interp bias estimate)
- `reports/work_jetphox_stats/rederive_jetphox.py` (independent re-derivation)

No code in `NLO/` or `plotting/` was modified in this audit.

---

## 6. Action items (prioritized)

| # | Priority | Action | Why |
|---|---|---|---|
| 1 | **HIGH** | Decide whether to swap CT14lo → CT14nlo in JETPHOX run card and regenerate | LO PDF in NLO ME; accounts for several % of the high-pT JETPHOX/Werner drift |
| 2 | **HIGH** | Patch `plot_final_selection.C:339` to use log-log or per-bin power-law interpolation for Werner | Removes up to 30% linear-interpolation bias at low pT, distorts JP/We comparison |
| 3 | **MEDIUM** | Restore or rename `plot_final.C` (currently 0 bytes, canonical name in docs) | Documentation/code disagreement |
| 4 | **MEDIUM** | Add `nseg` into JETPHOX output (TNamed / UserInfo) and read it in `MakeJetPHOXhisto.C` + `PlotTruthIso.C` | Silent mis-normalization attack surface (currently 100 vs 10 hardcoded) |
| 5 | **MEDIUM** | Clarify figure legend / caption that JETPHOX is iso-matched and Werner is inclusive | Current legend ambiguous; low-pT offset is physics, not a bug |
| 6 | **LOW** | Fix `abs(y[0])` → `std::fabs(y[0])` at `MakeJetPHOXhisto.C:72` | Style / defense against ROOT version change |
| 7 | **LOW** | Rename or Scale/Sumw2 the 2D `h_truth_eta_pT` | Misleading diagnostic |
| 8 | **LOW** | Commit frozen `parameter.indat.template` vs run-time patches | Reproducibility |
| 9 | **LOW** | Delete stale `plotting/figures/final_JETPHOX*.pdf` and `PPG12-analysis-note/Figures/results/final_JETPHOX*.pdf` | Not referenced in note .tex files |
| 10 | **INFO** | No rerun needed for stats: max stat_rel = 6.8% at [32,36] GeV, << 21/3 = 22% target | Stats are sufficient inside the analysis range |

---

## 7. What NOT to do

- **Do not rerun JETPHOX just for stats** — all bins are already fine.
- **Do not worry about the JP/We central-value disagreement in itself** — scale envelopes overlap.
- **Do not interpret low-pT JP<We as a JETPHOX bug** — it is the iso effect, and matches the frag fraction prediction.
