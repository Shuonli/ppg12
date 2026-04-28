# PPG12 Analysis Note — Publication-Readiness Audit

**Scope.** Five-axis review of `PPG12-analysis-note/` (entry point `main.tex` plus all `\input` chapters and the consolidated DI / vertex-rw appendices). The four review axes requested were (1) numerical self-consistency, (2) prose vs YAML configs, (3) prose/captions vs the actual figures, (4) publication tone (no jargon, no dev-process artefacts). A fifth axis covers physics correctness of formulas and methodology.

**Verdict: NEEDS REVISION before circulation.** The note has 11 blocker-level items (wrong cut values quoted in prose, total-systematic disagreement across chapters, stale narrative, broken cross-reference, plot/caption mismatches), ~30 major items (publication tone, code-name exposure, internal-narrative passages, captions describing the wrong content), and a long tail of minor polish.

---

## 1. Reviewer pass matrix

| Audit axis | Agent | Verdict | Notable findings |
|---|---|---|---|
| 1a Numerical self-consistency | `note-critic` | NEEDS-REVISION | Stale iso/BDT formulas in 4 places; total-syst disagreement (29.7 vs 31.97) |
| 1b Prose vs YAML configs | `data-explorer` | NEEDS-REVISION | Iso & BDT formula values stale; jet-sample table uses DI binning; 9 vs 25 features |
| 1c Prose vs plots | `plot-cosmetics-reviewer` | NEEDS-REVISION | ROC caption claims wrong nominal; CT14lo↔CT14nlo flip; iso-ET 6-panel mislabel; 16.6 pb⁻¹ baked into tower-acceptance plots |
| 1d Publication tone | `note-critic` | NEEDS-REVISION | Code names, "all-z", "PPG12", "_nom"/"_double", lab-notebook section titles |
| 1e Physics correctness | `physics-reviewer` | NEEDS-REVISION | Cut-value text errors (CRITICAL); cross-section formula doesn't factorize ε; ABCD R-factor not stated |

All five reports converge on the same top items, increasing confidence. No fixes have been applied; this report is findings-only.

---

## 2. Blocker findings (must fix before publication)

### B1. Stale isolation cut formula in three places

| Location | Note text | Config truth | Source |
|---|---|---|---|
| `reconstruction.tex:82` | `isoET < 0.502095 + 0.0433036·E_T_reco` | `0.490 + 0.037·E_T` | `config_bdt_nom.yaml:33-34` |
| `reconstruction.tex:86` | non-iso `> 0.8 + 0.502095 + 0.0433036·E_T_reco` | non-iso `> 0.8 + 0.490 + 0.037·E_T` | same |
| `systematics.tex:81` | identical stale formula | identical fix | same |

The plot `Figures/analysis/ETCut_FitResults.pdf` legend correctly shows `0.490 + 0.0370·pT` — i.e. the **plot is right and the prose is stale**. Fix direction: update text to `0.490 + 0.037·E_T`. The 0.502095/0.0433036 numbers appear to come from a previous iso fit and were never refreshed.

### B2. Wrong photon-ID BDT thresholds in Table 4

`reconstruction.tex:168, 173` (Table `tab:photonid`):

| Edge | Note quotes | Config truth |
|---|---|---|
| Tight `BDT >` | `0.8667 − 0.00667·E_T` | `0.8333 − 0.00333·E_T` |
| Non-tight upper | `< 0.8667 − 0.00667·E_T` | `< 0.6667 + 0.00333·E_T` (positive slope) |
| Non-tight lower | `> 0.9 − 0.02·E_T` | `> 0.7333 − 0.01333·E_T` |

These are not rounding errors. The note describes a different parametrization that crosses (E_T=10, BDT=0.8) but diverges from the production cut at higher E_T (e.g. at 30 GeV the note says `>0.6667`, the production says `>0.7333`). The non-tight upper edge slope sign is also wrong in the note (note: negative; production: positive). A reader reproducing the analysis from the note would apply different selections. **Source of truth**: `config_bdt_nom.yaml:173-174, 229-234`.

### B3. Three different total-systematic numbers across chapters

| Location | Quoted total at highest E_T |
|---|---|
| `systematics.tex:186` (prose) | `−29.7% / +40.6%` |
| `systematics.tex:209` (Table 3) | `29.70 / 40.59` |
| `conclusion.tex:5` (and abstract) | `−31.97% / +41.18%` |

The prose at `systematics.tex:186` itself states that the row-wise quadrature in Table 3 is *not* the per-bin total ("the row-wise quadrature of the per-source maxima in the table does not equal the per-bin total"), then quotes 29.7/40.6 as the per-bin total, which are exactly the row-wise quadrature numbers. The conclusion's 31.97/41.18 is yet a third value. Pick one set, propagate to abstract / conclusion / systematics intro / Table 3 caption, and explain explicitly in the systematics intro how the per-source-max table relates (or doesn't) to the per-bin total.

### B4. `noniso04` value misquoted

`systematics.tex:83` says `noniso04: isoET > reco_iso_max + 0.4 GeV`. The actual config `config_bdt_noniso04.yaml:45` has `reco_noniso_min_shift: 0.1` GeV. The variant name "04" refers to a historical 0.4 GeV shift; the production-time value is 0.1 GeV. The "looser" `noniso10` (=1.0 GeV) is correctly stated. Fix the quoted shift, and either drop the variant name or rename it to track the actual value.

### B5. Selection table 2 lists DI MC samples, not main-pipeline samples

`selection.tex:51-66` Table `table:MCdataset` lists `jet8 / jet12 / jet20 / jet30` and `photon5 / photon10 / photon20`. The main analysis pipeline uses `jet10 / jet15 / jet20 / jet30 / jet50` (per `CrossSectionWeights.h` and `MergeSim.C`); `jet8 / jet12 / jet20 / jet30 / jet40` is the DI study binning. The current Table 2 is a hybrid that matches neither set exactly. Replace the jet rows with the main-pipeline binning; the photon rows are correct.

### B6. "1.5 mrad nominal published period" stale narrative in DI appendix

`double_interaction_consolidated.tex:21`: `\textbf{1.5 mrad} (nominal published period)`. The current published nominal is the **all-range, all-z combined** sample with L = 64.37 pb⁻¹ (selection.tex:119, abstract, conclusion). The 1.5 mrad period alone is no longer the nominal. Drop "nominal published period" from this sentence; treat both periods symmetrically as inputs to the merged measurement.

### B7. Tower-acceptance appendix presents `mask_phisymm_tight` as a cross-check, but it is the production nominal

The nominal config `config_bdt_nom.yaml:106-108` switched to `tower_mask_on: 1`, `tower_mask_name: "mask_phisymm_tight"` on 2026-04-23 — the tight phi-symmetric mask is now applied symmetrically to data and MC at production time. The appendix (`appendix.tex:140-294`) still describes all four masks (preselect, common, tight, OR) as cross-check shifts vs an *unmasked* baseline `σ_nom = 1963 pb`. The integrated `1.37%` shift labelled `mask_phisymm_tight` is now **absorbed into the nominal**, not a delta to it. Rewrite the appendix so the "nominal" reference is the masked-tight result; the other three masks stay as cross-checks of the choice.

### B8. Broken cross-reference: `\ref{sec:sys:tower_acc}`

`appendix.tex:284` cites `Section~\ref{sec:sys:tower_acc}`. No such label exists anywhere in the note (the section is commented out at `systematics.tex:128-132`). The LaTeX compilation log already prints `Reference 'sec:sys:tower_acc' on page 65 undefined`. Either resurrect the section in `systematics.tex` or rewrite the appendix sentence to point to the appendix-internal context.

### B9. DI section narrative contradicts Table 4

`double_interaction_consolidated.tex:282`: "Mixing in 22.4% of DI MC closes the data–MC gap substantially for both at E_T≥14 GeV."

Table 4 (`di:tab:chi2_full`, lines ~380-386) shows for the 1.5 mrad period (which uses 7.9%, not 22.4%):

| variable, bin | nominal χ²/ndf | mixed-MC χ²/ndf |
|---|---|---|
| weta_cog 10–14 GeV | 9.72 | **42.31** (worse) |
| BDT score 10–14 GeV | 5.98 | **33.08** (worse) |

So at low E_T the mix actually worsens agreement at 1.5 mrad. Either the narrative is overgeneralizing from the 0 mrad row (where mixing helps), or the table values are wrong. Reconcile.

### B10. CT14lo vs CT14nlo flipped between caption and plot

`results.tex:26` (caption of Fig `final_bdt_nom`): "The JETPHOX NLO curve shown here uses the **CT14lo** PDF set". The plot legend on `Figures/results/final_bdt_nom.pdf` reads "**CT14nlo PDF / BFG II FF**". One of the two is the displayed truth. The appendix Fig `compare_3pdf.pdf` correctly uses CT14lo as the LO baseline, so the natural reading is that `final_bdt_nom.pdf` is the CT14**nlo** comparison and the caption needs to be corrected — but verify against the JETPHOX run that produced the file.

### B11. Tower-acceptance maps have stale 16.6 pb⁻¹ baked into the figures

`Figures/tower_acceptance/official_tower_*_mc.pdf` and `..._data.pdf` have an in-figure title that reads "16.6 pb⁻¹ (1.5 mrad), all-range merge". The note's nominal lumi is 64.37 pb⁻¹ (all-z fiducial). These plots predate the all-z migration; they need to be regenerated from the current pipeline output before they go in a publication-track note. The captions at `appendix.tex:218-232` don't quote a number, so the contradiction is purely between what the figure *shows* and the rest of the note — but it is visible to any reader.

---

## 3. Major findings (should fix)

### 3.1 Numerical / config consistency

- **9 vs 25 BDT features**: physics-reviewer noted that `reconstruction.tex:306` describes the **NPB BDT** with 25 features and `reconstruction.tex:356` describes the **photon-ID BDT** with 9 features. The two are different classifiers; both descriptions can be self-consistent. Verification needed: (a) the photon-ID BDT (`base_v3E` model in `bdt_et_bin_models`) does in fact use 9 features and not 25 at inference time; (b) `FunWithxgboost/config.yaml` `data.features` (25 entries) is the NPB or full-feature training input. If both are true the note is correct; if `apply_BDT.C` reads 25 features for the photon-ID BDT, the prose Table 5 is missing 16 features.
- **50/50 train/test split (prose) vs 70/10/20 (config)**: `reconstruction.tex:356` says "MC events are split 50/50". `FunWithxgboost/config.yaml:160-162` defines `train_size 0.7 / val 0.1 / test 0.2`. Reconcile.
- **Reported binning omits 32-36 GeV bin**: `analysis.tex:220` lists fiducial binning `[10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32]` — 10 bins. Config has 11 reco bins ending at 36 GeV. The reported range is 12–32 GeV, so the 32-36 bin is **computed but unreported**. Add a single sentence explaining this rather than silently dropping the bin.
- **Truth-pT prior reweight description**: `analysis.tex:222` says "fit range [8, 35] GeV" with "[10, 30] GeV clamp"; config has `xmin: 0`, `xmax: 100`, `clamp_pT_min: 10`, `clamp_pT_max: 30`. The clamp matches. The "8–35" is the range used to derive the rational (Padé) form, and is consistent with the comment at `config_bdt_nom.yaml:277-280`. Suggest stating "fit derived on E_T_truth ∈ [8, 35] GeV; clamped to [10, 30] GeV at evaluation" so the reader can map the prose to the YAML.
- **Run range**: `selection.tex:11` Table 1 quotes `47289–53880` while `config_bdt_nom.yaml:262-263` has `run_min: 47289 / run_max: 54000`. The 53880 may be the highest run with data; the wider config bound is harmless. Add a sentence to clarify rather than leave the discrepancy.

### 3.2 Plot/caption mismatches

- **ROC caption misnames the nominal isolation**: `reconstruction.tex:52` caption claims "the nominal 70 MeV tower-based definition (R = 0.3) is shown in blue"; the analysis nominal is **topo-cluster R = 0.4** (config `use_topo_iso: 2`). Body text on line 66 already states this correctly. Fix the caption.
- **Iso-ET 6-panel grid mislabel**: `analysis.tex:80-82` caption says "signal MC with tight ID (blue)" while the actual plot legend labels the blue series "PYTHIA Inclusive Jet (tight)". The caption mislabels inclusive jet MC as signal MC; it also omits the orange "PYTHIA Inclusive Jet (non-tight)" curve and the two exponential fits visible in the plot. Rewrite the caption to match what is shown.
- **S/B caption number**: `reconstruction.tex:103` says "S/B is about 15% at 10 GeV and rises to about 50% near 30 GeV"; the plot `analysis/SB.pdf` shows ~18% at 10 GeV and ~65–70% at 30 GeV. Update the rounding.
- **Energy-resolution systematic**: caption / Table 3 quote +21.24% high-side max for energy resolution; the plot `systematic/syst_rel_eres.pdf` last bin reaches roughly +30%. Either the table value is from a non-final bin, or the plot was regenerated and the number wasn't refreshed. Re-extract numerically and reconcile.
- **Debug labels visible in figures**: several plots have in-figure annotations like `bdt_nom`, `tune: bdt_nom`. Cosmetic but should not appear in a publication-track note. Plots flagged: `et_sbs.pdf`, `purity_bdt_nom.pdf`, `mc_purity_correction_bdt_nom.pdf`.
- **Python matplotlib renders mixed with ROOT plots**: vertex-rw figures, run-by-run iso-ET, and JETPHOX PDF compare are matplotlib renders that break the visual style of the rest of the note. Acceptable for an internal note; should be rerendered through the sPHENIX style helpers if these stay in the published version.
- **pT-range labels in shower-shape captions**: `pt0`/`pt1`/`pt2`/`pt3`/`pt4` in shower-shape figure macros expand to E_T ranges that don't consistently match the captions. E.g. `appendix.tex:19` caption says `12–14 GeV` for `pt1` while `reconstruction.tex:197` caption says `10–14 GeV` for `pt0` (with the same `pt`-to-GeV scheme). Re-derive the mapping from `plotcommon.h:ptRanges`, harmonize all captions, and remove the cut0/cut1/cut2 internal tags from prose entirely.
- **`26–35 GeV` figure-caption labels**: at least 4 captions (reconstruction.tex:215, :251, appendix.tex:37, :77, :135) write "26–35 GeV" for the high-pT showershape panels. The pT axis upper edge in the analysis is **36 GeV**, not 35. These should read "26–36 GeV" (or, more precisely, the GeV range corresponding to the actual figure x-range — possibly `pt2`/`pt3`/`pt4` map to bins higher than 26 GeV).

### 3.3 Physics-correctness items

- **Cross-section efficiency not factorized**: `results.tex:13-19` defines $d\sigma/dE_T^\gamma d\eta = (1/\mathcal{L}) \cdot Y^{rec} / (\mathcal{E} \cdot \Delta E_T \cdot \Delta\eta)$. The workflow lists $\varepsilon_{reco} \times \varepsilon_{id} \times \varepsilon_{iso} \times \varepsilon_{MBD}$ separately at `analysis.tex:185-191`. Spell out the factorization explicitly in the equation so the reader can map prose to formula.
- **ABCD R-factor not stated**: `analysis.tex:97-103` writes the leakage-free background as `N^B · N^C / N^D`. The user-spec form is `R · N^B · N^C / N^D` with R = 1 only when the tight/non-tight axis and iso/non-iso axis are independent in the *background*. The MC-purity closure correction `C_MC` (analysis.tex:172) absorbs the residual non-independence. State explicitly that the R = 1 assumption is being made and that `C_MC` corrects the residual deviation.
- **MC-purity closure narrative ("few-percent") inconsistent with 10.5% systematic**: `analysis.tex:68` says "the ABCD-extracted purity differs from the MC truth purity at the few-percent level"; `systematics.tex:108` quotes a max systematic of ±10.5%. Rephrase as "up to ~10%, treated as the dominant purity systematic (Section 4.x)".
- **L1 plateau systematic +0.40% flat may understate +1.5% bin-wise correction**: the L1 trigger correction itself reaches +1.5% at the lowest analysis bin (`analysis.tex:7`); treating the systematic as a flat 0.40% post-unfold underrepresents the bin-by-bin sensitivity at low E_T. Either propagate the fit covariance band per bin, or argue explicitly why a flat 0.40% bounds the bin-wise variation.
- **Truth-vertex reweight 1.5 mrad non-closure not propagated**: `vertex_rw.tex:111` shows `χ²/ndf = 13.7` (ndf 33) and `max|r-1| = 0.111` (i.e. 11% closure floor at the closure-window edge) for 1.5 mrad. There is no explicit cross-section systematic for this non-closure in `systematics.tex`. Confirm whether it is absorbed into the DI-fraction or response-reweight systematics, or add a separate envelope.
- **DI fraction $f_{best} = 0.000$ at 1.5 mrad**: `systematics.tex:173` reports the data $\chi^2$-best DI fraction at 1.5 mrad pulls to zero. This is a non-trivial physics statement (data prefers no DI MC at 1.5 mrad). Add a short sentence explaining whether this is interpreted as a real fluctuation (low DI rate close to zero is plausible at 1.5 mrad), a methodological limitation (BDT-axis insensitivity at low DI fraction), or a residual data-MC tension unrelated to DI. The DI consolidated appendix already discusses this in the caveat block — link the systematics text to that discussion.

### 3.4 Publication tone

The note carries substantial dev-process and pipeline-internal artifacts. Worst offenders, with suggested replacements:

| Where | Offending text | Replace with |
|---|---|---|
| `selection.tex:93` | `applied at fill time as a per-cluster weight 1/ε_L1 in \texttt{RecoEffCalculator\_TTreeReader.C}` | "applied as a per-cluster weight 1/ε_L1 at the analysis efficiency-correction stage" |
| `selection.tex:110` | `Cluster threshold is 5 GeV (slimtree minimum)` | "Cluster threshold is 5 GeV (analysis preselection minimum)" |
| `selection.tex:119` | `the beam-delivered (all-z) fiducial convention` | "the beam-delivered fiducial convention (no truth-vertex restriction in the efficiency denominator)" |
| `analysis.tex:5` | `Cluster EMCal towers using RawClusterBuilder without sub-cluster splitting` | "Cluster EMCal towers using the standard cluster builder with sub-cluster splitting disabled" |
| `analysis.tex:15` | `combined at histogram fill time / per-period MC histograms ... summed across periods to produce the all-range MC` | "the per-period MC histograms are luminosity-weighted and summed to produce the combined MC" |
| `analysis.tex:190` | `truth-vertex window set by \texttt{vertex\_cut\_truth} (the all-z fiducial in the nominal configuration)` | "an open truth-vertex window so the efficiency converts the beam-delivered luminosity into the analysis fiducial sample" |
| `analysis.tex:222, 229` | `nominal-baseline pipeline run / stored under analysis.unfold.truth\_reweight in the nominal config` | "preliminary unfold of the data spectrum / the rational fit applied photon-by-photon at response-matrix build time" |
| `systematics.tex:1, 13, 169, 282 etc.` | `response-driven (Resp) variations / Resp group / Resp envelope` | "response-driven variations" (drop the abbreviation; expand once and reuse) |
| `systematics.tex:80-84` | `\mathrm{reco\_iso\_max} / (noniso04) / (noniso10)` | drop config-field tag and variant labels; quote only the gap value |
| `appendix.tex:174-180, 260-263, 405-410` | `\texttt{preselect}, \texttt{common}, \texttt{tight}, \texttt{OR}, mask\_phisymm\_tight`, `\texttt{b2bjet\_pt5}`, `\texttt{b2bjet\_pt5\_npb0}` | "preselect mask", "tight mask", "OR mask"; "b2b $p_T^{jet} \ge 5$ GeV (NPB on)" etc. |
| `double_interaction_consolidated.tex:27, 56, 167, 411-417, 583-588` | `\texttt{photon10\_double}, \texttt{jet12\_double}, \*\_nom, \*\_double, mix\_weight, hadded, signal\_combined, jet\_inclusive\_combined, mixed\_0mrad\_vtxrw, compare\_energy\_response.C` | sample-name and filename mentions removed; describe by physics ("dedicated single-interaction signal-photon MC", "blended sample with the truth-vertex reweight applied", etc.) |
| `double_interaction_consolidated.tex:200, 207` | `the 20 pp gap between single (93.5%) and double (73.1%)` | "the 20 percentage-point gap" (current text overloads `pp` as collision system AND percentage points — visually identical to `\pp`) |
| `vertex_rw.tex:122, 301` | section/paragraph titles `Hero figure / Why 0 mrad is easy / Why 1.5 mrad is hard` | "Summary figure / Convergence at 0 mrad / Convergence at 1.5 mrad" |
| `vertex_rw.tex:64, 209-217, 244` | "we are trying to match", "stuck edge bins", "silently rejected upstream" | neutral phrasing |
| `vertex_rw.tex:337-340` | "extension to the full sample blend (including jet12) is anticipated to shift the closure metrics by O(1%)" | either drop or report the actual measured cross-check |
| `double_interaction_consolidated.tex:27` | `The PPG12 DI programme uses ...` | "This analysis uses dedicated double-interaction MC ..." (`PPG12` is internal collaboration tag) |
| `selection.tex:93` | reference to `(Section~\ref{sec:analysis} step~2)` | cross-reference by section name, not by step number |
| `introduction.tex:12` | `future PPG sPHENIX heavy-ion programs` | "future sPHENIX heavy-ion programs" (`PPG` is internal acronym) |
| `results.tex:26` | `NLO pQCD calculations by Werner Vogelsang` | cite the calculation (paper / private communication note); do not name an individual without a citation |
| `appendix.tex:105` | `common cuts (cut1: prob, e11/e33, NPB score > 0.5, weta_cogx bounds)` | "the pre-selection (e11/e33, NPB score, w_η bounds)" — drop `cut1` and the `_cogx` underscore |

The shared theme is to remove every config-field name, code-class name, file path, branch-suffix, `_nom`/`_double`/`_0rad`/`_1p5mrad` tag, condor-script reference, internal-narrative paragraph header, and "as we ran the pipeline" commentary. None of those carry physics meaning to an external reader; all of them break publication tone.

### 3.5 Notation / style

- `\GeV` macro used in some places (introduction, selection) and bare `GeV` in others (results.tex:21). Settle on `\GeV` everywhere.
- `\pb` macro vs `pb$^{-1}$` literal: defs.sty defines `\pb`; use it consistently.
- `non-physical background` (reconstruction.tex:197) vs `non-physics background` (everywhere else): pick `non-physics`.
- Section / subsection title capitalization uneven; standardize.
- Three-place decimal precision for percentages (e.g. `+9.13%` vs `9.1%`) — fine, but match across abstract / conclusion / systematics caption.
- Several `roughly`, `essentially`, `kind of` hedges — drop or replace with quantitative phrasing.
- `analysis.tex:3` "summarized as following:" → "summarized as follows:".
- `double_interaction_consolidated.tex:56-58`: missing closing parenthesis at end of sentence (LaTeX still compiles; reader notices).

---

## 4. Recommended fix order

1. **Numerical/config corrections (B1, B2, B3, B4, B5)** — these are direct factual errors that an external reviewer would call out immediately. Update the prose to match `config_bdt_nom.yaml`; if any prose value is intentionally different from the production cut, add a sentence explaining.
2. **B6/B7/B8/B11** — restructure the DI consolidated appendix and tower-acceptance appendix so the all-range allz nominal is the framing throughout, and replace the stale tower-acceptance figures with current-pipeline outputs.
3. **B9** — reconcile the DI mixed-MC narrative with Table 4 numerics (either edit narrative to acknowledge the 1.5 mrad low-E_T degradation, or re-extract the table with the corrected DI flow).
4. **B10** — verify which PDF set produced `final_bdt_nom.pdf` and either fix the caption or replot.
5. **Plot/caption pass (3.2)** — rerun captions against actual figures; remove debug labels; rerender Python plots in sPHENIX style.
6. **Physics-correctness wording (3.3)** — factorize the cross-section formula, state R=1 assumption explicitly, soften "few-percent" claim, decide on L1-plateau systematic propagation.
7. **Tone pass (3.4)** — global find-and-replace on code names, file names, jargon. The replacement table above covers the bulk.
8. **Notation / style polish (3.5)** — last.

After items 1–4 the note is technically defensible. After 5–6 it reads as a measurement document. After 7–8 it is a publication-track note.

---

## 5. Files referenced

Prose audited:
- `PPG12-analysis-note/{main, introduction, selection, reconstruction, analysis, systematics, results, conclusion, appendix, vertex_rw, double_interaction_consolidated}.tex`

Configs cross-checked:
- `efficiencytool/config_bdt_nom.yaml`
- `efficiencytool/config_bdt_nom_{0rad, 1p5mrad}.yaml`
- `efficiencytool/config_bdt_noniso{04, 10}.yaml`
- `FunWithxgboost/config.yaml`

Code cross-checked:
- `efficiencytool/RecoEffCalculator_TTreeReader.C`
- `efficiencytool/make_bdt_variations.py`
- `efficiencytool/CrossSectionWeights.h`

Figures inspected (sample of ~25 of the ~70 cited):
- `Figures/results/{final_bdt_nom, final_phenix_bdt_nom, final_all_bdt_nom}.pdf`
- `Figures/analysis/{ETCut_FitResults, SB, iso_ET_pt0, purity_bdt_nom, mc_purity_correction_bdt_nom, et_sbs, roc/roc_combined_pt0_10_15}.pdf`
- `Figures/tower_acceptance/{official_tower_preselect_mc, official_tower_preselect_data, official_tower_common_mc, official_tower_common_data}.pdf`
- `Figures/systematic/{syst_rel_eres, syst_sum_rel}.pdf`
- `Figures/Run_by_run/mean_isoET_vs_run.pdf`
- `Figures/double_interaction/vertex_reweight/vertex_reweight_iter_two_period.pdf`
- `Figures/etc/{h_maxEnergyClus_NewTriggerFilling_doNotScale_PhotonTurnOn (1), combine}.pdf`
- `Figures/jetphox/compare_3pdf.pdf`

LaTeX log warning confirming B8: `main.log:3413 — Reference 'sec:sys:tower_acc' on page 65 undefined`.
