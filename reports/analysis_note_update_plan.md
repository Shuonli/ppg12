# PPG12 Analysis Note — Update Plan

**Date:** 2026-04-23
**Scope:** Section-by-section punch list of plot regenerations and text edits needed to bring `PPG12-analysis-note/main.tex + sections/*.tex` up to date with the current pipeline as documented in `reports/` and codebase `wiki/`.
**Not in scope:** Executing any of these edits. Reviewing the DI appendix itself (already pulled in as a pre-built `\includepdf` on 2026-04-18 — its internal numbers were already audited in `reports/double_interaction_note`).

---

## 0. How this plan was built

- **Note baseline:** Main body last substantively edited on submodule commit **2026-03-10**; only the DI appendix (`\includepdf`) was added on 2026-04-18. Every non-DI-appendix section therefore reflects the 2026-03-10 state.
- **Comparison corpus:**
  - 25+ deduped reports under `reports/` (2026-04-08 → 2026-04-23), source `.md`/`.tex` read.
  - Codebase wiki under `wiki/pipeline/`, `wiki/concepts/`, `wiki/reference/`, `wiki/guides/` (all 24 articles refreshed 2026-04-20).
  - Rule files under `.claude/rules/` (yaml-config, double-interaction, root-macros, python-analysis).
  - Git log (2026-02-23 → 2026-04-23) for pipeline-path commits.
- **Number-level audit:** 10 specific cited numbers in the note were re-derived against live configs (`config_bdt_nom.yaml`), `wiki/reference/constants-sync.md`, memory notes, and commit diffs (agent `numerical-rederivation` — 7 FAIL, 2 PASS, 1 AMBIGUOUS).

Within each note section below, three explicit buckets:
1. **Stale figures** — figures that have been regenerated or superseded, or that are referenced by a path that does not exist on disk.
2. **Stale text** — outdated numbers, retired methods, or descriptions inconsistent with the current pipeline.
3. **Missing content** — subsections/figures that should exist but don't.

Priority: **HIGH** = physics-wrong or retired-method-as-current, **MED** = stale figure or outdated number that doesn't invalidate the physics, **LOW** = polish/cosmetic.

---

## 1. TL;DR — one-screen summary

The single most consequential gap: **the analysis note still describes a 1.5-mrad-only, legacy reco-vertex DI analysis with 16.6 pb⁻¹ luminosity and trigger bit 26. The production pipeline has moved to an all-run-range all-z nominal (64.3718 pb⁻¹) with lumi-weighted period merging, single-pass truth-vertex DI reweighting, bit-30 (Photon_4_GeV) L1 trigger, and a post-tower-mask integrated cross-section** (pre-mask σ_nom = 1963 pb; post-mask σ_nom = 1990.57 pb per `reports/tower_acceptance_2026-04-22.tex` L77 and L474 — both stated explicitly as all-z fiducial). The all-z nominal uses `vertex_cut_truth: 9999.0` to widen the truth-vertex denominator of the MBD-vertex efficiency, while `vertex_cut: 60` is retained for data and the MBD-efficiency numerator. The 60-cm fiducial setup (48.9309 pb⁻¹) is now a **cross-check only**. Sections 2, 4, 5, 6 all need substantive edits; the DI appendix (as `\includepdf`) is current but the methodology is nowhere in the main-body sections.

Secondary big-picture items:
- **Isolation** formula in the note (`1.08128 + 0.0299107·ET`) is the legacy FunWithxgboost training value; nominal analysis now uses `0.502095 + 0.0433036·ET`. The |vz|-dependent cone mis-centering bug (iso_cone_fix report, 2026-04-18) is not mentioned.
- **Systematics section** is missing: tower acceptance (+2.5% one-sided envelope), L1 trigger plateau low-bias (~0.4%), truth-vertex-reweight closure (~1–2%), donut isolation variants, NT-BDT non-closure.
- **Results section** cites values computed with pre-tower-mask, pre-full-range-merge conventions — needs a single full regeneration pass. σ_nom cascade affects §5 (total syst envelope), §6 (final number and figures), and §7 (conclusion).
- **`conclusion.tex` is empty** (0 lines, not 1) — must be written from scratch after §6 is updated.
- **`vertex_rw.tex` is orphaned** — 502-line file not `\input`-ed by `main.tex` or `appendix.tex`; it never compiles into the PDF. Lives in the repo as stale source; action is either delete or promote into the appendix.
- **Recommended starting point:** `reports/post_preliminary_updates.tex` is an explicit consolidated-methodology-delta report (BDT refresh, NPB, DI studies, shower shape comparisons) covering exactly the territory this plan maps. A postdoc should open that report alongside this plan.

Estimated effort: one postdoc × 2–3 weeks for the text edits + one coordinated pipeline rerun (all-range nominal, plots regenerated, updated JETPHOX curve when ready).

---

## 2. Cross-cutting themes

| # | Theme | Affected sections | Source reports / wiki |
|---|-------|-------------------|-----------------------|
| T1 | **Full-run-range all-z merge (Plan B)** is now the nominal — note still describes per-period (1.5 mrad) 60-cm-only analysis. Primary lumi is 64.3718 pb⁻¹ (all-z, `vertex_cut_truth: 9999.0`); 60-cm (48.9309 pb⁻¹) is a cross-check only | §2 (lumi), §4 (workflow, fiducial def), §5 (lumi syst), §6 (results) | `wiki/pipeline/full-run-range.md`, `.claude/rules/yaml-config.md`, commit `402104c`, `config_bdt_nom.yaml:261` (`lumi: 64.3718`) |
| T2 | **Single-pass truth-vertex reweight** replaces legacy two-pass reco-vertex reweight everywhere | §3 (reco), §4 (efficiency), §A (vertex_rw) | `wiki/concepts/double-interaction-efficiency.md`, `.claude/rules/double-interaction.md`, `reports/truth_vertex_closure_nominal.tex`, `reports/vtxreweight_impact.tex` |
| T3 | **L1 trigger** is bit 30 (Photon_4_GeV + MBD NS), not bit 26; plateau efficiency 99.58% not 99.1%; ~0.4% low bias if uncorrected | §2 (trigger), §5 (syst) | `reports/trigger_bit30_turnon.md`, memory `project_trigger_bit30_turnon.md`, `project_cut_efficiency_audit.md` |
| T4 | **Tower acceptance / dead-tower masks** — new systematic channel (+2.5% one-sided); dead-tower MC-vs-data asymmetry documented | §3 (selection), §5 (syst) | `reports/tower_acceptance_2026-04-22.tex`, memories `project_tower_acceptance_maps.md`, `project_tower_masking_variants.md` |
| T5 | **Isolation method update** — nominal parametric coefficients updated; |vz|-dependent cone mis-centering fix landed; donut variants (`_excl` family, inner-R scan) available as cross-checks | §3 (iso), §5 (iso syst) | `reports/iso_cone_fix.tex`, `reports/donut_iso_variants.tex`, `reports/isolation_ET_comparison.tex`, `wiki/concepts/isolation-cuts.md` |
| T6 | **DI blending in main-body methodology** — 8 SI/DI MC pairs blended at per-period cluster-weighted fractions (22.4% / 7.9%); deltaR matching replaced by track-ID-only; ET-dependent truth-pT stitching (photon10/20 at 22 GeV, not 30) | §2 (MC), §3 (reco), §4 (efficiency), §5 | `.claude/rules/double-interaction.md`, `reports/double_interaction_note.tex`, `wiki/reference/mc-samples.md`, commit `ce90372` |
| T7 | **BDT infrastructure audit** (LOW-priority reference, NOT a blocking gap) — photon-ID nominal `base_v1E` (9 features) unchanged; dual split/no-split support now exists but isn't physics-impacting for the main result. Mentioned here only so a postdoc doesn't think the BDT section needs a major rewrite. | §3 (photon ID BDT) — no edit needed beyond LOW-priority mention | `reports/bdt_training_audit.tex`, `reports/split_nosplit_dual_bdt.md`, `wiki/pipeline/02-bdt-training.md` |
| T8 | **JETPHOX rework pending** — CT14lo → CT14nlo / NNPDF3.1 planned; CT14lo shift −27% to +11% depending on pT. JETPHOX-vs-Werner drift explained | §1, §6 (theory comparison) | `reports/jetphox_audit.md`, `reports/jetphox_3pdf_comparison.md`, `reports/jetphox_nlo_rerun_setup.md` |
| T9 | **Purity / ABCD** — ET-dependent non-tight BDT bound is current; pure-background R-factor crosses 1.0 at pT≈14 GeV (genuine physics, not artifact); MC purity correction larger at low pT than previously documented | §4 (purity), §5 (purity syst) | `reports/mc_purity_correction_investigation.tex`, `reports/purity_nonclosure_ntbdt.tex`, `wiki/concepts/mc-purity-correction.md`, `wiki/concepts/abcd-method.md` |
| T10 | **Eta-edge migration / delta-R removal** — delta-R<0.1 geometric cut retired from production; eta-edge migration quantified (1.2% SI, 13.7% DI fakes flowing into response matrix) | §3, §4 (unfolding / response) | `wiki/concepts/eta-edge-migration.md`, `wiki/concepts/deltaR-truth-matching.md` |

---

## 3. Section-by-section punch list

### §1 Introduction (`introduction.tex`, ~27 lines)

**Stale figures** (all files exist on disk — they compile correctly; issue is *content*, not missing assets)
- [MED] `Figures/introduction/compton_scattering.png`, `annihilation.png`, `fragmentation.png` — present on disk. Stale status depends on whether the schematic is current; if the Feynman diagrams are standard textbook/external reference images they are likely fine. Postdoc action: inspect once and sign off, not regenerate.
- [MED] `Figures/introduction/pT_Spectrum_IsoET_4GeV_Cut.pdf` (direct/fragmentation fractions) — present on disk. Content may be stale if it was generated with pre-22-GeV stitch boundary and pre-DI-blend MC. Regenerate with the current MC configuration (post `ce90372` stitch fix, with DI blending if desired for consistency with §4) OR leave as is if the note is only illustrating qualitative behaviour.

**Stale text**
- [HIGH] Line 24: `integrated luminosity of 16.6~\pb` → verdict **FAIL**. Replace with the current nominal: **all-range allz 64.3718 pb⁻¹** (live in `config_bdt_nom.yaml:261`; all-range sum of 0 mrad 47.2076 + 1.5 mrad 17.1642). Mention the 60-cm fiducial value (48.9309 pb⁻¹) only as a cross-check in §5. See `.claude/rules/yaml-config.md` lines 26–34 and `project_nominal_allz_setup.md`.
- [LOW] `main.tex:31` abstract currently states "10 < etg < 26 GeV" kinematic range; if §4.5 / §6 resolve to a wider reporting range (e.g., 10–36 GeV post-update), update the abstract in lockstep.

**Missing content**
- [LOW] A one-line statement of scope of the measurement could mention that it is now a combined 0 mrad + 1.5 mrad result — this would make the shift in §4 less jarring.

---

### §2 Data and MC Samples (`selection.tex`, ~137 lines)

**Stale figures** (all present on disk; issue is *stale content*, not missing assets)
- [HIGH] `Figures/etc/10_over_05.pdf`, `20_over_10.pdf`, `combine.pdf` — exist on disk. Content likely stale: these illustrate the MC stitching between photon5/10/20 samples at the pre-commit-`ce90372` boundaries (14 GeV and 30 GeV). Regenerate against the updated stitch boundaries (14 GeV and **22 GeV**). `reports/sample_combining_check.tex` has the infrastructure.
- [HIGH] `Figures/etc/h_maxEnergyClus_NewTriggerFilling_doNotScale_Overlay (1).pdf` and `Photon_4_GeV_ConstantFit.pdf` — exist on disk (the filename-with-space is tolerated by `\includegraphics`; rename for hygiene at some point). Content is from pre-bit-30 measurement. Regenerate against bit-30 selection. `reports/trigger_bit30_turnon.md` has the numbers; a direct regeneration from the `trigger_bit30_turnon` macro is the cleanest path.

**Stale text** (all flagged by `numerical-rederivation` agent Claims #1–4)
- [HIGH] §2.3 Trigger, line 93: `the "Photon 4 GeV + NBD NS $\geq$ 1" trigger (trigger bit 26)` → should be `the "Photon 4 GeV" trigger (trigger bit 30)` (drop the "+ MBD NS ≥ 1" part — Shuonli confirmed 2026-04-23 that the text should just call it "Photon 4 GeV"). Bit 26 is a legacy reference from `makeHist_showerShapes_v6.C`; all nominal configs set `trigger_used: 30`.
- [HIGH] §2.3 line 93, §2.4 line 122: `the trigger efficiency is at a plateau about 99.1\%` → should be **99.58 ± 0.01%** (integrated ET ≥ 8 GeV; per-bin table in `reports/trigger_bit30_turnon.md`). **Number substitution is Phase 1; the logistic-fit figure AND surrounding write-up need regeneration against bit-30 selection — Phase 3a (task #7)**. Shuonli flagged this explicitly 2026-04-23.
- [HIGH] Line 45 figure caption: `photon 10~GeV(20~GeV) sample is fully efficient at 14~GeV(30~GeV) leading photon \etg` → the **20 → 22 GeV** boundary is the one currently in `CrossSectionWeights.h`. Both the caption and the conceptual description of sample coverage need the update.
- [HIGH] §2.2 MC Samples: the note enumerates `photon5/10/20` and `jet5/8/12/20/30`; it does not describe the **8 SI/DI pairs** blended into nominal simulation. Either add a subsection here or defer to §3.5/§4 — but the sample count quoted in the note needs to match what the production pipeline uses.
- [HIGH] §2.1 Data (selection.tex L12–13): the "Run 2024 dataset" table lists `ana509/ana521` rows. **Drop the ana509 row and keep ana521. Run range stays `47289–53880`** (Shuonli confirmed 2026-04-23 — the nominal uses ana521 only).
- [HIGH] §2.4 Event Selection: add a paragraph describing the all-z nominal setup — `|z_reco| < 60 cm` on data + MBD N&S ≥ 1 (numerator), with `vertex_cut_truth: 9999.0` widening the truth-vertex denominator in MC. This is the mechanism that makes the nominal lumi = 64.3718 pb⁻¹ rather than 48.9309 pb⁻¹. Source: `.claude/rules/yaml-config.md` lines 33–34, `project_nominal_allz_setup.md`.

**Missing content**
- [HIGH] Add a paragraph describing the **DI MC samples** (8 SI/DI pairs: photon5/10/20 ± double, jet8/12/20/30/40 ± double). Source: `.claude/rules/double-interaction.md` + `wiki/reference/mc-samples.md`. Mention the cluster-weighted DI fractions (22.4% at 0 mrad, 7.9% at 1.5 mrad, per-period blend).
- [MED] Add a paragraph on sample-combining validation — the `reports/sample_combining_check.tex` closure (machine-precision combining, cross-period additivity) gives confidence in the merge and is a natural bullet here.
- [LOW] Mention the `run28` production tag explicitly (currently no production-tag text beyond `ana509/ana521`).

---

### §3 Photon Reconstruction (`reconstruction.tex`, ~484 lines)

#### §3.1–3.2 Clustering & Energy Response
**Stale figures / text / missing:** PASS on first pass (figures exist; methodology unchanged — main change is the underlying ana521 tag, not the algorithm). LOW-priority polish only.

#### §3.3 Isolated Photons (Topo-cluster R=0.4)
**Stale text**
- [HIGH] Line 77: `\isoET < 1.08128 + 0.0299107 \cdot \etreco` → verdict **FAIL**. Current nominal uses **`0.502095 + 0.0433036·E_T`** (from `config_bdt_nom.yaml` lines 33–34, `reco_iso_max_b`/`reco_iso_max_s`). The `1.08128/0.0299107` numbers are the **legacy FunWithxgboost training-stage** isolation formula and appear only as commented-out lines in the current nominal config.
- [HIGH] Line 81 (non-iso sideband): `\isoET > 1 + 1.08128 + 0.0299107 \cdot \etreco` → should be `\isoET > 0.8 + 0.502095 + 0.0433036 \cdot \etreco`. **Note the shift changes from 1 → 0.8**: live `config_bdt_nom.yaml:45` has `reco_noniso_min_shift: 0.8`, not 1.0 as in the note.

**Missing content**
- [HIGH] Add a subsection on the **|vz|-dependent cone mis-centering fix** (2026-04-18, `reports/iso_cone_fix.tex`). The pre-fix symptom: `cluster_iso_005 ≈ −cluster_ET` for 100% of clusters at |vz|>80 cm (23.7% of Run 24 sample). Post-fix: 0.0% pathological cases. Even if the systematic impact is absorbed, this is a methodology fix that must be documented. Overlap with `reports/post_preliminary_updates.tex` — reuse text from that report.
- [HIGH] Add a paragraph on **available isolation variants** as a foundation for the corresponding systematics: `use_topo_iso = 2` (nominal topo R=0.4), `= 4` (donutFull), `= 5` (donutExcl), plus EMCal-only, inner-ring, `cluster_iso_01`. Source: `reports/donut_iso_variants.tex`, `reports/isolation_ET_comparison.tex`, `wiki/concepts/isolation-cuts.md`. (Elevated from MED to HIGH on critic review — this was dropped from Phase 2 in the prior draft.)

#### §3.4 Signal / Background Identification (Shower Shapes)
**Stale figures / text / missing:** PASS on first pass; figures exist, pre-selection values (`E11/E33 < 0.98`, `0.6 < et1 < 1.0`, `0.8 < E32/E35 < 1.0`, `weta_cogx < 0.6`, NPB > 0.5) match current pipeline.

#### §3.5 NPB BDT (lines 334–410)
**Stale figures / text / missing:** PASS on first pass.
- [LOW] Could note the dual split/no-split cluster node support (`reports/split_nosplit_dual_bdt.md`), but this is more operational than physics.

#### §3.6 Photon Identification BDT (lines 413–484)
**Stale text**
- [MED] The BDT-threshold parametric formula (`tight_bdt_min = 0.8667 − 0.00667·ET`, `non_tight_bdt_max` same, `non_tight_bdt_min = 0.9 − 0.02·ET`) is NOT explicitly stated in the lines reviewed. If the note describes flat thresholds anywhere, update. **Also flag: `constants-sync.md:89` lists the non-tight as `0.7333 − 0.01333·ET`, but live `config_bdt_nom.yaml:225-226` uses `0.9 − 0.02·ET` — pipeline-side audit item**.
- [LOW] Line 421 mentions `base_v1E` with 9 features — matches current (verdict **PASS**). Could mention the 11-variant family as systematic foundation.

**Missing content**
- [LOW] Brief note of the BDT training audit findings (ET/eta reweight flatness validated, `base_vr` dead variant, per-variant filename collision hygiene bug) — optional, and the note-level reader probably doesn't need it.

---

### §4 Analysis Procedure (`analysis.tex`, ~262 lines)

#### §4.1 Workflow Summary
**Missing content**
- [HIGH] The workflow figure and text describe per-period analysis. Add a paragraph on the **Plan B full-run-range merge**: per-period MC pre-scaled via `lumi_weight = lumi/lumi_target`, two feeder configs (`_0rad`, `_1p5mrad`), `merge_periods.sh` hadd across periods to give the nominal all-range MC. Source: `wiki/pipeline/full-run-range.md`, commit `402104c`. The workflow cartoon in Figure 4.1 (if present) likely needs a panel or arrow for this step.

#### §4.2 Truth-level fiducial definition
**Stale text / missing:** PASS on first pass (|η|<0.7, isoET_truth<4 GeV, ΔR=0.3 all match current).

#### §4.3 Purity / ABCD
**Stale text**
- [MED] The ABCD non-tight lower bound is parametric ET-dependent. Make sure the note does not describe a flat non-tight range. Source: `wiki/concepts/abcd-method.md`, `reports/mc_purity_correction_investigation.tex`.

**Missing content**
- [MED] A short paragraph acknowledging that the **R-factor crosses 1.0 at pT≈14 GeV** (pure-background MC, genuine physics, not artifact of non-tight dead zone). This shapes the expected sign of the purity correction and motivates the current MC purity correction prescription. Source: `reports/mc_purity_correction_investigation.tex`, `wiki/concepts/mc-purity-correction.md`.
- [MED] Mention the NT-BDT parametric bound scan as cross-check: mean |ABCD non-closure| = 0.062–0.073 core range; sets the non-closure systematic basis. Source: `reports/purity_nonclosure_ntbdt.tex`.

#### §4.4 Efficiency
**Missing content**
- [HIGH] Add paragraph on the **deltaR-to-track-ID matching switch** (deltaR<0.1 cut retired from production; track-ID-only matching replaces it). The motivation is the DI scenario where vertex displacement forces deltaR to fail geometrically but the track-ID still matches correctly. Source: `wiki/concepts/deltaR-truth-matching.md`, `wiki/concepts/efficiency-deltaR-comparison.md`.
- [MED] Mention the **DI efficiency decomposition** (94%→31% drop breaks into 88% matching artifact + 12% genuine loss; net relative reduction 9.4% on mixed-0 mrad photon-ID chain; cross-section impact +2.8% at 1.5 mrad). Source: `reports/double_interaction_note.tex`.

#### §4.5 Unfolding
**Stale text**
- [HIGH] Line 207: fiducial binning `[10, 12, 14, 16, 18, 20, 22, 24, 26, 35] GeV` → Shuonli confirmed 2026-04-23 the correct reporting binning is `[10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32]` (11 edges, 10 bins, 10–32 GeV range). This adds `28, 32` as edges and drops the 32–35 bin. Also remove/update the `[35, 40]` overflow-bin sentence; if kept, describe it as extending beyond the reporting range for unfolding tails only. The reco/analysis binning in `plotcommon.h` stays at 13 edges `[8, ..., 36]` as the finer internal binning for unfolding; the note's fiducial binning is a reporting aggregation and should say so explicitly.

**Missing content**
- [HIGH] A paragraph on **eta-edge migration** (inward truth→reco migration 1.2% SI, 13.7% DI; flows into response matrix as spurious truth pT; quantified contamination in region C up to 52% at high pT). This affects the unfolding closure interpretation. Source: `wiki/concepts/eta-edge-migration.md`. (Elevated from MED on critic review — was dropped from Phase 2.)

#### §4.6 Luminosity
**Stale text**
- [HIGH] `analysis.tex` L261: `The result is 16.6~pb⁻¹.` → replace with the all-z all-range nominal **64.37 pb⁻¹** (live in `config_bdt_nom.yaml:261` at full precision 64.3718). Per-period allz values: 0 mrad = 47.2076 pb⁻¹, 1.5 mrad = 17.1642 pb⁻¹. The 60-cm fiducial values (32.6574 / 16.2735 / 48.9309 pb⁻¹) appear only as a cross-check in §5. The second 16.6 pb⁻¹ occurrence is in `selection.tex:135` (`$16.6^{+1.4}_{-1.2}~\mathrm{pb}^{-1}$`) — update to `$64.37^{+5.4}_{-4.7}~\mathrm{pb}^{-1}$` (fractional uncertainty carried over from the 60-cm preliminary per Shuonli 2026-04-23).
- [HIGH] Explain in §4.6 text that the lumi is beam-delivered (all-z fiducial) rather than 60-cm-restricted, because the truth-vertex denominator of the MBD-vertex efficiency uses `vertex_cut_truth: 9999.0`. Data still has |z_reco|<60 cm (MBD-eff numerator). The second occurrence of 16.6 pb⁻¹ is in `selection.tex:135` (`$16.6^{+1.4}_{-1.2}~\mathrm{pb}^{-1}$` with asymmetric uncertainties) — update both in lockstep.

---

### §5 Systematic Uncertainties (`systematics.tex`, ~148 lines)

The note lists: energy scale (8–25%), energy resolution (<5%), tight ID (5–10%), non-tight selection (up to 15%), sideband isolation (4–20%), purity fit (10–20%), unfolding (<10%), luminosity (+15.1%/−11.2%), MBD trigger (<10%).

**Stale figures**
- [HIGH] All nine systematic figures (`syst_rel_escale.pdf` … `syst_sum_rel.pdf`) need regeneration against the current all-range nominal. The aggregation script `calc_syst_bdt.py` now operates with `lumi_weight` scaling; output values may shift ~1–5%.

**Stale text**
- [MED] Line 71: upper non-iso sideband variation `\isoET > 2 + ...` and line 73 lower variation `\isoET > 0.5 + ...` — these are parametric offsets on top of the nominal parametric iso cut. Consistent with current pipeline (confirmed in `make_bdt_variations.py`), but verify the iso intercept/slope are the updated ones (see §3.3).
- [LOW] Line 124 lumi syst `+15.1% / −11.2%`: source was the TGC lumi uncertainty model; cross-check against current `lumi/` note text and memory `project_lumi_list_provenance.md` (which flags provenance inconsistency in 60-cm pins vs Bit30Corr sums). Decide whether the asymmetric lumi systematic is still the best estimate.

**Missing content (new systematic channels)**
- [HIGH] **Tower acceptance / dead-tower mask systematic** — recommended envelope **+2.5% one-sided** on nominal cross-section. Four variants (preselect / common / tight phi-symm; OR mask) all give integrated shifts in the +1.4% to +2.3% range. Source: `reports/tower_acceptance_2026-04-22.tex`, memories `project_tower_acceptance_maps.md`, `project_tower_masking_variants.md`. Requires a new figure (comparison of 4 mask variants) and a subsection in §5.
- [HIGH] **L1 trigger plateau low-bias** — ~0.4% integrated low bias if ε_L1 = 1 is assumed (as the current pipeline does). Either add as a one-sided systematic OR apply `1/ε_L1(ET)` correction to data-side yield and drop the systematic. Source: `reports/trigger_bit30_turnon.md` (per-pT table), memory `project_trigger_bit30_turnon.md`.
- [HIGH] **Truth-vertex reweight closure systematic** — post-pipeline closure χ²/ndf = 1.31 (0 mrad), 14.79 (1.5 mrad non-convergence), 9.67 (all-range hadd). Cross-section impact O(1–2%). Source: `reports/truth_vertex_closure_nominal.tex`, memory `project_truth_vertex_reweight.md`.
- [MED] **Isolation cone fix / donut variants** — the iso-cone |vz|-fix (now nominal) removed a pathological class of clusters; the `_excl` and donut variants give a systematic handle on the iso-method choice. Even if not quoted as a systematic (absorbed into iso-sideband variation), should be described here. Source: `reports/iso_cone_fix.tex`, `reports/donut_iso_variants.tex`.
- [MED] **NT-BDT non-closure** — mean |NC| 0.062–0.073 across the core range; flat vs parametric non-tight comparisons. Should appear as a purity-systematic cross-check or be folded into the existing purity-fit variation (currently 10–20%). Source: `reports/purity_nonclosure_ntbdt.tex`.
- [LOW] **Split vs no-split cluster-node BDT** — 22-model family gives a cross-check systematic on clustering choice. Current pipeline ships dual-variant scores; a systematic that swaps cluster node is a natural next step. Source: `reports/split_nosplit_dual_bdt.md`.

---

### §6 Results (`results.tex`, ~48 lines)

**Stale figures**
- [HIGH] `Figures/results/final.pdf` and `final_phenix.pdf` — both need regeneration against:
  - All-range nominal (48.9309 pb⁻¹), not 1.5-mrad-only.
  - Post-tower-mask nominal (σ_nom shift +2.01% to +2.26% vs pre-mask).
  - New systematic envelope (which picks up the tower-acceptance, L1 turn-on, and truth-vertex-closure contributions from §5).
  - Current JETPHOX curve — ideally CT14nlo once production completes (`reports/jetphox_nlo_rerun_setup.md` has the staged pipeline), otherwise hold the CT14lo curve but flag it in the caption.

**Stale text** (all HIGH — σ_nom cascade affects every downstream number in §6/§7)
- [HIGH] Integrated cross-section quoted in the note should be updated to match the new nominal. Two reference numbers in `reports/tower_acceptance_2026-04-22.tex`, both **all-z fiducial**:
  - **Pre-mask σ_nom = 1963 pb** (L77) over ET ∈ [10, 36] GeV — baseline all-range, lumi 64.3718 pb⁻¹.
  - **Post-mask σ_nom = 1990.57 pb** (L474) with tight phi-symm mask applied.
  - Decide which is the primary result and state the tower-mask choice explicitly in the text. See also T4 for systematic treatment.
- [MED] Add a cross-check paragraph reporting the 60-cm fiducial cross-section (uses lumi 48.9309 pb⁻¹; different fiducial acceptance than the allz nominal). This gives a closure on the vertex-cut modelling.
- [HIGH] The PYTHIA / JETPHOX comparison paragraph should be updated:
  - CT14lo → CT14nlo once JETPHOX production lands (pending, per `reports/jetphox_nlo_rerun_setup.md` and `reports/jetphox_3pdf_comparison.md`).
  - Explain the CT14lo → CT14nlo shift (−27% to +11%, monotonic with ET) as a theory-side update in the caption.
  - The JETPHOX-vs-Werner drift (0.68 low-pT → 1.44 high-pT) has a three-driver explanation (iso cuts, Jensen-inequality bin averaging, CT14lo PDF); the note should reference this audit. Source: `reports/jetphox_audit.md`.
  - (Phase 3 only — gated on CT14nlo production completion.)

**Missing content**
- [MED] A short paragraph on the 60-cm fiducial cross-check (lumi 48.9309 pb⁻¹) vs the allz nominal (64.3718 pb⁻¹) — gives a closure-type check on the vertex-cut modelling. Previously listed as "allz cross-check" in an earlier draft; inverted here now that allz is nominal.

---

### §7 Conclusion (`conclusion.tex`, **0 lines — empty file**)

**Missing content**
- [HIGH] `conclusion.tex` is completely empty (not a single line as an earlier draft said). Write from scratch after §6 is finalized; include the integrated cross-section, total systematic uncertainty, and the comparison with JETPHOX / NLO pQCD. Also explicitly state the kinematic range (after §4.5 binning decision). **Sequencing:** this is a Phase 3 item — depends on §6's regenerated final plot and systematic envelope.

---

### Appendix (`appendix.tex`, `vertex_rw.tex`, `double_interaction_standalone.pdf`)

**§A.1–A.3 Shower Shape Variables, IsoET Distributions, Shower-Shape Cut Comparison**
- PASS on first pass. The figures exist; §A matches §3.4 pre-selection values.

**Double Interaction Consolidated (`\includepdf` of `double_interaction_standalone.pdf`)**
- Already updated (2026-04-18 commit `ab6f280`). The internal numbers were audited in `reports/double_interaction_note.tex`. **No action** on the appendix itself.
- However: once §3–§5 of the main body are updated to describe the DI blending inline, decide whether the standalone DI PDF is still needed as an appendix or can be trimmed to a pointer. **Recommend: keep appendix but add a cross-reference from the main body.**

**Vertex Reweight (`vertex_rw.tex`, 502 lines — ORPHANED — DELETE)**
- [HIGH] `vertex_rw.tex` is not `\input`-ed by `main.tex` or `appendix.tex`; it never compiles into the PDF. Truth-vertex reweight methodology is **already in the DI appendix PDF**: `double_interaction_consolidated.tex:25-31,70-75` describes "an iterative data-driven truth-vertex reweight w(z) factorized as w(z_h)·w(z_mb) ..." and that file compiles to `double_interaction_standalone.pdf`, which `appendix.tex:195` pulls in via `\includepdf`. No content migration needed. **Action: `git rm vertex_rw.tex`** (Shuonli confirmed 2026-04-23).

---

## 4. Recommended phasing

### Phase 1 (immediate, ~3–5 days) — Numbers and bits
Low-risk quick wins that fix stale physics statements without requiring any pipeline rerun. Decisions finalized with Shuonli 2026-04-23:
- **Trigger bit 26 → 30; name simplified to "Photon 4 GeV"** (drop the "+ MBD NS ≥ 1" part — selection.tex:93). Plot regen + surrounding logistic-fit write-up deferred to Phase 3a (task #7).
- **Trigger efficiency 99.1% → 99.58%** (selection.tex:93 & :122). Plot regen deferred to Phase 3a.
- **Luminosity values → allz nominal 64.37 pb⁻¹** with fractional uncertainty carried over from 60-cm preliminary: `$64.37^{+5.4}_{-4.7}~\mathrm{pb}^{-1}$` (introduction.tex:24, selection.tex:135, analysis.tex:261). Per-period allz values (0 mrad 47.2076, 1.5 mrad 17.1642) where needed.
- **Iso formula 1.08128/0.0299107 → 0.502095/0.0433036** (reconstruction.tex:77); **non-iso shift 1 → 0.8** (reconstruction.tex:81 — live `config_bdt_nom.yaml:45` `reco_noniso_min_shift: 0.8`, not 1).
- **Photon20 truth-pT boundary 30 → 22 GeV** in caption (selection.tex:45).
- **Data tag ana509 row dropped**; keep `ana521` only. Run range stays `47289–53880` (selection.tex:12–13 table).
- **Fiducial binning (analysis.tex:207): `[10,12,14,16,18,20,22,24,26,28,32]`** — 10 bins, 10–32 GeV range. Update/remove the `[35,40]` overflow sentence accordingly.
- **Abstract kinematic range (main.tex:31) → `10 < \etg < 32~\GeV`** to match new binning.
- **Delete orphaned `vertex_rw.tex`** — the truth-vertex reweight methodology is already in `double_interaction_consolidated.tex:25-31,70-75` which compiles into `double_interaction_standalone.pdf` and is `\includepdf`'d by `appendix.tex:195`. No content migration needed.
- Sanity-check that `cite.bib` has entries for every \cite key referenced in Phase 2 edits — add any missing entries (lumi IAN `ppg09IAN`, JETPHOX NLO refs, DI pileup fraction source).

### Phase 2 (1–2 weeks) — Text additions for methodology gaps
**Revised scope per Shuonli 2026-04-23.** Executed via 5 parallel `report-writer` agents + final `note-critic` sweep. Dropped: iso cone-fix subsection (#4), DI efficiency decomposition (#9). Several items moved to appendix rather than main body.

**Main body (short & focused):**
- #1 [HIGH] DI MC samples paragraph (§2.2) — 8 SI/DI pairs, cluster-weighted fractions 22.4%/7.9%. Source: `.claude/rules/double-interaction.md`, `reports/double_interaction_note.tex`
- #2 [HIGH] All-z nominal event-selection paragraph (§2.4) — `vertex_cut_truth=9999`, data retains 60 cm. Source: `.claude/rules/yaml-config.md`, `project_nominal_allz_setup.md`
- #3 [MED] Sample-combining validation one-line (§2.2). Source: `reports/sample_combining_check.tex`
- #5 [HIGH] **Topo-cluster isoET R=0.4 calculation methodology** (§3.3) — describes how nominal is computed; no variants (per Shuonli: variants don't work and nominal is topo R=0.4). Source: `wiki/concepts/isolation-cuts.md`
- #6 [HIGH] Full-run-range Plan B workflow (§4.1) — `lumi_weight`, `merge_periods`, all-z nominal. Source: `wiki/pipeline/full-run-range.md`
- #7 [HIGH] **Brief** truth-vertex reweight mention (2–3 sentences + reference to DI appendix) (§4). Source: DI appendix `double_interaction_consolidated.tex:25–31,70–75`
- #8 [HIGH] Delta-R retirement / track-ID matching (§4.4). Source: `wiki/concepts/deltaR-truth-matching.md`
- #11 [HIGH] MC purity closure correction as systematic — brief description + derivation + nominal correction figure (§4.3 or §5). Source: `reports/mc_purity_correction_investigation.tex`, `wiki/concepts/mc-purity-correction.md`
- #12a [HIGH] **Brief** tower-acceptance mention in §5 + reference to new appendix section (see #12b)
- #13 [HIGH] L1 trigger plateau systematic (§5) — ~0.4% low bias if ε_L1=1 assumed. Source: `reports/trigger_bit30_turnon.md`
- #14 [HIGH] Truth-vertex-reweight closure systematic (§5) — O(1–2%) cross-section impact. Source: `reports/truth_vertex_closure_nominal.tex`
- #15 ~~NT-BDT non-closure~~ **DROPPED (obsolete per Shuonli 2026-04-23)**
- #16a [MED] **Brief** JETPHOX PDF-comparison mention in §6 caption + reference to new appendix (see #16b)

**Appendix additions (new sections / fill broken references):**
- #10 ~~Eta-Edge Migration new DI appendix subsection~~ **HELD per Shuonli 2026-04-23** (broken forward ref at `double_interaction_consolidated.tex:88` remains unresolved — revisit later). Keep broken ref as a TODO marker.
- #12b [HIGH] **New appendix section on phi-symmetry tower-masking method** in `appendix.tex` — +2.5% envelope derivation, 4 mask variants (preselect/common/tight/OR). Source: `reports/tower_acceptance_2026-04-22.tex`
- #16b [MED] **New appendix section on JETPHOX PDF comparison** (CT14lo vs CT14nlo vs NNPDF3.1) in `appendix.tex`. Source: `reports/jetphox_3pdf_comparison.md`

### Phase 3 (coordinated pipeline rerun + final-plot regen)

**Phase 3a — independent of JETPHOX (~1 week):**
- Regenerate §5 figures (9 systematic figures) with updated syst groupings against the all-range nominal.
- Regenerate §6 final cross-section plot with all-range nominal + tower mask + new systematic envelope (keeping CT14lo JETPHOX curve if needed, with caption flag).
- **Regenerate §2 trigger efficiency figures against bit-30 selection AND update the surrounding logistic-fit write-up** (task #7). `Figures/etc/Photon_4_GeV_ConstantFit.pdf` + maxClusterEnergy overlay/turnon plots. Per-pT values from `reports/trigger_bit30_turnon.md`.
- Regenerate §2 MC-stitching supporting figures (photon ratio plots with 14/22 GeV boundaries).
- Regenerate §6 PHENIX comparison plot.
- **Write conclusion (§7)** once §6 is final.

**Phase 3b — gated on CT14nlo / NNPDF3.1 JETPHOX production landing:**
- Regenerate §6 JETPHOX comparison (primary curve CT14nlo; NNPDF3.1 as cross-check).
- Update §6 caption with CT14lo → CT14nlo shift explanation.
- Consider adding allz cross-check panel to §6.

**Optional / LOW priority:**
- Visually audit §1 Feynman diagram figures — no regeneration likely needed.
- Regenerate §1 direct/fragmentation fraction figure only if current version is pedagogically misleading against current MC blend.

---

## 5. Pipeline-side audit items (flagged along the way — not note edits, but blockers/adjacent fixes)

These are issues surfaced during this audit that are **not analysis-note edits** but will affect the content if not resolved:
1. `wiki/reference/constants-sync.md:89` says non-tight BDT min = `0.7333 − 0.01333·ET`, but live `config_bdt_nom.yaml:225-226` says `0.9 − 0.02·ET`. One of the two is stale. Resolve before documenting the parametric formula in the note.
2. 60-cm lumi pins vs Bit30Corr sums mismatch (memory `project_lumi_list_provenance.md`). If the 48.9309 pb⁻¹ number is to be cited in the note, its provenance needs a definitive audit.
3. JETPHOX production on CT14nlo / NNPDF3.1 is staged but not yet submitted (`reports/jetphox_nlo_rerun_setup.md`). §6 regeneration is blocked on this.
4. `NLO/MakeJetPHOXhisto.C` has uncommitted edits (git status). Resolve before a coordinated rerun.

---

## 6. Deliverables checklist

A tracking checklist for future sessions (one row per concrete edit/regeneration):

### Section 1
- [ ] Verify 3 Feynman diagram figures render correctly in compiled PDF (files exist); sign off, no regen needed unless content is misleading
- [ ] Verify direct/fragmentation fraction figure (`pT_Spectrum_IsoET_4GeV_Cut.pdf`) renders; regen only if it's PYTHIA-only with pre-22-GeV stitch
- [ ] Update lumi citation (introduction.tex:24) → 64.37 pb⁻¹
- [ ] Update abstract kinematic range (main.tex:31) → `10 < \etg < 32~\GeV` (follows new fiducial binning)

### Section 2
- [ ] Regenerate MC ratio stitching figures (`10_over_05.pdf`, `20_over_10.pdf`, `combine.pdf`) with 22 GeV boundary
- [ ] Regenerate trigger efficiency figures (bit 30)
- [ ] Update trigger bit 26 → 30
- [ ] Update trigger eff 99.1 → 99.58%
- [ ] Update photon20 "fully efficient at 30 GeV" → 22 GeV
- [ ] Add DI MC sample paragraph
- [ ] Update data run range to current ana521 / run28 production
- [ ] Add event-selection paragraph explaining the all-z nominal setup (|z_reco|<60 cm on data + `vertex_cut_truth: 9999.0` widening MC truth denominator)

### Section 3
- [ ] Update iso intercept/slope (reconstruction.tex:77) → `0.502095 + 0.0433036·ET`
- [ ] Update non-iso shift (reconstruction.tex:81) → 0.8 (not 1.0)
- [ ] Add iso cone-fix subsection
- [ ] Add iso variants paragraph (donut, `_excl`, inner-R, EMCal-only)
- [ ] (Optional) Note dual split/no-split BDT support

### Section 4
- [ ] Add full-run-range Plan B workflow paragraph
- [ ] Add truth-vertex reweight paragraph (either here or §3)
- [ ] Update analysis.tex:207 fiducial binning → `[10,12,14,16,18,20,22,24,26,28,32]` (10 bins, 10–32 GeV); remove/rephrase `[35,40]` overflow sentence
- [ ] Add deltaR-retirement / track-ID matching paragraph
- [ ] Add DI efficiency decomposition paragraph
- [ ] Add eta-edge migration paragraph
- [ ] Update lumi in §4.6 to 64.37 pb⁻¹ (all-z)

### Section 5
- [ ] Regenerate 9 systematic figures against current nominal
- [ ] Add tower-acceptance systematic section + figure
- [ ] Add L1 trigger plateau systematic (or correction)
- [ ] Add truth-vertex-reweight closure systematic
- [ ] Add NT-BDT non-closure cross-check (or fold into existing)
- [ ] Add iso-variant cross-check text

### Section 6
- [ ] Regenerate final cross-section plot (all-z, all-range, decide mask choice — pre-mask 1963 pb or post-mask 1990.57 pb)
- [ ] State integrated cross-section value explicitly with mask + fiducial definition (all-z nominal)
- [ ] Update JETPHOX curve to CT14nlo (Phase 3b — when production lands)
- [ ] Update PHENIX comparison plot
- [ ] Add 60-cm fiducial **cross-check** paragraph (lumi 48.9309 pb⁻¹) — now the cross-check, not the nominal

### Section 7
- [ ] Write conclusion from scratch (file is empty) with updated integrated cross-section and systematic envelope — Phase 3a after §6 finalized

### Appendix
- [ ] `git rm vertex_rw.tex` — truth-vertex reweight already in the DI appendix PDF via `double_interaction_standalone.pdf` (confirmed 2026-04-23)
- [ ] Add cross-reference from main body §4 to the DI appendix PDF

### Pipeline-side blockers (NOT note edits, but gating factors)
- [ ] Resolve `wiki/reference/constants-sync.md:89` (non-tight BDT `0.7333 − 0.01333·ET`) vs `config_bdt_nom.yaml:225-226` (`0.9 − 0.02·ET`) mismatch — one is stale
- [ ] Audit 60-cm lumi pin provenance (memory `project_lumi_list_provenance.md` flags inconsistency with Bit30Corr sums)
- [ ] Complete JETPHOX CT14nlo/NNPDF3.1 production
- [ ] Commit `NLO/MakeJetPHOXhisto.C` uncommitted changes

---

## 7. Reviewer matrix

| Artifact | 2a physics review | 2b cosmetics | 2c numerical re-derivation | 5 note-critic |
|----------|-------------------|--------------|----------------------------|----------------|
| Wave 1a note inventory | PASS | N/A | N/A | N/A |
| Wave 1b reports inventory | PASS | N/A | N/A | N/A |
| Wave 1c wiki inventory | PASS | N/A | N/A | N/A |
| Wave 1d git timeline | PASS | N/A | N/A | N/A |
| 10 cited-number claims in note | N/A | N/A | PASS (7 FAIL, 2 PASS, 1 AMBIGUOUS) | N/A |
| Plan itself (draft 1) | N/A | N/A | N/A | REVISE (4 substantive fixes) |
| Plan itself (draft 2 — current) | N/A | N/A | N/A | Expected PASS (all 4 REVISE items applied) |

---

## 8. Sources cited

### Reports (deduped set, from `reports/`)
- **`post_preliminary_updates.tex` — USE FIRST.** Consolidated post-preliminary methodology-delta report (BDT refresh, NPB, DI studies, shower shapes). Covers most Phase 2 content; treat as the starting template for new paragraphs.
- `tower_acceptance_2026-04-22.tex` — acceptance systematic, +2.5% envelope; pre-mask σ_nom = 1963 pb (L77), post-mask σ_nom = 1990.57 pb (L474)
- `truth_vertex_closure_nominal.tex` — reweight closure χ²/ndf (1.31 / 14.79 / 9.67 per period)
- `vtxreweight_impact.tex` — ~30–75% xsec shift without reweight (reweight is essential, not systematic)
- `vtxreweight_method_comparison.tex` — jet12 vs photon10 forward-model cross-check
- `sample_combining_check.tex` — combining arithmetic validation + jet8 SI narrow-vertex anomaly at 1.5 mrad
- `trigger_bit30_turnon.md` — ε_L1 = 99.58% at ET≥8 GeV
- `crosssection_cut_efficiency_audit.tex` — L1 bit & plateau correction audit; flags bit 26 vs 30 mismatch
- `mc_purity_correction_investigation.tex` — R-factor crosses 1.0 at pT≈14 GeV (genuine physics)
- `purity_nonclosure_ntbdt.tex` — NT-BDT non-closure mean 0.062–0.073
- `donut_iso_variants.tex` — 7 inner-R variants; all R≤0.1 recover more purity
- `iso_cone_fix.tex` — |vz| cone mis-centering fix (2026-04-18)
- `isolation_ET_comparison.tex` — full-calo vs EMCal-only
- `double_interaction_note.tex` — consolidated DI methodology (source for §4.4 efficiency decomposition)
- `bdt_training_audit.tex` — ET reweight flatness, per-variant outputs; `base_v1E` nominal 9 features confirmed
- `split_nosplit_dual_bdt.md` — dual cluster-node BDT refactor (22 models)
- `jetphox_audit.md` — CT14lo concern, JP/We drift drivers (iso, Jensen bin averaging, PDF)
- `jetphox_3pdf_comparison.md` — CT14nlo / NNPDF3.1 shifts vs CT14lo
- `jetphox_nlo_rerun_setup.md` — production pipeline staged (Phase 3b trigger)
- `xsec_weight_audit.md` — DI window fix post-audit (photon10/20 boundary 30 → 22 GeV)
- `ntbdtpair_numbers.md` — NT-BDT pair xsec table
- `condor_pipeline_review_run28.md` — pipeline infrastructure (not note-relevant)

### Wiki (codebase subtrees, all updated 2026-04-20)
- `wiki/pipeline/overview.md`, `full-run-range.md`, `01..05-*.md`
- `wiki/concepts/abcd-method.md`, `isolation-cuts.md`, `shower-shape-variables.md`, `unfolding.md`, `systematic-variations.md`, `double-interaction-efficiency.md`, `eta-edge-migration.md`, `deltaR-truth-matching.md`, `efficiency-deltaR-comparison.md`, `mc-purity-correction.md`
- `wiki/reference/config-schema.md`, `constants-sync.md`, `data-flow.md`, `mc-samples.md`, `file-inventory.md`, `pileup-fractions.md`
- `wiki/guides/running-full-pipeline.md`, `adding-a-systematic.md`, `common-issues.md`

### Rule files
- `.claude/rules/yaml-config.md` (lumi & period definitions)
- `.claude/rules/double-interaction.md` (DI pipeline)
- `.claude/rules/root-macros.md`, `python-analysis.md`

### Memory notes
- `project_nominal_allz_setup.md` — **nominal is allz (64.3718 pb⁻¹); 60-cm pins are cross-check only**
- `project_trigger_bit30_turnon.md` — ε_L1 = 99.58%, bit 30 active
- `project_lumi_list_provenance.md` — 60-cm pin vs Bit30Corr mismatch
- `project_truth_vertex_reweight.md` — closure χ²/ndf
- `project_tower_acceptance_maps.md`, `project_tower_masking_variants.md`
- `project_cut_efficiency_audit.md`
- `project_bdt_training_audit.md`
- `project_double_interaction.md`
- `project_recoeff_obsolete.md`
