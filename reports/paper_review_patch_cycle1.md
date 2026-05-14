# PPG12 Paper — Proposed Patch (Cycle 2, grouped by reviewer)

**Source PDF**: `sPHENIX_PPG12_Paper.pdf` (May 9 2026, 19 pages, line-numbered) — the version Jamie/Virginia/Justin reviewed.
**Target file**: `PPG12-Paper/main.tex` (currently 415 lines; will drift from PDF as edits land).
**Rules in force**: see memory `feedback_paper_revision_rules.md`. Cycle 2 of ≤3. Max 5 subagents per cycle. Only well-supported changes applied; uncertain ones left as `% TODO [reviewer X.N]`.

**Status of cycle 2**: PDF mapping + Justin ✓ (cycle 1); EMCal resolution + π⁰/η merger pT ✓ (cycle 2); Virginia, HERWIG+Emma, mass-peak files, critic audit still running in background.

**Shuonli's cycle-1 decisions (2026-05-14)**:
1. **Gluon-PDF wording (Jamie L20/L43/L368/L372 + Justin J2/J3/J14)**: option (b) — **soften to "depends on" everywhere**. Apply in abstract + intro + summary.
2. **J-L306 (Emma's data-driven MBD-eff writeup)**: leave as TODO item; do not block patch.
3. **J-Fig1 (pure-background histogram)**: TODO — Jamie likely wants a pure-background histogram in Fig. 1.
4. **J-L285 ("check with Emma")**: this comment is about **MBD efficiency** (L306 area), NOT the L285 escale shift sentence. Mapping corrected.

---

## Changelog (cycle 3 — applied 2026-05-14)

All ✅ CLEAR items from the cycle-2 patch applied to working tree (no git commit yet). Build verified: 19-page PDF, clean compile, all 36 changes rendered correctly. 11 UNRESOLVED items inserted as `\pendingTODO{...}` red markers with verbatim reviewer quotes.

**Files modified**:
- `PPG12-Paper/main.tex` — 28 text edits + 11 UNRESOLVED red TODO markers
- `PPG12-Paper/references.bib` — 6 author-field normalizations to "{XXX Collaboration}"; 2 new entries (`NNPDF:2021njg`, `sPHENIX:2025dET`)
- `PPG12-Paper/defs.sty` — `\pythia` macro (PYTHIA-8 → PYTHIA); `\RunOne` macro (Run 24 → RHIC Run 24)
- `plotting/paper/plot_paper_showershape.C` L199 — legend "Inclusive MC" → "Inclusive jet MC"
- `PPG12-analysis-note/selection.tex` L122–126 — stale 4% smearing prose replaced with canonical scheme summary

**Applied items by category**:

*Abstract*: full redraft to one acronym (sPHENIX), depends-on gluon-PDF wording, split run-on first sentence, uniform 64.4 pb⁻¹ precision.

*Intro (§I)*: STAR + ISR/FSR + truth-level UNRESOLVED markers; "first sPHENIX measurement"; cross-section hyphenation cleanup; gluon-PDF "depends directly on" + `Ichou:2010wc` citation; HARRISON ref; PHENIX "prompt-photon" wording (V1).

*Detector (§II)*: granularity 0.024 × 0.024 (matches note); two-dimensional projective geometry wording; "comprises" UNRESOLVED marker; "RHIC Run 24" via macro update.

*Trigger (§III)*: "online-calculated" UNRESOLVED marker; 8×8 vs 4×4 window flagged as TODO.

*Simulation (§III)*: HERWIG demoted to consistency check with 1.5%/3-6% agreement statement; smearing prose simplified (additive ET-dependent, no formula); truth-iso "excluding the photon itself" (drop false "excluding neutrinos"); muons UNRESOLVED marker; "simulation sample" wording; remove "hard, high-energy" → "high-energy"; explicit (1-f_PU)/f_PU pileup weighting.

*Reconstruction (§III)*: resolution σ_E/E definition; π⁰/η merger pT with red `\pendingTODO{11}`/`\pendingTODO{46}`; "centre" → "center"; "beam-pipe interactions"; "kinematic range of this analysis"; tight/non-tight BDT cuts updated to 0.8156 / -0.00156 / 0.6844 / +0.00156.

*Sideband (§III)*: Eq 2 wrapped in parens, comma not period; Eq 3 dot UNRESOLVED marker; redundant text trim; Padé equation; "a Padé fit" + no comma splice.

*Efficiency (§III)*: MBD-vertex eff direction reversed (59% → 52% across fiducial 12-32 GeV).

*Cross Section Determination (§III)*: subsection title capitalization; lumi uniform 64.4⁺⁵·⁹₋⁴·³ pb⁻¹.

*Systematics (§IV)*: escale prose expansion with co-dominance; eres simplified bracketing; "PPG10 ref TBD" UNRESOLVED marker for L321; three-treatment fix (eres → two-sided, unfolding-iteration → third-treatment max-symmetrized); L325 envelope sentence adds eres co-dominance; Fig 6 caption mentions lumi inclusion.

*Results (§V)*: middle/lower panel description added; truth-level fiducial UNRESOLVED marker on JETPHOX-comparison line; Fig 8 caption parens split.

*Summary (§VI)*: lumi 64.4 pb⁻¹; "depends on the gluon parton distribution function".

**11 UNRESOLVED items rendered as red `\pendingTODO` markers in the PDF**:
1. J-L25 ISR/FSR — intro
2. J-L37 STAR overlay — intro
3. J-L48 truth-level — intro (Jamie + Justin J5)
4. J-L67.window + J-J9 online — trigger sentence
5. J-J6 comprises — detector
6. J-L96 muons — simulation
7. J-J11 Eq 3 dot — sideband
8. J-L306 PPG10 ref — systematics
9. Justin J13 truth-level — results
10. J-L128 π⁰/η merger pT × 2 (red `\pendingTODO{11}`, `\pendingTODO{46}`)
11. J-Abs-1 x-range — abstract

**Build artifacts**: `PPG12-Paper/main.pdf` (19 pages, 498 KB, zero errors, one font warning). Build used `unsrt.bst` because `unsrturl.bst` is not installed in the local TeX-Live tree — final Overleaf build will use the correct style. The bibstyle is reverted to `unsrturl` in the committed `main.tex`.

**Not yet done in cycle 3**:
- Git commit (not requested).
- Overleaf push (not requested; needs OVERLEAF_TOKEN re-export).
- Plot cosmetic items (Fig 6 vertical lines, Fig 8 colors + stat band style): TODO for plot-rerun work outside the paper-text scope.
- JETPHOX NNPDF4.0 rerun for Fig 7 lower panel: separate condor-sweep task.

**Response letters drafted** at `reports/response_letters.md` (66 reviewer comments × verbatim quote + Addressed/Held/Not-changed response, ~5,300 words).

---

## Patch — grouped by reviewer comment

Each item lists: **PDF line(s)** → **current main.tex line(s)** → **proposed change** → **evidence status**.

Evidence status legend:
- ✅ **CLEAR** — change is directly supported by PDF/code/ROOT/bib evidence; safe to apply.
- ⚠️ **PARTIAL** — wording known, but tone/scope needs Shuonli confirm. Apply with the suggested wording; revisit if Shuonli objects.
- 🔴 **UNRESOLVED** — needs decision (held by Shuonli, or pending external input).
- ⏳ **CYCLE-1 PENDING** — agent still running; will be filled in cycle 2.

---

### JAMIE

#### J-Abs-1 — Abstract has too many defined acronyms
- **PDF L6–22**: abstract defines (Eγ_T), (p+p), (EMCal), (NLO), (pQCD), JETPHOX, MC, PYTHIA-8, PHENIX. (8 acronyms.)
- **Current main.tex L60–62**: same acronym density (still defines all of the above; PHENIX sentence already removed from current draft).
- **Critic fix**: reconcile with J-L20 gluon-PDF wording (use "depends on" everywhere per Shuonli's option (b)). "0.12 ≲ x ≲ 0.32" claim flagged for sign-off — the LO kinematic estimate x ≈ 2 ETγ/√s is standard but not currently stated in the paper.
- **Proposed change** (publication-ready, jargon stripped, gluon-PDF consistent):

```latex
\begin{abstract}
We report a measurement of the differential cross section of isolated prompt-photon production at midrapidity in $p+p$ collisions at $\sqrt{s}=200$~GeV. The data were recorded by the sPHENIX experiment at the Relativistic Heavy Ion Collider during the 2024 run, corresponding to an integrated luminosity of $64.4~\mathrm{pb}^{-1}$. The measurement spans photon transverse energies from 12 to 32~GeV. Photon candidates are reconstructed in the electromagnetic calorimeter and identified using a multivariate shower-shape classifier, with the residual jet background statistically subtracted using a sideband method. Photons are required to be isolated within $\Delta R = 0.3$, with the total additional transverse energy in the cone below 4~GeV at the truth level. The measured cross section is compared with predictions from the PYTHIA 8.307 Monte Carlo event generator and with two next-to-leading-order perturbative-QCD calculations. The result establishes a $p+p$ baseline for forthcoming sPHENIX measurements of isolated photons in heavy-ion collisions and depends on the gluon parton distribution function of the proton in the moderate-$x$ region.
\end{abstract}
```

- **Status**: ✅ CLEAR — folds J-J1 (split run-on first sentence) + J-L20 / J-J2 (gluon "depends on") + J-L18 (PHENIX comparison absent) + J-L10 (64.4) + Jamie-1 (acronyms: only sPHENIX remains undefined → defined in body). The "0.12 ≲ x ≲ 0.32" LO claim removed from abstract since it's not stated elsewhere in the paper; keep the kinematic-x discussion in §I intro after L78 instead, where the literature citation can support it. **TODO [Jamie-1.x-range]**: confirm with Shuonli whether the x range belongs in the abstract or only in the body.

#### J-L10 — Luminosity numeric precision consistency
- **PDF L10**: "corresponding to an integrated luminosity of 64.37 pb-1".
- **Current main.tex**:
  - Abstract L61: `$64.4\,\pb$` (1 dp, no uncertainty)
  - Intro L80: `$\mathscr{L}=64.4\,\pb$` (1 dp)
  - Lumi def L308: `$\mathscr{L}=64.37^{+5.88}_{-4.34}~\pb$` (2 dp central + 2 dp uncertainty)
  - Syst L321: `$\mathscr{L}=64.37^{+5.88}_{-4.34}~\mathrm{pb}^{-1}$` (2 dp)
  - Summary L371: `$\mathscr{L}=64.37~\pb$` (2 dp, no uncertainty)
- **Shuonli's decision (2026-05-14)**: use **64.4** everywhere; when the uncertainty is quoted, round to matched **1 dp** precision.
- **Precision arithmetic**: relative uncertainty is +9.13% / −6.75% (from the Vernier scan, fixed in earlier session). Applied to 64.4 pb⁻¹:
  - `+9.13% × 64.4 = +5.88` → rounds to **+5.9** (1 dp)
  - `−6.75% × 64.4 = −4.35` → rounds to **−4.3** (1 dp)
- **Proposed change** (uniform 1 dp everywhere):
  - L61 (abstract): keep `$64.4\,\pb$` — no change.
  - L80 (intro): keep `$\mathscr{L}=64.4\,\pb$` — no change.
  - L308 (lumi def): change `$\mathscr{L}=64.37^{+5.88}_{-4.34}~\pb$` → `$\mathscr{L}=64.4^{+5.9}_{-4.3}~\pb$`.
  - L321 (syst): change `$\mathscr{L}=64.37^{+5.88}_{-4.34}~\mathrm{pb}^{-1}$` → `$\mathscr{L}=64.4^{+5.9}_{-4.3}~\mathrm{pb}^{-1}$`.
  - L371 (summary): change `$\mathscr{L}=64.37~\pb$` → `$\mathscr{L}=64.4~\pb$`.
- **Status**: ✅ CLEAR — all five sites updated to uniform 1-dp precision.

#### J-L18 — Remove PHENIX comparison sentence from abstract
- **PDF L17–18**: *"A comparison with the previous PHENIX measurement of direct photons at the same collision energy is also presented."*
- **Current main.tex L60–62**: **already removed.**
- **Shuonli noted "Done"**.
- **Status**: ✅ CLEAR — no action.

#### J-L20 / J-L43 / J-L368 / J-L372 + Justin J2/J3/J14 — Soften gluon-PDF language to "depends on"
- **Shuonli's decision (2026-05-14)**: option (b) — soften to "depends on" everywhere. Drop "sensitivity to / sensitive to" framing across abstract, intro, summary.
- **Current main.tex L62 (abstract), L78 (intro), L376 (summary)**: all still use "sensitive to / sensitivity to".
- **Proposed changes** (3 locations, parallel wording):
  - **Abstract L62 (current)**: *"is sensitive to the gluon parton distribution in the proton, and establishes the \pp{} baseline for forthcoming sPHENIX measurements of isolated photons in heavy-ion collisions."*
    **→ After**: *"depends on the gluon parton distribution in the proton, and establishes the \pp{} baseline for forthcoming sPHENIX measurements of isolated photons in heavy-ion collisions."*
  - **Intro L78 (current)**: *"making the cross section particularly sensitive to the gluon parton distribution function (PDF) of the proton."*
    **→ After**: *"so that the cross section depends directly on the gluon parton distribution function (PDF) of the proton."*
  - **Summary L376 (current)**: *"sensitivity to the gluon parton distribution function"*
    **→ After**: *"a dependence on the gluon parton distribution function"*
- **Status**: ✅ CLEAR — supported by Shuonli's decision. Resolves Jamie L20/L43/L368/L372 and Justin J2/J3/J14 in one stroke.

#### J-L25 — ISR/FSR interference; fragmentation = FSR?
- **PDF L23–25**: *"Prompt photons refer to those produced either directly at the initial partonic hard-scattering, so-called direct photons, or from the collinear fragmentation of a final-state parton, so-called fragmentation photons."*
- **Current main.tex L70**: identical wording.
- **Shuonli noted "hold pending expert"**.
- **Status**: 🔴 UNRESOLVED — `% TODO [Jamie-L25]: hold for expert consultation on ISR/FSR phrasing; no edit applied`.

#### J-L30 — Reference [1] inside the period
- **PDF L29–30**: *"...collinear fragmentation of final-state partons into photons. [1]"* (citation appears AFTER the period in PDF render).
- **Current main.tex L73**: `photons~\cite{Aurenche:2006vj}.` — tex source has citation BEFORE the period.
- **Reading**: PDF rendering vs tex source mismatch. Bibstyle's `\bibpunct` or hyperref settings may push the citation past the period during compile. Need to compile and inspect.
- **Status**: ⚠️ PARTIAL — `% TODO [Jamie-L30]: tex source has \cite before period; check compiled PDF placement and adjust if hyperref renders citation post-period`.

#### J-L34 — "cross-section" vs "cross section" hyphenation
- **PDF L33–34**: *"Isolated prompt photon cross-sections..."* (hyphenated)
- **Current main.tex**: ~20 occurrences with no hyphen ("cross section"); subsection title at L289 reads "cross section Determination" (lowercase + no hyphen).
- **Proposed change**: APS house style is "cross section" (no hyphen). Sed pass through main.tex; rename subsection title at L289 to "Cross Section Determination" (capitalize Title Case).
- **Status**: ✅ CLEAR.

#### J-L37 — STAR 0912.3838 + newer measurements
- **PDF L35–37**: *"At the Relativistic Heavy-Ion Collider (RHIC) [13], the PHENIX experiment reported prompt photon cross-sections without an isolation requirement [14]."*
- **Current main.tex L75**: same wording, cites `HARRISON2003235` for RHIC and `PHENIX:2012jgx` for PHENIX.
- **Shuonli's decision (2026-05-14)**: **HOLD until the STAR points are overlaid on the existing PHENIX-comparison figure**. The citation should land together with the new plot, not as a standalone intro reference. Plotting work is needed in `plotting/paper/plot_paper_final.C` (the `final_phenix.pdf` macro) to read STAR data + overlay points with a third color/marker style.
- **Bibtex entry (staged but NOT yet applied)** — verified arxiv ID, DOI, journal. Will be added in the same commit as the plotting update:

```bibtex
@article{STAR:2009ojw,
    author        = "{STAR Collaboration}",
    collaboration = "STAR",
    title         = "{Inclusive $\pi^0$, $\eta$, and direct photon production at high transverse momentum in $p+p$ and $d+\mathrm{Au}$ collisions at $\sqrt{s_{NN}}=200$ GeV}",
    eprint        = "0912.3838",
    archivePrefix = "arXiv",
    primaryClass  = "hep-ex",
    doi           = "10.1103/PhysRevC.81.064904",
    journal       = "Phys. Rev. C",
    volume        = "81",
    pages         = "064904",
    year          = "2010"
}
```

- **Status**: 🔴 UNRESOLVED — bib entry verified, prose proposed, but commit blocked until Shuonli completes STAR-overlay plotting in `final_phenix.pdf`. **TODO [Jamie-L37]**: Shuonli to extract STAR data points from arXiv 0912.3838 (data table or HEPData), update `plot_paper_final.C` to overlay STAR alongside PHENIX, then land citation + intro-line update + bib entry + figure update in one commit.

#### J-L43 — Gluon PDF (intro)
- Folded into J-L20 above (Shuonli option (b)).
- **Status**: ✅ CLEAR.

#### J-L44 — "Article" not used in PRC/PRD; Justin overlaps with "first sPHENIX measurement"
- **PDF L43–46**: *"This Article reports the measurement of the differential cross section of isolated prompt-photon production..."*.
- **Current main.tex L80**: same.
- **Proposed change** (combines Jamie J-L44 and Justin J4):

```latex
This paper presents the first sPHENIX measurement of the differential cross section of isolated prompt-photon production as a function of \etg{} in \pp{} collisions at \comHEP{}, in the kinematic range $|\etag|<0.7$ and $12<\etg<32~\GeV$.
```

- **Justin's "first sPHENIX measurement" verified**: no prior sPHENIX isolated-photon paper on arxiv as of May 2026 (Justin's fact-check agent confirmed via BNL Newsroom, arxiv 2504.02242, arxiv 2604.15667).
- **Status**: ✅ CLEAR.

#### J-L48 — Define truth level à la ATLAS; Justin overlaps with J5
- **PDF L46–48**: *"Photons are required to be isolated at the truth level, defined by the total transverse energy of final-state particles within ∆R = 0.3 around the photon, excluding the photon itself, being below 4 GeV."*
- **Current main.tex L80**: same.
- **Shuonli's decision (2026-05-14)**: **HOLD pending expert consultation** on the canonical truth-level fiducial definition. This impacts J-L96 (Jamie's "muons?"), Justin J5 (truth-level jargon for data), and Justin J13 (truth-level in JETPHOX-comparison paragraph) — all four are linked through the same fiducial-definition prose.
- **Provisional proposed change (NOT to commit until expert sign-off)**:

```latex
The cross section is reported in a fiducial region defined at the truth level. Truth photons are stable photons originating directly from the hard scatter or from the parton-shower of an outgoing parton, excluding photons from hadron decays (such as $\pizero \to \gamma\gamma$ and $\eta\to\gamma\gamma$). A truth photon is declared isolated when the scalar sum of the transverse energies of all stable final-state particles within $\DR=0.3$ around the photon, excluding the photon itself, is below 4~\GeV. Experimentally, candidates passing a reconstruction-level isolation requirement (Sec.~\ref{ssec:iso}) are unfolded back to this truth-level fiducial.
```

- **Open issues to resolve with experts**:
  - ν+μ inclusion in the cone sum (per our HepMC-primary code, neutrinos are pre-decay so they don't enter; ATLAS convention explicitly excludes ν+μ at the stable-particle level).
  - Whether "stable" maps to PYTHIA status-1 or to GEANT-decayed stable particles.
  - Whether to give the JETPHOX parton-level vs particle-level equivalence sentence.
- **Status**: 🔴 UNRESOLVED — held pending expert consultation. **TODO [Jamie-L48]**: Shuonli to consult with theory + ATLAS-experienced collaborators on canonical truth-level fiducial wording, then commit J-L48/J-L96/J-J5/J-J13 together as one unified paragraph.

#### J-L67 + L70 — sPHENIX DAQ/trigger underweight, EMCal too brief
- **PDF L65–71**: detector + trigger description is two short paragraphs.
- **Current main.tex L87–103**: same brevity.
- **Cycle-2 fact-check findings**:
  - EMCal resolution **16%/√E ⊕ 5%** confirmed from `PPG12-analysis-note/reconstruction.tex` L4 citing `sPHENIX:2022emcal_testbeam` (NOT `Aidala:2020toz` as my v2 plan had — bibkey corrected).
  - Underlying testbeam paper (arXiv 2003.13685) measured (15.4%, 3.0%) and (13.3%, 3.5%) configurations. The "16% ⊕ 5%" form is a rounded-up wording the note adopts; the **constant term 5% likely reflects in-situ (not testbeam 3%)**. Need EMCal expert sign-off on the 5% in-situ.
  - **Trigger bit 30 is "Photon 4 GeV + MBD N&S coincidence"**, NOT a pure EMCal singles trigger. My v2 wording missed the MBD coincidence requirement (`selection.tex` L144, project memory `project_trigger_bit30_turnon.md`).
  - "Online-calculated" qualifier: NOT independently verified from local files. Likely correct per sPHENIX L1 firmware architecture, but flag as pending EMCal/L1 expert sign-off.

- **Proposed change** (corrected per cycle-2 findings):

```latex
The sPHENIX detector~\cite{PHENIX:2015siv, Belmont:2023fau} is a new central-rapidity detector at RHIC commissioned in 2024 and optimized for jet and heavy-flavor measurements in heavy-ion and \pp{} collisions. The detector is built around a 1.4~T superconducting solenoid and provides full azimuthal coverage for $|\eta|<1.1$ with a layered tracking system, an electromagnetic calorimeter (EMCal), and inner and outer hadronic calorimeters (IHCAL, OHCAL)~\cite{sPHENIX:2017lqb}. Forward minimum-bias triggering and primary-vertex reconstruction are provided by the Minimum-Bias Detector (MBD), a pair of Cherenkov-radiator quartz hodoscopes covering $3.51 < |\eta| < 4.61$ on each side of the interaction region.

The EMCal is a scintillating-fiber tungsten sampling calorimeter with two-dimensional projective geometry in $\eta$--$\phi$~\cite{Aidala:2020toz, sPHENIX:2017lqb} and approximately $0.024 \times 0.024$ $\eta$--$\phi$ segmentation. The towers are grouped into $2\times 2$ superblocks read out by silicon photomultipliers.

The level-1 trigger system uses hardware-summed EMCal energy in tower-window regions to form photon and jet trigger primitives. For the dataset reported here, photon candidates are selected on the EMCal photon trigger requiring an online-calculated tower-window sum above 4~GeV, with a sharp turn-on for cluster $\etg > 8.5$~\GeV. The minimum-bias trigger requires at least one charged-particle hit on each side of the MBD.
```

- **Shuonli's decisions (2026-05-14)**:
  1. Drop the single-particle energy-resolution sentence. The resolution is implied by the data-vs-MC smearing description in §III (J-L92 rewrite). Also clears the J-L67.5pct UNRESOLVED item.
  2. **The analysis trigger requirement is ONLY the online-calculated EMCal tower-window threshold — NOT in coincidence with an MBD hit.** Removed the "in coincidence with an MBD north-and-south hit" clause. The MBD minimum-bias trigger is described separately for detector-completeness, not as part of the photon-trigger selection.
- **Critic fixes applied**:
  - `sPHENIX:2022emcal_testbeam` NOT in `PPG12-Paper/references.bib`. Reverted to `Aidala:2020toz` (testbeam paper, arXiv 2003.13685) and `sPHENIX:2017lqb` (combined prototype) — both verified present in `references.bib`.
  - 4×4 vs 8×8 trigger window: agent says 4×4, current PDF says 8×8, code says 4×4. Pulled the specific number out of the proposed prose pending expert sign-off; using generic "tower-window sum" instead. **TODO [Jamie-L67.window]**: confirm trigger primitive size (4×4 or 8×8) with John Haggerty before final submission.
  - Granularity 0.024 (matches `reconstruction.tex` L4), not 0.025 — applied.
- **Evidence**: granularity 0.024 from `PPG12-analysis-note/reconstruction.tex` L4; trigger bit 30 = "Photon 4 GeV + MBD_NS" from `selection.tex` L144 and `config_bdt_nom*.yaml` (`trigger_used: 30`).
- **Open items requiring sPHENIX expert sign-off**:
  - "Online-calculated" qualifier on the L1 sum — ask John Haggerty / L1 firmware team.
  - 4×4 vs 8×8 trigger-window size — verify with sPHENIX trigger team.
- **Status**: ⚠️ PARTIAL — text drafted, two trigger-related expert-confirm items still open. Energy-resolution sentence DROPPED per Shuonli.

#### J-L77 — PYTHIA-8 8.307 redundant 8
- **PDF L76–77**: *"PYTHIA-8 8.307"*.
- **Current main.tex L107**: `\pythia~8.307` — the `\pythia` macro expands to "PYTHIA-8" (causing the duplicate "8").
- **Proposed change**: redefine the `\pythia` macro in `defs.sty` to just "PYTHIA" so `\pythia~8.307` renders as "PYTHIA 8.307". Or use `PYTHIA~8.307` directly.
- **Status**: ✅ CLEAR — `defs.sty` edit + body unchanged.

#### J-L92 — π⁰/η mass+width detail + figure; Virginia overlaps (V-L92)
- **PDF L90–92**: *"...convolved with an additional Gaussian of relative width σE/E = 4% [TODO: 4% pending the final EMCal calibration], applied as a multiplicative N(1, 0.04) factor on every MC cluster ET."*
- **Current main.tex L111**: same wording (STALE — code is now additive ET-dependent per commit `40d1c02`).
- **Proposed change** (publication-ready, tight prose per Shuonli — formula details dropped, no pending-calibration footnote):

```latex
The generated \pythia{} and \herwig{} events are propagated through the full sPHENIX detector using the \geant{} simulation package~\cite{AGOSTINELLI2003250} with detector-level noise included to match the data, and photons are reconstructed the same way as in data. To account for residual differences between data and simulation in the single-particle EMCal energy resolution, an additional ET-dependent Gaussian smearing is applied to the simulated cluster \et{} when filling the unfolding response matrix, with a width derived such that the simulated single-particle resolution matches the resolution measured in data. The unsmeared cluster \et{} is used for the photon-identification and isolation selections so that the additional resolution affects only the unfolding response.
```

- **Shuonli's decision (2026-05-14)**: drop the σ_extra formula, the (p_0, p_1, p_2) parametrization, the 17 GeV / 20 GeV threshold values, and the "pending the final EMCal calibration" footnote. The reader gets the physics intent (smearing tunes simulation to data); the detailed parametrization belongs in the analysis note.
- **Net effect**: clears the J-L92.thresholds TODO (no longer relevant since the thresholds aren't quoted).
- **Status**: ✅ CLEAR.

- **Evidence**: code at `efficiencytool/RecoEffCalculator_TTreeReader.C` L549, L557–577, L2125–2133, L2899–2907, L2960–2988; canonical note `PPG12-analysis-note/systematics.tex` L53–86; commit `40d1c02`.
- **π⁰/η mass+width figure** (Jamie's second part): **cycle-2 agent confirms NO mass-peak fit files exist inside PPG12** (`efficiencytool/results/`, `plotting/figures/`, `PPG12-analysis-note/Figures/` all empty for `pi0|eta_mass|peak|mgg|invmass`). External candidates exist in EMCal calibration WG areas (`/gpfs/mnt/gpfs02/sphenix/user/ecroft/Pi0Mass.pdf`, `/gpfs/mnt/gpfs02/sphenix/user/mrehman/merged_output/pi0_peak_*.pdf`) but are run-level peak monitoring, not data-vs-MC overlay. **TODO [Jamie-L92.fig]**: choose one of (a) defer figure to a later revision, (b) borrow an EMCal-WG plot (run-level only — partial answer to Jamie), or (c) run a new diphoton-pair invariant-mass study (new code in `BDTinput.C` / `MergeSim.C` / `ShowerShapeCheck.C`, none of which currently build diphoton pairs). Recommend (a) for this revision and add to revision-2 todo list.
- **Status**: ⚠️ PARTIAL — text rewrite is supported; figure pending file location.

#### J-L96 — Truth-level isolation definition: drop false "excluding neutrinos"
- **PDF L94–96**: *"At the truth level, signal photons are defined as prompt photons (including both direct and fragmentation photons) that satisfy the isolation requirement: the total transverse energy of all final state particles within ∆R = 0.3, excluding neutrinos and the photon itself, is below 4 GeV."*
- **Current main.tex L113**: same wording. The "excluding neutrinos" claim is **false** per code audit (`anatreemaker/source/CaloAna24.cc` L934–940 PID block commented out — no PID filter on the iso loop).
- **Shuonli's decision (2026-05-14)**: minimal truthful fix — just drop "excluding neutrinos and" → "excluding the photon itself". The fuller ATLAS-style fiducial definition (J-L48 bundle) remains held pending expert consultation, but the false claim must be corrected regardless.
- **Proposed change** at main.tex L113:

```latex
At the truth level, signal photons are defined as prompt photons (including both direct and fragmentation photons) that satisfy the isolation requirement: the total transverse energy of all final-state particles within $\DR=0.3$, excluding the photon itself, is below 4~\GeV. A reconstructed photon candidate is matched to a generated photon when that photon contributes the largest fraction of the simulated energy deposited in the candidate's EMCal cluster.
```

- **Net change**: 4-word deletion ("excluding neutrinos and"). Now consistent with the intro at PDF L46–48 / current main.tex L80, which already reads "excluding the photon itself" with no neutrino mention.
- **Evidence**: `anatreemaker/source/CaloAna24.cc` L874–959 (no PID filter applied to truth-iso sum); cycle-1 audit on 747,706 prompt photons confirms ν+μ contribution is <1 MeV per photon on average, so the practical impact of including them is negligible — no re-run needed.
- **Status**: ✅ CLEAR — minimal targeted edit. The broader J-L48 fiducial-definition expansion remains 🔴 UNRESOLVED.

#### J-L112 — Remove "hard"
- **PDF L110–112**: *"Second, fPU is the corresponding probability that an event producing a hard, high-energy cluster contains two or more p+p interactions."*
- **Current main.tex L117**: same.
- **Proposed change**: `hard, high-energy` → `high-energy`.
- **Status**: ✅ CLEAR.

#### J-L124 — "resolution" undefined
- **PDF L123–124**: *"Reconstructed photons have a resolution of approximately 6% across the entire Eγ_T range."*
- **Current main.tex L126**: same.
- **Proposed change**: *"Reconstructed photons have an energy resolution $\sigma_E/E$ of approximately 6\% across the entire \etg{} range, dominated by the EMCal stochastic and constant terms."*
- **Status**: ✅ CLEAR.

#### J-L128 — pT at which π⁰/η photons merge
- **PDF L126–128**: shower-shape introduction paragraph.
- **Current main.tex L128**: same.
- **Cycle-2 fact-check finding**: at midrapidity with EMCal tower Δθ ≈ 0.024 (from `reconstruction.tex` L4) and m_π⁰ = 135 MeV, m_η = 548 MeV (PDG):
  - π⁰ merger onset: p_T > 2m/Δθ = 2 × 0.135 / 0.024 = **11.2 GeV**
  - η merger onset: p_T > 2 × 0.548 / 0.024 = **45.7 GeV**
  - (No hard-coded merger pT in `anatreemaker/source/CaloAna24.cc`; the cluster splitter handles merging case-by-case via `cluster_merged_prob`.)
- **Shuonli's decision (2026-05-14)**: add the sentence but tag the numerical values as `\pendingTODO{}` red placeholders so they can be tightened later.
- **Proposed change** (append to L128):

```latex
The two photons from $\pizero \to \gamma\gamma$ become indistinguishable in a single EMCal cluster above $p_T \approx \pendingTODO{11}$~\GeV; for $\eta \to \gamma\gamma$ the corresponding threshold is $p_T \approx \pendingTODO{46}$~\GeV. The merging is gradual rather than a hard cut.
```

- **Evidence (used to derive the placeholder values)**: granularity 0.024 (`PPG12-analysis-note/reconstruction.tex` L4); EMCal inner radius 93.5 cm (`efficiencytool/DoubleInteractionCheck.C` L797); PDG masses. The geometric derivation can be tightened with shower-shape simulation if needed.
- **Status**: ✅ CLEAR — prose committed with placeholder numbers in red.

#### J-L132 — "centre" → "center"; hyphens
- **PDF L130–133**: *"centre-of-gravity"*, mixed with later hyphenated/non-hyphenated forms.
- **Current main.tex L128, L136, L137, L139, L143**: "centre"; L172: "center".
- **Proposed change**: sed-pass `centre` → `center` everywhere; keep `center-of-gravity` hyphenated (standard compound adjective).
- **Status**: ✅ CLEAR.

#### J-L138 — Upstream-muon physical insight (z-axis?)
- **PDF L134–138**: NCB paragraph.
- **Current main.tex L152**: same.
- **Shuonli's decision (2026-05-14)**: drop "beam-halo and" — just say "beam-pipe interactions". Also avoid semicolon per project writing preference.
- **Proposed change** (append to NCB paragraph):

```latex
These NCB clusters originate predominantly from beam--pipe interactions upstream of the interaction region. The resulting muons propagate approximately parallel to the beam direction and produce clusters with elongated $\eta$-like shower shapes and arrival times shifted earlier than the in-time beam crossing.
```

- **Status**: ✅ CLEAR — phrasing tightened to "beam-pipe interactions" only, matching what the analysis-note NCB study supports.

#### J-L144 — What other observables (NCB BDT features); Virginia overlaps (V-L144)
- **PDF L142–144**: *"using the shower-shape variables listed in Table 1 together with additional shower-shape observables."*
- **Current main.tex L154**: same.
- **Shuonli's note**: response-letter only, NOT in paper.
- **Proposed change**: no paper edit. Response letter will list the additional features pulled from `FunWithxgboost/config.yaml`.
- **Status**: ✅ CLEAR — no paper edit. **TODO [Jamie-L144]**: draft response-letter text enumerating 17 additional features.

#### J-Fig1 — Show background shape too (already shown, just mislabeled)
- **PDF Fig. 1 caption**: *"the blue histogram shows the inclusive jet background MC"* — caption already calls it background.
- **Current legend in `dis_*_cut1.pdf`**: just says "Inclusive MC" (verified by reading the rendered PDF text).
- **Code reality** (`plotting/paper/plot_paper_showershape.C` L66, L199): `proj_bkg` is read from `MC_efficiencyshower_shape_jet_inclusive_combined_*.root` — that IS the inclusive jet MC (all jet pT-hat bins combined, weighted by cross-section). The blue histogram **is** the pure background shape.
- **Shuonli's decision (2026-05-14)**: the background is already shown; Jamie's concern is the misleading legend label. Just relabel "Inclusive MC" → "Inclusive jet MC" (or "Background MC") so it's unambiguous.
- **Proposed change** (`plotting/paper/plot_paper_showershape.C` L199):

```diff
- leg->AddEntry(proj_bkg, "Inclusive MC", "l");
+ leg->AddEntry(proj_bkg, "Inclusive jet MC", "l");
```

And in main.tex Fig. 1 caption (L162–168), confirm the caption reads "inclusive jet background MC" to match. Currently already does.
- **Status**: ✅ CLEAR — single-line legend relabel + regenerate `dis_*_cut1.pdf`. No new histogram needed.

#### NEW — Tight / non-tight BDT cuts in paper are stale (audit-only finding)
- **Paper PDF L179–180 / current main.tex L179–181**: still has the old cuts
  - Tight: `ID BDT > 0.833 − 0.00333 ETγ`
  - Non-tight: `0.733 − 0.01333 ETγ < ID BDT < 0.667 + 0.00333 ETγ`
- **Analysis note `reconstruction.tex` L260, L265** (canonical):
  - Tight: `BDT score > 0.8156 − 0.00156 · ETγ`
  - Non-tight: `0.7333 − 0.01333 · ETγ < BDT score < 0.6844 + 0.00156 · ETγ`
- **Production config `efficiencytool/config_bdt_nom.yaml`** (L188–189, L244–245, L248–249):
  - Tight: intercept `0.815625`, slope `-0.0015625`
  - Non-tight upper: intercept `0.684375`, slope `+0.0015625`
  - Non-tight lower: intercept `0.7333…`, slope `-0.01333…`
- **Three numbers drifted** (only the non-tight lower bound matches):
  - Tight intercept: 0.833 (paper) → 0.8156 (note/code)
  - Tight slope: -0.00333 (paper) → -0.00156 (note/code)
  - Non-tight upper intercept: 0.667 (paper) → 0.6844 (note/code)
  - Non-tight upper slope: +0.00333 (paper) → +0.00156 (note/code)
- **Proposed change** at main.tex L179–181:

```latex
The ``tight'' photon identification (\gid) selection is defined by an \etg{}-dependent BDT-score threshold, optimized to maintain an approximately 80\% identification (ID) efficiency relative to the pre-selection across the analysis range. The tight selection requires $\mathrm{ID\ BDT} > 0.8156 - 0.00156\,\etg.$
A ``non-tight'' \gid{} selection is defined by $0.7333 - 0.01333\,\etg < \mathrm{ID\ BDT} < 0.6844 + 0.00156\,\etg,$
which enriches the background photon contribution and is used as the photon-ID sideband in the data-driven purity calculation described in Sec.~\ref{ssec:abcd}.
```

- **Also check**: `systematics.tex` L97–98 has commented-out variant cuts that reference the old (?) numbers — verify and confirm they correspond to the variation around the **new** nominal, not the old one. Quick check: L97 says `0.8156 − 0.00156` (matches new nominal), L102 says variation is `0.8156 - 0.00156·ETreco + 0.05` (a flat ±0.05 shift around the new nominal). Consistent with new cuts. ✓
- **Evidence**: `PPG12-analysis-note/reconstruction.tex` L260,L265; `efficiencytool/config_bdt_nom.yaml` L186–249.
- **Status**: ✅ CLEAR — purely a paper-prose update; production code and analysis note are already on the new cuts.

#### J-L179 — Show BDT distribution data vs MC; Virginia overlaps
- **PDF L177–179**: tight-γID definition.
- **Current main.tex L179**: same.
- **Proposed change**: new figure paralleling Fig. 1 showing the BDT score distribution (data + signal MC + inclusive MC) at one representative pT slice, with the tight-BDT threshold drawn.
- **Status**: ⚠️ PARTIAL — new plotting macro needed (not a text edit). **TODO [Jamie-L179]**: implement plot using existing `cut1` infrastructure.

#### J-Fig3 — Side-by-side schematic + 2D data heatmap
- **Status**: ⚠️ PARTIAL — new plotting macro. **TODO [Jamie-Fig3]**.

#### J-L236-237 — Restates figure, remove
- **PDF L234–237**: redundant phrasing about leakage correction.
- **Current main.tex L257**: same wording.
- **Proposed change**: trim L257. Keep only *"The signal leakage-corrected purity ranges from approximately 0.5 at 10 GeV to approximately 0.9 around 30 GeV."*
- **Status**: ✅ CLEAR.

#### J-Pade-eq — Add equation for P[1,1](x)
- **Proposed change**: insert *"a [1/1] Padé approximant, $P_{[1/1]}(x) = (a_0 + a_1 x)/(1 + b_1 x)$, is fitted to smooth statistical fluctuations."*
- **Status**: ✅ CLEAR.

#### J-L266 — MBD-vertex efficiency clarification; Shuonli: no new figure
- **PDF L263–266**: *"...ranges from approximately 55% to 60% across the measured ETγ range."*
- **Current main.tex L280**: same. **Number is wrong both in direction and approximate values per ROOT-file extraction**: g_mbd_eff = 58.68% at ETγ=13 GeV down to 51.80% at ETγ=30 GeV (verified).
- **Proposed change** (concise — text-only, no new figure):

```latex
The MBD-vertex efficiency, $\varepsilon_{\mathrm{vtx}}$, is defined as the fraction of truth-level isolated signal photon events in which the MBD records at least one hit on each side and yields a reconstructed vertex within $|z_\mathrm{reco}|<60$~cm. The denominator uses the all-$z$ truth-vertex fiducial, so this efficiency converts the beam-delivered luminosity into the analysis-fiducial vertex acceptance. The efficiency decreases from approximately $59\%$ at $\etg\!\simeq\!12$~\GeV{} to approximately $52\%$ at $\etg\!\simeq\!30$~\GeV. The photon trigger efficiency is at least $99.5\%$ for $\etg>8.5$~\GeV{} and is applied per cluster before the remaining efficiency corrections, as described in Sec.~\ref{ssec:eventsel}. All efficiencies are estimated from MC simulation and applied as bin-by-bin corrections to the unfolded \etg{} spectrum.
```

- **Shuonli's decisions (2026-05-14)**:
  - Dropped the "(covering $3.51<|\eta|<4.61$)" MBD pseudorapidity coverage parenthetical — already stated in the detector section.
  - Dropped the "reflecting the falling forward-rapidity particle density..." physics-trend clause — the numerical trend speaks for itself.
- **Status**: ✅ CLEAR — numbers verified against `efficiencytool/results/Photon_final_bdt_nom.root::g_mbd_eff`.

#### J-L285 — Escale prose expansion (mapping resolved: Shuonli's "check with Emma" is about MBD eff, lives in J-L306)
- **PDF L285–287**: *"The reconstructed EMCal cluster energy in MC is shifted by ±1.5% [TODO: Placeholder pending the final EMCal calibration]..."*. **(PDF version)**
- **Current main.tex L315**: *"...shifted by ±1.1% to account for the residual data-MC difference..."* (drift: ±1.5% → ±1.1% in current main.tex).
- **Shuonli's decision (2026-05-14)**: confirmed — "check with Emma" maps to **L306 MBD-eff**, NOT to this escale sentence. J-L285 is therefore unambiguously about escale-prose expansion.
- **Proposed change** for L315 escale (Justin J28-equivalent context):

```latex
The reconstructed EMCal cluster transverse energy in MC is varied by a flat multiplicative factor of $1 \pm 0.011$, bracketing the residual data--MC offset of the reconstructed \pizero{} and $\eta$ peak positions measured across the EMCal acceptance. The variation is applied per cluster, before any identification or isolation requirements, and is propagated through the full analysis chain to give a symmetric two-sided systematic, added in quadrature with the other contributions. This uncertainty is among the dominant systematic contributions at high \etg{}, together with the energy-resolution uncertainty.
```

- **Shuonli's decision (2026-05-14)**: changed "is the dominant systematic" → "is among the dominant systematic contributions, together with the energy-resolution uncertainty" — at high ETγ the escale and eres systematics are now co-dominant (per the post eres-max-symmetrize budget where eres reaches ±21.2% at high ETγ).
- **Evidence**: `make_bdt_variations.py` L213–216 (energyscale11up/down at clusterescale=1.011/0.989); `config_bdt_nom.yaml` L95 (`cluster_escale: 1.0` nominal); `PPG12-analysis-note/systematics.tex` L85 (eres max ±21.2% at high ETγ).
- **Status**: ✅ CLEAR.

#### J-L306 — Vague sentence (Shuonli: "check with Emma" maps HERE, not L285)
- **PDF L304–306**: *"The efficiency uncertainty is taken from variations of the empirical Eiso_T correction described in Sec. 3.3, together with a separate envelope assigned to the MBD-vertex efficiency based on the data-driven MBD-vertex efficiency study."*
- **Current main.tex L321**: same.
- **Shuonli's decision (2026-05-14)**: **HOLD**. The data-driven MBD-efficiency systematic and its source (PPG10 dijet analysis) need to be settled before this sentence can be tightened.
- **Reference context (preserved for cycle 3)**: there is no author "Emma" in the analysis note. The reference is to **PPG10** (sPHENIX dijet analysis), already cited in `PPG12-analysis-note/systematics.tex` L202–205 with text *"taken from the data-driven MBD-coincidence efficiency reported in [PPG10, ref TBD]"*. Currently in the budget table as **±6.00% placeholder** (`systematics.tex` L306).
- **PPG10 bibtex entry**: MISSING from `references.bib` and `PPG12-analysis-note/cite.bib`. Need Shuonli to obtain the proper PPG10 paper bibkey.
- **Provisional rewrite (NOT to commit until hold lifts)**:

```latex
A flat $\pm 6\%$ uncertainty on the MBD-vertex-finding efficiency is assigned, taken from the data-driven MBD-coincidence efficiency measurement in the sPHENIX dijet analysis~\cite{PPG10}. This systematic is transferred from the dijet to the prompt-photon analysis since the MBD-vertex efficiency is driven by forward-rapidity underlying-event activity and is therefore common between QCD-dijet and prompt-photon events at the same leading truth $p_T$. The MBD-vertex-found fraction in PYTHIA-8 photon and jet samples agrees to within 1\% across $5 < p_T^\mathrm{lead} < 35$~\GeV{}, justifying this transfer.
```

- **Status**: 🔴 UNRESOLVED — held. **TODO [Jamie-L306]**: Shuonli to provide PPG10 bibtex key; resolve the red placeholders `[PPG10, ref TBD]` and `[$\pm 6\%$ placeholder]` in `PPG12-analysis-note/systematics.tex` first.

#### J-PPG09 — Truth-level Gaussian smearing; closure test
- **Per v2 plan**: NOT needed — the code already smears at the response-matrix fill stage with width sized by `pT_truth`. The systematic envelope (zero data params vs wider data params) is the natural resolution-uncertainty closure.
- **Status**: ✅ CLEAR — response-letter explanation: *"Our current smear sizes σ_extra by truth pT, which is the PPG09-equivalent prescription; the response-matrix-only application avoids double-counting in the bin-by-bin efficiency."*

#### J-Fig6 — Vertical lines at x-axis edges
- **Proposed change**: call `SetErrorX(0)` in the syst-sum plot macro, OR rebin overflow to uniform 2-GeV bins. Regenerate `syst_sum_rel.pdf`.
- **Status**: ⚠️ PARTIAL — cosmetic, plotting macro edit. **TODO [Jamie-Fig6]**.

#### NEW (cycle-3.1 audit) — L325 envelope sentence: mention eres co-dominance
- **Current main.tex L325**: *"The total systematic uncertainty has a strong bin-to-bin correlation, dominated by the energy-scale uncertainty at high \etg{} and by the purity-closure uncertainty at low \etg{}. The bin-by-bin envelope reaches approximately $-31\%/+32\%$ in the highest \etg{} region."*
- **Audit finding**: contradicts the J-L285 update that explicitly states escale is "co-dominant with the energy-resolution uncertainty" at high ETγ. The L325 sentence must mirror this.
- **Proposed change** at L325:

```latex
The total systematic uncertainty has a strong bin-to-bin correlation, dominated by the energy-scale and energy-resolution uncertainties at high \etg{} and by the purity-closure uncertainty at low \etg{}. The bin-by-bin envelope reaches approximately $-31\%/+32\%$ in the highest \etg{} region.
```

- **Net change**: insert "and energy-resolution" so the L325 sentence is consistent with the J-L285 escale rewrite.
- **Status**: ✅ CLEAR.

#### J-L324 — HERWIG in systematics
- **PDF L323–324**: *"The bin-by-bin envelope reaches approximately ±25% in the highest Eγ_T region."* — HERWIG not mentioned in PDF systematics section either.
- **Current main.tex L325**: envelope is "-31%/+32%" (drift from PDF's "±25%").
- **Cycle-2 fact-check finding**: my v2 proposal *"agree within statistical precision"* is **WRONG**. The analysis note `systematics.tex` L191–194 + L305 actually finds:
  - ε_iso: HERWIG sits **3–6% higher** than PYTHIA → propagated as a **±5.3% max systematic in [12,14] GeV**, **5.83% in the budget table**
  - ε_ID: agrees within ±1.5% (this part is fine)
  - ε_MBD: HERWIG sits **4–7% above** PYTHIA (analysis-note appendix L590)
- **Proposed change** (now corrected, supported by `systematics.tex` L191–194):

```latex
Prompt-photon samples are also generated with \herwig~\cite{Bellm:2015jjp} as a generator cross-check on the photon-identification and isolation efficiencies. The photon-identification efficiency agrees with the \pythia{} estimate within 1.5\%, while the isolation efficiency differs by 3--6\%; this residual difference is propagated as a generator-modelling systematic uncertainty (Sec.~\ref{sec:systematics}).
```

- **Note**: This change also fixes Jamie's underlying complaint — HERWIG IS in the systematics; the budget table at `systematics.tex` L305 carries a 5.83% entry. The current paper §IV systematics description should also be updated to list "generator modelling (HERWIG vs PYTHIA)" as one of the variations.
- **Evidence**: `PPG12-analysis-note/systematics.tex` L191–194 ("ε_iso sits 3–6% higher in HERWIG", propagated as ±5.3% max in [12,14] GeV), L305 (5.83% in budget table); `PPG12-analysis-note/appendix.tex` L573–590 (full HERWIG appendix); HERWIG samples on disk at `/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon{5,10,20}_herwig/`; configs `efficiencytool/config_bdt_herwig_1p5mrad*.yaml`; results `efficiencytool/results/Photon_final_herwig*.root`; pipeline driver `efficiencytool/oneforall_tree_herwig.sh`.
- **Status**: ✅ CLEAR — text now matches the actual analysis-note quantification.

#### J-L338-mid + L338-bot — Middle panel + bottom panel of Fig 7
- **Proposed change** (middle): one sentence between current L340 and L346.

```latex
The middle panel shows the ratio of each prediction to the data, with the experimental total uncertainty drawn as a band around unity and theoretical scale plus PDF uncertainties shown as boxes around each prediction.
```

- **Proposed change** (bottom): paragraph about PDF set comparison. Depends on NNPDF4.0 rerun (#35 in v2) — if not yet done, keep NNPDF3.1 text.

```latex
The lower panel shows JETPHOX predictions divided by the data for four alternative PDF sets (CT18NLO, NNPDF3.1, CTEQ6.6, MSHT20NLO). The four sets bracket the data within experimental uncertainties across the measured \etg{} range.
```

- **Status**: ⚠️ PARTIAL — middle-panel sentence ✅ CLEAR; bottom-panel paragraph contingent on NNPDF4.0 update. **TODO [Jamie-L338]**: include NNPDF4.0 once rerun lands.

#### J-NNPDF4.0 — Update to NNPDF4.0; add Hessian band
- Add `NNPDF:2021njg` (verified arxiv 2109.02653, EPJC 82 428):

```bibtex
@article{NNPDF:2021njg,
    author        = "{NNPDF Collaboration}",
    collaboration = "NNPDF",
    title         = "{The path to proton structure at 1\% accuracy}",
    eprint        = "2109.02653",
    archivePrefix = "arXiv",
    primaryClass  = "hep-ph",
    doi           = "10.1140/epjc/s10052-022-10328-7",
    journal       = "Eur. Phys. J. C",
    volume        = "82",
    number        = "5",
    pages         = "428",
    year          = "2022"
}
```

- JETPHOX rerun with NNPDF4.0: ~3-hour wallclock per scale/PDF. Pre-resubmission.
- Hessian band on lower panel: **DEFERRED to revision 2**.
- **Status**: ✅ CLEAR for bib add; ⚠️ PARTIAL for JETPHOX rerun (separate work). **TODO [Jamie-NNPDF]**: queue JETPHOX NNPDF4.0 run.

#### J-Ichou — Gluon-PDF sensitivity literature citation
- Even with Shuonli's option (b) "depends on" softening, citing the established gluon-PDF / direct-photon sensitivity literature would strengthen the wording.
- **Critic fix**: bibkey `Ichou:2010wc` already exists at `references.bib:15` (same arXiv 1005.4529 / PRD 82 014015). DO NOT add a duplicate entry; use the existing key.
- **Cite at intro L78** after the gluon-PDF "depends on" sentence:

```latex
...so that the cross section depends directly on the gluon parton distribution function (PDF) of the proton~\cite{Ichou:2010wc}.
```

- **Status**: ✅ CLEAR — citation key verified to exist; no bibtex entry needs to be added.

#### J-Fig8-stat — Stat uncertainties around 1
- **Per audit**: stat band already drawn around y=1 but light blue and hidden by theory boxes.
- **Proposed change**: darken band, add explicit legend entry "sPHENIX (stat ⊕ syst)" in ratio panels.
- **Status**: ⚠️ PARTIAL — plotting cosmetics. **TODO [Jamie-Fig8.stat]**.

#### J-Fig8-parens — "(Details...)" parens
- **PDF Fig. 8 caption**: *"...is overlaid (Details of the correction procedure are described in the text.)."*
- **Current main.tex L358**: same.
- **Proposed change**: split into separate sentence:

```latex
\caption{... is overlaid. Details of the correction procedure are described in the text.}
```

- **Status**: ✅ CLEAR.

#### J-L368, L372 — Gluon PDF again
- Folded into J-L20 above (Shuonli option (b) — "depends on" in summary).
- **Status**: ✅ CLEAR.

#### J-Refs — Refs 14, 15 + 19, 24, 25 Collaboration consistency
- **Proposed bib changes** (verified line numbers in `references.bib`):

| bibkey | author@line | current | replace with |
|---|---|---|---|
| `PHENIX:2010vgy` | 30 | `"Adare, A. and others"` | `"{PHENIX Collaboration}"` |
| `ATLAS:2011ezy` | 218 | `"Aad, Georges and others"` | `"{ATLAS Collaboration}"` |
| `PHENIX:2015siv` | 807 | `"Adare, A. and others"` | `"{PHENIX Collaboration}"` |
| `Aune:2024ysr` | 859 | `"Aune, S. and others"` | `"{sPHENIX Collaboration}"` + add `collaboration = "sPHENIX",` |
| `Aidala:2020toz` | 871 | `"Aidala, C. A. and others"` | `"{sPHENIX Collaboration}"` + add `collaboration = "sPHENIX",` |
| `sPHENIX:2017lqb` | 884 | `"Aidala, C. A. and others"` | `"{sPHENIX Collaboration}"` |

- **Status**: ✅ CLEAR.

---

### VIRGINIA (cycle-2, fact-checked)

#### V1 — PHENIX is "prompt (direct + fragmentation)" not "direct"
- **Evidence**: PHENIX paper arXiv 1205.5533 (PRL 108 081402) p. 3–4: *"Direct photons are defined as photons that do not originate from hadronic decays."* and *"Beyond leading order, direct photons can be produced either by bremsstrahlung from any quark line in a 2-to-2 subprocess, e.g. g+q→g+q+γ, or in a parton shower from fragmentation..."*. The PHENIX inclusive measurement therefore includes fragmentation — which in PPG12 terminology is "prompt", not "direct".
- **Current main.tex inconsistency**: L75 already calls PHENIX "prompt"; L80 still calls it "direct-photon".
- **Proposed change** at L80:

```latex
...and to the previous inclusive prompt-photon measurement by PHENIX~\cite{PHENIX:2012jgx}.
```

- Optional parenthetical at L351 (PHENIX comparison paragraph):

```latex
The PHENIX measurement was performed without an isolation requirement (and so corresponds to inclusive prompt photons in the present terminology), in a narrower acceptance...
```

- **Status**: ✅ CLEAR.

#### V2 — Reference for calorimeter calibration details
- **Cycle-2 finding**: no sPHENIX paper with the "PPG03" tag exists; the closest published reference is arXiv 2504.02242 (sPH-BULK-2025-02, Phys. Rev. C 2025) which establishes the EMCal energy-scale calibration via iterative π⁰-mass fitting.
- **Proposed bib entry** (new, not in `references.bib`):

```bibtex
@article{sPHENIX:2025dET,
    author        = "{sPHENIX Collaboration}",
    collaboration = "sPHENIX",
    title         = "{Measurement of the transverse energy density in Au+Au collisions at $\sqrt{s_{NN}} = 200$ GeV with the sPHENIX detector}",
    eprint        = "2504.02242",
    archivePrefix = "arXiv",
    primaryClass  = "nucl-ex",
    journal       = "Phys. Rev. C",
    volume        = "112",
    pages         = "024908",
    year          = "2025",
    reportNumber  = "sPH-BULK-2025-02"
}
```

(Cycle-3.1 audit fix: dropped the non-standard DOI placeholder `10.1103/h8d5-swg6`; canonical PRC 112 024908 metadata included instead. DOI can be added back as `10.1103/PhysRevC.112.024908` if verified against InspireHEP.)

- **Proposed cite location**: insert in §III "MC simulation" near the smearing sentence, *"The EMCal energy scale and resolution are established with an iterative $\eta$-dependent $\pi^0$-mass calibration~\cite{sPHENIX:2025dET}."*
- **Status**: ✅ CLEAR — bibkey + arxiv ID verified.

#### V3 — Combine SI + DI samples explicit
- **Cycle-2 finding**: `efficiencytool/oneforall_tree_double.sh` L34–86 and `RecoEffCalculator_TTreeReader.C` L104 confirm: per-event weight = cross-section × mix_weight × vtx_reweight, with mix_weight = SINGLE_FRAC for `_nom` and DOUBLE_FRAC for `_double`. The current paper L117 sentence ("simulated photon yield is constructed from the luminosity-weighted combination of single-interaction and pile-up samples") uses awkward "photon yield" phrasing and is also ambiguous about the pileup-fraction weighting.
- **Shuonli's decision (2026-05-14)**: rephrase from "simulated photon yield" → "the simulation sample" (the construction is over MC samples, not directly over a yield).
- **Proposed change** at L117 (replace the existing sentence):

```latex
The simulation sample used in this analysis is constructed within each crossing-angle period by combining the single-interaction and pile-up samples weighted by $(1-f_{\mathrm{PU}})$ and $f_{\mathrm{PU}}$ respectively, and the two crossing-angle periods are then combined by their integrated luminosities.
```

- **Status**: ✅ CLEAR — makes the SI/DI combination explicit AND fixes the "photon yield" phrasing per Shuonli.

#### V4 — L137/L152 wording
- **Current main.tex L152**: *"...non-negligible relative to the physics signal at the analysis $\etg{}$."*
- **Proposed change**:

```latex
...non-negligible relative to the physics signal in the kinematic range of this analysis.
```

- **Status**: ✅ CLEAR.

#### V5 — Specify number of iteration variations
- **Cycle-2 finding**: `make_bdt_variations.py` L266–277 — variants are `unfold_iter1` (cross-check only, NOT in systematic), `unfold_iter3` (max envelope), `unfold_iter4` (max envelope). Nominal is iter=2.
- **Proposed change** at L321 (systematics paragraph):

```latex
The unfolding uncertainty is evaluated by repeating the unfolding without the prior reweighting and by varying the number of Bayesian iterations to three and four (the nominal is two).
```

- **Status**: ✅ CLEAR.

#### V6 — Three treatments listed (paper wording wrong: third treatment is unfolding-iteration max, not energy resolution)
- **Cycle-2 fact-check**: commit `049e7df` added a third treatment to main.tex L323 but mis-identified what the third treatment was. Compared against canonical `systematics.tex` L21:
- **Shuonli's clarification (2026-05-14)**: the third treatment is the **unfolding iteration max-deviation**, NOT energy resolution. Energy resolution belongs in the two-sided per-bin group.

| Group | Paper L323 (current — wrong) | Note `systematics.tex` L21 (canonical) |
|---|---|---|
| Two-sided per-bin asymmetric | energy scale, non-isolation boundary, NCB cut, purity-68% CI, luminosity | energy scale, **energy resolution**, non-isolated boundary, NPB cut, luminosity, purity-fit CI |
| One-sided symmetrized | tight γID, non-tight γID, fit form, MC closure, iso eff, unfolding reweighting, **unfolding iteration**, pile-up fraction | tight ID, non-tight ID, fit form, MC closure, iso eff, response reweighting, double-interaction blending fraction |
| Third treatment (max-symmetrized) | (incorrectly) "energy resolution" | **unfolding iteration scan: max{iter3, iter4} symmetrized** |

- **Proposed L323 rewrite** (align to note):

```latex
Each variation falls into one of three treatments. Two-sided sources (energy scale, energy resolution, non-isolation boundary, NCB cut, purity $68\%$ confidence interval, luminosity) yield asymmetric up and down deviations that are kept separate per \etg{} bin. One-sided sources (tight \gid{} threshold, non-tight \gid{} threshold, fit functional form, MC closure correction, isolation efficiency, unfolding reweighting, pile-up fraction) yield a single deviation that is symmetrized about the nominal. The unfolding iteration uncertainty is taken as the per-bin maximum of the iter-3 and iter-4 deviations, symmetrized about the nominal.
```

- **Net changes**:
  - Move **energy resolution** from third-treatment to two-sided.
  - Move **unfolding iteration** from one-sided to its own third-treatment max-symmetrized group.
  - Drop "energy resolution is treated as the bin-by-bin maximum of the two bracket variations, symmetrized" sentence.
  - Keep "pile-up fraction" naming in the paper (consistent with the pile-up introduction in §III; the note's "DI blending fraction" is internal-language).
- **Evidence**: `PPG12-analysis-note/systematics.tex` L21 (canonical taxonomy).
- **Status**: ✅ CLEAR — paper L323 rewrite proposed, aligned with note. No re-aggregation needed (this is paper prose, not the syst quadrature).

#### V7 — Total systematic includes/excludes luminosity
- **Cycle-2 finding**: `plotting/paper/plot_paper_systematics.py` L199–203 confirms `add_flat(*quadrature_sum(total_comps), h_nom, FLAT_SYSTS)` adds the asymmetric lumi (+9.13%/-6.75% per `make_bdt_variations.py` L487–488) to the Total envelope but NOT as a separate breakdown line.
- **Proposed change** at Fig. 6 caption (L331):

```latex
\caption{Breakdown of relative systematic uncertainties as a function of \etg. The Total envelope includes the global luminosity uncertainty (added in quadrature, bin-independent at $^{+9.1\%}_{-6.8\%}$) in addition to the seven component groups shown.}
```

- **Status**: ✅ CLEAR (caption-only fix). Optional plot regen: add a thin lumi line to the breakdown.

#### V8 — Fig. 8 colors hard to distinguish
- **Cycle-2 finding**: current `plot_paper_final.C` L29 + L1070 — sPHENIX = kAzure+2 (blue), PHENIX original = kPink+5 (pink, RGB 0.8,0.2,0.6), PHENIX corrected = kViolet+1 (purple, RGB 0.6,0.2,1.0). The two PHENIX colors sit too close in color space.
- **Proposed change** (colorblind-friendly, code diff in `plotting/paper/plot_paper_final.C`):

```diff
- const int col[]   = {kAzure + 2, kPink + 5, ...};
+ const int col[]   = {kAzure + 2, kRed - 4, ...};
- const Color_t kCorrColor = kViolet + 1;
+ const Color_t kCorrColor = kOrange + 7;
```

And caption text (L358):

```latex
Comparison of the present isolated prompt photon cross section measurement (blue) with the PHENIX measurement~\cite{PHENIX:2012jgx} (red). The corrected PHENIX measurement (orange), including the bin-width correction and the pseudorapidity-acceptance correction to $|\eta^\gamma|<0.7$, is overlaid. Details of the correction procedure are described in the text. (Folds J-Fig8-parens.)
```

- **Status**: ⚠️ PARTIAL — code change + figure regen needed. **TODO [Virginia-V8]**: apply color change, rerun `plot_paper_final.C`, verify visually.

#### V9 — Fig. 8 stat-error-bar style
- **Cycle-2 finding**: `plot_paper_final.C` L1166–1275 — currently draws the PHENIX and sPHENIX statistical errors as two vertical bars offset by ±0.25 GeV from the bin center. They look like a single fat bar at a glance.
- **Proposed change**: move sPHENIX stat to a thin band around y=1 (mirroring the systematic-band style); keep PHENIX stat as a vertical bar at the bin center (per-point).
- **Code change**: ~15 lines in `plot_paper_final.C` L1224–1268.
- **Caption update**:

```latex
Per-point vertical bars show the PHENIX statistical uncertainty; the inner dark band around unity shows the sPHENIX statistical uncertainty and the outer light band shows the sPHENIX systematic uncertainty.
```

- **Status**: ⚠️ PARTIAL — plot code change needed. **TODO [Virginia-V9]**.

#### V-Run24 — L46 "Run 24" → "RHIC Run 24"
- **Current**: macro `\RunOne` expands to "Run 24".
- **Proposed**: change macro definition (in `defs.sty`) so `\RunOne` → "RHIC Run 24". Single-line edit affects all occurrences.
- **Status**: ✅ CLEAR.

#### V-Bdt-fig — BDT distribution figure
- Overlaps Jamie J-L179. Same plan: new figure showing BDT score data vs MC at one representative pT slice, tight-BDT threshold drawn as a vertical line.
- **Status**: ⚠️ PARTIAL — overlaps with Jamie J-L179. Single new plotting macro covers both.

---

### JUSTIN

#### J-J1 — Split abstract first sentence at "200 GeV"
- **Proposed change**: split into two sentences. **Folded into J-Abs-1 abstract rewrite above.**
- **Status**: ✅ CLEAR.

#### J-J2 — Abstract "sensitive to" → "dependence on"
- **Conflict with Shuonli's "keep" stance**. See J-L20.
- **Status**: 🔴 UNRESOLVED.

#### J-J3 — Intro "particularly sensitive to"
- Same as J-J2. **Status**: 🔴 UNRESOLVED.

#### J-J4 — "Article" + "first sPHENIX measurement"
- **Folded into J-L44.** ✅ CLEAR.

#### J-J5 — L47 "truth-level" applied to data
- **Folded into J-L48.** ✅ CLEAR.

#### J-J6 — "comprises" → "consists of"
- **PDF L57 / main.tex L89**: *"sPHENIX comprises the following subsystems..."*.
- Note: "comprises" used actively is grammatically correct (Garner's, Strunk & White). Justin's "is comprised of" is **debatable**. Alternative: "consists of" is unambiguous.
- **Proposed change**: `comprises` → `consists of`.
- **Status**: ⚠️ PARTIAL — Justin's literal suggestion is grammatically debatable, but "consists of" satisfies the intent. **TODO [Justin-J6]**: confirm Shuonli prefers "consists of" over keeping "comprises".

#### J-J7 — L65 missing "a"
- **Proposed change**: *"The EMCal is **a** scintillating-fiber tungsten sampling calorimeter..."*.
- **Folded into J-J8 below.**
- **Status**: ✅ CLEAR.

#### J-J8 — L66/67 "η-φ two-dimensional projective geometry" + testbeam reference
- **Proposed change** (combines J7 + J8):

```latex
The EMCal is a scintillating-fiber tungsten sampling calorimeter with two-dimensional projective geometry in $\eta$--$\phi$~\cite{Aidala:2020toz, sPHENIX:2017lqb} and approximately $0.025 \times 0.025$ $\eta$-$\phi$ segmentation.
```

- **Verification**: `Aidala:2020toz` confirmed to be the "2-D Projective sPHENIX EMCal Prototype" paper (arxiv 2003.13685). Already cited at current L89.
- **Status**: ✅ CLEAR.

#### J-J9 — L70 "online-calculated" trigger threshold
- **PDF L69–70 / main.tex L103**: *"Events are selected by a single high-ET photon trigger that requires the total energy in any 8 × 8 EMCal tower window to exceed 4 GeV."*
- **Proposed change**: insert "online-calculated".

```latex
Events are selected by a single high-$E_T$ photon trigger that requires the total online-calculated energy in any $8\times 8$ EMCal tower window to exceed 4~\GeV.
```

- **Evidence**: per `selection.tex` L132–142 of the analysis note, bit-30 is an L1 hardware trigger on FPGA primitives. Trigger turn-on plateau at offline-cluster ETγ ≈ 8 GeV confirms the firmware sums are less-calibrated than offline (`reports/trigger_bit30_turnon.md`).
- **Status**: ✅ CLEAR.

#### J-J10 — Eq 2 parens + remove dot
- **Current main.tex L214–221**: ends with `.` after `\frac{N^C}{N^D}`.
- **Proposed change**:

```latex
\begin{equation}
N^{A}_{\text{signal}}
= N^{A}_{\text{raw}}
- N^{B}_{\text{raw}}
\left( \frac{N^{C}_{\text{raw}}}{N^{D}_{\text{raw}}} \right),
\label{eq:purity_noleakage}
\end{equation}
```

(period → comma, ratio wrapped in parens)
- **Status**: ✅ CLEAR.

#### J-J11 — Eq 3 remove dot
- **Current main.tex L225–246**: already ends with `,`. The `\cdot` between B-factor and (C/D)-ratio is what Justin likely means.
- **Proposed change**: remove the `\cdot` (parallel to J10 once ratio is parenthesized).
- **Status**: ⚠️ PARTIAL — interpretation. **TODO [Justin-J11]**: confirm with Justin whether he means the `\cdot` or something else.

#### J-J12 — Fig 4 caption "an Pade" → "a Pade"
- **Current main.tex L266**: *"fitted by an Pad\'e fit"*.
- **Proposed change**:

```latex
\caption{Purity as a function of \etg{} with and without signal leakage correction. The purity with leakage correction is fitted with a Pad\'e function, and the shaded area shows the 68.3\% confidence interval of the fit.}
```

- **Status**: ✅ CLEAR.

#### J-J13 — L358 truth-level in JETPHOX comparison
- **Folded into J-L48 / J-J5 fiducial-definition fix.**
- **Status**: ✅ CLEAR.

#### J-J14 — L372 conclusion "sensitivity to"
- Same as J-J2 / Shuonli's gluon-PDF stance. **Status**: 🔴 UNRESOLVED.

---

## Cycle-1 summary

### Items ready to apply (✅ CLEAR — strong evidence, Shuonli decisions applied)
1. J-L10 lumi consistency (already aligned)
2. J-L18 PHENIX abstract sentence (already removed)
3. J-L20/J-L43/J-L368/J-L372/J-J2/J-J3/J-J14 **gluon-PDF soften to "depends on"** (Shuonli's option (b))
4. J-L34 cross-section hyphenation (sed pass)
5. J-L37 STAR bibtex add + cite
6. J-L44 Article → paper + first sPHENIX
7. J-L48 truth-level definition expansion (folds J5 + J13 + L96 + L113)
8. J-L77 PYTHIA-8 8.307 redundancy (defs.sty macro)
9. J-L92 / L113 stale smearing + truth-iso rewrites
10. J-L96 truth-iso ν+μ clarification
11. J-L112 remove "hard"
12. J-L124 resolution definition
13. J-L128 π⁰/η merger pT (cycle-2 verified: 11 GeV / 46 GeV at Δθ=0.024)
14. J-L132 centre → center
15. J-L236-237 trim redundant text
16. J-L285 escale prose (Shuonli confirmed: this is escale; MBD-Emma is L306)
17. J-L324 HERWIG demotion (cycle-2 corrected: ε_iso 3–6% disagreement IS propagated as ±5.83% syst)
18. J-Pade equation
19. J-L266 MBD-vertex eff direction + numbers
20. J-Fig8-parens
21. J-Refs Collaboration consistency
22. J-NNPDF (bib add only)
23. J-Ichou bib add (optional, available if wanted)
24. J-J8 EMCal 2D projective + cite
25. J-J9 online-calculated (with EMCal-expert sign-off recommendation)
26. J-J10 Eq 2 parens + comma
27. J-J12 a Pade

### Items requiring decision/external input (🔴 UNRESOLVED)
- **J-L25** — ISR/FSR phrasing (Shuonli: hold for expert)
- **J-L306** — vague MBD-eff sentence; **TODO: get PPG10 bibtex key from Shuonli; text drafted but `\cite{PPG10}` is placeholder**
- **J-L67.window** — L1 photon-trigger window size (PDF says 8×8, code/memory says 4×4 — verify with sPHENIX trigger team)
- **J-J6** — "comprises" → "consists of" (Shuonli confirm exact wording)

### Items partially supported (⚠️ PARTIAL — apply with TODO marker)
- **J-L30** — ref [1] inside period (verify compiled rendering)
- **J-L67** — detector expansion (drafted; needs L1-window-size + 5% in-situ EMCal expert sign-off)
- **J-L92.fig** — π⁰/η mass-peak figure (cycle-2 agent still locating ROOT files)
- **J-L138** — upstream-muon physical insight (PPG12-specific verify pending)
- **J-L144** — NCB feature list (response-letter only)
- **J-L179** — BDT distribution figure (new plotting macro)
- **J-Fig1** — pure-background histogram (Shuonli's TODO: implement in plotting macro)
- **J-Fig3** — schematic + 2D heatmap (new plotting macro)
- **J-Fig6** — vertical lines (plotting macro)
- **J-L338** — middle + bottom panel (depends on NNPDF4.0)
- **J-NNPDF** — JETPHOX rerun (separate work)
- **J-Fig8.stat** — stat band cosmetic
- **J-J11** — eq 3 dot interpretation (confirm with Justin)

### Items still in flight
- All Virginia items (V1–V9, V-Run24, V-Bdt-fig) — ⏳ CYCLE-2 PENDING (background agent)
- π⁰/η mass-peak files locator — ⏳ CYCLE-2 PENDING (background agent)
- Critic audit of cycle-1/2 patch — ⏳ CYCLE-2 PENDING (background agent)

---

## Cycle 2 result summary

**Cycle 2 returned**: 5 agents (EMCal resolution + L1 trigger ✓; π⁰/η merger pT ✓; mass-peak files ✓; HERWIG+Emma+Ichou ✓; critic audit ✓) + Virginia (originally cycle-1 background, returned during cycle 2) ✓.

**Shuonli's cycle-2 confirmations (2026-05-14)**:
- L285 = escale shift (confirmed; matches J-L285 in patch).
- Publication-ready writing rule added — all proposed LaTeX text must be jargon-free and journal-grade.

**Critic-flagged blockers fixed in cycle 2**:
- `Ichou:2010er` → existing key `Ichou:2010wc` (no duplicate entry).
- `Harrison:2003sb` → existing key `HARRISON2003235`.
- `sPHENIX:2022emcal_testbeam` (not in paper's `references.bib`) → fall back to existing `Aidala:2020toz` + `sPHENIX:2017lqb`.
- L67 4×4 vs J-J9 8×8 trigger window — pulled the specific number out of the proposed prose; awaiting John-Haggerty / L1 expert sign-off.
- L67 0.024 vs J-J8 0.025 granularity — reconciled to 0.024 (matches `reconstruction.tex` L4).
- L372 lumi inconsistency: paper L371/372 still says 64.37; add to lumi-fix scope.
- J-Abs-1 abstract uses "depends on" gluon-PDF wording (per Shuonli's option (b)), matching J-L20/L43/L368/L372.
- Mass-peak figure (J-L92.fig) — no PPG12 internal files; recommend defer to revision 2 or borrow EMCal-WG plot.

## Final cycle-2 patch summary

### Items ready to apply (✅ CLEAR — strong evidence + Shuonli decisions applied + critic-audited)
Jamie: L10 (uniform 1-dp `64.4⁺⁵·⁹₋⁴·³` everywhere), L18 (already done), L20/L43/L368/L372 gluon "depends on", L34, L44 (folds J-J4), L77, L92 (smearing prose simplified), L96 (drop false "excluding neutrinos"), L112, L124, L128 (π⁰/η merger with pendingTODO), L132, L138 (beam-pipe interactions), L236-237, L266, L285 (escale prose), L324 (HERWIG with 3-6% disagreement), L325 envelope (add eres co-dominance), Fig1 (relabel "Inclusive MC" → "Inclusive jet MC"), Pade equation, Fig8-parens, Refs Coll., NNPDF bib add, Ichou cite, J-J8 EMCal 2D projective, J-J10 Eq 2 parens, J-J12 a Pade, Abstract redraft.
Virginia: V1 (direct→prompt), V2 (sPHENIX:2025dET bib add), V3 (simulation-sample wording + pileup-frac weighting), V4 (kinematic-range wording), V5 (iter 3,4), V6 (third treatment = unfolding iter max, not eres; eres moved to two-sided), V7 (lumi caption), V-Run24 (RHIC Run 24).
**NEW (audit-flagged)**: Tight/non-tight BDT cuts paper-vs-note drift — three numbers update at main.tex L179–181 (intercepts 0.833→0.8156 and 0.667→0.6844; slopes ∓0.00333→∓0.00156). ✅ CLEAR.

### Items requiring external input (🔴 UNRESOLVED — blocked until Shuonli or expert)
- **J-L25** ISR/FSR (Shuonli: hold for expert).
- **J-L37** STAR 0912.3838 — hold until STAR-overlay plotting lands in `final_phenix.pdf`; citation + bib + intro update lands in same commit.
- **J-L48 + J-J5 + J-J13 (truth-level fiducial bundle, minus L96)** — hold pending expert consultation on canonical truth-level fiducial wording (stable vs PYTHIA-status-1, JETPHOX equivalence statement, intro-vs-jargon framing for data). J-L96 has been cleared with a minimal truthful fix (drop false "excluding neutrinos"). Remaining items will commit as one unified paragraph once resolved.
- **J-L306 / J-Emma** PPG10 bibtex key needed.
- **J-L67.window** 4×4 vs 8×8 trigger primitive (verify with John Haggerty).
- **J-L67.online** "online-calculated" qualifier (sPHENIX L1 firmware team).
- **J-J6** "comprises" → "consists of" (Shuonli confirm exact wording).
- **J-1.x-range** "0.12 ≲ x ≲ 0.32" in abstract (Shuonli decide if it stays).
- **J-J11** Eq 3 dot interpretation (confirm with Justin).

### Items requiring follow-up plot/code work (⚠️ PARTIAL — separate work item)
- **J-L92.fig** π⁰/η mass-peak figure: defer to revision 2 (no PPG12 internal files; recommend new diphoton study or borrow EMCal-WG plot).
- **J-L92.thresholds** σ_extra(pT) thresholds 17 / 20 GeV — verify by direct extraction.
- **J-L179 / V-Bdt-fig** BDT distribution figure (new plotting macro).
- **J-Fig3** schematic + 2D data heatmap (new plotting macro).
- **J-Fig6** SetErrorX(0) cosmetic.
- **J-Fig8.stat / V8 / V9** Fig. 8 recolor + sPHENIX stat band (code change in `plot_paper_final.C`).
- **J-L338** middle + bottom panel discussion (depends on NNPDF4.0).
- **J-NNPDF** JETPHOX rerun with NNPDF4.0 (~3-hour wallclock, separate work).
- **J-L30** ref [1] inside period (verify post-compile rendering).
- **J-L144** NCB feature list — response-letter only.
- **J-L36** CT18NLO Hessian band — DEFERRED to revision 2.

---

## Response-letter format (per Shuonli, 2026-05-14)

For each reviewer comment in cycle 3.5, use this format:

```
> [reviewer's actual comment, verbatim]

[Our response: "Addressed: ..." OR "Held: ..." OR "Not changed because ..."]
```

One section per reviewer (Jamie, Virginia, Justin). Save to `reports/response_letters.md`.

## UNRESOLVED marker policy in main.tex (per Shuonli, 2026-05-14)

For every 🔴 UNRESOLVED item with a specific reviewer comment, insert a `\pendingTODO{...}` next to the relevant text in `main.tex` containing the **verbatim** reviewer comment. The macro renders in red so the unresolved items are visible in the compiled PDF.

Format:
```latex
\pendingTODO{Reviewer X (Jamie/Virginia/Justin): "verbatim comment"}
```

For audit-only items without a specific reviewer comment (e.g., x-range in abstract, J-L67.window), use a generic note:
```latex
\pendingTODO{TODO: brief description of what's blocked}
```

The 9 UNRESOLVED items + verbatim comments:

| # | Item | main.tex location | Verbatim reviewer comment |
|---|---|---|---|
| 1 | J-L25 ISR/FSR | L70 (prompt-photon definition) | Jamie: *"in field theory, I thought ISR and FSR were indistinguishable and can also interfere. Thus, are fragmentation photons somehow identifiable as exclusively FSR?"* |
| 2 | J-L37 STAR overlay | L75 (PHENIX/RHIC reference) | Jamie: *"I believe PHENIX has measurements in pp, dAu, AuAu as paralleled by LHC. Also, STAR has a measurement (https://arxiv.org/pdf/0912.3838), and not sure if there is something more recent. Please double check."* |
| 3 | J-L48 truth-level (intro) | L80 | Jamie: *"should we define 'truth level' here, a la the detailed definition used by ATLAS?"* + Justin J5: *"truth-level bad here if you are just describing what's done to the data"* |
| 4 | J-L48 truth-level (Sec III) | L113 | Jamie L96: *"see above point regarding truth level. Anything about muons here?"* |
| 5 | J-J13 truth-level (JETPHOX) | L340 | Justin: *"see previous comment about use of truth level, I don't like this if just meant as jargon — what does truth level mean for data?"* |
| 6 | J-L306 PPG10 / MBD-eff vague | L321 | Jamie: *"this sentence is quite vague"* + Shuonli note: check PPG10 ref |
| 7 | J-L67.window 4×4 vs 8×8 | L103 (trigger sentence) | (audit-flagged, no specific reviewer comment) |
| 8 | J-L67.online qualifier | L103 (trigger sentence) | Justin: *"[requires] total online-calculated energy to be above 4 GeV"* |
| 9 | J-J6 comprises/consists | L89 (detector) | Justin: *"change comprises to 'is comprised of'"* |
| 10 | J-1.x-range in abstract | L60–62 (abstract) | (no specific comment, optional) |
| 11 | J-J11 Eq 3 dot | L225 (eq:purity_leak) | Justin: *"L225 equation 3 remove the dot it's not needed"* |

These markers will be removed once the corresponding item is resolved (expert input, PPG10 ref, etc.).

## Cycle 3 plan (pending Shuonli's review of this cycle-2 patch)

1. Shuonli reviews this patch and provides:
   - Decisions on the 10 UNRESOLVED items (or marks them as defer-to-revision-2).
   - PPG10 bibtex key (for J-L306).
   - Expert sign-off (or wait) on EMCal resolution / trigger window / L1 calibration.
2. Apply all approved CLEAR items to `PPG12-Paper/main.tex` and `references.bib` in a single commit.
3. Build LaTeX: `cd PPG12-Paper && pdflatex main && bibtex main && pdflatex main && pdflatex main`.
4. Inspect compile log + rendered PDF for any new errors (broken citations, undefined macros, etc.).
5. Commit the changelog into `reports/paper_review_patch_cycle1.md` (this file).
6. If 3 cycles are exhausted and UNRESOLVED items remain: STOP, report exact remaining blockers to Shuonli.
