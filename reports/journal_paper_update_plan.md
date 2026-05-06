# PPG12 Conference Note → Journal Manuscript: Proposed Update

**Date:** 2026-05-05
**Source:** `PPG12-Paper/main.tex` (sPH-CONF-JET-2025-02 conf-note draft, 280 lines, ~2k words of prose)
**Target:** journal manuscript draft, **Phys. Rev. D** (settled — see §1)
**Status:** **NOTHING APPLIED.** This document collects all proposed text changes for the user to review. `PPG12-Paper/main.tex` is unchanged.

> Decisions: journal target = **PRD** (§1), section structure = **Option A / single §3 Analysis Procedure** (§2). Placeholders in §5.1 still pending calorimeter / lumi / vertex-reweight input.

---

## Changelog
- Pass 1 (2026-05-05): Applied reviewer-pass-1 fixes. CRITICALs: corrected the bit-10 coincidence definition to "MBD\,N\&S\,$\geq 1$" (one PMT each side, not two; this is a fix to how bit 10 is *described*, independent of its role in the analysis — see Pass 5), EMCal granularity to $0.024\times 0.024$, \etone definition to four-tower sum, BDT-training-bin count to 2, NPB-cut efficiency claim, low-\etg{} purity bin to $[10,12]$ (or $[12,14]$ TBD from ROOT), Padé[1/1] as nominal with erf as alternative, response-matrix-prior reweight input as iter-1 D'Agostini unfold ratio, $\varepsilon_\mathrm{trig}$ applied at fill time outside $\mathcal{E}$, unfolding envelope as iter-3/iter-4 only, DI-fraction systematic as $f_\mathrm{best}$ from BDT-score $\chi^2$ fit. WARNINGs: pPb refs restored, $L_{0\mathrm{mrad}}/L_{1.5\mathrm{mrad}}$ inlined ($47.21$, $17.16$ pb${}^{-1}$), MBD timing-veto cut to follow live config, noise-overlay sentence weakened, \etg{} reach reframed via $x_T$, purity systematic enumerated as five sources, retired trigger-plateau syst noted, tower-acceptance cross-check noted, PHENIX rescaling source corrected (truth simulation, not "$\eta$-flat assumption"). MISSINGs: bib keys to add (`ATLAS:2017topoclusters`, `ATLAS:2019isophoton13TeV`), new §5.1 entries for 1.5 mrad truth-vertex closure, bit-26-vs-bit-30, photon-trigger-plateau correction. Bloat: removed already-settled R=0.4-vs-R=0.3 TODO and consolidated §5.5 against inline TODOs. Style: stripped semicolons in proposed paragraphs.
- Pass 2 (2026-05-05): Reviewer signed off (PASS, no Task #4). Tightened §5.2 NPB description: the ${\approx}99\%/{\sim}99\%$ retention/purity numbers are now anchored to the analysis-note's $22<\etg<28~\GeV$ representative bin (matching `reconstruction.tex:162`), with a TODO to quote the lowest-\etg{} bin once that figure is regenerated. No other changes; reviewer's other Pass-2 note (lumi $64.4$ vs.\ $64.37$ pb${}^{-1}$ rounding) was acknowledged as acceptable.
- Pass 3 (2026-05-05): Locked in PRD as the journal target. The user's original "PRC/PRD" framing was an implicit ack that either is acceptable, and the precedent analysis converged on PRD. The §1 heading and prose are updated, and the trailing "Open question for user" line is removed. PRC remains a valid override if explicitly requested.
- Pass 4 (2026-05-05): Locked in Option B (split §3 into Data / MC / Photon-Reco-ID / Analysis-Procedure top-level sections) as the structural decision. The proposed prose was already written for Option B, so this is a labelling-only change to §2: the heading is now "Section structure: split (Option B, decided)" and the user-confirm/override line is removed. Reverting to flat §3 (Option A) is still cheap if requested later (label rename only).
- Pass 5 (2026-05-05): Corrected §3.2 trigger description per the user's note that the new analysis no longer uses the MBD trigger as an event-selection trigger. The "Two triggers are used in this analysis" framing is replaced. Bit~30 (Photon 4 GeV) is now described as the sole sample-selection trigger, with bit~10 (MBD\,N\&S\,$\geq 1$) used only as (i) the unbiased denominator for the photon-trigger turn-on calibration and (ii) the per-run live-scaler input to the luminosity calculation. A sentence is added that ties the $|\vz|<60$~cm cut and the all-$z$ lumi together via the MBD-vertex efficiency (consistent with `analysis.tex:303`).
- Pass 6 (2026-05-05): Applied deep-reviewer audit against the live `PPG12-analysis-note/`. Categories of changes:
  - **Lumi rounding unified to 64.37 pb⁻¹** (matches `analysis.tex:303`, `selection.tex:170`, `conclusion.tex:1`). Asymmetric uncertainty written as $^{+5.88}_{-4.34}$~pb⁻¹ everywhere (abstract, §3.1, §6.5, §9, TL;DR).
  - **Subscript fix `\sigma_{\mathrm{MB}}`→`\sigma_{\mathrm{MBD}}`** in §3.1, §7 lumi paragraph, §5.4 bib list, and TL;DR (matches `analysis.tex:303`, `systematics.tex:218`).
  - **Theory attribution "Gordon and Vogelsang"→"Vogelsang"** in abstract, §1, §6 intent line, §8 prose, §8 figure caption, and §9 summary (matches `results.tex:26`'s "Werner Vogelsang" only).
  - **MBD timing-spread veto removed from §3.2.** No `mbd_avg_sigma_max` exists in `config_bdt_nom*.yaml`. Replaced with a sentence stating the analysis has no event-level pileup veto and that DI is handled MC-side (matches `selection.tex:166`).
  - **§4.4 DI fraction derivation** changed from "MBD-timing-derived event-level pileup probability" → "Poisson decomposition of the MBD-coincidence rate" (matches `double_interaction_consolidated.tex:18-58`).
  - **§4.4 truth-vertex reweight** caveats added: photon10-only training subset (per memory `project_truth_vertex_reweight.md`) and 1.5 mrad fit non-convergence flagged inline.
  - **§4.4 truth-reco match** wording changed from vague "EMCal cluster–truth association" → "track-ID-based via `CaloRawClusterEval`" (matches `analysis.tex:221`).
  - **§4.4 PYTHIA framing** corrected: "separately for signal and background productions" → "both `HardQCD:all = on` and `PromptPhoton:all = on` enabled in every production, partitioned by truth-level event filter" (matches `selection.tex:32`).
  - **§4.4 MC-sample-table TODO** updated to mirror the analysis-note pair `table:MCdataset` (eight SI/DI rows including `jet40`/`jet40_double`) and `table:stitching` (truth leading-pT windows).
  - **§5.1 cluster-reconstruction prose** trimmed: dropped unsupported "0.5 GeV total / 0.2 GeV seed / 3×3 local-maximum" thresholds, kept only the 70 MeV tower threshold + cluster-splitting algorithm (matches `reconstruction.tex:8`).
  - **§5.1 in-situ resolution** "${\sim}6\%$ at $E\sim 10$~\GeV{}" → "${\sim}5\%$ across the analysis range" (matches `analysis.tex:246`); also corrected the §2 detector-text claim that data was smeared to MC (it is the other way around).
  - **§5.2 pre-selection list** dropped \wetacogx (it is a BDT input feature, not a rectangular pre-selection cut; `tab:photonid` lists only E11/E33, et1, E32/E35, NPB).
  - **§5.2 NPB BDT framework** "built with TMVA" → "trained in XGBoost and exported to TMVA" (matches `reconstruction.tex:148`).
  - **§5.3 et1 definition** "four towers surrounding the CoG" → "$2\times 2$ block including the CoG tower" (matches `reconstruction.tex:78-82`).
  - **§5.3 BDT training** dropped the unsupported "two ETγ bins ($8<\etg<15$ and $15<\etg<35$)" claim (analysis-note is silent on the binning); added the live $\sim 750$-tree, depth-5, 50/50 train/test parameters from `reconstruction.tex:325`.
  - **§6.1 isolation prose** corrected: dropped the unsupported "60 MeV tower seed" and "tower-adjacent exclusion in $\eta,\phi$" prescriptions; replaced with "seed-significance threshold" and "candidate's own reconstructed-cluster $E_T$ subtracted from the cone sum" (matches `reconstruction.tex:398`).
  - **§6.1 iso-fudge wording** added "to MC only; data are unmodified" (matches `analysis.tex:163`).
  - **§6.2 bootstrap-trial inline TODO** dropped — `analysis.tex:113` already states 10,000 trials.
  - **§6.2 figure caption** "fitted with an error function" → "fitted with an order-$[1/1]$ Pad\'e approximant" (the erf form is the alternative, not the nominal — matches `analysis.tex:134`, `systematics.tex:99`).
  - **§6.3 unfolding-iteration criterion** "minimise the quadrature sum of iter-to-iter change and bootstrap stat" → "successive-iteration changes fall below per-bin statistical uncertainty while iteration-bias remains subdominant" (matches `analysis.tex:259`).
  - **§7 energy-scale** "up to ${\sim}30\%$ above 26 GeV" → "$-22.6\%/+29.3\%$ at the high-\etg{} end" (matches `systematics.tex:28`).
  - **§7 energy-resolution** "$5$--$13\%$" → "$-17.1\%/+6.3\%$" (matches `systematics.tex:43`).
  - **§7 purity** restructured from five contributions → six: tight ID, non-tight ID, non-iso boundary, fit functional form, fit-CI ($\pm 1\sigma$ Pad\'e parameter posterior), and ABCD MC closure $C_\mathrm{MC}$ as its own enumerated entry (matches `systematics.tex:6-11`, sub-section 5.4 + 5.5 + 5.6 split).
  - **§7 efficiency** rewritten: dropped tight-threshold and BDT-model variations (not in the live efficiency budget per `systematics.tex:135-142` and memory `project_syst_redo_needed.md`); kept only the iso-pedestal $b\to 0$ variation reaching $\pm 4.5\%$, plus the $\pm 6\%$ MBD-vertex-eff placeholder split out as its own envelope.
  - **§7 truth-vertex-reweight closure** added as a new cross-check paragraph (matches `systematics.tex:213-214` and `double_interaction_consolidated.tex:580-583`).
  - **§7 totals** "$-30\%$ to $+35\%$" → "$-30.4\%/+35.1\%$" (matches `systematics.tex:236`).
  - **§8 results / §9 summary** lumi 64.4 → 64.37; total uncertainty range tightened to $-30.4\%/+35.1\%$.
  - **§8 $x_T$-scaling** TODO marked explicitly as new content not present in the analysis-note, requiring collaboration sign-off.
  - **§5.5 TODO list** pruned: removed entries for MBD-timing-veto cut, run-range upper bound, and bootstrap trial count (now resolved). Added an entry tracking the 1.5 mrad truth-vertex closure and an MC-sample-table pointer.
- Pass 6b (2026-05-05): Applied three deep-reviewer INFO follow-ups. §7 tower-mask wording aligned to "$\phi$-symmetry" (was "hard-dead"). §7 totals tightened to precise $-30.4\%/+35.1\%$ with `tab:syst_summary` reference (removing redundant soft phrasing). §8 PHENIX caption "assuming $\eta$-flat" replaced with "from the prompt-photon truth distribution in the \pythia{} signal simulation" to match the Pass-1-fixed prose.
- Pass 7 (2026-05-06): Flipped the structure decision from Option B (split) to **Option A (single §3 Analysis Procedure)** per user request. (i) §2 lock-in updated to Option A. (ii) Top-of-doc decisions box and TL;DR section list updated. (iii) §4.1 intro-outline paragraph rewritten to enumerate §3.1–§3.7 as subsections of a single §3, with §4 Systematics, §5 Results, §6 Summary as the remaining top-level sections. (iv) Plan §4.3–§4.6 consolidated into a single §4.3 "§3 Analysis Procedure" block. Each former top-level `\section{...}` (Data Sample / Monte Carlo / Photon Reco-ID) is demoted to `\subsection{...}` of §3. The former `\subsection{Double-interaction handling}` becomes a `\paragraph{...}` inside §3.2; the former §5 sub-blocks (Cluster reconstruction, Pre-selection, Photon-ID BDT) become `\paragraph{...}` blocks inside §3.3. The §6 Isolation + Signal Extraction merge into §3.4 Signal Extraction with iso content as a `\paragraph{...}` lead-in. (v) Plan §4.7 → §4.4 (Systematics), §4.8 → §4.5 (Results), §4.9 → §4.6 (Summary), §4.10 → §4.7 (trailing matter). (vi) Section labels updated: `sec:dataset` → `ssec:eventsel`, `sec:mc` → `ssec:mc`, `sec:reco` → `ssec:reco`. The single §3 keeps `\label{sec:analysis}` as the parent.
- Pass 8 (2026-05-06): Applied the user's seven inline comments + Wave-2-flagged fixes. Categories of changes:
  - Comment 1 (§4.1): dropped dir/frag motivation figure block entirely (PYTHIA-LO unreliable per user). Removed the commented `\begin{figure}` block, the `\pendingTODO{regenerate pT_Spectrum_IsoET_4GeV_Cut.pdf}`, and the 75–85% direct-photon retention sentence that referenced it.
  - Comment 2 (§4.2): reverted detector section to the conf-note original prose (sPHENIX-wide convention text). The Pass-6 expansion (coordinate system, MBD detail, in-situ resolution paragraph, etc.) is replaced with the unchanged conf-note paragraphs, with the standard EMCal/HCal calorimeter detector citations from the conf-note (`Aidala:2020toz`, `sPHENIX:2017lqb`) retained.
  - Comment 3 (§4.3.1): trimmed §3.1 from 13 sentences (3.3× conf-note) to 5 sentences at conf-note granularity. Removed all trigger bit numbers (no "bit 10" / "bit 30" / "Level-1 bit"). Removed internal-note cites (`ppg09IAN`, sPHENIX vdM-scan note). Dropped the dataset / triggers / vertex-selection paragraph splits and the dataset-summary table TODO.
  - Comment 4 (§4.3.2): replaced the multi-sentence truth-vertex reweight prose with the Wave-1E single sentence ("a per-event truth-vertex reweight is applied to the MC so that the MC reconstructed-vertex distribution matches the data after the standard analysis cuts"). Dropped the χ²/ndf=14.79 pendingTODO from §3.2 entirely.
  - Comment 5 (§4.3.3): verified §3.3 against conf-note granularity (1.6× per ledger, within tolerance). Inserted two shower-shape distribution figures (`dis_mixed_weta_cogx_eta0_pt2_cut1.pdf` and `dis_mixed_e32_to_e35_eta0_pt2_cut1.pdf` from `plotting/figures/showershape_1p5mrad/`) with caption + 2-sentence motivation paragraph illustrating signal/background separation at the NPB-pre-selection level.
  - Comment 6 (§4.3.4): inserted Wave-1A's topo-cluster paragraph in §3.4 Signal Extraction, citing `ATLAS:2017topoclusters` and naming the seed–neighbour–local-maximum-splitting algorithm with the live (4σ, 2σ, 1σ) significance steps drawn from `coresoftware/offline/packages/CaloReco/RawClusterBuilderTopo.{h,cc}` and `Fun4All_run_sim.C`.
  - Comment 7 (§4.5): updated PDF-panel description to four PDFs (CT18NLO, NNPDF3.1, CTEQ6.6, MSHT20NLO) per the live `plot_final_selection.C` third pad. Stated CT18NLO + truth-isolation requirement for both \jetphox{} and the Vogelsang calculation. Dropped the χ²/ndf TODO from the §5 ratio statement and the \jetphox{}-CT14NLO-vs-CT18NLO diff TODO. Dropped the $x_T$-scaling paragraph entirely. Updated the abstract and §6 summary to match (CT14NLO → CT18NLO).
  - User addition (cross-cutting): the seven inline `comment: ...` lines are removed from the plan body now that they are addressed.
- Pass 8b (2026-05-06): Five fixes from Wave 4 deep audit. (i) [CRITICAL] §3.1 EMCal granularity regression $0.025\times 0.025$ → $0.024\times 0.024$ (Pass 1 fix lost during Pass 8's conf-note revert, now restored). (ii) [WARNING] §3.4 noise-values phrasing: rewrote to clarify $5.3/35.1/68.4$ MeV are three standard deviations of per-cell pedestal noise (not $\sigma_\mathrm{noise}$), paired with $4\times/2\times/1\times$ significance multipliers from `RawClusterBuilderTopo`. (iii) [STYLE / user-decision option (b)] §3.3 light trim: dropped BDT training-parameters sentence, consolidated tight + non-tight parametric BDT equations into one labelled tight-cut equation plus inline non-tight band, with §3.3 granularity now ~1.8× vs the prior 2.4× while keeping the user-requested shower-shape figures and weta equation. (iv) [STYLE] Pass-8 changelog cosmetic: corrected the EMCal test-beam-paper retention claim. (v) [INFO] §5.4 added an `ATLAS:2016fta` status row.
- Pass 9 (2026-05-06): Aggressive trim + NPB→NCB rename + L346 follow-up. Categories of changes:
  - **NPB → "non-collision background (NCB)"**: 16 forward-looking body occurrences renamed across 12 plan lines (introduction at §3.3 first use defines the term, subsequent uses NCB). Changelog history left as "NPB" verbatim. §4 \paragraph{Background-rejection cut (NPB).} → \paragraph{Background-rejection cut (NCB).}. Definition extended at §3.3 first use to read "non-collision backgrounds (NCB, predominantly beam-induced clusters with anomalously broad $\eta$ profiles and out-of-time EMCal--MBD timing)" per `reconstruction.tex` evidence that NCB is tagged via EMCal--MBD timing, which is the signature of beam-induced (= non-collision) clusters. Label keys (`fig:npb_purity_yield`) and code identifiers untouched.
  - **§3.4 iso lead-in (per L346 inline comment)**: dropped explicit topo-cluster noise scales ($5.3/35.1/68.4$ MeV) and seed/grow/close significance multipliers ($4\times/2\times/1\times$). Kept the topo-cluster algorithm citation (`ATLAS:2017topoclusters`) and the EMCal/IHCal/OHCal three-layer adaptation. Also dropped the redundant "use cells from all three calorimeter layers and span the full $|\eta|<1.1$ central region" sentence since S2 already names the three-layer stack.
  - **§3.4 iso fudge (per L346)**: dropped explicit ($a=1.2$, $b=0.1$ \GeV) values from body. Replaced with "a small empirical correction is applied to the simulated $\isoET$" + systematic pointer to §4. The §4 \paragraph{Efficiency corrections} keeps the (a, b) values for the systematic-variation description.
  - **§3.4 purity Padé-vs-erf (per L346)**: dropped the body sentence comparing Padé and erf forms. Body now reads "fitted with a smooth functional form, and the fitted value at each bin centre is used to correct the residual background." Caption simplified to "fitted with a smooth functional form, with the shaded band showing the 68.3\% confidence interval of the fit." The §4 \paragraph{Purity} retains the Padé/erf comparison and ${\sim}0.3\%$ figure as a systematic-source detail.
  - **Aggressive trim (Wave-1A ledger)**: §4 intro paragraph dropped the redundant 2nd-sentence enumeration of source categories (already stated as \paragraph headers below). §4 photon-energy-resolution trimmed: $0\%/8\%/$twice-nominal trinity replaced with "varied between zero and twice the nominal value". §4 purity collapsed from six numbered contributions to three quadrature-summed groups (thresholds, fit form + posterior shift, MC closure), 7→2 sentences. §4 unfolding dropped the iter-1-not-converged footnote sentence (4→3 sentences). §4 double-interactions dropped the explicit $f_{\mathrm{DI},\mathrm{alt}}=0.290/0$ values, kept the data-driven $\chi^2$ derivation. §3.5 unfolding collapsed the iter-1 D'Agostini ratio + Padé[2/2] truth-prior reweight construction into a single sentence ("the response-matrix truth prior is reweighted to bring it into agreement with the data \etg{} spectrum"), removing 2 implementation-detail sentences. §3.6 efficiency tightened C1 sentence to drop the all-z-truth-vertex-denominator detail. §3.4 dropped the no-leakage-equation lead-in's explicit $N^B_{\mathrm{bkg}}/N^A_{\mathrm{bkg}}=N^D_{\mathrm{bkg}}/N^C_{\mathrm{bkg}}$ ratio, "Equation~\ref{eq:purity_leak} is solved iteratively." sentence, the 10,000-trial-bootstrap clause, the unfolding-buffer-underflow-bin parenthetical, and the explicit $0.51/0.72/0.90\pm$values (purity now stated as "$\mathcal{P}\sim 0.5$ at the lowest \etg{} bin to $\mathcal{P}\sim 0.9$ above $26~\GeV$"). §3.4 residual-correlation sentence shortened to drop "explicitly", "in the inclusive-jet simulation", "of the ABCD method". §3.3 cluster-reconstruction paragraph dropped the forward-pointing "The full \etg-binned values are reported in Sec.~\ref{sec:results}" sentence.
  - **L346 inline comment**: removed (resolved). Updated the §3.4 intent paragraph above the equations to describe the trimmed prose accurately (no explicit topo-cluster numerical thresholds, no MC iso-fudge $(a,b)$, Padé/erf moved to §4).

---

## TL;DR

The conference note is now obsolete in **every** quantitative respect:
- lumi 16.6 → 64.37 pb⁻¹ (3.9× larger)
- ETγ range 10–26 → 12–32 GeV
- vertex cut |z|<30 → |z|<60 cm
- isolation tower R=0.3 → topo R=0.4 with parametric `IsoET < 0.490 + 0.0370·ET^reco`
- photon ID: 5 rectangular shower-shape cuts → 11-feature BDT (`BDT > 0.8333 − 0.00333·ETγ`) + NCB BDT
- new physics content: full GEANT DI MC blending (f_DI = 22.4% / 7.9%), truth-vertex reweighting, all-z MBD-vertex efficiency, asymmetric +9.13/−6.75% lumi systematic from σ_MBD = 25.2 +2.3/−1.7 mb.

In addition the conf-note is missing structural elements expected by PRD:
- explicit shower-shape variable definitions (formula for `weta`)
- table of MC samples with p̂_T,min and event counts
- bin-by-bin systematic uncertainty table
- bin-by-bin cross-section value table (for HEPData)
- expanded detector subsection with calibration sources
- stronger consistency statement (mean ratio + χ²/ndf, not just "consistent")

The proposed update keeps the existing prose where possible and rewrites everything that's been superseded. The content below is grouped into:

1. Journal target (decided: PRD)
2. Section structure (decided: Option A / single §3 Analysis Procedure)
3. Header/abstract: side-by-side proposed change
4. Section-by-section proposed changes (full LaTeX inline)
5. Open items requiring user/collaboration input (placeholders + figure regeneration + bib + tables)
6. Style/voice fixes

---

## 1. Journal target: PRD (decided)

The journal target is **Phys. Rev. D**, settled per the user's original "PRC/PRD" framing combined with the precedent analysis below.

**PHENIX:2012jgx** (the cited prior measurement at the same √s = 200 GeV, same observable) was published in **Phys. Rev. D 86, 072008 (2012)**. ATLAS/CMS inclusive isolated-photon pp papers predominantly went to PRD (ATLAS:2013sdn = PRD 89, 052004; ATLAS:2011ezy similar). The conf-note text already foregrounds NLO pQCD comparisons and gluon-PDF sensitivity — both core PRD topics.

A PRC framing would be consistent with sPHENIX's heavy-ion mission, but until an Au+Au companion paper exists the pp-only result reads more naturally as a QCD/PDF measurement. A forward-looking paragraph in the introduction motivates the baseline role for upcoming sPHENIX heavy-ion photon measurements without requiring a PRC framing.

---

## 2. Section structure: single §3 Analysis Procedure (Option A, decided)

The conf-note's flat §3 layout is preserved. The journal paper has six top-level sections, mirroring the conf-note exactly:

- §1 Introduction
- §2 sPHENIX detector
- §3 Analysis Procedure
  - §3.1 Event Selection
  - §3.2 Monte Carlo Simulations (with double-interaction handling as an inline paragraph)
  - §3.3 Photon Reconstruction and Identification (cluster reconstruction, pre-selection / NCB BDT, photon-ID BDT as inline paragraphs)
  - §3.4 Signal Extraction (isolation + ABCD background subtraction)
  - §3.5 Unfolding
  - §3.6 Efficiency Correction
  - §3.7 Cross-Section Determination
- §4 Systematic Uncertainties
- §5 Results
- §6 Summary

The proposed text in §4 below uses LaTeX header levels `\section` / `\subsection` / `\paragraph` to match this layout. The single-section structure has lower diff churn against the conf-note, lets line-by-line comparison stay clean, and is consistent with prior sPHENIX paper drafts.

---

## 3. Header / abstract proposed changes

### 3.1 Header (lines 1–33)

**WAS:**
```latex
\documentclass[type=note,linenumbering=off]{sphenix}
\usepackage{amsmath}
\usepackage{lipsum}
\usepackage{cleveref}
\usepackage{mathrsfs}
\usepackage{xspace}
\crefrangeformat{figure}{Figs. #3#1#4--#5#2#6}

\usepackage{defs}

\def\GeV{...}
\def\MeV{...}

\begin{document}

\title{
{\bf{sPHENIX Paper}}
\ \\
\ \\
Measurement of isolated prompt photons in \pp{} collisions at $\sqrt{s}$ = 200~\GeV{} with the sPHENIX detector
}
\author{sPHENIX Collaboration}
\date{\today}

\doctag{sPH-CONF-JET-2025-02}
\docversion{1}
```

**PROPOSED:**
```latex
% sphenix.cls accepts type=note|proceedings|report — there is no `paper` option.
% TODO: switch to revtex4-2 with [aps,prd] for actual journal submission.
\documentclass[type=report,linenumbering=on]{sphenix}
\usepackage{amsmath}
\usepackage{amssymb}     % for \lesssim, \gtrsim
\usepackage{cleveref}
\usepackage{mathrsfs}
\usepackage{xspace}
\usepackage{xcolor}      % for \textcolor TODO markers; remove for final submission
\usepackage{booktabs}    % for systematics + cross-section tables
\crefrangeformat{figure}{Figs. #3#1#4--#5#2#6}

\usepackage{defs}

% TODO: remove \pendingTODO definition once all placeholders are resolved
\newcommand{\pendingTODO}[1]{\textcolor{red}{[\textbf{TODO:} #1]}}

\def\GeV{...}    % unchanged
\def\MeV{...}    % unchanged

\begin{document}

\title{Measurement of isolated prompt-photon production in \pp{} collisions at $\sqrt{s}=200$~\GeV{} with the sPHENIX detector}
\author{sPHENIX Collaboration}
\date{\today}

\doctag{\pendingTODO{paper-track tag, e.g.\ sPH-PAPER-2026-XXX}}
\docversion{0.1 (journal draft)}
```

**Notes:**
- `type=note` → `type=report`. sphenix.cls only accepts `note|proceedings|report` — there is **no** `paper` option (verified by grep). The actual journal submission needs revtex4-2.
- Removed `\usepackage{lipsum}` — unused, and broke compile on this host. `\lipsum{}` is not called anywhere in the body.
- Added `xcolor` + `\pendingTODO{...}` macro for visible red TODO markers.
- Added `amssymb` (for `\lesssim` used in proposed §7) and `booktabs` (for the systematics + cross-section tables).
- Removed the `{\bf{sPHENIX Paper}}` placeholder text from the title.

### 3.2 Abstract (lines 35–37)

**WAS:**
> This sPHENIX Conference Note reports the measurement of the isolated prompt photon cross-section as a function of transverse energy (\etg) in proton--proton collisions at $\sqrt{s}$ = 200~\GeV, using data collected in 2024 with the sPHENIX detector at the Relativistic Heavy Ion Collider with integrated luminosity of $16.6~\mathrm{pb}^{-1}$. Photons are measured within $|\etag|<$ 0.7 and 10 $<\etg<$ 26~\GeV{}. They are reconstructed using the Electromagnetic Calorimeter and identified using electromagnetic shower shapes. An isolation selection using both the Electromagnetic and Hadronic Calorimeters is applied to suppress both fragmentation photons and background photons mostly originating from neutral-meson decays. The \etg yield is corrected for purity and efficiency, and then unfolded for detector response. The \etg-differential cross-sections are compared with theoretical predictions from Monte Carlo generator \pythia as well as next-to-leading-order perturbative Quantum Chromodynamics calculations including \jetphox.

**PROPOSED:**
> The differential cross section of isolated prompt-photon production is measured as a function of photon transverse energy (\etg) in proton--proton collisions at \comHEP\ using data recorded in 2024 with the sPHENIX detector at the Relativistic Heavy Ion Collider, corresponding to an integrated luminosity of $64.37\,\pb$. Photons are measured in $|\etag|<0.7$ and $12<\etg<32~\GeV$. They are reconstructed in the Electromagnetic Calorimeter and identified using a multivariate analysis of the electromagnetic shower shape. The measurement reports the fiducial cross section defined by the truth-level isolation requirement that the total transverse energy of final-state particles in a cone of radius $\DR=0.3$ around the photon, excluding the photon itself, is below 4~\GeV. A tighter \etg-dependent reconstruction-level isolation requirement on the topo-cluster transverse energy in $\DR=0.4$ is applied to the candidates and corrected to the truth fiducial in the analysis. A data-driven sideband technique is used to subtract the residual background from neutral-meson decays, and the resulting yield is corrected for the detector response by Bayesian iterative unfolding. The measured cross section is compared to next-to-leading-order (NLO) perturbative QCD calculations from \jetphox{} and from Vogelsang, and to the prediction of the \pythia{} 8.307 Monte Carlo generator. A comparison with the previous PHENIX measurement of inclusive direct photons in the same collision system is also presented. The result establishes the pp baseline for forthcoming sPHENIX measurements of isolated photons in nuclear collisions.

**Diff highlights:**
- "This sPHENIX Conference Note reports" → "The differential cross section ... is measured" (journal voice).
- 16.6 → 64.37 pb⁻¹ (matches `analysis.tex:303`, `selection.tex:170`, `conclusion.tex:1`, `introduction.tex:37`).
- 10–26 → 12–32 GeV.
- "electromagnetic shower shapes" → "multivariate analysis of the electromagnetic shower shape" (journal version names the BDT method up front).
- New explicit isolation definition in the abstract (matches the analysis-note's "excluding the photon energy itself" wording, no "neutrinos" gloss).
- New PHENIX comparison + heavy-ion baseline outlook sentences.

---

## 4. Section-by-section proposed changes

For the long sections below the format is: short description of intent, then the full proposed LaTeX block. Where prose is preserved verbatim from the conf-note, that fact is noted.

### 4.1 §1 Introduction (lines 41–62)

Intent: replace the obsolete framing ("In this note", 16.6 pb⁻¹, 10–26 GeV), add x-reach context, add explicit Tevatron citations, add an outline paragraph, add the heavy-ion baseline motivation specific to RHIC. The dir/frag motivation figure is held back from the journal version because the underlying PYTHIA prediction is LO-only and is not a quantitatively meaningful direct/fragmentation decomposition (per user).

**PROPOSED:**
```latex
\section{Introduction}
\label{sec:introduction}

Prompt photons in hadronic collisions are produced either directly in the hard parton--parton scattering (\textit{direct} photons) or in the collinear fragmentation of a final-state parton (\textit{fragmentation} photons). At leading order in perturbative quantum chromodynamics (pQCD) the direct contribution is dominated by quark--gluon Compton scattering, $qg \to q\gamma$, and quark--antiquark annihilation, $q\bar{q} \to g\gamma$~\cite{Owens:1986mp}. At next-to-leading order (NLO), additional channels with collinear fragmentation of final-state partons into photons contribute.

Because they do not hadronize, prompt photons provide a clean experimental probe of the hard scattering and a stringent test of pQCD over a wide range of momentum transfers~\cite{Aurenche:2006vj,Gordon:1993qc}. The cross section is particularly sensitive to the gluon parton distribution function (PDF). The kinematic reach of this measurement, $12<\etg<32~\GeV$ at $\sqrt{s}=200~\GeV$, corresponds to $x_T \equiv 2\etg/\sqrt{s}$ in the range $0.12$--$0.32$, which at leading-order $2\to 2$ kinematics samples Bjorken-$x$ values of comparable magnitude in the proton, complementing the lower-$x$ coverage of LHC photon measurements~\cite{CMS:2018qao,ATLAS:2016fta,ATLAS:2013sdn,ATLAS:2011ezy,ALICE:2019rtd,ALICE:2024kgy}. The result is also a necessary baseline for forthcoming sPHENIX measurements of photons and photon--jet correlations in nuclear collisions, where photons emerge unmodified from the quark--gluon plasma and serve as a calibrated probe of partonic energy loss~\cite{CMS:2020oen,ATLAS:2015rlt,ALICE:2024yvg}.

The fragmentation contribution and background photons originating from hadron decays are suppressed by requiring the photon to be isolated, that is, by restricting the total transverse energy of final-state particles in a cone of fixed radius $\DR \equiv \sqrt{\deta^{2}+\dphi^{2}}$ around the photon direction.

Isolated prompt-photon cross sections in \pp{} collisions have been measured extensively at the LHC~\cite{CMS:2018qao,ATLAS:2016fta,ATLAS:2013sdn,ATLAS:2011ezy,ALICE:2019rtd,ALICE:2024kgy} and in \ppb{}~\cite{ALICE:2025bnc,ATLAS:2019ery} and \pbpb{} collisions~\cite{CMS:2020oen,ALICE:2024yvg,ATLAS:2015rlt}, and at the Tevatron~\pendingTODO{add D0/CDF isolated-photon refs to the bib, e.g.\ D0 PRL 87 251805 (2001) and CDF PRD 80 111106 (2009); confirm the InspireHEP keys before insertion}. At the Relativistic Heavy-Ion Collider (RHIC)~\cite{HARRISON2003235} the PHENIX experiment measured the inclusive direct-photon cross section in \pp{} collisions at \comHEP{}~\cite{PHENIX:2012jgx} without an isolation requirement, and STAR has reported related photon-tagged measurements~\pendingTODO{cite STAR direct-photon papers}.

This Article reports the measurement of the differential cross section of isolated prompt-photon production as a function of \etg{} in \pp{} collisions at \comHEP{}, in the kinematic range $|\etag|<0.7$ and $12<\etg<32~\GeV$. The measurement uses the full \RunOne\ \pp{} data sample collected with the sPHENIX detector, corresponding to an integrated luminosity of $64.37\,\pb$. Photons are required to be isolated at the truth level, defined by the total transverse energy of final-state particles within $\DR=0.3$ around the photon, excluding the photon itself, being below 4~\GeV. The cross section is compared to NLO pQCD calculations from \jetphox~\cite{Aurenche:2006vj} and from Vogelsang~\cite{Gordon:1993qc}, to the \pythia{}~8.307~\cite{Sjostrand:2014zea} Monte Carlo generator with the Detroit tune~\cite{Aguilar:2021sfa}, and to the previous inclusive direct-photon measurement by PHENIX~\cite{PHENIX:2012jgx}.

The remainder of this Article is organised as follows. Section~\ref{sec:sphenix_info} describes the sPHENIX detector. Section~\ref{sec:analysis} details the analysis procedure, including event selection and Monte Carlo simulation, photon reconstruction and identification, signal extraction with the data-driven sideband subtraction, unfolding, the efficiency corrections, and the cross-section determination. Systematic uncertainties are discussed in Section~\ref{sec:systematics}, and the results are presented in Section~\ref{sec:results}. A summary follows in Section~\ref{sec:summary}.
```

### 4.2 §2 sPHENIX detector (lines 65–75)

Intent: leave this section unchanged from the conf-note. The sPHENIX detector description is unified across collaboration papers and any wording change here would diverge from the convention. The Pass-6 expansions (coordinate system, EMCal granularity inline, in-situ resolution paragraph, MBD detail, solenoid promotion) are reverted. The journal version reproduces the conf-note prose verbatim.

**PROPOSED (verbatim copy of `PPG12-Paper/main.tex` §2, lines 64–76):**
```latex
\section{sPHENIX detector}
\label{sec:sphenix_info}

sPHENIX~\cite{PHENIX:2015siv, Belmont:2023fau} is a new detector designed to measure jet and heavy-flavor probes of the quark-gluon plasma (QGP) created in Au+Au collisions at the RHIC. A precision tracking system enables measurements of heavy-flavor and jet-substructure observables while the electromagnetic and hadronic calorimeter system is crucial for measuring the energy of jets and identifying direct photons and electrons.

Going outwards starting from the beam line, sPHENIX comprises the following subsystems~\cite{Sphenix:TDR}: the MAPS-based Vertex Detector (MVTX); the INTermediate Tracker (INTT); the Time Projection Chamber (TPC)~\cite{Klest:2020sdb}; the Time Projection Chamber Outer Tracker (TPOT)~\cite{Aune:2024ysr}; the Electromagnetic Calorimeter (EMCAL)~\cite{Aidala:2020toz,sPHENIX:2017lqb}; the Inner Hadronic Calorimeter (IHCAL)~\cite{sPHENIX:2017lqb}; the 1.4\,T superconducting solenoid magnet~\cite{OConnor:1998llb} and the Outer Hadronic Calorimeter (OHCAL)~\cite{sPHENIX:2017lqb}.
Except for TPOT, all detectors have full azimuthal coverage and span $\left|\eta\right| < 1.1$ in pseudorapidity. sPHENIX also includes a number of forward detectors, namely the Minimum Bias Detectors (MBD), the sPHENIX Event Plane Detectors (sEPD), and the Zero Degree Calorimeters (ZDC), that includes the Shower Maximum Detector (SMD).

sPHENIX began its commissioning process in RHIC Run-2023 with Au+Au collisions at \sqsn = 200~\GeV. During RHIC Run-2024, sPHENIX collected a large sample of transversely polarized \pp{} physics data at $\sqrt{s}$=200~\GeV{}, alongside a smaller sample of Au+Au data to complete its commissioning phase in that collision system.
```

### 4.3 §3 Analysis Procedure — single-section parent (replaces conf-note §3)

Intent: in Option A the journal §3 is a single section with seven subsections, mirroring the conf-note layout. The proposed LaTeX block opens with a single `\section{Analysis Procedure}` and `\label{sec:analysis}`. Subsections §4.3.1–§4.3.7 below are the seven `\subsection{...}` blocks of journal §3.

**PROPOSED (parent header):**
```latex
\section{Analysis Procedure}
\label{sec:analysis}
```

### 4.3.1 §3.1 Event Selection (replaces conf-note §3.1)

Intent: update lumi (16.6 → 64.37 \pb{}) and the vertex cut (30 → 60 cm). Keep the granularity at conf-note level (per user): no trigger bit numbers, no internal-note cites, no per-paragraph splits, no dataset-summary table. The conf-note presents triggers + vertex selection in a single short paragraph and the journal version preserves that structure.

**PROPOSED (replaces conf-note §3.1, lines 82–83):**
```latex
\subsection{Event Selection}
\label{ssec:eventsel}
The data analysed in this Article were recorded by sPHENIX during \RunOne{} at \comHEP{}, corresponding to an integrated luminosity of $\mathscr{L}=64.37^{+5.88}_{-4.34}\,\pb$. Events were selected by a photon trigger which requires the total energy in any 8$\times$8 EMCAL tower window (each tower covers a $0.024\times 0.024$ segment in $\eta$--$\phi$) to exceed 4~\GeV{}, in coincidence with a minimum-bias signal in the MBD requiring at least one photomultiplier firing on each side. The bin-by-bin efficiency of the photon trigger is measured \emph{in situ} from a minimum-bias reference sample as a function of the leading-cluster transverse energy and is corrected for at the level of the photon yield. Events are required to have a reconstructed MBD $z$-vertex in $|\vz|<60~\text{cm}$. The integrated luminosity is determined from the live count of the MBD coincidence and the in-situ measurement of the MBD inelastic cross section, $\sigma_{\mathrm{MBD}}=25.2^{+2.3}_{-1.7}$~mb, and the residual fraction of events with a vertex outside the analysis fiducial is restored at the cross-section level via the MBD-vertex efficiency described in Sec.~\ref{ssec:efficiency}.
```

### 4.3.2 §3.2 Monte Carlo Simulations (replaces conf-note §3.2)

Intent: standardise PYTHIA 8.307, add MC sample table, fold the double-interaction blending content into an inline `\paragraph{...}`, describe truth-vertex reweighting in a single sentence (per user — not a multi-paragraph derivation), state the truth signal definition explicitly.

**PROPOSED (replaces conf-note §3.2, lines 85–89):**
```latex
\subsection{Monte Carlo Simulations}
\label{ssec:mc}
Simulated events are generated with \pythia~8.307~\cite{Sjostrand:2014zea} using the Detroit tune~\cite{Aguilar:2021sfa}, with both \texttt{HardQCD:all = on} and \texttt{PromptPhoton:all = on} processes enabled in every production. The generated samples are partitioned into prompt-photon (``signal'') and inclusive-jet (``background'') productions by a truth-level event filter on the leading truth photon (signal) or leading truth jet (background) within $|\eta|<1.5$ above a per-sample $\hat{p}_{\mathrm{T}}$ threshold. The generated events are propagated through the full sPHENIX simulation, based on \geant{}~4~\cite{AGOSTINELLI2003250}, and are reconstructed with the same algorithms as the data. A noise contribution is included in the simulation that reproduces the per-tower pedestal RMS observed in minimum-bias data \pendingTODO{confirm the exact procedure (random-tower overlay vs.\ parametric Gaussian) with the EMCal calibration team}.
The set of generated samples and their truth-level leading-pT windows, used to tile the analysed \etg{} range without overlap, are listed in Table~\pendingTODO{add Table~\ref{tab:mc_samples} mirroring the analysis-note pair (\texttt{table:MCdataset} for the eight SI/DI generation rows including \texttt{jet40}/\texttt{jet40\_double}, and \texttt{table:stitching} for the truth leading-pT windows)}.

At the truth level, signal photons are defined as prompt photons (direct and fragmentation) that satisfy the truth-level isolation requirement: the total transverse energy $\isoETtruth$ of final-state particles in a cone $\DR=0.3$ around the photon direction, excluding the photon energy itself, is below 4~\GeV. Reconstructed photon candidates are matched to truth photons using the track-ID-based truth-cluster association provided by the \texttt{CaloRawClusterEval} module of the GEANT-based simulation, which assigns to each reconstructed cluster the truth primary particle whose simulated hits contribute the largest summed energy across the cluster's towers.

A per-event truth-vertex reweighting is applied to the simulated events so that the reconstructed-vertex distribution in the simulation matches the data after the standard analysis cuts.

\paragraph{Double-interaction handling.}
At the instantaneous luminosities delivered during \RunOne{}, a non-negligible fraction of bunch crossings contain more than one inelastic \pp{} interaction. Multiple interactions distort the reconstructed primary vertex, the cluster kinematics, and the shower-shape distributions of the photon candidates. To account for this effect, dedicated full-\geant{} double-interaction (DI) Monte Carlo samples are generated, in which each event combines a hard \pythia{} event with a minimum-bias \pythia{} event placed at an independent vertex along the beam line. The cluster-weighted DI fractions are $f_{\mathrm{DI}}^{0\,\mathrm{mrad}}=0.224$ and $f_{\mathrm{DI}}^{1.5\,\mathrm{mrad}}=0.079$ for the two crossing-angle configurations, derived from a Poisson decomposition of the measured MBD-coincidence rate into mean interactions per bunch crossing $\mu$, with triple and higher orders folded into the double-interaction yield, and validated against a data-driven $\chi^2$ fit of the cluster photon-ID BDT-score distribution at the preselection cut.
The simulated photon yield used in the analysis is a linear combination of the single-interaction (SI) and DI samples weighted by these fractions and by the corresponding integrated-luminosity weights of each crossing-angle configuration.
```

### 4.3.3 §3.3 Photon Reconstruction and Identification (replaces conf-note §3.3)

Intent: replace the "five shower-shape variables" hand-wave with the explicit `weta` formula and the parametric BDT cut. Fold cluster reconstruction, pre-selection (NCB BDT), and the photon-ID BDT into inline `\paragraph{...}` blocks within a single `\subsection{...}`. Granularity ledger (Wave-1C): 1.6× the conf-note baseline, within tolerance — no further trimming required. Per the user's mid-Wave-1 addition, two shower-shape distribution figures (`weta_cogx` and `e32_to_e35`, eta0/pt2/cut1, NCB-pre-selection level) are inserted with surrounding prose to motivate the BDT input variables.

**PROPOSED (replaces conf-note §3.3, lines 92–99):**
```latex
\subsection{Photon Reconstruction and Identification}
\label{ssec:reco}
\paragraph{Cluster reconstruction.}
Photon candidates are reconstructed from clusters of EMCal towers. Adjacent towers (sharing an edge) above a 70~\MeV{} threshold and passing tower quality cuts are aggregated into clusters by the standard sPHENIX \texttt{RawClusterBuilderTemplate} algorithm. The clusters grow until no adjacent towers are above the threshold, after which a cluster-splitting algorithm separates merged showers, and the resulting sub-clusters are the photon candidates used in the analysis.
The cluster energy and direction are corrected for the position of the EMCal shower in the calorimeter and for the longitudinal collision vertex.
The energy resolution of the reconstructed photon candidates, measured \emph{in situ} from the response matrix of matched truth photons, is approximately ${\sim}5\%$ across the analysis range. \pendingTODO{quote final $\sigma_E/E$ parametrisation, replacing the approximate number with the fit.}

\paragraph{Pre-selection.}
A pre-selection is applied to all photon candidates before the photon-ID BDT in order to remove non-collision backgrounds (NCB, predominantly beam-induced clusters with anomalously broad $\eta$ profiles and out-of-time EMCal--MBD timing) and unphysical reconstructed clusters. The pre-selection consists of rectangular cuts on the cluster shower-shape variables $E_{\mathrm{t}1}$, $E_{1\times 1}/E_{3\times 3}$, and $E_{3\times 2}/E_{3\times 5}$, together with a multivariate ``NCB BDT'' trained in \textsc{XGBoost} on 25 cluster-level features and exported to \textsc{TMVA} for application at the reconstruction level. The training uses a data-derived NCB sample (clusters with anomalously early EMCal--MBD timing) against an inclusive-MC physics-cluster sample. The NCB BDT score is required to exceed $0.5$, retaining ${\approx}99\%$ of physics-cluster signal at ${\sim}99\%$ purity, derived at the $22<\etg<28~\GeV$ representative bin from a data-driven template fit to the EMCal--MBD timing distribution \pendingTODO{quote the lowest-\etg{} bin (\texttt{pt0}, $10<\etg<14~\GeV$) signal-retention and purity values from the analysis-note Fig.~\ref{fig:npb_purity_yield} once that figure is regenerated}. The full pre-selection cut values are summarised in Table~\pendingTODO{add Table~\ref{tab:photonid}, mirror analysis-note Table~3}.

\paragraph{Photon-ID BDT.}
After pre-selection, photon candidates are classified using a Boosted Decision Tree (BDT) trained on the EMCal shower-shape variables. Photons and electromagnetic clusters produce narrower transverse showers in the EMCal than the merged cluster of a high-\pt{} neutral meson decay ($\pizero\to\gamma\gamma$). The dominant input variable is the energy-weighted second moment of the tower $\eta$ distribution computed in the cluster centre-of-gravity reference frame,
\begin{equation}
\wetacogx = \sqrt{\frac{\sum_i E_i\,(\eta_i-\bar\eta)^2}{\sum_i E_i}}, \qquad \bar\eta = \frac{\sum_i E_i\,\eta_i}{\sum_i E_i},
\label{eq:weta}
\end{equation}
where the sum is over towers $i$ in the cluster. An analogous variable $w_{\phi}^{\mathrm{COGX}}$ is computed in the $\phi$ direction. Additional variables include the energy fraction $\etone$ defined as the fraction of the cluster energy contained in the $2\times 2$ block of towers around the cluster centre-of-gravity (CoG) — including the CoG tower itself — and the energy ratios $E_{3\times 2}/E_{3\times 5}$ and $E_{1\times 1}/E_{3\times 3}$ which separate point-like electromagnetic showers from the broader hadronic profile~\pendingTODO{full list of the 11 photon-ID BDT input variables, aligned with analysis-note Table~\ref{tab:bdt\_features}}.

The discriminating power of the shower-shape variables is illustrated in Fig.~\ref{fig:showershape_dis}, which compares the distributions in data, prompt-photon signal simulation, and inclusive-jet background simulation at the NCB pre-selection cut for a representative $\etg$ bin. Signal photons populate the narrow-shower side of \wetacogx{} and the hadronic-profile-suppressed side of $E_{3\times 2}/E_{3\times 5}$, while the data are dominated by the inclusive-jet background at the pre-selection level.

\begin{figure}[tbp!]
    \centering
    \includegraphics[width=0.45\linewidth]{figures/dis_mixed_weta_cogx_eta0_pt2_cut1.pdf}
    \includegraphics[width=0.45\linewidth]{figures/dis_mixed_e32_to_e35_eta0_pt2_cut1.pdf}
    \caption{Distributions of the shower-shape variables \wetacogx{} (left) and $E_{3\times 2}/E_{3\times 5}$ (right) at the NCB pre-selection level in a representative \etg{} bin ($14<\etg<16~\GeV$, $|\etag|<0.7$). Black points are data, the red histogram is the prompt-photon signal simulation, and the blue histogram is the inclusive-jet background simulation, both normalised to the data integral. The signal/background separation in the two variables motivates their use as inputs to the photon-ID BDT. \pendingTODO{regenerate against the all-range allz nominal; choose representative \etg{} bin consistent with the rest of the section.}}
    \label{fig:showershape_dis}
\end{figure}

The shower-shape variables are combined into a single photon-ID score using the \textsc{XGBoost} implementation~\pendingTODO{cite the Chen-Guestrin XGBoost paper}, trained on \pythia{} prompt-photon signal against \pythia{} inclusive-jet background with 11 input features.
The ``tight'' photon-identification (\gid) selection is defined parametrically in $\etg$ as
\begin{equation}
\mathrm{BDT} \;>\; 0.8333 - 0.00333\,\etg/\GeV,
\label{eq:tight_bdt}
\end{equation}
chosen to maintain a roughly flat 80\% identification efficiency relative to the pre-selection across the analysis range. The non-tight \gid{} band, used as the photon-ID sideband in the data-driven purity calculation, is defined by the parametric range $0.7333 - 0.01333\,\etg/\GeV < \mathrm{BDT} < 0.6667 + 0.00333\,\etg/\GeV$.
The systematic uncertainty associated with the NCB cut is evaluated by varying the score threshold to $0.3$ and $0.7$, and an independent closure cross-check on a back-to-back-recoil-jet-enriched data sample with and without the NCB cut shows the unfolded cross section to agree to within $\sim 1\%$ across the analysis \etg{} range.
```

### 4.3.4 §3.4 Signal Extraction (replaces conf-note §3.4)

Intent: combine the parametric isolation requirement and the ABCD background-subtraction procedure into a single subsection mirroring the conf-note layout, with the isolation content as a `\paragraph{...}` lead-in. Re-anchor the parametric isolation `IsoET < 0.490 + 0.0370·ET^reco` (topo R = 0.4). Quote purity as rising from $\sim 0.5$ at the lowest \etg{} bin to $\sim 0.9$ above $26~\GeV$. The isolation lead-in cites the ATLAS topo-cluster algorithm (`ATLAS:2017topoclusters`) without restating the per-cell noise scales or the seed/grow/close significance multipliers, and the MC iso-fudge values are deferred to §4 systematics. The purity is fitted with a smooth functional form, with the Padé/erf comparison kept as a §4 systematic-source detail.

**PROPOSED (replaces conf-note §3.4, lines 101–170):**
```latex
\subsection{Signal Extraction}
\label{ssec:abcd}
\paragraph{Isolation.}
A photon-isolation requirement is imposed to suppress fragmentation photons and the residual background of high-\pT{} meson decays that survive the shower-shape selection. The reconstructed isolation transverse energy $\isoET^{\mathrm{reco}}$ is computed as the scalar sum of the transverse energies of topologically reconstructed clusters (``topo-clusters'') in the calorimeter system, following the algorithm originally developed by ATLAS~\cite{ATLAS:2017topoclusters} and adapted to the sPHENIX three-layer EMCal/IHCal/OHCal stack. Topo-clusters whose axes fall inside a cone of radius $\DR=0.4$ around the photon direction are summed, and the photon candidate's own reconstructed-cluster $E_T$ is subtracted from the cone sum to approximately remove the candidate's energy deposit from the topo-cluster total.
The reconstructed-photon isolation requirement is parameterised as
\begin{equation}
\isoET^{\mathrm{reco}} \;<\; 0.490\,\GeV \;+\; 0.0370 \,\etreco,
\label{eq:iso_cut}
\end{equation}
with coefficients tuned to maintain an 80\% isolation efficiency on simulated truth-isolated direct photons across the analysis range. Non-isolated candidates used as background in the ABCD method are required to satisfy
\begin{equation}
\isoET^{\mathrm{reco}} \;>\; 0.8\,\GeV \;+\; (0.490\,\GeV + 0.0370\,\etreco).
\label{eq:noniso_cut}
\end{equation}
A small empirical correction is applied to the simulated $\isoET$ to absorb residual data-MC mismatch, with the variation propagated as a systematic uncertainty (Sec.~\ref{sec:systematics}).

\paragraph{ABCD background subtraction.}
The residual background to tight, isolated photon candidates is subtracted using a data-driven sideband (``ABCD'') technique adapted from inclusive-photon measurements at the LHC~\cite{ATLAS:2016fta,ATLAS:2017nah}. Four populations of photon candidates are defined according to the photon-ID BDT score (tight or non-tight) and the isolation requirement (isolated or non-isolated): A (tight, isolated -- signal region), B (tight, non-isolated), C (non-tight, isolated), and D (non-tight, non-isolated), as illustrated in Fig.~\ref{fig:sideband_diagram}.

\begin{figure}[tbp!]
    \centering
    \includegraphics[width=0.45\linewidth]{figures/sideband_diagram.png}
    \caption{Schematic of the signal region (A) and sideband regions (B, C, D) used in the data-driven background subtraction.}
    \label{fig:sideband_diagram}
\end{figure}

Assuming the photon-ID and isolation observables are uncorrelated for background, the number of signal photons in region A in the absence of signal leakage is
\begin{equation}
N^{A}_{\mathrm{sig}}
= N^{A}_{\mathrm{raw}}
- N^{B}_{\mathrm{raw}}
\,\cdot\,
\frac{N^{C}_{\mathrm{raw}}}{N^{D}_{\mathrm{raw}}}.
\label{eq:purity_noleakage}
\end{equation}
Because of imperfect ID--iso decorrelation and finite signal leakage into the sideband regions, Eq.~\ref{eq:purity_noleakage} is generalised to
\begin{equation}
N^{A}_{\mathrm{sig}}
= N^{A}_{\mathrm{raw}}
-
  \bigl(
    N^{B}_{\mathrm{raw}}
    - f^{B,\mathrm{MC}}\,N^{A}_{\mathrm{sig}}
  \bigr)
  \,\cdot\,
  \frac{
    N^{C}_{\mathrm{raw}}
    - f^{C,\mathrm{MC}}\,N^{A}_{\mathrm{sig}}
  }{
    N^{D}_{\mathrm{raw}}
    - f^{D,\mathrm{MC}}\,N^{A}_{\mathrm{sig}}
  },
\label{eq:purity_leak}
\end{equation}
where the leakage fractions $f^{X,\mathrm{MC}}\equiv N^{X,\mathrm{sig}}_{\mathrm{MC}}/N^{A,\mathrm{sig}}_{\mathrm{MC}}$ are determined from truth-matched signal photons in the inclusive-prompt-photon simulation. The purity, defined as the signal fraction in the signal region,
\begin{equation}
\mathcal{P}(\etg) \;\equiv\; \frac{N^{A}_{\mathrm{sig}}(\etg)}{N^{A}_{\mathrm{raw}}(\etg)},
\label{eq:purity_def}
\end{equation}
is determined in the same $\etg$ bins as the cross-section measurement, with bootstrap statistical uncertainties. The leakage-corrected purity rises from $\mathcal{P}\sim 0.5$ at the lowest \etg{} bin to $\mathcal{P}\sim 0.9$ above $26~\GeV$. \pendingTODO{quote bin-by-bin purity from live \texttt{Photon\_final\_*.root} once the all-range allz nominal is finalised.}

The BDT and isolation requirements are optimised to minimise the residual ID--iso correlation, and the residual MC closure is propagated as a systematic uncertainty (Sec.~\ref{sec:systematics}).

Figure~\ref{fig:isoDist} illustrates the data-MC consistency by overlaying the $\isoET$ distribution in tight (signal) and non-tight (background) data, together with the simulation of tight signal photons. The non-tight data and tight signal MC distributions are stacked to demonstrate the mixture observed in the data.

\begin{figure}[tbp!]
    \centering
    \includegraphics[width=0.45\linewidth]{figures/h1D_iso_nom_6_8.pdf}
    \caption{$\isoET$ distributions of tight \gid{} candidates and non-tight \gid{} candidates in data, and of tight \gid{} signal photons in simulation, in a representative \etg{} bin. The non-tight-data and signal-MC histograms are stacked. \pendingTODO{regenerate this figure for the full $12<\etg<32~\GeV$ range, multi-panel by \etg{} bin if practical.}}
    \label{fig:isoDist}
\end{figure}

Figure~\ref{fig:purity} shows the resulting purity as a function of $\etg$. The signal-leakage correction is at the few-percent level, reflecting the high purity of the tight signal region. To reduce bin-to-bin statistical fluctuations, the leakage-corrected purity values are fitted with a smooth functional form, and the fitted value at each bin centre is used to correct the residual background. The functional-form spread is propagated as a systematic uncertainty (Sec.~\ref{sec:systematics}).

\begin{figure}[tbp!]
    \centering
    \includegraphics[width=0.45\linewidth]{figures/purity_nom.pdf}
    \caption{Purity as a function of \etg{} with and without the signal-leakage correction. The leakage-corrected purity is fitted with a smooth functional form, with the shaded band showing the 68.3\% confidence interval of the fit. \pendingTODO{regenerate against the all-range allz nominal.}}
    \label{fig:purity}
\end{figure}
```

### 4.3.5 §3.5 Unfolding (replaces conf-note §3.5)

**PROPOSED (replaces conf-note §3.5, lines 172–174):**
```latex
\subsection{Unfolding}
\label{ssec:unfolding}
The purity-corrected photon-yield distribution is unfolded for the residual detector-response effects (notably the EMCal energy scale and resolution and the bin-to-bin migrations) using the iterative D'Agostini Bayesian method~\cite{DAgostini:2010hil} as implemented in \textsc{RooUnfold}~3.0.5~\cite{Adye:2011gm}. The response matrix is constructed from the prompt-photon simulation, with reconstructed-\etg{} bin edges $\{10,12,14,16,18,20,22,24,26,28,32,36\}~\GeV$ and truth-\etg{} bin edges $\{8,10,12,14,16,18,20,22,24,26,28,32,36,45\}~\GeV$. The bins outside the reported $12<\etg<32~\GeV$ analysis range are retained as overflow buffers. The response-matrix truth prior is reweighted to bring it into agreement with the data \etg{} spectrum. The number of iterations is set to two, chosen so that successive-iteration changes fall below the per-bin statistical uncertainty while the iteration-dependent bias remains subdominant to the other systematics. The corresponding \etg-binned closure tests are summarised in Sec.~\ref{sec:systematics}.
```

### 4.3.6 §3.6 Efficiency Correction (replaces conf-note §3.6)

**PROPOSED (replaces conf-note §3.6, lines 176–183):**
```latex
\subsection{Efficiency Correction}
\label{ssec:efficiency}
The unfolded yield is corrected for the photon reconstruction (\effreco), identification (\effid), isolation (\effiso), and MBD-vertex ($\varepsilon_{\mathrm{vtx}}$) efficiencies, all of which are evaluated bin-by-bin in $\etgtruth$ from the prompt-photon simulation. The photon-trigger efficiency is included as a separate event-by-event weight at fill time as discussed in Sec.~\ref{ssec:eventsel}. The MBD-vertex efficiency accounts for the loss of events whose reconstructed vertex falls outside the $|\vz|<60$~cm fiducial.
The breakdown of the efficiencies as a function of $\etgtruth$ is shown in Fig.~\ref{fig:eff_photon}.

\begin{figure}[hbtp!]
    \centering
    \includegraphics[width=0.45\linewidth]{figures/eff_photon.pdf}
    \caption{Reconstruction (\effreco), identification (\effid), isolation (\effiso) and convolved (\efftot) efficiencies as a function of truth photon \etgtruth. \pendingTODO{regenerate against the all-range allz nominal; add MBD-vertex efficiency curve $\varepsilon_{\mathrm{vtx}}$.}}
    \label{fig:eff_photon}
\end{figure}
```

### 4.3.7 §3.7 Cross-Section Determination (replaces conf-note §3.7)

**PROPOSED (replaces conf-note §3.7, lines 185–202):**
```latex
\subsection{Cross-Section Determination}
\label{ssec:xsec}
The differential cross section per pseudorapidity unit and per photon transverse-energy unit is
\begin{equation}
\frac{\mathrm{d}^{2}\sigma}{\mathrm{d}\etg\,\mathrm{d}\etag}
= \frac{1}{\mathscr{L}}\;
\frac{Y^{\mathrm{rec}}(\etg)}{\mathcal{E}(\etg)\,\Delta\etg\,\Delta\etag},
\label{eq:cs}
\end{equation}
where $Y^{\mathrm{rec}}(\etg)$ is the unfolded, purity-corrected photon yield,
\begin{equation}
Y^{\mathrm{rec}}(\etg)
= \mathrm{Unfold}\!\Bigl[\,\Nsig(\etg)\;\mathcal{P}(\etg)\,\Bigr],
\label{eq:rec}
\end{equation}
$\Nsig(\etg)$ is the raw yield of tight, isolated photon candidates, $\mathcal{P}(\etg)$ is the purity, $\mathcal{E}(\etg) = \effreco\cdot\effid\cdot\effiso\cdot\varepsilon_{\mathrm{vtx}}$ is the combined efficiency, and $\Delta\etg$ and $\Delta\etag$ are the bin widths. The photon-trigger efficiency $\varepsilon_{\mathrm{trig}}(\etg)$ is applied per-cluster as the inverse weight $1/\varepsilon_{\mathrm{trig}}$ at fill time prior to the purity, unfolding, and combined-efficiency corrections, as described in Sec.~\ref{ssec:eventsel}, and is therefore not multiplied into $\mathcal{E}$.
The integrated luminosity, $\mathscr{L}=64.37^{+5.88}_{-4.34}\,\pb$, is determined as described in Sec.~\ref{ssec:eventsel}.
```

### 4.4 §4 Systematic Uncertainties (replaces conf-note §4)

Intent: replace the four-paragraph conf-note text with nine `\paragraph{...}` blocks matching the analysis-note's 14-source breakdown. Update lumi systematic. Add new sources: DI fraction, NCB cut, MC closure of ABCD.

**PROPOSED (replaces conf-note §4, lines 205–221):**
```latex
\section{Systematic Uncertainties}
\label{sec:systematics}
Each systematic source is evaluated by repeating the full analysis chain with the corresponding variation. For each \etg{} bin, the relative shifts of the cross section are taken as the systematic uncertainties, and added in quadrature to obtain the total. A bin-by-bin breakdown is given in Table~\pendingTODO{add Table~\ref{tab:syst_summary}, mirroring analysis-note Table 5} and visualised in Fig.~\ref{fig:syst_sum}.

\paragraph{Photon energy scale.}
The EMCal absolute energy scale is calibrated using \pizero{} and $\eta$ mass peaks reconstructed in the data and in simulation. The residual data-MC difference is propagated as a $\pm 2.6\%$ shift of the reconstructed cluster energy. This is the dominant systematic uncertainty in the high-\etg{} bins, reaching $-22.6\%$ on the down side and $+29.3\%$ on the up side at the high-\etg{} end of the analysis range. \pendingTODO{the 2.6\% shift remains a placeholder pending the final EMCal calibration; replace with the calibration paper's value or shrink-justify the variation.}

\paragraph{Photon energy resolution.}
A per-cluster Gaussian smearing of $\sigma_E/E=4\%$ is applied to the MC reconstructed cluster energy to absorb the residual data-MC mismatch in the photon-cluster energy resolution. The systematic uncertainty is derived by varying the smearing between zero and twice the nominal value, yielding maximum bin-by-bin shifts of $-17.1\%$ on the down side and $+6.3\%$ on the up side relative to the nominal. \pendingTODO{4\% is a placeholder pending the final calibration. Revisit jointly with the energy-scale entry above.}

\paragraph{Purity.}
The purity systematic comprises three contributions added in quadrature: variation of the tight and non-tight \gid{} thresholds and the non-isolation boundary, variation of the smooth functional form (Pad\'e[1/1] versus erf, ${\sim}0.3\%$) together with the $\pm 1\sigma$ shift of the nominal-fit parameter posterior, and an MC-closure correction $C_\mathrm{MC}(\etg)\equiv\mathcal{P}_\mathrm{truth}^\mathrm{MC}/\mathcal{P}_\mathrm{ABCD}^\mathrm{MC}$ derived from the inclusive-jet simulation. The closure term is the dominant component at low \etg{}, and the combined purity systematic ranges from ${\sim}8\%$ at low \etg{} to ${\sim}4\%$ above $22~\GeV$.

\paragraph{Efficiency corrections.}
The combined efficiency systematic is taken from the residual additive pedestal mismodelling of the simulated $\isoET$. The MC iso-fudge $\isoET^\mathrm{MC,corr}=a\cdot\isoET^\mathrm{MC}+b$ with $(a,b)=(1.2,0.1\,\GeV)$ is varied by setting $b\to 0$ while keeping $a=1.2$, and the resulting one-sided shift is symmetrised. The variation reaches $\pm 4.5\%$ relative to the nominal at its peak. The MBD-vertex efficiency carries a separate placeholder envelope of $\pm 6\%$, flat in \etg{}, taken from the PPG10 data-driven MBD-coincidence efficiency~\pendingTODO{replace with the PPG10-derived MBD-vertex eff systematic; the current $\pm 6\%$ is a placeholder. Memory \texttt{project\_mbd\_vertex\_eff\_pT\_dependence.md} documents a photon-vs-jet cross-check at the $1\%$ level that may justify shrinking.}.

\paragraph{Unfolding.}
The unfolding systematic comprises two contributions, taken in quadrature. The first is the prior-reweighting variation, in which the data-driven smoothed truth prior is replaced by the un-reweighted MC truth prior of the response matrix. The second is the iteration-count envelope, taken as the per-bin maximum of the deviations of the $N_{\mathrm{iter}}=3$ and $N_{\mathrm{iter}}=4$ unfolds from the nominal $N_{\mathrm{iter}}=2$, symmetrised. The combined unfolding systematic stays $\lesssim 6\%$ across the analysis range.

\paragraph{Double interactions.}
The DI admixture fractions $f_{\mathrm{DI}}^{0\,\mathrm{mrad}}$ and $f_{\mathrm{DI}}^{1.5\,\mathrm{mrad}}$ are varied to the alternative best-fit values from the data-driven $\chi^2$ fit of the cluster photon-ID BDT-score distribution at the preselection cut, and the analysis is rerun for each variation. The single one-sided shifts are symmetrised, and the associated cross-section shift stays below $2.5\%$ across the analysis range.

\paragraph{Background-rejection cut (NCB).}
The NCB BDT score threshold is varied between $0.3$ and $0.7$ about the nominal value of $0.5$, and the associated relative cross-section shift is asymmetric, $-1.5\%$ to $+3.6\%$, increasing with \etg.

\paragraph{Photon-trigger correction.}
The bin-by-bin photon-trigger correction described in Sec.~\ref{ssec:eventsel} absorbs the per-bin Gumbel-CDF turn-on uncertainty into the central value, and no separate plateau systematic is assigned. \pendingTODO{confirm this treatment in the journal version, and otherwise propagate the $1\sigma$ fit-covariance band as a per-bin systematic.}

\paragraph{Tower-acceptance cross-check.}
A cross-check based on the EMCal $\phi$-symmetry of the analysis-fiducial tower acceptance is performed by comparing the unfolded cross section under nominal and $\phi$-symmetry tower-mask variants. The variation is propagated as a cross-check rather than as a quadrature term, consistent with the analysis-note treatment.

\paragraph{Truth-vertex-reweight closure cross-check.}
The truth-vertex reweighting described in Sec.~\ref{ssec:mc} is closure-tested per crossing-angle period against the reconstructed-vertex profile in the data. The residual non-closure, ${\sim}1$--$2\%$ and concentrated at low \etg{}, is propagated as a cross-check rather than as a quadrature term, consistent with the analysis-note treatment.

\paragraph{Luminosity.}
The integrated luminosity is determined from the live MBD\,N\&S coincidence count and the in-situ van~der~Meer measurement of the MBD inelastic cross section, $\sigma_{\mathrm{MBD}}=25.2^{+2.3}_{-1.7}$~mb. The asymmetric uncertainty propagates to a normalisation systematic of $+9.13\%/-6.75\%$, which is fully correlated across all \etg{} bins.

The total systematic uncertainty is strongly dependent on \etg{} and is dominated by the energy-scale term at high \etg{} and by the purity-closure term at low \etg{}, with the bin-by-bin envelope reaching $-30.4\%$ on the down side and $+35.1\%$ on the up side near the high-\etg{} end of the analysis range (Table~\ref{tab:syst_summary}). Because the energy-scale uncertainty is the dominant contribution, the total systematic uncertainty has a strong bin-to-bin correlation and is shown as a band rather than as bin-uncorrelated points in the result figures of Sec.~\ref{sec:results}.

\begin{figure}[tbp!]
    \centering
    \includegraphics[width=0.8\linewidth]{figures/syst_sum_rel.pdf}
    \caption{Breakdown of relative systematic uncertainties as a function of \etg. \pendingTODO{regenerate after the systematic redo (\texttt{tightbdt50/70}, \texttt{novtxcut}, \texttt{ntbdtmin05/10} retired; SYST\_GROUPS["eff"] no longer contains \texttt{tight\_bdt}); update legend.}}
    \label{fig:syst_sum}
\end{figure}
```

### 4.5 §5 Results (replaces conf-note §5)

Intent: PYTHIA version standardised to 8.307. \jetphox{} and Vogelsang nominals are now both at CT18NLO with the truth-level isolation requirement applied (matching `plot_final_selection.C`). The final-plot lower panel adds a four-PDF comparison (CT18NLO, NNPDF3.1, CTEQ6.6, MSHT20NLO), per the live macro third pad. PHENIX rescaling is described explicitly. The $x_T$-scaling paragraph is held back for now (per user — not yet sign-off-ready). The data-vs-theory ratio is stated as "consistent within uncertainties" without a $\chi^2/\mathrm{ndf}$ for the moment.

**PROPOSED (replaces conf-note §5, lines 226–252):**
```latex
\section{Results}
\label{sec:results}

The measured differential cross section of isolated prompt-photon production as a function of \etg{} in $|\etag|<0.7$ is shown in Fig.~\ref{fig:final}, together with predictions from the \pythia{}~8.307 Monte Carlo with the Detroit tune~\cite{Aguilar:2021sfa}, the NLO pQCD calculation from \jetphox{}~v$1.3.1\_4$~\cite{Aurenche:2006vj}, and the NLO pQCD calculation by Vogelsang~\cite{Gordon:1993qc}.
Both NLO calculations use the \textsc{CT18NLO} parton distribution functions~\cite{Hou:2019efy} as their nominal PDF set, and the \jetphox{} prediction additionally includes the BFG~set~II~\cite{Bourhis:1997yu} parton-to-photon fragmentation functions for the fragmentation contribution. The renormalisation, factorisation, and (for \jetphox{}) fragmentation scales are set to $\mu_{R}=\mu_{F}=\mu_{f}=\etg$, and a scale uncertainty is obtained by simultaneously varying them by factors of $1/2$ and $2$.
The \pythia{}, \jetphox{}, and Vogelsang predictions all impose the same truth-level isolation requirement as the data, $\isoETtruth<4~\GeV{}$ in $\DR=0.3$, so the comparison is carried out at the same fiducial level. The measurement is consistent with the NLO pQCD predictions over the analysis range $12<\etg<32~\GeV$.

The measured cross section spans approximately two orders of magnitude over the reported range $12<\etg<32~\GeV$. The statistical uncertainty is sub-percent at low \etg{} and reaches $\sim 30\%$ in the highest reported bin. The total systematic uncertainty reaches $-30.4\%$ on the down side and $+35.1\%$ on the up side at the high-\etg{} end of the reported range, dominated by the EMCal energy-scale uncertainty at high \etg{} and by the ABCD MC-closure uncertainty at low \etg{}.

\begin{figure}[tbp!]
    \centering
    \includegraphics[width=0.85\linewidth]{figures/final.pdf}
    \caption{Differential cross section of isolated prompt-photon production as a function of \etg{} in \pp{} collisions at \comHEP{}, in the kinematic range $|\etag|<0.7$ and $12<\etg<32~\GeV$. The data are compared in the upper panel with predictions from \pythia{}~8.307 (Detroit tune), \jetphox{} (CT18NLO), and the NLO pQCD calculation by Vogelsang (CT18NLO), all with the same truth-level isolation requirement as the data. Statistical uncertainties are shown as vertical bars, and the total systematic uncertainty as shaded bands. The boxes on the \jetphox{} points show the scale-variation uncertainty from $\mu_{R,F,f}\to\etg/2,2\etg$. The middle panel shows the theory-to-data ratio with the experimental systematic uncertainty drawn as a band around unity, and the lower panel shows the \jetphox{} prediction at \comHEP{} with four alternative proton PDF sets (\textsc{CT18NLO}~\cite{Hou:2019efy}, \textsc{NNPDF3.1}~\cite{NNPDF:2017mvq}, \textsc{CTEQ6.6}~\cite{Nadolsky:2008zw}, and \textsc{MSHT20NLO}~\cite{Bailey:2020ooq}) divided by the data. \pendingTODO{regenerate against the all-range allz nominal at $\mathscr{L}=64.37~\pb$ with extended $\etg$ range.}}
    \label{fig:final}
\end{figure}

\paragraph{Comparison with PHENIX.} The cross section is compared with the previous direct-photon measurement by the PHENIX experiment in Fig.~\ref{fig:final_phenix}~\cite{PHENIX:2012jgx}. The PHENIX measurement was performed without an isolation requirement, in a narrower pseudorapidity acceptance ($|\eta^\gamma|<0.25$), and with each data point reporting the cross section at the bin centre rather than the bin-integrated average. To enable a direct comparison, the PHENIX data are rescaled in this Article by the bin-width-recovery factor and by the ratio of $\eta$-density between the two acceptances, where the $\eta$-density ratio is taken from the prompt-photon truth distribution in the \pythia{} signal simulation rather than from a flat-$\eta$ ansatz. The two measurements agree within the quoted experimental and theoretical uncertainties, and the present result extends the previous measurement to higher \etg{} and to the full barrel rapidity coverage.

\begin{figure}[tbp!]
    \centering
    \includegraphics[width=0.75\linewidth]{figures/final_phenix.pdf}
    \caption{Comparison of the present measurement of the isolated-prompt-photon cross section with the PHENIX direct-photon measurement~\cite{PHENIX:2012jgx} (rescaled by bin-width recovery and by the $\eta$-density ratio for $|\eta^\gamma|<0.7$ taken from the prompt-photon truth distribution in the \pythia{} signal simulation). Statistical uncertainties are shown as vertical bars and total systematic uncertainties as shaded bands. \pendingTODO{regenerate against the all-range allz nominal; show both the rescaled PHENIX and the original (non-rescaled) PHENIX points for transparency.}}
    \label{fig:final_phenix}
\end{figure}
```

### 4.6 §6 Summary (replaces conf-note §6)

Intent: update lumi/range; add outlook on heavy-ion programme; quantify total uncertainty envelope.

**PROPOSED (replaces conf-note §6, lines 258–262):**
```latex
\section{Summary}
\label{sec:summary}
The differential cross section of isolated prompt-photon production as a function of \etg{} has been measured in \pp{} collisions at \comHEP{}, using the full \RunOne{} sPHENIX dataset corresponding to an integrated luminosity of $\mathscr{L}=64.37^{+5.88}_{-4.34}\,\pb$. Photons are measured in the kinematic range $|\etag|<0.7$ and $12<\etg<32~\GeV{}$, with a fiducial truth-level isolation requirement of $\isoETtruth<4~\GeV{}$ within $\DR=0.3$. Photons are reconstructed in the EMCal and identified with a multivariate analysis of the electromagnetic shower shape, and the residual background from neutral-meson decays is subtracted by a four-region (ABCD) sideband technique. The yield is corrected for purity and unfolded for detector response using the iterative Bayesian method, then converted to a cross section with reconstruction, identification, isolation, MBD-vertex, and trigger efficiency corrections.
The measured cross section spans approximately two orders of magnitude over the reported range and has a total relative uncertainty reaching $-30.4\%$ on the down side and $+35.1\%$ on the up side at the high-\etg{} end, dominated by the EMCal energy-scale uncertainty at high \etg{} and by the ABCD MC-closure uncertainty at low \etg{}.
The result is consistent with NLO pQCD predictions from \jetphox{} and from Vogelsang, both with the \textsc{CT18NLO} PDF set and with the same truth-level isolation requirement as the data, with the \pythia{}~8.307 Monte Carlo using the Detroit tune, and, where comparable, with the previous PHENIX direct-photon measurement at the same collision energy~\cite{PHENIX:2012jgx}. The present measurement extends the existing $\sqrt{s}=200$~\GeV{} prompt-photon coverage to higher \etg{} and to the full barrel acceptance, and constitutes the pp baseline for forthcoming sPHENIX measurements of photon and photon--jet observables in nuclear collisions, where the high-precision pp reference is essential for the interpretation of any in-medium modification.
```

### 4.7 New trailing matter (Acknowledgments, Data availability, Appendix A)

**PROPOSED (replaces conf-note appendix block, lines 270–280):**
```latex
\section*{Acknowledgments}
\pendingTODO{standard sPHENIX collaboration acknowledgements + RHIC + CADR + funding agencies; copy template from a recent sPHENIX publication.}

\section*{Data availability}
\pendingTODO{tabulated cross-section values are provided in Tab.~\ref{tab:xsec_values} and uploaded to HEPData (DOI to be assigned).}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bibliographystyle{unsrturl}
\bibliography{references}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\clearpage
\appendix
\section{Tabulated cross-section values}
\label{app:xsec_values}
The bin-by-bin cross section is reported in Table~\ref{tab:xsec_values}. Statistical and systematic uncertainties are listed separately. The systematic uncertainty includes all the contributions described in Sec.~\ref{sec:systematics} added in quadrature. The luminosity uncertainty is fully correlated across bins and is not included in the per-bin systematic to avoid double counting.

\begin{table}[h!]
    \centering
    \caption{Differential cross section of isolated prompt-photon production, $\mathrm{d}^{2}\sigma/(\mathrm{d}\etg\,\mathrm{d}\eta)$ in $|\etag|<0.7$, in nine \etg{} bins. Statistical (stat) and bin-by-bin systematic (syst) uncertainties are quoted, and the luminosity uncertainty of $+9.13\%/-6.75\%$ is fully correlated and is not included in syst.}
    \label{tab:xsec_values}
    \pendingTODO{populate from \texttt{Photon\_final\_bdt\_allz.root} after the systematic redo. Suggested columns: \etg{} bin range; bin centre; $\mathrm{d}^{2}\sigma/(\mathrm{d}\etg\,\mathrm{d}\eta)$~[pb/\GeV]; stat unc; +syst; -syst.}
\end{table}
```

---

## 5. Open items requiring user / collaboration input

### 5.1 Pending physics inputs (BLOCKING)
1. **Energy-scale shift** (2.6%) — placeholder pending final EMCal calibration. Confirm or update value with calorimeter group.
2. **Energy-resolution smear** (4%) — placeholder. Same dependency.
3. **MBD-vertex eff systematic** (±6%) — placeholder pending PPG10 / dedicated derivation. Memory `project_mbd_vertex_eff_pT_dependence.md` documents the photon-vs-jet cross-check at 1% level — may justify shrinking.
4. **1.5 mrad truth-vertex reweight** — fit did not converge per memory `project_truth_vertex_reweight.md` (closure χ²/ndf = 14.79). Decide journal-text treatment: re-fit, footnote-level caveat, or absorb into the DI systematic envelope.
5. **Full systematic redo** — pending per memory `project_syst_redo_needed.md` (tightbdt50/70, novtxcut, ntbdtmin05/10 retired; SYST_GROUPS["eff"] no longer contains tight_bdt). Final quadrature must be redone before regenerating `syst_sum_rel.pdf`.
6. **Photon-trigger plateau correction** — per memory `project_cut_efficiency_audit.md`, the lack of an explicit Photon-4-GeV L1 plateau correction biases the cross section low by ~0.9%. Confirm the per-bin Gumbel-CDF correction subsumes this, or apply an explicit plateau-pull factor.
7. **Bit-26 vs bit-30 reconciliation** — per memory `project_cut_efficiency_audit.md` and `project_trigger_bit30_turnon.md`. Analysis-note and code agree on bit 30 being the Photon-4-GeV trigger; bit 26 is legacy. The journal text should be self-consistent on the bit number, and the lumi-counting code should be cross-checked against the analysis selection bit before submission.
8. **Bit-30 trigger-eff fit form** — confirm the Gumbel-CDF $p_0\exp[-\exp(-(\etg-\mu)/\beta)]$ choice for the journal version, or move to a more standard erfc/sigmoid for reviewer familiarity.

### 5.2 Pending figures
| Figure | Status | Depends on syst redo | Action |
|---|---|---|---|
| `final.pdf` | Existing, lumi-obsolete | Yes | Regenerate via `plotting/plot_final_selection.C` against `Photon_final_bdt_allz.root` (or refreshed nominal). \jetphox{} and Vogelsang now both use CT18NLO with the truth-isolation requirement; the lower panel adds the four-PDF comparison (CT18NLO/NNPDF3.1/CTEQ6.6/MSHT20NLO). |
| `final_phenix.pdf` | Existing, lumi-obsolete | Yes | Regenerate. Show both rescaled and original PHENIX points. |
| `syst_sum_rel.pdf` | Existing, syst-list-obsolete | Yes | Rerun `calc_syst_bdt.py` after the syst redo. |
| `purity_nom.pdf` | Existing, lumi-obsolete | No | Regenerate against allz/all-range. |
| `eff_photon.pdf` | Existing, lumi-obsolete | No | Regenerate. Add MBD-vertex efficiency curve. |
| `h1D_iso_nom_6_8.pdf` | Existing, lumi-obsolete | No | Regenerate. Consider multi-panel by ETγ bin. |
| `dis_mixed_weta_cogx_eta0_pt2_cut1.pdf` | Existing, in `plotting/figures/showershape_1p5mrad/` | No | Copy to `PPG12-Paper/figures/` for the §3.3 shower-shape distribution figure. |
| `dis_mixed_e32_to_e35_eta0_pt2_cut1.pdf` | Existing, in `plotting/figures/showershape_1p5mrad/` | No | Copy to `PPG12-Paper/figures/` for the §3.3 shower-shape distribution figure. |

### 5.3 Pending tables
| Table | Required by | Action |
|---|---|---|
| `tab:datasample` | §3.1 | Run period, crossing-angle, integrated luminosity, event statistics for the two subsamples. |
| `tab:mc_samples` | §3.2 | MC generator settings, p̂_T,min ranges, generated event counts for the 8 SI/DI sample pairs. |
| `tab:photonid` | §3.3 | Photon-identification preselection cut values (parametric thresholds). Mirror analysis-note Table 3. |
| `tab:bdt_features` | §3.3 | 11 input features of the photon-ID BDT. Mirror analysis-note Table 4. |
| `tab:syst_summary` | §4 | Bin-by-bin systematic uncertainty breakdown. Mirror analysis-note Table 5. |
| `tab:xsec_values` | App. A | Bin-by-bin cross-section values for HEPData. Required by PRD. |

### 5.4 Pending citations (need to add to `references.bib`)
- `Hou:2019efy` — CT18NLO PDFs (nominal PDF set for both \jetphox{} and Vogelsang)
- `NNPDF:2017mvq` — NNPDF3.1 PDFs (lower-panel PDF comparison)
- `Nadolsky:2008zw` — CTEQ6.6 PDFs (lower-panel PDF comparison)
- `Bailey:2020ooq` — MSHT20NLO PDFs (lower-panel PDF comparison)
- `Owens:1986mp` — Owens prompt-photon review (LO production processes)
- `ATLAS:2017topoclusters` — ATLAS topo-cluster algorithm (used in §3.4 isolation lead-in; available in `PPG12-analysis-note/cite.bib`, copy across)
- `ATLAS:2019isophoton13TeV` — ATLAS 13 TeV isolated-photon (referenced in MC-purity-closure description; available in `PPG12-analysis-note/cite.bib`, copy across)
- D0:2005wir, CDF:2009gd or equivalent — Tevatron isolated-photon refs (confirm InspireHEP keys before insertion)
- STAR direct-photon papers (e.g. STAR @ 200 GeV and 510 GeV)
- Chen-Guestrin XGBoost paper for the BDT framework
- ATLAS:2017nah is **already** in `PPG12-Paper/references.bib` (line 194); reuse the existing entry
- `ATLAS:2016fta` is **already** in `PPG12-Paper/references.bib` (line 144); reuse the existing entry (referenced from §1 LHC-iso-photon list and §3.4 ABCD-method paragraph)

### 5.5 Pending text content (cross-linked to inline `\pendingTODO{}` markers)

This list is a checklist against the inline `\pendingTODO{}` markers in §3 and §4; use it as a tracking summary, not as a duplicate task list. Each entry below points to the inline TODO that owns the resolution.

- **EMCal in-situ energy-resolution parametrisation** → §4.3.3 cluster-reconstruction paragraph. Quote final $a/\sqrt{E\,[\GeV]}\oplus b$ once the analysis-note in-situ fit is final.
- **Acknowledgments** → §4.7 (acknowledgments TODO). Copy template from a recent sPHENIX publication.
- **Noise-overlay procedure** → §4.3.2 (noise procedure inline TODO). Confirm random-tower overlay vs. parametric Gaussian with the EMCal calibration team.
- **Photon-trigger plateau-correction systematic** → §4.4 (photon-trigger paragraph TODO). Confirm the per-bin Gumbel-CDF correction is the only treatment, otherwise propagate the $1\sigma$ fit-covariance band.
- **MC sample table** → §4.3.2 (Table~`tab:mc_samples`). Mirror analysis-note `table:MCdataset` (note's table is missing the `jet40`/`jet40_double` row; use the eight-pair list from `CrossSectionWeights.h` to populate the journal version) plus the truth-leading-pT stitching windows from `table:stitching`.
- **1.5 mrad truth-vertex reweight closure** — handled internally by the closure systematic per memory `project_truth_vertex_reweight.md`. The journal text now reduces the truth-vertex reweighting to a single descriptive sentence per user, and the residual non-closure is propagated through the post-pipeline closure cross-check rather than quoted inline (§4.3.2 + §4.4 closure paragraph).
- **MBD-timing-veto removal cross-check** — confirmed: no `mbd_avg_sigma_max` field in `config_bdt_nom*.yaml`. The journal text reflects this absence, and the per-bunch-crossing pileup effect is handled at the MC level by the SI+DI blending of §4.3.2 rather than an event-level data veto.
- **Shower-shape distribution figures** (NEW) → §4.3.3 (Fig.~`fig:showershape_dis`). Two-panel figure showing \wetacogx{} and $E_{3\times 2}/E_{3\times 5}$ at the NCB pre-selection level for a representative \etg{} bin. Files exist at `plotting/figures/showershape_1p5mrad/dis_mixed_*_eta0_pt2_cut1.pdf`. Copy into `PPG12-Paper/figures/`.

---

## 6. Style/voice fixes embedded in proposed text

- "this Conference Note" / "in this note" → "this Article" / "the present analysis" (8+ instances).
- "we define photons to be isolated" → "Photons are required to be isolated" (passive voice, journal style).
- "approximately 80% to 90%" hand-waving → flagged as TODO for actual binned values.
- PYTHIA version standardised to "8.307" (was inconsistent in the conf-note: "PYTHIA.307" in §3.2 vs "PYTHIA 8.307" in §5).
- Cleaned `\etg^{...}` constructs (illegal double-superscript with the `\etg` macro that already has `^\gamma`). Use `\etreco` / `\etgtruth` macros from `defs.sty` and prose alternatives.
- Semicolons in the proposed body have been swept and replaced with periods or commas+conjunctions, per project convention (memory `feedback_no_semicolons.md`).
- Outline paragraph added at end of §1 (referencing all subsequent section labels).
- Reframed §3.3 (Photon Reconstruction and Identification) to a "Cluster reconstruction → Pre-selection (NCB BDT + rectangular cuts) → Photon-ID BDT" ordering inside a single subsection, matching the analysis-note grouping where the non-collision-background veto is part of the pre-selection chain rather than a separate post-ID step.
- §2 Detector text reverted to the conf-note original (Pass 8) — the description is unified across sPHENIX collaboration papers.
- §3.1 Event Selection trimmed to conf-note granularity (Pass 8) — no per-paragraph splits, no trigger bit numbers, no internal-note cites, no dataset-summary table.
- §3.2 truth-vertex reweighting reduced to a single descriptive sentence (Pass 8) — no Padé / iterative-fit / 1.5 mrad non-closure detail in the journal text.
- §4.5 \jetphox{} and Vogelsang both at CT18NLO with the truth-isolation requirement (Pass 8); the final-plot lower panel adds the four-PDF comparison (CT18NLO/NNPDF3.1/CTEQ6.6/MSHT20NLO) per the live `plot_final_selection.C` third pad. The $x_T$-scaling paragraph and the $\chi^2/\mathrm{ndf}$ data/theory TODO are held back from the journal version for now.

---

## 7. Suggested workflow for the user

1. **Read this document first.** Skim the section deltas in §4. The journal target (PRD) and the section structure (Option A, single §3 Analysis Procedure) are already settled — see §1 and §2.
2. **Resolve §5.1 placeholders** (energy scale 2.6%, EMCal smear 4%, MBD-vertex eff ±6%) with the calorimeter/lumi/vertex-reweight groups. These gate the abstract's headline uncertainty.
3. **Once placeholders are resolved**, the proposed text can be applied to `PPG12-Paper/main.tex` in a single editing session. The file currently is at conf-note state (untouched).
4. **Independently**, the figures and tables flagged in §5.2 and §5.3 can be regenerated/populated. Several depend on the systematic redo finishing first.

---

## 8. Files reviewed (read-only, not modified)

- `PPG12-Paper/main.tex` — the conf-note source (unchanged)
- `PPG12-Paper/references.bib` — bibliography
- `PPG12-Paper/defs.sty` — LaTeX macros
- `PPG12-Paper/figures/` — 11 figures inventory
- `PPG12-analysis-note/` — current internal note (sPH-ppg-2024-012)
- `wiki/physics/` and `wiki/raw/papers/` — journal-target precedent + comparison-paper context
- `efficiencytool/results/Photon_final_bdt_allz.root` (key listing only) — current nominal output

## 9. Files NOT modified

`PPG12-Paper/main.tex` is at git HEAD (commit a3090a4). All proposed edits are in this document only.
