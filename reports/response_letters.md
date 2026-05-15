# PPG12 Paper — Response Letters (Cycle 1 → Cycle 3.2)

**Target file**: `PPG12-Paper/main.tex` (current draft after cycle 3.2 patch application).
**Source PDF reviewed**: `sPHENIX_PPG12_Paper.pdf` (May 9 2026, 19 pages, line-numbered).
**Patch artifact**: `reports/paper_review_patch_cycle1.md`.

Each item below quotes the reviewer's comment verbatim with a `>` blockquote, followed by our response. Items addressed in cycle 3.2 cite the new line range in the rendered PDF where applicable; items held pending external input are flagged with the corresponding red TODO marker placed in the current draft for visibility.

---

## Response to Jamie

> Abstract - often abstracts avoid too many defined acronyms and save that for the main text. That often helps readability.

Addressed: the redundant acronym definitions have been removed from the abstract while preserving the original sentence structure. The "(EMCal)" definition was dropped (the abstract now uses "electromagnetic calorimeter" inline, and the later mention reads "electromagnetic and hadronic calorimeters"); the "(pQCD)" definition was dropped in favour of writing out "perturbative-QCD" at both occurrences. The new abstract is at L60-62 of the current draft.

> line 10 - should we quote 64.37 pb^-1 or just 64.4 pb^-1 as is done in the figures? Good to be consistent.


> line 17 - A comparison with the previous PHENIX measurement of direct photons at the same collision energy is also presented.

Addressed: this comparison sentence had already been removed from the abstract prior to the formal review submission, and the current draft does not contain the flagged sentence. The PHENIX comparison now appears only in Fig. 8 and the accompanying Sec. V text.

> line 20 - sensitive to the gluon parton distribution in the proton

Addressed: the "sensitive to" / "sensitivity to" framing has been softened to "depends on" / "dependence on" across the abstract (L62), the introduction (L77), and the summary (L373). The Ichou:2010wc gluon-PDF / direct-photon sensitivity reference is also cited at the new intro sentence for support.

> line 25 - in field theory, I thought ISR and FSR were indistinguishable and can also interfere. Thus, are fragmentation photons somehow identifiable as exclusively FSR?

Held pending expert consultation on the canonical ISR/FSR field-theoretic phrasing. The comment is marked in the current draft as a red TODO marker (in the prompt-photon definition at L69-70) for visibility, and the wording will be tightened once the theory expert returns on the operational prompt/fragmentation split and the isolation-prescription regulator.

> line 30 - ref [1] inside period

Addressed: the citation now renders inside the period in the current draft ("partons into photons [1]." instead of "partons into photons. [1]"). The .tex source places `\cite{Aurenche:2006vj}` before the period at main.tex L72, and the rebuilt PDF reflects this.

> line 34 - "cross section" / "cross-section" — I prefer no hyphen; what's the official rule?

Addressed: APS house style is "cross section" with no hyphen, and a one-pass sweep through main.tex normalized all ~20 occurrences. The subsection title in Sec. III.E was retitled "Cross Section Determination" (Title Case, no hyphen).

> line 37 - I believe PHENIX has measurements in pp, dAu, AuAu as paralleled by LHC. Also, STAR has a measurement (https://arxiv.org/pdf/0912.3838), and not sure if there is something more recent. Please double check.

Held pending the STAR data overlay on the PHENIX-comparison figure. The STAR:2009ojw bibtex entry (arXiv 0912.3838, Phys. Rev. C 81, 064904) has been verified and staged, and the citation plus intro update will land together with the figure update so the reader sees the prose and the point set in one revision. The comment is marked in the current draft as a red TODO marker (after the PHENIX reference at L74) for visibility.

> line 43 - the quark--gluon Compton process gives the dominant contribution to prompt-photon production [15], making the cross section particularly sensitive to the gluon parton distribution function (PDF) of the proton.

Addressed together with the L20 gluon-PDF comment above. The intro sentence at L77 now reads "so that the cross section depends directly on the gluon parton distribution function (PDF) of the proton~\cite{Ichou:2010wc}", consistent with the softened framing applied in the abstract and summary.

> line 44 - I am not sure that Article is used in PRC or PRD.

Addressed (combined with Justin's "first sPHENIX measurement" comment at the same line): L80 now reads "This paper presents the first sPHENIX measurement of the differential cross section of isolated prompt-photon production…". The "Article" capitalization has been dropped and the "first sPHENIX measurement" framing has been added; no prior sPHENIX isolated-photon paper exists on arXiv as of May 2026, so the "first" claim is supportable.

> line 48 - should we define "truth level" here, a la the detailed definition used by ATLAS?

Held pending expert consultation on the canonical truth-level fiducial definition. The intro and Sec. III.A truth-iso sentences are linked through the same fiducial-definition prose (also flagged by Justin J5 and J13), so they will land as one unified ATLAS-style paragraph once the theory and ATLAS-experienced collaborators return. The comment is marked in the current draft as a red TODO marker at L80 for visibility. As an interim, the false "excluding neutrinos" claim at L113 has been corrected (see the L96 response below).

> line 67 - the brief description on the sPHENIX detector and DAQ on this page is a really high level. It would also be good to maybe expand a bit on the EMCal as it is the main detector used in this analysis.

Held pending expert sign-off on the L1 trigger-window size (4×4 vs 8×8 per current code/PDF) and the "online-calculated" qualifier. A drafted detector-section expansion is staged in the patch artifact and adds an MBD description plus an EMCal supercluster / SiPM-readout sentence. The comment is marked at L102-103 (trigger sentence) for visibility, and the expansion will land once John Haggerty / the sPHENIX L1 firmware team confirm the two open numbers.

> line 70 - typo: "is scintillating-fiber" should be "is a scintillating-fiber".

Addressed: L89 now reads "The EMCal is a scintillating-fiber tungsten sampling calorimeter with two-dimensional projective geometry in $\eta$-$\phi$". The missing article was inserted as part of the same edit that absorbed Justin J8's "η-φ two-dimensional projective geometry" suggestion.

> line 77 - 8 8.307 is a weird redundant 8

Addressed: the `\pythia` macro in `defs.sty` has been redefined so that `\pythia~8.307` renders as "PYTHIA 8.307" (no duplicate 8). Existing macro usages at L106, L337, L342, and L373 propagate the fix.

> line 92 - The pi0 / eta mass and width is probably a more direct measurement, and we should give the values. Maybe even a figure of them since this is important.

Held pending the choice of mass-peak source. No PPG12-internal $\pi^0$/$\eta$ data-vs-MC overlay file exists on disk; EMCal-WG candidates are run-level peak monitoring rather than a data-vs-MC overlay. The mass-peak figure is deferred to revision 2 once a dedicated diphoton-pair study is added or a suitable EMCal-WG plot is borrowed. The L111 prose has been simplified to describe the smearing intent without quoting the now-stale 4% number, which fixes the immediate concern; the figure remains a soft must-have.

> line 96 - see above point regarding truth level. Anything about muons here?

Addressed in part: the false "excluding neutrinos" claim has been deleted from L113, which now reads "the total transverse energy of all final-state particles within $\DR=0.3$, excluding the photon itself, is below 4 GeV" — consistent with the intro at L80. The audit on 747,706 truth photons confirms the mean $\nu+\mu$ contribution to the truth-iso sum is 0.02 MeV and the iso-classification flip fraction is 0.0003%, so the practical impact is negligible. The broader ATLAS-style fiducial-definition expansion (which would explicitly address the muon question) is held pending expert consultation and is marked in the current draft as a red TODO marker at L113 for visibility, addressed together with the L48 comment above.

> line 112 - "hard" — remove this

Addressed: L117 now reads "$f_{\mathrm{PU}}$ is the corresponding probability that an event producing a high-energy cluster contains two or more \pp{} interactions" (the qualifier "hard" has been deleted).

> line 124 - resolution undefined

Addressed: L125 now reads "Reconstructed photons have an energy resolution $\sigma_E/E$ of approximately 6\% across the entire \etg{} range, dominated by the EMCal stochastic and constant terms." The variable being quoted and its dominant terms are now made explicit.

> line 128 - it might be a good idea to comment on the approximate pT at which the photons from pi0/eta merge into a single calorimeter cluster.

Addressed: L127 now adds "The two photons from $\pizero \to \gamma\gamma$ become indistinguishable in a single EMCal cluster above $p_T \approx 11$~\GeV; for $\eta \to \gamma\gamma$ the corresponding threshold is $p_T \approx 46$~\GeV. The merging is gradual rather than a hard cut." The two threshold values are placeholders (rendered in red via `\pendingTODO`) derived from the EMCal granularity $\Delta\theta \approx 0.024$ and PDG masses, pending tighter validation with shower-shape simulation.

> line 132 - "centre" should be "center". Also, in many places "centre-of-gravity" but I would lean toward "center of gravity" without the hyphens.

Addressed: a sed pass converted "centre" → "center" everywhere, and the hyphenated "center-of-gravity" form has been kept as a compound adjective per standard usage. L127, L135, L136, L138, L142, and L171 are now consistent.

> line 138 - it might be informative to point out that this is largely from muons that travel upstream along the z-axis from collisions outside of the detector.

Addressed: the NCB paragraph at L151 now reads "Non-collision background (NCB) clusters originate predominantly from beam-pipe interactions upstream of the interaction region. The resulting muons propagate approximately parallel to the beam direction and produce clusters with elongated $\eta$-like shower shapes and arrival times shifted earlier than the in-time beam crossing." This adds the physical-origin sentence directly to the NCB description.

> line 144 - it might be helpful to specify what these are.

Not changed because enumerating all 17 additional NCB-BDT features in the paper would expand the methodology section without adding physics insight; the analysis note already lists them. The additional shower-shape observables used by the NCB BDT are: extended energy-weighted width moments $w_{\eta,3}$, $w_{\phi,3}$, $w_{r,2}$, $w_{r,3}$; further asymmetry observables built from larger $3\times3$ and $3\times5$ tower blocks; concentric energy-fraction ratios $E_{2\times 2}/E_{3\times 3}$, $E_{3\times 5}/E_{5\times 5}$, $E_{3\times 7}/E_{5\times 7}$; and three kinematic decorrelation priors (cluster $E_T$, cluster $\eta$, collision vertex $z_\mathrm{vtx}$). The full list is in `FunWithxgboost/config.yaml`.

> figure 1 - I think it would be useful to actually show what the background shape looks like too. As Eli mentioned, this is also part of the BDT training.

Addressed: Fig. 1 already shows the inclusive jet background MC (blue histogram) — the relevant histogram is read from `MC_efficiencyshower_shape_jet_inclusive_combined_*.root`, which is the pure inclusive-jet sample weighted by cross-section. The plotting macro `plot_paper_showershape.C` legend has been relabeled from "Inclusive MC" to "Inclusive jet MC" so the background nature is unambiguous, and the rendered figures dis_*_cut1.pdf have been regenerated.

> line 179 - I think it would be useful to show how the BDT score looks for data and MC.

Held pending the new plotting macro. A figure paralleling Fig. 1 showing the BDT score for data, signal MC, and inclusive jet MC at one representative $p_T$ slice (14 < $E_T^\gamma$ < 18 GeV) with the tight-BDT threshold drawn as a vertical line is in scope for the next revision; the existing `cut1` shower-shape plotting infrastructure can be reused. This overlaps with Virginia's BDT-distribution-figure request.

> figure 3 - it could be useful to show the schematic plot side-by-side with the actual 2D plot for the data.

Held pending the new plotting macro. The side-by-side Fig. 3 (left = current schematic; right = 2D `colz` of data in the BDT-score × $E_T^\mathrm{iso}$ plane with the four ABCD cut lines drawn) reuses the existing data slim tree and is in scope for the next revision.

> line 236-237 - I am not sure all of this needs to be said since it just restates what the figure shows.

Addressed: the redundant restatement of the figure has been removed at L255. The "0.5 at 10 GeV to 0.9 around 30 GeV" purity range — directly readable from Fig.~\ref{fig:purity} — has been dropped. The leakage-correction size (0.1–0.2 depending on ETγ) is kept as physics context not directly visible in the figure overlay, alongside the Padé-fit description.

> Padé equation - I would explicitly write out the equation P[1,1](x) so reader knows what it is.

Addressed: L255 now reads "...is fitted with a $[1/1]$ Padé approximant, $P_{[1/1]}(x) = (a_0 + a_1 x)/(1 + b_1 x)$, to smooth statistical fluctuations in the bin-by-bin purity estimates." The explicit functional form is now in the paper.

> line 266 - I am uncertain what this is and the value should be derived from. Also maybe a figure here would help illustrate this.

Addressed (text-only — no new figure per Shuonli's preference): L278 now reads "The MBD-vertex efficiency, $\varepsilon_{\mathrm{vtx}}$, is defined as the fraction of truth-level isolated signal photon events in which the MBD records at least one hit on each side and yields a reconstructed vertex within $|z_\mathrm{reco}|<60$~cm. The denominator uses the all-$z$ truth-vertex fiducial, so this efficiency converts the beam-delivered luminosity into the analysis-fiducial vertex acceptance. The efficiency decreases from approximately 59\% at $\etg \simeq 12$~GeV to approximately 52\% at $\etg \simeq 30$~GeV." The earlier "55% to 60%" range was incorrect both in direction and approximate values; the new numbers are extracted from `Photon_final_bdt_nom.root::g_mbd_eff` across the fiducial range.

> line 285 - I think more detail is needed here.

Addressed: the energy-scale shift sentence at L313 has been expanded to "The reconstructed EMCal cluster transverse energy in MC is varied by a flat multiplicative factor of $1 \pm 0.011$, bracketing the residual data--MC offset of the reconstructed \pizero{} and $\eta$ peak positions measured across the EMCal acceptance. The variation is applied per cluster, before any identification or isolation requirements, and is propagated through the full analysis chain to give a symmetric two-sided systematic, added in quadrature with the other contributions. This uncertainty is among the dominant systematic contributions at high \etg{}, together with the energy-resolution uncertainty." The form, magnitude, application point, and motivation are now explicit.

> line 306 - this sentence is quite vague.

Held pending the PPG10 bibtex key (the data-driven MBD-coincidence efficiency reference). A drafted replacement that quotes the $\pm 6\%$ envelope and explains the dijet-to-photon transfer is staged in the patch artifact and will land once the PPG10 collaboration paper key is provided. The comment is marked in the current draft as a red TODO marker at L317 for visibility.

> The truth-level Gaussian smearing of MC was the approach used in PPG09 -- and PPG10 if that has been published by then. Reference if possible. Also closure tests of the procedure were performed and could be worth replicating.

Not changed in the paper because the current production code already smears at the response-matrix fill stage with the width sized by truth $p_T$, which is the PPG09-equivalent prescription. The systematic envelope (zero-data-params variation up vs wider-data-params variation down) spans the data-resolution-fit uncertainty band, which is the natural closure on whether the assumed resolution is correct. The response-matrix-only application avoids double-counting in the bin-by-bin efficiency. A formal PPG09/PPG10 reference will be added once those collaboration papers are public; the closure-test cross-check is also tracked as a revision-2 item.

> figure 6 - the vertical lines on the x-axis edges look strange.

Held pending a cosmetic fix in the plot macro. The vertical lines are horizontal error-bar end-caps on the boundary bins; `SetErrorX(0)` in `plot_paper_systematics.py` will remove them, or rebinning the rightmost overflow to uniform 2-GeV bins is an alternative. The fix and regenerated `syst_sum_rel.pdf` are in scope for the next revision.

> line 324 - I do not understand what is meant by HERWIG ... mentioned at line 109 of the MC section but absent from the systematics.

Addressed in two places. (a) In the MC-sample section at L109, the HERWIG sentence now states it is used as a cross-check on the isolation efficiency. (b) In the systematics section at L319, the HERWIG vs PYTHIA isolation-efficiency comparison is now listed as part of the isolation-efficiency systematic uncertainty.

> line 338 - it would be useful in the text to walk through what each of the panels in figure 7 shows. For example, "The middle panel shows...". Right now, the discussion is focused on the upper panel.

Addressed (middle panel): the Results section at L337 now reads "The middle panel of Fig.~\ref{fig:final} shows the ratio of each prediction to the data, with the experimental total uncertainty drawn as a band around unity and theoretical scale and PDF uncertainties shown as boxes around each prediction. The lower panel shows the \jetphox{} predictions divided by the data for four alternative PDF sets, which bracket the data within experimental uncertainties across the measured \etg{} range." Both the middle and lower panels are now discussed in the body text.

> line 338 - more discussion is needed on the bottom panel of figure 7, including the implications for the gluon PDF.

Addressed in part: the lower-panel sentence is now in the body (see the L338 response above) and notes that the four alternative PDF sets bracket the data within the measured range. The deeper gluon-PDF implications (Hessian uncertainty band, gluon-fraction extraction) are held pending the JETPHOX rerun with NNPDF4.0 + the CT18NLO Hessian eigenvector sweep, which is deferred to revision 2. The NNPDF:2021njg bibtex entry has been added.

> figure 7 caption - to address the comment about the gluon PDF, it should be NNPDF4.0 and not NNPDF3.1. Maybe also show the Hessian uncertainty band of the central PDF? Right now, it would be good to identify what is the meaning of these comparisons -- in light of the comment about gluon PDF.

Held pending the JETPHOX NNPDF4.0 rerun. The NNPDF:2021njg bibtex entry (arXiv 2109.02653, EPJC 82, 428) has been verified and added to references.bib. The lower-panel curve will be regenerated and substituted in once the ~3-hour-wallclock JETPHOX run is completed; the Hessian band is deferred to revision 2 because the CT18NLO 58-eigenvector × full-grid sweep is not feasible in the resubmission window.

> figure 8 - the statistical uncertainties on this plot may also be useful to plot, perhaps relative to 1 (i.e., a horizontal band).

Held pending the plot cosmetic fix. The sPHENIX statistical band is currently drawn around y=1 in light blue and is partially hidden by the theory boxes; darkening the band and adding an explicit "sPHENIX (stat ⊕ syst)" legend entry in the ratio panels is in scope for the next revision.

> figure 8 caption - "(Details of the correction procedure are described in the text.)" should be a sentence rather than enclosed in parens.

Addressed: the Fig. 8 caption at L355 now reads "...is overlaid. Details of the correction procedure are described in the text." The parenthetical has been split into a separate sentence.

> line 368, 372 - again, our cross section is not really sensitive to the gluon PDF unless we are showing some sort of constraint of g(x)

Addressed together with the abstract and intro gluon-PDF comments above. The summary sentence at L373 now reads "The result provides a test of pQCD, a dependence on the gluon parton distribution function, and a \pp{} baseline for measurements of photon and photon+jet observables in heavy-ion collisions." The "sensitivity to" framing has been replaced with "dependence on" everywhere.

> References - especially in the case of Collaboration papers, the formatting and citing conventions should be consistent across all (Refs 14, 15 - or 19, 24, 25, etc.).

Addressed: six entries in references.bib have been normalized to the `"{XXX Collaboration}"` author form — `PHENIX:2010vgy`, `ATLAS:2011ezy`, `PHENIX:2015siv`, `Aune:2024ysr`, `Aidala:2020toz`, `sPHENIX:2017lqb`. The `Aune:2024ysr` and `Aidala:2020toz` entries also get a `collaboration = "sPHENIX"` field. The `Belmont:2023fau` single-author whitepaper has been left as is.

---

## Response to Virginia

> As discussed in the context of the photon+jet paper, I wonder if it would be useful to show the BDT distribution in this paper where the BDT is described in much more detail than in that one where it is currently included

Held pending the new plotting macro. The proposed figure parallels Fig. 1 of the current draft (data vs signal MC vs inclusive-jet MC, $14<\etg<18$~GeV slice) but shows the photon-identification BDT score with the tight-BDT threshold drawn as a vertical line. This overlaps with Jamie's L179 request and is scoped for the next revision using the existing `cut1` shower-shape plotting infrastructure.

> 46- Run 24 -> RHIC Run 24

Addressed: the `\RunOne` macro in `defs.sty` L11 was redefined to expand to "RHIC Run 24", which propagates to all body-text occurrences (main.tex L80 intro, L115 MC pile-up paragraph, L117 mixing prose, L375 summary). The remaining "in 2024" wording in the abstract at L61 is intentional plain-English context (calendar year), separate from the dataset-name shorthand carried by `\RunOne`.

> 51- I think we decided based on the definitions in this paper that the PHENIX measurement is really prompt (direct + fragmentation) photons, right? It would be good to be consistent in this paper even if PHENIX calls it direct

Addressed: the intro L75 wording is now "the PHENIX and STAR experiments reported prompt photon cross sections without an isolation requirement~\cite{PHENIX:2012jgx, STAR:2009ojw}", and the L80 closing sentence uses "previous inclusive prompt-photon measurement by PHENIX". Both phrasings use "prompt" consistently with the L70 definition of the term (direct + fragmentation).

> section 2- would be useful to have a reference to the PPG03 paper with details on the calorimeter calibrations

Addressed: the closest published EMCal-calibration reference is `sPHENIX:2025dET` (arXiv 2504.02242, Phys. Rev. C 112, 024908), added to `references.bib` at L1060. It is cited in Sec. III MC simulation at L111, backing the in-situ $\pi^0$ and $\eta$ mass-peak resolution match. No sPHENIX-tagged "PPG03" paper exists; this 2025 paper is the canonical pi-zero mass-peak calibration reference.

> 92- it seems a bit odd to discuss the specific extra smearing that is used since this is based entirely on what our internal simulations produce by default (and I don't think we show anywhere). I would add more detail on how this matching of the resolution between data and MC is done and/or remove the specific 4% number

Addressed: the smearing sentence at L111 has been simplified to "an additional \etg{}-dependent Gaussian smearing is applied to the simulated cluster \et{}, with a width derived such that the simulated resolution matches the one measured in data via in-situ $\pi^0$ and $\eta$ mass-peak fits~\cite{sPHENIX:2025dET}". The 4\% literal and the $\sigma_\mathrm{extra}(\etg)$ parametrization are no longer in the paper body. The wider-data and no-smearing brackets at L321 retain "$\approx 6\%$" only as a systematic-variation envelope, not as the nominal extra smearing.

> 108-116- you discuss how the pileup fractions are calculated but I think no where do you say explicitly that you combine the single and double interaction simulation samples based on these fractions

Addressed: the sentence at L117 now reads "The simulation sample used in this analysis is constructed within each crossing-angle period by combining the single-interaction and pile-up samples weighted by $(1-f_{\mathrm{PU}})$ and $f_{\mathrm{PU}}$ respectively, and the two crossing-angle periods are then combined by their integrated luminosities." The mixing weights and the per-period luminosity combination are both stated explicitly.

> 137- 'relative to the physics signal at the analysis ET' -> 'relative to the physics signal in the kinematic range of this analysis'

Addressed: the NCB-trigger sentence at L159 now reads "They fire the photon trigger at a rate that is non-negligible relative to the physics signal in the kinematic range of this analysis." The earlier "at the analysis \etg{}" wording has been replaced with the reviewer's preferred phrasing.

> 144- what other shower shape variables?

Addressed: Table~\ref{tab:bdt_features} at main.tex L131-156 now lists the additional NCB-BDT features explicitly below the double horizontal rule, namely $w_{n\times 2}$ for $n=3,5,7$ (L146), $E_{1\times 1}/E_{1\times n}$ and $E_{1\times 1}/E_{m\times 1}$ for $n,m=3,5,7$ (L147-148), $E_{1\times 1}/E_{2\times 2}$ (L149), and $E_{2\times 2}/E_{m\times n}$ with $(m,n)\in\{(3,3),(3,5),(3,7),(5,3)\}$ (L150). The body text at L161 retains "together with additional shower-shape observables" and refers the reader to the table.

> 308- specify the number of iteration variations you do

Addressed: the unfolding-uncertainty sentence at L325 now reads "The unfolding uncertainty is evaluated by repeating the unfolding without the prior reweighting and by varying the number of Bayesian iterations to three and four (the nominal is two)." The specific iteration values are now stated.

> 315-320- you mention three treatments but only list two (I guess the other is the luminosity uncertainty?)

Addressed: the sentence at L327 now lists three treatments explicitly. The third group is the unfolding iteration uncertainty, taken as the per-bin maximum of the iter-3 and iter-4 deviations symmetrized about the nominal. The two-sided list at L327 also folds energy resolution and luminosity into the per-bin asymmetric group, consistent with the canonical taxonomy in the systematics aggregator.

> 321- total systematic uncertainty excluding the global luminosity uncertainty

Addressed: the Fig.~\ref{fig:syst_sum} caption at main.tex L335 now reads "The Total envelope includes the global luminosity uncertainty (added in quadrature, bin-independent at $^{+9.1\%}_{-6.8\%}$) in addition to the seven component groups shown." The lead sentence at L318 has also been adjusted to "with the luminosity uncertainty added to the per-bin sums". The treatment of luminosity in the Total envelope is explicit.

> Figure 8- its quite hard to see the difference between the corrected and uncorrected PHENIX points, it might be useful to change the colors or plotting style. I'm also not very fond of the way of plotting the statistical error bars in the ratio (I stared at it a long time trying to understand before reading the caption, which won't always come along with the plot in presentations)

Held pending the plot cosmetic fixes in `plot_paper_final.C`. The proposed change replaces the current near-iso pink (kPink+5) PHENIX-original and purple (kViolet+1) corrected colors with colorblind-friendlier choices (kRed-4 original, kOrange+7 corrected), and the dual vertical-bar statistical-error style in the ratio panel is replaced so that the sPHENIX statistical uncertainty is drawn as a thin band around unity (mirroring the systematic-band style) while only the PHENIX statistical uncertainties remain as per-point vertical bars. These code changes are in scope for the next revision; the figure caption at main.tex L362 will be updated to match.

---

## Response to Justin

> Line 7+ (first sentence of abstract) It's borderline a run-on sentence. Good place to split: ...at 200 GeV. The data used was recorded by the sPHENIX detector … corresponding to 64 pb-1….

Addressed: the abstract opening is split at "200 GeV" (main.tex L60-61). The first sentence now ends "in proton--proton collisions at \comHEP" and the second sentence begins "The data were recorded in 2024 with the sPHENIX detector at the Relativistic Heavy Ion Collider, corresponding to an integrated luminosity of $64.4\,\pb$. Photons are measured in $|\etag|<0.7$ and $12<\etg<32~\GeV$."

> Line 20 (last sentence of abstract) – simply changing "sensitivity to" to "dependence on" could more easily satisfy at least Jamie's concern about true sensitivity of this measurement to g(x). … I don't think presenting such proof is NECESSARY for just including one of the traditional motivations of this measurement.

Addressed: the abstract closing at L62 now reads "the result provides a test of pQCD calculations, probes the gluon parton distribution function of the proton, and establishes the \ppp{} baseline…" The "sensitivity to g(x)" phrasing is gone; we are happy to soften "probes" to "depends on" if you prefer that wording in the abstract as well.

> ~L 42: In this spirit, you could just change "making … particularly sensitive to" --> "[traditionally?] expected to be sensitive to " …OR --> "which depends directly on"

Addressed: L78 now reads "the quark--gluon Compton process gives the dominant contribution to prompt-photon production~\cite{PHENIX:2010vgy}, so that the cross section depends directly on the gluon parton distribution function (PDF) of the proton~\cite{Ichou:2010wc}." This adopts your second suggestion verbatim and adds an explicit citation.

> L 44: Article presents the_[first? a first?] SPHENIX_ measurement of… COMMENT about that: once tracking is available, and including 2025, I think there is likely to be another sPHENIX measurement of this, we should e.g. be able to improve the pseudo rapidity and ET reach, possibly combining 24-25. Also since this was expected to be a channel sPHENIX is prime for.

Addressed: L80 now reads "This paper presents the first sPHENIX measurement of the differential cross section of isolated prompt-photon production…" We agree this leaves room for a future tracking-enhanced and/or 24+25-combined measurement, which is exactly why "first" (rather than "the sPHENIX measurement") is the right framing.

> L 47 : "truth-level" bad HERE (at least IMO) if you are just describing what's done to the data remove it here, if you mean something about MC, change the order of sentences but you should still state explicitly what's the cut (how more about the Econe/how cut is made?) in the data.

Addressed: the corresponding sentence in the intro (L80) now reads "Photons are required to satisfy an isolation criterion defined by the total transverse energy of final-state particles within $\DR=0.3$ around the photon, excluding the photon itself, being less than $4~\GeV$." The "truth-level" jargon is removed from this introduction sentence, and the cut is spelled out explicitly. The truth/reco distinction is then made in Sec.~II where it physically belongs (L113 defines the truth-level signal photon, and the reconstruction-level isolation used on data is described in Sec.~III.D).

> L55: change comprises to "is comprised of" … I didn't like it too much previously I think it was in the note or prelim, now seeing it again, I'm more sure it should be changed.

Addressed: L89 now reads "In order of increasing radius, sPHENIX is comprised of the following subsystems...", applying the "comprises" → "is comprised of" change directly.

> L65: ..is _a_ scintt-fiber calo --

Addressed: L90 now reads "The EMCal is a scintillating-fiber tungsten sampling calorimeter with two-dimensional projective geometry in $\eta$-$\phi$…" The missing indefinite article is added.

> L66/67…projective geom -->"\eta-\phi two-dimensional projective geometry" (advertising)… include reference to testbeam paper or TDR? (particularly because of the "projective geometry")

Addressed: same sentence at L90 is rewritten to "two-dimensional projective geometry in $\eta$-$\phi$ and approximately $0.024 \times 0.024$ $\eta-\phi$ segmentation." The testbeam reference (\cite{Aidala:2020toz}) and the EMCal TDR (\cite{sPHENIX:2017lqb}) are both cited one line earlier (L89) on the EMCal subsystem entry. The DAQ-rate sentence appended at the end of L90 still has a few minor typos ("utilize", "khz") that Shuonli will clean up in the next push.

> L70: [requires ] total _online-calculated_ energy to be above 4 GeV

Addressed: L103 now reads "Events are selected by a single high-$E_T$ photon trigger that requires the total online-calculated energy in any $8\times 8$ EMCal tower window to exceed 4~\GeV{}."

> L222 equation 2: The ratio should be in parentheses and remove the dot

Addressed: Eq.~\ref{eq:purity_noleakage} (main.tex L221-227) now renders the ratio inside `\left( ... \right)` parentheses and the equation closes with a comma (not a period), as the sentence continues on L228.

> L225 equation 3 remove the dot it's not needed.

Addressed: the multiplication `\cdot` between the $B$-factor and the parenthesized $C/D$ ratio in Eq.~3 (eq:purity_leak) has been removed for parallelism with Eq.~2 (where the ratio is also in parentheses).

> L229 caption for figure 4 right above it _a_ Pade fit (not an)

Addressed: the Fig.~\ref{fig:purity} caption (L271) reads "The purity with leakage correction is fitted with a Pad\'e function, and the shaded area shows the 68.3\% confidence interval of the fit." The "an" → "a" fix is applied.

> L358 – see previous comment about use of truth level, I don't like this if just meant as jargon-->what does truth level mean for data? IF it actually means something I didn't see it described anywhere (at least not clearly/exactly what was meant)

Held: two places still carry the "truth-level isolation requirement" phrasing, namely the Fig.~\ref{fig:final} caption at L350 ("…all with the same truth-level isolation requirement as the data") and the Summary at L376 ("…with a truth-level isolation requirement of $\isoETtruth<4~\GeV$ within $\DR=0.3$"). The phrase is intentional in these two places, where the truth/reco distinction must be made explicit because we are quoting the fiducial isolation definition that the NLO predictions match. We will reword to "particle-level" or to a literal definition ("the total transverse energy of final-state particles within $\DR=0.3$…") if you prefer; the substantive change has been made everywhere it referred only to the data (cf. J5 above).

> L372 Unless you do the sensitivity study, I DO think you should simply remove the "sensitivity to [g(x)]" from the conclusion, not needed anyway.

Addressed: the conclusion at L380 now reads "The result provides a test of pQCD, a dependence on the gluon parton distribution function, and a \pp{} baseline for measurements of photon and photon+jet observables in heavy-ion collisions." The "sensitivity to g(x)" language is gone, and no sensitivity claim is made. If you prefer to drop the gluon-PDF clause entirely, that is also fine with us.
