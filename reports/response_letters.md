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

Addressed: a generator-modelling systematic is now listed explicitly in §IV. The main-sources list at L311 adds "generator modelling"; the efficiency-uncertainty paragraph at L319 includes the sentence "A generator-modelling envelope on the photon-identification and isolation efficiencies is taken from the comparison of HERWIG and PYTHIA cross-checks (Sec.~III.A)"; the three-treatments classification at L320 lists "generator modelling" among the one-sided symmetrized sources. The MC-section sentence at L109 pointing the reader to Sec.~IV now resolves to a concrete mention in §IV. The 5.83\% per-bin envelope used in the budget is carried by the analysis note (`systematics.tex` L191-194, L305).

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

> L51: PHENIX measurement is really prompt (direct + fragmentation), not "direct" — be consistent across the document.

Addressed: the previously inconsistent "direct-photon measurement by PHENIX" phrasing at L79 has been changed to "previous inclusive prompt-photon measurement by PHENIX", consistent with the L75 wording. The PHENIX paper arXiv 1205.5533 defines "direct photons" to include both quark-line bremsstrahlung and fragmentation contributions, which corresponds to "prompt" in the PPG12 terminology.

> Section 2 - reference to PPG03 paper with calorimeter calibration details would help.

Addressed: a new bibtex entry `sPHENIX:2025dET` (arXiv 2504.02242, Phys. Rev. C 112, 024908) has been added to references.bib for the sPHENIX measurement that establishes the EMCal energy-scale calibration via iterative $\pi^0$-mass fitting. The citation is inserted in Sec. III "MC simulation" near the energy-resolution smearing sentence. No paper with the literal "PPG03" tag exists; this is the closest published EMCal-calibration reference.

> L46: "Run 24" → "RHIC Run 24"

Addressed: the `\RunOne` macro in `defs.sty` has been redefined to expand to "RHIC Run 24". All occurrences in the body text (intro L80, MC section L114, summary L368) update automatically through the macro.

> L92: extra smearing 4% specific number odd; add detail on data/MC resolution matching or remove specific 4%.

Addressed together with Jamie's L92 comment: the 4% specific number has been removed in favor of a description of the physics intent. L110 now reads "To account for residual differences between data and simulation in the single-particle EMCal energy resolution, an additional \etg{}-dependent Gaussian smearing is applied to the simulated cluster \et{} when filling the unfolding response matrix, with a width derived such that the simulated single-particle resolution matches the resolution measured in data. The unsmeared cluster \et{} is used for the photon-identification and isolation selections so that the additional resolution affects only the unfolding response." The detailed $\sigma_\mathrm{extra}(p_T)$ parametrization is documented in the analysis note.

> L108-116: pileup fractions are calculated, but there is no explicit statement that you combine SI + DI samples using these fractions.

Addressed: L116 now ends with "The simulation sample used in this analysis is constructed within each crossing-angle period by combining the single-interaction and pile-up samples weighted by $(1-f_{\mathrm{PU}})$ and $f_{\mathrm{PU}}$ respectively, and the two crossing-angle periods are then combined by their integrated luminosities." The mixing weights and the per-period combination are now both explicit.

> L137: "relative to the physics signal at the analysis ET" → "relative to the physics signal in the kinematic range of this analysis"

Addressed: L151 now reads "They fire the photon trigger at a rate that is non-negligible relative to the physics signal in the kinematic range of this analysis." The earlier "at the analysis $\etg{}$" wording has been replaced with the reviewer's preferred phrasing.

> L144: what other shower-shape variables (same comment as Jamie)?

Addressed together with Jamie's L144 comment. The full list of 17 additional NCB-BDT features is given in the Jamie response above and in `FunWithxgboost/config.yaml`; the paper retains the compact "together with additional shower-shape observables" wording to keep the methodology section focused.

> BDT distribution figure — useful here vs the photon+jet paper?

Held pending the new plotting macro. The BDT score data-vs-MC distribution figure overlaps with Jamie's L179 request and will land as one new figure paralleling Fig. 1, with the tight-BDT threshold drawn as a vertical line at one representative $p_T$ slice.

> L308: specify the number of iteration variations you do.

Addressed: L318 now reads "The unfolding uncertainty is evaluated by repeating the unfolding without the prior reweighting and by varying the number of Bayesian iterations to three and four (the nominal is two)." The specific iteration values used in the variations are now stated.

> L315-320: three treatments mentioned, but only two listed. Is luminosity the third?

Addressed: the L323 sentence now lists three treatments explicitly, with the unfolding iteration uncertainty as the third (max-symmetrized) group. L320 reads "Each variation falls into one of three treatments. Two-sided sources (energy scale, energy resolution, non-isolation boundary, NCB cut, purity 68\% confidence interval, luminosity) yield asymmetric up and down deviations that are kept separate per \etg{} bin. One-sided sources (tight \gid{} threshold, non-tight \gid{} threshold, fit functional form, MC closure correction, isolation efficiency, unfolding reweighting, pile-up fraction) yield a single deviation that is symmetrized about the nominal. The unfolding iteration uncertainty is taken as the per-bin maximum of the iter-3 and iter-4 deviations, symmetrized about the nominal." Energy resolution has been moved into the two-sided per-bin group, consistent with the canonical taxonomy in `systematics.tex` L21.

> L321: is the total systematic uncertainty excluding the global luminosity uncertainty?

Addressed: the Fig. 6 caption at L327 now reads "Breakdown of relative systematic uncertainties as a function of \etg. The Total envelope includes the global luminosity uncertainty (added in quadrature, bin-independent at $^{+9.1\%}_{-6.8\%}$) in addition to the seven component groups shown." The treatment of the luminosity in the Total envelope is now explicit.

> Figure 8 — corrected vs uncorrected PHENIX colors hard to distinguish; statistical error-bar style hard to read.

Held pending the plot cosmetic fixes. The proposed change replaces the current pink (kPink+5) PHENIX-original color and purple (kViolet+1) corrected color with colorblind-friendlier choices (kRed-4 for original, kOrange+7 for corrected), and the dual vertical-bar stat-error style will be changed so that the sPHENIX statistical uncertainty is drawn as a thin band around unity (mirroring the systematic-band style) and only PHENIX statistical uncertainties remain as per-point vertical bars. These code changes in `plot_paper_final.C` are in scope for the next revision.

---

## Response to Justin

> First sentence of the abstract is borderline run-on. Suggest splitting after "200 GeV".

Addressed: the abstract has been rewritten. The first sentence now reads "We report a measurement of the differential cross section of isolated prompt photon production at midrapidity in $p+p$ collisions at $\sqrt{s}=200$ GeV." and a second sentence carries the dataset and luminosity information. The new abstract is at L60-62 of the current draft.

> Last sentence of the abstract: "sensitive to" → "dependence on" the gluon PDF, unless a sensitivity study has been performed.

Addressed: the abstract closing sentence now reads "...and depends on the gluon parton distribution function of the proton in the moderate-$x$ region." The "sensitive to" framing has been replaced with "depends on" throughout (abstract, intro, summary), consistent with the position that we cite the gluon-PDF / direct-photon literature without claiming a sensitivity-study constraint.

> Intro: "making the cross section particularly sensitive to" → softer language, e.g., "traditionally expected to be sensitive to" or "which depends directly on".

Addressed together with the abstract comment above: L77 now reads "so that the cross section depends directly on the gluon parton distribution function (PDF) of the proton~\cite{Ichou:2010wc}." The Ichou:2010wc reference supports the gluon-PDF / direct-photon sensitivity statement without claiming an in-paper constraint.

> L44: "Article presents the [first? a first?] sPHENIX measurement of..."

Addressed together with Jamie's L44 comment: L80 now reads "This paper presents the first sPHENIX measurement of the differential cross section of isolated prompt-photon production…". The "Article" capitalization has been replaced with "paper" and the "first sPHENIX measurement" framing has been added. No prior sPHENIX isolated-photon paper exists on arXiv as of May 2026.

> L47: "truth-level" bad here if just describing what's done to the data; if MC, change the word order; describe E_cone / cut on data.

Held pending expert consultation on the canonical truth-level fiducial definition, addressed together with Jamie's L48 comment. The intro, Sec. III.A, and the Results-section JETPHOX-comparison sentence are linked through the same fiducial-definition prose and will land as one unified ATLAS-style paragraph once the theory and ATLAS-experienced collaborators return. The comment is marked in the current draft as a red TODO marker at L80 for visibility.

> L55: "comprises" → "is comprised of"

Held pending Shuonli's preferred wording. "Comprises" used actively (whole-comprises-parts) is correct per Garner's and Strunk & White, while "is comprised of" is a common but grammatically debated form; "consists of" is unambiguous and would satisfy the intent. The comment is marked in the current draft as a red TODO marker at L88 for visibility, and the wording will be finalized in cycle 3.3.

> L65: "is a scintillating-fiber calo" — insert "a"

Addressed: L89 now reads "The EMCal is a scintillating-fiber tungsten sampling calorimeter with two-dimensional projective geometry in $\eta$-$\phi$ and approximately $0.024 \times 0.024$ $\eta$-$\phi$ segmentation. The towers are grouped into $2\times 2$ superblocks read out by silicon photomultipliers." The missing article was inserted as part of the same edit that absorbed Justin J8 below.

> L66/67: "projective geom" → "η-φ two-dimensional projective geometry"; include testbeam / TDR reference for the EMCal.

Addressed together with J7 above: L89 now uses "two-dimensional projective geometry in $\eta$-$\phi$" and cites `Aidala:2020toz` (the 2-D projective sPHENIX EMCal prototype paper, arXiv 2003.13685) plus `sPHENIX:2017lqb`. Both references are already in references.bib.

> L70: "[requires] total online-calculated energy to be above 4 GeV"

Held pending sPHENIX L1 firmware-team confirmation of the "online-calculated" qualifier. The phrase has been inserted at L102 as "Events are selected by a single high-$E_T$ photon trigger that requires the total online-calculated energy in any $8\times 8$ EMCal tower window to exceed 4 GeV", and the comment is marked in the current draft as a red TODO marker for visibility. The 4×4 vs 8×8 tower-window size is also flagged for John Haggerty to confirm before final submission.

> L222 equation 2: ratio should be in parentheses and remove the dot.

Addressed: Eq. 2 (eq:purity_noleakage) at L213-219 now reads
$$
N^{A}_{\text{signal}} = N^{A}_{\text{raw}} - N^{B}_{\text{raw}} \left( \frac{N^{C}_{\text{raw}}}{N^{D}_{\text{raw}}} \right),
$$
The $N^C/N^D$ ratio is wrapped in parentheses and the terminating period has been replaced with a comma.

> L225 equation 3: remove the dot, it's not needed.

Held pending interpretation. The current Eq. 3 (eq:purity_leak) already ends with a comma rather than a period, so the literal "dot" at the end of the equation has already been removed. If the reviewer's "dot" refers instead to the `\cdot` multiplication symbol between the $B$-factor and the parenthesized $C/D$ ratio inside the equation, that can also be deleted for parallelism with Eq. 2 (now that the ratio is in parentheses). The comment is marked in the current draft as a red TODO marker at L222 for visibility, and the intent will be confirmed with the reviewer in the next round.

> L229 / Figure 4 caption: "an Pade fit" → "a Pade fit"

Addressed: the Fig. 4 caption at L264 now reads "The purity with leakage correction is fitted with a Padé function, and the shaded area shows the 68.3\% confidence interval of the fit." The wording has been changed from "an Padé fit" to "a Padé function" — both fixes the article and replaces the duplicated "fit" word.

> L358: truth level — if jargon, what does it mean for data?

Held pending expert consultation, addressed together with Jamie's L48 / J5 truth-level comments. The intro and Sec. III.A truth-iso sentences and the Results-section JETPHOX-comparison sentence are linked through the same fiducial-definition prose, and the comment is marked in the current draft as a red TODO marker at L337 (the JETPHOX-comparison sentence in the Results section) for visibility.

> L372: "sensitivity to the gluon parton distribution function" in the conclusion — remove unless a sensitivity study has been done.

Addressed together with the abstract and intro gluon-PDF comments. The summary sentence at L373 now reads "The result provides a test of pQCD, a dependence on the gluon parton distribution function, and a \pp{} baseline for measurements of photon and photon+jet observables in heavy-ion collisions." The "sensitivity to" claim has been softened to "a dependence on" everywhere; the deeper gluon-PDF constraint study (Hessian band, gluon-fraction extraction) is deferred to revision 2.
