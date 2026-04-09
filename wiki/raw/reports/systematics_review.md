# PPG12 Systematic Uncertainties Review

Review of PPG12's systematic uncertainty budget against standard experimental practice (ATLAS, CMS, ALICE, PHENIX). Based on `systematics.tex`, `SYST_PIPELINE_README.md`, `make_bdt_variations.py`, and the physics wiki.

---

## 1. Systematic Categories: PPG12 vs Other Experiments

### PPG12 current categories (from `systematics.tex`)

| Category | Subsections | Typical size |
|----------|-------------|-------------|
| Energy scale | 2.6% cluster ET shift | 8--25% |
| Energy resolution | 5% extra smearing | < 5% |
| Efficiency (tight ID) | BDT threshold 0.70/0.90 vs nom 0.80 | 5--10% |
| Purity: non-tight selection | nt_bdt_min = 0.05, 0.10 | up to 15% |
| Purity: sideband isolation | noniso shift +0.1/+1.0 | 4--20% |
| Purity: NPB score cut | npb = 0.3/0.7 vs nom 0.5 | documented but no figure |
| Purity: MC closure | ABCD closure correction | one-sided |
| Purity: fitting function | Pade vs erf + fit parameter CI | 10--20% |
| Unfolding | response matrix reweighting on/off | < 10% |
| Luminosity | MBD cross-section 25.2 mb (+9.1%/-6.7%) | asymmetric, flat |
| MBD trigger efficiency | +/- 5% variation | < 10% |

### ATLAS standard categories (from 1908.02746, 13 TeV)

| Category | PPG12 equivalent | Status in PPG12 |
|----------|------------------|-----------------|
| Photon energy scale (76-parameter model) | Energy scale (2.6% shift) | Covered (simpler model) |
| Photon energy resolution | Energy resolution (5% smearing) | Covered |
| Photon identification efficiency | Tight ID variation | Covered |
| Isolation efficiency | Implicit in purity variations | Partially covered |
| Luminosity (~2%) | Luminosity (+9.1%/-6.7%) | Covered |
| Purity / background estimation | Purity group (5 sub-sources) | Covered |
| R_bg background correlation | **Not explicitly included** | **Gap** |
| Trigger efficiency | MBD trigger efficiency | Covered |
| Material budget (conversion probability) | **Not included** | **Gap** (see below) |
| Signal leakage MC dependence | MC purity closure | Partially covered |

### PHENIX standard categories (from 1205.5533, pp 200 GeV)

| Category | PPG12 equivalent | Status |
|----------|------------------|--------|
| Energy scale (1.5%) | Energy scale (2.6%) | Covered |
| pi0 efficiency / cocktail | Not needed (ABCD replaces R_gamma) | N/A |
| Acceptance | **Not included** | **See discussion** |
| Luminosity (BBC trigger xsec) | Luminosity (MBD trigger xsec) | Covered |
| Conversion correction | **Not included** | **Gap** |

### CMS categories (from 1807.00782, 13 TeV BDT)

| Category | PPG12 equivalent | Status |
|----------|------------------|--------|
| Energy scale | Energy scale | Covered |
| Pileup | **Not included as systematic** | **Gap** |
| BDT signal template (MC generator) | BDT model variation | Partially covered |
| BDT background template | Non-tight selection variation | Covered |
| Isolation efficiency SF | Implicit in purity variations | Partially covered |

---

## 2. Missing or Underrepresented Systematics

### 2a. R_bg (ABCD Background Correlation Factor) -- HIGH PRIORITY

ATLAS devotes significant effort to validating R_bg != 1 by subdividing control regions B and D into B'/B'' and D'/D'' with shifted isolation boundaries, finding 15--40% deviations from unity that propagate as < 1% to 3.6% on the cross section.

PPG12 uses R_bg = 1 (the `R = params[7]` in `myfunc()` of `CalculatePhotonYield.C`). The wiki ABCD article states: "PPG12 currently uses R = 1 (uncorrelated assumption) in the default analysis, with deviations propagated as a systematic uncertainty." However, **no R_bg variation appears in the VARIANTS list** or in `SYST_TYPES`/`SYST_GROUPS`. The non-tight and non-iso variations indirectly probe ABCD decorrelation sensitivity, but they do not directly test R_bg != 1.

**Recommendation:** Add an explicit R_bg systematic, either by:
- Data-driven B'/D' subdivision (ATLAS method)
- Varying R_bg = 0.9 / 1.1 (or wider, based on MC estimate)
- Computing R_bg from background MC and propagating the deviation

This is standard practice at ATLAS and ALICE (the alpha_MC correction).

### 2b. Pileup / Double Interaction -- MEDIUM PRIORITY

CMS includes pileup as a systematic (0--11% in pp at 5.02 TeV). PHENIX at 510 GeV applied rate-dependent pileup corrections. PPG12 has extensive double-interaction studies (`DoubleInteractionCheck.C`, `run_showershape_double.sh`, `pileup_mixing.tex`, `double_interaction_efficiency.tex`) but the analysis note currently states:

> "The pileup mixing is not applied as a correction to the nominal cross-section result, since the nominal analysis uses the 1.5 mrad period where the pileup fraction is small (13.3% cluster-weighted) and the effect on the BDT-based photon identification is within the existing systematic uncertainties."

While this is defensible for 1.5 mrad (7.2% double interaction fraction, 13.3% cluster-weighted), **no explicit pileup systematic uncertainty is propagated to the final budget**. The double interaction studies demonstrate sub-percent overall effects after weighting, but a formal uncertainty should be assigned, especially if the 0 mrad data (18.7% fraction) is ever included.

**Recommendation:** Either:
- Add a one-sided pileup systematic based on the difference between single-only and blended (single+double) efficiency/purity
- Or explicitly document in `systematics.tex` why the effect is negligible and covered by existing variations (e.g., it is smaller than the efficiency variation)

### 2c. Material Budget / Photon Conversion -- LOW-MEDIUM PRIORITY

ATLAS propagates a material budget systematic that accounts for photon conversions before the calorimeter, which reduce the photon reconstruction efficiency. PHENIX included conversion corrections. For sPHENIX, the material budget before the EMCal (beam pipe, MVTX, INTT, TPC, inner HCAL support) creates a non-negligible conversion probability.

PPG12 does not include a material systematic. If the simulation correctly models the material, this should cancel in the efficiency ratio. However, mismodeling of the material budget is a standard source of uncertainty.

**Recommendation:** Assess whether the GEANT4 material description of the sPHENIX inner detectors has been validated. If material budget uncertainty is available from the detector group, propagate it as a systematic on the reconstruction efficiency.

### 2d. MC Generator Dependence -- LOW-MEDIUM PRIORITY

ATLAS cross-checks signal leakage fractions between Pythia and Sherpa generators. CMS validates BDT signal templates against data control samples (Z -> ee). PPG12 currently uses only Pythia 8 for signal MC.

The `bdt_model` variation (using different BDT model architectures) partially probes this, and the MC purity closure correction addresses ABCD assumption violations. However, a comparison of efficiency and purity using an alternative generator (e.g., Herwig, or Pythia with different tunes) would strengthen the analysis.

**Recommendation:** If alternative generator samples are available, run the pipeline with them as a cross-check. If not, this should be noted as a limitation.

### 2e. Unfolding Method Variation -- LOW PRIORITY

PPG12 evaluates the unfolding systematic by toggling response matrix reweighting (on/off), finding < 10%. ATLAS uses SVD unfolding and varies the regularization parameter. CMS evaluates bin migration effects.

The current single-variation approach is adequate for a first measurement, but could be strengthened by:
- Varying the number of Bayesian iterations
- Comparing Bayesian vs bin-by-bin unfolding
- Varying the prior distribution

This is low priority since the current variation already captures < 10%.

### 2f. Isolation Efficiency (Signal-Side) -- LOW PRIORITY

ATLAS and CMS separately evaluate the efficiency of the isolation cut for signal photons (distinct from the background purity effect). PPG12 folds this into the overall efficiency from MC, but does not vary the isolation cut on the signal side to assess the MC modeling of the isolation efficiency.

The `mc_iso_scale` and `mc_iso_shift` variations exist in `make_bdt_variations.py` but are currently flagged as `syst_type=None` (cross-checks only, not propagated to the final uncertainty). This may be intentional if they are deemed negligible.

**Recommendation:** Evaluate whether the MC isolation scale/shift variations produce significant deviations. If so, promote them to a formal systematic.

---

## 3. Dominant Systematics

### PPG12 hierarchy (from `systematics.tex`)

At lower ET (8--12 GeV):
- **Energy scale: 20--25%** (dominant)
- Purity (sideband isolation): 10--20%
- Purity (fitting): 10--20%
- Non-tight selection: up to 15%
- Efficiency: 5--10%
- Energy resolution: < 5%
- Unfolding: < 10%
- MBD trigger: < 10%

At higher ET (>20 GeV):
- Energy scale remains significant (~8%)
- Purity variations decrease
- Isolation and tight variations contribute ~10%

### Comparison with other experiments

| Regime | ATLAS | ALICE | PHENIX | PPG12 |
|--------|-------|-------|--------|-------|
| Low ET | Purity (5--15%) | Isolation probability (20%) | R_gamma amplification (30--50%) | **Energy scale (20--25%)** |
| Mid ET | Energy scale (5--10%) | Energy scale (5--10%) | Energy scale (5--8%) | Energy scale + purity (~15%) |
| High ET | Energy scale (5--16%) | Energy scale (5--10%) | Energy scale (5--8%) | Energy scale (8--10%) |

**Finding: PPG12's energy scale systematic is unusually large.** The 2.6% energy scale uncertainty propagates to 8--25% on the cross section, driven by the steep falling spectrum (power law index n ~ 5--6). This is the expected physics: a 2.6% shift in ET moves events across bins on a steeply falling spectrum, so the fractional change in the cross section is approximately n * 0.026 = 13--16%. The 8--25% range in the analysis note is consistent with this.

For comparison:
- ATLAS uses a 1--2% energy scale uncertainty, propagating to 5--16% on the cross section
- PHENIX achieved 1.5% energy scale precision
- ALICE quotes ~3% energy scale uncertainty

The 2.6% scale uncertainty recommended by the sPHENIX Calorimeter Calibration Working Group is the key driver. **Reducing this from 2.6% to ~1.5% would roughly halve the dominant systematic.** This is a detector calibration issue, not an analysis methodology issue.

**The overall pattern is consistent with expectations:** energy scale dominates (especially at lower pT where the spectrum is steepest), with purity contributions significant at the lowest pT bins. This matches the ATLAS and ALICE experience.

---

## 4. Luminosity Systematic

### Current status

`systematics.tex` states:
> "The luminosity systematic uncertainty is propagated directly from the minimum-bias dimuon (MBD) trigger cross-section measurement of 25.2 +2.3/-1.7 mb \cite{MBDxsec}."

This gives an asymmetric uncertainty of +9.1%/-6.7%. The code implements this in `make_bdt_variations.py`:
```python
LUMI_SYST = {"down": (25.2 - 23.5) / 25.2, "up": (27.5 - 25.2) / 25.2}
```

### The `\cite{MBDxsec}` issue

The citation key `MBDxsec` likely refers to a forthcoming sPHENIX MBD cross-section measurement. If this paper is not yet published or available as a preprint, the systematic value is not independently verifiable. This should be resolved before publication:
- If a preliminary result exists, cite the sPHENIX internal note number
- If based on PHENIX BBC measurements scaled to MBD, cite the PHENIX papers and state the extrapolation method
- The PHENIX BBC trigger cross section was 21.8 +/- 2.1 mb (pp at 200 GeV), roughly 10% uncertainty

### Comparison with other experiments

| Experiment | Luminosity uncertainty | Method |
|------------|----------------------|--------|
| ATLAS 13 TeV | 1.7--2.1% | van der Meer scans |
| CMS 13 TeV | 2.3--2.5% | Pixel cluster counting |
| ALICE 13 TeV | 1.6--2.1% | van der Meer scans |
| PHENIX pp 200 GeV | ~10% | BBC trigger cross-section |
| **PPG12** | **+9.1%/-6.7%** | MBD trigger cross-section |

PPG12's luminosity uncertainty is 3--5x larger than LHC experiments. This is expected: RHIC does not perform van der Meer scans with the same precision as the LHC, and the MBD trigger cross-section carries inherent model dependence from the Glauber/BBC framework. However, at +9.1%, this is one of the largest single contributions to the systematic budget.

**Recommendations:**
1. Resolve the dangling `\cite{MBDxsec}` before submission -- either cite a published/submitted paper or an internal note
2. Clearly state whether the 25.2 mb value is from a sPHENIX-specific measurement or adopted from PHENIX with corrections
3. In the analysis note, add a sentence noting that the luminosity uncertainty is common to all sPHENIX cross-section measurements and will improve with dedicated van der Meer scan analyses

---

## 5. Correlation Treatment

### Current PPG12 treatment

The analysis note (`systematics.tex`) states:
> "Variations within each systematic group are combined in quadrature, and the group totals are then combined in quadrature to produce the final systematic uncertainty band."

All systematics appear to be treated as **bin-by-bin uncertainties** -- each pT bin gets its own fractional uncertainty from each source, and these are summed in quadrature. The systematic uncertainty band is drawn as a box around each data point in the final cross-section plot.

The code (`calc_syst_bdt.py`) confirms this: deviations are computed per bin, then aggregated.

**What is missing:** The analysis note does not explicitly classify each systematic as:
- **Type A (uncorrelated):** Statistical-like, varies independently bin to bin
- **Type B (correlated):** Moves all bins in the same direction (e.g., energy scale, luminosity)
- **Type C (global normalization):** Affects all bins identically (luminosity)

### What other experiments do

**ATLAS:** Classifies systematics by correlation type. Energy scale and luminosity are fully correlated across pT bins. Purity and background systematics may be partially correlated. When providing results for combination or reinterpretation, ATLAS publishes the full covariance matrix or lists of nuisance parameters.

**PHENIX:** Typically quotes three categories of uncertainties separately:
1. Statistical (uncorrelated)
2. Point-to-point systematic (partially correlated)
3. Overall normalization (fully correlated, e.g., energy scale + luminosity)

The final figure often shows statistical error bars, point-to-point systematic boxes, and a separate normalization uncertainty band or percentage.

**ALICE:** Similar to PHENIX -- separates uncorrelated, partially correlated, and fully correlated sources. The luminosity and energy scale are drawn as a separate global box at the bottom of the figure.

### PPG12 classification (suggested)

| Source | Correlation across pT bins | Presentation |
|--------|---------------------------|-------------|
| Energy scale | Fully correlated (same shift applied everywhere) | Separate box or percentage |
| Energy resolution | Fully correlated | With energy scale |
| Luminosity | Global normalization | Separate percentage or box |
| MBD trigger efficiency | Fully correlated | With global normalization |
| Purity (all sub-sources) | Partially correlated | Point-to-point boxes |
| Efficiency (tight ID) | Partially correlated | Point-to-point boxes |
| Unfolding | Partially correlated | Point-to-point boxes |

**Recommendation:** In the analysis note and final data tables, separate the uncertainties into at least two classes:
1. **Point-to-point (uncorrelated + partially correlated):** purity, efficiency, unfolding, NPB
2. **Global (fully correlated):** energy scale, energy resolution, luminosity, MBD trigger

This is standard practice for cross-section measurements that may be used in global fits or comparisons with theory.

---

## 6. Summary of Findings

### Coverage assessment

PPG12 covers 6 of the ~8 standard systematic categories for an isolated photon cross-section measurement. This is a reasonable budget for a first measurement.

### Action items by priority

| Priority | Item | Estimated impact |
|----------|------|-----------------|
| HIGH | Add R_bg (ABCD correlation factor) systematic | Potentially 1--5% based on ATLAS experience |
| HIGH | Resolve `\cite{MBDxsec}` citation | Documentation only |
| HIGH | Classify systematics by correlation type | Presentation only, but important for data reuse |
| MEDIUM | Assess pileup systematic formally | Likely < 2% for 1.5 mrad, but should be documented |
| MEDIUM | Promote `mc_iso_scale`/`mc_iso_shift` to systematics if non-negligible | Evaluate from existing cross-check results |
| LOW | Material budget systematic | Requires detector group input |
| LOW | MC generator dependence | Requires alternative MC samples |
| LOW | Additional unfolding variations (iterations, method) | Current variation gives < 10% |

### What PPG12 does well

1. **Automated pipeline:** The 3-step workflow (generate configs, run condor, aggregate) is reproducible and scalable -- superior to hand-written systematic macros.
2. **Comprehensive purity decomposition:** Five sub-sources for purity (non-tight, non-iso, NPB, MC closure, fitting) is more granular than typical ATLAS/ALICE treatments.
3. **Parametric cut variations:** ET-dependent cuts with intercept/slope variations properly capture the kinematic dependence.
4. **Asymmetric uncertainties:** Two-sided variations are preserved (up/down separately), not forced to be symmetric.
5. **Luminosity treatment:** Asymmetric luminosity uncertainty correctly propagated from MBD cross-section.

### What could strengthen the analysis

1. Explicit R_bg systematic (standard at ATLAS/ALICE)
2. Correlation classification (Type A/B/C) in the analysis note
3. Formal pileup uncertainty, even if small
4. Reduction of the energy scale uncertainty from 2.6% (detector calibration improvement)
5. Cross-check with alternative MC generator for signal leakage fractions

---

## Source Files Reviewed

- `PPG12-analysis-note/systematics.tex` -- systematic uncertainty chapter
- `SYST_PIPELINE_README.md` -- 3-step pipeline documentation
- `efficiencytool/make_bdt_variations.py` -- VARIANTS, SYST_TYPES, SYST_GROUPS definitions
- `plotting/calc_syst_bdt.py` -- aggregation and quadrature logic
- `efficiencytool/CalculatePhotonYield.C` -- ABCD solver (myfunc, R parameter)
- `PPG12-analysis-note/double_interaction_efficiency.tex` -- pileup efficiency study
- `PPG12-analysis-note/pileup_mixing.tex` -- pileup mixing cross-check
- `wiki/physics/measurements/lhc-photon-measurements.md` -- ATLAS/CMS/ALICE methods
- `wiki/physics/measurements/rhic-photon-measurements.md` -- PHENIX methods
- `wiki/physics/techniques/abcd-sideband-method.md` -- ABCD systematic treatments
