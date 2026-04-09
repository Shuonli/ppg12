# Theory Comparison Review: PPG12 NLO Predictions

**Date**: 2026-04-08
**Scope**: Review of JETPHOX setup, Vogelsang calculation, scale variations, isolation matching, PDF choice, and missing comparisons in the PPG12 final cross-section plot.

---

## 1. JETPHOX Setup

### What PPG12 Uses

Verified directly from the JETPHOX parameter file (`/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/working/parameter.indat`):

| Parameter | Value | Notes |
|-----------|-------|-------|
| **PDF** | `CT14lo` (line 80) | LO PDF set with NLO matrix elements -- see Finding 1 below |
| **FF** | BFG set II (`0200`, line 156) | Standard choice; BFG I vs II negligible after isolation |
| **Scale choice** | `4` = cm * pT_photon (line 169) | Multiplied by cm factor |
| **Scale cm (nominal)** | 1.0D0 (lines 173/177/181 for mu_F/mu_R/mu_ff) | i.e., mu = ET |
| **Scale cm (low)** | 0.5D0 | mu = ET/2, set via sed in `run_condor05_SL7.sh` |
| **Scale cm (high)** | Presumably 2.0D0 | mu = 2*ET, set via `run_condor20_SL7.sh` |
| **alpha_s** | From LHAPDF (option `0`, line 196) | Determined by the PDF set |
| **N_f** | 5 (line 200) | Standard |
| **alpha_em** | Constant, 0 loops (line 188) | i.e., alpha_em = 1/137 |
| **Contributions** | `dir,onef` at NLO (`TRUE,TRUE`) | Direct + fragmentation, both at NLO |
| **Box diagram** | Off (`0`, line 208) | Born only for direct; box is numerically tiny |
| **sqrt(s)** | 200 GeV (line 286) | Correct |
| **Rapidity** | [-0.7, 0.7] (lines 290/294) | Matches PPG12 eta acceptance |
| **Isolation** | Fixed cone, R = 0.3, ET_max = 4 GeV (lines 315/327/339) | See Finding 3 below |

### Are These Well-Motivated?

The overall setup is standard and well-motivated for an NLO isolated photon calculation at RHIC kinematics. BFG set II, 5 active flavors, constant alpha_em, and the direct+fragmentation decomposition at NLO are all standard choices matching the literature (Aurenche et al. hep-ph/0602133, Ichou & d'Enterria 1005.4529).

---

## Finding 1: CT14lo PDF with NLO Matrix Elements (CRITICAL)

### The Issue

The analysis note (`results.tex` line 21 and `ppg12_conf_note/main.tex` line 230) both explicitly state:

> "the CT14lo parton distribution functions are employed"

This is confirmed by the JETPHOX parameter file (line 80): `CT14lo`.

**Using an LO PDF set with NLO partonic matrix elements is formally inconsistent.** The standard practice for an NLO calculation is to use NLO PDFs (CT14nlo). The reason: PDFs are extracted from data using matrix elements at the same perturbative order. LO PDFs are fitted to data using LO matrix elements and compensate for missing NLO corrections by distorting the gluon and sea distributions. Combining LO PDFs with NLO matrix elements double-counts some effects and misses others, leading to predictions that are neither formally LO nor NLO.

### What Other Experiments Use

- **ATLAS** (1908.02746, 1807.00782): CT14nlo, MMHT2014nlo, NNPDF3.0nlo for NLO comparisons; NNPDF3.1nnlo and MMHT2014nnlo for NNLO comparisons.
- **CMS** (1906.01371): CT14 NLO PDFs with JETPHOX.
- **ALICE** (1903.02209): CT14 NLO PDFs with JETPHOX.
- **PHENIX** (1205.5533): CTEQ6M (NLO) PDFs with the Vogelsang NLO code.

No major experiment uses LO PDFs with NLO JETPHOX in their published results.

### Impact

The quantitative impact depends on the x range. At the PPG12 kinematics (x ~ 0.04-0.35, Q ~ 8-35 GeV), differences between CT14lo and CT14nlo gluon PDFs are typically 5-15%, with CT14lo having a harder gluon at moderate x. This is comparable to the scale uncertainty band itself. The effect would shift the central NLO prediction and could artificially improve or worsen data-theory agreement.

### Recommendation

**Switch to CT14nlo** (or at minimum CT14nlo, member 0) for the JETPHOX calculation. This requires only changing line 80 of `parameter.indat` from `CT14lo` to `CT14nlo` and rerunning. The alpha_s from LHAPDF will automatically update to the NLO running value. This is a one-line fix in the JETPHOX configuration.

The wiki paper note (`wiki/raw/papers/C_theory/1506.07443.md`, line 95) already flags this concern:

> "JETPHOX is an NLO calculation, so using CT14 NLO PDFs would be formally consistent."

---

## Finding 2: Scale Variation (ADEQUATE)

### What PPG12 Does

PPG12 runs JETPHOX at three scales:

| Label | mu value | File |
|-------|----------|------|
| "05" | mu = 0.5 * ET | `NLO/rootFiles/jetPHOX_05.root` |
| "10" (nominal) | mu = ET | `NLO/rootFiles/jetPHOX_10.root` |
| "20" | mu = 2.0 * ET | `NLO/rootFiles/jetPHOX_20.root` |

In `plot_final_selection.C`, the scale band is constructed as an asymmetric error: upper edge = mu=0.5*ET (larger cross section from larger alpha_s at lower scale), lower edge = mu=2*ET. This is standard three-point scale variation with mu_R = mu_F = mu_ff varied simultaneously.

### Is This Standard?

Yes. The simultaneous three-point variation of all scales from mu/2 to 2*mu is the standard NLO uncertainty prescription used by ATLAS, CMS, ALICE, and PHENIX. Some collaborations also perform a 7-point variation (varying mu_R and mu_F independently by factors of 1/2 and 2, excluding opposite extremes), which typically gives a slightly larger band. PPG12's choice is the more conservative (smaller band) option.

### Central Scale Choice

The wiki notes that mu = pT/2 is often preferred as the central value (motivated by NLL resummation by Aurenche et al. and empirically providing better fits to fixed-target and ISR data). PPG12 uses mu = ET as the central value with pT/2 and 2*pT as the band edges. This is a defensible choice -- ATLAS and CMS also use mu = ET as the nominal. The choice of central scale affects the central prediction but not the uncertainty band, so this is a matter of convention rather than error.

### Assessment

The scale variation implementation is correct and standard. No action needed.

---

## Finding 3: Isolation Matching (IMPORTANT DISCREPANCY)

### JETPHOX Isolation Parameters

From `parameter.indat`:
- **Cone radius**: R = 0.3 (line 327) -- matches experiment
- **ET_max**: 4.0 GeV (line 339) -- **does NOT match the experimental parametric cut**

### Experimental Isolation

PPG12 uses a parametric (ET-dependent) isolation threshold:

```
reco_iso_max = b + s * ET = 0.502 + 0.0433 * ET
```

At representative ET values:
- ET = 10 GeV: E_T^max = 0.93 GeV
- ET = 20 GeV: E_T^max = 1.37 GeV
- ET = 30 GeV: E_T^max = 1.80 GeV

The truth-level isolation in MC uses a flat cut: `iso_ET_truth < 4 GeV`.

### The Mismatch

JETPHOX uses a flat 4 GeV isolation threshold matching the truth-level cut, not the tighter parametric reco-level cut. At ET = 10 GeV, the JETPHOX isolation is 4x looser than the experimental reco isolation (4.0 vs 0.93 GeV). At ET = 30 GeV, it is still 2x looser (4.0 vs 1.80 GeV).

This is partially justified: the measurement is corrected for efficiency back to the truth-level fiducial region defined by iso_ET_truth < 4 GeV, so the JETPHOX calculation should match the truth-level definition. The key question is whether the efficiency correction fully accounts for the difference between truth iso < 4 GeV and reco iso < 0.502 + 0.0433*ET. If the unfolding and efficiency correction are performed correctly with respect to the truth-level fiducial definition, then the JETPHOX calculation at 4 GeV is the appropriate comparison.

However, the `NLO/PlotTruthIso.C` study shows that the isolated cross section varies significantly with the isolation threshold, especially for the fragmentation component. Going from iso = 4 GeV to iso = 1 GeV changes the fragmentation fraction measurably. PPG12 should verify that the fiducial cross section definition used for the efficiency correction is genuinely iso_ET_truth < 4 GeV and that this matches the JETPHOX run.

### Cone vs Smooth

The wiki's isolation-theory article warns about ~5% differences between standard (fixed) cone and Frixione smooth-cone isolation. PPG12 correctly uses standard cone in JETPHOX, matching the experimental procedure. This is not a concern.

### R = 0.3 Caveat

The wiki correctly notes that R = 0.3 is borderline for NLO reliability. Catani et al. (hep-ph/0204023) showed that R = 0.4 is safely in the NLO regime while R = 0.1 gives unphysical results. At R = 0.3, unresummed log(R) corrections are moderate but near the edge of applicability. This is worth noting as a caveat but is not easily addressed without changing the experimental analysis.

### Recommendation

- Verify that the efficiency correction truly uses iso_ET_truth < 4 GeV as the fiducial denominator, not the parametric reco cut.
- Consider running JETPHOX also at iso = 2 GeV as a cross-check to assess the sensitivity of the prediction to the isolation threshold.
- The R = 0.3 cone-size caveat should be acknowledged in the analysis note.

---

## Finding 4: Vogelsang Calculation (MULTIPLE ISSUES)

### What Is Shown

The Vogelsang NLO calculation is shown on the final plot with scale bands (sc05, sc1, sc2). Data files are in `plotting/sphenix_nlo/photons_newphenix_sc{05,1,2}.dat`. The plotting code reads columns: `pt, yield, fragyield, foo`, sums `yield + fragyield`, and integrates over PPG12 pT bins.

### Issues

**4a. No isolation applied.** Both the analysis note and the plot legend explicitly state this:

> "The calculation by Gordon and Vogelsang does not require a E_T^iso condition" (conf note line 233)

The plot legend annotation reads: `(no E_T^iso / mu_f = mu_F = mu_R = E_T^gamma)`.

This means the Vogelsang curve includes the full fragmentation contribution without isolation suppression. At PPG12 kinematics, isolation suppresses fragmentation by 60-80%, so the non-isolated cross section is ~10-20% higher than the isolated one. Comparing an isolated measurement to a non-isolated theory prediction is not a proper apples-to-apples comparison.

**4b. Unknown PDF set.** The wiki states CTEQ6M PDFs, but this cannot be verified from the data files alone. The `nlo_plot.C` helper (by J. Nagle) labels the curves as "NLO pQCD W. Vogelsang" but does not document the PDF set. The data files carry no metadata.

**4c. Unknown isolation, cone size, and ET binning parameters.** Without documentation of the Vogelsang calculation setup, it is impossible to verify that the rapidity range, sqrt(s), and other parameters match the PPG12 measurement.

**4d. Lack of citation detail.** The conf note cites "Gordon:1993qc" (PRD 48, 3136, 1993) but this paper is >30 years old and predates modern PDF sets. It is unclear whether Werner Vogelsang provided updated calculations specifically for sPHENIX or whether these are legacy files.

### Recommendation

- Clarify the Vogelsang calculation parameters (PDF set, isolation, rapidity range, sqrt(s)) in the analysis note.
- Either request an isolated version of the Vogelsang calculation, or clearly state in the note that the comparison is to a non-isolated NLO prediction and discuss the expected ~10-20% excess.
- Consider whether showing a non-isolated NLO prediction alongside an isolated measurement adds value or confusion.

---

## Finding 5: Missing Comparisons

### What PPG12 Currently Shows

1. Data (sPHENIX)
2. PYTHIA 8.307 (Detroit tune)
3. JETPHOX NLO (CT14lo, BFG II, three scales)
4. Vogelsang NLO (no isolation)
5. PHENIX direct photon data (separate plot)

### What Other Experiments Include

| Comparison | ATLAS | CMS | ALICE | PHENIX | PPG12 |
|------------|-------|-----|-------|--------|-------|
| JETPHOX NLO | Yes | Yes | Yes | No | Yes |
| NNLOJET NNLO | Yes (13 TeV) | No | No | No | No |
| Vogelsang NLO | No | No | No | Yes | Yes |
| PYTHIA 8 | Yes | Yes | Yes | No | Yes |
| Sherpa | Yes | Yes | No | No | No |
| PDF uncertainty band | Yes | Sometimes | Yes | No | No |
| Multiple PDF sets | Yes (4-5 sets) | Sometimes | Yes | No | No |
| xT scaling | Sometimes | No | No | Sometimes | No |
| Direct/frag decomposition | Sometimes | No | No | No | No (but available in PlotTruthIso.C) |

### Potentially Valuable Additions

**5a. PDF uncertainty band (RECOMMENDED).** CT14 provides 56 Hessian error sets at 90% CL. Running JETPHOX with all error sets would give a PDF uncertainty band of ~3-5%, which is smaller than the scale band but represents a different (and independent) source of uncertainty. ATLAS routinely shows this. Since JETPHOX already loops over PDF error sets when the "call all PDF error sets" flag is TRUE (it IS set to TRUE in the parameter file), the infrastructure exists.

**5b. NNLO comparison (DESIRABLE BUT NOT AVAILABLE).** NNLOJET NNLO calculations for isolated photon production at sqrt(s) = 200 GeV do not currently exist (available only at 13 TeV). This would be the gold standard but requires external theory support. Worth flagging as a future extension.

**5c. xT scaling (OPTIONAL).** Expressing the cross section as a function of xT = 2pT/sqrt(s) and overlaying measurements from different sqrt(s) values (ISR, SPS, Tevatron, LHC) tests the universality of the pQCD subprocess picture. This is a standard analysis in PHENIX/STAR prompt photon papers but is more of a complementary plot than a replacement.

**5d. Multiple PDF sets (OPTIONAL).** Running JETPHOX with NNPDF3.1nlo and MMHT2014nlo in addition to CT14nlo would demonstrate PDF robustness. The differences are typically 5-10% and are covered by the scale band, so this adds confirmation rather than new physics.

**5e. Direct/fragmentation decomposition plot (OPTIONAL).** The `PlotTruthIso.C` macro already produces this. Showing the fraction of direct vs fragmentation as a function of ET for the PPG12 isolation parameters would be pedagogically valuable and is straightforward.

---

## Summary of Findings and Priority

| # | Finding | Severity | Action |
|---|---------|----------|--------|
| 1 | CT14lo used with NLO matrix elements | **Critical** | Switch to CT14nlo in JETPHOX parameter.indat |
| 2 | Scale variation | OK | Standard 3-point variation, no change needed |
| 3 | JETPHOX isolation = 4 GeV flat, vs parametric reco cut | **Important** | Verify fiducial definition matches; consider iso sensitivity study |
| 4 | Vogelsang calculation has no isolation | **Important** | Clarify in note; request isolated version or explain discrepancy |
| 5a | No PDF uncertainty band | Recommended | Enable existing PDF error set infrastructure |
| 5b | No NNLO comparison | Desirable | Not available at 200 GeV; flag as future work |
| 5c | No xT scaling | Optional | Complementary physics plot |
| 5d | No multiple PDF sets | Optional | Would demonstrate PDF robustness |
| 5e | No direct/frag decomposition | Optional | Straightforward from existing PlotTruthIso.C |

---

## Files Examined

| File | Purpose |
|------|---------|
| `/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/working/parameter.indat` | JETPHOX configuration (PDF, FF, scales, isolation, kinematics) |
| `/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/run_condor_SL7.sh` | Nominal scale condor submission |
| `/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/run_condor05_SL7.sh` | mu=0.5 scale variation |
| `NLO/MakeJetPHOXhisto.C` | Processes JETPHOX trees into histograms |
| `NLO/PlotTruthIso.C` | Truth isolation threshold study |
| `NLO/rootFiles/jetPHOX_{05,10,20}.root` | Binned JETPHOX predictions |
| `plotting/plot_final_selection.C` | Final cross-section plot |
| `plotting/sphenix_nlo/photons_newphenix_sc{05,1,2}.dat` | Vogelsang tabulated predictions |
| `plotting/sphenix_nlo/nlo_plot.C` | J. Nagle's original NLO plotting helper |
| `PPG12-analysis-note/results.tex` | Theory comparison discussion in analysis note |
| `ppg12_conf_note/main.tex` | Conference note theory discussion |
| `wiki/physics/theory/nlo-calculations.md` | Wiki NLO theory article |
| `wiki/physics/theory/prompt-photon-production.md` | Wiki production mechanisms article |
| `wiki/physics/theory/fragmentation-functions.md` | Wiki BFG FF article |
| `wiki/physics/theory/isolation-theory.md` | Wiki isolation theory article |
