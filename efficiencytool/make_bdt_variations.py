#!/usr/bin/env python3
"""Generate BDT variation config files from a nominal base config.

Usage:
    python make_bdt_variations.py <base_config> [--outdir DIR]

Each entry in VARIANTS produces one output file named
config_bdt_{name}.yaml in the output directory (defaults to the same
directory as the base config).

Supported override keys and their YAML paths:
    bdt_model_name         -> input.bdt_model_name
    reco_noniso_min_shift  -> analysis.reco_noniso_min_shift
    mc_purity_correction   -> analysis.mc_purity_correction
    tight_bdt_min          -> analysis.tight.bdt_min
    nt_bdt_max             -> analysis.non_tight.bdt_max
    nt_bdt_min             -> analysis.non_tight.bdt_min
    npb_score_cut          -> analysis.common.npb_score_cut
    bdt_et_bin_edges       -> input.bdt_et_bin_edges
    bdt_et_bin_models      -> input.bdt_et_bin_models
"""

import argparse
import copy
import os
import sys

from ruamel.yaml import YAML

# ---------------------------------------------------------------------------
# Variant definitions
# Each dict must have a 'name' key.  All other keys are override parameters.
# ---------------------------------------------------------------------------
VARIANTS = [
    # Canonical all-range nominal. Auto-expanded by generate_variants() into
    # 3 on-disk configs: config_bdt_nom.yaml (all-range, identical to base),
    # config_bdt_nom_0rad.yaml (0mrad merge-feeder), config_bdt_nom_1p5mrad.yaml
    # (1.5mrad merge-feeder). Merge-feeders bake per-event lumi_weight =
    # lumi/lumi_target into MC fills so a plain hadd across periods
    # reproduces the all-range expectation.
    dict(name="nom",
         syst_type=None, syst_role=None),

    # NOTE (2026-04-22): the all-z configuration (beam-delivered lumi
    # 64.3718 pb^-1, vertex_cut_truth=9999) IS now the nominal. The
    # explicit `allz` / `allz_0rad` / `allz_1p5mrad` entries are removed
    # because they would duplicate nom. If a 60cm-fiducial cross-check is
    # needed, add a `cm60` variant with lumi=48.9309, vertex_cut_truth=60,
    # and period lumis 32.6574 / 16.2735.

    # No-vertex-cut cross-check: removes BOTH the reco (vertex_cut) AND truth
    # (vertex_cut_truth) fiducial cuts. Data uses the full beam-delivered
    # sample (|z_reco| unrestricted), and the MBD-eff denominator includes
    # all truth events. Pair with beam-delivered lumi (64.3718 pb^-1).
    # Tests the fiducial chain at its limit: if σ_novtxcut ≈ σ_nom, the
    # combined reco+truth fiducial selection is self-consistent.
    #iso scale and shift
    dict(name="mciso_noscaleshift", mc_iso_scale=1.0, mc_iso_shift=0.0,
         syst_type=None, syst_role=None),
    # Promoted 2026-04-25: probes additive-pedestal mismodeling between data
    # and MC iso ET (multiplicative scale stays at nominal 1.2). _no_scale
    # mixes pedestal+width with escale/eres, so we use _no_shift only.
    dict(name="mciso_no_shift",  mc_iso_shift=0.0,
         syst_type="iso_resolution", syst_role="one_sided"),
    # Generator-driven iso-eff systematic (HERWIG7 vs PYTHIA8 Detroit, 1.5 mrad).
    # The variant cross-section file Photon_final_bdt_iso_generator.root is NOT
    # produced by oneforall.sh; it is built post-hoc by
    # plotting/build_iso_generator_variant.C, which scales the nominal sigma by
    # eps_iso_P / eps_iso_H bin-by-bin. Hence no config keys are emitted here.
    dict(name="iso_generator",
         syst_type="iso_generator", syst_role="one_sided",
         aggregate_only=True),
    dict(name="mciso_no_scale",  mc_iso_scale=1.0,
         syst_type=None, syst_role=None),
    # Flat-threshold BDT partitions (cross-checks).  Tight is flat [T, 1.0];
    # non-tight is flat [nt_min, T] with no gap at T.  Explores the ABCD
    # sensitivity to both the tight cut T and the non-tight width.
    dict(name="flat_t90_nt50",
         tight_bdt_min_intercept=0.90, tight_bdt_min_slope=0.0,
         nt_bdt_max_intercept=0.90,    nt_bdt_max_slope=0.0,
         nt_bdt_min=0.50,
         syst_type=None, syst_role=None),
    dict(name="flat_t90_nt70",
         tight_bdt_min_intercept=0.90, tight_bdt_min_slope=0.0,
         nt_bdt_max_intercept=0.90,    nt_bdt_max_slope=0.0,
         nt_bdt_min=0.70,
         syst_type=None, syst_role=None),
    dict(name="flat_t85_nt50",
         tight_bdt_min_intercept=0.85, tight_bdt_min_slope=0.0,
         nt_bdt_max_intercept=0.85,    nt_bdt_max_slope=0.0,
         nt_bdt_min=0.50,
         syst_type=None, syst_role=None),
    dict(name="flat_t95_nt50",
         tight_bdt_min_intercept=0.95, tight_bdt_min_slope=0.0,
         nt_bdt_max_intercept=0.95,    nt_bdt_max_slope=0.0,
         nt_bdt_min=0.50,
         syst_type=None, syst_role=None),

    # Non-iso sideband boundary shift  →  syst: noniso (purity)
    dict(name="noniso04",     reco_noniso_min_shift=0.1,
         syst_type="noniso", syst_role="down"),   # tighter noniso window
    dict(name="noniso10",     reco_noniso_min_shift=1.0,
         syst_type="noniso", syst_role="up"),     # looser noniso window

    # Reco-iso parametric shift cross-check (2026-04-27): tighter intercept
    # 0.490 (vs nominal 0.502095) and slightly looser slope 0.037 (vs
    # 0.0433036). Probes joint sensitivity of purity and iso efficiency
    # to the parametric isolation cut shape; not a syst group entry.
    dict(name="iso_p49_s37",
         reco_iso_max_b=0.490, reco_iso_max_s=0.037,
         syst_type=None, syst_role=None),

    # NPB score cut  →  syst: npb_cut (two_sided). Looser/tighter cut
    # produces the expected anti-correlated yield shifts after the merge
    # fix (Pearson r=-0.998, post-fix mean +0.3%/-1.2%).
    dict(name="npb03",        npb_score_cut=0.3,
         syst_type="npb_cut", syst_role="down"),  # looser npb cut
    dict(name="npb07",        npb_score_cut=0.7,
         syst_type="npb_cut", syst_role="up"),    # tighter npb cut

    # Purity fit option  →  syst: purity_fit (purity)
    dict(name="purity_pade", fit_option=0,
         syst_type="purity_fit", syst_role="one_sided"),
    # Purity fit 1σ CI band  →  syst: purity_fit_ci (purity, two-sided).
    # Replaces the central f_purity_fit value with the upper/lower edge of
    # the 0.683 CL fitter confidence interval (TVirtualFitter::GetConfidenceIntervals)
    # in CalculatePhotonYield.C, controlled by analysis.fittingerror = +1/-1.
    dict(name="purity_fit_ci_up",   fittingerror=1,
         syst_type="purity_fit_ci", syst_role="up"),
    dict(name="purity_fit_ci_down", fittingerror=-1,
         syst_type="purity_fit_ci", syst_role="down"),
    # Use MC-driven purity correction ratio for data purity fit
    dict(name="mc_purity_correction", mc_purity_correction=1,
         syst_type="mc_purity_correction", syst_role="one_sided"),

    # Vertex reweighting off — kept as cross-check only. Demoted 2026-04-25:
    # under the all-z nominal, the legacy reco-vertex reweight is already off
    # (truth_vertex_reweight_on=1 is nominal), so this variant probes the
    # legacy path and overlaps the truth-vertex closure that the locked
    # budget intentionally drops.
    dict(name="vtxreweight0", vertex_reweight_on=0,
         syst_type=None, syst_role=None),

    # no unfolding reweighting — note: name says "unfolding", but the actual
    # effect is upstream in RecoEffCalculator (response prior is built with
    # f_reweight tilt; CalculatePhotonYield force-zeros reweight). Verified
    # non-null via Phase-1.A1 audit: variant builds an unweighted response
    # matrix and produces a real shift vs nominal.
    dict(name="no_unfolding_reweighting", reweight=0,
         syst_type="reweight", syst_role="one_sided"),

    # Photon-ID systematics — BDT cut placement (2026-04-25 locked).
    # tightup_p05: tight_bdt_min_intercept += 0.05 (slope unchanged).
    #   Nominal intercept = 0.8333..., new = 0.8833... Tightens the signal
    #   selection at all ET; effect is dominantly an efficiency loss in A.
    # ntdown_m10: nt_bdt_min_intercept -= 0.10 (slope unchanged).
    #   Nominal intercept = 0.7333..., new = 0.6333... Loosens the non-tight
    #   inner edge so C/D widen; tests bkg-template-shape stability.
    dict(name="tightup_p05",
         tight_bdt_min_intercept=0.8833333333333334,
         syst_type="photon_id_tight", syst_role="one_sided"),
    dict(name="ntdown_m10",
         nt_bdt_min_intercept=0.6333333333333333,
         syst_type="photon_id_nontight", syst_role="one_sided"),

    # BDT model cross-checks — single model across all ET bins. Kept as
    # cross-checks only (not part of the published systematic budget). The
    # nominal photon-ID BDT (base_v3E) uses 11 input features; base_v0E uses
    # only 7 of those, so the comparison entangles training-vintage effects
    # with feature-engineering changes and is not a clean closure variation.
    dict(name="bdtmodel_v0",  bdt_et_bin_edges=[8, 15, 35], bdt_et_bin_models=["base_v0",  "base_v0"],
         syst_type=None, syst_role=None),
    dict(name="bdtmodel_v0E", bdt_et_bin_edges=[8, 15, 35], bdt_et_bin_models=["base_v0E", "base_v0E"],
         syst_type=None, syst_role=None),

    # ET-binned BDT model cross-checks (etbin_v3E_v3E was retired 2026-04-24
    # when v3E became the nominal model, so it is now degenerate with nom).
    dict(name="etbin_E_E",      bdt_et_bin_edges=[8, 15, 35], bdt_et_bin_models=["base_E",   "base_E"],
         syst_type=None, syst_role=None),

    # b2bjet requirement — cross-check variants (syst_type=None). Require a
    # back-to-back jet above pt_min on top of nominal selection; enriches the
    # photon+jet component of the signal. Two thresholds kept: 5 GeV and 7 GeV.
    dict(name="b2bjet_pt5", common_b2bjet_cut=1, common_b2bjet_pt_min=5.0,
         syst_type=None, syst_role=None),
    dict(name="b2bjet_pt7", common_b2bjet_cut=1, common_b2bjet_pt_min=7.0,
         syst_type=None, syst_role=None),

    # b2bjet requirement with NPB common cut disabled — tests whether the
    # b2bjet photon+jet requirement alone controls purity without NPB.
    dict(name="b2bjet_pt5_npb0", common_b2bjet_cut=1, common_b2bjet_pt_min=5.0, npb_score_cut=0.0,
         syst_type=None, syst_role=None),
    dict(name="b2bjet_pt7_npb0", common_b2bjet_cut=1, common_b2bjet_pt_min=7.0, npb_score_cut=0.0,
         syst_type=None, syst_role=None),


    # Energy scale and resolution  →  syst: escale, eres
    # eres scheme (May 2026): additive ET-dependent smearing keyed off the
    # truth-photon pT, applied ONLY to the response-matrix fills in
    # RecoEffCalculator_TTreeReader.C. The sigma_extra(pT) curve is the
    # quadrature difference between an external data-resolution fit
    # (function_compare_wide.root) and the PPG12 PYTHIA sigma_MC fit
    # (q0=0.185, q2=0.04 from fit_sigma_mc.py). Nominal data fit is
    # f_energy_example_cE_0p05 → sigma_extra ~ 2% at high ET. Two-sided
    # eres systematic: down = cE_0p08 (~6%), up = no extra smearing (data
    # params all zero). The legacy multiplicative `cluster_eres` field is
    # retired but kept in field-map for backward compatibility.
    # Energy-scale envelope: placeholder at +/-1.1% pending the final
    # EMCal calibration. The +/-1.5% and +/-2.6% pairs are kept alongside
    # as cross-checks but the active journal-text systematic is taken from
    # the 1.1% pair (mapped to syst_type="escale"); the 1.5% and 2.6%
    # pairs have syst_type=None so they do not enter the syst aggregator.
    dict(name="energyscale11up",   clusterescale=1.011,
         syst_type="escale", syst_role="down"),
    dict(name="energyscale11down", clusterescale=0.989,
         syst_type="escale", syst_role="up"),
    dict(name="energyscale15up",   clusterescale=1.015,
         syst_type=None,     syst_role=None),
    dict(name="energyscale15down", clusterescale=0.985,
         syst_type=None,     syst_role=None),
    dict(name="energyscale26up",   clusterescale=1.026,
         syst_type=None,     syst_role=None),
    dict(name="energyscale26down", clusterescale=0.974,
         syst_type=None,     syst_role=None),
    # Role labels follow the *output* (cross-section) direction, matching the
    # escale convention above: less smearing → unfolded yield goes up,
    # more smearing → unfolded yield goes down. Without this convention the
    # breakdown plot would draw the eres band with sign flipped relative to
    # the per-type signed plot.
    # New additive-smearing eres variants. Nominal lives in the base config
    # (cluster_eres_data_p* = cE_0p05 params, cluster_eres_mc_p* = PYTHIA fit).
    # "up" variant (no extra smearing): zero out the data params; lambda
    # returns 0 when sigma_data^2 < sigma_MC^2.
    # "down" variant (cE_0p08 data fit): wider data → ~6% extra smearing.
    dict(name="eres_smear_none",
         cluster_eres_data_p0=0.0,  cluster_eres_data_p1=0.0,  cluster_eres_data_p2=0.0,
         syst_type="eres", syst_role="up"),
    dict(name="eres_smear_cE0p08",
         cluster_eres_data_p0=0.13, cluster_eres_data_p1=0.08, cluster_eres_data_p2=0.08,
         syst_type="eres", syst_role="down"),

    # ----------------------------------------------------------
    # T1: Unfolding regularization scan (Bayes iteration count).
    # The base CalculatePhotonYield writes h_unfold_sub_leak_{1..10}; these
    # variants are produced via postprocess_unfold_iter_scan.py which
    # promotes the iter-N histogram into h_unfold_sub_result so the
    # standard aggregator can pick them up. No new condor needed.
    # All-range pinned (run_min/max from nominal) so they generate as a
    # single bare-name config rather than auto-expanding into period
    # feeders — the iteration count is period-independent.
    # ----------------------------------------------------------
    # iter=1 is unconverged (Bayes regularization too strong), so it
    # inflates the iter-scan envelope. Demoted to cross-check only —
    # variant is still produced for the closure-overlay plot in the
    # unfolding appendix but is excluded from the systematic envelope.
    dict(name="unfold_iter1", resultit=1,
         run_min=47289, run_max=54000, lumi=64.3718, lumi_target=64.3718,
         vertex_cut_truth=9999.0, truth_vertex_reweight_on=1,
         syst_type=None, syst_role=None),
    dict(name="unfold_iter3", resultit=3,
         run_min=47289, run_max=54000, lumi=64.3718, lumi_target=64.3718,
         vertex_cut_truth=9999.0, truth_vertex_reweight_on=1,
         syst_type="unfold_iter", syst_role="max"),
    dict(name="unfold_iter4", resultit=4,
         run_min=47289, run_max=54000, lumi=64.3718, lumi_target=64.3718,
         vertex_cut_truth=9999.0, truth_vertex_reweight_on=1,
         syst_type="unfold_iter", syst_role="max"),

    # ----------------------------------------------------------
    # T2: Unfolding prior — flat-truth prior. Demoted to cross-check on
    # 2026-04-25: a fully-flat prior on the steeply-falling spectrum at
    # iter=2 produces an unphysical syst (>25000% at 32-36 GeV) because
    # the prior is too far from the true shape and 2 Bayes iterations
    # cannot recover. The iter scan (T1) covers regularization at 1-7%
    # physically; this flat-prior variant is kept only as a cross-check
    # to bound prior dependence in an extreme limit.
    # ----------------------------------------------------------
    dict(name="unfold_prior_flat", flat_prior=1,
         run_min=47289, run_max=54000, lumi=64.3718, lumi_target=64.3718,
         vertex_cut_truth=9999.0, truth_vertex_reweight_on=1,
         syst_type=None, syst_role=None),

    # ----------------------------------------------------------
    # T3: Double-interaction blending fraction ±20% per period.
    # DOUBLE_FRAC nominal (cluster-weighted, triple+ folded into double):
    #   0 mrad:   f_best = 0.290 (analytic 0.224, |Δf| = 0.066)
    #   1.5 mrad: f_best = 0.000 (analytic 0.079, |Δf| = 0.079)
    # Manual triplet (bare-name + 2 feeders) so each crossing carries its own
    # f_best while the all-range merge_periods.sh step produces a single
    # Photon_final_bdt_di_frac_fit.root for calc_syst_bdt.py. Only the
    # bare-name is tagged with syst_type so the aggregator picks up exactly
    # one all-range varied spectrum (not the two per-period pre-merge ones).
    # ----------------------------------------------------------
    dict(name="di_frac_fit",
         run_min=47289, run_max=54000, lumi=64.3718, lumi_target=64.3718,
         vertex_cut_truth=9999.0, truth_vertex_reweight_on=0,
         truth_vertex_reweight_file="/sphenix/user/shuhangli/ppg12/efficiencytool/truth_vertex_reweight/output/1p5mrad/reweight.root",
         syst_type="di_fraction", syst_role="one_sided"),
    dict(name="di_frac_fit_0rad", double_frac_override=0.290,
         run_min=47289, run_max=51274, lumi=47.2076, lumi_target=64.3718,
         vertex_cut_truth=9999.0, truth_vertex_reweight_on=1,
         truth_vertex_reweight_file="/sphenix/user/shuhangli/ppg12/efficiencytool/truth_vertex_reweight/output/0mrad/reweight.root",
         syst_type=None, syst_role=None),
    dict(name="di_frac_fit_1p5mrad", double_frac_override=0.000,
         run_min=51274, run_max=54000, lumi=17.1642, lumi_target=64.3718,
         vertex_cut_truth=9999.0, truth_vertex_reweight_on=1,
         truth_vertex_reweight_file="/sphenix/user/shuhangli/ppg12/efficiencytool/truth_vertex_reweight/output/1p5mrad/reweight.root",
         syst_type=None, syst_role=None),
]

# ---------------------------------------------------------------------------
# purity_nonclosure_ntbdt: ET-parametric tight/non-tight pair scan (cross-check)
# Each variant is defined by three anchors at ET=10 and ET=25 GeV:
#   tight_lo, nt_hi (non-tight upper), nt_lo (non-tight lower).
# Unlike the previous set, nt_hi is decoupled from tight_lo — the C' low
# anchor has nt_hi < tight_lo at ET=10, introducing a "gap" [nt_hi, tight_lo]
# that is excluded from both tight and non-tight regions. Not a systematic.
# ---------------------------------------------------------------------------

def _anchors_to_param(v10: float, v25: float):
    slope = (v25 - v10) / 15.0
    intercept = v10 - 10.0 * slope
    return intercept, slope

_NTBDT_LOW_ANCHORS = [
    # label, (tight_lo@10, nt_hi@10, nt_lo@10)
    ("A", (0.80, 0.80, 0.70)),
    ("B", (0.80, 0.80, 0.60)),
    ("C", (0.80, 0.70, 0.60)),  # GAP: nt_hi < tight_lo at ET=10
]
_NTBDT_HIGH_ANCHORS = [
    # label, (tight_lo@25, nt_hi@25, nt_lo@25)
    ("X", (0.70,  0.70,  0.40)),
    ("Y", (0.75,  0.75,  0.50)),
    ("Z", (0.75,  0.75,  0.40)),
]
for _low_label, (_t10, _h10, _l10) in _NTBDT_LOW_ANCHORS:
    for _high_label, (_t25, _h25, _l25) in _NTBDT_HIGH_ANCHORS:
        # C+Z = ntbdtpair_t80_75_h70_75_l60_40 is the current nominal cut
        # (promoted 2026-04-24). Kept here as a variant (degenerate with
        # nominal) so the config survives for easy swap if the nominal choice
        # changes later.
        _t_int, _t_slp = _anchors_to_param(_t10, _t25)
        _h_int, _h_slp = _anchors_to_param(_h10, _h25)
        _l_int, _l_slp = _anchors_to_param(_l10, _l25)
        _name = (f"ntbdtpair_t{int(round(_t10*100))}_{int(round(_t25*100))}"
                 f"_h{int(round(_h10*100))}_{int(round(_h25*100))}"
                 f"_l{int(round(_l10*100))}_{int(round(_l25*100))}")
        # 2026-04-25: All ntbdtpair anchor combos demoted to cross-check
        # (syst_type=None). The locked photon-ID coverage is now
        # bdt_tightup_p05 + bdt_ntdown_m10 (constant-shift variants on the
        # nominal C+Z parametric form), not the earlier A+X anchor pair.
        # Anchor scan kept for purity-non-closure cross-checks only.
        _syst_type, _syst_role = None, None
        VARIANTS.append(dict(
            name=_name,
            tight_bdt_min_intercept=_t_int, tight_bdt_min_slope=_t_slp,
            nt_bdt_max_intercept=_h_int,   nt_bdt_max_slope=_h_slp,
            nt_bdt_min_intercept=_l_int,   nt_bdt_min_slope=_l_slp,
            syst_type=_syst_type, syst_role=_syst_role,
        ))

# ---------------------------------------------------------------------------
# Tower-acceptance mask variants — data-driven phi-symmetry method.
#
# For each selection level (preselect, common, tight), an independent mask
# is built from the data-only phi-symmetry test: within each 2-eta-row band,
# fit Gaussian to the 256 per-phi data counts and flag phi bins with
# z < -2 (anomalously low relative to their band's phi-symmetric
# expectation). No MC assumption needed.
#
# The earlier log(R)/Poisson-based masks (mask_{common,tight,or}_*) have been
# retired; their configs are archived under archive_logR_masks/configs/ and
# their old results under archive_logR_masks/results/.
#
# Not counted in SYST_GROUPS quadrature; the cross-section shift vs nominal
# is the direct acceptance-systematic envelope.
# ---------------------------------------------------------------------------
_TOWER_MASK_FILE = "/sphenix/user/shuhangli/ppg12/efficiencytool/tower_masks_bdt_nom.root"
for _lvl in ("preselect", "common", "tight", "or"):
    # All four phi-symmetry mask variants (preselect, common, tight, or) are
    # kept as cross-checks only — none enter the systematic quadrature.
    # Phi-symmetry residuals are absorbed into the photon-ID response chain.
    VARIANTS.append(dict(
        name=f"mask_phisymm_{_lvl}",
        tower_mask_on=1,
        tower_mask_file=_TOWER_MASK_FILE,
        tower_mask_name=f"mask_phisymm_{_lvl}",
        syst_type=None, syst_role=None,
    ))

# ---------------------------------------------------------------------------
# Systematic type definitions
# mode: "two_sided" | "one_sided" | "max" | "placeholder"
# group: which SYST_GROUPS key this type contributes to
# ---------------------------------------------------------------------------
SYST_TYPES = {
    # ---- response/unfolding ----
    "reweight":     {"mode": "one_sided",   "group": "unfolding"},
    # ---- purity: tight/non-tight BDT placement, isolation sideband, fit form,
    #     MC closure correction (matches systematics.tex itemize) ----
    "photon_id_tight":    {"mode": "one_sided", "group": "purity"},
    "photon_id_nontight": {"mode": "one_sided", "group": "purity"},
    "noniso":       {"mode": "two_sided",   "group": "purity"},
    "purity_fit":   {"mode": "one_sided",   "group": "purity"},
    "purity_fit_ci": {"mode": "two_sided",  "group": "purity"},
    "mc_purity_correction": {"mode": "one_sided", "group": "purity"},
    # ---- NPB cut: standalone group (not part of purity) per analysis-note
    #     structure. Pearson r(npb03, npb07) = -0.998 (clean two-sided).
    "npb_cut":      {"mode": "two_sided",   "group": "npb"},
    # ---- detector iso resolution / pedestal: efficiency group ----
    "iso_resolution": {"mode": "one_sided", "group": "efficiency"},
    # ---- generator-driven iso-eff systematic: HERWIG7 vs PYTHIA8 Detroit
    #     1.5 mrad cross-check propagated as sigma_nom × eps_iso_P / eps_iso_H
    #     (built by plotting/build_iso_generator_variant.C). Documented in
    #     systematics.tex §5.4.2 and Appendix sec:appendix:herwig_xcheck.
    "iso_generator": {"mode": "one_sided", "group": "efficiency"},
    # ---- tower acceptance (phi-symmetry mask envelope) ----
    # Demoted to a cross-check (no quadrature contribution) on 2026-04-27;
    # phi-symmetry residuals are absorbed into the photon-ID response chain.
    # The mask_phisymm_or variant now carries syst_type=None.
    # ---- energy scale and resolution ----
    "escale":       {"mode": "two_sided",   "group": "escale"},
    "eres":         {"mode": "two_sided",   "group": "eres"},
    # ---- DI blending fraction: data-driven envelope from chi^2 fit at the
    #     preselection cut. Two period-pinned variants at f_best (0.290 at
    #     0 mrad, 0.000 at 1.5 mrad), each syst_role="one_sided"; calc_syst_bdt
    #     applies |Δσ| symmetrically per crossing and quadrature-sums.
    "di_fraction":  {"mode": "one_sided",   "group": "di_fraction"},
    "unfold_iter":  {"mode": "max",         "group": "unfolding"},
    # bdt_model retired 2026-04-29: bdtmodel_v0E swaps an alternative
    # photon-ID classifier with a different feature set (7 vs 11), so the
    # variation entangles training vintage with feature engineering and
    # was not a clean closure systematic. Kept as a cross-check only.
    # unfold_prior removed 2026-04-25 — flat-prior gave unphysical 25000%
    # syst at high pT (see unfold_prior_flat variant comment).
    # ---- retired keys kept declared so old VARIANTS dicts that still
    #     reference them load cleanly (variants below this line are all
    #     syst_type=None and contribute nothing to quadrature):
    #     tight_bdt, nt_bdt, vtx_reweight, b2bjet, timing,
    #     mbd, nor, bdt_model — removed 2026-04-25 / 2026-04-29.
}

# Quadrature grouping: group name -> list of syst_type names.
# Structure mirrors the analysis-note systematics itemize (purity merges
# tight/non-tight ID + noniso + fit form + MC closure; NPB and efficiency
# are stand-alone groups; unfolding combines reweight + iter scan).
SYST_GROUPS = {
    "escale":      ["escale"],
    "eres":        ["eres"],
    "purity":      ["photon_id_tight", "photon_id_nontight",
                    "noniso", "purity_fit", "purity_fit_ci",
                    "mc_purity_correction"],
    "efficiency":  ["iso_resolution", "iso_generator"],
    "unfolding":   ["reweight", "unfold_iter"],
    "di_fraction": ["di_fraction"],
    "npb":         ["npb_cut"],
}

# Groups included in the final total systematic quadrature sum (per-bin Resp).
# The lumi multiplicative-flat source is added via FLAT_SYSTS in
# calc_syst_bdt.py rather than here — it is not response-matrix-driven and is
# applied post-unfold.
FINAL_SYSTS = [
    "escale", "eres", "purity", "efficiency",
    "unfolding", "di_fraction", "npb",
]

# Multiplicative-flat sources (asymmetric fractional). Applied post-unfold by
# add_flat() in calc_syst_bdt.py instead of being summed inside any group.
#   lumi: from MBD inelastic xsec 25.2 +2.3/-1.7 mb (Joey/lumi memo, 2026).
#         Reconciled 2026-04-25 — note had stale +15.1/-11.2 from sigma=15.2.
# l1_plateau retired 2026-04-28: the bit-30 turn-on is corrected per-\etg
# bin in CalculatePhotonYield via the fitted L1 efficiency, so there is no
# residual systematic to assign.
LUMI_SYST   = {"down": (25.2 - 23.5) / 25.2, "up": (27.5 - 25.2) / 25.2}  # +9.13/-6.75%
FLAT_SYSTS  = {"lumi": LUMI_SYST}

# ---------------------------------------------------------------------------
# Mapping from flat override key -> (yaml_section_path, yaml_leaf_key)
# ---------------------------------------------------------------------------
OVERRIDE_MAP = {
    "bdt_model_name":        (["input"],                              "bdt_model_name"),
    "reco_noniso_min_shift": (["analysis"],                           "reco_noniso_min_shift"),
    "iso_emcalinnerr":       (["analysis"],                           "iso_emcalinnerr"),
    "tight_bdt_min":         (["analysis", "tight"],                  "bdt_min"),
    "nt_bdt_max":            (["analysis", "non_tight"],              "bdt_max"),
    "nt_bdt_min":            (["analysis", "non_tight"],              "bdt_min"),
    "tight_bdt_min_intercept": (["analysis", "tight"], "bdt_min_intercept"),
    "tight_bdt_min_slope": (["analysis", "tight"], "bdt_min_slope"),
    "nt_bdt_max_intercept": (["analysis", "non_tight"], "bdt_max_intercept"),
    "nt_bdt_max_slope": (["analysis", "non_tight"], "bdt_max_slope"),
    "nt_bdt_min_intercept": (["analysis", "non_tight"], "bdt_min_intercept"),
    "nt_bdt_min_slope": (["analysis", "non_tight"], "bdt_min_slope"),
    "npb_score_cut":         (["analysis", "common"],                 "npb_score_cut"),
    "common_wr_cogx_bound": (["analysis", "common"],                 "wr_cogx_bound"),
    "mc_iso_scale":          (["analysis"],                           "mc_iso_scale"),
    "reco_iso_max_b":        (["analysis"],                           "reco_iso_max_b"),
    "reco_iso_max_s":        (["analysis"],                           "reco_iso_max_s"),
    "apply_trigger_eff_correction": (["analysis"],                    "apply_trigger_eff_correction"),
    "trigger_eff_p0":        (["analysis"],                           "trigger_eff_p0"),
    "trigger_eff_mu":        (["analysis"],                           "trigger_eff_mu"),
    "trigger_eff_beta":      (["analysis"],                           "trigger_eff_beta"),
    "vertex_cut":            (["analysis"],                           "vertex_cut"),
    "vertex_cut_truth":      (["analysis"],                           "vertex_cut_truth"),
    "lumi_target":           (["analysis"],                           "lumi_target"),
    "vertex_reweight_on":    (["analysis"],                           "vertex_reweight_on"),
    "truth_vertex_reweight_on":   (["analysis"],                      "truth_vertex_reweight_on"),
    "truth_vertex_reweight_file": (["analysis"],                      "truth_vertex_reweight_file"),
    "cluster_node_name":     (["input"],                              "cluster_node_name"),
    "data_file":             (["input"],                              "data_file"),
    "photon_jet_file_branch_dir": (["input"],                       "photon_jet_file_branch_dir"),
    "common_b2bjet_cut":        (["analysis"],                 "common_b2bjet_cut"),
    "common_b2bjet_pt_min":     (["analysis"],                 "common_b2bjet_pt_min"),
    "cluster_mbd_time_min":        (["analysis"],                 "cluster_mbd_time_min"),
    "cluster_mbd_time_max":        (["analysis"],                 "cluster_mbd_time_max"),
    "use_topo_iso":        (["analysis"],                 "use_topo_iso"),
    "mc_iso_shift":        (["analysis"],                 "mc_iso_shift"),
    "fit_option":        (["analysis"],                 "fit_option"),
    "fittingerror":      (["analysis"],                 "fittingerror"),
    "mc_purity_correction": (["analysis"],              "mc_purity_correction"),
    "bdt_et_bin_edges":  (["input"],                              "bdt_et_bin_edges"),
    "bdt_et_bin_models": (["input"],                              "bdt_et_bin_models"),
    "run_min": (["analysis"],                              "run_min"),
    "run_max": (["analysis"],                              "run_max"),
    "lumi": (["analysis"],                              "lumi"),
    "clusterescale": (["analysis"],                              "cluster_escale"),
    "clustereres": (["analysis"],                              "cluster_eres"),
    # Additive ET-dependent smearing: 3 data params + 3 MC params
    "cluster_eres_data_p0": (["analysis"],                       "cluster_eres_data_p0"),
    "cluster_eres_data_p1": (["analysis"],                       "cluster_eres_data_p1"),
    "cluster_eres_data_p2": (["analysis"],                       "cluster_eres_data_p2"),
    "cluster_eres_mc_p0":   (["analysis"],                       "cluster_eres_mc_p0"),
    "cluster_eres_mc_p1":   (["analysis"],                       "cluster_eres_mc_p1"),
    "cluster_eres_mc_p2":   (["analysis"],                       "cluster_eres_mc_p2"),
    "mbd_avg_sigma_max": (["analysis"],                              "mbd_avg_sigma_max"),
    "mbd_avg_sigma_min": (["analysis"],                              "mbd_avg_sigma_min"),
    "reweight": (["analysis", "unfold"],                              "reweight"),
    "resultit": (["analysis", "unfold"],                              "resultit"),
    "flat_prior": (["analysis", "unfold"],                            "flat_prior"),
    "tower_mask_on":   (["analysis"],                                 "tower_mask_on"),
    "tower_mask_file": (["analysis"],                                 "tower_mask_file"),
    "tower_mask_name": (["analysis"],                                 "tower_mask_name"),
    # T3 (di_fraction): consumed by oneforall_tree_double_dispatch.sh once
    # the dispatcher learns to read this YAML field instead of computing
    # DOUBLE_FRAC from run_min. Configs exist now; condor-side support TODO.
    "double_frac_override": (["analysis"],                            "double_frac_override"),
}


_METADATA_KEYS = {"name", "syst_type", "syst_role", "aggregate_only"}
_ROLE_BUCKETS = {"up", "down", "one_sided", "max"}
_MODE_TO_ALLOWED_ROLES = {
    "two_sided": {"up", "down"},
    "one_sided": {"one_sided"},
    "max": {"max"},
    "placeholder": set(),
}

# ---------------------------------------------------------------------------
# Per-period merge-feeder override bundles (60cm fiducial family).
#
# Any VARIANTS entry that does NOT pin its own run range is auto-expanded
# into 3 on-disk configs by generate_variants(): the bare name (all-range,
# inherits the base config) plus two merge-feeders with these overrides
# applied on top of the variant's own overrides. The per-event lumi_weight
# = lumi/lumi_target pre-scales MC so a plain hadd across periods reproduces
# the all-range expectation.
# ---------------------------------------------------------------------------
PER_PERIOD_OVERRIDES = {
    # Nominal flipped to all-z on 2026-04-22: beam-delivered lumi with
    # vertex_cut_truth=9999 so the MBD-eff truth denominator covers the full
    # beam-delivered sample (no truth-vtx fiducial). The old 60cm-fiducial
    # values (32.6574, 16.2735, lumi_target=48.9309) are retained in git
    # history if a 60cm cross-check is needed.
    "0rad": {
        "run_min": 47289,
        "run_max": 51274,
        "lumi": 47.2076,
        "lumi_target": 64.3718,
        "vertex_cut_truth": 9999.0,
        "truth_vertex_reweight_on": 1,
        "truth_vertex_reweight_file": "/sphenix/user/shuhangli/ppg12/efficiencytool/truth_vertex_reweight/output/0mrad/reweight.root",
    },
    "1p5mrad": {
        "run_min": 51274,
        "run_max": 54000,
        "lumi": 17.1642,
        "lumi_target": 64.3718,
        "vertex_cut_truth": 9999.0,
        "truth_vertex_reweight_on": 1,
        "truth_vertex_reweight_file": "/sphenix/user/shuhangli/ppg12/efficiencytool/truth_vertex_reweight/output/1p5mrad/reweight.root",
    },
}


def validate_variants() -> None:
    """Validate systematic metadata so downstream aggregation can rely on it."""
    for variant in VARIANTS:
        name = variant["name"]
        syst_type = variant.get("syst_type")
        syst_role = variant.get("syst_role")

        if (syst_type is None) != (syst_role is None):
            raise ValueError(
                f"Variant '{name}' must set both syst_type and syst_role together"
            )

        if syst_type is None:
            continue

        if syst_type not in SYST_TYPES:
            raise ValueError(f"Variant '{name}' has unknown syst_type '{syst_type}'")

        if syst_role not in _ROLE_BUCKETS:
            raise ValueError(f"Variant '{name}' has unsupported syst_role '{syst_role}'")

        mode = SYST_TYPES[syst_type]["mode"]
        allowed_roles = _MODE_TO_ALLOWED_ROLES[mode]
        if syst_role not in allowed_roles:
            raise ValueError(
                f"Variant '{name}' uses role '{syst_role}' but syst_type "
                f"'{syst_type}' expects one of {sorted(allowed_roles)}"
            )


validate_variants()


def apply_overrides(doc, overrides: dict) -> None:
    """Mutate *doc* (ruamel.yaml CommentedMap) with the given overrides."""
    for key, value in overrides.items():
        if key in _METADATA_KEYS:
            continue
        if key not in OVERRIDE_MAP:
            raise ValueError(f"Unknown override key: '{key}'")
        section_path, leaf = OVERRIDE_MAP[key]
        node = doc
        for part in section_path:
            node = node[part]
        node[leaf] = value

    # Auto-default lumi_target to lumi unless explicitly overridden. This
    # prevents a merge-feeder lumi_target inherited from the base config from
    # silently propagating to systematic variants — every variant is treated
    # as a per-period standalone (lumi_weight = 1) by default; only entries
    # that explicitly set lumi_target become merge-feeders.
    if "lumi_target" not in overrides:
        doc["analysis"]["lumi_target"] = doc["analysis"]["lumi"]


def _write_config(yaml, base_doc, outdir, name, overrides):
    """Deep-copy base_doc, apply overrides, set var_type=bdt_{name}, write."""
    doc = copy.deepcopy(base_doc)
    apply_overrides(doc, overrides)
    doc["output"]["var_type"] = f"bdt_{name}"
    outfile = os.path.join(outdir, f"config_bdt_{name}.yaml")
    with open(outfile, "w") as fh:
        yaml.dump(doc, fh)
    summary = ", ".join(f"{k}={v}" for k, v in overrides.items()) if overrides else "(no overrides)"
    print(f"  Wrote {outfile}  [{summary}]")


def generate_variants(base_config: str, outdir: str) -> None:
    yaml = YAML()
    yaml.preserve_quotes = True

    with open(base_config, "r") as fh:
        base_doc = yaml.load(fh)

    for variant in VARIANTS:
        name = variant["name"]
        # aggregate_only variants are registered for calc_syst_bdt.py's variant
        # map but produce no config file (e.g. synthetic variants whose result
        # ROOT file is built post-hoc by a separate macro).
        if variant.get("aggregate_only", False):
            continue
        overrides = {k: v for k, v in variant.items() if k not in _METADATA_KEYS}

        # Always write the bare-name config (all-range if not period-pinned,
        # or whatever the variant's own overrides specify).
        _write_config(yaml, base_doc, outdir, name, overrides)

        # Auto-expand into per-period merge-feeders UNLESS the variant has
        # already pinned its own run range. This applies to regular
        # systematic variants (tightbdt50, noniso04, etc.) which inherit
        # the all-range base; it does NOT apply to the ntbdtpair scan or
        # the allz_{0rad,1p5mrad,all} entries that specify their own run
        # range, nor to the allz parent entry (all-range but has run_min
        # in its overrides to pin the lumi semantics).
        is_period_pinned = "run_min" in overrides or "run_max" in overrides
        if not is_period_pinned:
            for period, period_ovr in PER_PERIOD_OVERRIDES.items():
                expanded_overrides = {**overrides, **period_ovr}
                expanded_name = f"{name}_{period}"
                _write_config(yaml, base_doc, outdir, expanded_name, expanded_overrides)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("base_config", help="Path to the nominal base YAML config")
    parser.add_argument("--outdir", default=None,
                        help="Output directory (default: same directory as base_config)")
    args = parser.parse_args()

    if not os.path.isfile(args.base_config):
        sys.exit(f"ERROR: base config not found: {args.base_config}")

    outdir = args.outdir if args.outdir else os.path.dirname(os.path.abspath(args.base_config))
    os.makedirs(outdir, exist_ok=True)

    print(f"Base config : {args.base_config}")
    print(f"Output dir  : {outdir}")
    print(f"Variants    : {len(VARIANTS)}")
    print()
    generate_variants(args.base_config, outdir)
    print()
    print("Done.")


if __name__ == "__main__":
    main()
