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
    dict(name="mciso_no_shift",  mc_iso_shift=0.0,
         syst_type=None, syst_role=None),
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

    # NPB score cut  →  syst: npb_cut (efficiency)
    dict(name="npb03",        npb_score_cut=0.3,
         syst_type="npb_cut", syst_role="down"),  # looser npb cut
    dict(name="npb07",        npb_score_cut=0.7,
         syst_type="npb_cut", syst_role="up"),    # tighter npb cut

    # Purity fit option  →  syst: purity_fit (purity)
    dict(name="purity_pade", fit_option=0,
         syst_type="purity_fit", syst_role="one_sided"),
    # Use MC-driven purity correction ratio for data purity fit
    dict(name="mc_purity_correction", mc_purity_correction=1,
         syst_type="mc_purity_correction", syst_role="one_sided"),

    # Vertex reweighting off  →  syst: vtx_reweight (efficiency)
    dict(name="vtxreweight0", vertex_reweight_on=0,
         syst_type="vtx_reweight", syst_role="one_sided"),

    # no unfolding reweighting
    dict(name="no_unfolding_reweighting", reweight=0,
         syst_type="reweight", syst_role="one_sided"),

    # BDT model cross-checks — single model across all ET bins
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
    dict(name="energyscale26up",   clusterescale=1.026,
         syst_type="escale", syst_role="down"),
    dict(name="energyscale26down", clusterescale=0.974,
         syst_type="escale", syst_role="up"),
    dict(name="energyresolution5", clustereres=0.05,
         syst_type="eres", syst_role="max"),
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
        # A+X = t80_70_h80_70_l70_40 is the "symmetric-boundary" ABCD (nt_hi ==
        # tight_lo). Promoted to the sole nt_bdt one-sided systematic against the
        # all-range nominal. Other anchor combos remain cross-checks. Variants are
        # auto-expanded (bare + _0rad + _1p5mrad) via PER_PERIOD_OVERRIDES.
        if _low_label == "A" and _high_label == "X":
            _syst_type, _syst_role = "nt_bdt", "one_sided"
        else:
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
    "reweight":     {"mode": "one_sided",   "group": "unfold"},
    "tight_bdt":    {"mode": "two_sided",   "group": "eff"},
    "nt_bdt":       {"mode": "one_sided",   "group": "purity"},
    "noniso":       {"mode": "two_sided",   "group": "purity"},
    "npb_cut":      {"mode": "two_sided",   "group": "eff"},
    "purity_fit":   {"mode": "one_sided",   "group": "purity"},
    "mc_purity_correction": {"mode": "one_sided", "group": "purity"},
    "vtx_reweight": {"mode": "one_sided",   "group": "eff"},
    "bdt_model":    {"mode": "max",         "group": "eff"},
    "b2bjet":       {"mode": "one_sided",   "group": "eff"},
    "timing":       {"mode": "two_sided",   "group": "eff"},
    "escale":       {"mode": "two_sided",   "group": "escale"},
    "eres":         {"mode": "max",         "group": "eres"},
    "mbd":          {"mode": "placeholder", "group": "mbd"},   # TODO: add mbdeffup/down variants
    "nor":          {"mode": "placeholder", "group": "unfolding"},   # TODO: add nr variant
}

# Quadrature grouping: group name -> list of syst_type names
SYST_GROUPS = {
    "purity": ["noniso", "nt_bdt", "purity_fit", "mc_purity_correction"],
    "eff":    ["npb_cut"],   # tight_bdt removed 2026-04-24 with tightbdt50/70 cleanup — full syst needs redo
    "escale": ["escale"],
    "eres":   ["eres"],
    "mbd":    ["mbd"],
    "unfolding": ["reweight"],
}

# Groups included in the final total systematic quadrature sum
FINAL_SYSTS = ["purity", "eff", "escale", "eres", "mbd", "unfolding"]

# Luminosity uncertainty (asymmetric, fractional): central=25.2, down=23.5, up=27.5 pb^-1
LUMI_SYST = {"down": (25.2 - 23.5) / 25.2, "up": (27.5 - 25.2) / 25.2}

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
    "mc_purity_correction": (["analysis"],              "mc_purity_correction"),
    "bdt_et_bin_edges":  (["input"],                              "bdt_et_bin_edges"),
    "bdt_et_bin_models": (["input"],                              "bdt_et_bin_models"),
    "run_min": (["analysis"],                              "run_min"),
    "run_max": (["analysis"],                              "run_max"),
    "lumi": (["analysis"],                              "lumi"),
    "clusterescale": (["analysis"],                              "cluster_escale"),
    "clustereres": (["analysis"],                              "cluster_eres"),
    "mbd_avg_sigma_max": (["analysis"],                              "mbd_avg_sigma_max"),
    "mbd_avg_sigma_min": (["analysis"],                              "mbd_avg_sigma_min"),
    "reweight": (["analysis", "unfold"],                              "reweight"),
    "tower_mask_on":   (["analysis"],                                 "tower_mask_on"),
    "tower_mask_file": (["analysis"],                                 "tower_mask_file"),
    "tower_mask_name": (["analysis"],                                 "tower_mask_name"),
}


_METADATA_KEYS = {"name", "syst_type", "syst_role"}
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
