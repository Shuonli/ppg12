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
    # Sanity-check: identical to base, only var_type differs
    dict(name="nom",
         syst_type=None, syst_role=None),

    # Run range for data processing — cross-checks only, not systematics
    dict(name="0rad", run_min=47289, run_max=51274, lumi=32.6574,
         syst_type=None, syst_role=None),
    dict(name="all", run_min=47289, run_max=54000, lumi=49.562,
         syst_type=None, syst_role=None),
    dict(name="mbdrms01", mbd_avg_sigma_max=0.1,
        syst_type=None, syst_role=None),
    dict(name="mbdrms03", mbd_avg_sigma_max=0.3,
        syst_type=None, syst_role=None),
    dict(name="mbdrms05", mbd_avg_sigma_max=0.5,
        syst_type=None, syst_role=None),

    # Tight BDT cut variations  →  syst: tight_bdt (efficiency)
    dict(name="tightbdt50",   tight_bdt_min_intercept=0.70, tight_bdt_min_slope=-0.015, nt_bdt_max_intercept=0.70, nt_bdt_max_slope=-0.015,
         syst_type="tight_bdt", syst_role="down"),   # looser cut
    dict(name="tightbdt70",   tight_bdt_min_intercept=0.90, tight_bdt_min_slope=-0.015, nt_bdt_max_intercept=0.90, nt_bdt_max_slope=-0.015,
         syst_type="tight_bdt", syst_role="up"),     # tighter cut
    dict(name="tightbdtflat06",   tight_bdt_min_intercept=0.6, tight_bdt_min_slope=0, nt_bdt_max_intercept=0.6, nt_bdt_max_slope=0,
         syst_type=None, syst_role=None),

    #iso scale and shift
    dict(name="mciso_noscaleshift", mc_iso_scale=1.0, mc_iso_shift=0.0,
         syst_type=None, syst_role=None),
    dict(name="mciso_no_shift",  mc_iso_shift=0.0,
         syst_type=None, syst_role=None),
    dict(name="mciso_no_scale",  mc_iso_scale=1.0,
         syst_type=None, syst_role=None),
    # Non-tight BDT lower boundary  →  syst: nt_bdt (purity)
    dict(name="ntbdtmin05",   nt_bdt_min=0.05,
         syst_type="nt_bdt", syst_role="one_sided"),
    dict(name="ntbdtmin10",   nt_bdt_min=0.10,
         syst_type="nt_bdt", syst_role="one_sided"),

    #common cuts
    dict(name="common_wr03",   common_wr_cogx_bound=0.3,
         syst_type=None, syst_role=None),
    dict(name="common_wr05",   common_wr_cogx_bound=0.5,
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

    # Split cluster node — cross-check only
    dict(name="split",   cluster_node_name="CLUSTERINFO_CEMC", data_file="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/condorout/part_*_with_bdt_split.root",
          photon_jet_file_branch_dir="/bdt_split.root",
          syst_type=None, syst_role=None),

    # Vertex reweighting off  →  syst: vtx_reweight (efficiency)
    dict(name="vtxreweight0", vertex_reweight_on=0,
         syst_type="vtx_reweight", syst_role="one_sided"),

    # no unfolding reweighting
    dict(name="no_unfolding_reweighting", reweight=0,
         syst_type="reweight", syst_role="one_sided"),

    # ET-binned BDT model variants  →  syst: bdt_model (efficiency)
    dict(name="etbin_v3E_v3E",  bdt_et_bin_edges=[8, 15, 35], bdt_et_bin_models=["base_v3E", "base_v3E"],
         syst_type="bdt_model", syst_role="max"),
    dict(name="etbin_E_E",      bdt_et_bin_edges=[8, 15, 35], bdt_et_bin_models=["base_E",   "base_E"],
         syst_type=None, syst_role=None),
    dict(name="etbin_v1E_v3E",  bdt_et_bin_edges=[8, 15, 35], bdt_et_bin_models=["base_v1E", "base_v3E"],
         syst_type=None, syst_role=None),

    # Back-to-back jet cut  →  syst: b2bjet (efficiency)
    dict(name="b2bjet", common_b2bjet_cut=1,
         syst_type=None, syst_role=None),
     

    # Timing cut variations  →  syst: timing (efficiency)
    dict(name="timingcut_2", cluster_mbd_time_min=-2.0, cluster_mbd_time_max=2.0,
         syst_type=None, syst_role=None),   # tighter timing window
    dict(name="timingcut_5", cluster_mbd_time_min=-5.0, cluster_mbd_time_max=5.0,
         syst_type=None, syst_role=None),     # looser timing window

    # Energy scale and resolution  →  syst: escale, eres
    dict(name="energyscale26up",   clusterescale=1.026,
         syst_type="escale", syst_role="down"),
    dict(name="energyscale26down", clusterescale=0.974,
         syst_type="escale", syst_role="up"),
    dict(name="energyresolution7", clustereres=0.07,
         syst_type=None, syst_role=None),
    dict(name="energyresolution8", clustereres=0.08,
         syst_type=None, syst_role=None),
    dict(name="energyresolution5", clustereres=0.05,
         syst_type="eres", syst_role="max"),
]

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
    "eff":    ["tight_bdt", "npb_cut"],
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
    "npb_score_cut":         (["analysis", "common"],                 "npb_score_cut"),
    "common_wr_cogx_bound": (["analysis", "common"],                 "wr_cogx_bound"),
    "mc_iso_scale":          (["analysis"],                           "mc_iso_scale"),
    "reco_iso_max_b":        (["analysis"],                           "reco_iso_max_b"),
    "reco_iso_max_s":        (["analysis"],                           "reco_iso_max_s"),
    "vertex_cut":            (["analysis"],                           "vertex_cut"),
    "vertex_reweight_on":    (["analysis"],                           "vertex_reweight_on"),
    "cluster_node_name":     (["input"],                              "cluster_node_name"),
    "data_file":             (["input"],                              "data_file"),
    "photon_jet_file_branch_dir": (["input"],                       "photon_jet_file_branch_dir"),
    "common_b2bjet_cut":        (["analysis"],                 "common_b2bjet_cut"),
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
}


_METADATA_KEYS = {"name", "syst_type", "syst_role"}
_ROLE_BUCKETS = {"up", "down", "one_sided", "max"}
_MODE_TO_ALLOWED_ROLES = {
    "two_sided": {"up", "down"},
    "one_sided": {"one_sided"},
    "max": {"max"},
    "placeholder": set(),
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


def generate_variants(base_config: str, outdir: str) -> None:
    yaml = YAML()
    yaml.preserve_quotes = True

    with open(base_config, "r") as fh:
        base_doc = yaml.load(fh)

    for variant in VARIANTS:
        name = variant["name"]
        overrides = {k: v for k, v in variant.items() if k not in _METADATA_KEYS}

        # Deep-copy preserving ruamel.yaml comment objects
        doc = copy.deepcopy(base_doc)

        # Apply parameter overrides
        apply_overrides(doc, overrides)

        # Set output.var_type to the variant name
        doc["output"]["var_type"] = f"bdt_{name}"

        outfile = os.path.join(outdir, f"config_bdt_{name}.yaml")
        with open(outfile, "w") as fh:
            yaml.dump(doc, fh)

        changed = [k for k in overrides]
        summary = ", ".join(f"{k}={overrides[k]}" for k in changed) if changed else "(no overrides)"
        print(f"  Wrote {outfile}  [{summary}]")


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
