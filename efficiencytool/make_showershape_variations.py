#!/usr/bin/env python3
"""Generate shower-shape variation config files from a nominal base config.

Usage:
    python make_showershape_variations.py <base_config> [--outdir DIR]

Each entry in VARIANTS produces one output file named
config_showershape_{name}.yaml in the output directory (defaults to the same
directory as the base config).

Supported override keys and their YAML paths:
    vertex_cut             -> analysis.vertex_cut
    run_min                -> analysis.run_min
    run_max                -> analysis.run_max
    lumi                   -> analysis.lumi
    run_list_file          -> analysis.run_list_file
    cluster_mbd_time_min   -> analysis.cluster_mbd_time_min
    cluster_mbd_time_max   -> analysis.cluster_mbd_time_max
    mbd_avgsigma_cut_on    -> analysis.mbd_avgsigma_cut_on
    mbd_avgsigma_max       -> analysis.mbd_avgsigma_max
    mc_iso_scale           -> analysis.mc_iso_scale
    mc_iso_shift           -> analysis.mc_iso_shift
    vertex_reweight_on     -> analysis.vertex_reweight_on
    npb_score_cut          -> analysis.common.npb_score_cut
    reco_noniso_min_shift  -> analysis.reco_noniso_min_shift
    cluster_node_name      -> input.cluster_node_name
    data_file              -> input.data_file
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
    # Sanity-check: identical to base
    dict(name="nom"),

    # Run period: 0rad (run_min: 47289, run_max: 51274, lumi: 32.6574 pb^-1)
    dict(name="0rad",
         run_min=47289, run_max=51274, lumi=32.6574),

    # Run period: 1p5rad (run_min: 51274, run_max: 54000, lumi: 16.8588 pb^-1)
    dict(name="1p5rad",
         run_min=51274, run_max=54000, lumi=16.8588),

    # 0rad period with tight timing cut (±2 ns cluster-MBD window)
    dict(name="0radt2",
         run_min=47289, run_max=51274, lumi=32.6574,
         cluster_mbd_time_min=-4.0, cluster_mbd_time_max=0.0),
    
    #back to back jet cut
    dict(name="b2bjet", common_b2bjet_cut=1),

    # MBD avg-sigma pileup rejection (avgsigma < 0.1)
    dict(name="mbd01",
         mbd_avgsigma_cut_on=1, mbd_avgsigma_max=0.1),

    # Low-pileup run list with tight vertex cut
    dict(name="0radlowpu",
         run_list_file="/sphenix/user/shuhangli/ppg12/efficiencytool/low_pileup.list"),

    # Low-pileup run list (800 MeV threshold)
    dict(name="0radlowpu800",
         run_list_file="/sphenix/user/shuhangli/ppg12/efficiencytool/low_pileup_800.list"),
    dict(name="etbin_v3E_v3E",  bdt_et_bin_edges=[8, 15, 35], bdt_et_bin_models=["base_v3E", "base_v3E"]),
]

# ---------------------------------------------------------------------------
# Mapping from flat override key -> (yaml_section_path, yaml_leaf_key)
# ---------------------------------------------------------------------------
OVERRIDE_MAP = {
    "vertex_cut":            (["analysis"],          "vertex_cut"),
    "run_min":               (["analysis"],          "run_min"),
    "run_max":               (["analysis"],          "run_max"),
    "lumi":                  (["analysis"],          "lumi"),
    "run_list_file":         (["analysis"],          "run_list_file"),
    "cluster_mbd_time_min":  (["analysis"],          "cluster_mbd_time_min"),
    "cluster_mbd_time_max":  (["analysis"],          "cluster_mbd_time_max"),
    "mbd_avgsigma_cut_on":   (["analysis"],          "mbd_avgsigma_cut_on"),
    "mbd_avgsigma_max":      (["analysis"],          "mbd_avgsigma_max"),
    "bdt_et_bin_edges":      (["analysis"],          "bdt_et_bin_edges"),
    "bdt_et_bin_models":     (["analysis"],          "bdt_et_bin_models"),
    "mc_iso_scale":          (["analysis"],          "mc_iso_scale"),
    "mc_iso_shift":          (["analysis"],          "mc_iso_shift"),
    "vertex_reweight_on":    (["analysis"],          "vertex_reweight_on"),
    "reco_noniso_min_shift": (["analysis"],          "reco_noniso_min_shift"),
    "reco_iso_max_b":        (["analysis"],          "reco_iso_max_b"),
    "common_b2bjet_cut":     (["analysis"],          "common_b2bjet_cut"),
    "reco_iso_max_s":        (["analysis"],          "reco_iso_max_s"),
    "npb_score_cut":         (["analysis", "common"], "npb_score_cut"),
    "cluster_node_name":     (["input"],             "cluster_node_name"),
    "data_file":             (["input"],             "data_file"),
    "photon_jet_file_branch_dir": (["input"],        "photon_jet_file_branch_dir"),

}


def apply_overrides(doc, overrides: dict) -> None:
    """Mutate *doc* (ruamel.yaml CommentedMap) with the given overrides."""
    for key, value in overrides.items():
        if key == "name":
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
        overrides = {k: v for k, v in variant.items() if k != "name"}

        doc = copy.deepcopy(base_doc)

        apply_overrides(doc, overrides)

        # Set var_type: "bdt_showershape" for nom, "bdt_showershape_{name}" otherwise
        if name == "nom":
            doc["output"]["var_type"] = "bdt_showershape"
        else:
            doc["output"]["var_type"] = f"bdt_showershape_{name}"

        outfile = os.path.join(outdir, f"config_showershape_{name}.yaml")
        with open(outfile, "w") as fh:
            yaml.dump(doc, fh)

        changed = [k for k in overrides]
        summary = ", ".join(f"{k}={overrides[k]}" for k in changed) if changed else "(no overrides)"
        print(f"  Wrote {outfile}  [{summary}]")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("base_config", help="Path to the nominal base YAML config (config_showershape.yaml)")
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
    generate_variants(base_config=args.base_config, outdir=outdir)
    print()
    print("Done.")


if __name__ == "__main__":
    main()
