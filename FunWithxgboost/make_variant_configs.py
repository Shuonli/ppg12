"""
Generate variant_configs/config_{model}_{variant}.yaml for all 11 photon-ID
model variants. Supports both "split" and "nosplit" cluster-node inputs.

Run once locally before submitting Condor jobs, e.g.:
    python make_variant_configs.py --variant split
    python make_variant_configs.py --variant nosplit
"""
import argparse
import copy
import os
import yaml

BASE_CONFIG = "config.yaml"
OUT_DIR = "variant_configs"

# Feature sets ordered to match apply_BDT.C x_list exactly
MODEL_FEATURES = {
    "base": [
        "vertexz", "cluster_Eta", "e11_over_e33",
        "cluster_et1", "cluster_et2", "cluster_et3", "cluster_et4",
    ],
    "base_vr": [
        "vertexz", "cluster_Eta", "e11_over_e33",
        "cluster_et1", "cluster_et2", "cluster_et3", "cluster_et4",
    ],
    "base_v0": [
        "vertexz", "cluster_Eta", "e11_over_e33",
        "cluster_et2", "cluster_et3", "cluster_et4",
    ],
    "base_v1": [
        "cluster_weta_cogx", "vertexz", "cluster_Eta", "e11_over_e33",
        "cluster_et1", "cluster_et2", "cluster_et3", "cluster_et4",
    ],
    "base_v2": [
        "cluster_weta_cogx", "cluster_wphi_cogx", "vertexz", "cluster_Eta", "e11_over_e33",
        "cluster_et1", "cluster_et2", "cluster_et3", "cluster_et4",
    ],
    "base_v3": [
        "cluster_weta_cogx", "cluster_wphi_cogx", "vertexz", "cluster_Eta", "e11_over_e33",
        "cluster_et1", "cluster_et2", "cluster_et3", "cluster_et4", "e32_over_e35",
    ],
    "base_E": [
        "cluster_Et", "vertexz", "cluster_Eta", "e11_over_e33",
        "cluster_et1", "cluster_et2", "cluster_et3", "cluster_et4",
    ],
    "base_v0E": [
        "cluster_Et", "vertexz", "cluster_Eta", "e11_over_e33",
        "cluster_et2", "cluster_et3", "cluster_et4",
    ],
    "base_v1E": [
        "cluster_Et", "cluster_weta_cogx", "vertexz", "cluster_Eta", "e11_over_e33",
        "cluster_et1", "cluster_et2", "cluster_et3", "cluster_et4",
    ],
    "base_v2E": [
        "cluster_Et", "cluster_weta_cogx", "cluster_wphi_cogx", "vertexz", "cluster_Eta", "e11_over_e33",
        "cluster_et1", "cluster_et2", "cluster_et3", "cluster_et4",
    ],
    "base_v3E": [
        "cluster_Et", "cluster_weta_cogx", "cluster_wphi_cogx", "vertexz", "cluster_Eta", "e11_over_e33",
        "cluster_et1", "cluster_et2", "cluster_et3", "cluster_et4", "e32_over_e35",
    ],
}


def build_file_lists(variant: str):
    """Return (signal_files, background_files) for the given variant."""
    signal_files = [
        f"shapes_{variant}_photon20.txt",
        f"shapes_{variant}_photon10.txt",
        f"shapes_{variant}_photon5.txt",
    ]
    background_files = [
        f"shapes_{variant}_jet40.txt",
        f"shapes_{variant}_jet30.txt",
        f"shapes_{variant}_jet20.txt",
        f"shapes_{variant}_jet12.txt",
        f"shapes_{variant}_jet5.txt",
    ]
    return signal_files, background_files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--variant",
        choices=["split", "nosplit"],
        required=True,
        help="Cluster-node variant: 'split' (CLUSTERINFO_CEMC) or 'nosplit' (CLUSTERINFO_CEMC_NO_SPLIT)",
    )
    args = parser.parse_args()
    variant = args.variant

    with open(BASE_CONFIG) as f:
        base = yaml.safe_load(f)

    os.makedirs(OUT_DIR, exist_ok=True)

    signal_files, background_files = build_file_lists(variant)

    filenames = []
    for name, features in MODEL_FEATURES.items():
        cfg = copy.deepcopy(base)

        cfg["data"]["features"] = features
        cfg["data"]["single_file_set"]["signal"] = signal_files
        cfg["data"]["single_file_set"]["background"] = background_files
        cfg["data"]["use_single_file_set"] = True
        cfg["output"]["root_file_prefix"] = f"model_{name}_{variant}"

        fname = f"config_{name}_{variant}.yaml"
        out_path = os.path.join(OUT_DIR, fname)
        with open(out_path, "w") as f:
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        filenames.append(fname)
        print(f"Wrote {out_path}")

    filelist_path = os.path.join(OUT_DIR, f"filelist_{variant}.txt")
    with open(filelist_path, "w") as f:
        f.write("\n".join(filenames) + "\n")
    print(f"Wrote {filelist_path}")

    print(f"\nDone — {len(MODEL_FEATURES)} {variant} configs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
