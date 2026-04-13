"""
Generate split_configs/config_{model}_split.yaml for all 11 photon-ID model variants.
Run once locally before submitting Condor jobs.
"""
import copy
import os
import yaml

BASE_CONFIG = "config.yaml"
OUT_DIR = "split_configs"

SIGNAL_FILES = [
    "shapes_split_photon20.txt",
    "shapes_split_photon10.txt",
    "shapes_split_photon5.txt",
]

BACKGROUND_FILES = [
    "shapes_split_jet40.txt",
    "shapes_split_jet30.txt",
    "shapes_split_jet20.txt",
    "shapes_split_jet12.txt",
    "shapes_split_jet5.txt",
]

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


def main():
    with open(BASE_CONFIG) as f:
        base = yaml.safe_load(f)

    os.makedirs(OUT_DIR, exist_ok=True)

    filenames = []
    for name, features in MODEL_FEATURES.items():
        cfg = copy.deepcopy(base)

        cfg["data"]["features"] = features
        cfg["data"]["single_file_set"]["signal"] = SIGNAL_FILES
        cfg["data"]["single_file_set"]["background"] = BACKGROUND_FILES
        cfg["data"]["use_single_file_set"] = True
        cfg["output"]["root_file_prefix"] = f"model_{name}_split"

        fname = f"config_{name}_split.yaml"
        out_path = os.path.join(OUT_DIR, fname)
        with open(out_path, "w") as f:
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        filenames.append(fname)
        print(f"Wrote {out_path}")

    filelist_path = os.path.join(OUT_DIR, "filelist.txt")
    with open(filelist_path, "w") as f:
        f.write("\n".join(filenames) + "\n")
    print(f"Wrote {filelist_path}")

    print(f"\nDone — {len(MODEL_FEATURES)} configs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
