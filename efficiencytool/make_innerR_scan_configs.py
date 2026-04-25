#!/usr/bin/env python3
"""
Generate 8 inner-R scan configs from config_bdt_nom.yaml:
  4 inner R values (0.05, 0.075, 0.10, 0.20) x 2 reweight states (on/off) = 8 configs.

Each variant sets:
  - analysis.iso_topo_innerr : inner R value (new key, inserted under analysis:)
  - output.var_type          : unique suffix (innerR_005, innerR_005_nore, ...)
  - analysis.unfold.reweight : 0 for _nore variants (1 otherwise, matching nominal)
"""
import argparse
from pathlib import Path
from ruamel.yaml import YAML

VARIANTS = [
    # (suffix, inner R, reweight)
    ("innerR_005",       0.05,  1),
    ("innerR_005_nore",  0.05,  0),
    ("innerR_0075",      0.075, 1),
    ("innerR_0075_nore", 0.075, 0),
    ("innerR_01",        0.10,  1),
    ("innerR_01_nore",   0.10,  0),
    ("innerR_02",        0.20,  1),
    ("innerR_02_nore",   0.20,  0),
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nom",   default="config_bdt_nom.yaml")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.indent(mapping=2, sequence=4, offset=2)

    nom_path = Path(args.nom)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for suffix, inner_r, reweight in VARIANTS:
        with nom_path.open() as f:
            cfg = yaml.load(f)

        # Set inner R on topo path
        cfg["analysis"]["iso_topo_innerr"] = float(inner_r)

        # Toggle reweight in unfold block (nominal is 1)
        cfg["analysis"]["unfold"]["reweight"] = int(reweight)

        # Unique var_type
        cfg["output"]["var_type"] = f"bdt_{suffix}"

        out_path = outdir / f"config_bdt_{suffix}.yaml"
        with out_path.open("w") as f:
            yaml.dump(cfg, f)
        print(f"  wrote {out_path}  (iso_topo_innerr={inner_r}, reweight={reweight}, var_type=bdt_{suffix})")

if __name__ == "__main__":
    main()
