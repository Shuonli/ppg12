#!/usr/bin/env python3
"""Derive parametric iso cuts (intercept, slope) for each variant by invoking FindETCut.C."""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
DEFAULT_VARIANTS = [
    "donutFull_005", "donutFull_0075", "donutFull_02",
    "donutExcl_005", "donutExcl_0075", "donutExcl_01", "donutExcl_02",
]

LINE_RE = re.compile(
    r"DERIVED_ISO_CUT\s+var_type=(\S+)\s+target_eff=(\S+)\s+intercept=(\S+)\s+slope=(\S+)"
)


def run_findetcut(var_name: str, target_eff: float):
    cmd = ["root", "-l", "-b", "-q", f'FindETCut.C("{var_name}", {target_eff})']
    res = subprocess.run(cmd, cwd=str(HERE), capture_output=True, text=True, timeout=300)
    for line in res.stdout.splitlines():
        m = LINE_RE.search(line)
        if m and m.group(1) == var_name:
            return float(m.group(3)), float(m.group(4))
    tail = "\n".join((res.stdout + "\n" + res.stderr).splitlines()[-20:])
    print(f"[WARN] DERIVED_ISO_CUT line missing for {var_name}; last 20 lines:\n{tail}",
          file=sys.stderr)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", default=",".join(DEFAULT_VARIANTS),
                    help="Comma-separated variant names (without 'bdt_' prefix)")
    ap.add_argument("--target-eff", type=float, default=0.8)
    ap.add_argument("--outfile", default=str(HERE / "derived_iso_cuts.yaml"))
    args = ap.parse_args()

    variants = [v.strip() for v in args.variants.split(",") if v.strip()]
    results = {}
    for v in variants:
        infile = RESULTS / f"MC_efficiency_bdt_{v}.root"
        if not infile.exists():
            print(f"[WARN] skip {v}: {infile} not found", file=sys.stderr)
            continue
        print(f"[INFO] running FindETCut for {v} (target_eff={args.target_eff})")
        pair = run_findetcut(v, args.target_eff)
        if pair is None:
            continue
        results[v] = pair

    if not results:
        print("[ERROR] no variants produced results", file=sys.stderr)
        sys.exit(1)

    yaml = YAML()
    yaml.default_flow_style = False
    yaml.indent(mapping=2, sequence=4, offset=2)
    doc = CommentedMap()
    header = (
        f"Derived iso cuts (target efficiency {args.target_eff}) -- "
        f"see derive_iso_cuts.py"
    )
    for v, (b, s) in results.items():
        entry = CommentedMap()
        entry["reco_iso_max_b"] = float(b)
        entry["reco_iso_max_s"] = float(s)
        doc[v] = entry
    doc.yaml_set_start_comment(header)
    with open(args.outfile, "w") as f:
        yaml.dump(doc, f)
    print(f"[INFO] wrote {args.outfile}")

    print("\nSummary (target_eff = {:.2f}):".format(args.target_eff))
    print(f"{'variant':<22} {'intercept':>14} {'slope':>14}")
    for v, (b, s) in results.items():
        print(f"{v:<22} {b:>14.6f} {s:>14.6f}")

    print("\n# Paste into make_bdt_variations.py VARIANTS placeholders:")
    for v, (b, s) in results.items():
        print(f'# {v}:')
        print(f'#     "reco_iso_max_b": {b:.6f},')
        print(f'#     "reco_iso_max_s": {s:.6f},')


if __name__ == "__main__":
    main()
