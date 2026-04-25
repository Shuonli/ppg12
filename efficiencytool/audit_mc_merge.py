#!/usr/bin/env python3
"""Audit all-range merged MC efficiency files for silent merge failures.

The PPG12 pipeline merges per-period (_0rad, _1p5mrad) MC outputs into an
all-range file via merge_periods.C (TFileMerger). When TFileMerger hits a
transient gpfs I/O error on one of the input files, it can silently drop
that input's contents — the all-range file ends up with only one period's
events. The truth-level histogram h_all_cluster_signal_0 is the cleanest
indicator because it should equal (0rad + 1p5mrad) bin-for-bin, with no
NPB / iso / cut dependence.

Usage:
    python audit_mc_merge.py [--results DIR] [--fix] [--tolerance 0.01]

Default behaviour: print BROKEN variants and exit 0 if all are OK, exit 1
if any are broken. With --fix, automatically re-runs merge_periods.sh on
the broken variants (reads CFG names from the merged-file basename) and
schedules a Phase-2 condor re-submission script for CalculatePhotonYield.
"""
import argparse
import glob
import os
import subprocess
import sys

import uproot

TRUTH_KEY = "h_all_cluster_signal_0"  # NPB-/iso-independent, ideal mismatch detector


def audit_one(base_path: str, base_name: str, tolerance: float):
    """Return (status, exp_sum, actual_sum) where status is OK / BROKEN / MISSING_INPUT."""
    f0 = os.path.join(base_path, f"MC_efficiency_bdt_{base_name}_0rad.root")
    f1 = os.path.join(base_path, f"MC_efficiency_bdt_{base_name}_1p5mrad.root")
    fb = os.path.join(base_path, f"MC_efficiency_bdt_{base_name}.root")
    if not os.path.exists(f0) or not os.path.exists(f1) or not os.path.exists(fb):
        return ("MISSING_INPUT", 0.0, 0.0)
    try:
        s0 = uproot.open(f0)[TRUTH_KEY].values().sum()
        s1 = uproot.open(f1)[TRUTH_KEY].values().sum()
        sb = uproot.open(fb)[TRUTH_KEY].values().sum()
    except (KeyError, OSError) as e:
        return ("READ_ERROR", 0.0, 0.0)
    exp = s0 + s1
    if exp <= 0:
        return ("EMPTY", 0.0, sb)
    rel = (sb - exp) / exp
    return ("OK" if abs(rel) < tolerance else "BROKEN", exp, sb)


def find_variants(base_path: str):
    """Discover variants that have both _0rad and _1p5mrad MC efficiency files."""
    files0 = glob.glob(os.path.join(base_path, "MC_efficiency_bdt_*_0rad.root"))
    variants = []
    for f0 in files0:
        base = os.path.basename(f0).replace("MC_efficiency_bdt_", "").replace("_0rad.root", "")
        f1 = os.path.join(base_path, f"MC_efficiency_bdt_{base}_1p5mrad.root")
        fb = os.path.join(base_path, f"MC_efficiency_bdt_{base}.root")
        if os.path.exists(f1) and os.path.exists(fb):
            variants.append(base)
    return sorted(variants)


def fix_one(efficiencytool_dir: str, base_name: str):
    """Re-run merge_periods.sh for a broken variant."""
    cfg0 = f"config_bdt_{base_name}_0rad.yaml"
    cfg1 = f"config_bdt_{base_name}_1p5mrad.yaml"
    cfgb = f"config_bdt_{base_name}.yaml"
    cmd = ["bash", "merge_periods.sh", cfg0, cfg1, cfgb]
    print(f"  [fix] running: {' '.join(cmd)}")
    r = subprocess.run(cmd, cwd=efficiencytool_dir, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [fix] FAILED rc={r.returncode}: {r.stderr.splitlines()[-3:]}")
        return False
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", default="/sphenix/user/shuhangli/ppg12/efficiencytool/results",
                   help="Directory with MC_efficiency_bdt_*.root files")
    p.add_argument("--efficiencytool", default="/sphenix/user/shuhangli/ppg12/efficiencytool",
                   help="Directory with merge_periods.sh and configs")
    p.add_argument("--fix", action="store_true",
                   help="Re-run merge_periods.sh on broken variants. Does NOT re-run "
                        "CalculatePhotonYield — submit Phase-2 condor manually after.")
    p.add_argument("--tolerance", type=float, default=0.01,
                   help="Relative mismatch threshold (default: 1%%)")
    p.add_argument("--quiet", action="store_true", help="Print only broken variants")
    args = p.parse_args()

    variants = find_variants(args.results)
    print(f"Auditing {len(variants)} variants in {args.results} ...")
    broken, ok = [], []
    for v in variants:
        status, exp, actual = audit_one(args.results, v, args.tolerance)
        if status == "OK":
            ok.append(v)
            if not args.quiet:
                pass
        elif status == "BROKEN":
            rel = (actual - exp) / exp if exp > 0 else 0
            print(f"  BROKEN  {v:>40s}  exp={exp:.3e}  actual={actual:.3e}  rel={rel*100:+.1f}%")
            broken.append(v)
        else:
            print(f"  {status}  {v}")

    print(f"\nSummary: {len(ok)} OK, {len(broken)} BROKEN")

    if not broken:
        print("All merges consistent.")
        return 0

    if args.fix:
        print(f"\nRe-running merge_periods.sh on {len(broken)} broken variants...")
        fixed = 0
        for v in broken:
            if fix_one(args.efficiencytool, v):
                fixed += 1
        print(f"\n{fixed}/{len(broken)} re-merges succeeded.")
        print("\nNEXT STEP: re-run CalculatePhotonYield via condor for these variants.")
        print("Suggested submit-file body:")
        print("  queue cfg in (")
        for v in broken:
            print(f"      config_bdt_{v}.yaml")
        print("  )")
        return 0 if fixed == len(broken) else 2
    else:
        print("\nRun with --fix to auto-remerge. Then re-run CalculatePhotonYield via condor.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
