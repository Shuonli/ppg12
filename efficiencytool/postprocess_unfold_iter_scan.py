#!/usr/bin/env python3
"""Post-process the nominal Photon_final outputs into iter-scan variants.

CalculatePhotonYield.C writes h_unfold_sub_leak_{1..10} (and the un-leakage
versions) at every Bayes iteration. The "result" key h_unfold_sub_result is
just a copy of h_unfold_sub_leak_{resultit}. So for the unfold_iter
systematic scan (T1), we don't need to re-run the pipeline — we can just
read the existing nominal outputs and promote h_unfold_sub_leak_{N} into
the h_unfold_sub_result slot in new files named after the variant.

This is much cheaper than running the full pipeline 3x with different
resultit configs (the response matrix is identical across iterations).

Usage: python3 postprocess_unfold_iter_scan.py [results_dir]
"""
import sys
import os
import shutil
import ROOT

ROOT.gROOT.SetBatch(True)

NOMINAL = "bdt_nom"
ITERS   = [1, 3, 4]   # nominal is 2; scan these iter values as the syst envelope

# Histograms to re-point per iteration. Pairs of (target_key, source_key_template).
# target_key is what the new file will contain at top level; source_key_template
# is the key in the nominal file we copy from (with {N} replaced by the iter).
PROMOTIONS = [
    ("h_unfold_sub_result",       "h_unfold_sub_leak_{N}"),
    ("h_unfold_sub_result_woeff", "h_unfold_sub_{N}"),
]


def promote_iter(src_path: str, dst_path: str, iter_n: int) -> None:
    """Copy src_path -> dst_path, then replace target keys with iter_n versions."""
    if not os.path.exists(src_path):
        print(f"  [SKIP] {src_path} missing")
        return
    shutil.copy2(src_path, dst_path)

    f = ROOT.TFile.Open(dst_path, "UPDATE")
    if not f or f.IsZombie():
        print(f"  [FAIL] cannot open {dst_path} for UPDATE")
        return

    for target_key, source_template in PROMOTIONS:
        source_key = source_template.format(N=iter_n)
        h_src = f.Get(source_key)
        if not h_src:
            print(f"  [WARN] {dst_path}: missing source key {source_key}; "
                  f"target {target_key} unchanged")
            continue
        h_new = h_src.Clone(target_key)
        h_new.Write(target_key, ROOT.TObject.kOverwrite)
    f.Close()


def main():
    results_dir = sys.argv[1] if len(sys.argv) > 1 else \
        "/sphenix/user/shuhangli/ppg12/efficiencytool/results"
    if not os.path.isdir(results_dir):
        sys.exit(f"ERROR: results dir not found: {results_dir}")

    nom_data = os.path.join(results_dir, f"Photon_final_{NOMINAL}.root")
    nom_mc   = os.path.join(results_dir, f"Photon_final_{NOMINAL}_mc.root")
    for f in (nom_data, nom_mc):
        if not os.path.exists(f):
            sys.exit(f"ERROR: nominal file missing: {f}")

    print(f"results dir : {results_dir}")
    print(f"source data : {nom_data}")
    print(f"source mc   : {nom_mc}")
    print(f"iters       : {ITERS}")
    print()

    for n in ITERS:
        var = f"unfold_iter{n}"
        dst_data = os.path.join(results_dir, f"Photon_final_bdt_{var}.root")
        dst_mc   = os.path.join(results_dir, f"Photon_final_bdt_{var}_mc.root")
        print(f"  iter={n}: writing {dst_data} and {dst_mc}")
        promote_iter(nom_data, dst_data, n)
        promote_iter(nom_mc,   dst_mc,   n)

    print("\nDone.")


if __name__ == "__main__":
    main()
