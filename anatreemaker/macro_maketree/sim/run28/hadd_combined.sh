#!/usr/bin/bash

# hadd OutDir*/caloana.root -> combined.root per sample.
# Covers all 18 MC samples (SI + DI + new DI) after the iso-cone fix + 9-branch add.
# Serial merge (-j 1 -n 50) to avoid /tmp exhaustion on large samples.
#
# Usage:
#   ./hadd_combined.sh                   # all 18 samples
#   ./hadd_combined.sh si                # 9 SI
#   ./hadd_combined.sh di                # 9 DI (old + new)
#   ./hadd_combined.sh di_new            # 7 new DI

base_dir="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28"
hadd_tmp="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/hadd_tmp"

SI=("photon5" "photon10" "photon20" "jet5" "jet8" "jet12" "jet20" "jet30" "jet40")
DI_OLD=("photon10_double" "jet12_double")
DI_NEW=("photon5_double" "photon20_double" "jet8_double" "jet20_double" "jet30_double" "jet40_double" "jet50_double")

mode="${1:-all}"
case "$mode" in
    si)     dirs=("${SI[@]}") ;;
    di)     dirs=("${DI_OLD[@]}" "${DI_NEW[@]}") ;;
    di_new) dirs=("${DI_NEW[@]}") ;;
    all)    dirs=("${SI[@]}" "${DI_OLD[@]}" "${DI_NEW[@]}") ;;
    *) echo "Unknown mode '$mode' (use: si | di | di_new | all)"; exit 1 ;;
esac

echo "=========================================="
echo "Mode: $mode   Samples: ${#dirs[@]}"
echo "Base: $base_dir   hadd_tmp: $hadd_tmp"
echo "=========================================="

mkdir -p "$hadd_tmp"

ok=0; failed=0; skipped=0
date
for dir in "${dirs[@]}"; do
    echo "=== $dir ==="
    cond="$base_dir/$dir/condorout"
    if [ ! -d "$cond" ]; then
        echo "  SKIP: $cond missing"; skipped=$((skipped+1)); continue
    fi
    cd "$cond" || { echo "  cd failed"; skipped=$((skipped+1)); continue; }

    files=(OutDir*/caloana.root)
    if [ ${#files[@]} -eq 0 ] || [ ! -f "${files[0]}" ]; then
        echo "  SKIP: no caloana.root"; skipped=$((skipped+1)); continue
    fi

    echo "  merging ${#files[@]} files (-j 1 -n 50) -> combined.root"
    rm -f combined.root
    hadd -f -j 1 -n 50 -d "$hadd_tmp" combined.root "${files[@]}" 2>&1 | tail -3
    if [ -f combined.root ]; then
        size=$(du -h combined.root | cut -f1)
        echo "  OK  -> combined.root ($size)"
        ok=$((ok+1))
    else
        echo "  FAIL no combined.root produced"
        failed=$((failed+1))
    fi
done

echo "=========================================="
echo "Summary: ok=$ok  failed=$failed  skipped=$skipped  (mode=$mode)"
date
echo "=========================================="
