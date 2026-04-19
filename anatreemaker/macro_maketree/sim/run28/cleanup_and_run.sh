#!/usr/bin/bash

# Clean + relaunch condor tree-making for ALL PPG12 MC samples
# (single-interaction + old DI + new DI expansion after iso cone fix + 9-branch add).
#
# Per-sample run_condor.sh: total_line=10000, lines_per_job=10 => 1000 jobs / sample.
# Sample groups:
#   SI  (9)  photon5, photon10, photon20, jet5, jet8, jet12, jet20, jet30, jet40
#   DI  (2)  photon10_double, jet12_double                 (existing, MDC2 cluster+mbd DST inputs)
#   DI+ (7)  photon5_double, photon20_double, jet8_double,
#            jet20_double, jet30_double, jet40_double, jet50_double
#            (new — MDC2-style macro, reads g4hits + truth_jet only)
#
# Usage:
#   ./cleanup_and_run.sh               # all 18 samples
#   ./cleanup_and_run.sh si            # only SI (9)
#   ./cleanup_and_run.sh di            # only DI — old + new (9)
#   ./cleanup_and_run.sh di_new        # only new DI (7)
#   ./cleanup_and_run.sh dry           # list what would run without submitting

base_dir="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28"

SI=("photon5" "photon10" "photon20" "jet5" "jet8" "jet12" "jet20" "jet30" "jet40")
DI_OLD=("photon10_double" "jet12_double")
DI_NEW=("photon5_double" "photon20_double" "jet8_double" "jet20_double" "jet30_double" "jet40_double" "jet50_double")

mode="${1:-all}"
case "$mode" in
    si)     dirs=("${SI[@]}") ;;
    di)     dirs=("${DI_OLD[@]}" "${DI_NEW[@]}") ;;
    di_new) dirs=("${DI_NEW[@]}") ;;
    all|dry) dirs=("${SI[@]}" "${DI_OLD[@]}" "${DI_NEW[@]}") ;;
    *) echo "Unknown mode '$mode' (use: si | di | di_new | all | dry)"; exit 1 ;;
esac

echo "=========================================="
echo "Mode: $mode       Samples: ${#dirs[@]}       (~$((${#dirs[@]}*1000)) jobs)"
echo "Base: $base_dir"
echo "=========================================="
echo ""

submitted=0
skipped=0
for dir in "${dirs[@]}"; do
    echo "Processing: $dir"

    if [ ! -d "$base_dir/$dir" ]; then
        echo "  WARNING: Directory $dir does not exist, skipping..."
        skipped=$((skipped+1)); continue
    fi

    cd "$base_dir/$dir" || { echo "  ERROR: cd $dir failed"; skipped=$((skipped+1)); continue; }

    # Input guard: every sample must have non-empty g4hits.list AND either
    # dst_calo_cluster.list (SI/old-DI path) or dst_truth_jet.list (new-DI MDC2 path).
    n_g4=$(wc -l < g4hits.list 2>/dev/null || echo 0)
    n_dst=$(wc -l < dst_calo_cluster.list 2>/dev/null || echo 0)
    n_tj=$(wc -l < dst_truth_jet.list 2>/dev/null || echo 0)
    if [ "$n_g4" -eq 0 ] || { [ "$n_dst" -eq 0 ] && [ "$n_tj" -eq 0 ]; }; then
        echo "  SKIP: empty input lists (g4=$n_g4 calo=$n_dst truthjet=$n_tj)"
        skipped=$((skipped+1)); continue
    fi

    if [ ! -f "run_condor.sh" ]; then
        echo "  SKIP: run_condor.sh not found"
        skipped=$((skipped+1)); continue
    fi

    if [ "$mode" = "dry" ]; then
        echo "  DRY: would clean condorout/OutDir* and submit run_condor.sh  (g4=$n_g4 calo=$n_dst truthjet=$n_tj)"
        submitted=$((submitted+1)); continue
    fi

    if [ -d "condorout" ]; then
        echo "  Removing condorout/OutDir* and stale ff.sub..."
        rm -rf condorout/OutDir*
        rm -f condorout/ff.sub
    fi

    echo "  Submitting run_condor.sh  (g4=$n_g4 calo=$n_dst truthjet=$n_tj)..."
    ./run_condor.sh
    submitted=$((submitted+1))
    echo "  Done with $dir"
    echo ""
done

echo "=========================================="
echo "Summary: submitted=$submitted  skipped=$skipped  (mode=$mode)"
echo "=========================================="
