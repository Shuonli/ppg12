#!/usr/bin/bash

# Script to clean up condor output directories and rerun condor jobs
# Iterates through jet and photon directories, removes OutDir* folders, and runs condor scripts

base_dir="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28"
# Full re-production after iso cone fix (commit 89941f1) — all PPG12 samples
dirs=("photon5" "photon10" "photon20" "photon10_double" "jet5" "jet8" "jet10" "jet12" "jet12_double" "jet15" "jet20" "jet30" "jet40" "jet50")

echo "=========================================="
echo "Starting cleanup and condor job submission"
echo "Base directory: $base_dir"
echo "=========================================="
echo ""

for dir in "${dirs[@]}"; do
    echo "Processing: $dir"

    # Check if directory exists
    if [ ! -d "$base_dir/$dir" ]; then
        echo "  WARNING: Directory $dir does not exist, skipping..."
        continue
    fi

    # Change to the target directory
    cd "$base_dir/$dir" || {
        echo "  ERROR: Failed to change to directory $dir"
        continue
    }

    # Check if condorout directory exists
    if [ ! -d "condorout" ]; then
        echo "  WARNING: condorout directory does not exist in $dir, skipping cleanup..."
    else
        # Remove all OutDir* directories and any stale shared submit file
        echo "  Removing OutDir* directories and ff.sub from condorout/..."
        rm -rf condorout/OutDir*
        rm -f condorout/ff.sub
        echo "  Cleanup complete for $dir"
    fi

    # Check if run_condor.sh exists
    if [ ! -f "run_condor.sh" ]; then
        echo "  WARNING: run_condor.sh not found in $dir, skipping job submission..."
    else
        # Run the condor script
        echo "  Running run_condor.sh for $dir..."
        ./run_condor.sh
        echo "  Job submission complete for $dir"
    fi

    echo "  Done with $dir"
    echo ""
done

echo "=========================================="
echo "All directories processed!"
echo "=========================================="
