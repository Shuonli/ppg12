#!/usr/bin/bash

# Script to clean up condor output directories
# Iterates through jet and photon directories and removes OutDir* folders

base_dir="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28"
dirs=("jet5" "jet12" "jet20" "jet30" "jet40" "photon10" "photon20" "photon5")

echo "=========================================="
echo "Starting cleanup"
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

    # Check if condorout directory exists
    if [ ! -d "$base_dir/$dir/condorout" ]; then
        echo "  WARNING: condorout directory does not exist in $dir, skipping..."
    else
        # Remove all OutDir* directories
        echo "  Removing OutDir* directories from condorout/..."
        rm -rf "$base_dir/$dir/condorout/OutDir"*
        echo "  Cleanup complete for $dir"
    fi

    echo "  Done with $dir"
    echo ""
done

echo "=========================================="
echo "All directories processed!"
echo "=========================================="
