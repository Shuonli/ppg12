#!/usr/bin/bash

# Script to hadd (merge) ROOT files from condor output directories
# Combines all OutDir*/caloana.root files into a single combined.root file

base_dir="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28"
hadd_tmp="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/hadd_tmp"
dirs=("photon5" "photon10" "photon20" "photon10_double" "jet5" "jet8" "jet10" "jet12" "jet12_double" "jet15" "jet20" "jet30" "jet40" "jet50")
echo "=========================================="
echo "Starting hadd for all directories"
echo "Base directory: $base_dir"
echo "Hadd temp directory: $hadd_tmp"
echo "=========================================="
echo ""

# Create hadd_tmp directory if it doesn't exist
if [ ! -d "$hadd_tmp" ]; then
    echo "Creating hadd_tmp directory: $hadd_tmp"
    mkdir -p "$hadd_tmp"
fi

for dir in "${dirs[@]}"; do
    echo "Processing: $dir"

    # Check if directory exists
    if [ ! -d "$base_dir/$dir/condorout" ]; then
        echo "  WARNING: condorout directory does not exist in $dir, skipping..."
        continue
    fi

    # Change to the condorout directory
    cd "$base_dir/$dir/condorout" || {
        echo "  ERROR: Failed to change to condorout directory in $dir"
        continue
    }

    # Check if there are any OutDir* directories with caloana.root files
    if ! ls OutDir*/caloana.root &>/dev/null; then
        echo "  WARNING: No caloana.root files found in OutDir* directories, skipping..."
        cd "$base_dir"
        continue
    fi

    # Count how many files will be merged
    file_count=$(ls OutDir*/caloana.root 2>/dev/null | wc -l)
    echo "  Found $file_count caloana.root files to merge"

    # Run hadd to combine all ROOT files
    echo "  Running hadd to create combined.root..."
    hadd -j 14 -k -f -d "$hadd_tmp" combined.root OutDir*/caloana.root

    if [ $? -eq 0 ]; then
        echo "  Successfully created combined.root for $dir"
        # Check the size of the output file
        if [ -f "combined.root" ]; then
            size=$(du -h combined.root | cut -f1)
            echo "  Output file size: $size"
        fi
    else
        echo "  ERROR: hadd failed for $dir"
    fi

    echo "  Done with $dir"
    echo ""
done

echo "=========================================="
echo "All directories processed!"
echo "=========================================="
