#!/bin/bash

# Script to run PlotSpectra macro on multiple simulation samples
# Usage: ./run_plot_spectra.sh [mode]
#   mode: "all" (default), "photon", "jet", "nofilter"

MODE=${1:-"all"}

echo "=========================================="
echo "Running PlotSpectra on simulation samples"
echo "Mode: $MODE"
echo "=========================================="

# Define sample-specific pT ranges based on RecoEffCalculator_TTreeReader.C
# These are the default ranges for each sample

run_sample() {
    local sample=$1
    local min_photon=$2
    local max_photon=$3
    local min_jet=$4
    local max_jet=$5

    echo ""
    echo "=========================================="
    echo "Processing: $sample"
    echo "Photon filter: $min_photon - $max_photon GeV"
    echo "Jet filter: $min_jet - $max_jet GeV"
    echo "=========================================="

    root -l -b -q "PlotSpectra.C(\"$sample\", $min_photon, $max_photon, $min_jet, $max_jet)"

    if [ $? -eq 0 ]; then
        echo "SUCCESS: $sample completed"
    else
        echo "ERROR: $sample failed"
        return 1
    fi
}

# Photon samples
run_photon_samples() {
    echo ""
    echo "===== Processing Photon Samples ====="

    # photon5: 0-14 GeV (use 5-12 for cleaner separation)
    run_sample "photon5" 5 12 -1 -1

    # photon10: 14-30 GeV (use 12-25)
    run_sample "photon10" 12 25 -1 -1

    # photon20: 30-200 GeV (use 25-100)
    run_sample "photon20" 25 100 -1 -1
}

# Jet samples
run_jet_samples() {
    echo ""
    echo "===== Processing Jet Samples ====="

    # jet10: 10-15 GeV
    run_sample "jet10" -1 -1 10 15

    # jet15: 15-20 GeV
    run_sample "jet15" -1 -1 15 20

    # jet20: 20-30 GeV
    run_sample "jet20" -1 -1 20 30

    # jet30: 30-50 GeV
    run_sample "jet30" -1 -1 30 50

    # jet50: 50-100 GeV
    run_sample "jet50" -1 -1 50 100
}

# Run without filters (all events passing vertex cut)
run_nofilter_samples() {
    echo ""
    echo "===== Processing All Samples (No pT Filters) ====="

    for sample in photon5 photon10 photon20 jet10 jet15 jet20 jet30 jet50; do
        run_sample "$sample" -1 -1 -1 -1
    done
}

# Main execution based on mode
case $MODE in
    "photon")
        run_photon_samples
        ;;
    "jet")
        run_jet_samples
        ;;
    "nofilter")
        run_nofilter_samples
        ;;
    "all")
        run_photon_samples
        run_jet_samples
        ;;
    *)
        echo "ERROR: Unknown mode '$MODE'"
        echo "Valid modes: all, photon, jet, nofilter"
        exit 1
        ;;
esac

echo ""
echo "=========================================="
echo "All jobs completed!"
echo "Output files: PlotSpectra_*.root"
echo "=========================================="

# List generated files
echo ""
echo "Generated files:"
ls -lh PlotSpectra_*.root 2>/dev/null || echo "No output files found"
