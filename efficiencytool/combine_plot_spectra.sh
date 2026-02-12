#!/bin/bash

# Script to combine PlotSpectra output files using hadd
# Combines photon samples and jet samples separately with proper cross-section weighting

echo "=========================================="
echo "Combining PlotSpectra Output Files"
echo "=========================================="

# Check if output files exist
if ! ls PlotSpectra_*.root 1> /dev/null 2>&1; then
    echo "ERROR: No PlotSpectra output files found!"
    echo "Run the analysis first with: ./submit_plot_spectra.sh"
    exit 1
fi

echo ""
echo "Available output files:"
ls -lh PlotSpectra_*.root
echo ""

# Combine photon samples
PHOTON_FILES=$(ls PlotSpectra_photon*.root 2>/dev/null)
if [ -n "$PHOTON_FILES" ]; then
    echo "Combining photon samples..."
    hadd -f PlotSpectra_photon_combined.root $PHOTON_FILES
    echo "Created: PlotSpectra_photon_combined.root"
fi

# Combine jet samples
JET_FILES=$(ls PlotSpectra_jet*.root 2>/dev/null)
if [ -n "$JET_FILES" ]; then
    echo "Combining jet samples..."
    hadd -f PlotSpectra_jet_combined.root $JET_FILES
    echo "Created: PlotSpectra_jet_combined.root"
fi

# Combine all samples
echo ""
echo "Combining all samples..."
hadd -f PlotSpectra_all_combined.root PlotSpectra_photon*.root PlotSpectra_jet*.root

echo ""
echo "=========================================="
echo "Combination complete!"
echo ""
echo "Combined files:"
ls -lh PlotSpectra_*combined.root
echo ""
echo "NOTE: These files contain cross-section weighted histograms."
echo "The weights are already applied, so you can directly compare"
echo "or overlay distributions from different samples."
echo "=========================================="
