#!/bin/bash

# Condor worker script for PlotSpectra
# This script runs on condor nodes
# Arguments: $1=sample, $2=min_photon, $3=max_photon, $4=min_jet, $5=max_jet

SAMPLE=$1
MIN_PHOTON=$2
MAX_PHOTON=$3
MIN_JET=$4
MAX_JET=$5

echo "=========================================="
echo "Condor Job for PlotSpectra"
echo "Sample: $SAMPLE"
echo "Photon filter: $MIN_PHOTON - $MAX_PHOTON GeV"
echo "Jet filter: $MIN_JET - $MAX_JET GeV"
echo "Node: $(hostname)"
echo "Time: $(date)"
echo "=========================================="

# Set up ROOT environment (adjust path if needed)
source /opt/sphenix/core/bin/sphenix_setup.sh -n

# Navigate to work directory
cd /gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool

# Run the macro
root -l -b -q "PlotSpectra.C(\"$SAMPLE\", $MIN_PHOTON, $MAX_PHOTON, $MIN_JET, $MAX_JET)"

EXIT_CODE=$?

echo ""
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: Job completed"
    ls -lh PlotSpectra_${SAMPLE}*.root
else
    echo "ERROR: Job failed with exit code $EXIT_CODE"
fi
echo "=========================================="

exit $EXIT_CODE
