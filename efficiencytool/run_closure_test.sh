#!/bin/bash

# Script to run closure tests for unfolding validation
# Usage: ./run_closure_test.sh [config_file] [max_iterations] [save_all]
#   config_file: YAML config file (default: config_bdt_none.yaml)
#   max_iterations: Maximum number of iterations to test (default: 10)
#   save_all: Save all iteration histograms to ROOT file (default: false)

CONFIG=${1:-"config_bdt_none.yaml"}
MAX_ITER=${2:-10}
SAVE_ALL=${3:-false}

echo "========================================="
echo "Running Closure Test"
echo "========================================="
echo "Config file: $CONFIG"
echo "Max iterations: $MAX_ITER"
echo "Save all iterations: $SAVE_ALL"
echo ""

# Check if config file exists
if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config file $CONFIG not found!"
    exit 1
fi

# Check if Closure.C exists
if [ ! -f "Closure.C" ]; then
    echo "ERROR: Closure.C not found!"
    exit 1
fi

# Run the closure test
echo "Running ROOT macro..."
root -l -b -q "Closure.C(\"$CONFIG\", $MAX_ITER, $SAVE_ALL)"

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "Closure test completed successfully!"
    echo "Check output files:"
    echo "  - ClosureTest_${CONFIG%.yaml}.root"
    echo "  - figure/FullClosure_eta*_iter*.png"
    echo "  - figure/HalfClosure_eta*_iter*.png"
    echo "  - figure/FullClosure_chi2_eta*.png"
    echo "  - figure/HalfClosure_chi2_eta*.png"
    echo "  - figure/ClosureComparison_eta*.png"
    echo "========================================="
else
    echo ""
    echo "ERROR: Closure test failed!"
    exit 1
fi
