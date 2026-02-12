#!/bin/bash

# Wrapper script to submit PlotSpectra jobs to Condor
# Usage: ./submit_plot_spectra.sh [method]
#   method: "condor" (default) or "local"

METHOD=${1:-"condor"}

echo "=========================================="
echo "PlotSpectra Job Submission Script"
echo "Method: $METHOD"
echo "=========================================="

# Create logs directory if it doesn't exist
if [ ! -d "logs" ]; then
    echo "Creating logs directory..."
    mkdir -p logs
fi

case $METHOD in
    "condor")
        echo ""
        echo "Submitting jobs to Condor..."
        echo ""

        condor_submit run_plot_spectra.sub

        if [ $? -eq 0 ]; then
            echo ""
            echo "=========================================="
            echo "Jobs submitted successfully!"
            echo "Monitor with: condor_q"
            echo "Check logs in: logs/plot_spectra_*.{out,err,log}"
            echo "=========================================="
        else
            echo "ERROR: Job submission failed"
            exit 1
        fi
        ;;

    "local")
        echo ""
        echo "Running locally (sequential)..."
        echo "This will take a long time!"
        echo ""
        read -p "Continue? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            ./run_plot_spectra.sh all
        else
            echo "Cancelled."
            exit 0
        fi
        ;;

    *)
        echo "ERROR: Unknown method '$METHOD'"
        echo "Valid methods: condor, local"
        exit 1
        ;;
esac
