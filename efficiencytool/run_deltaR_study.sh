#!/bin/bash
# Run the deltaR truth-matching study on single-interaction MC
# Usage: bash run_deltaR_study.sh [maxEvents]
#   maxEvents: optional, default -1 (all events)

set -e

MAXEVENTS=${1:--1}
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
PLOTDIR="$(cd "$SCRIPTDIR/../plotting" && pwd)"

echo "=== Phase 1: Running analysis macro ==="
cd "$SCRIPTDIR"
root -l -b -q "study_deltaR_single.C(${MAXEVENTS})"

echo "=== Phase 2: Making plots ==="
cd "$PLOTDIR"
root -l -b -q 'plot_deltaR_single_study.C'

echo "=== Done. Plots in plotting/figures/deltaR_study/ ==="
