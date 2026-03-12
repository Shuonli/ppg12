#!/usr/bin/env bash
# run_truthjet_windows.sh
# Step 1: fills histograms for all ET windows, one ROOT job per MC sample (parallel).
# Step 2: hadd per-sample files into a single combined file.
# Step 3: plots from the combined ROOT file.
#
# Usage: ./run_truthjet_windows.sh [config_file]
# Example: ./run_truthjet_windows.sh config_bdt_nom.yaml

source /sphenix/u/shuhang98/setup.sh

CONFIGNAME=${1:-config_bdt_nom.yaml}
LOGDIR="logs/truthjet"
OUTDIR="results"
PLOTDIR="truthjet_plots"
COMBINED="${OUTDIR}/truthjet_ETwindows_all.root"
mkdir -p "${LOGDIR}" "${OUTDIR}" "${PLOTDIR}"

echo "Using config: ${CONFIGNAME}"

echo "--- Step 1: filling histograms (jet samples in parallel) ---"
root -l -b -q 'TruthJetInETWindow.C("'"${CONFIGNAME}"'","jet5")'     > "${LOGDIR}/jet5.log"     2>&1 &
root -l -b -q 'TruthJetInETWindow.C("'"${CONFIGNAME}"'","jet12")'    > "${LOGDIR}/jet12.log"    2>&1 &
root -l -b -q 'TruthJetInETWindow.C("'"${CONFIGNAME}"'","jet20")'    > "${LOGDIR}/jet20.log"    2>&1 &
root -l -b -q 'TruthJetInETWindow.C("'"${CONFIGNAME}"'","jet30")'    > "${LOGDIR}/jet30.log"    2>&1 &
root -l -b -q 'TruthJetInETWindow.C("'"${CONFIGNAME}"'","jet40")'    > "${LOGDIR}/jet40.log"    2>&1

wait
echo "All samples done."

# Check that all expected output files exist
SAMPLES=(jet5 jet12 jet20 jet30 jet40)
MISSING=0
INFILES=()
for s in "${SAMPLES[@]}"; do
    f="${OUTDIR}/truthjet_ETwindows_${s}.root"
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing output file $f (check ${LOGDIR}/${s}.log)"
        MISSING=$((MISSING + 1))
    else
        INFILES+=("$f")
    fi
done
if [[ $MISSING -gt 0 ]]; then
    echo "Aborting: $MISSING sample(s) failed. Fix errors and re-run."
    exit 1
fi

echo "--- Step 2: hadding ---"
hadd -f "${COMBINED}" "${INFILES[@]}"

if [[ $? -ne 0 ]]; then
    echo "ERROR: hadd failed. Aborting."
    exit 1
fi

echo "--- Step 3: plotting ---"
root -l -b -q "PlotTruthJetWindows.C(\"${COMBINED}\",\"${PLOTDIR}\")" \
    2>&1 | tee "${LOGDIR}/plot.log"

echo "Done. Plots in ${PLOTDIR}/"
