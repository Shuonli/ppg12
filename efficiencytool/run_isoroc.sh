#!/usr/bin/env bash
# run_isoroc.sh
# Runs IsoROC_calculator.C across all jet and photon MC samples in parallel.
# Usage: ./run_isoroc.sh [config_file]
# Example: ./run_isoroc.sh config_bdt_nom.yaml

source /sphenix/u/shuhang98/setup.sh

CONFIGNAME=${1:-config_bdt_nom.yaml}

echo "Using config: ${CONFIGNAME}"
echo "Starting all samples..."

root -l -b -q 'IsoROC_calculator.C("'"${CONFIGNAME}"'","photon5")'  &
root -l -b -q 'IsoROC_calculator.C("'"${CONFIGNAME}"'","photon10")' &
root -l -b -q 'IsoROC_calculator.C("'"${CONFIGNAME}"'","photon20")' &
root -l -b -q 'IsoROC_calculator.C("'"${CONFIGNAME}"'","jet5")'     &
root -l -b -q 'IsoROC_calculator.C("'"${CONFIGNAME}"'","jet12")'    &
root -l -b -q 'IsoROC_calculator.C("'"${CONFIGNAME}"'","jet20")'    &
root -l -b -q 'IsoROC_calculator.C("'"${CONFIGNAME}"'","jet30")'    &
root -l -b -q 'IsoROC_calculator.C("'"${CONFIGNAME}"'","jet40")'

wait
echo "All samples done."
