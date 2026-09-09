#!/bin/bash
# Data-only RecoEff pass for one config (used when a change affects the data
# run-range filter but not the MC, e.g. the 2026-09-09 period-boundary fix).
set -eo pipefail
source /sphenix/u/shuhang98/setup.sh
cd /sphenix/user/shuhangli/ppg12/efficiencytool
CONFIG=${1:?config required}
echo "[data-only] ${CONFIG}"
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIG}"'", "data")'
