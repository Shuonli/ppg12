#!/usr/bin/env bash
# Runs RecoEffCalculator_TTreeReader for jet8 only, two passes.
# Used to recover after jet8/bdt_split.root was regenerated post-Stage-C.
set -e
source /sphenix/u/shuhang98/setup.sh
cd /sphenix/user/shuhangli/ppg12/efficiencytool
CONFIG=${1:-config_bdt_nosplit.yaml}
echo "[jet8 pass 1 vertex scan] config=$CONFIG"
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIG}"'", "jet8", true)'
echo "[jet8 pass 2 full analysis] config=$CONFIG"
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIG}"'", "jet8")'
echo "[jet8 done]"
