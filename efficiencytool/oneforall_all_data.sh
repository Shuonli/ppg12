#!/usr/bin/env bash
# Data-side regenerator for all-range configs.
#
# All-range configs (config_bdt_*_all.yaml) skip the per-sample TTree-reader
# stage in oneforall_tree_double.sh because their MC comes from a hadd of
# per-period merge-feeder MC. But the DATA side (filetype="data") still needs
# RecoEffCalculator_TTreeReader.C to be run on the all-range data file with
# the all-range config's cuts and run filter — otherwise data_histo_{var_type}.root
# is either missing (CalculatePhotonYield segfaults) or stale (wrong cuts).
#
# This script does that one-shot per config:
#   1. RecoEffCalculator_TTreeReader.C(config, "data")  -> data_histo
#   2. CalculatePhotonYield.C(config, false)            -> Photon_final (data)
#
# The MC closure (CalculatePhotonYield with isMC=true) was already produced
# by oneforall.sh in the merge_periods stage and is not redone here.

set -eo pipefail
source /sphenix/u/shuhang98/setup.sh
set -u

CONFIGNAME=${1:?config file required}
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPTDIR}"

echo "[oneforall_all_data] config: ${CONFIGNAME}"
echo "[oneforall_all_data] === Step 1: RecoEffCalculator data ==="
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "data")'

echo "[oneforall_all_data] === Step 2: CalculatePhotonYield (isMC=false) ==="
root -l -b -q 'CalculatePhotonYield.C("'"${CONFIGNAME}"'", false)'

echo "[oneforall_all_data] done: ${CONFIGNAME}"
