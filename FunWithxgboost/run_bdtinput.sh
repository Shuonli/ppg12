#!/usr/bin/env bash
# Condor executable for one BDTinput.C job.
# Args: $1 = YAML config, $2 = filetype (photon5, jet12, data, ...)
set -e
source /sphenix/u/shuhang98/setup.sh
WORKDIR="/sphenix/user/shuhangli/ppg12/FunWithxgboost"
cd "$WORKDIR"
echo "[run_bdtinput] config=$1 filetype=$2"
root -l -b -q 'BDTinput.C("'"$1"'", "'"$2"'")'
