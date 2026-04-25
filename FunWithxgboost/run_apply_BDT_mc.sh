#!/usr/bin/env bash
# Condor executable for one apply_BDT.C MC job.
# Arg: $1 = filetype. Valid values: photon{5,10,20}, photon{5,10,20}_double,
#      jet{5,8,12,20,30,40}, jet{8,12,20,30,40}_double, jet50_double.
#      See mc_filetypes.list for the full 18-sample list consumed by
#      submit_apply_BDT_mc.sub.
set -e
source /sphenix/u/shuhang98/setup.sh
WORKDIR="/sphenix/user/shuhangli/ppg12/FunWithxgboost"
cd "$WORKDIR"
echo "[run_apply_BDT_mc] filetype=$1"
root -l -b -q 'apply_BDT.C("config_nom.yaml", "'"$1"'")'
