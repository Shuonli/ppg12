#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/sphenix/user/shuhangli/ppg12/FunWithxgboost"
MACRO="$WORKDIR/apply_BDT.C"
CONFIG="config_nom.yaml"

filetypes=(
  #photon5
  #photon10
  #photon20
  #jet5
  #jet8
  #jet12
  #jet10
  #jet15
  #jet20
  #jet30
  #jet40
  photon10_double
  jet12_double
)

mkdir -p "$WORKDIR/logs"

pids=()
for ft in "${filetypes[@]}"; do
  logfile="$WORKDIR/logs/apply_BDT_${ft}.log"
  echo "Starting apply_BDT for $ft (logging to $logfile)"
  root -l -b -q "$MACRO(\"$CONFIG\",\"$ft\")" >"$logfile" 2>&1 &
  pids+=("$!")
done

exit_code=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    exit_code=1
  fi
done

exit "$exit_code"




