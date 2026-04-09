#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/sphenix/user/shuhangli/ppg12/FunWithxgboost"
MACRO="$WORKDIR/export_cluster_arrays_csv.C"
CONFIG="config_nom.yaml"
OUTPREFIX="clusters"

# How many ROOT jobs to run at once.
# Override like: MAX_JOBS=4 ./run_export_clusters_parallel.sh
MAX_JOBS="${MAX_JOBS:-8}"

# Optional: cap events per filetype (useful for quick tests). -1 means all events.
# Override like: MAX_EVENTS=10000 ./run_export_clusters_parallel.sh
MAX_EVENTS="${MAX_EVENTS:- -1}"

filetypes=(
  photon5
  photon10
  photon20
  jet10
  jet20
  jet15
  jet30
  jet50
)

mkdir -p "$WORKDIR/logs"

running=0
pids=()
for ft in "${filetypes[@]}"; do
  logfile="$WORKDIR/logs/export_clusters_${ft}.log"
  echo "Starting export_cluster_arrays_csv for $ft (logging to $logfile)"

  # Produces: ${OUTPREFIX}_${ft}.csv (e.g. clusters_photon5.csv)
  root -l -b -q "$MACRO(\"$CONFIG\",\"$ft\",\"$OUTPREFIX\",$MAX_EVENTS,49)" >"$logfile" 2>&1 &
  pids+=("$!")
  ((running+=1))

  # simple concurrency limit
  if (( running >= MAX_JOBS )); then
    if ! wait -n; then
      echo "A job failed; check logs under $WORKDIR/logs" >&2
      exit 1
    fi
    ((running-=1))
  fi
done

exit_code=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    exit_code=1
  fi
done

exit "$exit_code"


