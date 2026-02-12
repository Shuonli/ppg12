#!/usr/bin/env bash

source /sphenix/u/shuhang98/setup.sh

# Check if a config file name was provided
if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file> [use_leading_tower_time] [apply_leading_time_corr]"
  echo "  use_leading_tower_time: true/false or 1/0 (default: true)"
  echo "  apply_leading_time_corr: true/false or 1/0 (default: false; only meaningful if use_leading_tower_time=true)"
  exit 1
fi

CONFIGNAME=$1
USE_LEADING=${2:-true}
APPLY_CORR=${3:-false}

# Normalize common values
case "${USE_LEADING}" in
  1|true|TRUE|True) USE_LEADING=true ;;
  0|false|FALSE|False) USE_LEADING=false ;;
  *)
    echo "ERROR: invalid use_leading_tower_time='${USE_LEADING}'. Use true/false or 1/0."
    exit 2
    ;;
esac

case "${APPLY_CORR}" in
  1|true|TRUE|True) APPLY_CORR=true ;;
  0|false|FALSE|False) APPLY_CORR=false ;;
  *)
    echo "ERROR: invalid apply_leading_time_corr='${APPLY_CORR}'. Use true/false or 1/0."
    exit 2
    ;;
esac

# Run plot_cluster_time for photon samples
root -l -b -q 'plot_cluster_time.C("'"${CONFIGNAME}"'", "photon5", '"${USE_LEADING}"', '"${APPLY_CORR}"')' &
root -l -b -q 'plot_cluster_time.C("'"${CONFIGNAME}"'", "photon10", '"${USE_LEADING}"', '"${APPLY_CORR}"')' &
root -l -b -q 'plot_cluster_time.C("'"${CONFIGNAME}"'", "photon20", '"${USE_LEADING}"', '"${APPLY_CORR}"')' &

# Run plot_cluster_time for jet samples
root -l -b -q 'plot_cluster_time.C("'"${CONFIGNAME}"'", "jet10", '"${USE_LEADING}"', '"${APPLY_CORR}"')' &
root -l -b -q 'plot_cluster_time.C("'"${CONFIGNAME}"'", "jet15", '"${USE_LEADING}"', '"${APPLY_CORR}"')' &
root -l -b -q 'plot_cluster_time.C("'"${CONFIGNAME}"'", "jet20", '"${USE_LEADING}"', '"${APPLY_CORR}"')' &
root -l -b -q 'plot_cluster_time.C("'"${CONFIGNAME}"'", "jet30", '"${USE_LEADING}"', '"${APPLY_CORR}"')' &
root -l -b -q 'plot_cluster_time.C("'"${CONFIGNAME}"'", "jet50", '"${USE_LEADING}"', '"${APPLY_CORR}"')'

# Wait for all background jobs to finish
wait

echo "All plot_cluster_time jobs completed"

# Run merge script (match the same cluster-time definition)
root -l -b -q 'merge_cluster_time.C('"${USE_LEADING}"', '"${APPLY_CORR}"')'

echo "Merging completed"
