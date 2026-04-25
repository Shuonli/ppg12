#!/usr/bin/env bash
# Dispatch wrapper for make_selection_plots.sh — accepts a config filename
# and extracts the BDT suffix (bdt_<var_type>) for use by the plotting macros.
# Designed for condor submission (one job per config_bdt_*.yaml).

set -eo pipefail
source /sphenix/u/shuhang98/setup.sh
set -u

CONFIGNAME=${1:?config file required}
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPTDIR}"

# Strip path prefix and .yaml suffix; strip leading "config_" to get e.g. "bdt_nom".
basename=$(basename "${CONFIGNAME}" .yaml)
suffix="${basename#config_}"

echo "[make_selection_plots_dispatch] config: ${CONFIGNAME}"
echo "[make_selection_plots_dispatch] suffix: ${suffix}"

bash make_selection_plots.sh "${suffix}"
