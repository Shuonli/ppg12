#!/usr/bin/env bash
# run_saturation_study.sh
# Runs SaturationStudy.C for all file types in parallel, then hadds the outputs.
#
# Usage: ./run_saturation_study.sh <config_file>
# Example: ./run_saturation_study.sh config_showershape.yaml

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <config_file>"
    exit 1
fi

CONFIGNAME=$1

# ---------------------------------------------------------------------------
# Extract output path info from the YAML config so we can build hadd targets.
# The values in the config look like:
#   eff_outfile:  "/path/to/results/MC_efficiency"
#   data_outfile: "/path/to/results/data_histo"
#   var_type:     "showershape"
# ---------------------------------------------------------------------------
strip_quotes() { echo "$1" | sed 's/["\x27 ]//g'; }

EFF_OUTFILE=$(strip_quotes "$(grep 'eff_outfile' "$CONFIGNAME"  | head -1 | cut -d: -f2-)")
DATA_OUTFILE=$(strip_quotes "$(grep 'data_outfile' "$CONFIGNAME" | head -1 | cut -d: -f2-)")
VAR_TYPE=$(strip_quotes "$(grep 'var_type' "$CONFIGNAME"          | head -1 | cut -d: -f2-)")

echo "eff_outfile  : ${EFF_OUTFILE}"
echo "data_outfile : ${DATA_OUTFILE}"
echo "var_type     : ${VAR_TYPE}"

# ---------------------------------------------------------------------------
# Per-filetype output files produced by SaturationStudy.C
# ---------------------------------------------------------------------------
OUT_PHOTON5="${EFF_OUTFILE}_saturation_photon5_${VAR_TYPE}.root"
OUT_PHOTON10="${EFF_OUTFILE}_saturation_photon10_${VAR_TYPE}.root"
OUT_PHOTON20="${EFF_OUTFILE}_saturation_photon20_${VAR_TYPE}.root"
OUT_JET10="${EFF_OUTFILE}_saturation_jet10_${VAR_TYPE}.root"
OUT_JET15="${EFF_OUTFILE}_saturation_jet15_${VAR_TYPE}.root"
OUT_JET20="${EFF_OUTFILE}_saturation_jet20_${VAR_TYPE}.root"
OUT_JET30="${EFF_OUTFILE}_saturation_jet30_${VAR_TYPE}.root"
OUT_JET50="${EFF_OUTFILE}_saturation_jet50_${VAR_TYPE}.root"
OUT_DATA="${DATA_OUTFILE}_saturation_${VAR_TYPE}.root"

# hadd targets
HADD_PHOTON="${EFF_OUTFILE}_saturation_photon_${VAR_TYPE}.root"
HADD_JET="${EFF_OUTFILE}_saturation_jet_${VAR_TYPE}.root"

# ---------------------------------------------------------------------------
# Run all file types in parallel
# ---------------------------------------------------------------------------
echo "========================================"
echo "Running SaturationStudy.C for all types"
echo "========================================"

root -l -b -q "SaturationStudy.C(\"${CONFIGNAME}\", \"photon5\")"  &
root -l -b -q "SaturationStudy.C(\"${CONFIGNAME}\", \"photon10\")" &
root -l -b -q "SaturationStudy.C(\"${CONFIGNAME}\", \"photon20\")" &
root -l -b -q "SaturationStudy.C(\"${CONFIGNAME}\", \"jet10\")"    &
root -l -b -q "SaturationStudy.C(\"${CONFIGNAME}\", \"jet15\")"    &
root -l -b -q "SaturationStudy.C(\"${CONFIGNAME}\", \"jet20\")"    &
root -l -b -q "SaturationStudy.C(\"${CONFIGNAME}\", \"jet30\")"    &
root -l -b -q "SaturationStudy.C(\"${CONFIGNAME}\", \"jet50\")"    &
root -l -b -q "SaturationStudy.C(\"${CONFIGNAME}\", \"data\")"     &

wait
echo "All file types done."

# ---------------------------------------------------------------------------
# hadd
# ---------------------------------------------------------------------------
echo "========================================"
echo "hadding photon samples -> ${HADD_PHOTON}"
echo "========================================"
hadd -f "${HADD_PHOTON}" "${OUT_PHOTON5}" "${OUT_PHOTON10}" "${OUT_PHOTON20}"

echo "========================================"
echo "hadding jet samples    -> ${HADD_JET}"
echo "========================================"
hadd -f "${HADD_JET}" "${OUT_JET10}" "${OUT_JET15}" "${OUT_JET20}" "${OUT_JET30}" "${OUT_JET50}"

echo "========================================"
echo "Done. Outputs:"
echo "  Signal MC : ${HADD_PHOTON}"
echo "  Jet MC    : ${HADD_JET}"
echo "  Data      : ${OUT_DATA}"
echo "========================================"
