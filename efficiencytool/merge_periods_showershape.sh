#!/usr/bin/env bash
#
# merge_periods_showershape.sh -- combine two per-crossing-angle showershape
# outputs (e.g. showershape_0rad + showershape_1p5mrad) into a single
# all-range showershape MC.
#
# Per-period MC is assumed to already be pre-scaled to the merge-target lumi
# via cross_weight *= lumi/lumi_target inside ShowerShapeCheck.C (single-pass
# truth-vertex reweight path). With the scaling baked into per-event Fill
# weights, plain hadd correctly combines TH1/TH2/TH3/TProfile.
#
# Usage:
#   bash merge_periods_showershape.sh [cfg_0rad] [cfg_1p5mrad] [cfg_combined]
# Defaults:
#   config_showershape_0rad.yaml + config_showershape_1p5mrad.yaml -> config_showershape.yaml
#
# Also optionally hadd's the per-period data histos as a cross-check against
# the dedicated all-range data ShowerShapeCheck job. The dedicated run is the
# authoritative all-range data (because it uses the all-range run_min/run_max
# directly); the hadd cross-check is written with the ".merged" suffix so it
# can be diffed manually.
set -e

CFG0=${1:-config_showershape_0rad.yaml}
CFG1=${2:-config_showershape_1p5mrad.yaml}
CFGO=${3:-config_showershape.yaml}

SUFFIX0=${CFG0#config_}; SUFFIX0=${SUFFIX0%.yaml}
SUFFIX1=${CFG1#config_}; SUFFIX1=${SUFFIX1%.yaml}
SUFFIXO=${CFGO#config_}; SUFFIXO=${SUFFIXO%.yaml}

RESULTS_DIR=${RESULTS_DIR:-results}

echo "[merge_periods_showershape] period0 suffix = ${SUFFIX0}"
echo "[merge_periods_showershape] period1 suffix = ${SUFFIX1}"
echo "[merge_periods_showershape] all-range suffix = ${SUFFIXO}"
echo "[merge_periods_showershape] results dir = ${RESULTS_DIR}"

SIGNAL_OUT="${RESULTS_DIR}/MC_efficiencyshower_shape_signal_combined_${SUFFIXO}.root"
JET_OUT="${RESULTS_DIR}/MC_efficiencyshower_shape_jet_inclusive_combined_${SUFFIXO}.root"
DATA_MERGED="${RESULTS_DIR}/data_histoshower_shape_${SUFFIXO}.merged.root"

echo
echo "--- signal_combined ---"
hadd -f "${SIGNAL_OUT}" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_signal_combined_${SUFFIX0}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_signal_combined_${SUFFIX1}.root"

echo
echo "--- jet_inclusive_combined ---"
hadd -f "${JET_OUT}" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet_inclusive_combined_${SUFFIX0}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet_inclusive_combined_${SUFFIX1}.root"

echo
echo "--- data hadd cross-check (.merged) ---"
hadd -f "${DATA_MERGED}" \
    "${RESULTS_DIR}/data_histoshower_shape_${SUFFIX0}.root" \
    "${RESULTS_DIR}/data_histoshower_shape_${SUFFIX1}.root" || echo "[warn] per-period data hadd skipped (files missing)"

echo
echo "--- signal SI-only (non-mixed reference) ---"
hadd -f "${RESULTS_DIR}/MC_efficiencyshower_shape_signal_${SUFFIXO}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_signal_${SUFFIX0}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_signal_${SUFFIX1}.root"

echo
echo "--- jet_inclusive SI-only (non-mixed reference) ---"
hadd -f "${RESULTS_DIR}/MC_efficiencyshower_shape_jet_inclusive_${SUFFIXO}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet_inclusive_${SUFFIX0}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet_inclusive_${SUFFIX1}.root"

echo
echo "[merge_periods_showershape] wrote ${SIGNAL_OUT}"
echo "[merge_periods_showershape] wrote ${JET_OUT}"
echo "[merge_periods_showershape] wrote ${DATA_MERGED} (cross-check; authoritative all-range data = data_histoshower_shape_${SUFFIXO}.root from dedicated run)"
echo "[merge_periods_showershape] wrote SI-only all-range: MC_efficiencyshower_shape_{signal,jet_inclusive}_${SUFFIXO}.root"
