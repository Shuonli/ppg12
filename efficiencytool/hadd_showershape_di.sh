#!/usr/bin/env bash
# Downstream hadd for the single-pass DI showershape pipeline.
# After submit_showershape_di.sub has produced all per-sample outputs, this
# script merges them into the two combined files consumed downstream.
#
# Per-event cross-section weights are already applied inside ShowerShapeCheck.C
# (via CrossSectionWeights.h) and multiplied by mix_weight = SINGLE_FRAC or
# DOUBLE_FRAC. Therefore plain hadd produces the correctly blended
# signal / jet_inclusive combined histograms.
#
# Mirrors the final hadd block of the retired two-pass reco-vertex pipeline
# (run_showershape_double_reco_legacy.sh lines 54-62), extended from the
# 2-sample (photon10, jet12) scope to the full 8-pair scope.
#
# Usage:
#   bash hadd_showershape_di.sh [config_showershape_0rad.yaml]
set -e

CONFIGNAME=${1:-config_showershape_0rad.yaml}
SUFFIX=${CONFIGNAME#config_}; SUFFIX=${SUFFIX%.yaml}
RESULTS_DIR=${RESULTS_DIR:-results}

echo "[hadd_showershape_di] config=${CONFIGNAME} suffix=${SUFFIX} results=${RESULTS_DIR}"

SIGNAL_OUT="${RESULTS_DIR}/MC_efficiencyshower_shape_signal_combined_${SUFFIX}.root"
JET_OUT="${RESULTS_DIR}/MC_efficiencyshower_shape_jet_inclusive_combined_${SUFFIX}.root"

# Signal = 3 photon pT windows x {nom, double} = 6 files.
hadd -f "${SIGNAL_OUT}" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_photon5_nom_${SUFFIX}.root"     \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_photon5_double_${SUFFIX}.root"  \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_photon10_nom_${SUFFIX}.root"    \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_photon10_double_${SUFFIX}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_photon20_nom_${SUFFIX}.root"    \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_photon20_double_${SUFFIX}.root"

# Jet inclusive = 5 jet pT windows x {nom, double} = 10 files.
hadd -f "${JET_OUT}" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet8_nom_inclusive_${SUFFIX}.root"     \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet8_double_inclusive_${SUFFIX}.root"  \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet12_nom_inclusive_${SUFFIX}.root"    \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet12_double_inclusive_${SUFFIX}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet20_nom_inclusive_${SUFFIX}.root"    \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet20_double_inclusive_${SUFFIX}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet30_nom_inclusive_${SUFFIX}.root"    \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet30_double_inclusive_${SUFFIX}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet40_nom_inclusive_${SUFFIX}.root"    \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet40_double_inclusive_${SUFFIX}.root"

echo "[hadd_showershape_di] wrote ${SIGNAL_OUT}"
echo "[hadd_showershape_di] wrote ${JET_OUT}"

# Background-only jet hadd: jet samples filled with doinclusive=false (truth-
# matched signal photons removed via the background_set veto in
# ShowerShapeCheck.C). Inputs are the per-sample files WITHOUT the _inclusive
# infix (produced by the showershape_di_jobs_*_bkgonly.list condor jobs).
# Used to overlay a "Background only" curve on Fig. 1 and Fig. 2 of the paper.
JET_BKGONLY="${RESULTS_DIR}/MC_efficiencyshower_shape_jet_background_only_combined_${SUFFIX}.root"
if ls "${RESULTS_DIR}/MC_efficiencyshower_shape_jet8_nom_${SUFFIX}.root" >/dev/null 2>&1; then
    hadd -f "${JET_BKGONLY}" \
        "${RESULTS_DIR}/MC_efficiencyshower_shape_jet8_nom_${SUFFIX}.root"     \
        "${RESULTS_DIR}/MC_efficiencyshower_shape_jet8_double_${SUFFIX}.root"  \
        "${RESULTS_DIR}/MC_efficiencyshower_shape_jet12_nom_${SUFFIX}.root"    \
        "${RESULTS_DIR}/MC_efficiencyshower_shape_jet12_double_${SUFFIX}.root" \
        "${RESULTS_DIR}/MC_efficiencyshower_shape_jet20_nom_${SUFFIX}.root"    \
        "${RESULTS_DIR}/MC_efficiencyshower_shape_jet20_double_${SUFFIX}.root" \
        "${RESULTS_DIR}/MC_efficiencyshower_shape_jet30_nom_${SUFFIX}.root"    \
        "${RESULTS_DIR}/MC_efficiencyshower_shape_jet30_double_${SUFFIX}.root" \
        "${RESULTS_DIR}/MC_efficiencyshower_shape_jet40_nom_${SUFFIX}.root"    \
        "${RESULTS_DIR}/MC_efficiencyshower_shape_jet40_double_${SUFFIX}.root"
    echo "[hadd_showershape_di] wrote ${JET_BKGONLY}"
else
    echo "[hadd_showershape_di] skipping background-only hadd (per-sample inputs not found; run showershape_di_jobs_${SUFFIX#showershape_}_bkgonly.list first)"
fi

# SI-only reference hadds (no DI blending) — consumed by
# plot_showershapes_variations.C when use_mixed=false (the `dis_*` prefix).
# Each _nom sample still carries mix_weight = SINGLE_FRAC, so per-bin shape
# normalization is unchanged from a true SI-only; only the overall scale
# differs by a constant, which cancels in the unit-area normalization used
# for shower-shape comparisons.
SIGNAL_SI="${RESULTS_DIR}/MC_efficiencyshower_shape_signal_${SUFFIX}.root"
JET_SI="${RESULTS_DIR}/MC_efficiencyshower_shape_jet_inclusive_${SUFFIX}.root"

hadd -f "${SIGNAL_SI}" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_photon5_nom_${SUFFIX}.root"  \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_photon10_nom_${SUFFIX}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_photon20_nom_${SUFFIX}.root"

hadd -f "${JET_SI}" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet8_nom_inclusive_${SUFFIX}.root"  \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet12_nom_inclusive_${SUFFIX}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet20_nom_inclusive_${SUFFIX}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet30_nom_inclusive_${SUFFIX}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet40_nom_inclusive_${SUFFIX}.root"

echo "[hadd_showershape_di] wrote ${SIGNAL_SI} (SI-only reference)"
echo "[hadd_showershape_di] wrote ${JET_SI} (SI-only reference)"
