#!/usr/bin/env bash
# LEGACY: two-pass reco-vertex showershape blending. The primary single-pass
# truth-vertex showershape pipeline is now submit_showershape_di.sub +
# showershape_di_jobs_{0rad,1p5rad}.list. Kept as a reference / cross-check.
# Renamed from run_showershape_double.sh.
#
# Usage: bash run_showershape_double_reco_legacy.sh [config_name] [double_frac]
#   config_name  — YAML config file name (default: config_showershape.yaml)
#   double_frac  — cluster-weighted double-interaction fraction (default: 0.224 = 22.4%, 0 mrad)
#                  Use 0.079 for 1.5 mrad beam runs (7.9%)
#                  Computed by calc_pileup_range.C: w(n)=n*f(n)/sum(k*f(k)), triple+ folded into double

source /sphenix/u/shuhang98/setup.sh
# Use config_showershape.yaml by default (optional first argument)
CONFIGNAME=${1:-config_showershape.yaml}
SUFFIX=${CONFIGNAME#config_}; SUFFIX=${SUFFIX%.yaml}

DOUBLE_FRAC=${2:-0.224}
SINGLE_FRAC=$(awk "BEGIN{printf \"%.6f\", 1.0 - ${DOUBLE_FRAC}}")
echo "[run_showershape_double] DOUBLE_FRAC=${DOUBLE_FRAC}  SINGLE_FRAC=${SINGLE_FRAC}"
RESULTS_DIR=results

# Pass 1: vertex scan for all samples (fast — fills h_vertexz only)
# MC samples are weighted by their physical interaction fractions so that
# hadd(nom_vtxscan, double_vtxscan) = blended MC vertex distribution.
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon10_double", false, true, '"${DOUBLE_FRAC}"')' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet12_double",    true,  true, '"${DOUBLE_FRAC}"')' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon10_nom",    false, true, '"${SINGLE_FRAC}"')' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet12_nom",       true,  true, '"${SINGLE_FRAC}"')' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "data",            false, true)'

wait

# Merge vtxscan files into combined per-sample-type vtxscan.
# Each input was filled with mix_weight, so hadd = 0.776*nom + 0.224*double.
PHOTON10_COMBINED_VTXSCAN="${RESULTS_DIR}/MC_efficiencyshower_shape_photon10_combined_${SUFFIX}_vtxscan.root"
JET12_COMBINED_VTXSCAN="${RESULTS_DIR}/MC_efficiencyshower_shape_jet12_combined_inclusive_${SUFFIX}_vtxscan.root"

hadd -f "${PHOTON10_COMBINED_VTXSCAN}" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_photon10_nom_${SUFFIX}_vtxscan.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_photon10_double_${SUFFIX}_vtxscan.root"

hadd -f "${JET12_COMBINED_VTXSCAN}" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet12_nom_inclusive_${SUFFIX}_vtxscan.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet12_double_inclusive_${SUFFIX}_vtxscan.root"

# Pass 2: full analysis for all samples.
# Both nom and double use the SAME combined vtxscan as sim input so that the
# vertex reweighting ratio is computed against the blended MC distribution.
# Each sample applies its own mix_weight to all histogram fills.
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon10_double", false, false, '"${DOUBLE_FRAC}"', "'"${PHOTON10_COMBINED_VTXSCAN}"'")' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet12_double",    true,  false, '"${DOUBLE_FRAC}"', "'"${JET12_COMBINED_VTXSCAN}"'")' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon10_nom",    false, false, '"${SINGLE_FRAC}"', "'"${PHOTON10_COMBINED_VTXSCAN}"'")' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet12_nom",       true,  false, '"${SINGLE_FRAC}"', "'"${JET12_COMBINED_VTXSCAN}"'")' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "data",            false)'

wait

# Merge final Pass 2 outputs.
# mix_weight was already embedded per-event, so plain hadd is correct.
hadd -f "${RESULTS_DIR}/MC_efficiencyshower_shape_photon10_combined_${SUFFIX}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_photon10_nom_${SUFFIX}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_photon10_double_${SUFFIX}.root"

hadd -f "${RESULTS_DIR}/MC_efficiencyshower_shape_jet12_combined_inclusive_${SUFFIX}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet12_nom_inclusive_${SUFFIX}.root" \
    "${RESULTS_DIR}/MC_efficiencyshower_shape_jet12_double_inclusive_${SUFFIX}.root"
