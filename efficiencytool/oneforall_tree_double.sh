#!/usr/bin/env bash
# Two-pass double-interaction MC blending for the cross-section pipeline.
# Calls RecoEffCalculator_TTreeReader.C with mix_weight fractions, then
# runs MergeSim.C and CalculatePhotonYield.C downstream.
#
# Usage: bash oneforall_tree_double.sh <config_file> [double_frac]
#   config_file  — YAML config file name
#   double_frac  — fraction of double-interaction clusters (default: 0.187 = 18.7%)
#                  Use 0.072 for 1.5 mrad beam runs (7.2% double-interaction rate)

source /sphenix/u/shuhang98/setup.sh

if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file> [double_frac]"
  exit 1
fi

CONFIGNAME=$1
DOUBLE_FRAC=${2:-0.187}
SINGLE_FRAC=$(awk "BEGIN{printf \"%.6f\", 1.0 - ${DOUBLE_FRAC}}")

echo "============================================================"
echo "[oneforall_tree_double] Config       : ${CONFIGNAME}"
echo "[oneforall_tree_double] DOUBLE_FRAC  : ${DOUBLE_FRAC}"
echo "[oneforall_tree_double] SINGLE_FRAC  : ${SINGLE_FRAC}"
echo "============================================================"

# Extract output keys from the YAML config
VAR_TYPE=$(python3 -c "import yaml; c=yaml.safe_load(open('${CONFIGNAME}')); print(c['output']['var_type'])")
EFF_OUTFILE=$(python3 -c "import yaml; c=yaml.safe_load(open('${CONFIGNAME}')); print(c['output']['eff_outfile'])")
RESP_OUTFILE=$(python3 -c "import yaml; c=yaml.safe_load(open('${CONFIGNAME}')); print(c['output']['response_outfile'])")

echo "[oneforall_tree_double] VAR_TYPE     : ${VAR_TYPE}"
echo "[oneforall_tree_double] EFF_OUTFILE  : ${EFF_OUTFILE}"
echo "[oneforall_tree_double] RESP_OUTFILE : ${RESP_OUTFILE}"
echo "============================================================"

# ====================================================================
# Pass 1: vertex scan (fast — fills h_vertexz only)
# MC samples are weighted by their physical interaction fractions so
# that hadd(nom_vtxscan, double_vtxscan) = blended MC vertex distribution.
# ====================================================================
echo ""
echo "[oneforall_tree_double] === Pass 1: vertex scan ==="
echo ""

# Non-double samples (no mix_weight needed, default 1.0)
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon5",  true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon20", true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet5",     true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet8",     true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet20",    true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet30",    true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet40",    true)' &

# Double-paired samples with mix_weight
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon10_nom",    true, '"${SINGLE_FRAC}"')' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon10_double", true, '"${DOUBLE_FRAC}"')' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet12_nom",       true, '"${SINGLE_FRAC}"')' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet12_double",    true, '"${DOUBLE_FRAC}"')' &

# Data (no mix_weight)
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "data", true)'

wait

# ====================================================================
# hadd vtxscans: nom + double → combined per sample type
# Each input was filled with mix_weight, so hadd = SINGLE_FRAC*nom + DOUBLE_FRAC*double.
# ====================================================================
echo ""
echo "[oneforall_tree_double] === hadd vtxscans ==="
echo ""

PHOTON10_COMBINED_VTXSCAN="${EFF_OUTFILE}_photon10_combined_${VAR_TYPE}_vtxscan.root"
JET12_COMBINED_VTXSCAN="${EFF_OUTFILE}_jet12_combined_${VAR_TYPE}_vtxscan.root"

hadd -f "${PHOTON10_COMBINED_VTXSCAN}" \
    "${EFF_OUTFILE}_photon10_nom_${VAR_TYPE}_vtxscan.root" \
    "${EFF_OUTFILE}_photon10_double_${VAR_TYPE}_vtxscan.root"

hadd -f "${JET12_COMBINED_VTXSCAN}" \
    "${EFF_OUTFILE}_jet12_nom_${VAR_TYPE}_vtxscan.root" \
    "${EFF_OUTFILE}_jet12_double_${VAR_TYPE}_vtxscan.root"

# ====================================================================
# Pass 2: full analysis
# Both nom and double use the SAME combined vtxscan as sim input so that
# vertex reweighting is computed against the blended MC distribution.
# Each sample applies its own mix_weight to all histogram fills.
# ====================================================================
echo ""
echo "[oneforall_tree_double] === Pass 2: full analysis ==="
echo ""

# Non-double samples (default args — no mix_weight, no vtxscan override)
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon5")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon20")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet5")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet8")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet20")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet30")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet40")' &

# Double-paired samples with mix_weight and combined vtxscan override
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon10_nom",    false, '"${SINGLE_FRAC}"', "'"${PHOTON10_COMBINED_VTXSCAN}"'")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon10_double", false, '"${DOUBLE_FRAC}"', "'"${PHOTON10_COMBINED_VTXSCAN}"'")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet12_nom",       false, '"${SINGLE_FRAC}"', "'"${JET12_COMBINED_VTXSCAN}"'")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet12_double",    false, '"${DOUBLE_FRAC}"', "'"${JET12_COMBINED_VTXSCAN}"'")' &

# Data (default args)
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "data")'

wait

# ====================================================================
# hadd final outputs: merge nom + double into standard filenames
# that MergeSim.C expects (photon10, jet12 without _nom/_double suffix).
# mix_weight was already embedded per-event, so plain hadd is correct.
# ====================================================================
echo ""
echo "[oneforall_tree_double] === hadd final outputs ==="
echo ""

# Efficiency files
hadd -f "${EFF_OUTFILE}_photon10_${VAR_TYPE}.root" \
    "${EFF_OUTFILE}_photon10_nom_${VAR_TYPE}.root" \
    "${EFF_OUTFILE}_photon10_double_${VAR_TYPE}.root"

hadd -f "${EFF_OUTFILE}_jet12_${VAR_TYPE}.root" \
    "${EFF_OUTFILE}_jet12_nom_${VAR_TYPE}.root" \
    "${EFF_OUTFILE}_jet12_double_${VAR_TYPE}.root"

# Response matrix files
hadd -f "${RESP_OUTFILE}_photon10_${VAR_TYPE}.root" \
    "${RESP_OUTFILE}_photon10_nom_${VAR_TYPE}.root" \
    "${RESP_OUTFILE}_photon10_double_${VAR_TYPE}.root"

hadd -f "${RESP_OUTFILE}_jet12_${VAR_TYPE}.root" \
    "${RESP_OUTFILE}_jet12_nom_${VAR_TYPE}.root" \
    "${RESP_OUTFILE}_jet12_double_${VAR_TYPE}.root"

# ====================================================================
# Downstream: MergeSim + CalculatePhotonYield
# ====================================================================
echo ""
echo "[oneforall_tree_double] === Downstream: MergeSim + CalculatePhotonYield ==="
echo ""

root -l -b -q 'MergeSim.C("'"${CONFIGNAME}"'")'
root -l -b -q 'CalculatePhotonYield.C("'"${CONFIGNAME}"'", true)'   # isMC = true
root -l -b -q 'CalculatePhotonYield.C("'"${CONFIGNAME}"'", false)'  # isMC = false

echo ""
echo "============================================================"
echo "[oneforall_tree_double] Done. Config: ${CONFIGNAME}, DOUBLE_FRAC: ${DOUBLE_FRAC}"
echo "============================================================"
