#!/usr/bin/env bash
# Two-pass double-interaction MC blending for the cross-section pipeline.
# Calls RecoEffCalculator_TTreeReader.C with mix_weight fractions on all 8
# single-interaction (SI) / double-interaction (DI) MC pairs, merges them into
# the single per-sample outputs that MergeSim.C expects, and finally runs
# MergeSim.C and CalculatePhotonYield.C downstream.
#
# Pass 1 (vertex scan) is automatically skipped when the config sets
# analysis.truth_vertex_reweight_on=1 — under truth-vertex reweighting the
# reco-vertex reweight branch is force-disabled in RecoEffCalculator and the
# Pass-1 vtxscan output is never read, so the ~15-20 min Pass 1 compute is
# wasted. For legacy configs with truth_vertex_reweight_on=0, both passes
# run and Pass 2 reads the hadded vtxscan as before.
#
# Sample coverage (8 pairs, 16 MC jobs + 1 data job per pass):
#   signal : photon5, photon10, photon20
#   bg     : jet8, jet12, jet20, jet30, jet40
#   skipped: jet5 (no DI partner), jet50_double (no SI partner)
#
# Usage: bash oneforall_tree_double.sh <config_file> [double_frac]
#   config_file  — YAML config file name
#   double_frac  — cluster-weighted double-interaction fraction (default: 0.224 = 22.4%, 0 mrad)
#                  Use 0.079 for 1.5 mrad beam runs (7.9%)
#                  Computed by calc_pileup_range.C: w(n)=n*f(n)/sum(k*f(k)), triple+ folded into double

source /sphenix/u/shuhang98/setup.sh

if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file> [double_frac]"
  exit 1
fi

CONFIGNAME=$1
DOUBLE_FRAC=${2:-0.224}
SINGLE_FRAC=$(awk "BEGIN{printf \"%.6f\", 1.0 - ${DOUBLE_FRAC}}")

# Single source of truth for the SI/DI pair list. Every entry here participates
# in Pass 1 (vtxscan), Pass 2 (full analysis), and the final hadd that produces
# the single per-sample file that MergeSim.C consumes.
SAMPLE_PAIRS=(photon5 photon10 photon20 jet8 jet12 jet20 jet30 jet40)

echo "============================================================"
echo "[oneforall_tree_double] Config       : ${CONFIGNAME}"
echo "[oneforall_tree_double] DOUBLE_FRAC  : ${DOUBLE_FRAC}"
echo "[oneforall_tree_double] SINGLE_FRAC  : ${SINGLE_FRAC}"
echo "[oneforall_tree_double] Sample pairs : ${SAMPLE_PAIRS[*]}"
echo "============================================================"

# Extract output keys from the YAML config
VAR_TYPE=$(python3 -c "import yaml; c=yaml.safe_load(open('${CONFIGNAME}')); print(c['output']['var_type'])")
EFF_OUTFILE=$(python3 -c "import yaml; c=yaml.safe_load(open('${CONFIGNAME}')); print(c['output']['eff_outfile'])")
RESP_OUTFILE=$(python3 -c "import yaml; c=yaml.safe_load(open('${CONFIGNAME}')); print(c['output']['response_outfile'])")
TRUTH_VTX_REWEIGHT_ON=$(python3 -c "import yaml; c=yaml.safe_load(open('${CONFIGNAME}')); print(c['analysis'].get('truth_vertex_reweight_on', 0))")

echo "[oneforall_tree_double] VAR_TYPE                : ${VAR_TYPE}"
echo "[oneforall_tree_double] EFF_OUTFILE             : ${EFF_OUTFILE}"
echo "[oneforall_tree_double] RESP_OUTFILE            : ${RESP_OUTFILE}"
echo "[oneforall_tree_double] truth_vertex_reweight_on: ${TRUTH_VTX_REWEIGHT_ON}"
echo "============================================================"

# Pass 1 is wasted work under truth-vertex reweighting — skip it.
if [[ "${TRUTH_VTX_REWEIGHT_ON}" == "1" ]]; then
    SKIP_PASS1=1
    echo "[oneforall_tree_double] truth_vertex_reweight_on=1 → Pass 1 (vertex scan) SKIPPED"
    echo "[oneforall_tree_double] Pass 2 will pass vtxscan_sim_override=\"\" (reco-vertex reweight is force-disabled in RecoEff anyway)"
else
    SKIP_PASS1=0
    echo "[oneforall_tree_double] truth_vertex_reweight_on=0 → running both Pass 1 and Pass 2 (legacy reco-vertex reweight)"
fi
echo "============================================================"

declare -A COMBINED_VTXSCAN

if [[ "${SKIP_PASS1}" == "0" ]]; then
    # ================================================================
    # Pass 1: vertex scan (fast — fills h_vertexz only).
    # MC samples are weighted by their physical interaction fractions so
    # that hadd(nom_vtxscan, double_vtxscan) = blended MC vertex shape.
    # ================================================================
    echo ""
    echo "[oneforall_tree_double] === Pass 1: vertex scan ==="
    echo ""

    for SAMPLE in "${SAMPLE_PAIRS[@]}"; do
        root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "'"${SAMPLE}"'_nom",    true, '"${SINGLE_FRAC}"')' &
        root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "'"${SAMPLE}"'_double", true, '"${DOUBLE_FRAC}"')' &
    done

    # Data (no mix_weight)
    root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "data", true)'

    wait

    # ================================================================
    # hadd vtxscans: nom + double → combined per sample type.
    # Each input was filled with mix_weight, so hadd = SINGLE_FRAC*nom +
    # DOUBLE_FRAC*double.
    # ================================================================
    echo ""
    echo "[oneforall_tree_double] === hadd vtxscans ==="
    echo ""

    for SAMPLE in "${SAMPLE_PAIRS[@]}"; do
        COMBINED_VTXSCAN[${SAMPLE}]="${EFF_OUTFILE}_${SAMPLE}_combined_${VAR_TYPE}_vtxscan.root"
        hadd -f "${COMBINED_VTXSCAN[${SAMPLE}]}" \
            "${EFF_OUTFILE}_${SAMPLE}_nom_${VAR_TYPE}_vtxscan.root" \
            "${EFF_OUTFILE}_${SAMPLE}_double_${VAR_TYPE}_vtxscan.root"
    done
else
    # Empty vtxscan override → RecoEff falls back to vertex_reweight_file in
    # the config (which under truth_vertex_reweight_on=1 is itself bypassed
    # by the truth-vertex weight branch in the event loop).
    for SAMPLE in "${SAMPLE_PAIRS[@]}"; do
        COMBINED_VTXSCAN[${SAMPLE}]=""
    done
fi

# ====================================================================
# Pass 2: full analysis. Each sample applies its own mix_weight to all
# histogram fills. When Pass 1 ran, vtxscan_sim_override is the hadded
# blended vertex shape; otherwise it's empty.
# ====================================================================
echo ""
echo "[oneforall_tree_double] === Pass 2: full analysis ==="
echo ""

for SAMPLE in "${SAMPLE_PAIRS[@]}"; do
    VTXSCAN="${COMBINED_VTXSCAN[${SAMPLE}]}"
    root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "'"${SAMPLE}"'_nom",    false, '"${SINGLE_FRAC}"', "'"${VTXSCAN}"'")' &
    root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "'"${SAMPLE}"'_double", false, '"${DOUBLE_FRAC}"', "'"${VTXSCAN}"'")' &
done

# Data (default args)
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "data")'

wait

# ====================================================================
# hadd final outputs: merge nom + double into the single per-sample
# filename that MergeSim.C expects (e.g. ${EFF_OUTFILE}_photon10_${VAR_TYPE}.root).
# mix_weight was already embedded per-event, so plain hadd is correct.
# ====================================================================
echo ""
echo "[oneforall_tree_double] === hadd final outputs ==="
echo ""

for SAMPLE in "${SAMPLE_PAIRS[@]}"; do
    # Efficiency file consumed by MergeSim.C
    hadd -f "${EFF_OUTFILE}_${SAMPLE}_${VAR_TYPE}.root" \
        "${EFF_OUTFILE}_${SAMPLE}_nom_${VAR_TYPE}.root" \
        "${EFF_OUTFILE}_${SAMPLE}_double_${VAR_TYPE}.root"

    # Response-matrix file (signal samples only: photon5/10/20)
    if [[ "${SAMPLE}" == photon* ]]; then
        hadd -f "${RESP_OUTFILE}_${SAMPLE}_${VAR_TYPE}.root" \
            "${RESP_OUTFILE}_${SAMPLE}_nom_${VAR_TYPE}.root" \
            "${RESP_OUTFILE}_${SAMPLE}_double_${VAR_TYPE}.root"
    fi
done

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
