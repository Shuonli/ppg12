#!/usr/bin/env bash
# Auto-wrapper for oneforall_tree_double.sh:
#   1) Compute f_double from per-run luminosity data (calc_pileup_range.C)
#   2) Run the two-pass cross-section blending pipeline (oneforall_tree_double.sh)
#
# Usage:
#   bash oneforall_tree_double_auto.sh 0mrad
#   bash oneforall_tree_double_auto.sh 1.5mrad
#   bash oneforall_tree_double_auto.sh 0mrad config_bdt_nom.yaml   # explicit config
#
# The crossing-angle argument selects the run range and default config:
#   0mrad   → runs 47289–51274, config_bdt_nom.yaml
#   1.5mrad → runs 51274–54000, config_bdt_nom.yaml

source /sphenix/u/shuhang98/setup.sh

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPTDIR}"

ANGLE="${1:-}"
if [[ -z "${ANGLE}" ]]; then
    echo "Usage: bash oneforall_tree_double_auto.sh <0mrad|1.5mrad> [config.yaml]"
    exit 1
fi

case "${ANGLE}" in
    0mrad|0rad)
        RN_MIN=47289
        RN_MAX=51274
        DEFAULT_CONFIG="config_bdt_nom.yaml"
        ;;
    1.5mrad|1p5mrad)
        RN_MIN=51274
        RN_MAX=54000
        DEFAULT_CONFIG="config_bdt_nom.yaml"
        ;;
    *)
        echo "ERROR: unrecognized crossing angle '${ANGLE}'. Use 0mrad or 1.5mrad."
        exit 1
        ;;
esac

CONFIGNAME="${2:-${DEFAULT_CONFIG}}"

echo "============================================================"
echo "[auto] Crossing angle : ${ANGLE}"
echo "[auto] Run range      : ${RN_MIN} – ${RN_MAX}"
echo "[auto] Config         : ${CONFIGNAME}"
echo "============================================================"

# --- Step 1: compute f_double from per-run luminosity data ---
echo ""
echo "[auto] Step 1: Running calc_pileup_range.C(${RN_MIN}, ${RN_MAX}) ..."
echo ""

PILEUP_LOG=$(mktemp /tmp/pileup_XXXXXX.log)
root -l -b -q "calc_pileup_range.C(${RN_MIN}, ${RN_MAX})" 2>&1 | tee "${PILEUP_LOG}"

# Parse cluster-weighted double fraction from the output line:
#   cw_double = 0.224000
CW_DOUBLE=$(grep 'cw_double' "${PILEUP_LOG}" | grep -oE '[0-9]+\.[0-9]+$')

if [[ -z "${CW_DOUBLE}" ]]; then
    echo "ERROR: could not parse cw_double from calc_pileup_range.C output."
    echo "Check ${PILEUP_LOG} for details."
    exit 1
fi

echo ""
echo "============================================================"
echo "[auto] Extracted cw_double = ${CW_DOUBLE}  (cluster-weighted, triple+ folded in)"
echo "============================================================"
echo ""

rm -f "${PILEUP_LOG}"

# --- Step 2: run the two-pass cross-section blending pipeline ---
echo "[auto] Step 2: Running oneforall_tree_double.sh ${CONFIGNAME} ${CW_DOUBLE} ..."
echo ""

bash oneforall_tree_double.sh "${CONFIGNAME}" "${CW_DOUBLE}"

echo ""
echo "============================================================"
echo "[auto] Done. Crossing angle: ${ANGLE}, cw_double: ${CW_DOUBLE}"
echo "[auto] Config: ${CONFIGNAME}"
echo "============================================================"
