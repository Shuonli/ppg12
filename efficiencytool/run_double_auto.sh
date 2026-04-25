#!/usr/bin/env bash
# Unified auto-wrapper for double-interaction MC blending pipelines.
# Computes f_double from per-run luminosity data, then runs the appropriate
# two-pass blending pipeline (showershape or cross-section).
#
# Usage:
#   bash run_double_auto.sh --mode showershape 0mrad
#   bash run_double_auto.sh --mode showershape 1.5mrad
#   bash run_double_auto.sh --mode crosssection 0mrad
#   bash run_double_auto.sh --mode crosssection 1.5mrad config_bdt_custom.yaml
#
# The crossing-angle argument selects the run range and default config:
#   0mrad   -> runs 47289-51274
#   1.5mrad -> runs 51274-54000
#
# Default configs per mode:
#   showershape:  0mrad -> config_showershape_0rad.yaml
#                 1.5mrad -> config_showershape_1p5rad.yaml
#   crosssection: config_bdt_nom.yaml (both angles)

source /sphenix/u/shuhang98/setup.sh

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPTDIR}"

# --- Parse arguments ---
MODE=""
ANGLE=""
CONFIGNAME=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)
            MODE="$2"; shift 2 ;;
        0mrad|0rad|1.5mrad|1p5mrad)
            ANGLE="$1"; shift ;;
        *)
            CONFIGNAME="$1"; shift ;;
    esac
done

if [[ -z "${MODE}" || -z "${ANGLE}" ]]; then
    echo "Usage: bash run_double_auto.sh --mode <showershape|crosssection> <0mrad|1.5mrad> [config.yaml]"
    exit 1
fi

# --- Crossing-angle -> run range ---
case "${ANGLE}" in
    0mrad|0rad)
        RN_MIN=47289
        RN_MAX=51274
        ;;
    1.5mrad|1p5mrad)
        RN_MIN=51274
        RN_MAX=54000
        ;;
    *)
        echo "ERROR: unrecognized crossing angle '${ANGLE}'. Use 0mrad or 1.5mrad."
        exit 1
        ;;
esac

# --- Mode -> default config and inner script ---
case "${MODE}" in
    showershape)
        if [[ -z "${CONFIGNAME}" ]]; then
            case "${ANGLE}" in
                0mrad|0rad)     CONFIGNAME="config_showershape_0rad.yaml" ;;
                1.5mrad|1p5mrad) CONFIGNAME="config_showershape_1p5rad.yaml" ;;
            esac
        fi
        # NOTE: This legacy two-pass reco-vertex script was renamed to
        # run_showershape_double_reco_legacy.sh; the primary single-pass
        # truth-vertex showershape flow is now submit_showershape_di.sub.
        INNER_SCRIPT="run_showershape_double_reco_legacy.sh"
        ;;
    crosssection)
        CONFIGNAME="${CONFIGNAME:-config_bdt_nom.yaml}"
        INNER_SCRIPT="oneforall_tree_double.sh"
        ;;
    *)
        echo "ERROR: unrecognized mode '${MODE}'. Use showershape or crosssection."
        exit 1
        ;;
esac

echo "============================================================"
echo "[auto] Mode            : ${MODE}"
echo "[auto] Crossing angle  : ${ANGLE}"
echo "[auto] Run range       : ${RN_MIN} - ${RN_MAX}"
echo "[auto] Config          : ${CONFIGNAME}"
echo "[auto] Inner script    : ${INNER_SCRIPT}"
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

# --- Step 2: run the two-pass blending pipeline ---
echo "[auto] Step 2: Running ${INNER_SCRIPT} ${CONFIGNAME} ${CW_DOUBLE} ..."
echo ""

bash "${INNER_SCRIPT}" "${CONFIGNAME}" "${CW_DOUBLE}"

echo ""
echo "============================================================"
echo "[auto] Done. Mode: ${MODE}, Crossing angle: ${ANGLE}, cw_double: ${CW_DOUBLE}"
echo "[auto] Config: ${CONFIGNAME}"
echo "============================================================"
