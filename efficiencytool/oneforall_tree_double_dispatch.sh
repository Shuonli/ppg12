#!/usr/bin/env bash
# oneforall_tree_double_dispatch.sh
#
# Per-config dispatch for the DI cross-section pipeline. Derives DOUBLE_FRAC
# from the config's analysis.run_min (no calc_pileup_range.C call — saves
# ~30s per config), then invokes oneforall_tree_double.sh which does the
# full per-sample TTree-reader → MergeSim → CalculatePhotonYield pipeline.
#
# DOUBLE_FRAC convention (from .claude/rules/double-interaction.md, cluster
# weighted, triple+ folded into double):
#   0 mrad   (run_min < 51274): 0.224
#   1.5 mrad (run_min >= 51274): 0.079
#
# Dispatch convention (post Plan B restructure):
#   config_bdt_X_0rad.yaml    → 0mrad merge-feeder → run DI TTree pipeline
#   config_bdt_X_1p5mrad.yaml → 1.5mrad merge-feeder → run DI TTree pipeline
#   config_bdt_X.yaml (bare)  → all-range parent → SKIP (oneforall.sh handles
#                               it via merge_periods.sh after Phase 1 finishes)
#   config_bdt_X_0mrad.yaml   → legacy ntbdtpair scan per-period (note the
#                               different "0mrad" vs "0rad" suffix) → run
#                               DI TTree pipeline
#
# A bare name may also be a true standalone (no companion merge-feeders);
# in that case we run the TTree pipeline for it too.

set -eo pipefail
source /sphenix/u/shuhang98/setup.sh
set -u

CONFIGNAME=${1:?config file required}
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPTDIR}"

# Skip bare all-range parent configs whose merge-feeder siblings exist.
case "${CONFIGNAME}" in
  *_0rad.yaml|*_1p5mrad.yaml|*_0mrad.yaml)
    : # per-period config — proceed with the TTree pipeline
    ;;
  *)
    base=${CONFIGNAME%.yaml}
    if [ -f "${base}_0rad.yaml" ] && [ -f "${base}_1p5mrad.yaml" ]; then
      echo "[dispatch] ${CONFIGNAME} is an all-range parent (merge-feeder siblings exist); skip TTree pipeline"
      echo "[dispatch] (handle via oneforall.sh after merge-feeders finish)"
      exit 0
    fi
    # Standalone (no per-period siblings): proceed with the TTree pipeline.
    ;;
esac

RUN_MIN=$(python3 -c "import yaml; print(yaml.safe_load(open('${CONFIGNAME}'))['analysis']['run_min'])")
# T3 (di_fraction systematic): allow per-config override of DOUBLE_FRAC via
# analysis.double_frac_override. Falls through to run_min-based default when
# the field is absent or null. ruamel emits "None" for null; we treat any
# non-numeric value as absent.
DOUBLE_FRAC_OVERRIDE=$(python3 -c "
import yaml
y = yaml.safe_load(open('${CONFIGNAME}'))
v = y.get('analysis', {}).get('double_frac_override', None)
print('' if v is None else v)
" 2>/dev/null)

if [[ -n "${DOUBLE_FRAC_OVERRIDE}" ]]; then
  DOUBLE_FRAC="${DOUBLE_FRAC_OVERRIDE}"
  CROSSING="override (run_min=${RUN_MIN})"
elif [[ "${RUN_MIN}" -ge 0 && "${RUN_MIN}" -lt 51274 ]]; then
  DOUBLE_FRAC=0.224
  CROSSING="0mrad"
elif [[ "${RUN_MIN}" -ge 51274 ]]; then
  DOUBLE_FRAC=0.079
  CROSSING="1.5mrad"
else
  # run_min == -1 (no filter): treat as 1.5 mrad nominal default
  DOUBLE_FRAC=0.079
  CROSSING="1.5mrad (default, run_min=${RUN_MIN})"
fi

echo "[dispatch] config       : ${CONFIGNAME}"
echo "[dispatch] run_min      : ${RUN_MIN}"
echo "[dispatch] crossing     : ${CROSSING}"
echo "[dispatch] DOUBLE_FRAC  : ${DOUBLE_FRAC}"

bash oneforall_tree_double.sh "${CONFIGNAME}" "${DOUBLE_FRAC}"
