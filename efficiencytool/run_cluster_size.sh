#!/bin/bash
# Run ClusterSizeStudy.C over data + nominal MC + DI MC samples.
#
# Two modes:
#   (default)          baseline (no cuts, no reweighting, all runs)
#                      → results/cluster_size_{sample}.root
#   PERIOD={1p5mrad|0mrad}
#                      per-period run with common/tight cuts applied,
#                      truth-vertex reweighting, run-range filter on data
#                      → results/cluster_size_cuts_{PERIOD}_{sample}.root
#
# Periods use the showershape configs for selection bounds:
#   1p5mrad → efficiencytool/config_showershape_1p5rad.yaml
#   0mrad   → efficiencytool/config_showershape_0rad.yaml
#
# Env vars:
#   PARALLEL=N    parallel sample jobs (default 4)
#   PERIOD=NAME   run in per-period cut mode; otherwise baseline mode

set -u

repo_root="/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12"
cd "${repo_root}/efficiencytool" || exit 1

mkdir -p results logs

SAMPLES=(
  data
  photon5 photon10 photon20
  photon10_nom photon10_double
  jet5 jet8 jet12 jet20 jet30 jet40
  jet12_nom jet12_double
)

PARALLEL="${PARALLEL:-4}"
PERIOD="${PERIOD:-}"

# Resolve config + output prefix based on period mode
CONFIGNAME=""
OUTPREFIX="results/cluster_size"
LOGPREFIX="logs/cluster_size"
if [ -n "${PERIOD}" ]; then
  case "${PERIOD}" in
    1p5mrad) CONFIGNAME="${repo_root}/efficiencytool/config_showershape_1p5rad.yaml" ;;
    0mrad)   CONFIGNAME="${repo_root}/efficiencytool/config_showershape_0rad.yaml"   ;;
    *)
      echo "unknown PERIOD=${PERIOD}  (must be 1p5mrad or 0mrad)"
      exit 2
      ;;
  esac
  if [ ! -f "${CONFIGNAME}" ]; then
    echo "config file missing: ${CONFIGNAME}"
    exit 2
  fi
  OUTPREFIX="results/cluster_size_cuts_${PERIOD}"
  LOGPREFIX="logs/cluster_size_cuts_${PERIOD}"
  echo "=== period=${PERIOD}  config=${CONFIGNAME} ==="
fi

run_one() {
  local sample="$1"
  local logfile="${LOGPREFIX}_${sample}.log"
  local outfile="${OUTPREFIX}_${sample}.root"
  echo "[$(date +%H:%M:%S)] starting ${sample}  -> ${outfile}"
  # Forward args: filetype, outpath_override, mix_weight, configname
  root -l -b -q "ClusterSizeStudy.C+(\"${sample}\",\"${outfile}\",1.0,\"${CONFIGNAME}\")" \
      &> "${logfile}"
  local rc=$?
  if [ ${rc} -eq 0 ] && [ -s "${outfile}" ]; then
    echo "[$(date +%H:%M:%S)] done     ${sample}  -> ${outfile}"
  else
    echo "[$(date +%H:%M:%S)] FAILED   ${sample}  (rc=${rc}, see ${logfile})"
  fi
}

export -f run_one
export LOGPREFIX OUTPREFIX CONFIGNAME

printf '%s\n' "${SAMPLES[@]}" | \
    xargs -n1 -P"${PARALLEL}" -I{} bash -c 'run_one "$@"' _ {}

echo ""
echo "Output summary:"
ls -lh "${OUTPREFIX}"_*.root 2>/dev/null | awk '{print $NF, $5}'
