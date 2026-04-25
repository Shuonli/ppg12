#!/usr/bin/env bash
# Condor executable for a single ShowerShapeCheck.C job in the DI (double-
# interaction) single-pass truth-vertex-reweight pipeline.
#
# Args:
#   $1 = config      (e.g. config_showershape_0rad.yaml)
#   $2 = sample      (e.g. photon10_nom, photon10_double, data)
#   $3 = doinclusive ("true" / "false")
#   $4 = mix_weight  (float; SINGLE_FRAC / DOUBLE_FRAC / 1.0 for data)
#
# The macro is called with do_vertex_scan=false and no vtxscan_sim_override:
# vertex reweighting is provided by TruthVertexReweightLoader.h reading the
# truth_vertex_reweight file referenced from the YAML (sample-invariant).
set -e
source /sphenix/u/shuhang98/setup.sh

WORKDIR="/sphenix/user/shuhangli/ppg12/efficiencytool"
cd "$WORKDIR"

CONFIG=${1:?missing config}
SAMPLE=${2:?missing sample}
DOINCLUSIVE=${3:?missing doinclusive}
MIX_WEIGHT=${4:?missing mix_weight}

echo "[run_showershape_di_job] config=${CONFIG} sample=${SAMPLE} doinclusive=${DOINCLUSIVE} mix_weight=${MIX_WEIGHT}"

root -l -b -q 'ShowerShapeCheck.C("'"${CONFIG}"'", "'"${SAMPLE}"'", '"${DOINCLUSIVE}"', false, '"${MIX_WEIGHT}"')'
