#!/usr/bin/env bash
# Submit nominal Phase 1 (per-period TTree jobs), wait for them to clear
# the condor queue, then submit Phase 2 (merge_periods + RecoEffCalculator
# on data + CalculatePhotonYield for the merged MC). Logs go to logs/.
#
# Triggered after the BDT-cut overlap fix (ET=42 crossover, intercept
# 0.815625 / 0.684375, slope +-0.0015625). See config_bdt_nom*.yaml.

set -eo pipefail
cd /sphenix/user/shuhangli/ppg12/efficiencytool

echo "[chain] $(date -Is) Submitting Phase 1 (per-period TTree jobs)..."
P1_OUT=$(condor_submit oneforall_tree_double_nom_only.sub)
echo "${P1_OUT}"
CLUSTER=$(echo "${P1_OUT}" | grep -oE "cluster [0-9]+" | grep -oE "[0-9]+" | tail -1)
echo "[chain] $(date -Is) Phase 1 cluster: ${CLUSTER}"

# Poll condor_q every 2 minutes for the cluster's jobs to clear.
echo "[chain] $(date -Is) Waiting for Phase 1 jobs (cluster ${CLUSTER}) to finish..."
while condor_q -submitter "${USER}" 2>/dev/null | awk -v c="${CLUSTER}" '$1 ~ "^"c"\\." {found=1} END{exit !found}'; do
    sleep 120
done

echo "[chain] $(date -Is) Phase 1 done. Submitting Phase 2 (merge_periods + final yield)..."
condor_submit oneforall_nom_only.sub
echo "[chain] $(date -Is) Phase 2 submitted. Pipeline chain complete."
