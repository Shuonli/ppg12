#!/usr/bin/env bash
# Poll condor for cluster 1274917 (Phase 1 nominal TTree jobs) to finish,
# then submit oneforall_nom_only.sub (Phase 2). Recovery from the broken
# chain script that prematurely submitted Phase 2 due to a condor_q
# -submitter filter mismatch.

set -eo pipefail
cd /sphenix/user/shuhangli/ppg12/efficiencytool

CLUSTER=1274917
echo "[wait] $(date -Is) waiting for cluster ${CLUSTER} (Phase 1) to finish..."

# Poll with condor_q <cluster> -nobatch which does NOT require a submitter
# filter. Loop exits when no rows starting with "<cluster>." remain.
while condor_q "${CLUSTER}" -nobatch 2>/dev/null | grep -qE "^${CLUSTER}\."; do
    sleep 120
done

echo "[wait] $(date -Is) cluster ${CLUSTER} cleared from queue. Checking history for success..."
# Quick history sanity: should show 2 completed jobs (.0 and .1).
# condor_history takes -constraint, not a positional cluster ID, so use that
# (or just `condor_history -limit 5` for the most recent). Don't let a
# history-query failure abort the script — disable -e for this call only.
set +e
condor_history -constraint "ClusterId==${CLUSTER}" 2>/dev/null | head -5
set -e

echo "[wait] $(date -Is) submitting Phase 2 (oneforall_nom_only.sub)..."
condor_submit oneforall_nom_only.sub
echo "[wait] $(date -Is) Phase 2 submitted. Chain complete."
