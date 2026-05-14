#!/usr/bin/env bash
# Variant pipeline chain (post 2026-05-13 BDT-cut fix, ET=42 crossover):
#   Phase 1 (oneforall_tree_double.sub, ~166 cfgs, 16 GB/job)
#   → wait
#   Phase 2 (oneforall.sub, ~166 cfgs, 2 GB/job)
#   → wait
#   → audit_merges.py
#
# Uses the patched polling pattern (condor_q <CLUSTER> -nobatch grep) — the
# earlier chain script that used `condor_q -submitter ${USER}` raced and
# exited immediately because the FQDN-submitter filter didn't match.
# `condor_history -constraint "ClusterId==..."` is wrapped in set +e so a
# history-query failure cannot abort the chain.

set -eo pipefail
cd /sphenix/user/shuhangli/ppg12/efficiencytool

poll_cluster() {
    local cluster="$1"
    local sleep_s="${2:-120}"
    # Give condor a few seconds to register the new submit before polling.
    sleep 5
    while condor_q "${cluster}" -nobatch 2>/dev/null | grep -qE "^${cluster}\."; do
        sleep "${sleep_s}"
    done
    set +e
    condor_history -constraint "ClusterId==${cluster}" 2>/dev/null | head -3
    set -e
}

# ---- Phase 1 ----
echo "[chain] $(date -Is) submitting Phase 1 (oneforall_tree_double.sub)..."
P1_OUT=$(condor_submit oneforall_tree_double.sub)
echo "${P1_OUT}"
P1_CLUSTER=$(echo "${P1_OUT}" | grep -oE "cluster [0-9]+" | grep -oE "[0-9]+" | tail -1)
echo "[chain] $(date -Is) Phase 1 cluster: ${P1_CLUSTER}"
echo "[chain] $(date -Is) polling for Phase 1 to clear (every 120s)..."
poll_cluster "${P1_CLUSTER}" 120
echo "[chain] $(date -Is) Phase 1 cleared."

# ---- Phase 2 ----
echo "[chain] $(date -Is) submitting Phase 2 (oneforall.sub)..."
P2_OUT=$(condor_submit oneforall.sub)
echo "${P2_OUT}"
P2_CLUSTER=$(echo "${P2_OUT}" | grep -oE "cluster [0-9]+" | grep -oE "[0-9]+" | tail -1)
echo "[chain] $(date -Is) Phase 2 cluster: ${P2_CLUSTER}"
echo "[chain] $(date -Is) polling for Phase 2 to clear (every 60s)..."
poll_cluster "${P2_CLUSTER}" 60
echo "[chain] $(date -Is) Phase 2 cleared."

# ---- Audit ----
echo "[chain] $(date -Is) running merge audit..."
set +e
python3 audit_merges.py > logs/audit_merges_post_cut_fix.log 2>&1
AUDIT_RC=$?
set -e
echo "[chain] $(date -Is) audit_merges.py exit=${AUDIT_RC}; log at logs/audit_merges_post_cut_fix.log"
echo "[chain] $(date -Is) variant pipeline + audit complete."
