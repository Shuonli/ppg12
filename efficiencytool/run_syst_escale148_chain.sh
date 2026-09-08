#!/usr/bin/env bash
# run_syst_escale148_chain.sh  (2026-06-08)
#
# Scoped systematic re-run on the new (waveform-sim) trees with the energy
# scale envelope changed to +/-1.48% (energyscale148up/down -> syst_type
# escale). Runs ONLY the syst_type variations + nominal listed in
# syst_configs.list (59 configs), skipping cross-check scans and nomold.
#
#   Phase 1  oneforall_tree_double_syst.sub  (16 GB/job, per-sample TTree)
#     -> wait
#   Phase 2  oneforall_syst.sub              (2 GB/job, merge + xsec)
#     -> wait
#   -> audit_merges.py
#   -> rebuild iso_generator variant (scales fresh nominal sigma by herwig factor)
#   -> calc_syst_bdt.py (aggregate escale-148 into the syst breakdown)
#
# Polling pattern: condor_q <CLUSTER> -nobatch | grep (the submitter-filter
# pattern raced and exited early).

set -eo pipefail
cd /sphenix/user/shuhangli/ppg12/efficiencytool
source /sphenix/u/shuhang98/setup.sh
set -u

LOG=logs/syst_escale148_chain.log
exec > >(tee -a "${LOG}") 2>&1

poll_cluster() {
    local cluster="$1"
    local sleep_s="${2:-120}"
    sleep 10
    while condor_q "${cluster}" -nobatch 2>/dev/null | grep -qE "^${cluster}\."; do
        sleep "${sleep_s}"
    done
    set +e
    condor_history -constraint "ClusterId==${cluster}" 2>/dev/null | head -3
    set -e
}

echo "=========================================================="
echo "[chain] $(date -Is) START syst escale148 re-run"
echo "[chain] configs: $(wc -l < syst_configs.list) (syst_configs.list)"
echo "=========================================================="

# ---- Phase 1 ----
echo "[chain] $(date -Is) submitting Phase 1 (oneforall_tree_double_syst.sub)..."
P1_OUT=$(condor_submit oneforall_tree_double_syst.sub)
echo "${P1_OUT}"
P1_CLUSTER=$(echo "${P1_OUT}" | grep -oE "cluster [0-9]+" | grep -oE "[0-9]+" | tail -1)
echo "[chain] $(date -Is) Phase 1 cluster: ${P1_CLUSTER}; polling every 120s..."
poll_cluster "${P1_CLUSTER}" 120
echo "[chain] $(date -Is) Phase 1 cleared."

# ---- Phase 2 ----
echo "[chain] $(date -Is) submitting Phase 2 (oneforall_syst.sub)..."
P2_OUT=$(condor_submit oneforall_syst.sub)
echo "${P2_OUT}"
P2_CLUSTER=$(echo "${P2_OUT}" | grep -oE "cluster [0-9]+" | grep -oE "[0-9]+" | tail -1)
echo "[chain] $(date -Is) Phase 2 cluster: ${P2_CLUSTER}; polling every 60s..."
poll_cluster "${P2_CLUSTER}" 60
echo "[chain] $(date -Is) Phase 2 cleared."

# ---- Audit ----
echo "[chain] $(date -Is) running merge audit..."
set +e
python3 audit_merges.py > logs/audit_merges_escale148.log 2>&1
echo "[chain] audit_merges.py exit=$? (log: logs/audit_merges_escale148.log)"

# ---- iso_generator variant (post-processed: scales fresh nominal sigma) ----
echo "[chain] $(date -Is) rebuilding iso_generator variant..."
( cd ../plotting && root -l -b -q build_iso_generator_variant.C ) \
    > logs/build_iso_generator_escale148.log 2>&1
echo "[chain] build_iso_generator_variant.C exit=$? (log: logs/build_iso_generator_escale148.log)"

# ---- Aggregate systematics ----
echo "[chain] $(date -Is) aggregating systematics (calc_syst_bdt.py)..."
( cd ../plotting && python calc_syst_bdt.py --results ../efficiencytool/results \
    --outdir rootFiles --figdir figures --skip-missing ) \
    > logs/calc_syst_bdt_escale148.log 2>&1
echo "[chain] calc_syst_bdt.py exit=$? (log: logs/calc_syst_bdt_escale148.log)"
set -e

echo "=========================================================="
echo "[chain] $(date -Is) DONE syst escale148 re-run (P1=${P1_CLUSTER} P2=${P2_CLUSTER})"
echo "=========================================================="
