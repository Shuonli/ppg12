#!/bin/bash
# Appendix cross-checks re-run with the current setup (2026-09-10): the tower-mask
# variants against the mask-off reference, and the back-to-back-jet NPB check.
# Phase 1 = per-period feeders, Phase 2 = all-range parents, then the merge audit.
# Usage: bash run_appendix_xcheck_chain.sh [feeder_list] [parent_list] [tag]
# The 2026-09-11 blending-fraction re-run used difrac_feeders_20260911.list,
# difrac_parents_20260911.list and tag difrac_fit.
set -eo pipefail
cd /sphenix/user/shuhangli/ppg12/efficiencytool
source /sphenix/u/shuhang98/setup.sh
set -u
FEEDERS=${1:-xcheck_feeders_20260910.list}
PARENTS=${2:-xcheck_parents_20260910.list}
TAG=${3:-appendix_xcheck}
LOG=logs/${TAG}_chain_$(date +%Y%m%d_%H%M).log
exec > >(tee -a "${LOG}") 2>&1
poll(){ local c="$1" s="${2:-120}"; sleep 10; while condor_q "${c}" -nobatch 2>/dev/null | grep -qE "^${c}\."; do sleep "${s}"; done; }
echo "[chain] $(date -Is) START feeders ($(wc -l < "${FEEDERS}") configs from ${FEEDERS})"
P1=$(condor_submit list_file="${FEEDERS}" oneforall_tree_double_xcheck.sub | grep -oE 'cluster [0-9]+' | grep -oE '[0-9]+' | tail -1)
[ -n "${P1}" ] || { echo "[chain] feeder submission failed, stopping"; exit 1; }
echo "[chain] feeders cluster ${P1}"; poll "${P1}" 120; echo "[chain] $(date -Is) feeders cleared"
condor_history -constraint "ClusterId==${P1}" -af ExitCode 2>/dev/null | sort | uniq -c | sed 's/^/[chain] feeder exit codes: /'
BAD=$(condor_history -constraint "ClusterId==${P1} && ExitCode!=0" -af ClusterId 2>/dev/null | wc -l)
[ "${BAD}" -eq 0 ] || { echo "[chain] ${BAD} feeder job(s) failed, not running the parents"; exit 1; }
P2=$(condor_submit list_file="${PARENTS}" oneforall_parents_xcheck.sub | grep -oE 'cluster [0-9]+' | grep -oE '[0-9]+' | tail -1)
[ -n "${P2}" ] || { echo "[chain] parent submission failed, stopping"; exit 1; }
echo "[chain] parents cluster ${P2}"; poll "${P2}" 60; echo "[chain] $(date -Is) parents cleared"
condor_history -constraint "ClusterId==${P2}" -af ExitCode 2>/dev/null | sort | uniq -c | sed 's/^/[chain] parent exit codes: /'
set +e
python3 audit_merges.py > logs/audit_merges_${TAG}.log 2>&1; echo "[chain] audit_merges.py exit=$?"
for v in $(sed 's/^config_bdt_//; s/\.yaml$//' "${PARENTS}"); do printf "[chain] %-26s %s\n" "$v" "$(ls -l --time-style=+%m-%d_%H:%M results/Photon_final_bdt_${v}.root 2>/dev/null | awk '{print $6}')"; done
set -e
echo "[chain] $(date -Is) DONE (feeders=${P1} parents=${P2})"
