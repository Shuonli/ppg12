#!/bin/bash
# Appendix cross-checks re-run with the current setup (2026-09-10): the tower-mask
# variants against the mask-off reference, and the back-to-back-jet NPB check.
# Phase 1 = per-period feeders, Phase 2 = all-range parents, then the merge audit.
set -eo pipefail
cd /sphenix/user/shuhangli/ppg12/efficiencytool
source /sphenix/u/shuhang98/setup.sh
set -u
LOG=logs/appendix_xcheck_chain_$(date +%Y%m%d_%H%M).log
exec > >(tee -a "${LOG}") 2>&1
poll(){ local c="$1" s="${2:-120}"; sleep 10; while condor_q "${c}" -nobatch 2>/dev/null | grep -qE "^${c}\."; do sleep "${s}"; done; }
echo "[chain] $(date -Is) START feeders ($(wc -l < xcheck_feeders_20260910.list) configs)"
P1=$(condor_submit oneforall_tree_double_xcheck.sub | grep -oE 'cluster [0-9]+' | grep -oE '[0-9]+' | tail -1)
echo "[chain] feeders cluster ${P1}"; poll "${P1}" 120; echo "[chain] $(date -Is) feeders cleared"
condor_history -constraint "ClusterId==${P1}" -af ExitCode 2>/dev/null | sort | uniq -c | sed 's/^/[chain] feeder exit codes: /'
P2=$(condor_submit oneforall_parents_xcheck.sub | grep -oE 'cluster [0-9]+' | grep -oE '[0-9]+' | tail -1)
echo "[chain] parents cluster ${P2}"; poll "${P2}" 60; echo "[chain] $(date -Is) parents cleared"
condor_history -constraint "ClusterId==${P2}" -af ExitCode 2>/dev/null | sort | uniq -c | sed 's/^/[chain] parent exit codes: /'
set +e
python3 audit_merges.py > logs/audit_merges_appendix_xcheck.log 2>&1; echo "[chain] audit_merges.py exit=$?"
for v in $(sed 's/^config_bdt_//; s/\.yaml$//' xcheck_parents_20260910.list); do printf "[chain] %-26s %s\n" "$v" "$(ls -l --time-style=+%m-%d_%H:%M results/Photon_final_bdt_${v}.root 2>/dev/null | awk '{print $6}')"; done
set -e
echo "[chain] $(date -Is) DONE (feeders=${P1} parents=${P2})"
