#!/bin/bash
# Full MC re-run after restoring the truth-ET, response-arm-only smearing (commit 40d1c02 scheme).
# Phase 1: per-period TTree pipeline for every feeder config (parents are skipped by the dispatcher).
# Phase 2: parents only (merge_periods + data hadd + CalculatePhotonYield), never the feeders again
#          (re-running feeder MergeSim in Phase 2 races the parent merges).
# Then: merge audit, iso_generator rebuild, unfold_iter post-processing, systematics aggregation.
set -eo pipefail
cd /sphenix/user/shuhangli/ppg12/efficiencytool
source /sphenix/u/shuhang98/setup.sh
set -u
LOG=logs/smear_revert_chain_$(date +%Y%m%d_%H%M).log
mkdir -p logs
exec > >(tee -a "${LOG}") 2>&1
poll_cluster() { local c="$1" s="${2:-120}"; sleep 10; while condor_q "${c}" -nobatch 2>/dev/null | grep -qE "^${c}\."; do sleep "${s}"; done; }
echo "[chain] $(date -Is) START (Phase 1: $(wc -l < syst_configs_rerun.list) configs)"
P1=$(condor_submit oneforall_tree_double_rerun.sub | grep -oE 'cluster [0-9]+' | grep -oE '[0-9]+' | tail -1)
echo "[chain] Phase 1 cluster ${P1}"; poll_cluster "${P1}" 120; echo "[chain] $(date -Is) Phase 1 cleared"
condor_history -constraint "ClusterId==${P1}" -af ExitCode 2>/dev/null | sort | uniq -c | sed 's/^/[chain] P1 exit codes: /'
P2=$(condor_submit oneforall_parents_rerun.sub | grep -oE 'cluster [0-9]+' | grep -oE '[0-9]+' | tail -1)
echo "[chain] Phase 2 (parents) cluster ${P2}"; poll_cluster "${P2}" 60; echo "[chain] $(date -Is) Phase 2 cleared"
condor_history -constraint "ClusterId==${P2}" -af ExitCode 2>/dev/null | sort | uniq -c | sed 's/^/[chain] P2 exit codes: /'
set +e
python3 audit_merges.py > logs/audit_merges_smear_revert.log 2>&1; echo "[chain] audit_merges.py exit=$?"
( cd ../plotting && root -l -b -q build_iso_generator_variant.C ) > logs/build_iso_generator_smear_revert.log 2>&1; echo "[chain] build_iso_generator_variant.C exit=$?"
python3 postprocess_unfold_iter_scan.py results > logs/postprocess_unfold_iter_smear_revert.log 2>&1; echo "[chain] postprocess_unfold_iter_scan.py exit=$?"
( cd ../plotting && python calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures --skip-missing ) > logs/calc_syst_bdt_smear_revert.log 2>&1; echo "[chain] calc_syst_bdt.py exit=$?"
set -e
echo "[chain] $(date -Is) DONE (P1=${P1} P2=${P2})"
