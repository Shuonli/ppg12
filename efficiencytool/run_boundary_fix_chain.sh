#!/bin/bash
# Period-boundary fix (0 mrad run_max 51274 -> 51273): only the DATA histograms
# change, so re-run the data pass for the 0 mrad feeders, then the parents
# (data hadd + CalculatePhotonYield), then the post-steps and aggregation.
set -eo pipefail
cd /sphenix/user/shuhangli/ppg12/efficiencytool
source /sphenix/u/shuhang98/setup.sh
set -u
LOG=logs/boundary_fix_chain_$(date +%Y%m%d_%H%M).log
exec > >(tee -a "${LOG}") 2>&1
poll(){ local c="$1" s="${2:-60}"; sleep 10; while condor_q "${c}" -nobatch 2>/dev/null | grep -qE "^${c}\."; do sleep "${s}"; done; }
echo "[chain] $(date -Is) START data-only 0rad ($(wc -l < data0rad_rerun.list) configs)"
P1=$(condor_submit oneforall_data0rad.sub | grep -oE 'cluster [0-9]+' | grep -oE '[0-9]+' | tail -1)
echo "[chain] data cluster ${P1}"; poll "${P1}" 60; echo "[chain] $(date -Is) data pass cleared"
condor_history -constraint "ClusterId==${P1}" -af ExitCode 2>/dev/null | sort | uniq -c | sed 's/^/[chain] data exit codes: /'
P2=$(condor_submit oneforall_parents_rerun.sub | grep -oE 'cluster [0-9]+' | grep -oE '[0-9]+' | tail -1)
echo "[chain] parents cluster ${P2}"; poll "${P2}" 45; echo "[chain] $(date -Is) parents cleared"
condor_history -constraint "ClusterId==${P2}" -af ExitCode 2>/dev/null | sort | uniq -c | sed 's/^/[chain] parent exit codes: /'
set +e
python3 audit_merges.py > logs/audit_merges_boundary_fix.log 2>&1; echo "[chain] audit_merges.py exit=$?"
( cd ../plotting && root -l -b -q build_iso_generator_variant.C ) > logs/build_iso_generator_boundary_fix.log 2>&1; echo "[chain] build_iso_generator_variant.C exit=$?"
python3 postprocess_unfold_iter_scan.py results > logs/postprocess_unfold_iter_boundary_fix.log 2>&1; echo "[chain] postprocess_unfold_iter_scan.py exit=$?"
( cd ../plotting && python calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures --skip-missing ) > logs/calc_syst_bdt_boundary_fix.log 2>&1; echo "[chain] calc_syst_bdt.py exit=$?"
set -e
echo "[chain] $(date -Is) DONE (data=${P1} parents=${P2})"
