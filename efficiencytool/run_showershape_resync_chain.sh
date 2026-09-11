#!/bin/bash
# Shower-shape DI pipeline re-run after config_showershape_*.yaml were synced to the
# nominal selection on 2026-09-10 (BDT model base_v3E, tight/non-tight thresholds,
# isolation cut, MC iso shift, tower mask). Every sample, both periods, then the
# per-period and all-range merges and the note figures built from them.
set -eo pipefail
cd /sphenix/user/shuhangli/ppg12/efficiencytool
source /sphenix/u/shuhang98/setup.sh
set -u
LOG=logs/showershape_resync_chain_$(date +%Y%m%d_%H%M).log
exec > >(tee -a "${LOG}") 2>&1
poll(){ local c="$1" s="${2:-120}"; sleep 10; while condor_q "${c}" -nobatch 2>/dev/null | grep -qE "^${c}\."; do sleep "${s}"; done; }
echo "[chain] $(date -Is) START showershape jobs ($(wc -l < showershape_di_jobs_resync_20260910.list) rows)"
P=$(condor_submit list_file=showershape_di_jobs_resync_20260910.list submit_showershape_di.sub | grep -oE 'cluster [0-9]+' | grep -oE '[0-9]+' | tail -1)
[ -n "${P}" ] || { echo "[chain] submission failed, stopping"; exit 1; }
echo "[chain] showershape cluster ${P}"; poll "${P}" 120; echo "[chain] $(date -Is) showershape jobs cleared"
condor_history -constraint "ClusterId==${P}" -af ExitCode 2>/dev/null | sort | uniq -c | sed 's/^/[chain] exit codes: /'
set +e
bash hadd_showershape_di.sh config_showershape_0rad.yaml;    echo "[chain] hadd 0rad exit=$?"
bash hadd_showershape_di.sh config_showershape_1p5mrad.yaml; echo "[chain] hadd 1p5mrad exit=$?"
bash merge_periods_showershape.sh config_showershape_0rad.yaml config_showershape_1p5mrad.yaml config_showershape.yaml; echo "[chain] all-range merge exit=$?"
cd ../plotting
root -l -b -q 'plot_showershapes_selections.C("config_showershape_0rad.yaml")' > ../efficiencytool/logs/plot_showershapes_selections_resync.log 2>&1; echo "[chain] grids exit=$?"
root -l -b -q plot_SB.C > ../efficiencytool/logs/plot_SB_resync.log 2>&1; echo "[chain] S/B exit=$?"
root -l -b -q plot_npb_time_bkgsub.C > ../efficiencytool/logs/plot_npb_resync.log 2>&1; echo "[chain] NPB exit=$?"
set -e
echo "[chain] $(date -Is) DONE (showershape=${P})"
