#!/usr/bin/env bash
# Combined: DI-dispatcher tree stage (per-sample SI+DI with mix_weight) +
# yield stage (via oneforall.sh: merge_periods.sh + RecoEff(data) +
# CalcPhotonYield x2 for bare-name; or per-period MergeSim for feeders).
# Used by the full-variant condor submission.

set -u
CONFIGNAME="$1"
cd /sphenix/user/shuhangli/ppg12/efficiencytool

stamp() { date '+%F %T'; }
echo "[$(stamp)] START $CONFIGNAME"

echo "[$(stamp)] tree stage begin (DI dispatcher)"
bash oneforall_tree_double_dispatch.sh "$CONFIGNAME"
tree_rc=$?
echo "[$(stamp)] tree stage rc=$tree_rc"

if [ "$tree_rc" -ne 0 ]; then
    echo "[$(stamp)] ABORT: tree stage failed; skipping yield stage"
    exit "$tree_rc"
fi

echo "[$(stamp)] yield stage begin"
bash oneforall.sh "$CONFIGNAME"
yield_rc=$?
echo "[$(stamp)] yield stage rc=$yield_rc"
echo "[$(stamp)] DONE $CONFIGNAME"
exit "$yield_rc"
