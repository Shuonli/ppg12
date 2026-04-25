#!/usr/bin/env bash
# Combined: RecoEffCalculator_TTreeReader (tree stage) + MergeSim + CalculatePhotonYield (yield stage).
# One-stop wrapper used by the inner-R scan condor submission.

set -u
CONFIGNAME="$1"
cd /sphenix/user/shuhangli/ppg12/efficiencytool

stamp() { date '+%F %T'; }
echo "[$(stamp)] START $CONFIGNAME"

echo "[$(stamp)] tree stage begin"
bash oneforall_tree.sh "$CONFIGNAME"
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
