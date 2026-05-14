#!/usr/bin/env bash
# Recovery from 2026-05-13 variant Phase 2 TFileMerger-via-hadd race that
# left 43 merged files silently broken. Steps:
#   1. Re-hadd all 43 broken MC files from the proposed hadd commands in
#      the audit log. Sequential so concurrent disk I/O can't re-break.
#   2. Re-run CalculatePhotonYield × 2 (isMC=true + isMC=false) for the 36
#      distinct affected configs. Sequential for safety.
#   3. Re-audit.
#
# Inputs:
#   /tmp/broken_configs.list — 36 unique affected configs
#   logs/audit_merges_post_cut_fix.log — audit log with proposed hadd cmds

set -eo pipefail
cd /sphenix/user/shuhangli/ppg12/efficiencytool
source /sphenix/u/shuhang98/setup.sh
set -u

AUDIT=logs/audit_merges_post_cut_fix.log
CONFIGS=/tmp/broken_configs.list
LOG=logs/recovery_post_cut_fix.log

# ---- Step 1: hadd recovery ----
echo "[recovery] $(date -Is) Step 1: hadd recovery for 43 broken merged files"
HADD_OK=0
HADD_FAIL=0
# Extract `hadd -f ...` lines (indented in audit log) and run each.
while IFS= read -r cmd; do
    echo "  $cmd"
    if eval "$cmd" > /dev/null 2>&1; then
        HADD_OK=$((HADD_OK + 1))
    else
        echo "  [WARN] hadd failed for: $cmd"
        HADD_FAIL=$((HADD_FAIL + 1))
    fi
done < <(grep -oE "^\s+hadd -f .*" "$AUDIT" | sed 's/^\s*//')

echo "[recovery] $(date -Is) hadd done: ok=${HADD_OK}, fail=${HADD_FAIL}"

# ---- Step 2: re-run CalculatePhotonYield × 2 for affected configs ----
echo "[recovery] $(date -Is) Step 2: re-run CalculatePhotonYield for $(wc -l < $CONFIGS) configs"
CPY_OK=0
CPY_FAIL=0
while IFS= read -r cfg; do
    echo "[recovery] $(date -Is)   ${cfg}: isMC=true"
    if root -l -b -q "CalculatePhotonYield.C(\"${cfg}\", true)" > /dev/null 2>&1; then
        :
    else
        echo "[recovery]   [WARN] CalcPhotonYield(${cfg}, true) failed"
        CPY_FAIL=$((CPY_FAIL + 1))
        continue
    fi
    echo "[recovery] $(date -Is)   ${cfg}: isMC=false"
    if root -l -b -q "CalculatePhotonYield.C(\"${cfg}\", false)" > /dev/null 2>&1; then
        CPY_OK=$((CPY_OK + 1))
    else
        echo "[recovery]   [WARN] CalcPhotonYield(${cfg}, false) failed"
        CPY_FAIL=$((CPY_FAIL + 1))
    fi
done < "$CONFIGS"

echo "[recovery] $(date -Is) CalcPhotonYield done: ok=${CPY_OK}, fail=${CPY_FAIL}"

# ---- Step 3: re-audit ----
echo "[recovery] $(date -Is) Step 3: re-audit"
python3 audit_merges.py > logs/audit_merges_post_recovery.log 2>&1
AUDIT_RC=$?

# Pull headline from new audit
BROKEN_AFTER=$(grep -oE "[0-9]+ broken" logs/audit_merges_post_recovery.log | head -1 || echo "?")
echo "[recovery] $(date -Is) re-audit exit=${AUDIT_RC}, ${BROKEN_AFTER} after recovery"
echo "[recovery] $(date -Is) Recovery complete."
