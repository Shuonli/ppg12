#!/usr/bin/env bash
# Round-2 recovery: hadd all-range data_histo from fresh per-period siblings,
# then re-run CalcPhotonYield × 2. The previous recovery rebuilt MC merges
# and re-ran CalcPhotonYield, but missed that the all-range data_histo was
# also stale (pre 5/13 BDT cut fix) because oneforall.sh aborted before its
# data hadd step in the original Phase 2.

set -eo pipefail
cd /sphenix/user/shuhangli/ppg12/efficiencytool
source /sphenix/u/shuhang98/setup.sh
set -u

CONFIGS=/tmp/broken_configs.list
LOG=logs/recovery_round2_data_histo.log

echo "[round2] $(date -Is) Starting data-histo + CalcPhotonYield refresh for $(wc -l < $CONFIGS) configs"

OK_HADD=0
FAIL_HADD=0
OK_CPY=0
FAIL_CPY=0

while IFS= read -r cfg; do
    base=${cfg%.yaml}
    var=${base#config_bdt_}
    me=results/data_histo_bdt_${var}.root
    p0=results/data_histo_bdt_${var}_0rad.root
    p1=results/data_histo_bdt_${var}_1p5mrad.root

    if [ ! -f "$p0" ] || [ ! -f "$p1" ]; then
        echo "  [SKIP] $var: per-period siblings missing"
        continue
    fi

    echo "[round2] $(date -Is)   hadd data_histo for $var"
    if hadd -f "$me" "$p0" "$p1" > /dev/null 2>&1; then
        OK_HADD=$((OK_HADD + 1))
    else
        echo "[round2]   [WARN] hadd failed for $var"
        FAIL_HADD=$((FAIL_HADD + 1))
        continue
    fi

    echo "[round2] $(date -Is)   CalcPhotonYield($var, true)"
    if root -l -b -q "CalculatePhotonYield.C(\"$cfg\", true)" > /dev/null 2>&1; then
        :
    else
        echo "[round2]   [WARN] CalcPhotonYield($var, true) failed"
        FAIL_CPY=$((FAIL_CPY + 1))
        continue
    fi
    echo "[round2] $(date -Is)   CalcPhotonYield($var, false)"
    if root -l -b -q "CalculatePhotonYield.C(\"$cfg\", false)" > /dev/null 2>&1; then
        OK_CPY=$((OK_CPY + 1))
    else
        echo "[round2]   [WARN] CalcPhotonYield($var, false) failed"
        FAIL_CPY=$((FAIL_CPY + 1))
    fi
done < "$CONFIGS"

echo "[round2] $(date -Is) Done. hadd ok=${OK_HADD}/fail=${FAIL_HADD}, CPY ok=${OK_CPY}/fail=${FAIL_CPY}"

# Re-aggregate systematics
echo "[round2] $(date -Is) Re-running calc_syst_bdt.py"
cd /sphenix/user/shuhangli/ppg12/plotting
python3 calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures > /tmp/calc_syst_round2.log 2>&1
echo "[round2] $(date -Is) Complete."
