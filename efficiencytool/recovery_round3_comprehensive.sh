#!/usr/bin/env bash
# Round-3 comprehensive recovery: re-hadd ALL stale merged files across
# 4 file families (data_histo, MC_efficiency, MC_efficiency_jet, MC_response)
# detected by the audit, then re-run CalcPhotonYield × 2 for every affected
# bare-name config, then re-aggregate.
#
# Stale = merged file's mtime < max(per-period siblings' mtime) − 60s.

set -eo pipefail
cd /sphenix/user/shuhangli/ppg12/efficiencytool
source /sphenix/u/shuhang98/setup.sh
set -u

LOG=logs/recovery_round3.log

echo "[round3] $(date -Is) Building stale-file list..."

python3 << 'PYEOF' > /tmp/round3_stale_configs.list
import os
results = "/sphenix/user/shuhangli/ppg12/efficiencytool/results"
families = ["data_histo_bdt_", "MC_efficiency_bdt_", "MC_efficiency_jet_bdt_", "MC_response_bdt_"]
bare_names = set()
for fn in os.listdir(results):
    if fn.endswith("_0rad.root") and fn.startswith("data_histo_bdt_"):
        base = fn[len("data_histo_bdt_"):-len("_0rad.root")]
        if os.path.exists(os.path.join(results, f"data_histo_bdt_{base}_1p5mrad.root")):
            bare_names.add(base)

stale_configs = set()
for base in bare_names:
    for prefix in families:
        m = os.path.join(results, f"{prefix}{base}.root")
        p0 = os.path.join(results, f"{prefix}{base}_0rad.root")
        p1 = os.path.join(results, f"{prefix}{base}_1p5mrad.root")
        if not (os.path.exists(m) and os.path.exists(p0) and os.path.exists(p1)):
            continue
        if os.path.getmtime(m) < max(os.path.getmtime(p0), os.path.getmtime(p1)) - 60:
            stale_configs.add(base)
            break  # one stale family is enough to flag the config

for base in sorted(stale_configs):
    print(base)
PYEOF
N=$(wc -l < /tmp/round3_stale_configs.list)
echo "[round3] $(date -Is) Found ${N} configs with at least one stale merged file."

OK_HADD=0
FAIL_HADD=0
OK_CPY=0
FAIL_CPY=0

while IFS= read -r base; do
    families="data_histo_bdt_ MC_efficiency_bdt_ MC_efficiency_jet_bdt_ MC_response_bdt_"
    for prefix in $families; do
        m=results/${prefix}${base}.root
        p0=results/${prefix}${base}_0rad.root
        p1=results/${prefix}${base}_1p5mrad.root
        if [ ! -f "$p0" ] || [ ! -f "$p1" ]; then continue; fi
        # Re-hadd unconditionally (safer than checking which family is stale)
        if hadd -f "$m" "$p0" "$p1" > /dev/null 2>&1; then
            OK_HADD=$((OK_HADD + 1))
        else
            FAIL_HADD=$((FAIL_HADD + 1))
        fi
    done

    cfg=config_bdt_${base}.yaml
    if [ ! -f "$cfg" ]; then
        echo "[round3]   [SKIP] $base: $cfg not found"
        continue
    fi

    echo "[round3] $(date -Is)   ${base}: CalcPhotonYield × 2"
    if root -l -b -q "CalculatePhotonYield.C(\"${cfg}\", true)" > /dev/null 2>&1; then
        :
    else
        echo "[round3]   [WARN] ${base}: isMC=true failed"
        FAIL_CPY=$((FAIL_CPY + 1))
        continue
    fi
    if root -l -b -q "CalculatePhotonYield.C(\"${cfg}\", false)" > /dev/null 2>&1; then
        OK_CPY=$((OK_CPY + 1))
    else
        echo "[round3]   [WARN] ${base}: isMC=false failed"
        FAIL_CPY=$((FAIL_CPY + 1))
    fi
done < /tmp/round3_stale_configs.list

echo "[round3] $(date -Is) Done. hadd ok=${OK_HADD}/fail=${FAIL_HADD}, CPY ok=${OK_CPY}/fail=${FAIL_CPY}"

echo "[round3] $(date -Is) Re-running calc_syst_bdt.py"
cd /sphenix/user/shuhangli/ppg12/plotting
python3 calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures > /tmp/calc_syst_round3.log 2>&1
echo "[round3] $(date -Is) Complete."
