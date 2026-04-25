#!/usr/bin/bash
# Post-tree-making driver:
#   A) once MC hadd (monitor_mc_data.sh) finishes, sanity-check combined.root
#      and submit MC BDT apply via condor
#   B) in parallel, once data condor (DONE_lines>=15298) finishes, run
#      hadd_data_parts.sh and submit data BDT apply via condor
#
# This script complements monitor_mc_data.sh — it doesn't run hadd itself,
# just watches for the existing monitor to complete then kicks off downstream.

set -u
POLL=${1:-180}
base="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree"
sim_base="$base/sim/run28"
data_ana="$base/data/ana521"
data_cond="$data_ana/condorout"
hadd_log="$base/hadd_combined.log"
driver_log="$base/finish_reprocess.log"
fwx="/sphenix/user/shuhangli/ppg12/FunWithxgboost"

MC_SAMPLES=(photon5 photon10 photon20 jet5 jet8 jet12 jet20 jet30 jet40
            photon10_double jet12_double
            photon5_double photon20_double jet8_double jet20_double
            jet30_double jet40_double jet50_double)
DATA_TARGET=${DATA_TARGET:-15298}

log() { echo "[$(date +'%m-%d %H:%M:%S')] $*" | tee -a "$driver_log"; }

mc_hadd_done() {
    # hadd_combined.sh emits "Summary: ok=N ..." at end
    tail -10 "$hadd_log" 2>/dev/null | grep -q "^Summary: ok="
}

mc_sanity_check() {
    local ok=0 bad=0
    log "--- MC sanity check ---"
    for s in "${MC_SAMPLES[@]}"; do
        local f="$sim_base/$s/condorout/combined.root"
        if [ ! -s "$f" ]; then
            log "  $s: MISSING or empty"; bad=$((bad+1)); continue
        fi
        local n
        n=$(root -l -b -q 2>/dev/null <<EOF
TFile f("$f","READ");
TTree *t=(TTree*)f.Get("slimtree");
std::cout<<"ENT="<<(t?t->GetEntries():-1)<<std::endl;
EOF
)
        local entries=$(echo "$n" | awk -F= '/^ENT=/{print $2}')
        if [ -z "$entries" ] || [ "$entries" -lt 10 ]; then
            log "  $s: FAIL (entries='$entries')"; bad=$((bad+1))
        else
            local size=$(du -h "$f" | cut -f1)
            log "  $s: OK  entries=$entries  size=$size"
            ok=$((ok+1))
        fi
    done
    log "MC sanity summary: ok=$ok  bad=$bad"
    [ "$bad" = 0 ]
}

submit_mc_bdt() {
    log "--- submit MC BDT apply (18 jobs) ---"
    cd "$fwx"
    mkdir -p logs
    condor_submit submit_apply_BDT_mc.sub 2>&1 | tee -a "$driver_log"
}

data_done_lines() {
    find "$data_cond/OutDir"* -name "DONE.txt" -print0 2>/dev/null \
        | xargs -0 cat 2>/dev/null | wc -l
}

data_hadd() {
    log "--- run hadd_data_parts.sh ---"
    cd "$data_ana"
    mkdir -p /tmp/ppg12_pipeline_logs
    bash hadd_data_parts.sh 2>&1 | tee -a "$driver_log" | tail -200 > /dev/null
    local nparts
    nparts=$(ls "$data_cond"/part_*.root 2>/dev/null | wc -l)
    log "data hadd produced $nparts part_*.root files"
}

submit_data_bdt() {
    log "--- build datalist.list + submit data BDT apply ---"
    cd "$fwx"
    ls "$data_cond"/part_*.root > datalist.list
    local n=$(wc -l < datalist.list)
    log "datalist.list has $n entries"
    mkdir -p logs
    condor_submit apply_BDT_condor.sub 2>&1 | tee -a "$driver_log"
}

log "==== finish_reprocess start  poll=${POLL}s  target=$DATA_TARGET ===="

mc_bdt_submitted=0
data_bdt_submitted=0

while true; do
    # ---------- MC side ----------
    if [ "$mc_bdt_submitted" = 0 ]; then
        if mc_hadd_done; then
            log "MC hadd finished (saw Summary line)."
            if mc_sanity_check; then
                submit_mc_bdt
                mc_bdt_submitted=1
            else
                log "MC sanity FAILED -- not submitting BDT. Investigate $hadd_log."
                # don't spin forever on a bad hadd
                mc_bdt_submitted=-1
            fi
        fi
    fi

    # ---------- Data side ----------
    if [ "$data_bdt_submitted" = 0 ]; then
        nd=$(data_done_lines)
        if [ "$nd" -ge "$DATA_TARGET" ]; then
            log "Data condor finished ($nd / $DATA_TARGET)"
            data_hadd
            submit_data_bdt
            data_bdt_submitted=1
        fi
    fi

    # ---------- exit when both downstream flows launched ----------
    if [ "$mc_bdt_submitted" != 0 ] && [ "$data_bdt_submitted" != 0 ]; then
        log "Both downstream flows launched (MC=$mc_bdt_submitted DATA=$data_bdt_submitted). Exiting."
        exit 0
    fi

    sleep "$POLL"
done
