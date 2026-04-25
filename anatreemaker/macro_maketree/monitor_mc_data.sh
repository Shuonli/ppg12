#!/usr/bin/bash
# Monitor MC (18-sample condor) + data (ana521) reprocess.
# When MC is fully done, auto-trigger hadd_combined.sh for all 18 samples.
#
# Usage: ./monitor_mc_data.sh [poll_sec] [no_hadd]
#   poll_sec: seconds between checks (default 300 = 5 min)
#   no_hadd:  pass "no_hadd" as 2nd arg to skip auto-hadd

POLL=${1:-300}
NO_HADD=${2:-}
base="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree"
sim_base="$base/sim/run28"
data_cond="$base/data/ana521/condorout"
hadd_log="$base/hadd_combined.log"
mon_log="$base/monitor.log"

SI=("photon5" "photon10" "photon20" "jet5" "jet8" "jet12" "jet20" "jet30" "jet40")
DI_OLD=("photon10_double" "jet12_double")
DI_NEW=("photon5_double" "photon20_double" "jet8_double" "jet20_double" "jet30_double" "jet40_double" "jet50_double")
ALL=("${SI[@]}" "${DI_OLD[@]}" "${DI_NEW[@]}")

log() { echo "[$(date +'%H:%M:%S')] $*" | tee -a "$mon_log"; }

mc_done_check() {
    local total_outdir=0 total_calo=0 total_ok=0
    local all_ok=1
    local per_sample=""
    for s in "${ALL[@]}"; do
        local nout=$(ls -d "$sim_base/$s/condorout/OutDir"*/ 2>/dev/null | wc -l)
        # Count caloana.root >= 1MB (exclude truncated/quota-failed 8MB-limit files are still fine,
        # but zero-size is excluded; the real disk-quota truncations were at 8388608 exactly
        # — here we just require non-zero file that opened at least one TKey).
        local nroot=$(find "$sim_base/$s/condorout/OutDir"* -maxdepth 1 -name caloana.root -size +1c 2>/dev/null | wc -l)
        total_outdir=$((total_outdir + nout))
        total_calo=$((total_calo + nroot))
        per_sample+="$(printf '%s=%d/%d' "$s" "$nroot" "$nout") "
        [ "$nroot" -lt "$nout" ] && all_ok=0
    done
    # MC queue (informational only; many schedds show zombie "completed" jobs after caloana.root is written)
    local qact
    qact=$(condor_q -global shuhang98 -const 'JobStatus==2' -af ClusterId 2>/dev/null | wc -l)
    log "MC caloana=$total_calo / outdir=$total_outdir   running(JS=2)=$qact"
    printf '  ' >> "$mon_log"; echo "$per_sample" >> "$mon_log"
    if [ "$all_ok" = 1 ]; then return 0; fi
    return 1
}

# run.sh writes "DONE" to each OutDir's DONE.txt once per sub-job (~10 per OutDir).
# Target is the total sub-job count across runList.txt (run.sh hardcodes ~15298 for ana521).
DATA_TARGET=${DATA_TARGET:-15298}

data_status() {
    local nout=$(ls -d "$data_cond/OutDir"* 2>/dev/null | wc -l)
    local ndone_lines=$(find "$data_cond/OutDir"* -name "DONE.txt" -print0 2>/dev/null | xargs -0 cat 2>/dev/null | wc -l)
    local nroot=$(find "$data_cond" -maxdepth 2 -name "OUTTREE_*.root" 2>/dev/null | wc -l)
    local fresh=$(find "$data_cond" -maxdepth 2 -name "DONE.txt" -mmin -${POLL_MIN:-10} 2>/dev/null | wc -l)
    local pct=$(( (ndone_lines * 100) / (DATA_TARGET > 0 ? DATA_TARGET : 1) ))
    log "DATA  DONE_lines=$ndone_lines / target=$DATA_TARGET (${pct}%)  slimtrees=$nroot  outdirs=$nout  fresh<${POLL_MIN:-10}min: $fresh"
}

hadd_done=0
log "==== monitor start  poll=${POLL}s  no_hadd=${NO_HADD:-no} ===="

while true; do
    POLL_MIN=$((POLL/60 + 1))
    data_status
    if mc_done_check; then
        if [ "$hadd_done" = 0 ] && [ -z "$NO_HADD" ]; then
            log "MC fully done — triggering hadd_combined.sh all"
            "$sim_base/hadd_combined.sh" all >> "$hadd_log" 2>&1
            log "hadd finished — see $hadd_log ; last lines:"
            tail -8 "$hadd_log" | tee -a "$mon_log"
            hadd_done=1
        else
            log "MC done; NO_HADD set or already ran -> nothing to do"
        fi
    fi
    # exit condition: MC done (hadd ran) AND data reached target (or stalled)
    if [ "$hadd_done" = 1 ] || [ -n "$NO_HADD" ]; then
        ndone_lines=$(find "$data_cond/OutDir"* -name "DONE.txt" -print0 2>/dev/null | xargs -0 cat 2>/dev/null | wc -l)
        if [ "$ndone_lines" -ge "$DATA_TARGET" ]; then
            log "DATA reached target ($ndone_lines >= $DATA_TARGET).  Exiting."
            exit 0
        fi
        n_fresh=$(find "$data_cond" -maxdepth 2 -name "DONE.txt" -mmin -$((POLL_MIN*3)) 2>/dev/null | wc -l)
        if [ "$n_fresh" = 0 ]; then
            log "DATA stalled: no fresh DONE.txt in last $((POLL_MIN*3)) min (current=$ndone_lines/$DATA_TARGET).  Exiting."
            exit 0
        fi
    fi
    sleep "$POLL"
done
