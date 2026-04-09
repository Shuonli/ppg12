#!/usr/bin/env bash
#
# monitor_condor.sh — Monitor running condor jobs for PPG12 efficiency pipeline
#
# Usage:
#   bash monitor_condor.sh                  # one-shot status report
#   bash monitor_condor.sh --loop [SEC]     # repeat every SEC seconds (default 120)
#   bash monitor_condor.sh --watch CONFIG   # detailed view of one config's logs
#

LOGDIR="/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/logs"
RESULTDIR="/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results"
CONFIGDIR="/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool"

# Minimum file size (bytes) for a ROOT file to be considered "complete"
# ROOT RECREATE creates ~700 byte headers; completed files are >50KB
MIN_ROOT_SIZE=10000

# Expected samples per pipeline
TREE_SAMPLES="photon5 photon10 photon20 jet5 jet8 jet12 jet20 jet30 jet40 data"
JET_SAMPLES="jet5 jet8 jet12 jet20 jet30 jet40"

# Auto-detect cluster IDs from condor_q (updated each run)
detect_clusters() {
    local qout="$1"
    CLUSTER_TREE=$(echo "$qout" | grep "oneforall_tree\.sh " | head -1 | awk -F'.' '{print $1}')
    CLUSTER_JET=$(echo "$qout" | grep "oneforall_tree_jet\.sh" | head -1 | awk -F'.' '{print $1}')
    # Fallback: detect from most recent condor log events if no jobs in queue
    if [ -z "$CLUSTER_TREE" ] && [ -f "$LOGDIR/config_bdt_nom.yaml.log" ]; then
        CLUSTER_TREE=$(grep "^001 " "$LOGDIR/config_bdt_nom.yaml.log" 2>/dev/null | tail -2 | head -1 | grep -oP '\(\K[0-9]+')
        CLUSTER_JET=$(grep "^001 " "$LOGDIR/config_bdt_nom.yaml.log" 2>/dev/null | tail -1 | grep -oP '\(\K[0-9]+')
    fi
}

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

separator() {
    printf '%*s\n' 80 '' | tr ' ' '='
}

thin_sep() {
    printf '%*s\n' 80 '' | tr ' ' '-'
}

print_header() {
    echo ""
    separator
    echo -e "${BOLD}  PPG12 Condor Job Monitor — $(date '+%Y-%m-%d %H:%M:%S')${NC}"
    separator
}

# ──────────────────────────────────────────────────────────────────────────────
# Section 1: condor_q summary
# ──────────────────────────────────────────────────────────────────────────────
condor_summary() {
    echo -e "\n${CYAN}[1] CONDOR QUEUE STATUS${NC}"
    thin_sep

    local qout
    qout=$(condor_q -nobatch 2>&1)
    detect_clusters "$qout"

    # Count by status
    local total running idle held completed
    total=$(echo "$qout" | grep "^[0-9]" | wc -l)
    running=$(echo "$qout" | grep "^[0-9]" | awk '$6=="R"' | wc -l)
    idle=$(echo "$qout" | grep "^[0-9]" | awk '$6=="I"' | wc -l)
    held=$(echo "$qout" | grep "^[0-9]" | awk '$6=="H"' | wc -l)

    # Cluster breakdown
    local clusters
    clusters=$(echo "$qout" | grep "^[0-9]" | awk -F'.' '{print $1}' | sort -u)

    echo -e "  Total jobs:   ${BOLD}${total}${NC}"
    echo -e "  Running:      ${GREEN}${running}${NC}"
    echo -e "  Idle:         ${YELLOW}${idle}${NC}"
    echo -e "  Held:         ${RED}${held}${NC}"
    echo ""

    for cid in $clusters; do
        local cname ccount crun cidle cheld
        ccount=$(echo "$qout" | grep "^${cid}\." | wc -l)
        crun=$(echo "$qout" | grep "^${cid}\." | awk '$6=="R"' | wc -l)
        cidle=$(echo "$qout" | grep "^${cid}\." | awk '$6=="I"' | wc -l)
        cheld=$(echo "$qout" | grep "^${cid}\." | awk '$6=="H"' | wc -l)
        # Identify pipeline type from first job
        cname=$(echo "$qout" | grep "^${cid}\.0 " | grep -o 'oneforall_tree[^ ]*\.sh')
        echo -e "  Cluster ${BOLD}${cid}${NC} (${cname:-unknown}): ${ccount} jobs  [R:${crun} I:${cidle} H:${cheld}]"
    done

    # Show held jobs with reasons
    if [ "$held" -gt 0 ]; then
        echo ""
        echo -e "  ${RED}HELD JOBS:${NC}"
        echo "$qout" | grep "^[0-9]" | awk '$6=="H" {print "    " $0}'
        echo ""
        echo "  Hold reasons:"
        condor_q -held -af ClusterId ProcId HoldReason 2>/dev/null | while read cid pid reason; do
            echo "    ${cid}.${pid}: ${reason}"
        done
    fi
}

# ──────────────────────────────────────────────────────────────────────────────
# Section 2: Log file analysis
# ──────────────────────────────────────────────────────────────────────────────
log_analysis() {
    echo -e "\n${CYAN}[2] LOG FILE ANALYSIS${NC}"
    thin_sep

    local n_errors=0
    local n_fatal=0
    local error_configs=""

    for cfg in "$CONFIGDIR"/config_bdt_*.yaml; do
        local cfgname
        cfgname=$(basename "$cfg")
        local errfile="${LOGDIR}/${cfgname}.err"
        local outfile="${LOGDIR}/${cfgname}.out"

        if [ ! -f "$errfile" ]; then
            continue
        fi

        # Check for real errors (skip ROOT Info/Warning noise and benign SysErrors)
        local real_errors
        real_errors=$(grep -iE "Error|Fatal|Segmentation|segfault|Abort|SIGSEGV|SIGABRT|core dump|bad_alloc|terminate" "$errfile" 2>/dev/null \
            | grep -viE "^Info in|^Warning in|TEfficiency|SetUseWeightedEvents|TTreeReader::SetEntryBase|could not delete.*No such file" \
            | head -5)

        if [ -n "$real_errors" ]; then
            n_errors=$((n_errors + 1))
            error_configs="${error_configs}  ${RED}${cfgname}${NC}\n"
            while IFS= read -r line; do
                error_configs="${error_configs}    ${line}\n"
            done <<< "$real_errors"
        fi

        # Check for fatal in stdout (ROOT crash)
        if [ -f "$outfile" ]; then
            local out_fatal
            out_fatal=$(grep -iE "Fatal|Abort|Segmentation|core dump|libc.*abort" "$outfile" 2>/dev/null | head -3)
            if [ -n "$out_fatal" ]; then
                n_fatal=$((n_fatal + 1))
                error_configs="${error_configs}  ${RED}${cfgname} (stdout crash)${NC}\n"
                while IFS= read -r line; do
                    error_configs="${error_configs}    ${line}\n"
                done <<< "$out_fatal"
            fi
        fi
    done

    if [ "$n_errors" -eq 0 ] && [ "$n_fatal" -eq 0 ]; then
        echo -e "  ${GREEN}No errors or crashes detected in log files.${NC}"
    else
        echo -e "  ${RED}Errors found in ${n_errors} .err files, ${n_fatal} stdout crashes:${NC}"
        echo -e "$error_configs"
    fi
}

# ──────────────────────────────────────────────────────────────────────────────
# Section 3: Per-config progress (vtxscan + full analysis completion)
# ──────────────────────────────────────────────────────────────────────────────
progress_check() {
    echo -e "\n${CYAN}[3] PER-CONFIG PROGRESS${NC}"
    thin_sep

    local completed_tree=0
    local completed_jet=0
    local partial_tree=0
    local partial_jet=0
    local total_configs=0
    local details=""

    for cfg in "$CONFIGDIR"/config_bdt_*.yaml; do
        local cfgname vartype
        cfgname=$(basename "$cfg")
        vartype=$(grep 'var_type' "$cfg" 2>/dev/null | head -1 | sed 's/.*: *"\?\([^"]*\)"\?.*/\1/')
        total_configs=$((total_configs + 1))

        if [ -z "$vartype" ]; then
            continue
        fi

        # Helper: check if file exists and is large enough to be a real result
        is_complete() { [ -f "$1" ] && [ "$(stat -c%s "$1" 2>/dev/null)" -gt "$MIN_ROOT_SIZE" ]; }

        # --- Condor termination events (definitive completion marker) ---
        local logfile="${LOGDIR}/${cfgname}.log"
        local tree_terminated=false jet_terminated=false
        if [ -f "$logfile" ]; then
            [ -n "$CLUSTER_TREE" ] && grep -q "^005 (${CLUSTER_TREE}\." "$logfile" 2>/dev/null && tree_terminated=true
            [ -n "$CLUSTER_JET" ] && grep -q "^005 (${CLUSTER_JET}\." "$logfile" 2>/dev/null && jet_terminated=true
        fi

        # --- oneforall_tree.sh outputs (photon + jet + data) ---
        local tree_vtxscan=0 tree_full=0 tree_vtxscan_total=0 tree_full_total=0
        for sample in $TREE_SAMPLES; do
            tree_vtxscan_total=$((tree_vtxscan_total + 1))
            tree_full_total=$((tree_full_total + 1))

            if [ "$sample" = "data" ]; then
                local vtx_file="${RESULTDIR}/data_histo_${vartype}_vtxscan.root"
                local full_file="${RESULTDIR}/data_histo_${vartype}.root"
            else
                local vtx_file="${RESULTDIR}/MC_efficiency_${sample}_${vartype}_vtxscan.root"
                local full_file="${RESULTDIR}/MC_efficiency_${sample}_${vartype}.root"
            fi

            is_complete "$vtx_file" && tree_vtxscan=$((tree_vtxscan + 1))
            is_complete "$full_file" && tree_full=$((tree_full + 1))
        done

        # --- oneforall_tree_jet.sh outputs (jet only, no vtxscan) ---
        local jet_full=0 jet_full_total=0
        for sample in $JET_SAMPLES; do
            jet_full_total=$((jet_full_total + 1))
            local jfile="${RESULTDIR}/MC_efficiency_${sample}_${vartype}.root"
            is_complete "$jfile" && jet_full=$((jet_full + 1))
        done

        # Determine tree pipeline status
        local tree_status
        if $tree_terminated; then
            tree_status="${GREEN}DONE${NC}"
            completed_tree=$((completed_tree + 1))
        elif [ "$tree_full" -eq "$tree_full_total" ] && [ "$tree_full" -gt 0 ]; then
            tree_status="${GREEN}DONE${NC}"
            completed_tree=$((completed_tree + 1))
        elif [ "$tree_vtxscan" -gt 0 ] || [ "$tree_full" -gt 0 ]; then
            tree_status="${YELLOW}vtx:${tree_vtxscan}/${tree_vtxscan_total} full:${tree_full}/${tree_full_total}${NC}"
            partial_tree=$((partial_tree + 1))
        else
            tree_status="waiting"
        fi

        # Determine jet pipeline status
        local jet_status
        if $jet_terminated; then
            jet_status="${GREEN}DONE${NC}"
            completed_jet=$((completed_jet + 1))
        elif [ "$jet_full" -eq "$jet_full_total" ] && [ "$jet_full" -gt 0 ]; then
            jet_status="${GREEN}DONE${NC}"
            completed_jet=$((completed_jet + 1))
        elif [ "$jet_full" -gt 0 ]; then
            jet_status="${YELLOW}${jet_full}/${jet_full_total}${NC}"
            partial_jet=$((partial_jet + 1))
        else
            jet_status="waiting"
        fi

        # Estimate progress from stdout (last "Processing entry X / Y" line)
        local pct=""
        local outfile="${LOGDIR}/${cfgname}.out"
        if [ -f "$outfile" ]; then
            local last_progress
            last_progress=$(tail -200 "$outfile" 2>/dev/null | grep "Processing entry" | tail -1)
            if [ -n "$last_progress" ]; then
                local current total_entries
                current=$(echo "$last_progress" | awk '{print $3}')
                total_entries=$(echo "$last_progress" | awk '{print $5}')
                if [ "$total_entries" -gt 0 ] 2>/dev/null; then
                    pct=" (~$(( current * 100 / total_entries ))%)"
                fi
            fi
        fi

        # Only show non-DONE configs in detail
        if ! $tree_terminated || ! $jet_terminated; then
            details="${details}  $(printf '%-38s' "$vartype") tree: ${tree_status}${pct}  |  jet: ${jet_status}\n"
        fi
    done

    echo -e "  Configs: ${BOLD}${total_configs}${NC}"
    echo -e "  Tree pipeline (photon+jet+data):  ${GREEN}${completed_tree}${NC} done, ${YELLOW}${partial_tree}${NC} in progress, $((total_configs - completed_tree - partial_tree)) waiting"
    echo -e "  Jet pipeline (jet only):          ${GREEN}${completed_jet}${NC} done, ${YELLOW}${partial_jet}${NC} in progress, $((total_configs - completed_jet - partial_jet)) waiting"

    if [ -n "$details" ]; then
        echo ""
        echo -e "  ${BOLD}In-progress / waiting configs:${NC}"
        echo -e "$details"
    else
        echo -e "\n  ${GREEN}All configs completed!${NC}"
    fi
}

# ──────────────────────────────────────────────────────────────────────────────
# Section 4: Resource usage
# ──────────────────────────────────────────────────────────────────────────────
resource_usage() {
    echo -e "\n${CYAN}[4] RESOURCE USAGE${NC}"
    thin_sep

    # Memory usage from condor log image size updates
    local max_mem=0
    local max_mem_cfg=""
    for logfile in "$LOGDIR"/config_bdt_*.yaml.log; do
        local mem
        mem=$(grep "MemoryUsage of job" "$logfile" 2>/dev/null | tail -1 | awk '{print $1}')
        if [ -n "$mem" ] && [ "$mem" -gt "$max_mem" ] 2>/dev/null; then
            max_mem=$mem
            max_mem_cfg=$(basename "$logfile" .log)
        fi
    done

    if [ "$max_mem" -gt 0 ]; then
        echo -e "  Peak memory: ${BOLD}${max_mem} MB${NC} (${max_mem_cfg})"
    fi

    # Runtime from condor_q
    local max_runtime
    max_runtime=$(condor_q -nobatch 2>/dev/null | grep "^[0-9]" | awk '{print $5}' | sort | tail -1)
    echo -e "  Max wall time: ${BOLD}${max_runtime:-N/A}${NC}"

    # Disk usage of results
    local results_size
    results_size=$(du -sh "$RESULTDIR" 2>/dev/null | awk '{print $1}')
    echo -e "  Results dir size: ${BOLD}${results_size:-N/A}${NC}"
}

# ──────────────────────────────────────────────────────────────────────────────
# Section 5: Completion estimate
# ──────────────────────────────────────────────────────────────────────────────
completion_estimate() {
    echo -e "\n${CYAN}[5] RECENTLY COMPLETED${NC}"
    thin_sep

    # Check for condor termination events (event 005) across all logs
    local completed_jobs=""
    for logfile in "$LOGDIR"/config_bdt_*.yaml.log; do
        local cfgname
        cfgname=$(basename "$logfile" .log)
        local terminations
        terminations=$(grep "^005 " "$logfile" 2>/dev/null)
        if [ -n "$terminations" ]; then
            while IFS= read -r line; do
                local jobid timestamp retval
                jobid=$(echo "$line" | grep -o '([0-9.]*)')
                timestamp=$(echo "$line" | awk '{print $3, $4}')
                # Check return value from next few lines
                retval=$(grep -A5 "^005 " "$logfile" 2>/dev/null | grep "return value" | tail -1 | awk '{print $NF}')
                if [ "$retval" = "0" ]; then
                    completed_jobs="${completed_jobs}  ${GREEN}✓${NC} ${cfgname} ${jobid} @ ${timestamp} (exit 0)\n"
                else
                    completed_jobs="${completed_jobs}  ${RED}✗${NC} ${cfgname} ${jobid} @ ${timestamp} (exit ${retval:-?})\n"
                fi
            done <<< "$terminations"
        fi
    done

    if [ -n "$completed_jobs" ]; then
        echo -e "$completed_jobs"
    else
        echo -e "  No jobs have terminated yet."
    fi
}

# ──────────────────────────────────────────────────────────────────────────────
# Detailed watch mode for a single config
# ──────────────────────────────────────────────────────────────────────────────
watch_config() {
    local cfg="$1"
    echo -e "\n${BOLD}Detailed view: ${cfg}${NC}"
    separator

    local outfile="${LOGDIR}/${cfg}.out"
    local errfile="${LOGDIR}/${cfg}.err"
    local logfile="${LOGDIR}/${cfg}.log"

    echo -e "\n${CYAN}Last 20 lines of stdout:${NC}"
    tail -20 "$outfile" 2>/dev/null || echo "  (no stdout yet)"

    echo -e "\n${CYAN}Last 10 lines of stderr (filtered):${NC}"
    grep -viE "^Info in|^Warning in|TEfficiency|SetUseWeightedEvents|TTreeReader" "$errfile" 2>/dev/null | tail -10
    echo "(no real errors)" 2>/dev/null

    echo -e "\n${CYAN}Condor events:${NC}"
    grep "^[0-9][0-9][0-9] " "$logfile" 2>/dev/null | tail -10

    echo -e "\n${CYAN}Memory usage:${NC}"
    grep "MemoryUsage" "$logfile" 2>/dev/null | tail -3
}

# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────
run_report() {
    print_header
    condor_summary
    log_analysis
    progress_check
    resource_usage
    completion_estimate
    echo ""
    separator
}

MONITOR_LOG="${LOGDIR}/monitor_condor.log"

case "${1:-}" in
    --loop)
        interval="${2:-120}"
        echo "Monitoring every ${interval}s. Log: ${MONITOR_LOG}"
        echo "Press Ctrl+C to stop."
        while true; do
            {
                run_report
                echo ""
                echo "--- Next check in ${interval}s ---"
            } 2>&1 | tee "$MONITOR_LOG"
            # Exit if no jobs remain in queue
            remaining=$(condor_q -nobatch 2>/dev/null | grep "^[0-9]" | wc -l)
            if [ "$remaining" -eq 0 ]; then
                echo ""
                echo -e "${GREEN}All condor jobs have left the queue. Final report above.${NC}" | tee -a "$MONITOR_LOG"
                break
            fi
            sleep "$interval"
        done
        ;;
    --watch)
        if [ -z "${2:-}" ]; then
            echo "Usage: $0 --watch config_bdt_nom.yaml"
            exit 1
        fi
        watch_config "$2"
        ;;
    *)
        run_report
        ;;
esac
