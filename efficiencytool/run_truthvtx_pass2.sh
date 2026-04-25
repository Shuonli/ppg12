#!/usr/bin/env bash
# Pass 2 only (vtxscan files already symlinked from bdt_vtxreweight0), then yield stage.
# Limited to 3 parallel jobs to respect 10GB free memory budget.

source /sphenix/u/shuhang98/setup.sh
set -u
CONFIGNAME="$1"
cd /sphenix/user/shuhangli/ppg12/efficiencytool

stamp() { date '+%F %T'; }
log() { echo "[$(stamp)] $*"; }

log "START Pass2 $CONFIGNAME"

run_sample() {
    local sample="$1"
    log "sample $sample begin"
    root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "'"${sample}"'")' \
        > "logs/recoff_${sample}_$(basename ${CONFIGNAME} .yaml).log" 2>&1
    log "sample $sample rc=$?"
}
export -f run_sample stamp log
export CONFIGNAME

mkdir -p logs

# 3 waves of 3-4 parallel jobs to stay within memory budget
log "wave A: photon5, photon10, photon20"
( run_sample photon5 ) &
( run_sample photon10 ) &
( run_sample photon20 ) &
wait
log "wave A done"

log "wave B: jet5, jet8, jet12"
( run_sample jet5 ) &
( run_sample jet8 ) &
( run_sample jet12 ) &
wait
log "wave B done"

log "wave C: jet20, jet30, jet40"
( run_sample jet20 ) &
( run_sample jet30 ) &
( run_sample jet40 ) &
wait
log "wave C done"

log "wave D: data"
run_sample data
log "wave D done"

log "yield stage begin"
bash oneforall.sh "$CONFIGNAME" > "logs/yield_$(basename ${CONFIGNAME} .yaml).log" 2>&1
log "yield stage rc=$?"
log "DONE $CONFIGNAME"
