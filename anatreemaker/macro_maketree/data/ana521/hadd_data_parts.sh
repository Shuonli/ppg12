#!/usr/bin/bash
# Merge per-job OUTTREE_DST_Jet_*.root files into part_N.root files for BDT
# apply input. 200 files per part, serial hadd (-j 1, -n 50) to avoid /tmp
# exhaustion.
#
# Safe to run while some condor jobs are still writing:
#   - Canonical full sort of all OUTTREE files determines the part assignment.
#     Stable across reruns as long as no new files appear.
#   - A part is hadd'd only if none of its input files live in an OutDir with
#     an active condor job (JobStatus==1/idle or 2/running). Parts that still
#     touch active OutDirs are deferred and reported.
#   - An already-produced part_N.root is left alone if it is newer than every
#     one of its inputs (no redo).
#   - Re-run after more jobs finish to sweep up the deferred parts.
#
# Matching is done by OutDir basename (e.g. "OutDir10141"), not full path, to
# avoid /sphenix/... vs /gpfs/mnt/gpfs02/sphenix/... symlink mismatches between
# $(pwd) and condor's Iwd.
#
# Note: hadd itself skips corrupt/zombie ROOT inputs with an error and still
# produces the merged output from the valid inputs. Grep the log for
# "Zombie\|cannot open file" to audit missing files after a run.

set -e
cd /sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/condorout
FILES_PER_PART=200
LIST=/tmp/ppg12_pipeline_logs/outtree_files.list
ACTIVE=/tmp/ppg12_pipeline_logs/hadd_active_iwd.txt
ACTIVE_NAMES=/tmp/ppg12_pipeline_logs/hadd_active_outdir_names.txt
SKIPLOG=/tmp/ppg12_pipeline_logs/hadd_skipped.list
: > "$SKIPLOG"

echo "Collecting OutDirs with active condor jobs (idle or running)..."
condor_q "$(whoami)" -constraint 'JobStatus==1 || JobStatus==2' -af Iwd 2>/dev/null | sort -u > "$ACTIVE"
# Extract just the OutDirNNN basename from each Iwd path.
awk -F'/' '{print $NF}' "$ACTIVE" | sort -u > "$ACTIVE_NAMES"
n_active=$(wc -l < "$ACTIVE_NAMES")
echo "Active OutDirs: $n_active"

echo "Finding OUTTREE files..."
find . -maxdepth 2 -name "OUTTREE_DST_Jet*.root" | sort > "$LIST"
N=$(wc -l < "$LIST")
echo "Total OUTTREE files: $N"
NPARTS=$(( (N + FILES_PER_PART - 1) / FILES_PER_PART ))
echo "Planning $NPARTS parts of ${FILES_PER_PART} files each (last may have fewer)"
date

done_parts=0
skipped_parts=0

for ((p=0; p<NPARTS; p++)); do
    start=$(( p * FILES_PER_PART + 1 ))
    end=$(( start + FILES_PER_PART - 1 ))
    chunk=/tmp/ppg12_pipeline_logs/outtree_chunk_$p.list
    sed -n "${start},${end}p" "$LIST" > "$chunk"
    nc=$(wc -l < "$chunk")
    out="part_${p}.root"

    # Skip already-fresh part.
    if [ -f "$out" ]; then
        out_mtime=$(stat -c %Y "$out")
        newest_input=$(stat -c %Y $(cat "$chunk") 2>/dev/null | sort -n | tail -1)
        if [ -n "$newest_input" ] && [ "$out_mtime" -gt "$newest_input" ]; then
            echo "=== part $p up-to-date ($out newer than all inputs), skipping ==="
            rm -f "$chunk"
            done_parts=$((done_parts+1))
            continue
        fi
    fi

    # Extract the OutDir basename of each file in the chunk and check against
    # the active set. Input lines look like "./OutDirNNN/OUTTREE_...root", so
    # the 2nd /-separated field is the OutDir basename.
    active_hit=""
    if [ "$n_active" -gt 0 ]; then
        active_hit=$(awk -F'/' '{print $2}' "$chunk" | sort -u | grep -Fxf "$ACTIVE_NAMES" | head -1 || true)
    fi

    if [ -n "$active_hit" ]; then
        echo "=== part $p: has file(s) in active $active_hit, deferring ==="
        echo "$p" >> "$SKIPLOG"
        rm -f "$chunk"
        skipped_parts=$((skipped_parts+1))
        continue
    fi

    echo "=== part $p ($nc files) -> $out ==="
    rm -f "$out"
    hadd -f -j 1 -n 50 "$out" $(cat "$chunk") 2>&1 | tail -3
    ls -la "$out" 2>/dev/null | awk '{print "  ->", $NF, $5, "bytes"}'
    rm -f "$chunk"
    done_parts=$((done_parts+1))
done

echo "------------------------------------------------------------"
echo "hadd summary: $done_parts / $NPARTS parts done, $skipped_parts deferred"
if [ "$skipped_parts" -gt 0 ]; then
    echo "Deferred parts (rerun this script after their jobs finish):"
    paste -sd',' "$SKIPLOG"
fi
date
