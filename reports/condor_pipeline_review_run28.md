# run28 HTCondor Pipeline Review

**Path:** `anatreemaker/macro_maketree/sim/run28/`
**Scope:** 14 sample subdirs (photon5/10/20, photon10_double, jet5/8/10/12/12_double/15/20/30/40/50) producing slimtrees from G4 DSTs via `Fun4All_run_sim.C` → `CaloAna24` → `caloana.root`
**Date:** 2026-04-17

## Executive summary

The pipeline works but is structurally fragile. Fourteen near-duplicated copies of `run_condor.sh` (only `dstpath`/`TargetDir` change) and one identical `CondorRunSim_new.sh` across all samples make any future fix a 14-place edit. Two variants of `run_condor.sh` already drifted (`jet5/jet8` auto-derive job count; the rest hardcode `total_line=10000`, which produces empty-input jobs on samples with <10000 DST files, e.g. `photon10_double` = 9998, `jet12_double` = 9997). Five issues are correctness- or debuggability-critical (empty caloana.root merged by `hadd_combined.sh`, interactive ROOT hang risk, condor log written to scheduler `/tmp`, `ff.sub` appended instead of overwritten, 10k separate `condor_submit` RPCs).

## Finding table

| # | Severity | Finding | File/Line | Impact |
|---|----------|---------|-----------|--------|
| 1 | CRITICAL | `Log = /tmp/sli_$i.log` writes to scheduler's `/tmp` | `*/run_condor.sh:49` | Logs lost on scheduler reboot/cleanup; collisions with any other user index; no breadcrumb for debugging held/failed jobs |
| 2 | CRITICAL | `cat >>ff.sub<<EOF` appends rather than overwrites | `*/run_condor.sh:38` | On resubmit without full cleanup (partial `OutDir*` left) subsequent runs concatenate stanzas → invalid submit file silently rejected |
| 3 | CRITICAL | `root "Fun4All_run_sim.C"` — no `-b -l -q` | `*/CondorRunSim_new.sh:17` | On ROOT error, job hangs at interactive prompt until walltime; wastes slot and delays debug |
| 4 | CRITICAL | No ROOT exit-code check; `echo "JOB COMPLETE!"` always runs | `*/CondorRunSim_new.sh:17-19` | Failed job exits 0; `hadd_combined.sh` blindly merges zombie/truncated `caloana.root` |
| 5 | CRITICAL | `condor_submit` inside the 1000-iteration loop | `*/run_condor.sh:55` | ~1000 sequential scheduler RPCs per sample × 14 samples; sub01 load spike, slow submit, unnecessary auth overhead |
| 6 | HIGH | `PeriodicHold = (NumJobStarts>=2 && JobStatus == 1)` without matching `PeriodicRelease` | `*/run_condor.sh:43` | Failed jobs accumulate in HELD state forever; requires manual `condor_release`; no `on_exit_hold_reason` for triage |
| 7 | HIGH | `total_line=10000` hardcoded in most variants | `*/run_condor.sh:10` | `photon10_double` (9998 lines), `jet12_double` (9997) → last 2-3 jobs read empty filelists, produce empty output merged into `combined.root`. `jet5`/`jet8` already fixed this with `wc -l` |
| 8 | HIGH | Two diverging `run_condor.sh` variants (A: multi-list hardcoded total_line; B: single test.list dynamic) | `jet5/jet8` vs rest | Future maintenance must track both; jet8 uses `test.list` only (no dst_calo_cluster/truth_jet/g4hits/mbd_epd lists) — check it still emits the same slimtree schema |
| 9 | HIGH | `+JobFlavour = "workday"` is CERN LXBATCH syntax | `*/run_condor.sh:39` | Silently ignored by BNL SDCC HTCondor; misleading when debugging wall-time issues |
| 10 | HIGH | No `request_disk` / `request_cpus` | all submit files | Scheduler defaults can evict; silent eviction with no message |
| 11 | HIGH | `source /sphenix/u/shuhang98/setup.sh` hardcoded | `*/CondorRunSim_new.sh:7` | Anyone else (or a renamed account) can't re-run these jobs; setup.sh is outside the repo so no commit history |
| 12 | HIGH | 14 duplicated `run_condor.sh` (only paths change) + 14 identical `CondorRunSim_new.sh` | every sample dir | Drift inevitable; any policy change (e.g. memory bump, new retry rule) is a 14-place edit |
| 13 | MED | Per-job shell copy `CondorRunTC$i.sh` | `*/run_condor.sh:28` | 1000 identical shell copies per sample; wastes quota and inodes on shared GPFS |
| 14 | MED | `Notify_user = sl4859@columbia.edu` paired with `Notification = Never` | `*/run_condor.sh:45,50` | `Notify_user` is inert; misleading |
| 15 | MED | `Priority = +80` with no documented rationale | `*/run_condor.sh:46` | SDCC policy may cap/reject; if capped, creates false confidence that jobs are prioritised |
| 16 | MED | `CondorRunSim_new.sh` missing `set -euo pipefail` | `*/CondorRunSim_new.sh:1` | If `source setup.sh` fails, job still proceeds and succeeds vacuously |
| 17 | MED | `cleanup.sh` / `cleanup_and_run.sh` / `hadd_combined.sh` each hardcode different sample lists | top-level | `cleanup_and_run.sh` lacks jet8/10/15/50/doubles; `hadd_combined.sh` has only the two doubles (main list commented out). Drift-prone, error-prone |
| 18 | MED | `rm -rf condorout/OutDir*` unconditional | `cleanup*.sh` | Any experimental dir under condorout starting with `OutDir` is blown away without confirmation |
| 19 | MED | `jet12_double` Fun4All macro adds `TruthJetInput::add_embedding_flag(2)` (pythia-in-hijing embed) | `photon10_double/Fun4All_run_sim.C:458`, `jet12_double/Fun4All_run_sim.C:458`, `jet8/Fun4All_run_sim.C:458` | Physics-motivated (embedded double-interaction samples), but buried in 3 of 14 macro copies — worth an explicit header comment + a single parameterised macro |
| 20 | LOW | No custom classads (`+Sample`, `+ProductionTag`) | submit files | Can't filter `condor_q` by sample; harder to identify orphans |
| 21 | LOW | `test.out` / `test.err` / `ff.sub` / `jobtime.root` filenames | submit files | Hard to grep; prefer `caloana_$(ClusterId).$(ProcId).out` or fixed `job.{out,err,log}` |

## Recommended refactor (minimal)

Consolidate to one top-level shared template + thin per-sample wrapper. This removes ~40×14 = 560 lines of duplicated bash and fixes all CRITICAL + most HIGH in one pass.

### 1. Shared `submit.sub` template at `run28/submit.sub.template`

```condor
Universe                = vanilla
Executable              = $(BaseDir)/CondorRunSim_new.sh
getenv                  = False
request_memory          = 4GB
request_disk            = 4GB
request_cpus            = 1
initialdir              = $(OutDir)
Output                  = job.out
Error                   = job.err
Log                     = $(OutDir)/job.log

Notification            = Never
max_retries             = 2
on_exit_hold            = (ExitCode =!= 0)
on_exit_hold_reason     = strcat("ROOT exit ", string(ExitCode))

+Sample                 = "$(SAMPLE)"
+ProductionTag          = "run28"

Queue OutDir matching dirs $(BaseDir)/condorout/OutDir*
```

One `condor_submit` per sample instead of 1000-10000.

### 2. Shared `CondorRunSim_new.sh` at `run28/CondorRunSim_new.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail

# Portable sPHENIX env
source /sphenix/user/shuhangli/ppg12/setup_ppg12.sh    # checked-in replacement for ~shuhang98/setup.sh

cd "${_CONDOR_SCRATCH_DIR:-$PWD}"  # condor scratch if staged; else WorkDir on GPFS
echo "[run] CWD=$PWD  HOST=$(hostname)  JOB=${_CONDOR_JOB_AD:-local}"

root -b -l -q Fun4All_run_sim.C
status=$?
if [[ $status -ne 0 ]]; then
    echo "[run] ROOT failed with exit $status"
    exit $status
fi

# Validate output
python3 -c '
import sys, ROOT
f = ROOT.TFile.Open("caloana.root")
if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
    print("caloana.root invalid"); sys.exit(2)
t = f.Get("slimtree")
if not t or t.GetEntries() == 0:
    print("slimtree empty"); sys.exit(3)
print(f"[run] slimtree entries = {t.GetEntries()}")
' || exit $?

echo "[run] JOB COMPLETE"
```

### 3. Shared `run_condor.sh` at `run28/run_condor.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
sample="${1:-}"
lines_per_job="${2:-10}"
[[ -z "$sample" ]] && { echo "Usage: $0 <sample> [lines_per_job]" >&2; exit 1; }

base=/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28
dir="$base/$sample"
[[ -d "$dir" ]] || { echo "no such sample: $sample" >&2; exit 1; }

cd "$dir"
total_line=$(wc -l < dst_calo_cluster.list)
total_jobs=$(( (total_line + lines_per_job - 1) / lines_per_job ))

echo "[submit] $sample: total_line=$total_line jobs=$total_jobs lines/job=$lines_per_job"

target="$dir/condorout"
mkdir -p "$target"

for ((i=0; i<total_jobs; i++)); do
    start=$(( i * lines_per_job + 1 ))
    end=$(( start + lines_per_job - 1 ))
    (( start > total_line )) && break
    workdir="$target/OutDir$i"
    mkdir -p "$workdir"
    for f in dst_calo_cluster dst_truth_jet g4hits dst_mbd_epd; do
        [[ -f "$f.list" ]] || continue
        sed -n "${start},${end}p" "$f.list" > "$workdir/$f.list"
    done
    # Symlink (not copy) the macro so updates propagate
    ln -sf "$dir/Fun4All_run_sim.C" "$workdir/Fun4All_run_sim.C"
done

condor_submit \
    -a "BaseDir=$base" \
    -a "SAMPLE=$sample" \
    "$base/submit.sub.template"
```

### 4. Per-sample wrapper (optional, for convenience)

Each sample dir keeps a trivial one-liner:
```bash
# photon10/run_condor.sh
exec /sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/run_condor.sh photon10 "$@"
```

### 5. Rewrite top-level orchestrators to use one source-of-truth sample list

```bash
# run28/samples.txt  (single source of truth)
photon5
photon10
photon20
photon10_double
jet5
jet8
jet10
jet12
jet12_double
jet15
jet20
jet30
jet40
jet50
```

```bash
# run28/cleanup_and_run.sh
#!/usr/bin/env bash
set -euo pipefail
base=/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28
while read -r s; do
    [[ -z "$s" || "$s" =~ ^# ]] && continue
    echo "--- $s ---"
    rm -rf "$base/$s/condorout/OutDir"*
    "$base/run_condor.sh" "$s"
done < "$base/samples.txt"
```

## Priority and effort

| Priority | Fix | LOC | Payoff |
|----------|-----|-----|--------|
| P0 (today) | Add `root -b -l -q` + exit-code check in CondorRunSim_new.sh | 3 | Stops interactive hangs; flags broken jobs |
| P0 (today) | Change `Log = /tmp/sli_$i.log` → `Log = $(OutDir)/job.log` | 1 | Persistent, per-job logs |
| P0 (today) | Change `cat >>ff.sub` → `cat >ff.sub` | 1 | Kills silent resubmit corruption |
| P1 (this week) | Replace hardcoded `total_line=10000` with `wc -l` (already done in jet5/jet8) | 2 | Fixes empty-output tail jobs in doubles |
| P1 (this week) | Replace `PeriodicHold` with `max_retries` + `on_exit_hold` | 3 | Clear retry policy, automatic triage |
| P2 (next week) | Consolidate to one shared `submit.sub` + `run_condor.sh` using `Queue OutDir matching` | ~100 net | Removes ~560 lines of duplication; drift-proof |
| P2 (next week) | Single `samples.txt` + rewrite cleanup/hadd_combined to read it | ~50 | Fixes missing-sample bugs in hadd_combined and cleanup_and_run |
| P3 | Add `+Sample` / `+ProductionTag` classads | 2 | Easier `condor_q` triage |

## Verification checklist after any fix

- [ ] Run one sample end-to-end (e.g. `photon5`), count `OutDir*/caloana.root` with `slimtree->GetEntries() > 0`
- [ ] Confirm held jobs release automatically (or surface `HoldReason`)
- [ ] Confirm no job writes to scheduler `/tmp`
- [ ] `hadd_combined.sh` on the same sample produces `combined.root` with entries ≈ 1000 × per-job entries
- [ ] Re-run a second time without cleanup — submit should not corrupt `ff.sub`/`submit.sub`

## Files surveyed

```
run28/
├── cleanup_and_run.sh                # hardcoded sample list (incomplete)
├── cleanup.sh                         # hardcoded sample list (diverged)
├── hadd_combined.sh                   # hardcoded sample list (only doubles)
├── {photon5,photon10,photon20,photon10_double}/
│   ├── run_condor.sh                  # 14 near-duplicates
│   ├── CondorRunSim_new.sh            # 14 identical (md5 f4f9f3d4...)
│   └── Fun4All_run_sim.C              # 11 identical + 3 (doubles + jet8) with TruthJetInput embed_flag=2
└── {jet5,jet8,jet10,jet12,jet12_double,jet15,jet20,jet30,jet40,jet50}/
    ├── ...same structure...
    └── (jet8 uses single test.list; jet5 has dynamic lines_per_job variant)
```

## Not touched by this review

- `jet50/condor_{0,plus8,minus8,time}/` — look like pedestal/time-width cross checks; not part of main production; worth archiving or documenting
- `jet12{,_double}/test_secondvtx/` — single-job test directories
- `mb/` — empty
- `jet70/` — has `Fun4All_run_sim.C` but no `run_condor.sh` / `CondorRunSim_new.sh`; looks abandoned
