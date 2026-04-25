---
description: Audit efficiencytool/results merged ROOT files for TFileMerger silent partial-merge regressions
---

Scan all merged MC ROOT files in `efficiencytool/results/` and verify their integrals
equal the sum of their per-period (`_0rad` + `_1p5mrad`) inputs across ALL TH1/TH2 keys.
Catches the TFileMerger silent partial-merge bug (some histograms drop one period while
others merge correctly, defeating single-probe validators).

Arguments: $ARGUMENTS (optional results directory override; default: `efficiencytool/results/`)

## Steps

1. Glob the candidate merged files in the results dir matching:
   - `MC_efficiency_bdt_*.root` (signal MC merge)
   - `MC_efficiency_jet_bdt_*.root` (jet MC merge)
   - `MC_response_bdt_*.root` (response merge)
   Exclude files ending in `_0rad.root`, `_1p5mrad.root`, `_0mrad.root`,
   `.broken_TFileMerger`, `.pre_hadd_backup`, or any other `.root` suffix
   already used as a snapshot.

2. For each candidate `<base>.root`, check that BOTH per-period siblings exist on disk:
   - `<base>_0rad.root`
   - `<base>_1p5mrad.root`
   If either is missing, skip silently — that variant is standalone (not from `merge_periods.sh`).

3. For each variant pair, open all three files via `uproot` and iterate ALL keys.
   For each key whose object is a TH1 or TH2 (anything with `.values()`):
   - Compute `i0 = sum(0rad.values)`, `i1 = sum(1p5mrad.values)`, `im = sum(merged.values)`.
   - Skip if `i0 + i1 == 0` (empty histogram, not a useful probe).
   - Compute `rel = abs(im - (i0+i1)) / (i0+i1)`.
   - If `rel > 0.05` (>5%), record as BROKEN: `(variant, key, i0, i1, im, rel)`.

4. Build the audit report at `reports/merge_audit_{YYYY-MM-DD}.md` with:
   - Summary line: "Audited N variants, M broken (X distinct keys affected)".
   - Per broken variant: a table of broken keys + integrals.
   - Remediation block at the bottom with the exact `hadd -f ...` commands needed:
     ```bash
     cd efficiencytool
     hadd -f results/<base>.root results/<base>_0rad.root results/<base>_1p5mrad.root
     ```
   - At the top: "**No breakages found.**" if the broken list is empty.

5. Compare against the most recent prior audit in `reports/merge_audit_*.md`:
   - List variants that are newly BROKEN this run (not in the prior audit).
   - List variants that were BROKEN before but are now OK (resolved).
   - If the prior audit doesn't exist, say so and skip the diff.

6. Print a one-line summary to the console: `Audit complete: <N> ok, <M> broken. See reports/merge_audit_<date>.md`.

7. Print the proposed `hadd` commands to the console for any broken variants — the user
   can copy-paste to remediate immediately.

## Notes

- This is a **read-only** audit (until the user explicitly runs the proposed `hadd`).
- Probe histograms whose integrals are exactly equal across all three files (degenerate
  zero-content) are skipped — they don't tell us anything.
- Tolerance 5% is the action threshold. Smaller deviations (1e-3 to 5%) are real but
  often arise from floating-point summation order differences on weighted fills; ignore
  them.
- Background context: `merge_periods.C` switched from `TFileMerger` to `gSystem->Exec("hadd")`
  on 2026-04-25 after a real-world incident in which 5 variants had ~50% of their cluster
  histograms silently lose the 0rad period. This audit is the regression check.
