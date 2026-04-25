---
description: Fast version of /ask-shuonli for simple lookup-style PPG12 tasks — single-pass, no wave pipeline
---

Fast variant of `/ask-shuonli` for lightweight tasks. Skip the wave pipeline, skip prompt refinement, single-pass execution. Escalate to `/ask-shuonli` (full) the moment the task grows beyond a quick lookup.

Task: $ARGUMENTS

## When to use FAST vs FULL

**FAST is appropriate for:**
- Single-file lookups ("what's in this ROOT file", "where is X defined", "what does this config set")
- Simple status checks ("did the jobs finish", "is this branch up to date", "is config Y consistent")
- Single config validation
- Quick Q&A about pipeline mechanics, paths, conventions
- "Show me the cross section weights" / "list the syst variations" type queries
- Inspecting an existing plot PDF (no regeneration)

**STOP and escalate to `/ask-shuonli` (full) if any of these apply:**
- The task compares multiple variations or sweeps
- Numbers will be cited in a report, paper, or plot caption
- A plot needs to be generated for collaborator/publication use
- The task touches more than 2 pipeline stages
- A code change is required (FAST may diagnose; FULL fixes via `code-writer` + `fix-validator`)
- Any cross-section, efficiency, purity, or systematic value will be reported externally

When escalating, tell the user one line:
> This needs the full pipeline — re-run as `/ask-shuonli <task>`.

## Principles (same as full, lighter execution)

1. **Physics first** — frame in physics terms; skip tutorial-level explanations
2. **Wiki only when relevant** — check `wiki/index.md` only if the task actually needs domain context
3. **Constants sync awareness** — flag risk if you spot a sync issue, even in passing
4. **Data-MC interpretation** — cut1, compare to inclusive MC, never "between signal and background"
5. **One pass, ≤1 supporting agent** — main agent does the work directly. At most ONE supporting agent (e.g., one `data-explorer` for tree inspection, or one `plot-cosmetics-reviewer` for an existing plot).

## Step 1 — Quick Classify

Read the request:
- Fits FAST? → proceed to Step 2
- Doesn't fit FAST? → STOP, recommend `/ask-shuonli`
- Ambiguous? → ask the user ONE clarifying question (don't guess, don't run prompt-refiner)

## Step 2 — Execute (single pass)

Do the work directly in the main agent. Use at most ONE supporting agent. Allowed agents in FAST mode:
- `data-explorer` — for one ROOT/config inspection task
- `plot-cosmetics-reviewer` — for one PDF cosmetic check
- `physics-reviewer` — for a one-file read-only review

Anything more than one of these → escalate to FULL.

If you cite a number from a ROOT file: read the ROOT file directly (uproot or `/inspect-root`) inline. Cite the file + histogram/branch + bin in your answer. No separate `numerical-rederivation` agent needed in FAST mode — but the source must be visible.

If you generate a plot in FAST mode (rare — usually you should escalate): launch `plot-cosmetics-reviewer` on the output PDF before presenting. This is non-negotiable.

## Step 3 — Present (concise)

2–4 sentences for the answer, plus one of:
- A direct value with file:line or `file.root:histogram[bin]` reference
- A short list (≤5 items)
- A one-line recommendation

No findings table, no reviewer pass matrix, no LaTeX report. If the user wants a structured report, suggest re-running as `/ask-shuonli`.

## Plotting Style (still required if a plot is involved)

Same standards as `/ask-shuonli`: `plotcommon.h`, `init_plot()`, shared frames/legends, pT bins `{8,10,12,14,16,18,20,22,24,26,28,32,36}`, save to `plotting/figures/` or `PPG12-analysis-note/Figures/`. `plot-cosmetics-reviewer` still gates any plot output.

## Domain Knowledge Quick Reference (same as full skill)

- pT: 12 bins, 8–36 GeV (`plotcommon.h:ptRanges[]`)
- Eta: `[-0.7, 0.7]`
- BDT/iso: parametric, never flat
- ABCD: A/B/C/D = tight·iso / tight·noniso / nontight·iso / nontight·noniso
- Lumi: 16.6 pb⁻¹ (1.5 mrad), 32.7 pb⁻¹ (0 mrad), 49.6 pb⁻¹ (all)
- var_type unique per config; Cluster node `CLUSTERINFO_CEMC`; seed 42

## Anti-Patterns

- Don't run the wave pipeline — that's `/ask-shuonli`
- Don't launch more than one supporting agent — escalate instead
- Don't run `prompt-refiner` — if ambiguous, ask the user ONE question
- Don't write LaTeX reports or modify the analysis note
- Don't deploy to GitHub Pages
- Don't claim a fix works — FAST mode does NOT run `fix-validator`; diagnose only, fixes belong in FULL
- Don't skip `plot-cosmetics-reviewer` even in FAST mode if you produce or inspect a plot
- Don't cite a number without showing the ROOT file + histogram source
- Don't try to handle a complex task by squeezing it into FAST — escalate cleanly
