---
description: Write a summary report of work done in this session or for a pipeline step
---

Generate a concise summary report of work done in this conversation session.

Arguments: $ARGUMENTS (optional topic or pipeline name, e.g., `training`, `systematics`, `efficiency bdt_nom`)

## Report structure

```markdown
## Session Report — [today's date] — [topic or "general"]

### What was done
- Bullet list of actions taken (files modified, scripts run, configs changed)

### Key results
- Numbers, plots generated, models trained, jobs submitted
- Cross-section values, systematic uncertainties if relevant

### Issues found
- Any errors, warnings, inconsistencies discovered during the session

### Next steps
- What should be done next in the pipeline
```

## Pipeline-specific reports

If a specific pipeline is named:

- **training**: Report features used, ET bins, model performance (AUC), output model files
- **systematics**: Which step completed (1=configs, 2=condor, 3=aggregation), how many variations ran, key deviation values
- **efficiency <var_type>**: Config used, output files, cross-section values per pT bin
- **plotting**: Which plots were generated, output directory, any missing inputs

## Steps

1. Review the conversation history to identify what was accomplished.

2. Write the report following the structure above.

3. Create the `reports/` directory if it doesn't exist:
   ```
   mkdir -p /sphenix/user/shuhangli/ppg12/reports
   ```

4. Save the report to `reports/session_YYYY-MM-DD_{topic}.md`.

5. If the user asks, also generate a LaTeX snippet suitable for the analysis note appendix using `\begin{itemize}` formatting and custom macros from `defs.sty`.
