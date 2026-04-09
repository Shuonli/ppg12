---
description: Check HTCondor job status and cross-reference with expected results
---

Check HTCondor job status and result completeness for the efficiency pipeline.

Arguments: $ARGUMENTS
- (empty) → summary view
- `detailed` → also inspect log files for errors

## Steps

1. Run `condor_q -submitter shuhang98` and summarize: total jobs, running, idle, held. If held jobs exist, report their hold reasons.

2. List all `config_bdt_*.yaml` files in `/sphenix/user/shuhangli/ppg12/efficiencytool/`. For each, extract the `var_type` and check if `results/Photon_final_{var_type}.root` AND `results/Photon_final_{var_type}_mc.root` exist.

3. Categorize each config as: COMPLETE, MISSING_DATA, MISSING_MC, MISSING_BOTH.

4. For any MISSING results, check if a condor log exists in `logs/`:
   - If the log exists, tail the last 20 lines to identify errors.
   - Report likely failure reasons (segfault, file not found, timeout).

5. Print a summary table:
   ```
   | Status       | Count | Configs |
   |-------------|-------|---------|
   | COMPLETE    |   N   | ...     |
   | MISSING     |   N   | ...     |
   ```

6. If `$ARGUMENTS` contains `detailed`: also show the last 5 lines of every non-empty `.err` file in `logs/`.

This is a **read-only** inspection. Do not modify any files.
