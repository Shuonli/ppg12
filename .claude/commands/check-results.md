---
description: Validate that all systematic variation results are complete and healthy
---

Check which systematic variation results exist and flag issues.

Arguments: $ARGUMENTS (optional results directory override, default: `efficiencytool/results/`)

## Steps

1. Read `efficiencytool/make_bdt_variations.py` and extract every variant `name` from the VARIANTS list. Each should produce `Photon_final_bdt_{name}.root` and `Photon_final_bdt_{name}_mc.root`.

2. List all files matching `Photon_final_*.root` in the results directory.

3. For each expected variant name, check:
   - Does `Photon_final_bdt_{name}.root` exist? Report size and modification date.
   - Does `Photon_final_bdt_{name}_mc.root` exist?
   - Flag files smaller than 10 KB as likely failed/empty.
   - Flag files older than 30 days as potentially stale.

4. List any unexpected `Photon_final_*` files that do NOT correspond to a VARIANTS entry (orphaned results).

5. Check for stale results: any result file older than the corresponding config file means it needs re-running.

6. Print a summary table:
   ```
   | Variant | Data | MC | Size | Modified | Status |
   ```
   With status: OK, MISSING, SMALL, STALE

This is a **read-only** inspection. Do not modify any files.
