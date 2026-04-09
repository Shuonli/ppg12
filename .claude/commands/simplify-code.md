---
description: Review and simplify analysis code with PPG12-specific quality checks
---

Review analysis code for simplification opportunities and PPG12-specific quality issues.

Arguments: $ARGUMENTS (file path or directory to review)

## Steps

1. Read the target file(s). If a directory is given, focus on `.C` and `.py` files.

2. For each file, check for these PPG12-specific issues:

### Hardcoded constants
- Cross-section weights (photon5cross, jet50cross, etc.) that should come from YAML config or a shared header
- Cut values (eta bounds, pT thresholds, isolation parameters) that should be config-driven
- Path strings that should use config `input.*` fields
- pT bin edges that should reference `plotcommon.h` or config

### Duplicated logic
- Repeated ABCD region code (tight_iso, tight_noniso, nontight_iso, nontight_noniso) that could use a loop or helper
- Copy-pasted histogram booking/filling across regions
- Repeated file opening patterns

### Dead code
- Commented-out blocks (large sections, not explanatory comments)
- Unused variables or includes
- Unreachable branches behind always-true/false conditions

### Missing conventions
- Output filenames without `var_type` suffix (silent overwrite risk)
- Flat BDT/isolation thresholds where parametric (ET-dependent) cuts should be used
- Missing `SaveYamlToRoot()` pattern for config archival

### Redundant I/O
- Opening the same ROOT file multiple times
- Re-reading config that was already parsed

3. For each finding, propose a concrete simplification with before/after code snippets.

4. Ask the user before applying any changes. If confirmed, make the edits.

5. Follow conventions from `.claude/rules/root-macros.md` and `.claude/rules/python-analysis.md`.
