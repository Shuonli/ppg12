---
description: Compare two YAML analysis configs and explain the physics differences
---

Compare two analysis config YAML files and explain the physics differences.

Arguments: $ARGUMENTS (two config file paths or variant shortnames separated by space, e.g., `bdt_nom bdt_tightbdt50`)

## Steps

1. Parse the two arguments. If they look like shortnames (no path separator), expand to `efficiencytool/config_bdt_{name}.yaml`.

2. Read both config files.

3. Compare every key-value pair recursively. For each difference, show:
   - YAML path (e.g., `analysis.tight.bdt_min_intercept`)
   - Value in file A vs file B
   - Physical meaning:
     - `bdt_min*` / `bdt_max*`: BDT classification threshold (tight/non-tight boundary)
     - `reco_iso_max_b/s`: isolation cone energy cut (intercept + slope x ET)
     - `reco_noniso_min_shift`: non-isolated sideband lower boundary offset
     - `npb_score_cut`: non-prompt background rejection threshold
     - `cluster_escale`: energy scale factor
     - `cluster_eres`: energy resolution smearing
     - `vertex_reweight_on`: MC vertex z reweighting toggle
     - `fit_option`: purity fit function (1=default, 0=Pade)
     - `bdt_et_bin_models`: which trained BDT model per ET bin

4. Ignore `var_type` differences (expected to differ).

5. For parametric cuts, evaluate at reference ET=15 GeV to show concrete threshold values.

6. Identify which systematic group this variation belongs to by checking `SYST_TYPES` in `make_bdt_variations.py`.

7. Summarize: "This variation tests [description] and belongs to the [group] systematic group."

Format as a clean table:
```
| Parameter | Config A | Config B |
```

This is a **read-only** comparison. Do not modify any files.
