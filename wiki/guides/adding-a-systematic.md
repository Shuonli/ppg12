# Adding a New Systematic Variation

## Overview

Adding a systematic requires three things:
1. Define the variation in `make_bdt_variations.py`
2. Run the efficiency pipeline for the new config
3. Add the variation to the aggregation taxonomy

## Step 1: Define the Variation

Edit `efficiencytool/make_bdt_variations.py`.

### 1a. Add to VARIANTS list

`VARIANTS` is a **list of dicts**. Each dict has a `name` key, override keys as flat top-level entries (not nested under an `"overrides"` key), plus `syst_type` and `syst_role`:

```python
VARIANTS = [
    # ... existing variants ...

    # Two-sided example (up/down pair):
    dict(name="my_var_up",   my_config_key=new_value_up,
         syst_type="my_type", syst_role="up"),
    dict(name="my_var_down", my_config_key=new_value_down,
         syst_type="my_type", syst_role="down"),

    # One-sided example (symmetric):
    dict(name="my_var_alt",  my_config_key=alt_value,
         syst_type="my_type", syst_role="one_sided"),

    # Cross-check only (not a systematic):
    dict(name="my_crosscheck", my_config_key=check_value,
         syst_type=None, syst_role=None),
]
```

The variant `name` becomes the suffix in `config_bdt_{name}.yaml` and `var_type` is set to `bdt_{name}`.

### 1b. Add flat key to OVERRIDE_MAP (if needed)

If your override key is not already in the OVERRIDE_MAP, add it. Each entry maps a flat key to a tuple of `(yaml_section_path_list, yaml_leaf_key)`:

```python
OVERRIDE_MAP = {
    # ... existing mappings ...
    "my_config_key": (["analysis"],          "my_config_key"),
    # Or for a nested section:
    "my_nested_key": (["analysis", "tight"], "my_leaf_key"),
}
```

Existing mappings include keys like `tight_bdt_min_intercept` -> `(["analysis", "tight"], "bdt_min_intercept")`, `npb_score_cut` -> `(["analysis", "common"], "npb_score_cut")`, `clusterescale` -> `(["analysis"], "cluster_escale")`, etc.

### 1c. Add to SYST_TYPES

Each syst_type has a `mode` and a `group` assignment:

```python
SYST_TYPES = {
    # ... existing types ...
    "my_type": {"mode": "two_sided", "group": "my_group"},
}
```

Modes:
- `two_sided`: expects variants with roles `up` and `down`; takes each direction separately
- `one_sided`: expects role `one_sided`; symmetric `low = high = |delta|`
- `max`: expects role `max`; bin-wise maximum `|delta|` over all variants of this type
- `placeholder`: reserved for types with no active variants yet

The `group` field determines which SYST_GROUPS entry this type contributes to.

### 1d. Add to SYST_GROUPS (if new group)

```python
SYST_GROUPS = {
    # ... existing groups ...
    "my_group": ["my_type"],
}
```

Or append to an existing group's type list.

### 1e. Add group to FINAL_SYSTS (if new group)

```python
FINAL_SYSTS = ["purity", "eff", "escale", "eres", "mbd", "unfolding", "my_group"]
```

## Step 2: Generate and Run

```bash
cd /sphenix/user/shuhangli/ppg12/efficiencytool

# Regenerate all variation configs
python make_bdt_variations.py config_bdt_nom.yaml

# Verify the new config was created
ls config_bdt_my_var_up.yaml config_bdt_my_var_down.yaml

# Run the new variations
bash oneforall_tree.sh config_bdt_my_var_up.yaml
bash oneforall.sh config_bdt_my_var_up.yaml
bash oneforall_tree.sh config_bdt_my_var_down.yaml
bash oneforall.sh config_bdt_my_var_down.yaml

# Or run all variations on condor
condor_submit oneforall.sub
```

## Step 3: Aggregate

```bash
cd /sphenix/user/shuhangli/ppg12/plotting
python calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures
```

This will automatically pick up the new type and group from the imported `make_bdt_variations.py` metadata.

## Step 4: Verify

Check the output:
```bash
ls rootFiles/syst_bdt_my_type.root       # Per-type result
ls rootFiles/syst_bdt_my_group.root      # Per-group result
ls figures/syst_bdt_rel_my_type.pdf      # Visualization
```

## Example: Adding an Isolation Cone Size Systematic

```python
# In make_bdt_variations.py VARIANTS list:
dict(name="cone_r03",        use_topo_iso=0,
     syst_type="cone", syst_role="max"),
dict(name="cone_r04_tower",  use_topo_iso=0, cone_size=4,
     syst_type="cone", syst_role="max"),

# OVERRIDE_MAP (use_topo_iso already exists; add cone_size if needed):
"cone_size": (["analysis"], "cone_size"),

# SYST_TYPES:
"cone": {"mode": "max", "group": "purity"},

# SYST_GROUPS -- add to existing group:
SYST_GROUPS["purity"].append("cone")
```

## Checklist

- [ ] Variant dict(s) added to `VARIANTS` list with unique `name`
- [ ] Flat override key mapped in `OVERRIDE_MAP` (if new key)
- [ ] `syst_type` added to `SYST_TYPES` with correct `mode` and `group`
- [ ] Type's group exists in `SYST_GROUPS` (or create new group)
- [ ] Group added to `FINAL_SYSTS` (if new group)
- [ ] Validation passes: `python -c "from make_bdt_variations import validate_variants"`
- [ ] Config generated: `python make_bdt_variations.py config_bdt_nom.yaml`
- [ ] Pipeline run: `bash oneforall_tree.sh` + `bash oneforall.sh`
- [ ] Aggregation run: `python calc_syst_bdt.py`
- [ ] Output verified: per-type ROOT file, per-group ROOT file, visualization

## Important Notes

- `VARIANTS` is a **list** of dicts, not a dict of dicts -- order determines generation order
- Override keys go directly in the dict (e.g., `dict(name="x", my_key=val)`), not nested under `"overrides"`
- Uses `ruamel.yaml` (not PyYAML) to preserve comments and field ordering
- The `var_type` in each config is set automatically to `bdt_{name}` -- must be unique
- `syst_type` and `syst_role` must both be set or both be None
- Valid roles per mode: `two_sided` -> up/down, `one_sided` -> one_sided, `max` -> max
- `validate_variants()` runs automatically at import time and will raise errors for invalid metadata
- For cross-checks (not systematic uncertainties), set `syst_type=None, syst_role=None`
