# Systematic Variations

## Overview

Systematic uncertainties are evaluated by re-running the full efficiency/yield pipeline with modified analysis configs. Each variation changes one or two config parameters from the nominal. Results are aggregated into groups via quadrature sums.

## Three-Step Pipeline

```bash
# Step 1: Generate 30+ config variants from nominal
cd efficiencytool
python make_bdt_variations.py config_bdt_nom.yaml

# Step 2: Run all variants in parallel on HTCondor
condor_submit oneforall.sub

# Step 3: Aggregate systematics
cd ../plotting
python calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures
```

## Taxonomy

Each variation has a `syst_type` (what physics it probes) and a `syst_role` (`up`/`down`/`one_sided`/`max`/`None`).

### Variation Modes

| Mode | How delta is computed |
|------|----------------------|
| `two_sided` | Separate up/down variants; mirrors if one direction is missing |
| `one_sided` | Symmetric: `low = high = |delta|` |
| `max` | Bin-wise maximum of `|delta|` over all `max`-role variants |
| `placeholder` | Reserved type with no active variants yet |

## Systematic Types (SYST_TYPES)

From `efficiencytool/make_bdt_variations.py` lines 148-164:

| syst_type | Mode | Group | Physics |
|-----------|------|-------|---------|
| `noniso` | two_sided | purity | Non-iso sideband boundary |
| `nt_bdt` | one_sided | purity | Non-tight BDT lower boundary |
| `purity_fit` | one_sided | purity | Purity fit function choice |
| `mc_purity_correction` | one_sided | purity | MC-driven purity correction |
| `tight_bdt` | two_sided | eff | Tight BDT threshold |
| `npb_cut` | two_sided | eff | NPB score cut |
| `vtx_reweight` | one_sided | eff | Vertex reweighting toggle |
| `bdt_model` | max | eff | BDT model variant |
| `b2bjet` | one_sided | eff | Back-to-back jet cut |
| `timing` | two_sided | eff | Cluster timing window |

> **Note:** `vtx_reweight`, `bdt_model`, `b2bjet`, and `timing` declare `group="eff"` in SYST_TYPES metadata, but their variants currently have `syst_type=None` (cross-checks only). They are NOT included in `SYST_GROUPS["eff"]`, which only aggregates `tight_bdt` and `npb_cut`. If promoted to real systematics, add them to SYST_GROUPS.
| `escale` | two_sided | escale | Calorimeter energy scale |
| `eres` | max | eres | Energy resolution smearing |
| `reweight` | one_sided | unfold | Unfolding response reweighting |
| `mbd` | placeholder | mbd | MBD efficiency (no active variants) |
| `nor` | placeholder | unfolding | Normalization (no active variants) |

## Systematic Groups (SYST_GROUPS)

From `make_bdt_variations.py` lines 167-174:

| Group | Member types | Physics |
|-------|-------------|---------|
| **purity** | noniso, nt_bdt, purity_fit, mc_purity_correction | Purity extraction |
| **eff** | tight_bdt, npb_cut | Identification efficiency |
| **escale** | escale | Calorimeter energy scale |
| **eres** | eres | Energy resolution |
| **mbd** | mbd | MBD trigger efficiency (placeholder) |
| **unfolding** | reweight | Unfolding prior |

### Aggregation

- **Per-group:** Quadrature sum of all member type deltas, bin-by-bin.
- **Total:** Quadrature sum of all `FINAL_SYSTS` groups (`["purity", "eff", "escale", "eres", "mbd", "unfolding"]`) + flat luminosity uncertainty.
- **Luminosity:** Asymmetric fractional uncertainty: down = (25.2 - 23.5) / 25.2, up = (27.5 - 25.2) / 25.2.

## All Variant Definitions

From `make_bdt_variations.py` VARIANTS list (lines 34-141). Variants with `syst_type=None` are cross-checks, not systematics.

### Purity-Related Variants

| Variant name | Config override | syst_type | syst_role |
|--------------|----------------|-----------|-----------|
| `noniso04` | `reco_noniso_min_shift=0.1` | noniso | down |
| `noniso10` | `reco_noniso_min_shift=1.0` | noniso | up |
| `ntbdtmin05` | `nt_bdt_min=0.05` | nt_bdt | one_sided |
| `ntbdtmin10` | `nt_bdt_min=0.10` | nt_bdt | one_sided |
| `purity_pade` | `fit_option=0` | purity_fit | one_sided |

> **Naming note:** `purity_pade` is misleadingly named — it sets `fit_option=0` which selects **erf** (not Pade). The nominal config uses `fit_option=1` (Pade), so this variant tests switching FROM Pade TO erf.
| `mc_purity_correction` | `mc_purity_correction=1` | mc_purity_correction | one_sided |

### Efficiency-Related Variants

| Variant name | Config override | syst_type | syst_role |
|--------------|----------------|-----------|-----------|
| `tightbdt50` | `tight_bdt_min_intercept=0.70, nt_bdt_max_intercept=0.70` | tight_bdt | down |
| `tightbdt70` | `tight_bdt_min_intercept=0.90, nt_bdt_max_intercept=0.90` | tight_bdt | up |
| `npb03` | `npb_score_cut=0.3` | npb_cut | down |
| `npb07` | `npb_score_cut=0.7` | npb_cut | up |
| `vtxreweight0` | `vertex_reweight_on=0` | vtx_reweight | one_sided |
| `etbin_v3E_v3E` | `bdt_et_bin_models=["base_v3E","base_v3E"]` | bdt_model | max |

### Energy Scale Variants

| Variant name | Config override | syst_type | syst_role |
|--------------|----------------|-----------|-----------|
| `energyscale26up` | `clusterescale=1.026` | escale | down |
| `energyscale26down` | `clusterescale=0.974` | escale | up |

Note: the up/down naming refers to the scale factor direction, but the syst_role is inverted (higher scale -> smaller cross-section -> labeled "down" in the systematics).

### Energy Resolution Variants

| Variant name | Config override | syst_type | syst_role |
|--------------|----------------|-----------|-----------|
| `energyresolution5` | `clustereres=0.05` | eres | max |

### Unfolding Variants

| Variant name | Config override | syst_type | syst_role |
|--------------|----------------|-----------|-----------|
| `no_unfolding_reweighting` | `reweight=0` | reweight | one_sided |

### Cross-Check Variants (syst_type=None)

| Variant name | Config override | Purpose |
|--------------|----------------|---------|
| `nom` | (none) | Sanity check: identical to nominal |
| `0rad` | `run_min=47289, run_max=51274, lumi=32.6574` | 0 mrad crossing angle data |
| `all` | `run_min=47289, run_max=54000, lumi=49.562` | All runs combined |
| `mbdrms01` / `mbdrms03` / `mbdrms05` | `mbd_avg_sigma_max=0.1/0.3/0.5` | MBD pileup cut scans |
| `tightbdtflat06` | `tight_bdt_min_intercept=0.6, slope=0` | Flat BDT threshold |
| `mciso_noscaleshift` | `mc_iso_scale=1.0, mc_iso_shift=0.0` | No MC iso correction |
| `mciso_no_shift` | `mc_iso_shift=0.0` | MC iso shift only |
| `mciso_no_scale` | `mc_iso_scale=1.0` | MC iso scale only |
| `common_wr03` / `common_wr05` | `common_wr_cogx_bound=0.3/0.5` | Radial width preselection |
| `split` | Different cluster node / file paths | Split cluster node check |
| `etbin_E_E` | `bdt_et_bin_models=["base_E","base_E"]` | BDT model cross-check |
| `etbin_v1E_v3E` | `bdt_et_bin_models=["base_v1E","base_v3E"]` | Mixed BDT models |
| `b2bjet` | `common_b2bjet_cut=1` | Back-to-back jet veto |
| `timingcut_2` / `timingcut_5` | `cluster_mbd_time_min/max=+/-2/5` | Timing cut scan |
| `energyresolution7` / `energyresolution8` | `clustereres=0.07/0.08` | Extra resolution smearing |

## Config Generation Details

`VARIANTS` is a **list of dicts** (not a dict of dicts). Each dict has a `name` key plus flat override keys. The `OVERRIDE_MAP` translates flat keys to YAML paths as `(section_path_list, leaf_key)` tuples:

```python
OVERRIDE_MAP = {
    "tight_bdt_min_intercept": (["analysis", "tight"], "bdt_min_intercept"),
    "nt_bdt_max":              (["analysis", "non_tight"], "bdt_max"),
    "npb_score_cut":           (["analysis", "common"], "npb_score_cut"),
    "clusterescale":           (["analysis"], "cluster_escale"),
    "reweight":                (["analysis", "unfold"], "reweight"),
    # ... etc
}
```

Uses `ruamel.yaml` to preserve comments and field ordering. Output: `config_bdt_{name}.yaml` with `var_type` automatically set to `bdt_{name}`.

## Adding a New Systematic

See [Adding a Systematic](../guides/adding-a-systematic.md) for step-by-step instructions.

## Output Files

### Per-type

`plotting/rootFiles/syst_bdt_{type}.root` with histograms:
- `h_dev_low`, `h_dev_high` -- absolute deviations
- `h_dev_rel_low`, `h_dev_rel_high` -- relative deviations

### Per-group

`plotting/rootFiles/syst_bdt_{group}.root` (same histogram structure)

### Total

`plotting/rootFiles/syst_bdt_total.root` and `syst_sum.root`

The `syst_sum.root` file (with `h_sum_low/high`, `h_sum_rel_low/high`) is consumed by `plot_final_selection.C` for the systematic band on the final plot.

## Visualization

- `figures/syst_bdt_rel_{type}.pdf` -- per-type relative deviation (up=red, down=blue)
- `figures/syst_bdt_breakdown.pdf` -- all groups + total on one canvas
- Group colors: purity (Violet+1), eff (Azure+7), escale (Green-2), eres (Orange+1), mbd (Violet-1), unfolding (Magenta+2)
