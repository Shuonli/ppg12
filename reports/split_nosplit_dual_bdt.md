# Split / No-Split Dual BDT Refactor

Date: 2026-04-13

## Summary

Refactored the BDT training and application pipeline to support **both** `CLUSTERINFO_CEMC` (split) and `CLUSTERINFO_CEMC_NO_SPLIT` (no-split) cluster nodes as first-class citizens. Before this change, the pipeline was hardcoded to the split variant via fixed filenames (`shapes_split_*.txt`, `bdt_split.root`, `model_*_split_*.root`). Now:

1. `BDTinput.C` emits `shapes_{split,nosplit}_{filetype}.txt`, picking the tag from `cluster_node_name` in the config.
2. Training infrastructure is variant-parametric via `make_variant_configs.py --variant {split,nosplit}`.
3. `apply_BDT.C` iterates over a list of `(cluster_node, model_suffix)` pairs from config and writes BOTH sets of BDT score branches to a single output file.
4. Efficiencytool already parametric on `cluster_node_name` — zero code change downstream.

Backward compatibility: every existing `config_bdt_*.yaml` keeps working because `apply_BDT.C` synthesizes a single-variant fallback when `cluster_variants` is absent.

## Files changed

### `FunWithxgboost/BDTinput.C`
- Added `cluster_tag` derivation from `clusternodename` (lines 174–184). Unknown cluster nodes fail loudly.
- Line 487: output filename now `"shapes_" + cluster_tag + "_" + filetype + ".txt"`.
- Line 491: data NPB output is `"shapes_" + cluster_tag + "_data_npb.txt"`. NPB *tagging* is cluster-variant agnostic, but the shower-shape features written per tagged cluster come from the configured cluster node, so the file is per-variant.

### `FunWithxgboost/apply_BDT.C`
- Full refactor (492 → 608 lines).
- New `VariantState` struct holds per-variant input pointers, output arrays, RBDT lists, NPB model. Instances are heap-allocated via `std::unique_ptr<VariantState>` stored in `std::vector`.
- New `parse_variants(config)` reads `input.cluster_variants` list; falls back to single-variant synthesis from legacy `cluster_node_name` + `use_split_bdt` fields.
- Single event loop fills both variants' output branches before calling `outtree->Fill()` exactly once per event.
- Output filenames kept as `{filetype}/bdt_split.root` and `*_with_bdt_split.root` for backward compatibility — TODO comments added.
- All 11 photon-ID `x_list` feature vectors and the 25-entry `x_npb` vector are byte-identical to the pre-refactor version — the training/apply contract is preserved.

### Training infrastructure (renames via `git mv`)
- `make_split_configs.py` → `make_variant_configs.py` — adds `--variant {split,nosplit}` argparse flag, writes to unified `variant_configs/config_{name}_{variant}.yaml`, `root_file_prefix = model_{name}_{variant}`.
- `train_split_model.py` → `train_variant_model.py` — variant-agnostic, docstring updates only.
- `run_split_training.sh` → `run_variant_training.sh` — updated executable reference.
- `run_split_training_local.sh` → `run_variant_training_local.sh` — now takes `{split|nosplit}` CLI arg, reads per-variant filelist.
- `submit_split_training.sub` → `submit_variant_training.sub` — default `VARIANT=split`, override via `condor_submit -append "VARIANT=nosplit"`.
- `README_split_training.md` → `README_variant_training.md` — rewritten for both workflows.

### New file
- `FunWithxgboost/config_npb_training_nosplit.yaml` — parallel to `_split.yaml`, points at `shapes_nosplit_jet{5,12,20,30,40}.txt` + `shapes_nosplit_data_npb.txt`, output model `npb_score_nosplit`. The split config was also updated so its NPB data file is now `shapes_split_data_npb.txt`.

### Wiki docs
- `wiki/reference/file-inventory.md` — updated FunWithxgboost file list.
- `wiki/pipeline/02-bdt-training.md` — updated training orchestration references and BDTinput output-file format description.

## How to run

### 1. Generate training input files (both variants)

```bash
cd /sphenix/user/shuhangli/ppg12/FunWithxgboost

# Split variant (current nominal — config_nom.yaml already points to CLUSTERINFO_CEMC)
for s in photon5 photon10 photon20 jet5 jet12 jet20 jet30 jet40; do
    root -l -b -q 'BDTinput.C("config_nom.yaml", "'$s'")'
done
root -l -b -q 'BDTinput.C("config_nom.yaml", "data")'   # NPB data

# No-split variant — need a config_nom_nosplit.yaml with cluster_node_name: CLUSTERINFO_CEMC_NO_SPLIT
# (one-line edit from config_nom.yaml)
for s in photon5 photon10 photon20 jet5 jet12 jet20 jet30 jet40; do
    root -l -b -q 'BDTinput.C("config_nom_nosplit.yaml", "'$s'")'
done
root -l -b -q 'BDTinput.C("config_nom_nosplit.yaml", "data")'
```

Output:
- Split: `shapes_split_{sample}.txt`, `shapes_split_data_npb.txt`
- No-split: `shapes_nosplit_{sample}.txt`, `shapes_nosplit_data_npb.txt`

### 2. Train both variants

```bash
# Generate per-model configs
python make_variant_configs.py --variant split
python make_variant_configs.py --variant nosplit

# Submit 11 photon-ID models per variant via Condor
mkdir -p logs
condor_submit submit_variant_training.sub                              # split (default)
condor_submit submit_variant_training.sub -append "VARIANT=nosplit"    # nosplit

# Train NPB scores
python train_npb_score.py --config config_npb_training_split.yaml
python train_npb_score.py --config config_npb_training_nosplit.yaml
```

Outputs:
- `binned_models/model_{name}_{variant}_single_tmva.root` (11 × 2 = 22 files)
- `npb_models/npb_score_{variant}_tmva.root` (2 files)

### 3. Apply both BDTs in one pass

Add to `FunWithxgboost/config_nom.yaml` (the apply_BDT config):

```yaml
input:
  cluster_variants:
    - node: CLUSTERINFO_CEMC
      model_suffix: "_split"
    - node: CLUSTERINFO_CEMC_NO_SPLIT
      model_suffix: "_nosplit"
```

Then:
```bash
for s in photon5 photon10 photon20 jet5 jet12 jet20 jet30 jet40; do
    root -l -b -q 'apply_BDT.C("config_nom.yaml", "'$s'")'
done
# + data loop over part_*.root
```

Output: `FunWithxgboost/{sample}/bdt_split.root` now contains both `cluster_bdt_CLUSTERINFO_CEMC_*` AND `cluster_bdt_CLUSTERINFO_CEMC_NO_SPLIT_*` branches (plus matching NPB score branches).

### 4. Run efficiency on either variant

`efficiencytool/config_bdt_nom.yaml` (split, existing) — works unchanged.

For no-split, create `efficiencytool/config_bdt_nosplit.yaml`:
```yaml
input:
  cluster_node_name: "CLUSTERINFO_CEMC_NO_SPLIT"
  # ... everything else same as config_bdt_nom.yaml ...
output:
  var_type: "bdt_nosplit"
```
Then `bash oneforall.sh config_bdt_nosplit.yaml`.

**Caveat:** `reco_iso_max_b/s`, `bdt_min_intercept/slope`, and `mc_iso_shift/scale` in the efficiencytool config were fit against split-cluster distributions. They will need re-calibration for the no-split variant before the cross-section numbers are physically meaningful.

## Review findings (Wave 2)

| Reviewer | Verdict |
|---|---|
| cross-file consistency | 6 PASS, 1 WARN, 1 FAIL (fixed) |
| backward compat | 6 PASS, 2 WARN, 0 FAIL |
| physics correctness | 14 PASS, 1 WARN, 0 FAIL |

### FAIL (fixed)
- Filename/config mismatch: `BDTinput.C` writes `shapes_{tag}_data_npb.txt` but the NPB configs originally referenced `shapes_data_npb.txt` (unprefixed). **Fix:** updated `config_npb_training_split.yaml` to `shapes_split_data_npb.txt` and `config_npb_training_nosplit.yaml` to `shapes_nosplit_data_npb.txt`. The prefixed naming is correct because the NPB data file contains shower-shape features extracted from the configured cluster node.

### WARNs
- `submit_variant_training.sub` requires `FunWithxgboost/logs/` to pre-exist — `mkdir -p logs` before first submit.
- `apply_BDT.C` output filename kept as `bdt_split.root` for backward compat — may be confusing when only no-split BDTs are written. Rename deferred.
- NPB no-split config reuses split hyperparameters verbatim — may benefit from a re-tune if no-split cluster distributions are meaningfully different, but not required for a first pass.

## Smoke tests passed

- `BDTinput.C` loads in Cling.
- `apply_BDT.C` loads in Cling (608 lines).
- `make_variant_configs.py` and `train_variant_model.py` parse as valid Python.
- Both NPB configs parse as valid YAML with matching 25-feature lists, correct sample paths, variant-specific `shapes_{split,nosplit}_data_npb.txt` backgrounds.
- All 11 photon-ID feature vectors verified byte-identical between `MODEL_FEATURES` (py) and `x_list` (C++).
- Legacy `config_bdt_nom.yaml` (split, `cluster_node_name: CLUSTERINFO_CEMC`) → `parse_variants` produces `[(CLUSTERINFO_CEMC, "")]` → reads `model_{name}_single_tmva.root` → no regression.

## Open items for the user

1. Decide whether to delete the empty `FunWithxgboost/split_configs/` directory (only contains untracked `filelist.txt`).
2. Create `FunWithxgboost/config_nom_nosplit.yaml` (one-line diff from `config_nom.yaml`) before running BDTinput.C for the no-split variant.
3. Add the `cluster_variants` block to `config_nom.yaml` before running `apply_BDT.C` in dual-variant mode.
4. Re-calibrate iso/BDT cut parameters in the no-split efficiencytool config once the no-split BDT-scored trees exist.
5. Stage and commit the changes with a descriptive message.

## Uncommitted changes

22 files touched:
- 2 modified (`BDTinput.C`, `apply_BDT.C`)
- 5 renamed via `git mv` (make/train/run/submit/README)
- 11 deleted (`split_configs/config_*_split.yaml` — stale)
- 1 created (`config_npb_training_nosplit.yaml`)
- 2 modified wiki files
- 1 new report (this file)
