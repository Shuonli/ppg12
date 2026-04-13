# Variant Training Pipeline

Trains all 11 photon-ID BDT variants and the NPB score model on either the
`_split` or `_nosplit` input files, in parallel on Condor. Outputs are loaded
by `apply_BDT.C` according to `use_split_bdt` in `config_nom.yaml`:

- `use_split_bdt: 1` -> `model_{name}_split_single_tmva.root` + `npb_score_split_tmva.root`
- `use_split_bdt: 0` -> `model_{name}_nosplit_single_tmva.root` + `npb_score_nosplit_tmva.root`

The cluster-node variants correspond to sPHENIX's `CLUSTERINFO_CEMC`
(splitting ON, nominal) and `CLUSTERINFO_CEMC_NO_SPLIT` (splitting OFF).

---

## Prerequisites

The following txt files must exist in the working directory (produced by
`BDTinput.C` in the appropriate configuration):

| Role | Split files | Nosplit files |
|------|-------------|---------------|
| Signal (photon-ID) | `shapes_split_photon{5,10,20}.txt` | `shapes_nosplit_photon{5,10,20}.txt` |
| Background (photon-ID) | `shapes_split_jet{5,12,20,30,40}.txt` | `shapes_nosplit_jet{5,12,20,30,40}.txt` |
| Signal (NPB) | same jet files as above | same jet files as above |
| Background (NPB) | `shapes_data_npb.txt` (shared — data has no split/nosplit distinction) |

The `logs/` directory must also exist:
```bash
mkdir -p logs
```

---

## Running split training

### Pipeline A: Photon-ID BDT (11 models, Condor)

```bash
# 1. Generate per-model configs into variant_configs/
python make_variant_configs.py --variant split

# 2. Submit (VARIANT=split is the default)
condor_submit submit_variant_training.sub

# 3. Monitor
condor_q

# 4. Verify outputs
ls binned_models/model_*_split_single_tmva.root
```

### Pipeline B: NPB Score Model (single job, local)

```bash
python train_npb_score.py --config config_npb_training_split.yaml
# -> npb_models/npb_score_split_tmva.root
```

---

## Running no-split training

### Pipeline A: Photon-ID BDT (11 models, Condor)

```bash
# 1. Generate per-model configs into variant_configs/
python make_variant_configs.py --variant nosplit

# 2. Submit with the VARIANT macro override
condor_submit submit_variant_training.sub -append "VARIANT=nosplit"

# 3. Monitor
condor_q

# 4. Verify outputs
ls binned_models/model_*_nosplit_single_tmva.root
```

### Pipeline B: NPB Score Model

```bash
python train_npb_score.py --config config_npb_training_nosplit.yaml
# -> npb_models/npb_score_nosplit_tmva.root
```

---

## Applying the trained models

In `config_nom.yaml`, set:
```yaml
use_split_bdt: 1   # or 0 for nosplit
```
`apply_BDT.C` will then load the matching `model_{name}_{variant}_single_tmva.root`
and `npb_score_{variant}_tmva.root` automatically.

---

## Troubleshooting

**Check a failed Condor job:**
```bash
cat logs/variant_train_<variant>_<Cluster>.<Process>.err
```

**Re-run a single model locally:**
```bash
python train_variant_model.py variant_configs/config_base_v3E_split.yaml
python train_variant_model.py variant_configs/config_base_v3E_nosplit.yaml
```

**Background file has 0 entries** — the txt file exists but is empty or missing.
Check that `BDTinput.C` was run for that jet pT sample with the correct
`cluster_node_name`.

**Wrong VARIANT picked up by Condor** — you forgot to run
`make_variant_configs.py --variant nosplit` before submitting with
`VARIANT=nosplit`, so `variant_configs/filelist_nosplit.txt` does not exist.

---

## File Reference

| File | Description |
|------|-------------|
| `make_variant_configs.py` | Generates `variant_configs/` YAMLs and `filelist_{variant}.txt`; takes `--variant {split,nosplit}` |
| `train_variant_model.py` | CLI wrapper: runs `BinnedTrainingPipeline.run_training()` + `save_results()` on one config |
| `run_variant_training.sh` | Condor executable — sources env, runs `train_variant_model.py` |
| `submit_variant_training.sub` | Condor submission file — defaults to `VARIANT=split`, override with `-append "VARIANT=nosplit"` |
| `config_npb_training_split.yaml` | NPB score training config using split jet files |
| `config_npb_training_nosplit.yaml` | NPB score training config using nosplit jet files |
| `variant_configs/` | Generated per-model configs and `filelist_{variant}.txt` |
| `binned_models/` | Photon-ID TMVA ROOT outputs (`model_{name}_{variant}_single_tmva.root`) |
| `npb_models/` | NPB score TMVA ROOT outputs (`npb_score_{variant}_tmva.root`) |
| `logs/` | Condor stdout/stderr (`variant_train_{variant}_<Cluster>.<Process>.out/err`) |
