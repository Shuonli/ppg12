---
description: Run BDT training pipeline via main_training.py
---

Run BDT training pipeline.

Config: $ARGUMENTS (default: config.yaml). Append `--dry-run` to only preview what would be trained.

Working directory: `/sphenix/user/shuhangli/ppg12/FunWithxgboost`

## Steps

1. Read the config file and report:
   - Number of features (from `data.features`)
   - ET bin edges (from `binning.pt_edges`)
   - Number of models that will be trained
   - Output directory

2. Verify training data files exist (`NOSPLIT/shapes_*.txt` or paths in config).

3. If `$ARGUMENTS` contains `--dry-run`: stop here and report what would be trained. Do not run.

4. Run:
   ```
   source /sphenix/u/shuhang98/setup.sh && cd /sphenix/user/shuhangli/ppg12/FunWithxgboost && python main_training.py --config {config}
   ```

5. After completion:
   - List new model files in `binned_models/`
   - Report training metrics if printed to stdout (AUC, accuracy)
   - Remind user to run `/apply-bdt all` next to apply the new model to samples
