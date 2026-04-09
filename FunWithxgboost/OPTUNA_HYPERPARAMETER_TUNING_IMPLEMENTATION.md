# Optuna Hyperparameter Tuning Implementation

## Overview

Successfully implemented Optuna-based hyperparameter tuning with k-fold cross-validation for the XGBoost BDT training pipeline. This is an **optional feature** that can be enabled via configuration to automatically find optimal hyperparameters.

## Implementation Summary

### Files Modified

1. **requirements.txt** - Added 3 dependencies:
   - `optuna>=3.0.0` - Hyperparameter optimization framework
   - `plotly>=5.0.0` - Interactive visualization library
   - `kaleido>=0.2.0` - Static image export for plots

2. **config.yaml** - Added hyperparameter tuning configuration:
   - Training section: Added 4 tuning control parameters
   - New tuning section: Added 7 hyperparameter search ranges

3. **model_builder.py** - Added parameter update capability:
   - `update_params()` method to override model hyperparameters

4. **main_training.py** - Core implementation (4 new methods + 2 modifications):
   - `_tune_hyperparameters()` - Main tuning orchestration
   - `_create_optuna_objective()` - Defines search space and objective
   - `_stratified_kfold_cv()` - K-fold CV with sample weights
   - `_save_tuning_results()` - Saves tuning results to files
   - Modified `_train_single_model()` - Integrated tuning
   - Modified `_train_per_bin_models()` - Integrated global tuning

### Code Statistics

- **New code**: ~270 lines
- **Modified code**: ~50 lines
- **Total impact**: 4 files, 320 lines

## How It Works

### Two-Phase Training Architecture

```
Phase 1 (Optional): Hyperparameter Tuning
├─ Enabled by config flag: enable_hyperparameter_tuning: true
├─ Combines train (70%) + val (10%) → 80% for tuning
├─ Runs Optuna with 5-fold stratified cross-validation
├─ Tests hyperparameter combinations to maximize AUC
├─ Finds best hyperparameters
└─ Test set (20%) NEVER touched

Phase 2: Final Training
├─ If tuning enabled: Train with best hyperparameters on 80% data
├─ If tuning disabled: Train with config defaults on 70% data
└─ Evaluate on test set (20%)
```

### Hyperparameters Tuned

7 extended set of hyperparameters optimized via Bayesian optimization:

1. **learning_rate** (log-uniform): [0.001, 0.3]
2. **max_depth** (integer): [3, 10]
3. **n_estimators** (integer, step 50): [100, 1000]
4. **subsample** (uniform): [0.5, 1.0]
5. **reg_alpha** (log-uniform): [1e-8, 10.0]
6. **reg_lambda** (log-uniform): [1e-8, 10.0]
7. **colsample_bytree** (uniform): [0.5, 1.0]

### K-Fold Cross-Validation

- **5-fold** stratified cross-validation
- Maintains **signal/background balance** in each fold
- Preserves **sample weights** during training and evaluation
- Returns mean AUC across folds as optimization metric

### Global Tuning Strategy

For both training modes (single model and per-bin):
- Tuning runs **once globally** on all available data
- Same **best hyperparameters applied to all bins**
- Computational efficiency: O(n_trials) instead of O(n_trials × n_bins)

## Usage

### Enable Hyperparameter Tuning

Edit `config.yaml`:

```yaml
training:
  enable_hyperparameter_tuning: true  # Enable tuning
  n_trials: 100                        # Number of Optuna trials
  n_cv_folds: 5                        # K-fold splits
  tuning_timeout_seconds: null         # Optional timeout
  tuning_data_fraction: 1.0            # Use 100% of data (or 0.1 for 10%)
```

Run training:

```python
from main_training import BinnedTrainingPipeline

pipeline = BinnedTrainingPipeline("config.yaml")
pipeline.run_training()
pipeline.evaluate_models()
pipeline.generate_plots()
pipeline.save_results()
```

Or from command line:

```bash
cd /sphenix/user/shuhangli/ppg12/FunWithxgboost
python main_training.py
```

### Expected Output

```
=== Starting Binned Training Pipeline ===
Using vertex binning with 3 bins: ['-100_-30', '-30_30', '30_100']
--- Loading and reweighting all data globally ---

=== Phase 1: Hyperparameter Tuning with Optuna ===
[I 2025-12-03 14:00:00] Trial 0 finished with value: 0.8521
[I 2025-12-03 14:01:45] Trial 1 finished with value: 0.8567
[I 2025-12-03 14:03:22] Trial 2 finished with value: 0.8489
...
[I 2025-12-03 15:30:15] Trial 99 finished with value: 0.8891

Best hyperparameters: {
  'learning_rate': 0.042,
  'max_depth': 6,
  'n_estimators': 650,
  'subsample': 0.85,
  'reg_alpha': 0.001,
  'reg_lambda': 1.2,
  'colsample_bytree': 0.9
}
Best cross-validation AUC: 0.8891

Saved best hyperparameters: binned_models/best_hyperparameters.json
Saved all trials: binned_models/optuna_trials.csv
Saved optimization history: binned_models/optuna_optimization_history.html

=== Phase 2: Training with Best Hyperparameters ===
Data split: 40000 train, 5000 validation, 10000 test samples
Training complete!

=== Evaluating Models on Test Set ===
Bin -100_-30: Test AUC = 0.8923
Bin -30_30: Test AUC = 0.9012
Bin 30_100: Test AUC = 0.8856
```

### Output Files

When tuning is enabled, additional files are created:

1. **best_hyperparameters.json** - Best hyperparameters and CV score
   ```json
   {
     "best_params": {
       "learning_rate": 0.042,
       "max_depth": 6,
       ...
     },
     "best_value": 0.8891,
     "n_trials": 100
   }
   ```

2. **optuna_trials.csv** - All trials with parameters and scores
   - Columns: trial number, parameters, AUC, duration, state

3. **optuna_optimization_history.html** - Interactive plot
   - Shows convergence of hyperparameter search
   - View in web browser

### Disable Tuning (Default)

```yaml
training:
  enable_hyperparameter_tuning: false  # Default
```

Result: Original pipeline behavior, uses static hyperparameters from config.

## Configuration Reference

### Training Section

```yaml
training:
  # Existing parameters
  global_seed: 42
  stratify: true
  train_single_model: true
  train_size: 0.7
  val_size: 0.1
  test_size: 0.2

  # NEW: Hyperparameter tuning control
  enable_hyperparameter_tuning: false  # Set true to enable
  n_trials: 50                          # Number of Optuna trials
  n_cv_folds: 5                         # K-fold splits
  tuning_timeout_seconds: null          # Optional timeout (seconds)
  tuning_data_fraction: 1.0             # Fraction of data for tuning (0.1 = 10%, 1.0 = 100%)
```

### Tuning Section (NEW)

```yaml
tuning:
  # Fallback threshold
  baseline_auc_threshold: 0.50  # If CV AUC < this, use defaults

  # Learning rate: log-uniform distribution
  learning_rate_min: 0.001
  learning_rate_max: 0.3

  # Tree depth: integer
  max_depth_min: 3
  max_depth_max: 10

  # Number of trees: integer with step size
  n_estimators_min: 100
  n_estimators_max: 1000
  n_estimators_step: 50

  # Data subsampling: uniform [0.5, 1.0]
  subsample_min: 0.5
  subsample_max: 1.0

  # L1 regularization: log-uniform
  reg_alpha_min: 1.0e-8
  reg_alpha_max: 10.0

  # L2 regularization: log-uniform
  reg_lambda_min: 1.0e-8
  reg_lambda_max: 10.0

  # Feature subsampling: uniform [0.5, 1.0]
  colsample_bytree_min: 0.5
  colsample_bytree_max: 1.0
```

## Performance Considerations

### Computational Cost

With default settings (50 trials, 5 folds):
- **Total model training**: 50 trials × 5 folds = **250 model fits**
- **Runtime increase**: ~50x-100x compared to single training
- **Typical duration**: 30 minutes to 2 hours (depends on dataset size)

### Memory Usage

- K-fold creates k model copies in memory
- Each fold: ~same memory as single model
- Total increase: ~10-20% over baseline

### Optimization Strategies

1. **Reduce trials for quick iteration**:
   ```yaml
   n_trials: 20  # Faster, less thorough
   ```

2. **Reduce folds for large datasets**:
   ```yaml
   n_cv_folds: 3  # Faster, less reliable
   ```

3. **Set timeout for time-constrained runs**:
   ```yaml
   tuning_timeout_seconds: 3600  # 1 hour max
   ```

4. **Use subset of data for tuning** (NEW - fastest speedup!):
   ```yaml
   tuning_data_fraction: 0.1  # Use only 10% of data for tuning
   ```
   - **Speed**: ~10x faster with 0.1 (10%), ~5x faster with 0.2 (20%)
   - **Quality**: Good hyperparameters often found with small subset
   - **Strategy**: Tune on subset, then final model trains on full data
   - **Maintains balance**: Uses stratified sampling to preserve signal/background ratio
   - **Example**: For 1M samples → tune on 100k → 10x speedup

5. **Use smaller hyperparameter ranges**:
   ```yaml
   n_estimators_min: 200
   n_estimators_max: 600  # Narrower range
   ```

## Key Implementation Details

### Sample Weights Preservation

Critical: Sample weights are properly handled throughout:

```python
# During k-fold CV training
pipeline.fit(X_train_fold, y_train_fold, clf__sample_weight=w_train_fold)

# During k-fold CV evaluation
fold_auc = roc_auc_score(y_val_fold, y_pred, sample_weight=w_val_fold)

# Final training with tuned params
pipe.fit(X_train_final, y_train_final, clf__sample_weight=w_train_final)
```

### Stratification

Maintains signal/background balance:

```python
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=global_seed)
```

### Test Set Isolation

Test set is NEVER used during tuning:
- Tuning: train (70%) + val (10%) = 80%
- Final evaluation: test (20%)
- Prevents data leakage and overfitting

### Reproducibility

All random operations are seeded:
- Optuna TPE sampler: `seed=global_seed`
- StratifiedKFold: `random_state=global_seed`
- XGBClassifier: `random_state=global_seed`
- Background subsampling: Uses configured seed

## Troubleshooting

### Issue: Tuning takes too long

**Solutions**:
1. Reduce `n_trials` (e.g., 20-30)
2. Reduce `n_cv_folds` (e.g., 3)
3. Set `tuning_timeout_seconds`
4. Narrow hyperparameter ranges

### Issue: CV AUC below baseline

**Behavior**: Automatically falls back to config defaults

```
WARNING: Best CV AUC (0.4821) below baseline (0.50)
Falling back to config defaults
```

**Solutions**:
1. Check hyperparameter ranges are reasonable
2. Increase `n_trials` for better exploration
3. Lower `baseline_auc_threshold` if appropriate

### Issue: Memory errors during tuning

**Solutions**:
1. Reduce `n_cv_folds` from 5 to 3
2. Enable background subsampling in config
3. Reduce `n_estimators_max`

### Issue: Import errors

```
ImportError: No module named 'optuna'
```

**Solution**: Install dependencies

```bash
pip install optuna>=3.0.0 plotly>=5.0.0 kaleido>=0.2.0
```

Or install all requirements:

```bash
cd /sphenix/user/shuhangli/ppg12/FunWithxgboost
pip install -r requirements.txt
```

## Testing Recommendations

### 1. Baseline Test (No Regression)

```yaml
training:
  enable_hyperparameter_tuning: false
```

Run pipeline and verify original behavior unchanged.

### 2. Quick Tuning Test

```yaml
training:
  enable_hyperparameter_tuning: true
  n_trials: 5
  n_cv_folds: 3
```

Verify tuning runs without errors (~5 minutes).

### 3. Full Tuning Test

```yaml
training:
  enable_hyperparameter_tuning: true
  n_trials: 50
  n_cv_folds: 5
```

Verify performance improvement over defaults (~1 hour).

### 4. Per-Bin Mode Test

```yaml
training:
  train_single_model: false  # Per-bin mode
  enable_hyperparameter_tuning: true
  n_trials: 20
```

Verify global tuning applies to all bins.

## Benefits

✅ **Automated hyperparameter optimization** - No manual grid search
✅ **Bayesian optimization** - Intelligent exploration with TPE sampler
✅ **Unbiased evaluation** - Test set isolation prevents overfitting
✅ **Reproducible** - Seeded random operations
✅ **Sample weight aware** - Preserves physics reweighting
✅ **Stratified CV** - Maintains signal/background balance
✅ **Comprehensive logging** - Saves all trials and visualizations
✅ **Backward compatible** - Opt-in feature, zero impact when disabled
✅ **Configurable** - All tuning parameters in config file

## Future Enhancements

Potential improvements for future versions:

1. **Distributed tuning** - Parallel trial execution via Optuna storage
2. **Early stopping** - Pruning unpromising trials
3. **Multi-objective optimization** - Optimize AUC + training time
4. **Per-bin tuning** - Separate hyperparameters for each bin
5. **Warm start** - Resume from previous tuning runs
6. **Custom search spaces** - User-defined hyperparameter distributions

## References

- **Optuna Documentation**: https://optuna.readthedocs.io/
- **XGBoost Documentation**: https://xgboost.readthedocs.io/
- **Scikit-learn Cross-Validation**: https://scikit-learn.org/stable/modules/cross_validation.html

## Date

Implementation completed: 2025-12-03

## Contact

For questions or issues, please refer to the main pipeline documentation or contact the development team.
