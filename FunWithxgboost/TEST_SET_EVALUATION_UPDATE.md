# Test Set Evaluation Update

## Summary
Updated the BDT training pipeline to use a proper **train/validation/test split** with final evaluation on the **test set** for unbiased performance metrics.

## What Changed

### Previous Behavior
- Data was split into **train (80%)** and **validation (20%)** only
- Final evaluation was performed on the **validation set**
- No separate holdout test set for unbiased evaluation
- Risk of overfitting to validation set during model development

### New Behavior
- Data is now split into **train (70%)**, **validation (10%)**, and **test (20%)**
- Final evaluation is performed on the **test set**
- Test set remains completely unseen until final evaluation
- Provides unbiased performance estimates

## Modified Files

### 1. `main_training.py`

#### `_split_train_val_test_with_features()` (new method)
- Replaces `_split_train_val_with_features()`
- Performs two-stage splitting:
  1. Split off test set (20% by default)
  2. Split remaining data into train/val (70%/10% by default)
- Returns: `X_train, X_val, X_test, y_train, y_val, y_test, w_train, w_val, w_test`

#### `_train_single_model()`
- Updated to use new 3-way split
- Stores test data instead of validation data in `trained_pipelines`
- Prints split sizes for transparency
- Keys changed: `X_val` → `X_test`, `y_val` → `y_test`, `w_val` → `w_test`, `val_index` → `test_index`
- Added `n_val` and `n_test` to metrics

#### `_train_per_bin_models()`
- Same updates as `_train_single_model()` but for per-bin training mode

#### `evaluate_models()`
- Now evaluates on **test set** instead of validation set
- Updated docstring: "Evaluate all trained models on the test set"
- Prints test AUC for each bin
- Results stored with both `auc_test` and `auc_val` keys (backward compatibility)

#### `_save_metrics()`
- Added `auc_test` and `n_test` columns to CSV output
- Maintains `auc_val` for backward compatibility (now contains test AUC)

### 2. `plotting.py`

#### `plot_bdt_score_distributions()`
- Updated to plot BDT scores from **test set**
- Docstring updated: "Plot BDT score distributions for each bin on test set"

#### `plot_isoet_bdt_correlation_global()`
- Updated to compute correlations on **test set**
- Docstring updated: "Global isoET vs BDT score correlation (all bins combined) on test set"

#### `plot_auc_bar()`
- Updated labels: "Validation AUC" → "Test AUC"
- Title: "AUC by pT bin" → "Test AUC by pT bin"

### 3. `config.yaml`

```yaml
training:
  train_size: 0.7  # 70% for training
  val_size: 0.1    # 10% for validation (for future hyperparameter tuning)
  test_size: 0.2   # 20% for final unbiased evaluation
```

## Configuration

The split ratios are fully configurable in `config.yaml`:

```yaml
training:
  train_size: 0.7   # Training set fraction
  val_size: 0.1     # Validation set fraction (for hyperparameter tuning)
  test_size: 0.2    # Test set fraction (for final evaluation)
```

**Important:** `train_size + val_size + test_size` must equal 1.0

### Recommended Configurations

**Default (70/10/20):**
- Good balance for most use cases
- Enough training data
- Validation set for future hyperparameter optimization
- Test set large enough for reliable estimates

**More Training Data (80/0/20):**
```yaml
train_size: 0.8
val_size: 0.0
test_size: 0.2
```
- Use when not doing hyperparameter tuning
- Maximum training data

**Strict Separation (60/20/20):**
```yaml
train_size: 0.6
val_size: 0.2
test_size: 0.2
```
- Use for hyperparameter optimization
- Larger validation set for better tuning

## Benefits

1. **Unbiased Evaluation**: Test set provides unbiased performance estimates
2. **Better Science**: Follows machine learning best practices
3. **Future-Proof**: Validation set available for hyperparameter tuning
4. **Transparency**: Split sizes printed during training
5. **Reproducibility**: Same random seed ensures consistent splits
6. **Backward Compatible**: Old `auc_val` key still available in results

## Output Changes

### Metrics CSV (`per_bin_metrics.csv`)

**New columns:**
- `n_val`: Number of validation samples
- `n_test`: Number of test samples
- `auc_test`: Test set AUC (unbiased performance)

**Note:** `auc_val` now contains the same value as `auc_test` for backward compatibility.

### Console Output

During training:
```
Data split: 35000 train, 5000 validation, 10000 test samples
```

During evaluation:
```
=== Evaluating Models on Test Set ===
Bin 6_10: Test AUC = 0.9234
Bin 10_15: Test AUC = 0.9456
...
```

### Plots

- All plot titles and labels now reference "Test" instead of "Validation"
- Plots use test set data only
- BDT score distributions on test set
- Correlations computed on test set

## Validation

To verify the changes are working correctly:

1. **Check split sizes:**
   - Look for "Data split: X train, Y validation, Z test samples" in output
   - Verify ratios match config (70/10/20 by default)

2. **Check test AUC:**
   - Look for "Test AUC = ..." messages during evaluation
   - Check `per_bin_metrics.csv` for `auc_test` column

3. **Verify plots:**
   - BDT score distribution plots should show test set data
   - Plot titles should reference "Test" not "Validation"

## Migration Notes

- **No action required** for users - backward compatible
- Existing code using `auc_val` will still work
- New code should use `auc_test` for clarity
- Plot labels automatically updated

## Technical Details

### Split Strategy

1. **First split:** Separate test set (stratified by class labels)
   - Random seed: `global_seed`
   - Test size: `test_size`

2. **Second split:** Separate train/validation (stratified by class labels)
   - Random seed: `global_seed + 1` (different from first split)
   - Validation size: `val_size / (train_size + val_size)`

### Stratification

- Both splits maintain class balance (signal/background ratio)
- Controlled by `stratify: true` in config
- Ensures representative samples in all sets

### Random Seed

- First split: Uses `global_seed` (default: 42)
- Second split: Uses `global_seed + 1` (default: 43)
- Ensures reproducibility while maintaining independence

## Date
2025-12-03
