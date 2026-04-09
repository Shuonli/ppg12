#!/usr/bin/env python3
"""
NPB Score Training Script

Trains a BDT model to distinguish NPB (detector artifacts) from real clusters.
- Positive (label=1): Real clusters from Jet MC
- Negative (label=0): NPB-tagged clusters from data

Usage:
    python3 train_npb_score.py --config config_npb_training.yaml
"""

import os
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import itertools
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    roc_curve, auc, classification_report,
    confusion_matrix, roc_auc_score
)
from scipy.interpolate import UnivariateSpline
import xgboost as xgb
import joblib

# Optional ROOT support for TMVA export
HAS_ROOT = False
try:
    import ROOT
    HAS_ROOT = True
except ImportError:
    print("WARNING: ROOT not available. Cannot save .root format.")


class NPBTrainer:
    """NPB rejection score trainer"""

    def __init__(self, config_path):
        """Load configuration"""
        with open(config_path) as f:
            self.config = yaml.safe_load(f)

        self.output_dir = Path(self.config['training']['output_dir'])
        self.output_dir.mkdir(exist_ok=True, parents=True)

        self.feature_cols = self.config['features']['feature_list']
        self.model = None
        self.best_params_ = None
        self.best_val_auc_ = None

        # Per-bin analysis config and storage
        self.per_bin_config = self.config.get('per_bin_analysis', {})
        self.per_bin_models = {}
        self.per_bin_auc_single = {}
        self.per_bin_auc_binned = {}

    def load_data(self):
        """Load signal (Jet MC) and background (NPB data) files"""
        print("\n[1/7] Loading data files...")

        # Load signal (real clusters from Jet MC)
        signal_dfs = []
        for fname in self.config['data']['signal_files']:
            if not os.path.exists(fname):
                print(f"  WARNING: {fname} not found, skipping")
                continue
            df = pd.read_csv(fname, sep=r'\s+')
            signal_dfs.append(df)
            print(f"  Loaded {fname}: {len(df)} entries")

        df_signal = pd.concat(signal_dfs, ignore_index=True) if signal_dfs else pd.DataFrame()
        df_signal['label'] = 1  # Real clusters

        # Load background (NPB from data)
        bg_dfs = []
        for fname in self.config['data']['background_files']:
            if not os.path.exists(fname):
                print(f"  WARNING: {fname} not found, skipping")
                continue
            df = pd.read_csv(fname, sep=r'\s+')
            bg_dfs.append(df)
            print(f"  Loaded {fname}: {len(df)} entries")

        df_background = pd.concat(bg_dfs, ignore_index=True) if bg_dfs else pd.DataFrame()
        df_background['label'] = 0  # NPB artifacts

        # Apply selections
        sig_sel = self.config['data']['signal_selection']
        bg_sel = self.config['data']['background_selection']
        df_signal = df_signal.query(sig_sel) if sig_sel else df_signal
        df_background = df_background.query(bg_sel) if bg_sel else df_background

        print(f"  Signal after selection: {len(df_signal)}")
        print(f"  Background after selection: {len(df_background)}")

        # Combine
        df_all = pd.concat([df_signal, df_background], ignore_index=True)

        # Apply kinematic cuts
        et_min = self.config['data']['et_min']
        et_max = self.config['data']['et_max']
        eta_max = self.config['data']['eta_max']

        df_all = df_all.query(f"{et_min} <= cluster_Et <= {et_max}")
        df_all = df_all.query(f"abs(cluster_Eta) <= {eta_max}")

        print(f"  Total after kinematic cuts: {len(df_all)}")
        print(f"    Signal: {(df_all['label'] == 1).sum()}")
        print(f"    Background: {(df_all['label'] == 0).sum()}")

        # ------------------------------------------------------------------
        # Robustness: drop non-finite feature rows (inf / -inf / NaN)
        #
        # XGBoost will hard-fail on inf values (common when ratio features have
        # zero denominators in the upstream txt export).
        # We ensure training is stable even if inputs contain such rows.
        # ------------------------------------------------------------------
        if self.feature_cols:
            feat_df = df_all[self.feature_cols].copy()
            # Replace inf with NaN then drop
            feat_df = feat_df.replace([np.inf, -np.inf], np.nan)
            mask_ok = feat_df.notna().all(axis=1)
            n_bad = int((~mask_ok).sum())
            if n_bad > 0:
                print(f"  Dropping non-finite feature rows: {n_bad}")
                df_all = df_all.loc[mask_ok].copy()
                print(f"  Total after non-finite cleanup: {len(df_all)}")
                print(f"    Signal: {(df_all['label'] == 1).sum()}")
                print(f"    Background: {(df_all['label'] == 0).sum()}")

        # Downsample signal peak region before splitting/reweighting
        if self.config['training'].get('downsample_signal_peak', False):
            df_all = self._downsample_signal_peak(df_all)

        return df_all

    def _downsample_signal_peak(self, df):
        """Downsample signal events in the low-ET peak region.

        The inclusive MC signal spectrum is steeply falling, producing a huge
        peak at low ET (e.g. 6-15 GeV).  This makes inverse-PDF reweighting
        weights very large in the tail, even with capping.

        Strategy: compute a histogram of signal ET.  Use the bin density at
        the upper edge of the peak region as a reference.  In each bin inside
        the peak region, randomly drop events so the bin count matches that
        reference level.  Events outside the peak region and all background
        events are kept untouched.
        """
        training_cfg = self.config['training']
        seed = training_cfg['random_seed']
        peak_lo = float(training_cfg.get('downsample_et_min', 6.0))
        peak_hi = float(training_cfg.get('downsample_et_max', 15.0))
        target_frac = training_cfg.get('downsample_target_fraction', None)
        if target_frac is not None:
            target_frac = float(target_frac)

        sig_mask = (df['label'] == 1)
        n_sig_before = int(sig_mask.sum())

        if n_sig_before == 0:
            return df

        print(f"\n  Downsampling signal peak in ET [{peak_lo}, {peak_hi}] GeV...")

        rng = np.random.RandomState(seed)

        if target_frac is not None:
            # Simple mode: keep a fixed fraction of signal in the peak region
            in_peak = sig_mask & (df['cluster_Et'] >= peak_lo) & (df['cluster_Et'] < peak_hi)
            n_in_peak = int(in_peak.sum())
            n_keep = int(n_in_peak * target_frac)

            peak_indices = df.index[in_peak].values
            keep_indices = rng.choice(peak_indices, size=n_keep, replace=False)
            drop_indices = np.setdiff1d(peak_indices, keep_indices)

            df = df.drop(index=drop_indices).reset_index(drop=True)

            print(f"    Fixed fraction mode: kept {n_keep}/{n_in_peak} "
                  f"signal events in peak (fraction={target_frac})")
        else:
            # Auto mode: use reference bin density from upper edge of peak
            et_vals = df.loc[sig_mask, 'cluster_Et'].values
            et_min_data = float(self.config['data']['et_min'])
            et_max_data = float(self.config['data']['et_max'])

            # Build histogram with 1 GeV bins
            bin_width = 1.0
            hist_bins = np.arange(et_min_data, et_max_data + bin_width, bin_width)
            hist_counts, bin_edges = np.histogram(et_vals, bins=hist_bins)

            # Find the reference bin: the first full bin at or above peak_hi
            ref_bin_idx = None
            for j in range(len(bin_edges) - 1):
                if bin_edges[j] >= peak_hi:
                    ref_bin_idx = j
                    break
            if ref_bin_idx is None:
                print("    WARNING: could not find reference bin, skipping downsampling")
                return df

            ref_count = int(hist_counts[ref_bin_idx])
            print(f"    Reference bin [{bin_edges[ref_bin_idx]:.0f}, "
                  f"{bin_edges[ref_bin_idx+1]:.0f}] GeV: {ref_count} signal events")

            # For each bin in the peak region, downsample to ref_count
            all_drop = []
            for j in range(len(bin_edges) - 1):
                blo, bhi = bin_edges[j], bin_edges[j + 1]
                if blo >= peak_hi or bhi <= peak_lo:
                    continue  # outside peak region

                in_bin = sig_mask & (df['cluster_Et'] >= blo) & (df['cluster_Et'] < bhi)
                n_in_bin = int(in_bin.sum())

                if n_in_bin <= ref_count:
                    continue  # already at or below target

                bin_indices = df.index[in_bin].values
                keep = rng.choice(bin_indices, size=ref_count, replace=False)
                drop = np.setdiff1d(bin_indices, keep)
                all_drop.append(drop)

                print(f"    Bin [{blo:.0f}, {bhi:.0f}) GeV: {n_in_bin} -> {ref_count} signal events")

            if all_drop:
                drop_indices = np.concatenate(all_drop)
                df = df.drop(index=drop_indices).reset_index(drop=True)

        n_sig_after = int((df['label'] == 1).sum())
        print(f"    Signal: {n_sig_before} -> {n_sig_after} "
              f"(dropped {n_sig_before - n_sig_after})")

        return df

    def split_data(self, df):
        """Split into train/val/test sets, returning ET values for reweighting"""
        print("\n[2/7] Splitting data...")

        train_size = self.config['training']['train_size']
        val_size = self.config['training']['val_size']
        test_size = self.config['training']['test_size']
        seed = self.config['training']['random_seed']

        # Split DataFrame first (stratified by label)
        df_temp, df_test = train_test_split(
            df, test_size=test_size, random_state=seed, stratify=df['label']
        )

        val_frac = val_size / (train_size + val_size)
        df_train, df_val = train_test_split(
            df_temp, test_size=val_frac, random_state=seed, stratify=df_temp['label']
        )

        # Extract feature arrays, labels, and ET values
        X_train = df_train[self.feature_cols].values
        X_val = df_val[self.feature_cols].values
        X_test = df_test[self.feature_cols].values

        y_train = df_train['label'].values
        y_val = df_val['label'].values
        y_test = df_test['label'].values

        et_train = df_train['cluster_Et'].values
        et_val = df_val['cluster_Et'].values
        et_test = df_test['cluster_Et'].values

        print(f"  Train: {len(X_train)} ({np.sum(y_train==1)} signal, {np.sum(y_train==0)} bkg)")
        print(f"  Val:   {len(X_val)} ({np.sum(y_val==1)} signal, {np.sum(y_val==0)} bkg)")
        print(f"  Test:  {len(X_test)} ({np.sum(y_test==1)} signal, {np.sum(y_test==0)} bkg)")

        return X_train, X_val, X_test, y_train, y_val, y_test, et_train, et_val, et_test

    def _compute_et_weights(self, et_values, label_mask):
        """Compute ET reweighting factors using inverse-PDF method (per-class).

        Algorithm (following reweighting.py pattern):
        1. Compute histogram density estimate
        2. Fit spline for smooth interpolation
        3. Evaluate PDF at each sample's ET value
        4. Compute weights = 1 / PDF
        5. Cap maximum weight and normalize

        Args:
            et_values: numpy array of ET values for all samples
            label_mask: boolean mask selecting samples for this class

        Returns:
            numpy array of weights (same length as et_values)
        """
        training_cfg = self.config['training']
        n_bins = int(training_cfg.get('et_reweight_bins', 20))
        weight_cap = float(training_cfg.get('et_reweight_cap', 10.0))
        min_pdf = float(training_cfg.get('et_reweight_min_pdf', 1e-3))

        # Initialize all weights to 1
        weights = np.ones(len(et_values))

        # Get ET values for this class
        et_class = et_values[label_mask]
        if len(et_class) == 0:
            return weights

        # Compute histogram density
        et_min, et_max = et_class.min(), et_class.max()
        bins = np.linspace(et_min, et_max, n_bins + 1)
        hist, bin_edges = np.histogram(et_class, bins=bins, density=True)
        # Scale histogram so integral ~1 over the bin width
        hist = hist * n_bins
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        # Fit spline for smooth interpolation
        spline = UnivariateSpline(bin_centers, hist, s=0.0, ext=3)

        # Evaluate PDF at each sample's ET
        pdf_vals = spline(et_class)

        # Clip minimum PDF to prevent division by ~0
        pdf_vals = np.clip(pdf_vals, a_min=min_pdf, a_max=None)

        # Compute inverse-PDF weights
        class_weights = 1.0 / pdf_vals

        # Cap maximum weight
        class_weights = np.clip(class_weights, a_min=None, a_max=weight_cap)

        # Normalize so mean = 1
        class_weights /= np.mean(class_weights)

        # Assign to output array
        weights[label_mask] = class_weights

        return weights

    def compute_weights(self, y, et_values=None):
        """Compute class-balanced sample weights, optionally with ET reweighting"""
        print("\n[3/7] Computing sample weights...")

        if not self.config['training']['balance_classes']:
            weights = np.ones(len(y))
            print("  No class balancing (equal weights)")
        else:
            # Class weights: make signal and background contribute equally
            n_signal = np.sum(y == 1)
            n_background = np.sum(y == 0)
            total = n_signal + n_background

            w_signal = total / (2.0 * n_signal) if n_signal > 0 else 1.0
            w_background = total / (2.0 * n_background) if n_background > 0 else 1.0

            weights = np.where(y == 1, w_signal, w_background)

            print(f"  Class balancing:")
            print(f"    Signal weight: {w_signal:.4f}")
            print(f"    Background weight: {w_background:.4f}")

        # Apply ET reweighting if enabled
        if self.config['training'].get('reweight_et', False) and et_values is not None:
            print("  ET reweighting (inverse-PDF):")

            # Compute ET weights for signal (y==1)
            sig_mask = (y == 1)
            et_weights_sig = self._compute_et_weights(et_values, sig_mask)

            # Compute ET weights for background (y==0)
            bkg_mask = (y == 0)
            et_weights_bkg = self._compute_et_weights(et_values, bkg_mask)

            # Combine: each sample gets its class-specific ET weight
            et_weights = np.where(y == 1, et_weights_sig, et_weights_bkg)

            # Compound weighting: multiply class weights by ET weights
            weights = weights * et_weights

            print(f"    Signal ET weights: min={et_weights[sig_mask].min():.3f}, "
                  f"max={et_weights[sig_mask].max():.3f}, mean={et_weights[sig_mask].mean():.3f}")
            print(f"    Background ET weights: min={et_weights[bkg_mask].min():.3f}, "
                  f"max={et_weights[bkg_mask].max():.3f}, mean={et_weights[bkg_mask].mean():.3f}")
            print(f"    Final weights: min={weights.min():.3f}, max={weights.max():.3f}, "
                  f"mean={weights.mean():.3f}")

        return weights

    def _get_tuning_config(self):
        """Return hyperparameter tuning config (backward compatible)."""
        training_cfg = self.config.get('training', {})

        # New preferred schema:
        # training:
        #   hyperparameter_tuning:
        #     enabled: true
        #     max_evals: 10
        #     param_grid: { max_depth: [4,6], learning_rate: [0.05,0.1], ... }
        tuning_cfg = training_cfg.get('hyperparameter_tuning', {})

        # Backward compatible alias (if user adds a single flat key)
        enabled = bool(
            tuning_cfg.get('enabled', False) or training_cfg.get('enable_hyperparameter_tuning', False)
        )

        max_evals = int(tuning_cfg.get('max_evals', 12))
        metric = str(tuning_cfg.get('metric', 'auc')).lower()

        # If no grid provided, use a minimal one around common defaults
        default_grid = {
            'max_depth': [4, 6, 8],
            'learning_rate': [0.05, 0.1],
            'subsample': [0.7, 0.9],
            'colsample_bytree': [0.7, 0.9],
            'reg_alpha': [0.0, 1.0],
            'reg_lambda': [1.0],
            'n_estimators': [300, 600, 900],
        }
        param_grid = tuning_cfg.get('param_grid', default_grid)

        # Safety: only keep keys that exist in base params (or are known XGBClassifier args)
        base_params = self.config.get('model', {}).get('params', {}).copy()
        allowed = set(base_params.keys()) | {
            'n_estimators', 'max_depth', 'learning_rate', 'subsample', 'colsample_bytree',
            'min_child_weight', 'gamma', 'reg_alpha', 'reg_lambda'
        }
        filtered_grid = {k: v for k, v in param_grid.items() if k in allowed}

        return {
            'enabled': enabled,
            'max_evals': max_evals,
            'metric': metric,
            'param_grid': filtered_grid,
        }

    def _iter_param_combinations(self, param_grid):
        """Yield dicts for each parameter combination (deterministic order)."""
        if not param_grid:
            yield {}
            return

        keys = list(param_grid.keys())
        values = []
        for k in keys:
            v = param_grid[k]
            if isinstance(v, (list, tuple, np.ndarray)):
                values.append(list(v))
            else:
                values.append([v])

        for combo in itertools.product(*values):
            yield dict(zip(keys, combo))

    def tune_hyperparameters(self, X_train, y_train, X_val, y_val, sample_weights):
        """Minimal hyperparameter tuning: try a small grid, pick best by val metric."""
        tuning_cfg = self._get_tuning_config()
        if not tuning_cfg['enabled']:
            return None, None

        metric = tuning_cfg['metric']
        if metric != 'auc':
            print(f"  WARNING: tuning metric '{metric}' not supported; using AUC")
            metric = 'auc'

        base_params = self.config['model']['params'].copy()
        seed = self.config['training']['random_seed']

        combos = list(self._iter_param_combinations(tuning_cfg['param_grid']))
        if len(combos) == 0:
            combos = [{}]

        # Cap the number of evaluations (shuffle deterministically for fairness)
        rs = np.random.RandomState(seed)
        rs.shuffle(combos)
        combos = combos[: max(1, tuning_cfg['max_evals'])]

        print("\n[4/7] Hyperparameter tuning (minimal grid search)...")
        print(f"  Trials: {len(combos)} (cap={tuning_cfg['max_evals']})")

        best_auc = -np.inf
        best_params = None

        for i, delta in enumerate(combos, start=1):
            trial_params = base_params.copy()
            trial_params.update(delta)

            model = xgb.XGBClassifier(
                **trial_params,
                random_state=seed,
                use_label_encoder=False
            )

            model.fit(
                X_train, y_train,
                sample_weight=sample_weights,
                eval_set=[(X_train, y_train), (X_val, y_val)],
                verbose=False
            )

            val_pred = model.predict_proba(X_val)[:, 1]
            val_auc = roc_auc_score(y_val, val_pred)

            print(f"  Trial {i:02d}/{len(combos)}: val_auc={val_auc:.5f} params={delta}")

            if val_auc > best_auc:
                best_auc = val_auc
                best_params = trial_params

        self.best_params_ = best_params
        self.best_val_auc_ = float(best_auc) if np.isfinite(best_auc) else None

        print(f"\n  Best val AUC: {best_auc:.5f}")
        print(f"  Best params delta vs base: { {k: best_params[k] for k in tuning_cfg['param_grid'].keys() if k in best_params} }")

        return best_params, best_auc

    def train_model(self, X_train, y_train, X_val, y_val, sample_weights):
        """Train XGBoost model"""
        print("\n[4/7] Training XGBoost model...")

        params = self.config['model']['params'].copy()
        seed = self.config['training']['random_seed']

        # Optional minimal hyperparameter tuning on the validation split
        best_params, _best_auc = self.tune_hyperparameters(X_train, y_train, X_val, y_val, sample_weights)
        if best_params is not None:
            params = best_params

        self.model = xgb.XGBClassifier(
            **params,
            random_state=seed,
            use_label_encoder=False
        )

        # If we tuned, retrain on (train+val) for better final performance before testing.
        if best_params is not None:
            X_train2 = np.concatenate([X_train, X_val], axis=0)
            y_train2 = np.concatenate([y_train, y_val], axis=0)
            # Combine ET values if available (for ET reweighting)
            et_train2 = None
            if hasattr(self, 'et_train') and hasattr(self, 'et_val'):
                et_train2 = np.concatenate([self.et_train, self.et_val], axis=0)
            w_train2 = self.compute_weights(y_train2, et_values=et_train2)
            eval_set = [(X_train2, y_train2)]
            self.model.fit(
                X_train2, y_train2,
                sample_weight=w_train2,
                eval_set=eval_set,
                verbose=50
            )
        else:
            self.model.fit(
                X_train, y_train,
                sample_weight=sample_weights,
                eval_set=[(X_train, y_train), (X_val, y_val)],
                verbose=50
            )

        # Training metrics
        train_pred = self.model.predict_proba(X_train)[:, 1]
        val_pred = self.model.predict_proba(X_val)[:, 1]

        train_auc = roc_auc_score(y_train, train_pred)
        val_auc = roc_auc_score(y_val, val_pred)

        print(f"\n  Train AUC: {train_auc:.4f}")
        print(f"  Val AUC:   {val_auc:.4f}")
        if self.best_val_auc_ is not None:
            print(f"  Best tuned Val AUC: {self.best_val_auc_:.4f}")

    def evaluate_model(self, X_test, y_test):
        """Evaluate on test set"""
        print("\n[5/7] Evaluating on test set...")

        y_pred_proba = self.model.predict_proba(X_test)[:, 1]
        y_pred = (y_pred_proba > 0.5).astype(int)

        # ROC curve
        fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
        test_auc = auc(fpr, tpr)

        print(f"  Test AUC: {test_auc:.4f}")
        print("\nClassification Report:")
        print(classification_report(y_test, y_pred, target_names=['NPB', 'Real']))

        return test_auc, fpr, tpr, y_pred_proba, y_pred

    def generate_plots(self, X_test, y_test, test_auc, fpr, tpr, y_pred_proba, y_pred):
        """Generate evaluation plots"""
        if not self.config['evaluation']['generate_plots']:
            return

        print("\n[6/7] Generating plots...")

        # 1. ROC Curve
        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, label=f'NPB Score (AUC = {test_auc:.4f})', lw=2)
        plt.plot([0, 1], [0, 1], 'k--', label='Random', lw=1)
        plt.xlabel('False Positive Rate (NPB → Real)', fontsize=12)
        plt.ylabel('True Positive Rate (Real → Real)', fontsize=12)
        plt.title('NPB Score ROC Curve', fontsize=14)
        plt.legend(fontsize=11)
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'roc_curve.pdf', bbox_inches='tight')
        plt.close()
        print(f"  Saved: roc_curve.pdf")

        # 2. Score Distributions
        plt.figure(figsize=(10, 6))
        scores_real = y_pred_proba[y_test == 1]
        scores_npb = y_pred_proba[y_test == 0]

        plt.hist(scores_npb, bins=50, alpha=0.6, label=f'NPB (n={len(scores_npb)})',
                 density=True, color='red')
        plt.hist(scores_real, bins=50, alpha=0.6, label=f'Real (n={len(scores_real)})',
                 density=True, color='blue')
        plt.xlabel('NPB Score', fontsize=12)
        plt.ylabel('Normalized Entries', fontsize=12)
        plt.title('NPB Score Distributions', fontsize=14)
        plt.legend(fontsize=11)
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'score_distributions.pdf', bbox_inches='tight')
        plt.close()
        print(f"  Saved: score_distributions.pdf")

        # 3. Feature Importance
        importance = self.model.feature_importances_
        indices = np.argsort(importance)[::-1][:20]  # Top 20

        plt.figure(figsize=(10, 10))
        plt.barh(range(len(indices)), importance[indices], color='steelblue')
        plt.yticks(range(len(indices)), [self.feature_cols[i] for i in indices])
        plt.xlabel('Feature Importance (Gain)', fontsize=12)
        plt.title('Top 20 Features for NPB Identification', fontsize=14)
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(self.output_dir / 'feature_importance.pdf', bbox_inches='tight')
        plt.close()
        print(f"  Saved: feature_importance.pdf")

        # 4. Confusion Matrix
        cm = confusion_matrix(y_test, y_pred)
        plt.figure(figsize=(7, 6))
        plt.imshow(cm, interpolation='nearest', cmap='Blues')
        plt.title('Confusion Matrix (threshold=0.5)', fontsize=14)
        plt.colorbar()
        tick_marks = np.arange(2)
        plt.xticks(tick_marks, ['NPB', 'Real'], fontsize=11)
        plt.yticks(tick_marks, ['NPB', 'Real'], fontsize=11)

        # Add text annotations
        thresh = cm.max() / 2.
        for i in range(2):
            for j in range(2):
                plt.text(j, i, format(cm[i, j], 'd'),
                        ha="center", va="center",
                        color="white" if cm[i, j] > thresh else "black",
                        fontsize=14)

        plt.ylabel('True Label', fontsize=12)
        plt.xlabel('Predicted Label', fontsize=12)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'confusion_matrix.pdf', bbox_inches='tight')
        plt.close()
        print(f"  Saved: confusion_matrix.pdf")

    def plot_et_distributions(self, df):
        """Plot overlaid histograms of signal/background ET distributions"""
        if not self.per_bin_config.get('plot_et_distributions', True):
            return

        print("\n  Plotting ET distributions...")

        df_signal = df[df['label'] == 1]
        df_background = df[df['label'] == 0]

        n_bins = self.per_bin_config.get('et_bins_for_plot', 40)
        et_min = self.config['data']['et_min']
        et_max = self.config['data']['et_max']

        plt.figure(figsize=(10, 6))
        plt.hist(df_signal['cluster_Et'], bins=n_bins, range=(et_min, et_max),
                 alpha=0.6, label=f'Signal (n={len(df_signal)})',
                 density=True, color='blue')
        plt.hist(df_background['cluster_Et'], bins=n_bins, range=(et_min, et_max),
                 alpha=0.6, label=f'Background (n={len(df_background)})',
                 density=True, color='red')

        plt.xlabel('cluster_Et [GeV]', fontsize=12)
        plt.ylabel('Normalized Entries', fontsize=12)
        plt.title('ET Distributions: Signal vs Background', fontsize=14)
        plt.legend(fontsize=11)
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'et_distributions.pdf', bbox_inches='tight')
        plt.close()
        print(f"  Saved: et_distributions.pdf")

    def plot_et_reweighting_qa(self, et_values, y, weights):
        """Plot ET distributions before and after reweighting for QA.

        Creates side-by-side plots showing:
        - Left: Original ET distributions for signal and background
        - Right: Weighted ET distributions (should be flatter)
        """
        print("\n  Plotting ET reweighting QA...")

        et_min = self.config['data']['et_min']
        et_max = self.config['data']['et_max']
        n_bins = 30

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Signal and background masks
        sig_mask = (y == 1)
        bkg_mask = (y == 0)

        # Left: Unweighted distributions
        ax = axes[0]
        ax.hist(et_values[sig_mask], bins=n_bins, range=(et_min, et_max),
                alpha=0.6, label=f'Signal (n={sig_mask.sum()})', density=True, color='blue')
        ax.hist(et_values[bkg_mask], bins=n_bins, range=(et_min, et_max),
                alpha=0.6, label=f'Background (n={bkg_mask.sum()})', density=True, color='red')
        ax.set_xlabel('cluster_Et [GeV]', fontsize=12)
        ax.set_ylabel('Normalized Entries', fontsize=12)
        ax.set_title('Before ET Reweighting', fontsize=14)
        ax.legend(fontsize=10)
        ax.grid(alpha=0.3)

        # Right: Weighted distributions
        ax = axes[1]
        ax.hist(et_values[sig_mask], bins=n_bins, range=(et_min, et_max),
                weights=weights[sig_mask], alpha=0.6, label='Signal (weighted)',
                density=True, color='blue')
        ax.hist(et_values[bkg_mask], bins=n_bins, range=(et_min, et_max),
                weights=weights[bkg_mask], alpha=0.6, label='Background (weighted)',
                density=True, color='red')
        ax.set_xlabel('cluster_Et [GeV]', fontsize=12)
        ax.set_ylabel('Normalized Entries', fontsize=12)
        ax.set_title('After ET Reweighting', fontsize=14)
        ax.legend(fontsize=10)
        ax.grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(self.output_dir / 'et_reweighting_qa.pdf', bbox_inches='tight')
        plt.close()
        print(f"  Saved: et_reweighting_qa.pdf")

    def calculate_per_bin_auc(self, X_test, y_test, y_pred_proba):
        """Calculate AUC per ET bin using the single trained model"""
        print("\n  Calculating per-bin AUC (single model)...")

        et_bin_edges = self.per_bin_config.get('et_bin_edges', [5, 10, 15, 20, 25, 30, 35, 40])
        et_bin_labels = self.per_bin_config.get('et_bin_labels',
                                                 ["5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40"])

        # cluster_Et is at index 0 in feature list
        et_idx = self.feature_cols.index('cluster_Et')
        et_values = X_test[:, et_idx]

        per_bin_auc = {}

        for i, label in enumerate(et_bin_labels):
            if i >= len(et_bin_edges) - 1:
                break
            et_low = et_bin_edges[i]
            et_high = et_bin_edges[i + 1]

            # Select entries in this ET bin
            mask = (et_values >= et_low) & (et_values < et_high)
            y_bin = y_test[mask]
            pred_bin = y_pred_proba[mask]

            if len(y_bin) < 10 or len(np.unique(y_bin)) < 2:
                print(f"    {label} GeV: insufficient data (n={len(y_bin)})")
                per_bin_auc[label] = np.nan
                continue

            bin_auc = roc_auc_score(y_bin, pred_bin)
            per_bin_auc[label] = bin_auc
            n_sig = np.sum(y_bin == 1)
            n_bkg = np.sum(y_bin == 0)
            print(f"    {label} GeV: AUC = {bin_auc:.4f} (n_sig={n_sig}, n_bkg={n_bkg})")

        return per_bin_auc

    def plot_per_bin_auc(self, per_bin_auc_dict, title_suffix=""):
        """Plot bar chart of AUC vs ET bin"""
        print(f"\n  Plotting per-bin AUC {title_suffix}...")

        labels = list(per_bin_auc_dict.keys())
        auc_values = [per_bin_auc_dict[l] for l in labels]

        plt.figure(figsize=(10, 6))
        x_pos = np.arange(len(labels))
        colors = ['steelblue' if not np.isnan(v) else 'gray' for v in auc_values]

        bars = plt.bar(x_pos, [v if not np.isnan(v) else 0 for v in auc_values],
                       color=colors, edgecolor='black', linewidth=1.2)

        # Add value labels on bars
        for i, (bar, val) in enumerate(zip(bars, auc_values)):
            if not np.isnan(val):
                plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                         f'{val:.3f}', ha='center', va='bottom', fontsize=10)

        plt.xticks(x_pos, labels, fontsize=11)
        plt.xlabel('ET Bin [GeV]', fontsize=12)
        plt.ylabel('AUC', fontsize=12)
        plt.title(f'Per-bin AUC {title_suffix}', fontsize=14)
        plt.ylim(0, 1.1)
        plt.grid(alpha=0.3, axis='y')
        plt.tight_layout()

        # Determine filename based on title
        if 'Single' in title_suffix:
            fname = 'per_bin_auc_single_model.pdf'
        elif 'Per-bin' in title_suffix or 'Binned' in title_suffix:
            fname = 'per_bin_auc_per-bin_models.pdf'
        else:
            fname = 'per_bin_auc.pdf'

        plt.savefig(self.output_dir / fname, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {fname}")

    def train_per_bin_models(self, df):
        """Train separate XGBoost models per ET bin"""
        print("\n  Training per-bin models...")

        et_bin_edges = self.per_bin_config.get('et_bin_edges', [5, 10, 15, 20, 25, 30, 35, 40])
        et_bin_labels = self.per_bin_config.get('et_bin_labels',
                                                 ["5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40"])

        params = self.config['model']['params'].copy()
        seed = self.config['training']['random_seed']
        test_size = self.config['training']['test_size']

        for i, label in enumerate(et_bin_labels):
            if i >= len(et_bin_edges) - 1:
                break
            et_low = et_bin_edges[i]
            et_high = et_bin_edges[i + 1]

            # Select data in this ET bin
            df_bin = df[(df['cluster_Et'] >= et_low) & (df['cluster_Et'] < et_high)]

            if len(df_bin) < 50:
                print(f"    {label} GeV: insufficient data (n={len(df_bin)}), skipping")
                self.per_bin_auc_binned[label] = np.nan
                continue

            n_sig = (df_bin['label'] == 1).sum()
            n_bkg = (df_bin['label'] == 0).sum()

            if n_sig < 10 or n_bkg < 10:
                print(f"    {label} GeV: insufficient signal/background (sig={n_sig}, bkg={n_bkg}), skipping")
                self.per_bin_auc_binned[label] = np.nan
                continue

            X_bin = df_bin[self.feature_cols].values
            y_bin = df_bin['label'].values

            # Split into train/test
            X_train, X_test, y_train, y_test = train_test_split(
                X_bin, y_bin, test_size=test_size, random_state=seed, stratify=y_bin
            )

            # Compute weights
            w_signal = len(y_train) / (2.0 * np.sum(y_train == 1)) if np.sum(y_train == 1) > 0 else 1.0
            w_background = len(y_train) / (2.0 * np.sum(y_train == 0)) if np.sum(y_train == 0) > 0 else 1.0
            weights = np.where(y_train == 1, w_signal, w_background)

            # Train model
            model = xgb.XGBClassifier(
                **params,
                random_state=seed,
                use_label_encoder=False
            )
            model.fit(X_train, y_train, sample_weight=weights, verbose=False)

            # Evaluate
            y_pred = model.predict_proba(X_test)[:, 1]
            bin_auc = roc_auc_score(y_test, y_pred)

            self.per_bin_models[label] = model
            self.per_bin_auc_binned[label] = bin_auc

            print(f"    {label} GeV: AUC = {bin_auc:.4f} (train={len(X_train)}, test={len(X_test)})")

        # Plot per-bin model AUC
        self.plot_per_bin_auc(self.per_bin_auc_binned, "(Per-bin Models)")

    def plot_auc_comparison(self, single_model_auc, per_bin_model_auc):
        """Plot comparison of single model vs per-bin models AUC"""
        print("\n  Plotting AUC comparison...")

        et_bin_edges = self.per_bin_config.get('et_bin_edges', [5, 10, 15, 20, 25, 30, 35, 40])
        et_bin_labels = self.per_bin_config.get('et_bin_labels',
                                                 ["5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40"])

        # Calculate bin centers
        bin_centers = []
        for i in range(len(et_bin_labels)):
            if i < len(et_bin_edges) - 1:
                bin_centers.append((et_bin_edges[i] + et_bin_edges[i + 1]) / 2.0)

        # Extract AUC values
        single_auc = [single_model_auc.get(l, np.nan) for l in et_bin_labels]
        binned_auc = [per_bin_model_auc.get(l, np.nan) for l in et_bin_labels]

        plt.figure(figsize=(10, 6))

        # Plot lines with markers
        valid_single = [(c, a) for c, a in zip(bin_centers, single_auc) if not np.isnan(a)]
        valid_binned = [(c, a) for c, a in zip(bin_centers, binned_auc) if not np.isnan(a)]

        if valid_single:
            x_s, y_s = zip(*valid_single)
            plt.plot(x_s, y_s, 'o-', color='blue', markersize=8, linewidth=2,
                     label='Single Model')

        if valid_binned:
            x_b, y_b = zip(*valid_binned)
            plt.plot(x_b, y_b, 's--', color='red', markersize=8, linewidth=2,
                     label='Per-bin Models')

        plt.xlabel('ET Bin Center [GeV]', fontsize=12)
        plt.ylabel('AUC', fontsize=12)
        plt.title('AUC Comparison: Single Model vs Per-bin Models', fontsize=14)
        plt.legend(fontsize=11)
        plt.grid(alpha=0.3)
        plt.ylim(0.5, 1.05)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'auc_comparison.pdf', bbox_inches='tight')
        plt.close()
        print(f"  Saved: auc_comparison.pdf")

    def save_model(self):
        """Save model in multiple formats"""
        print("\n[7/7] Saving model...")

        model_name = self.config['training']['model_name']

        # 1. Joblib (Python)
        if self.config['model']['save_joblib']:
            joblib_path = self.output_dir / f"{model_name}.joblib"
            joblib.dump(self.model, joblib_path)
            print(f"  Saved joblib: {joblib_path}")

        # 2. ROOT/TMVA (C++)
        if self.config['model']['save_root'] and HAS_ROOT:
            # Follow the same naming convention as the other training pipeline in this repo
            # (see `main_training.py`): write a TMVA-compatible ROOT file.
            root_path = self.output_dir / f"{model_name}_tmva.root"
            try:
                # Ensure deterministic feature naming expected by TMVA (f0, f1, ...)
                booster = self.model.get_booster()
                booster.feature_names = [f"f{i}" for i in range(len(self.feature_cols))]

                # Match the other training code: pass the XGBClassifier directly
                tmva_method_name = self.config.get('model', {}).get('tmva_method_name', 'myBDT')
                ROOT.TMVA.Experimental.SaveXGBoost(
                    self.model,
                    tmva_method_name,
                    str(root_path),
                    num_inputs=len(self.feature_cols)
                )
                print(f"  Saved ROOT: {root_path}")
            except (AttributeError, TypeError) as e:
                print(f"  WARNING: Could not save TMVA ROOT via SaveXGBoost(XGBClassifier,...): {e}")
                print("  Trying fallback SaveXGBoost(booster,...) signature...")
                try:
                    # Some ROOT builds expose a booster-based signature
                    ROOT.TMVA.Experimental.SaveXGBoost(
                        booster,
                        tmva_method_name,
                        str(root_path),
                        num_inputs=len(self.feature_cols)
                    )
                    print(f"  Saved ROOT (fallback): {root_path}")
                except Exception as e2:
                    print(f"  WARNING: Fallback SaveXGBoost(booster,...) also failed: {e2}")
                    print("  Saving XGBoost JSON as a last-resort fallback for C++ inference...")
                    try:
                        json_path = self.output_dir / f"{model_name}_xgb.json"
                        self.model.save_model(str(json_path))
                        print(f"  Saved XGBoost JSON: {json_path}")
                    except Exception as e3:
                        print(f"  ERROR: Could not save XGBoost JSON either: {e3}")
        elif self.config['model']['save_root'] and not HAS_ROOT:
            print("  WARNING: ROOT not available, cannot save .root format")

        # 3. Metadata
        metadata = {
            'model_name': model_name,
            'n_features': len(self.feature_cols),
            'features': self.feature_cols,
            'config': self.config,
            'best_params': self.best_params_,
            'best_val_auc': self.best_val_auc_,
            'per_bin_auc_single_model': self.per_bin_auc_single,
            'per_bin_auc_binned_models': self.per_bin_auc_binned,
            'et_bin_edges': self.per_bin_config.get('et_bin_edges', []),
        }
        metadata_path = self.output_dir / f"{model_name}_metadata.yaml"
        with open(metadata_path, 'w') as f:
            yaml.dump(metadata, f, default_flow_style=False)
        print(f"  Saved metadata: {metadata_path}")

    def run(self):
        """Execute full training pipeline"""
        print("=" * 70)
        print("NPB Score Training Pipeline")
        print("=" * 70)

        df = self.load_data()
        X_train, X_val, X_test, y_train, y_val, y_test, et_train, et_val, et_test = self.split_data(df)

        # Store ET values for use in train_model (hyperparameter tuning retrain)
        self.et_train = et_train
        self.et_val = et_val

        sample_weights = self.compute_weights(y_train, et_values=et_train)
        self.train_model(X_train, y_train, X_val, y_val, sample_weights)
        test_auc, fpr, tpr, y_pred_proba, y_pred = self.evaluate_model(X_test, y_test)
        self.generate_plots(X_test, y_test, test_auc, fpr, tpr, y_pred_proba, y_pred)

        # Store df for per-bin analysis
        self.df = df

        # Plot ET distributions
        if self.per_bin_config.get('plot_et_distributions', True):
            self.plot_et_distributions(df)

        # Plot ET reweighting QA if enabled
        if self.config['training'].get('reweight_et', False):
            self.plot_et_reweighting_qa(et_train, y_train, sample_weights)

        # Per-bin AUC calculation (single model)
        if self.per_bin_config.get('calculate_per_bin_auc', True):
            self.per_bin_auc_single = self.calculate_per_bin_auc(X_test, y_test, y_pred_proba)
            self.plot_per_bin_auc(self.per_bin_auc_single, "(Single Model)")

        # Per-bin model training (optional)
        if self.per_bin_config.get('train_per_bin_models', False):
            self.train_per_bin_models(df)
            self.plot_auc_comparison(self.per_bin_auc_single, self.per_bin_auc_binned)

        self.save_model()

        print("\n" + "=" * 70)
        print("Training Complete!")
        print("=" * 70)
        print(f"Output directory: {self.output_dir}/")
        print(f"Test AUC: {test_auc:.4f}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Train NPB rejection score BDT")
    parser.add_argument("--config", default="config_npb_training.yaml",
                       help="Path to configuration file")
    parser.add_argument("--per-bin-training", action="store_true",
                       help="Enable per-bin model training for comparison")
    args = parser.parse_args()

    trainer = NPBTrainer(args.config)
    if args.per_bin_training:
        trainer.per_bin_config['train_per_bin_models'] = True
    trainer.run()
