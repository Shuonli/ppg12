"""
Main training pipeline for binned training.
"""
import os
import json
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, roc_curve
from scipy import stats
import joblib

from data_loader import DataLoader
from reweighting import KinematicReweighter
from model_builder import ModelBuilder
from plotting import BinnedTrainingPlotter


class BinnedTrainingPipeline:
    """Main class for binned training pipeline."""
    
    def __init__(self, config_path: str):
        self.config = self._load_config(config_path)
        self._validate_config()
        
        # Initialize components
        self.data_loader = DataLoader(self.config)
        self.reweighter = KinematicReweighter(self.config)
        self.model_builder = ModelBuilder(self.config)
        
        # Training state
        self.trained_pipelines = {}
        self.train_metrics = {}
        self.val_results = {}
        self.correlation_results = {}
        
        # Setup output directory
        self.output_dir = Path(self.config['output']['output_dir'])
        self.output_dir.mkdir(exist_ok=True)
        
        # Set random seed
        np.random.seed(self.config['training']['global_seed'])
    
    def _load_config(self, config_path: str) -> Dict:
        """Load configuration from YAML file."""
        import yaml
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def _validate_config(self) -> None:
        """Validate configuration parameters."""
        required_keys = ['training', 'model', 'binning', 'data', 'reweighting', 'output']
        for key in required_keys:
            if key not in self.config:
                raise ValueError(f"Missing required config key: {key}")
        
        # Validate binning configuration based on mode
        binning_mode = self.config['binning'].get('mode', 'pt')
        
        if binning_mode == 'pt':
            # Use pt_edges/pt_labels if available, otherwise fall back to edges/labels
            if 'pt_edges' in self.config['binning'] and 'pt_labels' in self.config['binning']:
                edges = self.config['binning']['pt_edges']
                labels = self.config['binning']['pt_labels']
            else:
                edges = self.config['binning']['edges']
                labels = self.config['binning']['labels']
            
            if len(edges) < 2:
                raise ValueError("pT binning edges must have at least 2 values")
            if len(edges) != len(labels) + 1:
                raise ValueError("Number of pT bin labels must be one less than number of edges")
                
        elif binning_mode == 'vertex':
            edges = self.config['binning']['vertex_edges']
            labels = self.config['binning']['vertex_labels']
            
            if len(edges) < 2:
                raise ValueError("Vertex binning edges must have at least 2 values")
            if len(edges) != len(labels) + 1:
                raise ValueError("Number of vertex bin labels must be one less than number of edges")
        else:
            raise ValueError(f"Invalid binning mode: {binning_mode}. Must be 'pt' or 'vertex'")
    
    def _get_binning_config(self) -> Tuple[List[float], List[str], str, str]:
        """Get binning configuration based on mode."""
        binning_mode = self.config['binning'].get('mode', 'pt')
        
        if binning_mode == 'pt':
            # Use pt_edges/pt_labels if available, otherwise fall back to edges/labels
            if 'pt_edges' in self.config['binning'] and 'pt_labels' in self.config['binning']:
                edges = self.config['binning']['pt_edges']
                labels = self.config['binning']['pt_labels']
            else:
                edges = self.config['binning']['edges']
                labels = self.config['binning']['labels']
            bin_column = self.config['data']['pt_column']
            
        elif binning_mode == 'vertex':
            edges = self.config['binning']['vertex_edges']
            labels = self.config['binning']['vertex_labels']
            bin_column = 'vertexz'
            
        return edges, labels, binning_mode, bin_column
    
    def _split_train_val_test_with_features(self, df: pd.DataFrame, features: List[str]) -> Tuple:
        """Split data into train/validation/test sets with features."""
        X = df[features].copy()
        y = df[self.config['data']['label_column']].astype(int).values
        w = df["weight"].to_numpy(dtype=np.float32)

        strat = y if self.config['training']['stratify'] else None

        # Get split ratios from config
        train_size = self.config['training']['train_size']
        val_size = self.config['training'].get('val_size', 0.1)  # Default 10% validation
        test_size = self.config['training'].get('test_size', 0.2)  # Default 20% test

        # Validate split ratios
        total_size = train_size + val_size + test_size
        if not np.isclose(total_size, 1.0):
            raise ValueError(f"train_size + val_size + test_size must equal 1.0, got {total_size}")

        # First split: separate out test set
        X_trainval, X_test, y_trainval, y_test, w_trainval, w_test = train_test_split(
            X, y, w,
            test_size=test_size,
            random_state=self.config['training']['global_seed'],
            stratify=strat
        )

        # Second split: separate train and validation from the remaining data
        # Adjust validation size relative to trainval size
        val_size_adjusted = val_size / (train_size + val_size)
        strat_trainval = y_trainval if self.config['training']['stratify'] else None

        X_train, X_val, y_train, y_val, w_train, w_val = train_test_split(
            X_trainval, y_trainval, w_trainval,
            test_size=val_size_adjusted,
            random_state=self.config['training']['global_seed'] + 1,  # Different seed
            stratify=strat_trainval
        )

        return X_train, X_val, X_test, y_train, y_val, y_test, w_train, w_val, w_test

    def _tune_hyperparameters(self, X_trainval, y_trainval, w_trainval, features):
        """Run Optuna hyperparameter tuning with k-fold cross-validation.

        Args:
            X_trainval: Combined training and validation features (80% of data)
            y_trainval: Combined labels
            w_trainval: Combined sample weights
            features: List of feature names

        Returns:
            dict: Best hyperparameters or None (fallback to config defaults)
        """
        import optuna
        from optuna.samplers import TPESampler

        print("\n=== Phase 1: Hyperparameter Tuning with Optuna ===")

        # Configuration
        n_trials = self.config['training'].get('n_trials', 50)
        n_folds = self.config['training'].get('n_cv_folds', 5)
        timeout = self.config['training'].get('tuning_timeout_seconds', None)
        tuning_fraction = self.config['training'].get('tuning_data_fraction', 1.0)

        # Subsample data for faster tuning if requested
        if tuning_fraction < 1.0:
            from sklearn.model_selection import train_test_split

            n_total = len(y_trainval)
            n_tuning = int(n_total * tuning_fraction)

            print(f"Using {tuning_fraction*100:.1f}% of data for tuning ({n_tuning}/{n_total} samples)")

            # Stratified sampling to maintain class balance
            X_tune, _, y_tune, _, w_tune, _ = train_test_split(
                X_trainval, y_trainval, w_trainval,
                train_size=tuning_fraction,
                stratify=y_trainval,
                random_state=self.config['training']['global_seed']
            )

            X_trainval, y_trainval, w_trainval = X_tune, y_tune, w_tune
            print(f"Tuning data: Signal={np.sum(y_tune==1)}, Background={np.sum(y_tune==0)}")
        else:
            print(f"Using 100% of data for tuning ({len(y_trainval)} samples)")

        # Create Optuna study
        study = optuna.create_study(
            direction='maximize',
            sampler=TPESampler(seed=self.config['training']['global_seed'])
        )

        # Create objective function
        objective = self._create_optuna_objective(
            X_trainval, y_trainval, w_trainval, features, n_folds
        )

        # Optimize
        study.optimize(objective, n_trials=n_trials, timeout=timeout, show_progress_bar=True)

        # Get best params
        best_params = study.best_params
        best_cv_auc = study.best_value

        print(f"\nBest hyperparameters: {best_params}")
        print(f"Best cross-validation AUC: {best_cv_auc:.4f}")

        # Save results
        self._save_tuning_results(study)

        # Fallback check: if tuning performs worse than baseline, return None
        baseline_auc = self.config.get('tuning', {}).get('baseline_auc_threshold', 0.50)
        if best_cv_auc < baseline_auc:
            print(f"WARNING: Best CV AUC ({best_cv_auc:.4f}) below baseline ({baseline_auc:.4f})")
            print("Falling back to config defaults")
            return None

        return best_params

    def _create_optuna_objective(self, X, y, w, features, n_folds):
        """Create objective function for Optuna optimization.

        Returns:
            callable: Objective function that returns mean CV AUC
        """
        def objective(trial):
            # Suggest hyperparameters from configured ranges
            tuning_config = self.config.get('tuning', {})

            params = {
                'learning_rate': trial.suggest_float(
                    'learning_rate',
                    tuning_config.get('learning_rate_min', 0.001),
                    tuning_config.get('learning_rate_max', 0.3),
                    log=True
                ),
                'max_depth': trial.suggest_int(
                    'max_depth',
                    tuning_config.get('max_depth_min', 3),
                    tuning_config.get('max_depth_max', 10)
                ),
                'n_estimators': trial.suggest_int(
                    'n_estimators',
                    tuning_config.get('n_estimators_min', 100),
                    tuning_config.get('n_estimators_max', 1000),
                    step=tuning_config.get('n_estimators_step', 50)
                ),
                'subsample': trial.suggest_float(
                    'subsample',
                    tuning_config.get('subsample_min', 0.5),
                    tuning_config.get('subsample_max', 1.0)
                ),
                'reg_alpha': trial.suggest_float(
                    'reg_alpha',
                    tuning_config.get('reg_alpha_min', 1e-8),
                    tuning_config.get('reg_alpha_max', 10.0),
                    log=True
                ),
                'reg_lambda': trial.suggest_float(
                    'reg_lambda',
                    tuning_config.get('reg_lambda_min', 1e-8),
                    tuning_config.get('reg_lambda_max', 10.0),
                    log=True
                ),
                'colsample_bytree': trial.suggest_float(
                    'colsample_bytree',
                    tuning_config.get('colsample_bytree_min', 0.5),
                    tuning_config.get('colsample_bytree_max', 1.0)
                ),
            }

            # Run k-fold cross-validation
            cv_scores = self._stratified_kfold_cv(X, y, w, features, params, n_folds)

            return np.mean(cv_scores)

        return objective

    def _stratified_kfold_cv(self, X, y, w, features, params, n_folds):
        """Perform stratified k-fold cross-validation with sample weights.

        Returns:
            list: AUC scores for each fold
        """
        from sklearn.model_selection import StratifiedKFold

        skf = StratifiedKFold(
            n_splits=n_folds,
            shuffle=True,
            random_state=self.config['training']['global_seed']
        )

        cv_scores = []

        for fold_idx, (train_idx, val_idx) in enumerate(skf.split(X, y)):
            # Split fold
            X_train_fold = X.iloc[train_idx] if hasattr(X, 'iloc') else X[train_idx]
            X_val_fold = X.iloc[val_idx] if hasattr(X, 'iloc') else X[val_idx]
            y_train_fold = y[train_idx]
            y_val_fold = y[val_idx]
            w_train_fold = w[train_idx]
            w_val_fold = w[val_idx]

            # Create model builder with trial hyperparameters
            temp_model_builder = ModelBuilder(self.config)
            temp_model_builder.update_params(params)

            # Build and train pipeline
            pipeline = temp_model_builder.build_pipeline()
            pipeline.fit(X_train_fold, y_train_fold, clf__sample_weight=w_train_fold)

            # Evaluate on validation fold
            if hasattr(pipeline.named_steps["clf"], "predict_proba"):
                y_pred = pipeline.predict_proba(X_val_fold)[:, 1]
            else:
                y_pred = pipeline.decision_function(X_val_fold)

            # Calculate weighted AUC
            fold_auc = roc_auc_score(y_val_fold, y_pred, sample_weight=w_val_fold)
            cv_scores.append(fold_auc)

        return cv_scores

    def _save_tuning_results(self, study):
        """Save Optuna study results to files.

        Saves:
            - Best hyperparameters to JSON
            - All trials to CSV
            - Optimization history plot
        """
        import json
        from pathlib import Path

        output_dir = Path(self.output_dir)

        # Save best params as JSON
        best_params_path = output_dir / "best_hyperparameters.json"
        with open(best_params_path, 'w') as f:
            json.dump({
                'best_params': study.best_params,
                'best_value': study.best_value,
                'n_trials': len(study.trials),
            }, f, indent=2)
        print(f"Saved best hyperparameters: {best_params_path}")

        # Save all trials as CSV
        trials_df = study.trials_dataframe()
        trials_path = output_dir / "optuna_trials.csv"
        trials_df.to_csv(trials_path, index=False)
        print(f"Saved all trials: {trials_path}")

        # Save optimization history plot
        try:
            import plotly
            fig = optuna.visualization.plot_optimization_history(study)
            plot_path = output_dir / "optuna_optimization_history.html"
            fig.write_html(str(plot_path))
            print(f"Saved optimization history: {plot_path}")
        except Exception as e:
            print(f"Could not save optimization plot: {e}")

    def _calculate_correlations(self, iso_val: np.ndarray, y_proba: np.ndarray, 
                               w_val: np.ndarray, mask: np.ndarray) -> Dict:
        """Calculate correlations between isoET and BDT score."""
        iso_c = iso_val[mask]
        proba_c = y_proba[mask]
        w_c = w_val[mask]
        
        # Weighted Pearson correlation
        def weighted_pearson(x, y, w):
            mx = np.average(x, weights=w)
            my = np.average(y, weights=w)
            cov = np.average((x - mx) * (y - my), weights=w)
            sx = np.sqrt(np.average((x - mx)**2, weights=w))
            sy = np.sqrt(np.average((y - my)**2, weights=w))
            if sx * sy == 0:
                return np.nan
            return cov / (sx * sy)
        
        pearson_r = weighted_pearson(iso_c, proba_c, w_c) if len(iso_c) > 1 else np.nan
        spearman_rho, spearman_p = stats.spearmanr(iso_c, proba_c) if len(iso_c) > 1 else (np.nan, np.nan)
        
        return {
            "pearson_r": float(pearson_r) if np.isfinite(pearson_r) else np.nan,
            "pearson_p": np.nan,  # p-value for weighted correlation is complex
            "spearman_rho": float(spearman_rho) if np.isfinite(spearman_rho) else np.nan,
            "spearman_p": float(spearman_p) if np.isfinite(spearman_p) else np.nan,
            "n_pairs": int(np.sum(mask)),
        }
    
    def analyze_feature_correlations(self, df_all: pd.DataFrame) -> None:
        """Analyze correlation between isoET and training features for background samples."""
        print("\\n=== Analyzing Feature Correlations (Background Samples Only) ===")
        
        # Filter to background samples only
        background_mask = df_all[self.config['data']['label_column']] == 0
        df_background = df_all[background_mask]
        
        if len(df_background) == 0:
            print("Warning: No background samples found for correlation analysis")
            return
        
        print(f"Analyzing {len(df_background):,} background samples")
        
        # Get training features
        features = self.data_loader.autodetect_features(df_background)
        iso_col = self.config['data']['iso_column']
        
        if iso_col not in df_background.columns:
            print(f"Warning: {iso_col} column not found in data")
            return
        
        if not features:
            print("Warning: No training features found")
            return
        
        # Calculate correlations between isoET and each training feature
        correlations = {}
        iso_values = df_background[iso_col].values
        
        print(f"\\nCorrelations between {iso_col} and training features:")
        print("-" * 60)
        print(f"{'Feature':<25} {'Pearson r':<12} {'P-value':<12} {'Samples':<10}")
        print("-" * 60)
        
        for feature in features:
            if feature == iso_col:
                continue  # Skip self-correlation
                
            if feature not in df_background.columns:
                continue
                
            feature_values = df_background[feature].values
            
            # Remove any rows with NaN values
            valid_mask = ~(np.isnan(iso_values) | np.isnan(feature_values))
            if np.sum(valid_mask) < 10:  # Need at least 10 points
                correlations[feature] = {'r': np.nan, 'p': np.nan, 'n': 0}
                continue
            
            iso_clean = iso_values[valid_mask]
            feature_clean = feature_values[valid_mask]
            
            # Calculate Pearson correlation
            try:
                r, p = stats.pearsonr(iso_clean, feature_clean)
                correlations[feature] = {'r': r, 'p': p, 'n': len(iso_clean)}
                
                # Format and print
                r_str = f"{r:.4f}" if not np.isnan(r) else "N/A"
                p_str = f"{p:.2e}" if not np.isnan(p) and p > 0 else "< 1e-16" if not np.isnan(p) else "N/A"
                print(f"{feature:<25} {r_str:<12} {p_str:<12} {len(iso_clean):<10}")
                
            except Exception as e:
                print(f"{feature:<25} Error: {str(e)[:30]}")
                correlations[feature] = {'r': np.nan, 'p': np.nan, 'n': 0}
        
        print("-" * 60)
        
        # Store correlation results for later use
        self.feature_correlations = correlations
        
        # Generate visualization
        self._plot_feature_correlations(correlations, df_background, features, iso_col)
        
        # Save correlation results
        self._save_feature_correlations(correlations, iso_col)
    
    def _plot_feature_correlations(self, correlations: Dict, df_background: pd.DataFrame, 
                                  features: List[str], iso_col: str) -> None:
        """Create visualization of feature correlations."""
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            from scipy import stats
            
            # Prepare data for plotting
            features_with_corr = [f for f in features if f != iso_col and f in correlations]
            correlations_values = [correlations[f]['r'] for f in features_with_corr if not np.isnan(correlations[f]['r'])]
            features_clean = [f for f in features_with_corr if not np.isnan(correlations[f]['r'])]
            
            if not correlations_values:
                print("No valid correlations to plot")
                return
            
            # Create figure with three subplots
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 6))
            
            # Plot 1: Bar chart of correlations
            abs_corr = [abs(r) for r in correlations_values]
            colors = ['red' if r < 0 else 'blue' for r in correlations_values]
            
            bars = ax1.barh(range(len(features_clean)), correlations_values, color=colors, alpha=0.7)
            ax1.set_yticks(range(len(features_clean)))
            ax1.set_yticklabels(features_clean)
            ax1.set_xlabel('Pearson Correlation Coefficient')
            ax1.set_title(f'Background: {iso_col} vs Training Features')
            ax1.grid(True, alpha=0.3, axis='x')
            ax1.axvline(x=0, color='black', linestyle='-', alpha=0.3)
            
            # Add correlation values on bars
            for i, (bar, r) in enumerate(zip(bars, correlations_values)):
                ax1.text(r + 0.01*np.sign(r), i, f'{r:.3f}', 
                        va='center', ha='left' if r >= 0 else 'right')
            
            # Plot 2: Correlation matrix heatmap for top correlated features
            top_features = sorted(features_clean, key=lambda f: abs(correlations[f]['r']), reverse=True)[:15]
            
            if len(top_features) > 1:
                # Create correlation matrix for top features + isoET
                features_for_matrix = [iso_col] + top_features
                corr_matrix = df_background[features_for_matrix].corr()
                
                # Plot heatmap
                mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
                sns.heatmap(corr_matrix, mask=mask, annot=True, fmt='.3f', 
                           cmap='RdBu_r', center=0, ax=ax2, cbar_kws={'shrink': 0.8})
                ax2.set_title(f'Correlation Matrix: {iso_col} + Top Features')
            else:
                ax2.text(0.5, 0.5, 'Insufficient features for matrix', 
                        ha='center', va='center', transform=ax2.transAxes)
                ax2.set_title('Correlation Matrix')
            
            # Plot 3: 2D correlation plot with profile overlay for most correlated feature
            if top_features:
                most_correlated = top_features[0]
                # Check if y-axis zoom is enabled in config
                zoom_y = self.config.get('plotting', {}).get('zoom_y_axis', True)
                self._plot_2d_correlation_with_profile(ax3, df_background, most_correlated, iso_col, 
                                                     correlations[most_correlated]['r'], zoom_y=zoom_y)
            
            plt.tight_layout()
            
            # Save plot
            if self.config['output']['save_plots']:
                plot_path = self.output_dir / "feature_correlations_background.pdf"
                plt.savefig(plot_path, dpi=300, bbox_inches='tight')
                print(f"Saved correlation plot: {plot_path}")
            
            plt.show()
            
            # Create additional 2D plots for top correlated features
            self._create_2d_correlation_grid(df_background, top_features[:6], iso_col, correlations)
            
        except ImportError:
            print("Matplotlib/Seaborn not available - skipping correlation plots")
        except Exception as e:
            print(f"Error creating correlation plots: {e}")
    
    def _plot_2d_correlation_with_profile(self, ax, df: pd.DataFrame, feature: str, iso_col: str, 
                                        correlation: float, zoom_y: bool = True) -> None:
        """Plot 2D correlation with profile overlay."""
        try:
            # Remove NaN values
            valid_mask = ~(df[feature].isna() | df[iso_col].isna())
            x_data = df.loc[valid_mask, feature].values
            y_data = df.loc[valid_mask, iso_col].values
            
            if len(x_data) < 10:
                ax.text(0.5, 0.5, f'Insufficient data for {feature}', 
                       ha='center', va='center', transform=ax.transAxes)
                return
            
            # Set y-axis limits for better visualization
            if zoom_y:
                # Focus on the most relevant range of isoET values
                y_percentiles = np.percentile(y_data, [0.5, 99.5])  # Use 0.5-99.5% range
                y_min, y_max = y_percentiles
                
                # Add some padding (5% of range)
                y_range = y_max - y_min
                y_padding = y_range * 0.05
                y_min = max(0, y_min - y_padding)  # Don't go below 0 for isoET
                y_max = y_max + y_padding
                
                # Filter data to the zoomed range for better visualization
                zoom_mask = (y_data >= y_min) & (y_data <= y_max)
                x_data_zoom = x_data[zoom_mask]
                y_data_zoom = y_data[zoom_mask]
                
                print(f"  Y-axis zoom: {iso_col} range [{y_min:.2f}, {y_max:.2f}] ({np.sum(zoom_mask)}/{len(y_data)} points)")
            else:
                x_data_zoom = x_data
                y_data_zoom = y_data
                y_min, y_max = None, None
            
            # Create 2D histogram
            ax.hexbin(x_data_zoom, y_data_zoom, gridsize=50, cmap='Blues', alpha=0.7, mincnt=1)
            
            # Set y-axis limits if zooming
            if zoom_y and y_min is not None and y_max is not None:
                ax.set_ylim(y_min, y_max)
            
            # Calculate and plot profile (mean and std) using original data for better statistics
            n_bins = 20
            x_min, x_max = np.percentile(x_data, [1, 99])
            bin_edges = np.linspace(x_min, x_max, n_bins + 1)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            
            profile_means = []
            profile_stds = []
            profile_centers = []
            
            for i in range(n_bins):
                mask = (x_data >= bin_edges[i]) & (x_data < bin_edges[i + 1])
                if i == n_bins - 1:  # Include the last edge
                    mask = (x_data >= bin_edges[i]) & (x_data <= bin_edges[i + 1])
                
                if np.sum(mask) > 5:  # Need at least 5 points per bin
                    y_bin = y_data[mask]
                    profile_means.append(np.mean(y_bin))
                    profile_stds.append(np.std(y_bin))
                    profile_centers.append(bin_centers[i])
            
            if profile_centers:
                profile_centers = np.array(profile_centers)
                profile_means = np.array(profile_means)
                profile_stds = np.array(profile_stds)
                
                # Plot profile with error bars
                ax.errorbar(profile_centers, profile_means, yerr=profile_stds, 
                           color='red', marker='o', markersize=4, linewidth=2, 
                           capsize=3, capthick=1, label='Profile ± 1σ')
                
                # Add correlation coefficient to plot
                ax.text(0.05, 0.95, f'r = {correlation:.3f}', 
                       transform=ax.transAxes, fontsize=12, 
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                # Add zoom info if applicable
                if zoom_y and y_min is not None and y_max is not None:
                    ax.text(0.05, 0.85, f'Y: [{y_min:.1f}, {y_max:.1f}]', 
                           transform=ax.transAxes, fontsize=10, 
                           bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
            
            ax.set_xlabel(feature)
            ax.set_ylabel(iso_col)
            ax.set_title(f'2D Correlation: {feature} vs {iso_col}')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
        except Exception as e:
            print(f"Error creating 2D correlation plot for {feature}: {e}")
    
    def _create_2d_correlation_grid(self, df: pd.DataFrame, top_features: List[str], 
                                  iso_col: str, correlations: Dict) -> None:
        """Create a grid of 2D correlation plots for top features."""
        try:
            import matplotlib.pyplot as plt
            
            if not top_features:
                return
            
            n_features = len(top_features)
            n_cols = min(3, n_features)
            n_rows = (n_features + n_cols - 1) // n_cols
            
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
            if n_features == 1:
                axes = [axes]
            elif n_rows == 1:
                axes = axes.reshape(1, -1)
            
            # Check if y-axis zoom is enabled in config
            zoom_y = self.config.get('plotting', {}).get('zoom_y_axis', True)
            
            for i, feature in enumerate(top_features):
                row = i // n_cols
                col = i % n_cols
                ax = axes[row, col] if n_rows > 1 else axes[col]
                
                correlation = correlations[feature]['r']
                self._plot_2d_correlation_with_profile(ax, df, feature, iso_col, correlation, zoom_y=zoom_y)
            
            # Hide empty subplots
            for i in range(n_features, n_rows * n_cols):
                row = i // n_cols
                col = i % n_cols
                ax = axes[row, col] if n_rows > 1 else axes[col]
                ax.set_visible(False)
            
            plt.tight_layout()
            
            # Save plot
            if self.config['output']['save_plots']:
                plot_path = self.output_dir / "2d_correlations_background.pdf"
                plt.savefig(plot_path, dpi=300, bbox_inches='tight')
                print(f"Saved 2D correlation grid: {plot_path}")
            
            plt.show()
            
        except Exception as e:
            print(f"Error creating 2D correlation grid: {e}")
    
    def _save_feature_correlations(self, correlations: Dict, iso_col: str) -> None:
        """Save feature correlation results to CSV."""
        try:
            # Create DataFrame from correlations
            corr_data = []
            for feature, data in correlations.items():
                corr_data.append({
                    'feature': feature,
                    'iso_column': iso_col,
                    'pearson_r': data['r'],
                    'p_value': data['p'],
                    'n_samples': data['n']
                })
            
            df_corr = pd.DataFrame(corr_data)
            df_corr = df_corr.sort_values('pearson_r', key=abs, ascending=False)
            
            # Save to CSV
            corr_path = self.output_dir / "feature_correlations_background.csv"
            df_corr.to_csv(corr_path, index=False)
            print(f"Saved correlation results: {corr_path}")
            
        except Exception as e:
            print(f"Error saving correlation results: {e}")
    
    def run_training(self) -> None:
        """Run the complete training pipeline."""
        print("=== Starting Binned Training Pipeline ===")
        
        # Get binning configuration
        bin_edges, bin_labels, binning_mode, bin_column = self._get_binning_config()
        print(f"Using {binning_mode} binning with {len(bin_labels)} bins: {bin_labels}")
        
        # Load and reweight all data globally
        print("--- Loading and reweighting all data globally ---")
        
        all_dfs_unweighted = self.data_loader.get_all_bin_data(bin_edges, bin_labels)
        # Handle single file set mode (returns single DataFrame) vs per-bin mode (returns list)
        if len(all_dfs_unweighted) == 1 and self.config['data']['use_single_file_set']:
            df_all_unweighted = all_dfs_unweighted[0]
        else:
            df_all_unweighted = pd.concat(all_dfs_unweighted, ignore_index=True)
        
        # Apply reweighting BEFORE breaking into bins
        df_all_reweighted = self.reweighter.apply_all_reweighting(df_all_unweighted)
        # Keep globally for plotting later
        self.df_all_unweighted = df_all_unweighted
        self.df_all_reweighted = df_all_reweighted
        
        # DEBUG: Check vertex reweighting configuration and data
        print("--- Global reweighting complete ---")
        print(f"DEBUG: Vertex reweighting enabled: {self.config['reweighting'].get('vertex_reweight', False)}")
        print(f"DEBUG: Data columns available: {list(df_all_unweighted.columns)}")
        print(f"DEBUG: 'vertexz' in data: {'vertexz' in df_all_unweighted.columns}")
        if 'vertexz' in df_all_unweighted.columns:
            print(f"DEBUG: Vertex range: {df_all_unweighted['vertexz'].min():.2f} to {df_all_unweighted['vertexz'].max():.2f} cm")
        
        # Analyze feature correlations before training
        self.analyze_feature_correlations(df_all_reweighted)
        
        if self.config['training']['train_single_model']:
            self._train_single_model(df_all_reweighted, df_all_unweighted, bin_edges, bin_labels, binning_mode, bin_column)
        else:
            self._train_per_bin_models(df_all_reweighted, df_all_unweighted, bin_edges, bin_labels, binning_mode, bin_column)
        
        print("=== Training Complete ===")
    
    def _train_single_model(self, df_all_reweighted: pd.DataFrame,
                           df_all_unweighted: pd.DataFrame,
                           bin_edges: List[float], bin_labels: List[str],
                           binning_mode: str, bin_column: str) -> None:
        """Train a single model for all bins."""
        print(f"\n=== Training a single model for all {binning_mode} bins ===")

        # Get features and split
        features = self.data_loader.autodetect_features(df_all_reweighted)
        if not features:
            raise ValueError("No features found for combined dataset.")

        X_train, X_val, X_test, y_train, y_val, y_test, w_train, w_val, w_test = self._split_train_val_test_with_features(
            df_all_reweighted, features
        )

        # Keep separate references to the validation split for plotting outputs
        X_val_original = X_val.copy()
        y_val_original = y_val.copy()
        w_val_original = w_val.copy()

        print(f"Data split: {len(y_train)} train, {len(y_val)} validation, {len(y_test)} test samples")

        # Phase 1: Hyperparameter tuning (optional)
        if self.config['training'].get('enable_hyperparameter_tuning', False):
            # Combine train and val for tuning
            X_trainval = pd.concat([X_train, X_val], axis=0)
            y_trainval = np.concatenate([y_train, y_val])
            w_trainval = np.concatenate([w_train, w_val])

            best_params = self._tune_hyperparameters(X_trainval, y_trainval, w_trainval, features)

            if best_params is not None:
                print("\n=== Phase 2: Training with Best Hyperparameters ===")
                self.model_builder.update_params(best_params)
                # Retrain on combined train+val data
                X_train, y_train, w_train = X_trainval, y_trainval, w_trainval
            else:
                print("\n=== Phase 2: Training with Default Hyperparameters ===")
        else:
            print("\n=== Training with Default Hyperparameters (tuning disabled) ===")

        # Phase 2: Build and train final pipeline
        pipe = self.model_builder.build_pipeline()
        t0 = time.time()
        pipe.fit(X_train, y_train, clf__sample_weight=w_train)
        t1 = time.time()

        # Calculate training AUC
        if hasattr(pipe.named_steps["clf"], "predict_proba"):
            y_train_proba = pipe.predict_proba(X_train)[:, 1]
        else:
            y_train_proba = pipe.decision_function(X_train)

        auc_train = roc_auc_score(y_train, y_train_proba, sample_weight=w_train) if len(np.unique(y_train)) > 1 else np.nan

        # Store results for each bin using TEST set
        df_val_all = df_all_reweighted.loc[X_val_original.index] if len(X_val_original) else pd.DataFrame(columns=df_all_reweighted.columns)
        df_test_all = df_all_reweighted.loc[X_test.index]

        for i in range(len(bin_edges) - 1):
            bin_min, bin_max = bin_edges[i], bin_edges[i+1]
            bin_label = bin_labels[i]

            # Filter validation and test data for this bin
            val_mask = (df_val_all[bin_column] >= bin_min) & \
                       (df_val_all[bin_column] < bin_max)
            test_mask = (df_test_all[bin_column] >= bin_min) & \
                       (df_test_all[bin_column] < bin_max)
            df_val_bin = df_val_all[val_mask]
            df_test_bin = df_test_all[test_mask]

            if df_test_bin.empty:
                continue

            X_val_bin = df_val_bin[features] if not df_val_bin.empty else pd.DataFrame(columns=features)
            y_val_bin = df_val_bin[self.config['data']['label_column']].astype(int).values if not df_val_bin.empty else np.array([])
            w_val_bin = df_val_bin["weight"].to_numpy(dtype=np.float32) if not df_val_bin.empty else np.array([], dtype=np.float32)

            X_test_bin = df_test_bin[features]
            y_test_bin = df_test_bin[self.config['data']['label_column']].astype(int).values
            w_test_bin = df_test_bin["weight"].to_numpy(dtype=np.float32)

            self.trained_pipelines[bin_label] = {
                "pipeline": pipe,
                "features": features,
                "X_val": X_val_bin,
                "y_val": y_val_bin,
                "w_val": w_val_bin,
                "val_index": X_val_bin.index,
                "X_test": X_test_bin,
                "y_test": y_test_bin,
                "w_test": w_test_bin,
                "test_index": X_test_bin.index,
                "df_bin": df_all_reweighted,
                "df_bin_unweighted": df_all_unweighted,
            }

            # Store metrics with appropriate labels
            metrics_entry = {
                "train_time_sec": round(t1 - t0, 3),
                "n_train": int(len(y_train)),
                "n_val": int(len(y_val)),
                "n_test": int(len(y_test_bin)),
                "auc_train": float(auc_train) if not np.isnan(auc_train) else np.nan,
                "bin_min": float(bin_min),
                "bin_max": float(bin_max),
            }

            # Add legacy field names for backward compatibility
            if binning_mode == 'pt':
                metrics_entry["et_min"] = float(bin_min)
                metrics_entry["et_max"] = float(bin_max)

            self.train_metrics[bin_label] = metrics_entry
    
    def _train_per_bin_models(self, df_all_reweighted: pd.DataFrame,
                             df_all_unweighted: pd.DataFrame,
                             bin_edges: List[float], bin_labels: List[str],
                             binning_mode: str, bin_column: str) -> None:
        """Train separate models for each bin."""

        # Global hyperparameter tuning (optional) - runs once before binning
        best_params_global = None
        if self.config['training'].get('enable_hyperparameter_tuning', False):
            print("\n=== Phase 1: Global Hyperparameter Tuning (applies to all bins) ===")

            # Use all data for global tuning
            features = self.data_loader.autodetect_features(df_all_reweighted)
            if features:
                X_train, X_val, X_test, y_train, y_val, y_test, w_train, w_val, w_test = \
                    self._split_train_val_test_with_features(df_all_reweighted, features)

                # Combine train and val for tuning
                X_trainval = pd.concat([X_train, X_val], axis=0)
                y_trainval = np.concatenate([y_train, y_val])
                w_trainval = np.concatenate([w_train, w_val])

                best_params_global = self._tune_hyperparameters(X_trainval, y_trainval, w_trainval, features)

                if best_params_global is not None:
                    print(f"\nGlobal tuning complete. Best hyperparameters will be applied to all bins.")
                else:
                    print(f"\nGlobal tuning fallback - using default hyperparameters for all bins.")
            else:
                print("Warning: No features found for global tuning. Using default hyperparameters.")

        print(f"\n=== Phase 2: Training Models Per Bin ===")

        for i in range(len(bin_edges) - 1):
            bin_min, bin_max = bin_edges[i], bin_edges[i+1]
            bin_label = bin_labels[i]
            unit = "GeV" if binning_mode == "pt" else "cm"
            print(f"\n=== Training bin {bin_label} [{bin_min}, {bin_max}) {unit} ===")

            # Filter data for current bin
            bin_mask = (df_all_reweighted[bin_column] >= bin_min) & \
                      (df_all_reweighted[bin_column] < bin_max)
            df_bin = df_all_reweighted[bin_mask]

            unweighted_bin_mask = (df_all_unweighted[bin_column] >= bin_min) & \
                                 (df_all_unweighted[bin_column] < bin_max)
            df_bin_unweighted = df_all_unweighted[unweighted_bin_mask]

            if df_bin.empty:
                print(f"Skipping bin {bin_label} as it is empty after filtering.")
                continue

            # Get features and split
            features = self.data_loader.autodetect_features(df_bin)
            if not features:
                print(f"Skipping bin {bin_label} as no features found.")
                continue

            X_train, X_val, X_test, y_train, y_val, y_test, w_train, w_val, w_test = self._split_train_val_test_with_features(
                df_bin, features
            )

            print(f"Data split: {len(y_train)} train, {len(y_val)} validation, {len(y_test)} test samples")

            # Apply global tuned hyperparameters if available
            if best_params_global is not None:
                self.model_builder.update_params(best_params_global)
                # Combine train and val for training with tuned params
                X_train_final = pd.concat([X_train, X_val], axis=0)
                y_train_final = np.concatenate([y_train, y_val])
                w_train_final = np.concatenate([w_train, w_val])
            else:
                X_train_final = X_train
                y_train_final = y_train
                w_train_final = w_train

            # Build and train pipeline
            pipe = self.model_builder.build_pipeline()
            t0 = time.time()
            pipe.fit(X_train_final, y_train_final, clf__sample_weight=w_train_final)
            t1 = time.time()

            # Calculate training AUC
            if hasattr(pipe.named_steps["clf"], "predict_proba"):
                y_train_proba = pipe.predict_proba(X_train_final)[:, 1]
            else:
                y_train_proba = pipe.decision_function(X_train_final)

            auc_train = roc_auc_score(y_train_final, y_train_proba, sample_weight=w_train_final) if len(np.unique(y_train_final)) > 1 else np.nan

            # Store results using TEST set
            self.trained_pipelines[bin_label] = {
                "pipeline": pipe,
                "features": features,
                "X_val": X_val,
                "y_val": y_val,
                "w_val": w_val,
                "val_index": X_val.index,
                "X_test": X_test,
                "y_test": y_test,
                "w_test": w_test,
                "test_index": X_test.index,
                "df_bin": df_bin,
                "df_bin_unweighted": df_bin_unweighted,
            }

            # Store metrics with appropriate labels
            metrics_entry = {
                "train_time_sec": round(t1 - t0, 3),
                "n_train": int(len(y_train_final)),
                "n_val": int(len(y_val)),
                "n_test": int(len(y_test)),
                "auc_train": float(auc_train) if not np.isnan(auc_train) else np.nan,
                "bin_min": float(bin_min),
                "bin_max": float(bin_max),
            }

            # Add legacy field names for backward compatibility
            if binning_mode == 'pt':
                metrics_entry["et_min"] = float(bin_min)
                metrics_entry["et_max"] = float(bin_max)

            self.train_metrics[bin_label] = metrics_entry
    
    def evaluate_models(self) -> None:
        """Evaluate all trained models on the test set."""
        print("=== Evaluating Models on Test Set ===")

        for bin_label in self.trained_pipelines:
            entry = self.trained_pipelines[bin_label]
            pipe = entry["pipeline"]
            X_test = entry["X_test"]
            y_test = entry["y_test"]
            w_test = entry["w_test"]

            # Get predictions on test set
            if hasattr(pipe.named_steps["clf"], "predict_proba"):
                y_proba = pipe.predict_proba(X_test)[:, 1]
            else:
                y_proba = pipe.decision_function(X_test)

            # Calculate test AUC (unweighted for fair comparison) and ROC curve (weighted)
            if len(np.unique(y_test)) > 1:
                auc_test = roc_auc_score(y_test, y_proba)
                fpr, tpr, thr = roc_curve(y_test, y_proba, sample_weight=w_test)
            else:
                auc_test, fpr, tpr, thr = np.nan, None, None, None

            # Store test results (keeping old key name for backward compatibility)
            self.val_results[bin_label] = {
                "auc_val": float(auc_test) if not np.isnan(auc_test) else np.nan,
                "auc_test": float(auc_test) if not np.isnan(auc_test) else np.nan,
                "fpr": fpr,
                "tpr": tpr,
                "thresholds": thr,
                "n_val": int(len(y_test)),
                "n_test": int(len(y_test)),
            }

            # Calculate correlations on test set
            iso_test = entry["df_bin"].loc[entry["test_index"], self.config['data']['iso_column']].values
            mask = np.isfinite(iso_test) & np.isfinite(y_proba) & np.isfinite(w_test)

            self.correlation_results[bin_label] = self._calculate_correlations(
                iso_test, y_proba, w_test, mask
            )

            print(f"Bin {bin_label}: Test AUC = {auc_test:.4f}" if not np.isnan(auc_test) else f"Bin {bin_label}: Test AUC = N/A")
    
    def generate_plots(self) -> None:
        """Generate all plots."""
        print("=== Generating Plots ===")
        
        # Initialize plotter
        self.plotter = BinnedTrainingPlotter(self.trained_pipelines, self.config)
        
        # DEBUG: Check plotter configuration
        print(f"DEBUG: Plotter initialized with vertex_reweight = {self.config['reweighting'].get('vertex_reweight', False)}")
        
        # Global reweighting QA (since reweighting is done globally)
        if hasattr(self, 'df_all_unweighted') and hasattr(self, 'df_all_reweighted'):
            print("DEBUG: Calling plot_global_reweighting_qa...")
            print(f"DEBUG: df_all_unweighted shape: {self.df_all_unweighted.shape}")
            print(f"DEBUG: df_all_reweighted shape: {self.df_all_reweighted.shape}")
            self.plotter.plot_global_reweighting_qa(self.df_all_unweighted, self.df_all_reweighted)
        else:
            print("DEBUG: Missing df_all_unweighted or df_all_reweighted for QA plots")

        # Generate plots for each bin (score distributions only)
        for bin_label in self.trained_pipelines:
            entry = self.trained_pipelines[bin_label]
            # no per-bin 2D feature profile plots
        
        # Plot BDT score distributions
        self.plotter.plot_bdt_score_distributions()

        # Global summary plots to match original notebook
        self.plotter.plot_auc_bar(self.val_results)
        self.plotter.plot_combined_roc(self.val_results)
        self.plotter.plot_tpr_at_fixed_fpr(self.val_results, target_fpr=0.2)
        
        # Feature importance plots for BDT evaluation
        self.plotter.plot_feature_importance(save_plots=self.config['output']['save_plots'])
        
        # Global correlation plot between isoET and BDT score
        self.plotter.plot_isoet_bdt_correlation_global(save_plots=self.config['output']['save_plots'])

        # Overlay profiles for pt and iso
        self.plotter.overlay_feature_profiles(self.config['data']['pt_column'])
        self.plotter.overlay_feature_profiles(self.config['data']['iso_column'])
        # Also produce global 2D showershape vs BDT score across all bins
        # for pt, iso, and all training features
        features_for_global = [self.config['data']['pt_column'], self.config['data']['iso_column']]
        # include configured features if available
        cfg_feats = self.config['data'].get('features') or []
        for f in cfg_feats:
            if f not in features_for_global:
                features_for_global.append(f)
        for feat in features_for_global:
            self.plotter.plot_feature_profile_global(feat)
    
    def save_results(self) -> None:
        """Save all results."""
        print("=== Saving Results ===")
        
        # Save models
        if self.config['output']['save_models']:
            for bin_label, entry in self.trained_pipelines.items():
                model_path = self.output_dir / f"model_{bin_label}.joblib"
                joblib.dump(entry["pipeline"], model_path)
                print(f"Saved model: {model_path}")
        
        # Save TMVA/ROOT format models
        if self.config['output'].get('save_tmva', False):
            self._save_tmva_models()
        
        # Save metrics
        if self.config['output']['save_metrics']:
            self._save_metrics()
        
        # Save metadata
        metadata = {
            "config": self.config,
            "training_metrics": self.train_metrics,
            "validation_results": self.val_results,
            "correlation_results": self.correlation_results,
        }
        
        metadata_path = self.output_dir / "training_metadata.json"
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2, default=str)
        print(f"Saved metadata: {metadata_path}")
    
    def _save_tmva_models(self) -> None:
        """Save models in TMVA/ROOT format for C++ compatibility."""
        print("=== Saving TMVA/ROOT Models ===")
        
        # Check if ROOT is available
        try:
            import ROOT
            TMVA_AVAILABLE = True
        except ImportError:
            print("Warning: ROOT/PyROOT not available. Skipping TMVA model saving.")
            print("To enable TMVA saving, install ROOT with PyROOT support.")
            return
        
        # Check if we're using XGBoost
        if self.config['model']['classifier'] != 'xgb':
            print(f"Warning: TMVA saving only supported for XGBoost models, not {self.config['model']['classifier']}")
            return
        
        # Save TMVA models for each bin
        for bin_label, entry in self.trained_pipelines.items():
            try:
                pipeline = entry["pipeline"]
                features = entry["features"]
                
                # Extract XGBoost classifier from pipeline
                if hasattr(pipeline, 'named_steps') and 'clf' in pipeline.named_steps:
                    xgb_model = pipeline.named_steps['clf']
                else:
                    print(f"Warning: Could not extract XGBoost classifier from pipeline for bin {bin_label}")
                    continue
                
                # Get the booster and set feature names
                booster = xgb_model.get_booster()
                
                # Set feature names to match expected format (f0, f1, f2, ...)
                feature_names = [f"f{i}" for i in range(len(features))]
                booster.feature_names = feature_names
                
                # Determine model name and path using configurable prefix
                root_prefix = self.config['output'].get('root_file_prefix', 'model')
                
                if self.config['training']['train_single_model']:
                    # For single model, save it once with a generic name
                    model_name = "myBDT"
                    tmva_path = self.output_dir / f"{root_prefix}_single_tmva.root"
                else:
                    # For per-bin models, include bin label
                    model_name = f"myBDT_{bin_label}"
                    tmva_path = self.output_dir / f"{root_prefix}_{bin_label}_tmva.root"
                
                # Save TMVA model
                ROOT.TMVA.Experimental.SaveXGBoost(
                    xgb_model, model_name, str(tmva_path), num_inputs=len(features)
                )
                
                print(f"Saved TMVA model: {tmva_path}")
                
                # If using single model, only save once
                if self.config['training']['train_single_model']:
                    break
                    
            except Exception as e:
                print(f"Error saving TMVA model for bin {bin_label}: {e}")
                continue
    
    def _save_metrics(self) -> None:
        """Save metrics to CSV."""
        rows = []
        bin_labels = self.config['binning']['labels']
        
        for bin_label in bin_labels:
            row = {"bin": bin_label}
            
            # Training metrics
            row.update(self.train_metrics.get(bin_label, {}))
            
            # Test AUC (keeping auc_val key for backward compatibility, but also adding auc_test)
            if bin_label in self.val_results:
                row["auc_test"] = self.val_results[bin_label]["auc_test"]
                row["auc_val"] = self.val_results[bin_label]["auc_val"]  # For backward compatibility
                # n_test is already in train_metrics
            else:
                row["auc_test"] = np.nan
                row["auc_val"] = np.nan
            
            # Correlations
            if bin_label in self.correlation_results:
                row.update(self.correlation_results[bin_label])
            else:
                row.update({
                    "pearson_r": np.nan,
                    "pearson_p": np.nan,
                    "spearman_rho": np.nan,
                    "spearman_p": np.nan,
                    "n_pairs": 0,
                })
            
            rows.append(row)
        
        metrics_df = pd.DataFrame(rows)
        
        # Handle different binning modes (pt vs vertex)
        binning_mode = self.config['binning'].get('mode', 'pt')
        
        # Define columns based on binning mode
        if binning_mode == 'vertex':
            # For vertex binning, use generic bin_min/bin_max or fallback to et_min/et_max if available
            if 'bin_min' in metrics_df.columns and 'bin_max' in metrics_df.columns:
                min_col, max_col = 'bin_min', 'bin_max'
            else:
                min_col, max_col = 'et_min', 'et_max'  # Fallback for backward compatibility
        else:
            # For pT binning, use et_min/et_max or fallback to bin_min/bin_max
            if 'et_min' in metrics_df.columns and 'et_max' in metrics_df.columns:
                min_col, max_col = 'et_min', 'et_max'
            else:
                min_col, max_col = 'bin_min', 'bin_max'  # Fallback
        
        # Select columns that exist in the dataframe
        base_columns = ["bin", min_col, max_col, "n_train", "n_val", "n_test", "auc_train", "auc_test", "auc_val",
                       "pearson_r", "pearson_p", "spearman_rho", "spearman_p", "n_pairs", "train_time_sec"]
        
        # Only include columns that actually exist
        available_columns = [col for col in base_columns if col in metrics_df.columns]
        metrics_df = metrics_df[available_columns].sort_values("bin")
        
        metrics_path = self.output_dir / "per_bin_metrics.csv"
        metrics_df.to_csv(metrics_path, index=False)
        print(f"Saved metrics: {metrics_path}")
        
        return metrics_df


def main():
    """Main function to run the training pipeline."""
    config_path = "config.yaml"
    
    if not os.path.exists(config_path):
        print(f"Configuration file not found: {config_path}")
        return
    
    # Create and run pipeline
    pipeline = BinnedTrainingPipeline(config_path)
    
    try:
        pipeline.run_training()
        pipeline.evaluate_models()
        pipeline.generate_plots()
        pipeline.save_results()
        print("=== Pipeline completed successfully ===")
    except Exception as e:
        print(f"Pipeline failed with error: {e}")
        raise


if __name__ == "__main__":
    main()
