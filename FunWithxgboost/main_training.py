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
        
        if len(self.config['binning']['edges']) < 2:
            raise ValueError("Binning edges must have at least 2 values")
        
        if len(self.config['binning']['edges']) != len(self.config['binning']['labels']) + 1:
            raise ValueError("Number of bin labels must be one less than number of edges")
    
    def _split_train_val_with_features(self, df: pd.DataFrame, features: List[str]) -> Tuple:
        """Split data into train/validation sets with features."""
        X = df[features].copy()
        y = df[self.config['data']['label_column']].astype(int).values
        w = df["weight"].to_numpy(dtype=np.float32)
        
        strat = y if self.config['training']['stratify'] else None
        
        return train_test_split(
            X, y, w, 
            train_size=self.config['training']['train_size'], 
            random_state=self.config['training']['global_seed'], 
            stratify=strat
        )
    
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
    
    def run_training(self) -> None:
        """Run the complete training pipeline."""
        print("=== Starting Binned Training Pipeline ===")
        
        # Load and reweight all data globally
        print("--- Loading and reweighting all data globally ---")
        bin_edges = self.config['binning']['edges']
        bin_labels = self.config['binning']['labels']
        
        all_dfs_unweighted = self.data_loader.get_all_bin_data(bin_edges, bin_labels)
        df_all_unweighted = pd.concat(all_dfs_unweighted, ignore_index=True)
        
        # Apply reweighting
        df_all_reweighted = self.reweighter.apply_all_reweighting(df_all_unweighted)
        print("--- Global reweighting complete ---")
        
        if self.config['training']['train_single_model']:
            self._train_single_model(df_all_reweighted, df_all_unweighted, bin_edges, bin_labels)
        else:
            self._train_per_bin_models(df_all_reweighted, df_all_unweighted, bin_edges, bin_labels)
        
        print("=== Training Complete ===")
    
    def _train_single_model(self, df_all_reweighted: pd.DataFrame, 
                           df_all_unweighted: pd.DataFrame, 
                           bin_edges: List[float], bin_labels: List[str]) -> None:
        """Train a single model for all bins."""
        print("\n=== Training a single model for all bins ===")
        
        # Get features and split
        features = self.data_loader.autodetect_features(df_all_reweighted)
        if not features:
            raise ValueError("No features found for combined dataset.")
        
        X_train, X_val, y_train, y_val, w_train, w_val = self._split_train_val_with_features(
            df_all_reweighted, features
        )
        
        # Build and train pipeline
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
        
        # Store results for each bin
        df_val_all = df_all_reweighted.loc[X_val.index]
        
        for i in range(len(bin_edges) - 1):
            et_min, et_max = bin_edges[i], bin_edges[i+1]
            bin_label = bin_labels[i]
            
            # Filter validation data for this bin
            val_mask = (df_val_all[self.config['data']['pt_column']] >= et_min) & \
                      (df_val_all[self.config['data']['pt_column']] < et_max)
            df_val_bin = df_val_all[val_mask]
            
            if df_val_bin.empty:
                continue
            
            X_val_bin = df_val_bin[features]
            y_val_bin = df_val_bin[self.config['data']['label_column']].astype(int).values
            w_val_bin = df_val_bin["weight"].to_numpy(dtype=np.float32)
            
            self.trained_pipelines[bin_label] = {
                "pipeline": pipe,
                "features": features,
                "X_val": X_val_bin,
                "y_val": y_val_bin,
                "w_val": w_val_bin,
                "val_index": X_val_bin.index,
                "df_bin": df_all_reweighted,
                "df_bin_unweighted": df_all_unweighted,
            }
            
            self.train_metrics[bin_label] = {
                "train_time_sec": round(t1 - t0, 3),
                "n_train": int(len(y_train)),
                "n_val": int(len(y_val_bin)),
                "auc_train": float(auc_train) if not np.isnan(auc_train) else np.nan,
                "et_min": float(et_min),
                "et_max": float(et_max),
            }
    
    def _train_per_bin_models(self, df_all_reweighted: pd.DataFrame, 
                             df_all_unweighted: pd.DataFrame, 
                             bin_edges: List[float], bin_labels: List[str]) -> None:
        """Train separate models for each bin."""
        for i in range(len(bin_edges) - 1):
            et_min, et_max = bin_edges[i], bin_edges[i+1]
            bin_label = bin_labels[i]
            print(f"\n=== Training bin {bin_label} [{et_min}, {et_max}) GeV ===")
            
            # Filter data for current bin
            bin_mask = (df_all_reweighted[self.config['data']['pt_column']] >= et_min) & \
                      (df_all_reweighted[self.config['data']['pt_column']] < et_max)
            df_bin = df_all_reweighted[bin_mask]
            
            unweighted_bin_mask = (df_all_unweighted[self.config['data']['pt_column']] >= et_min) & \
                                 (df_all_unweighted[self.config['data']['pt_column']] < et_max)
            df_bin_unweighted = df_all_unweighted[unweighted_bin_mask]
            
            if df_bin.empty:
                print(f"Skipping bin {bin_label} as it is empty after filtering.")
                continue
            
            # Get features and split
            features = self.data_loader.autodetect_features(df_bin)
            if not features:
                print(f"Skipping bin {bin_label} as no features found.")
                continue
            
            X_train, X_val, y_train, y_val, w_train, w_val = self._split_train_val_with_features(
                df_bin, features
            )
            
            # Build and train pipeline
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
            
            # Store results
            self.trained_pipelines[bin_label] = {
                "pipeline": pipe,
                "features": features,
                "X_val": X_val,
                "y_val": y_val,
                "w_val": w_val,
                "val_index": X_val.index,
                "df_bin": df_bin,
                "df_bin_unweighted": df_bin_unweighted,
            }
            
            self.train_metrics[bin_label] = {
                "train_time_sec": round(t1 - t0, 3),
                "n_train": int(len(y_train)),
                "n_val": int(len(y_val)),
                "auc_train": float(auc_train) if not np.isnan(auc_train) else np.nan,
                "et_min": float(et_min),
                "et_max": float(et_max),
            }
    
    def evaluate_models(self) -> None:
        """Evaluate all trained models."""
        print("=== Evaluating Models ===")
        
        for bin_label in self.trained_pipelines:
            entry = self.trained_pipelines[bin_label]
            pipe = entry["pipeline"]
            X_val = entry["X_val"]
            y_val = entry["y_val"]
            w_val = entry["w_val"]
            
            # Get predictions
            if hasattr(pipe.named_steps["clf"], "predict_proba"):
                y_proba = pipe.predict_proba(X_val)[:, 1]
            else:
                y_proba = pipe.decision_function(X_val)
            
            # Calculate validation AUC (unweighted to match original) and ROC curve (weighted)
            if len(np.unique(y_val)) > 1:
                auc_val = roc_auc_score(y_val, y_proba)
                fpr, tpr, thr = roc_curve(y_val, y_proba, sample_weight=w_val)
            else:
                auc_val, fpr, tpr, thr = np.nan, None, None, None
            
            self.val_results[bin_label] = {
                "auc_val": float(auc_val) if not np.isnan(auc_val) else np.nan,
                "fpr": fpr,
                "tpr": tpr,
                "thresholds": thr,
                "n_val": int(len(y_val)),
            }
            
            # Calculate correlations
            iso_val = entry["df_bin"].loc[entry["val_index"], self.config['data']['iso_column']].values
            mask = np.isfinite(iso_val) & np.isfinite(y_proba) & np.isfinite(w_val)
            
            self.correlation_results[bin_label] = self._calculate_correlations(
                iso_val, y_proba, w_val, mask
            )
    
    def generate_plots(self) -> None:
        """Generate all plots."""
        print("=== Generating Plots ===")
        
        # Initialize plotter
        self.plotter = BinnedTrainingPlotter(self.trained_pipelines, self.config)
        
        # Generate plots for each bin
        for bin_label in self.trained_pipelines:
            entry = self.trained_pipelines[bin_label]
            df_bin = entry["df_bin"]
            df_bin_unweighted = entry["df_bin_unweighted"]
            
            # Plot reweighting QA
            self.plotter.plot_reweighting_qa(
                df_bin_unweighted, df_bin, bin_label
            )
            
            # Plot feature profiles
            pipe = entry["pipeline"]
            X_val = entry["X_val"]
            w_val = entry["w_val"]
            
            if hasattr(pipe.named_steps["clf"], "predict_proba"):
                y_proba = pipe.predict_proba(X_val)[:, 1]
            else:
                y_proba = pipe.decision_function(X_val)
            
            # Plot profiles for key features
            for feature in [self.config['data']['pt_column'], self.config['data']['iso_column']]:
                if feature in X_val.columns:
                    self.plotter.plot_feature_profile(
                        df_bin, entry["val_index"], y_proba, feature, bin_label, w_val
                    )

            # Plot profiles for all training features (excluding pt and iso)
            for feature in entry.get("features", []):
                if feature in (self.config['data']['pt_column'], self.config['data']['iso_column']):
                    continue
                if feature in X_val.columns:
                    self.plotter.plot_feature_profile(
                        df_bin, entry["val_index"], y_proba, feature, bin_label, w_val
                    )
        
        # Plot BDT score distributions
        self.plotter.plot_bdt_score_distributions()

        # Global summary plots to match original notebook
        self.plotter.plot_auc_bar(self.val_results)
        self.plotter.plot_combined_roc(self.val_results)
        self.plotter.plot_tpr_at_fixed_fpr(self.val_results, target_fpr=0.2)

        # Overlay profiles for pt and iso
        self.plotter.overlay_feature_profiles(self.config['data']['pt_column'])
        self.plotter.overlay_feature_profiles(self.config['data']['iso_column'])
    
    def save_results(self) -> None:
        """Save all results."""
        print("=== Saving Results ===")
        
        # Save models
        if self.config['output']['save_models']:
            for bin_label, entry in self.trained_pipelines.items():
                model_path = self.output_dir / f"model_{bin_label}.joblib"
                joblib.dump(entry["pipeline"], model_path)
                print(f"Saved model: {model_path}")
        
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
    
    def _save_metrics(self) -> None:
        """Save metrics to CSV."""
        rows = []
        bin_labels = self.config['binning']['labels']
        
        for bin_label in bin_labels:
            row = {"bin": bin_label}
            
            # Training metrics
            row.update(self.train_metrics.get(bin_label, {}))
            
            # Validation AUC
            if bin_label in self.val_results:
                row["auc_val"] = self.val_results[bin_label]["auc_val"]
                row["n_val"] = self.val_results[bin_label]["n_val"]
            else:
                row["auc_val"] = np.nan
                row["n_val"] = 0
            
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
        metrics_df = metrics_df[[
            "bin", "et_min", "et_max", "n_train", "n_val", "auc_train", "auc_val",
            "pearson_r", "pearson_p", "spearman_rho", "spearman_p", "n_pairs", "train_time_sec"
        ]].sort_values("bin")
        
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
