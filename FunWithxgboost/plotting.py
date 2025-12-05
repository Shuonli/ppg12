"""
Plotting utilities for binned training.
"""
import pandas as pd
from typing import Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

try:
    import xgboost as xgb
    XGB_AVAILABLE = True
except ImportError:
    xgb = None
    XGB_AVAILABLE = False


class BinnedTrainingPlotter:
    """Handles all plotting functionality for binned training."""
    
    def __init__(self, trained_pipelines: Dict, config: Dict):
        self.trained_pipelines = trained_pipelines
        self.config = config
        self.plot_config = config['plotting']
        
        # Use correct labels based on binning mode
        binning_mode = config['binning'].get('mode', 'pt')
        if binning_mode == 'vertex':
            self.bin_labels = config['binning']['vertex_labels']
        else:
            self.bin_labels = config['binning']['labels']
            
            
        self.label_col = config['data']['label_column']
        self.pt_col = config['data']['pt_column']
        self.iso_col = config['data']['iso_column']
        self.save_plots = config['output'].get('save_plots', False)
        self.output_dir = config['output']['output_dir']
        
        # Get model prefix for plot annotations
        self.root_file_prefix = config['output'].get('root_file_prefix', 'model')
        self.train_single_model = config['training'].get('train_single_model', True)
        
        # Set matplotlib defaults
        plt.rcParams["figure.figsize"] = tuple(self.plot_config['figure_size'])
        plt.rcParams["axes.grid"] = self.plot_config['grid']
        sns.set_context(self.plot_config['context'])
    
    def _add_sphenix_header(self, ax_or_fig, bin_label: str = None) -> None:
        """Add sPHENIX simulation text and prefix to plot."""
        # Add sPHENIX simulation text and prefix at the top
        sphenix_text = r"$\mathbf{\mathit{sPHENIX}}$ Simulation Internal"
        prefix_text = self.root_file_prefix
        
        # Check if this is a figure or axes object
        if hasattr(ax_or_fig, 'transAxes'):  # It's an axes object
            ax = ax_or_fig
            # Add text in top left of plot area
            ax.text(0.02, 0.98, sphenix_text, transform=ax.transAxes, 
                   fontsize=12, weight='bold', va='top', ha='left')
            ax.text(0.02, 0.93, prefix_text, transform=ax.transAxes, 
                   fontsize=9, va='top', ha='left')
        else:  # It's a figure object
            fig = ax_or_fig
            # Add text at top of figure
            fig.text(0.02, 0.98, sphenix_text, fontsize=12, weight='bold', va='top', ha='left')
            fig.text(0.02, 0.94, prefix_text, fontsize=9, va='top', ha='left')
    
    def plot_reweighting_qa(self, df_before: pd.DataFrame, df_after: pd.DataFrame, 
                           bin_label: str) -> None:
        """Plot reweighting QA plots."""
        fig = plt.figure(figsize=(14, 15))
        gs = fig.add_gridspec(3, 2)
        
        # Eta Before Reweighting
        ax1 = fig.add_subplot(gs[0, 0])
        for label_val, name in zip([0, 1], ["Background", "Signal"]):
            mask = df_before[self.label_col] == label_val
            if mask.any():
                ax1.hist(df_before.loc[mask, "cluster_Eta"], bins=40, range=(-0.7, 0.7), 
                         histtype='step', label=name, density=True)
        ax1.set_title(f"Bin {bin_label}: Eta Before Reweighting")
        ax1.set_xlabel("cluster_Eta")
        ax1.set_ylabel("Density")
        ax1.legend()
        ax1.grid(True)
        
        # Eta After Reweighting
        ax2 = fig.add_subplot(gs[0, 1])
        for label_val, name in zip([0, 1], ["Background", "Signal"]):
            mask = df_after[self.label_col] == label_val
            if mask.any():
                ax2.hist(df_after.loc[mask, "cluster_Eta"], bins=40, range=(-0.7, 0.7),
                         weights=df_after.loc[mask, "weight"],
                         histtype='step', label=name, density=True)
        ax2.set_title(f"Bin {bin_label}: Eta After Reweighting")
        ax2.set_xlabel("cluster_Eta")
        ax2.legend()
        ax2.grid(True)
        
        # Et Before Reweighting
        ax_et1 = fig.add_subplot(gs[1, 0])
        et_min, et_max = df_before[self.pt_col].min(), df_before[self.pt_col].max()
        for label_val, name in zip([0, 1], ["Background", "Signal"]):
            mask = df_before[self.label_col] == label_val
            if mask.any():
                ax_et1.hist(df_before.loc[mask, self.pt_col], bins=40, range=(et_min, et_max),
                            histtype='step', label=name, density=True)
        ax_et1.set_title(f"Bin {bin_label}: Et Before Reweighting")
        ax_et1.set_xlabel(self.pt_col)
        ax_et1.set_ylabel("Density")
        ax_et1.legend()
        ax_et1.grid(True)
        
        # Et After Reweighting
        ax_et2 = fig.add_subplot(gs[1, 1])
        for label_val, name in zip([0, 1], ["Background", "Signal"]):
            mask = df_after[self.label_col] == label_val
            if mask.any():
                ax_et2.hist(df_after.loc[mask, self.pt_col], bins=40, range=(et_min, et_max),
                            weights=df_after.loc[mask, "weight"],
                            histtype='step', label=name, density=True)
        ax_et2.set_title(f"Bin {bin_label}: Et After Reweighting")
        ax_et2.set_xlabel(self.pt_col)
        ax_et2.legend()
        ax_et2.grid(True)
        
        # Class Balance Before (counts)
        ax3 = fig.add_subplot(gs[2, 0])
        counts_before = df_before[self.label_col].value_counts()
        n_sig_before = counts_before.get(1, 0)
        n_bkg_before = counts_before.get(0, 0)
        ax3.bar(["Background", "Signal"], [n_bkg_before, n_sig_before], color=['blue', 'orange'])
        ax3.set_title("Class Counts Before Reweighting")
        ax3.set_ylabel("Number of Events")
        ax3.grid(axis='y')
        for i, v in enumerate([n_bkg_before, n_sig_before]):
            ax3.text(i, v, str(v), ha='center', va='bottom')

        # Class Balance After (weighted yields)
        ax4 = fig.add_subplot(gs[2, 1])
        w_sig_after = df_after.loc[df_after[self.label_col] == 1, "weight"].sum() if "weight" in df_after.columns else float((df_after[self.label_col] == 1).sum())
        w_bkg_after = df_after.loc[df_after[self.label_col] == 0, "weight"].sum() if "weight" in df_after.columns else float((df_after[self.label_col] == 0).sum())
        ax4.bar(["Background", "Signal"], [w_bkg_after, w_sig_after], color=['blue', 'orange'])
        ax4.set_title("Class Yield After Reweighting")
        ax4.set_ylabel("Sum of Weights")
        ax4.grid(axis='y')
        for i, v in enumerate([w_bkg_after, w_sig_after]):
            ax4.text(i, v, f"{v:.1f}", ha='center', va='bottom')
        
        plt.tight_layout()
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig, bin_label)
        
        if self.save_plots:
            fig.savefig(f"{self.output_dir}/{self.root_file_prefix}_reweighting_qa_{bin_label}.pdf", dpi=300, bbox_inches='tight')
        plt.show()

    def plot_global_reweighting_qa(self, df_before: pd.DataFrame, df_after: pd.DataFrame) -> None:
        """Global reweighting QA (eta, ET, vertex, class balance) across full dataset."""
        print("DEBUG: Starting plot_global_reweighting_qa...")
        
        # Check if vertex reweighting is enabled
        vertex_reweight_enabled = self.config.get('reweighting', {}).get('vertex_reweight', False)
        print(f"DEBUG: vertex_reweight_enabled = {vertex_reweight_enabled}")
        print(f"DEBUG: df_before columns = {list(df_before.columns)}")
        print(f"DEBUG: 'vertexz' in df_before = {'vertexz' in df_before.columns}")
        print(f"DEBUG: df_before shape = {df_before.shape}")
        print(f"DEBUG: df_after shape = {df_after.shape}")
        
        # Adjust figure size and grid based on whether vertex reweighting is enabled
        if vertex_reweight_enabled:
            fig = plt.figure(figsize=(14, 20))
            gs = fig.add_gridspec(4, 2)  # Add row for vertex plots
        else:
            fig = plt.figure(figsize=(14, 15))
            gs = fig.add_gridspec(3, 2)

        # Eta Before Reweighting
        ax1 = fig.add_subplot(gs[0, 0])
        for label_val, name in zip([0, 1], ["Background", "Signal"]):
            mask = df_before[self.label_col] == label_val
            if mask.any():
                ax1.hist(df_before.loc[mask, "cluster_Eta"], bins=40, range=(-0.7, 0.7),
                         histtype='step', label=name, density=True)
        ax1.set_title("Global: Eta Before Reweighting")
        ax1.set_xlabel("cluster_Eta"); ax1.set_ylabel("Density"); ax1.legend(); ax1.grid(True)

        # Eta After Reweighting
        ax2 = fig.add_subplot(gs[0, 1])
        for label_val, name in zip([0, 1], ["Background", "Signal"]):
            mask = df_after[self.label_col] == label_val
            if mask.any():
                ax2.hist(df_after.loc[mask, "cluster_Eta"], bins=40, range=(-0.7, 0.7),
                         weights=df_after.loc[mask, "weight"],
                         histtype='step', label=name, density=True)
        ax2.set_title("Global: Eta After Reweighting")
        ax2.set_xlabel("cluster_Eta"); ax2.legend(); ax2.grid(True)

        # Et Before Reweighting
        ax_et1 = fig.add_subplot(gs[1, 0])
        et_min, et_max = df_before[self.pt_col].min(), df_before[self.pt_col].max()
        for label_val, name in zip([0, 1], ["Background", "Signal"]):
            mask = df_before[self.label_col] == label_val
            if mask.any():
                ax_et1.hist(df_before.loc[mask, self.pt_col], bins=40, range=(et_min, et_max),
                            histtype='step', label=name, density=True)
        ax_et1.set_title("Global: Et Before Reweighting")
        ax_et1.set_xlabel(self.pt_col); ax_et1.set_ylabel("Density"); ax_et1.legend(); ax_et1.grid(True)

        # Et After Reweighting
        ax_et2 = fig.add_subplot(gs[1, 1])
        for label_val, name in zip([0, 1], ["Background", "Signal"]):
            mask = df_after[self.label_col] == label_val
            if mask.any():
                ax_et2.hist(df_after.loc[mask, self.pt_col], bins=40, range=(et_min, et_max),
                            weights=df_after.loc[mask, "weight"],
                            histtype='step', label=name, density=True)
        ax_et2.set_title("Global: Et After Reweighting")
        ax_et2.set_xlabel(self.pt_col); ax_et2.legend(); ax_et2.grid(True)

        # Vertex Z plots (only if vertex reweighting is enabled)
        print(f"DEBUG: About to check vertex condition...")
        print(f"DEBUG: vertex_reweight_enabled = {vertex_reweight_enabled}")
        print(f"DEBUG: 'vertexz' in df_before.columns = {'vertexz' in df_before.columns}")
        
        if vertex_reweight_enabled and 'vertexz' in df_before.columns:
            print("DEBUG: ✓ Creating vertex plots!")
            print(f"DEBUG: Grid shape will be 4x2")
            # Vertex Z Before Reweighting
            ax_vtx1 = fig.add_subplot(gs[2, 0])
            vtx_min, vtx_max = df_before['vertexz'].min(), df_before['vertexz'].max()
            vtx_range = (vtx_min, min(vtx_max, 100))  # Cap at 100 cm for visualization
            for label_val, name in zip([0, 1], ["Background", "Signal"]):
                mask = df_before[self.label_col] == label_val
                if mask.any():
                    ax_vtx1.hist(df_before.loc[mask, 'vertexz'], bins=50, range=vtx_range,
                                histtype='step', label=name, density=True)
            ax_vtx1.set_title("Global: Vertex Z Before Reweighting")
            ax_vtx1.set_xlabel("Vertex Z [cm]"); ax_vtx1.set_ylabel("Density"); ax_vtx1.legend(); ax_vtx1.grid(True)

            # Vertex Z After Reweighting
            ax_vtx2 = fig.add_subplot(gs[2, 1])
            for label_val, name in zip([0, 1], ["Background", "Signal"]):
                mask = df_after[self.label_col] == label_val
                if mask.any():
                    ax_vtx2.hist(df_after.loc[mask, 'vertexz'], bins=50, range=vtx_range,
                                weights=df_after.loc[mask, "weight"],
                                histtype='step', label=name, density=True)
            ax_vtx2.set_title("Global: Vertex Z After Reweighting")
            ax_vtx2.set_xlabel("Vertex Z [cm]"); ax_vtx2.legend(); ax_vtx2.grid(True)

        # Class Balance Before (counts)
        class_row = 3 if vertex_reweight_enabled else 2
        ax3 = fig.add_subplot(gs[class_row, 0])
        counts_before = df_before[self.label_col].value_counts()
        n_sig_before = counts_before.get(1, 0)
        n_bkg_before = counts_before.get(0, 0)
        ax3.bar(["Background", "Signal"], [n_bkg_before, n_sig_before], color=['blue', 'orange'])
        ax3.set_title("Global Class Counts Before Reweighting")
        ax3.set_ylabel("Number of Events"); ax3.grid(axis='y')
        for i, v in enumerate([n_bkg_before, n_sig_before]):
            ax3.text(i, v, str(v), ha='center', va='bottom')

        # Class Balance After (weighted yields)
        ax4 = fig.add_subplot(gs[class_row, 1])
        w_sig_after = df_after.loc[df_after[self.label_col] == 1, "weight"].sum() if "weight" in df_after.columns else float((df_after[self.label_col] == 1).sum())
        w_bkg_after = df_after.loc[df_after[self.label_col] == 0, "weight"].sum() if "weight" in df_after.columns else float((df_after[self.label_col] == 0).sum())
        ax4.bar(["Background", "Signal"], [w_bkg_after, w_sig_after], color=['blue', 'orange'])
        ax4.set_title("Global Class Yield After Reweighting")
        ax4.set_ylabel("Sum of Weights"); ax4.grid(axis='y')
        for i, v in enumerate([w_bkg_after, w_sig_after]):
            ax4.text(i, v, f"{v:.1f}", ha='center', va='bottom')

        plt.tight_layout()
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig)
        
        if self.save_plots:
            fig.savefig(f"{self.output_dir}/{self.root_file_prefix}_reweighting_qa_global.pdf", dpi=300, bbox_inches='tight')
        plt.show()

    def plot_isoet_bdt_correlation_global(self, save_plots: bool = True) -> None:
        """Global isoET vs BDT score correlation (all bins combined) on test set."""
        if not self.trained_pipelines:
            print("No trained pipelines available for correlation plots")
            return

        # Concatenate test slices across bins
        dfs = []
        for bin_label in self.bin_labels:
            if bin_label not in self.trained_pipelines:
                continue
            entry = self.trained_pipelines[bin_label]
            pipe = entry["pipeline"]
            X_test = entry["X_test"]
            df_bin = entry["df_bin"].loc[entry["test_index"]].copy()
            if hasattr(pipe.named_steps["clf"], "predict_proba"):
                y_proba = pipe.predict_proba(X_test)[:, 1]
            else:
                y_proba = pipe.decision_function(X_test)
            df_bin["bdt_score"] = y_proba
            dfs.append(df_bin[[self.iso_col, self.label_col, "bdt_score"]])

        if not dfs:
            print("No data for global correlation")
            return

        df_all = pd.concat(dfs, axis=0, ignore_index=True)

        # Scatter plot
        fig, ax = plt.subplots(figsize=(7, 6))
        # Sample for visualization
        max_points = 15000
        plot_df = df_all.sample(n=min(max_points, len(df_all)), random_state=42) if len(df_all) > max_points else df_all
        bkg = plot_df[plot_df[self.label_col] == 0]
        sig = plot_df[plot_df[self.label_col] == 1]
        if len(bkg) > 0:
            ax.scatter(bkg[self.iso_col], bkg["bdt_score"], alpha=0.25, s=1, color='red', label='Background')
        if len(sig) > 0:
            ax.scatter(sig[self.iso_col], sig["bdt_score"], alpha=0.35, s=1, color='blue', label='Signal')
        ax.set_xlabel(f"{self.iso_col} (GeV)"); ax.set_ylabel('BDT Score'); ax.grid(True, alpha=0.3)
        ax.set_title("Global: isoET vs BDT Score")
        ax.legend()

        # Correlations
        corr_overall = np.corrcoef(df_all[self.iso_col], df_all["bdt_score"])[0, 1] if len(df_all) > 10 else np.nan
        textstr = f"Overall: {corr_overall:.3f}"
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)

        plt.tight_layout()
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig)
        
        if self.save_plots and save_plots:
            fig.savefig(f"{self.output_dir}/{self.root_file_prefix}_isoet_bdt_correlation_global.pdf", dpi=300, bbox_inches='tight')
            print(f"Saved: {self.output_dir}/{self.root_file_prefix}_isoet_bdt_correlation_global.pdf")
        plt.show()
    
    def plot_feature_profile(self, df_bin: pd.DataFrame, val_index: pd.Index, 
                           y_proba: np.ndarray, feature_name: str, bin_label: str, 
                           w_val: np.ndarray) -> None:
        """Generate feature profile plots."""
        val_df = df_bin.loc[val_index].copy()
        val_df['bdt_score'] = y_proba
        val_df['weight'] = w_val
        
        sig_df = val_df[val_df[self.label_col] == 1]
        bkg_df = val_df[val_df[self.label_col] == 0]
        
        # Axis rule: feature on y-axis, except for isoET -> swap axes
        invert_axes = not (feature_name == self.iso_col)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), sharey=invert_axes)
        title_xy = f"BDT Score vs. {feature_name}" if invert_axes else f"{feature_name} vs. BDT Score"
        fig.suptitle(f"Profile of {title_xy} (Bin: {bin_label})", fontsize=16)
        
        score_bins = np.linspace(0, 1, 41)
        
        # Plot signal and background
        self._plot_single_profile(ax1, sig_df, feature_name, score_bins, "Signal", invert_axes)
        self._plot_single_profile(ax2, bkg_df, feature_name, score_bins, "Background", invert_axes)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig, bin_label)
        
        if self.save_plots:
            fig.savefig(f"{self.output_dir}/{self.root_file_prefix}_profile_{feature_name}_{bin_label}.pdf", dpi=300, bbox_inches='tight')
        plt.show()
    
    def _plot_single_profile(self, ax: plt.Axes, df: pd.DataFrame, feature_name: str, 
                           score_bins: np.ndarray, title: str, invert_axes: bool) -> None:
        """Helper to plot a single profile."""
        if invert_axes:
            x_vals = df[feature_name]
            y_vals = df['bdt_score']
            w_vals = df['weight']
            
            if len(x_vals) > 0:
                x_min, x_max = np.nanpercentile(x_vals, [1, 99])
                if not np.isfinite(x_min) or not np.isfinite(x_max) or x_min == x_max:
                    x_min, x_max = np.nanmin(x_vals), np.nanmax(x_vals)
            else:
                x_min, x_max = 0.0, 1.0
            
            feat_bins = np.linspace(x_min, x_max, 51)
            ax.hist2d(x_vals, y_vals, bins=[feat_bins, score_bins], cmin=1, 
                     weights=w_vals, cmap='viridis')
            
            # Profile (mean of y in x-bins)
            bin_centers = 0.5 * (feat_bins[:-1] + feat_bins[1:])
            bin_means, _, _ = stats.binned_statistic(x_vals, y_vals, statistic='mean', bins=feat_bins)
            ax.plot(bin_centers, bin_means, 'r-o', lw=2, label='Profile')
            ax.set_xlabel(feature_name)
            ax.set_ylabel("BDT Score")
        else:
            x_vals = df['bdt_score']
            y_vals = df[feature_name]
            w_vals = df['weight']
            
            ax.hist2d(x_vals, y_vals, bins=[score_bins, 50], cmin=1, weights=w_vals, cmap='viridis')
            
            # Profile (mean of y in BDT-score bins)
            bin_centers = (score_bins[:-1] + score_bins[1:]) / 2
            bin_means, _, _ = stats.binned_statistic(x_vals, y_vals, statistic='mean', bins=score_bins)
            ax.plot(bin_centers, bin_means, 'r-o', lw=2, label='Profile')
            ax.set_xlabel("BDT Score")
            ax.set_ylabel(feature_name)
        
        r_val = self._weighted_pearson(x_vals, y_vals, w_vals)
        ax.set_title(f"{title} (r={r_val:.3f})")
        ax.legend()
        ax.grid(True)
    
    def _weighted_pearson(self, x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
        """Calculate weighted Pearson correlation."""
        mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(w)
        x, y, w = x[mask], y[mask], w[mask]
        
        if x.size < 2:
            return np.nan
        
        w_sum = np.sum(w)
        mx = np.sum(w * x) / w_sum
        my = np.sum(w * y) / w_sum
        cov = np.sum(w * (x - mx) * (y - my)) / w_sum
        sx = np.sqrt(np.sum(w * (x - mx) ** 2) / w_sum)
        sy = np.sqrt(np.sum(w * (y - my) ** 2) / w_sum)
        
        if sx == 0 or sy == 0:
            return np.nan
        
        return cov / (sx * sy)
    
    def plot_bdt_score_distributions(self, density: bool = True) -> None:
        """Plot BDT score distributions for each bin on test set."""
        for bin_label in self.bin_labels:
            if bin_label not in self.trained_pipelines:
                continue

            entry = self.trained_pipelines[bin_label]
            pipe = entry["pipeline"]
            X_test = entry["X_test"]
            df_test = entry["df_bin"].loc[entry["test_index"]].copy()
            w_test = entry["w_test"]

            # Get BDT scores on test set
            if hasattr(pipe.named_steps["clf"], "predict_proba"):
                scores = pipe.predict_proba(X_test)[:, 1]
                score_label = "BDT Score"
                bin_edges = np.linspace(0, 1, self.plot_config['n_bins'] + 1)
            else:
                scores = pipe.decision_function(X_test)
                score_label = "BDT Score (decision_function)"
                p1, p99 = np.nanpercentile(scores[np.isfinite(scores)], [1, 99])
                if not np.isfinite(p1) or not np.isfinite(p99) or p1 == p99:
                    p1, p99 = np.nanmin(scores), np.nanmax(scores)
                    if p1 == p99:
                        p1, p99 = p1 - 0.5, p1 + 0.5
                bin_edges = np.linspace(p1, p99, self.plot_config['n_bins'] + 1)

            df_test["score"] = scores
            df_test["w"] = w_test

            sig = df_test[df_test[self.label_col] == 1]
            bkg = df_test[df_test[self.label_col] == 0]
            
            # Build the plot
            fig = plt.figure(figsize=(7, 5))
            plt.hist(bkg["score"].values, bins=bin_edges, weights=bkg["w"].values,
                     histtype="stepfilled", alpha=0.35, label="Background", density=density)
            plt.hist(sig["score"].values, bins=bin_edges, weights=sig["w"].values,
                     histtype="step", lw=2, label="Signal", density=density)
            
            plt.xlabel(score_label)
            plt.ylabel("Density" if density else "Weighted Counts")
            plt.title(f"BDT Score Distribution - Bin {bin_label}")
            plt.legend()
            plt.grid(True)
            
            # Add sPHENIX header to the figure
            self._add_sphenix_header(fig, bin_label)
            
            if self.save_plots:
                fig.savefig(f"{self.output_dir}/{self.root_file_prefix}_bdt_score_{bin_label}.pdf", dpi=300, bbox_inches='tight')
            plt.show()

    def plot_auc_bar(self, val_results: Dict) -> None:
        """Plot test AUC by bin as a bar chart."""
        auc_data = [(b, val_results[b]["auc_val"]) for b in self.bin_labels if b in val_results]
        if not auc_data:
            return
        bins_, aucs_ = zip(*auc_data)
        fig = plt.figure()
        sns.barplot(x=list(bins_), y=list(aucs_), color="steelblue")
        plt.ylabel("Test AUC")
        plt.xlabel("cluster_Et bin [GeV]")
        plt.title("Test AUC by pT bin")
        plt.ylim(0.5, 1.0)
        for i, v in enumerate(aucs_):
            plt.text(i, v + 0.01, f"{v:.3f}", ha="center")
        plt.tight_layout()
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig)
        
        if self.save_plots:
            fig.savefig(f"{self.output_dir}/{self.root_file_prefix}_auc_bar.pdf", dpi=300, bbox_inches='tight')
        plt.show()

    def plot_combined_roc(self, val_results: Dict) -> None:
        """Plot combined ROC curves across bins."""
        fig = plt.figure(figsize=(8, 8))
        for bin_label in self.bin_labels:
            if bin_label not in val_results:
                continue
            r = val_results[bin_label]
            if r.get("fpr") is not None and r.get("tpr") is not None:
                plt.plot(r["fpr"], r["tpr"], label=f"{bin_label} (AUC={r['auc_val']:.3f})")
        plt.plot([0, 1], [0, 1], "k--", alpha=0.5, label="Chance")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title("ROC Curves for all pT Bins")
        plt.legend()
        plt.grid(True)
        plt.axis('square')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig)
        
        if self.save_plots:
            fig.savefig(f"{self.output_dir}/{self.root_file_prefix}_combined_roc.pdf", dpi=300, bbox_inches='tight')
        plt.show()

    def plot_tpr_at_fixed_fpr(self, val_results: Dict, target_fpr: float = 0.2) -> None:
        """Bar plot of TPR at a fixed FPR (default 0.2) per bin."""
        eff_points = []
        for bin_label in self.bin_labels:
            if bin_label not in val_results:
                continue
            r = val_results[bin_label]
            fpr = r.get("fpr")
            tpr = r.get("tpr")
            if fpr is None or tpr is None:
                continue
            fpr = np.asarray(fpr)
            tpr = np.asarray(tpr)
            if fpr.size == 0:
                continue
            order = np.argsort(fpr)
            fpr_sorted = fpr[order]
            tpr_sorted = tpr[order]
            if target_fpr <= fpr_sorted[0]:
                tpr_at = tpr_sorted[0]
            elif target_fpr >= fpr_sorted[-1]:
                tpr_at = tpr_sorted[-1]
            else:
                idx = np.searchsorted(fpr_sorted, target_fpr)
                f1, f2 = fpr_sorted[idx-1], fpr_sorted[idx]
                t1, t2 = tpr_sorted[idx-1], tpr_sorted[idx]
                w = (target_fpr - f1)/(f2 - f1) if f2 > f1 else 0.0
                tpr_at = t1 + w*(t2 - t1)
            eff_points.append((bin_label, float(tpr_at)))
        if not eff_points:
            return
        bins_eff, effs = zip(*eff_points)
        fig = plt.figure()
        sns.barplot(x=list(bins_eff), y=list(effs), color="darkorange")
        plt.ylabel(f"TPR @ FPR={target_fpr}")
        plt.xlabel("cluster_Et bin [GeV]")
        plt.title("Efficiency at 20% FPR by pT bin")
        plt.ylim(0.0, 1.0)
        for i, v in enumerate(effs):
            plt.text(i, min(v + 0.025, 0.97), f"{v:.3f}", ha="center", fontsize=9)
        plt.axhline(0.5, ls='--', color='gray', alpha=0.4)
        plt.tight_layout()
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig)
        
        if self.save_plots:
            fpr_tag = str(target_fpr).replace('.', 'p')
            fig.savefig(f"{self.output_dir}/{self.root_file_prefix}_tpr_at_fpr_{fpr_tag}.pdf", dpi=300, bbox_inches='tight')
        plt.show()

    def overlay_feature_profiles(self, feature_name: str) -> None:
        """Overlay profile curves across bins for a feature, matching original notebook behavior."""
        xs_sig_all, xs_bkg_all = [], []
        profiles_sig, profiles_bkg = {}, {}
        invert_axes = not (feature_name == self.iso_col)
        score_bins = np.linspace(0, 1, 41)

        # Gather values to define common x-bins
        for bin_label in self.bin_labels:
            if bin_label not in self.trained_pipelines:
                continue
            entry = self.trained_pipelines[bin_label]
            pipe = entry["pipeline"]
            X_val = entry["X_val"]
            df_bin = entry["df_bin"].loc[entry["val_index"]].copy()
            if hasattr(pipe.named_steps["clf"], "predict_proba"):
                y_proba = pipe.predict_proba(X_val)[:, 1]
            else:
                y_proba = pipe.decision_function(X_val)
            df_bin["bdt_score"] = y_proba
            sig = df_bin[df_bin[self.label_col] == 1]
            bkg = df_bin[df_bin[self.label_col] == 0]
            if invert_axes:
                xs_sig_all.append(sig[feature_name].to_numpy())
                xs_bkg_all.append(bkg[feature_name].to_numpy())
            else:
                xs_sig_all.append(sig["bdt_score"].to_numpy())
                xs_bkg_all.append(bkg["bdt_score"].to_numpy())

        def robust_bins(arr_list, nbins=50):
            if len(arr_list) == 0:
                return np.linspace(0, 1, nbins+1)
            x_all = np.concatenate([a[np.isfinite(a)] for a in arr_list if a is not None])
            if x_all.size == 0:
                return np.linspace(0, 1, nbins+1)
            x1, x99 = np.nanpercentile(x_all, [1, 99])
            if not np.isfinite(x1) or not np.isfinite(x99) or x1 == x99:
                x1, x99 = np.nanmin(x_all), np.nanmax(x_all)
                if x1 == x99:
                    x1, x99 = x1 - 0.5, x1 + 0.5
            return np.linspace(x1, x99, nbins+1)

        if invert_axes:
            xbins_sig = robust_bins(xs_sig_all, nbins=50)
            xbins_bkg = robust_bins(xs_bkg_all, nbins=50)
        else:
            xbins_sig = score_bins
            xbins_bkg = score_bins

        def binned_mean(x, y, bins):
            inds = np.digitize(x, bins) - 1
            centers = 0.5 * (bins[:-1] + bins[1:])
            means = np.full(len(centers), np.nan)
            for bi in range(len(centers)):
                m = inds == bi
                if np.any(m):
                    means[bi] = np.nanmean(y[m])
            return centers, means

        # Compute profiles per bin
        for bin_label in self.bin_labels:
            if bin_label not in self.trained_pipelines:
                continue
            entry = self.trained_pipelines[bin_label]
            pipe = entry["pipeline"]
            X_val = entry["X_val"]
            df_bin = entry["df_bin"].loc[entry["val_index"]].copy()
            if hasattr(pipe.named_steps["clf"], "predict_proba"):
                y_proba = pipe.predict_proba(X_val)[:, 1]
            else:
                y_proba = pipe.decision_function(X_val)
            df_bin["bdt_score"] = y_proba
            sig = df_bin[df_bin[self.label_col] == 1]
            bkg = df_bin[df_bin[self.label_col] == 0]
            if invert_axes:
                x_s, y_s = sig[feature_name].to_numpy(), sig["bdt_score"].to_numpy()
                x_b, y_b = bkg[feature_name].to_numpy(), bkg["bdt_score"].to_numpy()
            else:
                x_s, y_s = sig["bdt_score"].to_numpy(), sig[feature_name].to_numpy()
                x_b, y_b = bkg["bdt_score"].to_numpy(), bkg[feature_name].to_numpy()
            cs, ms = binned_mean(x_s, y_s, xbins_sig)
            cb, mb = binned_mean(x_b, y_b, xbins_bkg)
            profiles_sig[bin_label] = (cs, ms)
            profiles_bkg[bin_label] = (cb, mb)

        # Plot overlays
        fig, (axS, axB) = plt.subplots(1, 2, figsize=(16, 6), sharey=not invert_axes)
        if invert_axes:
            axS.set_xlabel(feature_name); axS.set_ylabel("BDT Score")
            axB.set_xlabel(feature_name); axB.set_ylabel("BDT Score")
            supt = f"Overlayed Profiles: BDT vs {feature_name}"
        else:
            axS.set_xlabel("BDT Score");   axS.set_ylabel(feature_name)
            axB.set_xlabel("BDT Score");   axB.set_ylabel(feature_name)
            supt = f"Overlayed Profiles: {feature_name} vs BDT"
        for bin_label, (c, m) in profiles_sig.items():
            axS.plot(c, m, marker='o', lw=2, label=str(bin_label))
        axS.set_title("Signal"); axS.grid(True); axS.legend()
        for bin_label, (c, m) in profiles_bkg.items():
            axB.plot(c, m, marker='o', lw=2, label=str(bin_label))
        axB.set_title("Background"); axB.grid(True); axB.legend()
        fig.suptitle(supt, fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig)
        
        if self.save_plots:
            fig.savefig(f"{self.output_dir}/{self.root_file_prefix}_overlay_profile_{feature_name}.pdf", dpi=300, bbox_inches='tight')
        plt.show()
    
    def plot_feature_profile_global(self, feature_name: str) -> None:
        """Global 2D histogram/profile: showershape vs BDT score across all bins.
        Uses validation slices from each bin and weights.
        """
        # Gather validation data across bins
        dfs = []
        for bin_label in self.bin_labels:
            if bin_label not in self.trained_pipelines:
                continue
            entry = self.trained_pipelines[bin_label]
            pipe = entry["pipeline"]
            X_val = entry["X_val"]
            df_bin = entry["df_bin"].loc[entry["val_index"]].copy()
            w_val = entry["w_val"]
            if feature_name not in df_bin.columns:
                continue
            # Scores
            if hasattr(pipe.named_steps["clf"], "predict_proba"):
                y_proba = pipe.predict_proba(X_val)[:, 1]
            else:
                y_proba = pipe.decision_function(X_val)
            df_bin["bdt_score"] = y_proba
            df_bin["weight"] = w_val
            dfs.append(df_bin[[feature_name, self.label_col, "bdt_score", "weight"]])

        if not dfs:
            print(f"No data available for global feature profile: {feature_name}")
            return

        df_all = pd.concat(dfs, axis=0, ignore_index=True)

        # Axis rule consistent with per-bin profiles
        invert_axes = not (feature_name == self.iso_col)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), sharey=invert_axes)
        title_xy = f"BDT Score vs. {feature_name}" if invert_axes else f"{feature_name} vs. BDT Score"
        fig.suptitle(f"Global Profile of {title_xy}", fontsize=16)

        score_bins = np.linspace(0, 1, 41)

        def plot_one(ax, df, title):
            if invert_axes:
                x_vals = df[feature_name].to_numpy()
                y_vals = df["bdt_score"].to_numpy()
                w_vals = df["weight"].to_numpy()
                # robust x-bins
                if x_vals.size > 0:
                    x_min, x_max = np.nanpercentile(x_vals, [1, 99])
                    if not np.isfinite(x_min) or not np.isfinite(x_max) or x_min == x_max:
                        x_min, x_max = np.nanmin(x_vals), np.nanmax(x_vals)
                else:
                    x_min, x_max = 0.0, 1.0
                feat_bins = np.linspace(x_min, x_max, 51)
                h = ax.hist2d(x_vals, y_vals, bins=[feat_bins, score_bins], cmin=1, weights=w_vals, cmap='viridis')
                # profile of mean
                bin_centers = 0.5 * (feat_bins[:-1] + feat_bins[1:])
                means, _, _ = stats.binned_statistic(x_vals, y_vals, statistic='mean', bins=feat_bins)
                ax.plot(bin_centers, means, 'r-o', lw=2, label='Profile')
                ax.set_xlabel(feature_name)
                ax.set_ylabel("BDT Score")
            else:
                x_vals = df["bdt_score"].to_numpy()
                y_vals = df[feature_name].to_numpy()
                w_vals = df["weight"].to_numpy()
                h = ax.hist2d(x_vals, y_vals, bins=[score_bins, 50], cmin=1, weights=w_vals, cmap='viridis')
                # profile of mean
                bin_centers = (score_bins[:-1] + score_bins[1:]) / 2
                means, _, _ = stats.binned_statistic(x_vals, y_vals, statistic='mean', bins=score_bins)
                ax.plot(bin_centers, means, 'r-o', lw=2, label='Profile')
                ax.set_xlabel("BDT Score")
                ax.set_ylabel(feature_name)
            # weighted Pearson
            r_val = self._weighted_pearson(x_vals, y_vals, w_vals)
            ax.set_title(f"{title} (r={r_val:.3f})")
            ax.legend(); ax.grid(True)
            return h

        sig_df = df_all[df_all[self.label_col] == 1]
        bkg_df = df_all[df_all[self.label_col] == 0]
        h1 = plot_one(ax1, sig_df, "Signal")
        h2 = plot_one(ax2, bkg_df, "Background")

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig)
        
        if self.save_plots:
            fig.savefig(f"{self.output_dir}/{self.root_file_prefix}_profile_global_{feature_name}.pdf", dpi=300, bbox_inches='tight')
        plt.show()
    
    def plot_isoet_bdt_correlation(self, save_plots: bool = True) -> None:
        """Plot correlation between isoET and BDT score with correlation values."""
        print("Generating isoET vs BDT score correlation plots...")
        
        # Check if we have the required data
        if not self.trained_pipelines:
            print("No trained pipelines available for correlation plots")
            return
        
        # Collect correlation data for all bins
        correlation_data = {}
        
        # Create individual plots for each bin
        n_bins = len(self.bin_labels)
        fig, axes = plt.subplots(1, n_bins, figsize=(6*n_bins, 6))
        
        if n_bins == 1:
            axes = [axes]
        
        for idx, bin_label in enumerate(self.bin_labels):
            if bin_label not in self.trained_pipelines:
                continue
                
            entry = self.trained_pipelines[bin_label]
            pipe = entry["pipeline"]
            X_val = entry["X_val"]
            df_bin = entry["df_bin"].loc[entry["val_index"]].copy()
            
            # Get BDT scores
            if hasattr(pipe.named_steps["clf"], "predict_proba"):
                y_proba = pipe.predict_proba(X_val)[:, 1]
            else:
                y_proba = pipe.decision_function(X_val)
            
            df_bin["bdt_score"] = y_proba
            
            # Check if isoET column exists
            if self.iso_col not in df_bin.columns:
                print(f"Warning: {self.iso_col} not found in data for bin {bin_label}")
                continue
            
            # Separate signal and background
            sig_data = df_bin[df_bin[self.label_col] == 1]
            bkg_data = df_bin[df_bin[self.label_col] == 0]
            
            # Calculate correlations
            correlations = {}
            if len(sig_data) > 10:  # Need sufficient points for correlation
                sig_corr = np.corrcoef(sig_data[self.iso_col], sig_data["bdt_score"])[0, 1]
                correlations['signal'] = sig_corr
            else:
                correlations['signal'] = np.nan
                
            if len(bkg_data) > 10:
                bkg_corr = np.corrcoef(bkg_data[self.iso_col], bkg_data["bdt_score"])[0, 1]
                correlations['background'] = bkg_corr
            else:
                correlations['background'] = np.nan
            
            # Overall correlation
            if len(df_bin) > 10:
                overall_corr = np.corrcoef(df_bin[self.iso_col], df_bin["bdt_score"])[0, 1]
                correlations['overall'] = overall_corr
            else:
                correlations['overall'] = np.nan
            
            correlation_data[bin_label] = correlations
            
            # Create scatter plot
            ax = axes[idx]
            
            # Sample data for visualization if too many points
            max_points = 5000
            if len(sig_data) > max_points:
                sig_sample = sig_data.sample(n=max_points, random_state=42)
            else:
                sig_sample = sig_data
                
            if len(bkg_data) > max_points:
                bkg_sample = bkg_data.sample(n=max_points, random_state=42)
            else:
                bkg_sample = bkg_data
            
            # Plot points
            if len(bkg_sample) > 0:
                ax.scatter(bkg_sample[self.iso_col], bkg_sample["bdt_score"], 
                          alpha=0.3, s=1, color='red', label='Background')
            
            if len(sig_sample) > 0:
                ax.scatter(sig_sample[self.iso_col], sig_sample["bdt_score"], 
                          alpha=0.5, s=1, color='blue', label='Signal')
            
            # Add correlation text
            textstr = f'Correlations:\\n'
            if not np.isnan(correlations['overall']):
                textstr += f'Overall: {correlations["overall"]:.3f}\\n'
            if not np.isnan(correlations['signal']):
                textstr += f'Signal: {correlations["signal"]:.3f}\\n'
            if not np.isnan(correlations['background']):
                textstr += f'Background: {correlations["background"]:.3f}'
            
            # Position text box
            props = dict(boxstyle='round', facecolor='white', alpha=0.8)
            ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                   verticalalignment='top', bbox=props)
            
            ax.set_xlabel(f'{self.iso_col} (GeV)')
            ax.set_ylabel('BDT Score')
            ax.set_title(f'Bin {bin_label}: {self.iso_col} vs BDT Score')
            ax.grid(True, alpha=0.3)
            ax.legend()
        
        plt.suptitle(f'{self.iso_col} vs BDT Score Correlation', fontsize=16)
        plt.tight_layout()
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig)
        
        if save_plots:
            plt.savefig(f"{self.config['output']['output_dir']}/{self.root_file_prefix}_isoet_bdt_correlation.pdf", 
                       dpi=300, bbox_inches='tight')
            print(f"Saved: {self.config['output']['output_dir']}/{self.root_file_prefix}_isoet_bdt_correlation.pdf")
        
        plt.show()
        
        # Create summary correlation plot
        self._plot_correlation_summary(correlation_data, save_plots)
        
        # Print correlation summary
        print("\\n=== isoET vs BDT Score Correlation Summary ===")
        for bin_label, corrs in correlation_data.items():
            print(f"Bin {bin_label}:")
            for corr_type, value in corrs.items():
                if not np.isnan(value):
                    print(f"  {corr_type.capitalize()}: {value:.4f}")
                else:
                    print(f"  {corr_type.capitalize()}: N/A")
        
        return correlation_data
    
    def _plot_correlation_summary(self, correlation_data: Dict, save_plots: bool) -> None:
        """Plot summary of correlations across bins."""
        # Prepare data for plotting
        bins = list(correlation_data.keys())
        overall_corrs = [correlation_data[bin_label]['overall'] for bin_label in bins]
        signal_corrs = [correlation_data[bin_label]['signal'] for bin_label in bins]
        background_corrs = [correlation_data[bin_label]['background'] for bin_label in bins]
        
        # Create bar plot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        x = np.arange(len(bins))
        width = 0.25
        
        # Remove NaN values for plotting
        overall_clean = [v if not np.isnan(v) else 0 for v in overall_corrs]
        signal_clean = [v if not np.isnan(v) else 0 for v in signal_corrs]
        background_clean = [v if not np.isnan(v) else 0 for v in background_corrs]
        
        bars1 = ax.bar(x - width, overall_clean, width, label='Overall', alpha=0.8)
        bars2 = ax.bar(x, signal_clean, width, label='Signal', alpha=0.8)
        bars3 = ax.bar(x + width, background_clean, width, label='Background', alpha=0.8)
        
        # Add value labels on bars
        def add_value_labels(bars, values):
            for bar, val in zip(bars, values):
                if val != 0:  # Don't label zero values (NaN replacements)
                    height = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width()/2., height + 0.01 * np.sign(height),
                           f'{val:.3f}', ha='center', va='bottom' if height >= 0 else 'top')
        
        add_value_labels(bars1, [v if not np.isnan(v) else 0 for v in overall_corrs])
        add_value_labels(bars2, [v if not np.isnan(v) else 0 for v in signal_corrs])
        add_value_labels(bars3, [v if not np.isnan(v) else 0 for v in background_corrs])
        
        ax.set_xlabel('ET Bins')
        ax.set_ylabel('Pearson Correlation Coefficient')
        ax.set_title(f'{self.iso_col} vs BDT Score - Correlation Summary')
        ax.set_xticks(x)
        ax.set_xticklabels(bins)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        
        plt.tight_layout()
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig)
        
        if save_plots:
            plt.savefig(f"{self.config['output']['output_dir']}/{self.root_file_prefix}_isoet_bdt_correlation_summary.pdf", 
                       dpi=300, bbox_inches='tight')
            print(f"Saved: {self.config['output']['output_dir']}/{self.root_file_prefix}_isoet_bdt_correlation_summary.pdf")
        
        plt.show()
    
    def plot_feature_importance(self, save_plots: bool = True) -> None:
        """Plot feature importance for XGBoost models using gain, weight, and cover metrics."""
        if not XGB_AVAILABLE:
            print("XGBoost not available - skipping feature importance plots")
            return
        
        print("Generating feature importance plots...")
        
        # Check if we have XGBoost models
        has_xgb_models = False
        for bin_label, entry in self.trained_pipelines.items():
            pipeline = entry.get("pipeline")
            if pipeline and hasattr(pipeline.named_steps.get("clf"), "get_booster"):
                has_xgb_models = True
                break
        
        if not has_xgb_models:
            print("No XGBoost models found - skipping feature importance plots")
            return
        
        # Extract feature importance for each bin
        importance_data = {}
        feature_names = None
        
        for bin_label, entry in self.trained_pipelines.items():
            pipeline = entry.get("pipeline")
            features = entry.get("features", [])
            
            if not pipeline or not hasattr(pipeline.named_steps.get("clf"), "get_booster"):
                continue
                
            classifier = pipeline.named_steps["clf"]
            booster = classifier.get_booster()
            
            # Get feature importance scores
            try:
                gain_scores = booster.get_score(importance_type='gain')
                weight_scores = booster.get_score(importance_type='weight')
                cover_scores = booster.get_score(importance_type='cover')
                
                if feature_names is None:
                    feature_names = features
                
                # Map feature indices back to feature names
                gain_dict = {}
                weight_dict = {}
                cover_dict = {}
                
                for i, feature_name in enumerate(features):
                    f_key = f"f{i}"  # XGBoost uses f0, f1, f2, etc.
                    gain_dict[feature_name] = gain_scores.get(f_key, 0.0)
                    weight_dict[feature_name] = weight_scores.get(f_key, 0.0)
                    cover_dict[feature_name] = cover_scores.get(f_key, 0.0)
                
                importance_data[bin_label] = {
                    'gain': gain_dict,
                    'weight': weight_dict,
                    'cover': cover_dict
                }
                
            except Exception as e:
                print(f"Warning: Could not extract feature importance for bin {bin_label}: {e}")
                continue
        
        if not importance_data:
            print("No feature importance data extracted")
            return
        
        # Create comprehensive feature importance plots
        self._plot_feature_importance_comparison(importance_data, feature_names, save_plots)
        self._plot_feature_importance_heatmap(importance_data, feature_names, save_plots)
        self._plot_top_features_per_bin(importance_data, feature_names, save_plots)
    
    def _plot_feature_importance_comparison(self, importance_data: Dict, feature_names: List[str], save_plots: bool) -> None:
        """Plot feature importance comparison across bins for each importance type."""
        importance_types = ['gain', 'weight', 'cover']
        
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))
        
        for idx, imp_type in enumerate(importance_types):
            ax = axes[idx]
            
            # Prepare data for plotting
            plot_data = []
            for bin_label, data in importance_data.items():
                for feature, score in data[imp_type].items():
                    plot_data.append({
                        'Feature': feature,
                        'Importance': score,
                        'Bin': bin_label,
                        'Type': imp_type
                    })
            
            if not plot_data:
                continue
                
            df_plot = pd.DataFrame(plot_data)
            
            # Create grouped bar plot
            pivot_data = df_plot.pivot(index='Feature', columns='Bin', values='Importance').fillna(0)
            
            # Sort by mean importance across bins
            feature_order = pivot_data.mean(axis=1).sort_values(ascending=False).index
            pivot_data = pivot_data.loc[feature_order]
            
            # Plot with better spacing
            pivot_data.plot(kind='bar', ax=ax, width=0.8)
            ax.set_title(f'Feature Importance - {imp_type.title()}')
            ax.set_xlabel('Features')
            ax.set_ylabel(f'{imp_type.title()} Score')
            ax.legend(title='ET Bin', bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.grid(True, alpha=0.3)
            
            # Rotate x-axis labels for better readability
            ax.tick_params(axis='x', rotation=45)
        
        plt.suptitle('Feature Importance Comparison Across ET Bins', fontsize=16, y=1.02)
        plt.tight_layout()
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig)
        
        if save_plots:
            plt.savefig(f"{self.config['output']['output_dir']}/{self.root_file_prefix}_feature_importance_comparison.pdf", 
                       dpi=300, bbox_inches='tight')
            print(f"Saved: {self.config['output']['output_dir']}/{self.root_file_prefix}_feature_importance_comparison.pdf")
        
        plt.show()
    
    def _plot_feature_importance_heatmap(self, importance_data: Dict, feature_names: List[str], save_plots: bool) -> None:
        """Plot feature importance heatmap showing gain importance across bins."""
        # Focus on gain importance for heatmap
        heatmap_data = []
        
        for bin_label, data in importance_data.items():
            row_data = []
            for feature in feature_names:
                row_data.append(data['gain'].get(feature, 0.0))
            heatmap_data.append(row_data)
        
        # Create DataFrame for heatmap
        df_heatmap = pd.DataFrame(heatmap_data, 
                                 index=list(importance_data.keys()), 
                                 columns=feature_names)
        
        # Normalize by row (each bin) to show relative importance
        df_normalized = df_heatmap.div(df_heatmap.sum(axis=1), axis=0).fillna(0)
        
        # Create heatmap
        fig, axes = plt.subplots(1, 2, figsize=(20, 6))
        
        # Raw values
        sns.heatmap(df_heatmap, annot=True, fmt='.0f', cmap='YlOrRd', ax=axes[0])
        axes[0].set_title('Feature Importance (Gain) - Raw Values')
        axes[0].set_xlabel('Features')
        axes[0].set_ylabel('ET Bins')
        
        # Normalized values (percentages)
        sns.heatmap(df_normalized, annot=True, fmt='.2f', cmap='YlOrRd', ax=axes[1])
        axes[1].set_title('Feature Importance (Gain) - Relative per Bin')
        axes[1].set_xlabel('Features')
        axes[1].set_ylabel('ET Bins')
        
        plt.suptitle('Feature Importance Heatmaps', fontsize=16)
        plt.tight_layout()
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig)
        
        if save_plots:
            plt.savefig(f"{self.config['output']['output_dir']}/{self.root_file_prefix}_feature_importance_heatmap.pdf", 
                       dpi=300, bbox_inches='tight')
            print(f"Saved: {self.config['output']['output_dir']}/{self.root_file_prefix}_feature_importance_heatmap.pdf")
        
        plt.show()
    
    def _plot_top_features_per_bin(self, importance_data: Dict, feature_names: List[str], save_plots: bool, top_n: int = 15) -> None:
        """Plot top N most important features for each bin."""
        n_bins = len(importance_data)
        fig, axes = plt.subplots(1, n_bins, figsize=(6*n_bins, 6))
        
        if n_bins == 1:
            axes = [axes]
        
        for idx, (bin_label, data) in enumerate(importance_data.items()):
            ax = axes[idx]
            
            # Get top features by gain
            gain_scores = data['gain']
            sorted_features = sorted(gain_scores.items(), key=lambda x: x[1], reverse=True)[:top_n]
            
            if not sorted_features:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'Bin {bin_label}')
                continue
            
            features, scores = zip(*sorted_features)
            
            # Create horizontal bar plot
            y_pos = np.arange(len(features))
            bars = ax.barh(y_pos, scores, color=plt.cm.viridis(np.linspace(0, 1, len(features))))
            
            ax.set_yticks(y_pos)
            ax.set_yticklabels(features)
            ax.set_xlabel('Gain Score')
            ax.set_title(f'Top {len(features)} Features - Bin {bin_label}')
            ax.grid(True, alpha=0.3, axis='x')
            
            # Add value labels on bars
            for i, (bar, score) in enumerate(zip(bars, scores)):
                ax.text(score + max(scores)*0.01, i, f'{score:.0f}', 
                       va='center', ha='left', fontsize=9)
            
            # Invert y-axis to show most important at top
            ax.invert_yaxis()
        
        plt.suptitle(f'Top {top_n} Most Important Features by Gain (Per ET Bin)', fontsize=16)
        plt.tight_layout()
        
        # Add sPHENIX header to the figure
        self._add_sphenix_header(fig)
        
        if save_plots:
            plt.savefig(f"{self.config['output']['output_dir']}/{self.root_file_prefix}_top_features_per_bin.pdf", 
                       dpi=300, bbox_inches='tight')
            print(f"Saved: {self.config['output']['output_dir']}/{self.root_file_prefix}_top_features_per_bin.pdf")
        
        plt.show()
