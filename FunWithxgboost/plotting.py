"""
Plotting utilities for binned training.
"""
import pandas as pd
from typing import Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats


class BinnedTrainingPlotter:
    """Handles all plotting functionality for binned training."""
    
    def __init__(self, trained_pipelines: Dict, config: Dict):
        self.trained_pipelines = trained_pipelines
        self.config = config
        self.plot_config = config['plotting']
        self.bin_labels = config['binning']['labels']
        self.label_col = config['data']['label_column']
        self.pt_col = config['data']['pt_column']
        self.iso_col = config['data']['iso_column']
        
        # Set matplotlib defaults
        plt.rcParams["figure.figsize"] = tuple(self.plot_config['figure_size'])
        plt.rcParams["axes.grid"] = self.plot_config['grid']
        sns.set_context(self.plot_config['context'])
    
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
        """Plot BDT score distributions for each bin."""
        for bin_label in self.bin_labels:
            if bin_label not in self.trained_pipelines:
                continue
            
            entry = self.trained_pipelines[bin_label]
            pipe = entry["pipeline"]
            X_val = entry["X_val"]
            df_val = entry["df_bin"].loc[entry["val_index"]].copy()
            w_val = entry["w_val"]
            
            # Get BDT scores
            if hasattr(pipe.named_steps["clf"], "predict_proba"):
                scores = pipe.predict_proba(X_val)[:, 1]
                score_label = "BDT Score"
                bin_edges = np.linspace(0, 1, self.plot_config['n_bins'] + 1)
            else:
                scores = pipe.decision_function(X_val)
                score_label = "BDT Score (decision_function)"
                p1, p99 = np.nanpercentile(scores[np.isfinite(scores)], [1, 99])
                if not np.isfinite(p1) or not np.isfinite(p99) or p1 == p99:
                    p1, p99 = np.nanmin(scores), np.nanmax(scores)
                    if p1 == p99:
                        p1, p99 = p1 - 0.5, p1 + 0.5
                bin_edges = np.linspace(p1, p99, self.plot_config['n_bins'] + 1)
            
            df_val["score"] = scores
            df_val["w"] = w_val
            
            sig = df_val[df_val[self.label_col] == 1]
            bkg = df_val[df_val[self.label_col] == 0]
            
            # Build the plot
            plt.figure(figsize=(7, 5))
            plt.hist(bkg["score"].values, bins=bin_edges, weights=bkg["w"].values,
                     histtype="stepfilled", alpha=0.35, label="Background", density=density)
            plt.hist(sig["score"].values, bins=bin_edges, weights=sig["w"].values,
                     histtype="step", lw=2, label="Signal", density=density)
            
            plt.xlabel(score_label)
            plt.ylabel("Density" if density else "Weighted Counts")
            plt.title(f"BDT Score Distribution - Bin {bin_label}")
            plt.legend()
            plt.grid(True)
            plt.show()

    def plot_auc_bar(self, val_results: Dict) -> None:
        """Plot AUC by bin as a bar chart."""
        auc_data = [(b, val_results[b]["auc_val"]) for b in self.bin_labels if b in val_results]
        if not auc_data:
            return
        bins_, aucs_ = zip(*auc_data)
        plt.figure()
        sns.barplot(x=list(bins_), y=list(aucs_), color="steelblue")
        plt.ylabel("Validation AUC")
        plt.xlabel("cluster_Et bin [GeV]")
        plt.title("AUC by pT bin")
        plt.ylim(0.5, 1.0)
        for i, v in enumerate(aucs_):
            plt.text(i, v + 0.01, f"{v:.3f}", ha="center")
        plt.tight_layout()
        plt.show()

    def plot_combined_roc(self, val_results: Dict) -> None:
        """Plot combined ROC curves across bins."""
        plt.figure(figsize=(8, 8))
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
        plt.figure()
        sns.barplot(x=list(bins_eff), y=list(effs), color="darkorange")
        plt.ylabel(f"TPR @ FPR={target_fpr}")
        plt.xlabel("cluster_Et bin [GeV]")
        plt.title("Efficiency at 20% FPR by pT bin")
        plt.ylim(0.0, 1.0)
        for i, v in enumerate(effs):
            plt.text(i, min(v + 0.025, 0.97), f"{v:.3f}", ha="center", fontsize=9)
        plt.axhline(0.5, ls='--', color='gray', alpha=0.4)
        plt.tight_layout()
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
        plt.show()
