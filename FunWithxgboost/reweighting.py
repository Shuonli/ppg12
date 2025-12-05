"""
Kinematic reweighting utilities for binned training.
"""
from typing import Dict
import pandas as pd
import numpy as np
from scipy.interpolate import UnivariateSpline


class KinematicReweighter:
    """Handles kinematic reweighting for signal and background samples."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.reweight_config = config['reweighting']
        self.pt_col = config['data']['pt_column']
        self.label_col = config['data']['label_column']
    
    def apply_class_reweighting(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply class imbalance reweighting."""
        df_out = df.copy()
        
        if self.reweight_config.get('class_reweight', False):
            n_sig = (df_out[self.label_col] == 1).sum()
            n_bkg = (df_out[self.label_col] == 0).sum()
            
            if n_sig > 0 and n_bkg > 0:
                weight_sig = (n_sig + n_bkg) / (2.0 * n_sig)
                weight_bkg = (n_sig + n_bkg) / (2.0 * n_bkg)
                df_out['class_weight'] = np.where(df_out[self.label_col] == 1, weight_sig, weight_bkg)
            else:
                df_out['class_weight'] = 1.0
        else:
            df_out['class_weight'] = 1.0
        
        df_out["weight"] = df_out['class_weight']
        return df_out
    
    def apply_eta_reweighting(self, df: pd.DataFrame, label: int) -> pd.DataFrame:
        """Apply eta reweighting for a specific class."""
        df_out = df.copy()
        mask = df_out[self.label_col] == label
        
        if not mask.any():
            return df_out
        
        # Eta Reweighting
        eta_vals = df_out.loc[mask, "cluster_Eta"].values
        eta_min, eta_max = -0.7, 0.7
        bins_eta = np.linspace(eta_min, eta_max, 21)
        
        hist_eta, bin_edges_eta = np.histogram(eta_vals, bins=bins_eta, density=True)
        hist_eta = hist_eta * 20
        x_centers_eta = 0.5 * (bin_edges_eta[:-1] + bin_edges_eta[1:])
        
        spline_eta = UnivariateSpline(x_centers_eta, hist_eta, s=0.0)
        pdf_eta_vals = spline_eta(eta_vals)
        pdf_eta_vals = np.clip(pdf_eta_vals, a_min=1e-3, a_max=None)
        
        weights_eta = 1.0 / pdf_eta_vals
        weights_eta /= np.mean(weights_eta)
        df_out.loc[mask, "weight"] *= weights_eta
        
        return df_out
    
    def apply_et_reweighting(self, df: pd.DataFrame, label: int) -> pd.DataFrame:
        """Apply Et reweighting for a specific class."""
        df_out = df.copy()
        mask = df_out[self.label_col] == label
        
        if not mask.any():
            return df_out
        
        if not self.reweight_config.get('et_reweight', False):
            return df_out
        
        # Et Reweighting
        et_vals = df_out.loc[mask, self.pt_col].values
        if len(et_vals) == 0:
            return df_out
        
        # Use the actual ET range from the data for binning
        et_min, et_max = et_vals.min(), et_vals.max()
        n_bins = self.reweight_config.get('et_reweight_bins', 20)
        
        et_bins = np.linspace(et_min, et_max, n_bins + 1)
        hist_et, bin_edges_et = np.histogram(et_vals, bins=et_bins, density=True)
        hist_et = hist_et * n_bins
        x_centers_et = 0.5 * (bin_edges_et[:-1] + bin_edges_et[1:])
        
        spline_et = UnivariateSpline(x_centers_et, hist_et, s=0.0)
        pdf_et_vals = spline_et(et_vals)
        pdf_et_vals = np.clip(pdf_et_vals, a_min=1e-3, a_max=None)
        
        weights_et = 1.0 / pdf_et_vals
        # Only use et_reweight_max to cap weights, not to define bin edges
        weight_cap = self.reweight_config.get("et_reweight_max", None)
        if weight_cap is not None:
            weights_et = np.clip(weights_et, a_min=None, a_max=weight_cap)
        weights_et /= np.mean(weights_et)
        df_out.loc[mask, "weight"] *= weights_et
        
        return df_out
    
    def apply_vertex_reweighting(self, df: pd.DataFrame, label: int) -> pd.DataFrame:
        """Apply vertex z reweighting for a specific class to flatten the distribution."""
        print(f"DEBUG: apply_vertex_reweighting called for label {label}")
        print(f"DEBUG: vertex_reweight config = {self.reweight_config.get('vertex_reweight', False)}")
        print(f"DEBUG: 'vertexz' in df = {'vertexz' in df.columns}")
        
        df_out = df.copy()
        mask = df_out[self.label_col] == label
        
        if not mask.any():
            print(f"DEBUG: No samples for label {label}")
            return df_out
        
        if not self.reweight_config.get('vertex_reweight', False):
            print(f"DEBUG: Vertex reweighting disabled, skipping")
            return df_out
        
        # Vertex z Reweighting
        vertex_vals = df_out.loc[mask, "vertexz"].values
        if len(vertex_vals) == 0:
            print(f"DEBUG: No vertex values for label {label}")
            return df_out
        
        print(f"DEBUG: ✓ Applying vertex reweighting for label {label}")
        print(f"DEBUG: Vertex range: {vertex_vals.min():.2f} to {vertex_vals.max():.2f} cm")
        print(f"DEBUG: Number of samples: {len(vertex_vals)}")
        
        # Use the actual vertex range from the data for binning
        vertex_min, vertex_max = vertex_vals.min(), vertex_vals.max()
        n_bins = self.reweight_config.get('vertex_reweight_bins', 20)
        
        vertex_bins = np.linspace(vertex_min, vertex_max, n_bins + 1)
        hist_vertex, bin_edges_vertex = np.histogram(vertex_vals, bins=vertex_bins, density=True)
        hist_vertex = hist_vertex * n_bins
        x_centers_vertex = 0.5 * (bin_edges_vertex[:-1] + bin_edges_vertex[1:])
        
        spline_vertex = UnivariateSpline(x_centers_vertex, hist_vertex, s=0.0)
        pdf_vertex_vals = spline_vertex(vertex_vals)
        pdf_vertex_vals = np.clip(pdf_vertex_vals, a_min=1e-3, a_max=None)
        
        weights_vertex = 1.0 / pdf_vertex_vals
        # Apply weight cap if specified
        weight_cap = self.reweight_config.get("vertex_reweight_max", None)
        if weight_cap is not None:
            weights_vertex = np.clip(weights_vertex, a_min=None, a_max=weight_cap)
        weights_vertex /= np.mean(weights_vertex)
        df_out.loc[mask, "weight"] *= weights_vertex
        
        return df_out
    
    def apply_all_reweighting(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply all reweighting steps."""
        df_out = self.apply_class_reweighting(df)
        
        # Per-class kinematic reweighting
        for label in [0, 1]:
            df_out = self.apply_eta_reweighting(df_out, label)
            df_out = self.apply_et_reweighting(df_out, label)
            df_out = self.apply_vertex_reweighting(df_out, label)
        
        return df_out
