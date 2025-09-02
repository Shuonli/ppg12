"""
Data loading utilities for binned training.
"""
import os
from typing import Dict, List, Optional, Union
import pandas as pd
import numpy as np

try:
    import uproot  # for ROOT files
except ImportError:
    uproot = None


class DataLoader:
    """Handles data loading from various file formats and binning."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.pt_col = config['data']['pt_column']
        self.iso_col = config['data']['iso_column']
        self.label_col = config['data']['label_column']
        self.configured_features = config['data'].get('features')
        self.use_single_file_set = config['data']['use_single_file_set']
        self._cached_single_df: Optional[pd.DataFrame] = None
        
        # Background subsampling configuration
        self.subsampling_config = config['data'].get('background_subsampling', {})
        self.subsampling_enabled = self.subsampling_config.get('enabled', False)
        self.et_threshold = self.subsampling_config.get('et_threshold', 10.0)
        self.subsample_fraction = self.subsampling_config.get('subsample_fraction', 0.5)  # Used as fallback for old method
        self.random_seed = self.subsampling_config.get('random_seed', 42)
        self.flatten_distribution = self.subsampling_config.get('flatten_distribution', True)
        self.n_et_bins = self.subsampling_config.get('n_et_bins', 20)  # Number of bins for flattening
        
        if self.use_single_file_set:
            self.file_map = config['data']['single_file_set']
        else:
            self.file_map = config['data']['per_bin_file_sets']
    
    def _flatten_et_distribution(self, df_low_et: pd.DataFrame) -> pd.DataFrame:
        """
        Flatten the ET distribution by subsampling each ET bin to have equal representation.
        
        Args:
            df_low_et: DataFrame containing low ET background samples
            
        Returns:
            DataFrame with flattened ET distribution
        """
        if len(df_low_et) == 0:
            return df_low_et
        
        # Create ET bins for the low ET region
        et_min = df_low_et[self.pt_col].min()
        et_max = self.et_threshold
        et_bins = np.linspace(et_min, et_max, self.n_et_bins + 1)
        
        # Assign bin indices to each sample
        df_low_et = df_low_et.copy()
        df_low_et['et_bin'] = pd.cut(df_low_et[self.pt_col], bins=et_bins, 
                                     labels=False, include_lowest=True)
        
        # Count samples in each bin
        bin_counts = df_low_et['et_bin'].value_counts().sort_index()
        
        # Find the minimum count (excluding empty bins)
        min_count = bin_counts[bin_counts > 0].min() if len(bin_counts[bin_counts > 0]) > 0 else 0
        
        if min_count == 0:
            print(f"Warning: No samples found in ET bins for flattening")
            return df_low_et.drop('et_bin', axis=1)
        
        # Subsample each bin to have min_count samples
        np.random.seed(self.random_seed)
        flattened_dfs = []
        
        for bin_idx in range(self.n_et_bins):
            bin_data = df_low_et[df_low_et['et_bin'] == bin_idx]
            if len(bin_data) > 0:
                if len(bin_data) >= min_count:
                    # Subsample to min_count
                    sampled_indices = np.random.choice(
                        bin_data.index, 
                        size=min_count, 
                        replace=False
                    )
                    flattened_dfs.append(bin_data.loc[sampled_indices])
                else:
                    # Keep all samples if fewer than min_count
                    flattened_dfs.append(bin_data)
        
        if flattened_dfs:
            result = pd.concat(flattened_dfs, ignore_index=True)
            result = result.drop('et_bin', axis=1)
            
            original_count = len(df_low_et)
            flattened_count = len(result)
            print(f"ET distribution flattening: {flattened_count}/{original_count} samples kept "
                  f"(<{self.et_threshold} GeV, {self.n_et_bins} bins, min_count={min_count})")
            
            return result
        else:
            return df_low_et.drop('et_bin', axis=1)
    
    def load_single(self, path: str, names: Optional[List[str]] = None) -> pd.DataFrame:
        """Load data from a single file of various formats."""
        ext = os.path.splitext(path)[1].lower()
        
        if ext in [".csv", ".txt"]:
            return pd.read_csv(path, sep=r"\s+", header=0, names=names)
        elif ext in [".pq", ".parquet"]:
            return pd.read_parquet(path)
        elif ext in [".root"]:
            if uproot is None:
                raise ImportError("uproot is required to read ROOT files. pip install uproot")
            
            # Heuristic: read first tree and all branches
            with uproot.open(path) as f:
                trees = [k for k in f.keys() if isinstance(f[k], uproot.behaviors.TTree.TTree)]
                if not trees:
                    trees = [k for k, v in f.items() if hasattr(v, "arrays")]
                key = trees[0] if trees else list(f.keys())[0]
                arrs = f[key].arrays(library="pd")
                return arrs
        else:
            raise ValueError(f"Unsupported file type: {ext}")
    
    def _load_single_file_set_once(self) -> pd.DataFrame:
        """Load and cache the entire single file set once, with labels applied.
        Subsequent calls reuse the cached DataFrame.
        """
        if self._cached_single_df is not None:
            return self._cached_single_df
        
        dfs: List[pd.DataFrame] = []
        file_map = self.file_map  # type: ignore[assignment]
        for file_type, paths in file_map.items():
            # Normalize to list
            if isinstance(paths, str):
                paths = [paths]
            for path in paths:
                if not os.path.exists(path):
                    raise FileNotFoundError(f"Configured file not found: {path}")
                df_i = self.load_single(path, names=None)
                if file_type == "signal":
                    df_i[self.label_col] = 1
                    df_i = df_i[df_i["pid"].isin([1, 2])]  # photon only
                elif file_type == "background":
                    df_i[self.label_col] = 0
                    df_i = df_i[~df_i["pid"].isin([1, 2])]  # reject electrons
                    
                    # Apply background subsampling if enabled
                    if self.subsampling_enabled and self.pt_col in df_i.columns:
                        # Split background into low ET and high ET based on threshold
                        low_et_mask = df_i[self.pt_col] < self.et_threshold
                        high_et_mask = df_i[self.pt_col] >= self.et_threshold
                        
                        df_low_et = df_i[low_et_mask]
                        df_high_et = df_i[high_et_mask]
                        
                        # Apply flattening or fixed-fraction subsampling to low ET background
                        if len(df_low_et) > 0:
                            if self.flatten_distribution:
                                # Flatten the ET distribution in the low ET region
                                df_low_et_processed = self._flatten_et_distribution(df_low_et)
                            else:
                                # Use fixed-fraction subsampling (legacy method)
                                np.random.seed(self.random_seed)
                                subsample_size = int(len(df_low_et) * self.subsample_fraction)
                                subsample_indices = np.random.choice(
                                    df_low_et.index, 
                                    size=subsample_size, 
                                    replace=False
                                )
                                df_low_et_processed = df_low_et.loc[subsample_indices]
                                print(f"Background subsampling: kept {subsample_size}/{len(df_low_et)} "
                                      f"low ET samples (<{self.et_threshold} GeV), fraction={self.subsample_fraction}")
                            
                            # Combine processed low ET with all high ET background
                            df_i = pd.concat([df_low_et_processed, df_high_et], ignore_index=True)
                        else:
                            # No low ET samples to process
                            df_i = df_high_et
                else:
                    # Unknown type, skip
                    continue
                dfs.append(df_i)
        if not dfs:
            raise RuntimeError("No data loaded from single_file_set")
        df_all = pd.concat(dfs, ignore_index=True)
        
        # Validate required columns
        base_required = {self.pt_col, self.iso_col, self.label_col}
        missing = sorted(base_required - set(df_all.columns))
        if missing:
            raise KeyError(f"Missing required columns in single_file_set: {missing}")
        
        # Drop rows with missing label and key columns
        df_all = df_all.dropna(subset=[self.label_col])
        
        # Ensure numeric dtypes for key columns
        for c in [self.pt_col, self.iso_col, self.label_col]:
            if c in df_all.columns:
                df_all[c] = pd.to_numeric(df_all[c], errors="coerce")
        # Drop rows with NA on key columns after coercion
        df_all = df_all.dropna(subset=[self.pt_col, self.iso_col, self.label_col])
        
        self._cached_single_df = df_all
        return self._cached_single_df
    
    def load_bin_data(self, bin_label: str, et_min: float, et_max: float) -> pd.DataFrame:
        """Load data for a specific energy bin."""
        if self.use_single_file_set:
            # Slice from cached full dataset
            df_all = self._load_single_file_set_once()
            df = df_all[(df_all[self.pt_col] >= et_min) & (df_all[self.pt_col] < et_max)]
            return df
        
        file_map = self.file_map.get(bin_label, {})
        if not file_map:
            raise FileNotFoundError(f"No input files configured for bin {bin_label}")
        
        dfs = []
        for file_type, paths in file_map.items():
            # Normalize to list
            if isinstance(paths, str):
                paths = [paths]
            
            for path in paths:
                if not os.path.exists(path):
                    raise FileNotFoundError(f"Configured file not found: {path}")
                
                df_i = self.load_single(path, names=None)
                
                if file_type == "signal":
                    df_i[self.label_col] = 1
                    df_i = df_i[df_i["pid"].isin([1, 2])]  # photon only
                elif file_type == "background":
                    df_i[self.label_col] = 0
                    df_i = df_i[~df_i["pid"].isin([1, 2])]  # reject electrons
                    
                    # Apply background subsampling if enabled
                    if self.subsampling_enabled and self.pt_col in df_i.columns:
                        # Split background into low ET and high ET based on threshold
                        low_et_mask = df_i[self.pt_col] < self.et_threshold
                        high_et_mask = df_i[self.pt_col] >= self.et_threshold
                        
                        df_low_et = df_i[low_et_mask]
                        df_high_et = df_i[high_et_mask]
                        
                        # Apply flattening or fixed-fraction subsampling to low ET background
                        if len(df_low_et) > 0:
                            if self.flatten_distribution:
                                # Flatten the ET distribution in the low ET region
                                df_low_et_processed = self._flatten_et_distribution(df_low_et)
                            else:
                                # Use fixed-fraction subsampling (legacy method)
                                np.random.seed(self.random_seed)
                                subsample_size = int(len(df_low_et) * self.subsample_fraction)
                                subsample_indices = np.random.choice(
                                    df_low_et.index, 
                                    size=subsample_size, 
                                    replace=False
                                )
                                df_low_et_processed = df_low_et.loc[subsample_indices]
                                print(f"Background subsampling (bin {bin_label}): kept {subsample_size}/{len(df_low_et)} "
                                      f"low ET samples (<{self.et_threshold} GeV), fraction={self.subsample_fraction}")
                            
                            # Combine processed low ET with all high ET background
                            df_i = pd.concat([df_low_et_processed, df_high_et], ignore_index=True)
                        else:
                            # No low ET samples to process
                            df_i = df_high_et
                else:
                    # Unknown type, skip
                    continue
                
                dfs.append(df_i)
        
        if not dfs:
            raise RuntimeError(f"No data loaded for bin {bin_label}")
        
        df = pd.concat(dfs, ignore_index=True)
        
        # Validate required columns
        base_required = {self.pt_col, self.iso_col, self.label_col}
        missing = sorted(base_required - set(df.columns))
        if missing:
            raise KeyError(f"Missing required columns for bin {bin_label}: {missing}")
        
        # Drop rows with missing label and key columns
        df = df.dropna(subset=[self.label_col])
        
        # Filter by Et bin
        df = df[(df[self.pt_col] >= et_min) & (df[self.pt_col] < et_max)]
        
        # Ensure numeric dtypes for key columns
        for c in [self.pt_col, self.iso_col, self.label_col]:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")
        
        # Drop rows with NA on key columns after coercion
        df = df.dropna(subset=[self.pt_col, self.iso_col, self.label_col])
        
        return df
    
    def autodetect_features(self, df: pd.DataFrame) -> List[str]:
        """Return configured features if provided; otherwise auto-detect numeric features excluding PT/ISO/LABEL."""
        if self.configured_features:
            # Keep only features present in the dataframe, preserve order
            return [f for f in self.configured_features if f in df.columns]
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        exclude = {self.pt_col, self.iso_col, self.label_col}
        feats = [c for c in numeric_cols if c not in exclude]
        return feats
    
    def get_all_bin_data(self, bin_edges: List[float], bin_labels: List[str]) -> List[pd.DataFrame]:
        """Load data for all bins."""
        all_dfs = []
        if self.use_single_file_set:
            df_all = self._load_single_file_set_once()
            for i in range(len(bin_edges) - 1):
                et_min, et_max = bin_edges[i], bin_edges[i+1]
                bin_label = bin_labels[i]
                print(f"Loading data for bin {bin_label} (slicing cached single_file_set)...")
                df_bin = df_all[(df_all[self.pt_col] >= et_min) & (df_all[self.pt_col] < et_max)]
                all_dfs.append(df_bin)
            return all_dfs
        
        for i in range(len(bin_edges) - 1):
            et_min, et_max = bin_edges[i], bin_edges[i+1]
            bin_label = bin_labels[i]
            print(f"Loading data for bin {bin_label}...")
            df_bin = self.load_bin_data(bin_label, et_min, et_max)
            all_dfs.append(df_bin)
        return all_dfs
