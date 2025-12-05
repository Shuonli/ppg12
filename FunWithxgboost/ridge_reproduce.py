#!/usr/bin/env python
"""
Reproduce Kelly, Malamud, Zhou (2024) "The Virtue of Complexity in Return Prediction"
Journal of Finance, Vol. LXXIX, No. 1

This script implements:
1. Random Fourier Features (RFF) as in Rahimi & Recht (2007)
2. Ridge/ridgeless regression for market return prediction
3. Out-of-sample recursive forecasting with rolling windows
4. Performance metrics: R², Sharpe ratio, alpha, information ratio

Based on Section V of the paper (Empirical Analysis).
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Tuple, Dict, Optional, Sequence
import math
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')


class RandomFourierFeatures:
    """
    Random Fourier Features as in Rahimi & Recht (2007).
    
    Transforms input predictors G_t into features:
    S_{i,t} = [sin(γ * ω_i' * G_t), cos(γ * ω_i' * G_t)]'
    
    where ω_i ~ N(0, I) are random weights.
    """
    
    def __init__(self, n_features: int, n_predictors: int, gamma: float = 2.0, seed: int = 42):
        """
        Args:
            n_features: Number of RFF pairs to generate (P/2 in paper notation)
            n_predictors: Dimension of input predictors (15 in paper)
            gamma: Bandwidth parameter (γ in paper, default=2)
            seed: Random seed for reproducibility
        """
        self.n_features = n_features
        self.n_predictors = n_predictors
        self.gamma = gamma
        
        # Generate random weights ω_i ~ N(0, I)
        rng = np.random.RandomState(seed)
        self.omega = rng.randn(n_features, n_predictors)
    
    def transform(self, X: np.ndarray) -> np.ndarray:
        """
        Transform predictors into RFF.
        
        Args:
            X: (n_samples, n_predictors) array of predictors
            
        Returns:
            (n_samples, 2*n_features) array of RFF
        """
        # Compute γ * ω' * X'  (shape: n_features x n_samples)
        projections = self.gamma * self.omega @ X.T
        
        # Stack [sin, cos] features
        sin_features = np.sin(projections).T
        cos_features = np.cos(projections).T
        
        return np.hstack([sin_features, cos_features])


class RidgeRegression:
    """Ridge regression with optional ridgeless limit."""
    
    def __init__(self, alpha: float = 1e-3):
        """
        Args:
            alpha: Ridge penalty (z in paper). Use alpha=1e-10 for ridgeless.
        """
        self.alpha = alpha
        self.beta = None
    
    def fit(self, X: np.ndarray, y: np.ndarray):
        """
        Fit ridge regression: β = (α*I + X'X)^{-1} X'y
        
        Args:
            X: (n_samples, n_features) feature matrix
            y: (n_samples,) target vector
        """
        n_samples, n_features = X.shape
        
        # Compute X'X + α*I
        XtX = X.T @ X / n_samples
        ridge_matrix = self.alpha * np.eye(n_features) + XtX
        
        # Solve for β
        Xty = X.T @ y / n_samples
        self.beta = np.linalg.solve(ridge_matrix, Xty)
        
        return self
    
    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict y = X @ β"""
        return X @ self.beta


def standardize_features(X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Standardize features by their training sample standard deviation.
    
    Returns:
        X_std: Standardized features
        std: Standard deviations used
    """
    std = np.std(X, axis=0, ddof=1)
    std[std == 0] = 1.0  # Avoid division by zero
    return X / std, std


def compute_oos_metrics(y_true: np.ndarray, y_pred: np.ndarray, 
                        positions: np.ndarray) -> Dict[str, float]:
    """
    Compute out-of-sample performance metrics.
    
    Args:
        y_true: Realized returns
        y_pred: Predicted returns
        positions: Market timing positions (forecasts)
        
    Returns:
        Dictionary of performance metrics
    """
    # Forecast errors
    errors = y_true - y_pred
    mse = np.mean(errors**2)
    
    # Out-of-sample R²
    var_y = np.var(y_true, ddof=1)
    oos_r2 = 1 - mse / var_y
    
    # Market timing strategy returns
    strategy_returns = positions * y_true
    
    # Strategy performance
    mean_return = np.mean(strategy_returns)
    std_return = np.std(strategy_returns, ddof=1)
    sharpe_ratio = mean_return / std_return if std_return > 0 else 0
    
    # Alpha vs. buy-and-hold (untimed market)
    # Regress strategy returns on market returns
    X_reg = np.column_stack([np.ones(len(y_true)), y_true])
    try:
        coefs = np.linalg.lstsq(X_reg, strategy_returns, rcond=None)[0]
        alpha = coefs[0]
        beta_mkt = coefs[1]
        
        # Residuals for information ratio
        fitted = X_reg @ coefs
        residuals = strategy_returns - fitted
        tracking_error = np.std(residuals, ddof=1)
        information_ratio = alpha / tracking_error if tracking_error > 0 else 0
    except:
        alpha = np.nan
        information_ratio = np.nan
    
    return {
        'oos_r2': oos_r2,
        'mean_return': mean_return,
        'volatility': std_return,
        'sharpe_ratio': sharpe_ratio,
        'alpha': alpha,
        'information_ratio': information_ratio,
        'beta_norm': np.linalg.norm(positions) / np.sqrt(len(positions))
    }


def recursive_oos_forecast(df: pd.DataFrame, 
                          n_rff: int = 100,
                          train_window: int = 12,
                          ridge_alpha: float = 1e-3,
                          gamma: float = 2.0,
                          seed: int = 42,
                          phi_matrix: Optional[np.ndarray] = None) -> Dict[str, np.ndarray]:
    """
    Perform recursive out-of-sample forecasting as in the paper.
    
    Args:
        df: DataFrame with standardized predictors and returns
        n_rff: Number of RFF pairs (P/2)
        train_window: Rolling training window in months (T)
        ridge_alpha: Ridge penalty (z in paper)
        gamma: RFF bandwidth parameter
        seed: Random seed
        
    Returns:
        Dictionary with forecasts, realized returns, and positions
    """
    # Extract predictor columns (all *_std columns)
    pred_cols = [col for col in df.columns if col.endswith('_std')]
    X_raw = df[pred_cols].values
    y_raw = df['Mkt-RF_std'].values
    
    n_samples = len(df)
    n_predictors = len(pred_cols)
    
    # Initialize or reuse RFF features
    if phi_matrix is None:
        rff = RandomFourierFeatures(n_rff, n_predictors, gamma=gamma, seed=seed)
        phi_matrix = rff.transform(X_raw)
    else:
        if phi_matrix.shape[0] != n_samples:
            raise ValueError("Precomputed RFF matrix row count does not match sample size.")
        if phi_matrix.shape[1] != 2 * n_rff:
            raise ValueError("Precomputed RFF matrix feature count does not match n_rff.")
    
    # Storage for out-of-sample results
    forecasts = []
    realized = []
    positions = []
    
    # Recursive forecasting
    for t in tqdm(range(train_window, n_samples),
                  total=n_samples - train_window,
                  desc="Rolling forecast",
                  leave=False):
        # Training window: [t-train_window, t)
        train_idx = slice(t - train_window, t)
        y_train = y_raw[train_idx]
        
        # Test point: t
        y_test = y_raw[t]
        
        # Use cached RFF features
        X_train = phi_matrix[train_idx]
        X_test = phi_matrix[t:t+1]
        
        # Standardize features by training sample std
        X_train_std, train_std = standardize_features(X_train)
        X_test_std = X_test / train_std
        
        # Fit ridge regression
        model = RidgeRegression(alpha=ridge_alpha)
        model.fit(X_train_std, y_train)
        
        # Forecast
        y_pred = model.predict(X_test_std)[0]
        
        # Store results
        forecasts.append(y_pred)
        realized.append(y_test)
        positions.append(y_pred)  # Market timing position = forecast
    
    return {
        'forecasts': np.array(forecasts),
        'realized': np.array(realized),
        'positions': np.array(positions),
        'dates': df.index[train_window:].values
    }


def volatility_timed_momentum(df: pd.DataFrame,
                              momentum_window: int = 12,
                              vol_window: int = 12,
                              target_vol: float = 1.0) -> Dict[str, np.ndarray]:
    """
    Volatility-managed momentum strategy (Moreira & Muir style).
    
    Uses past momentum scaled by inverse variance to decide positions.
    """
    returns = df['Mkt-RF_std'].values
    dates = df.index.values
    
    momentum = pd.Series(returns).rolling(momentum_window).mean().shift(1)
    volatility = pd.Series(returns).rolling(vol_window).std(ddof=1).shift(1)
    
    start = max(momentum_window, vol_window)
    momentum = momentum.iloc[start:].to_numpy()
    volatility = volatility.iloc[start:].to_numpy()
    realized = returns[start:]
    signal = np.zeros_like(realized)
    
    valid = (~np.isnan(momentum)) & (~np.isnan(volatility)) & (volatility > 0)
    denom = volatility[valid] ** 2 + 1e-8
    signal[valid] = momentum[valid] / denom
    
    # Scale to target volatility
    strat_returns = signal[valid] * realized[valid] if np.any(valid) else np.array([])
    current_vol = np.std(strat_returns, ddof=1) if strat_returns.size > 1 else 0
    if current_vol > 0 and target_vol > 0:
        scale = target_vol / current_vol
        signal *= scale
    
    return {
        'forecasts': signal,
        'realized': realized,
        'positions': signal,
        'dates': dates[start:]
    }


def gaussian_kernel_regression(df: pd.DataFrame,
                               train_window: int = 12,
                               bandwidth: float = 1.0) -> Dict[str, np.ndarray]:
    """
    Gaussian kernel regression baseline using standardized predictors.
    """
    if bandwidth <= 0:
        raise ValueError("bandwidth must be positive.")
    
    pred_cols = [col for col in df.columns if col.endswith('_std')]
    X = df[pred_cols].values
    y = df['Mkt-RF_std'].values
    
    forecasts = []
    realized = []
    positions = []
    
    iterator = tqdm(
        range(train_window, len(df)),
        total=len(df) - train_window,
        desc="MLP Baseline",
        leave=False
    )
    for t in iterator:
        train_idx = slice(t - train_window, t)
        X_train = X[train_idx]
        y_train = y[train_idx]
        x_test = X[t:t+1]
        
        X_train_std, train_std = standardize_features(X_train)
        x_test_std = x_test / train_std
        
        diffs = X_train_std - x_test_std
        dist_sq = np.sum(diffs**2, axis=1)
        weights = np.exp(-0.5 * dist_sq / (bandwidth ** 2))
        weight_sum = np.sum(weights)
        
        if weight_sum <= 0:
            y_pred = 0.0
        else:
            y_pred = np.dot(weights, y_train) / weight_sum
        
        forecasts.append(y_pred)
        realized.append(y[t])
        positions.append(y_pred)
    
    return {
        'forecasts': np.array(forecasts),
        'realized': np.array(realized),
        'positions': np.array(positions),
        'dates': df.index[train_window:].values
    }


def _require_torch():
    try:
        import torch  # type: ignore
        from torch import nn  # noqa: F401
    except ImportError as exc:
        raise ImportError("PyTorch is required for MLP/SIREN baselines. "
                          "Install via `pip install torch`.") from exc
    return torch


def _standardize_window_features(features: np.ndarray, train_idx: slice, test_idx: int):
    X_train = features[train_idx]
    X_test = features[test_idx:test_idx + 1]
    X_train_std, train_std = standardize_features(X_train)
    X_test_std = X_test / train_std
    return X_train_std, X_test_std


def _prepare_features(df: pd.DataFrame,
                      use_rff: bool,
                      n_rff: int,
                      gamma: float,
                      seed: int) -> np.ndarray:
    pred_cols = [col for col in df.columns if col.endswith('_std')]
    X_raw = df[pred_cols].values
    if use_rff:
        rff = RandomFourierFeatures(n_rff, len(pred_cols), gamma=gamma, seed=seed)
        return rff.transform(X_raw)
    return X_raw


def mlp_recursive_baseline(df: pd.DataFrame,
                           train_window: int = 12,
                           hidden_sizes: Sequence[int] = (64, 32),
                           epochs: int = 200,
                           lr: float = 1e-3,
                           use_rff: bool = False,
                           n_rff: int = 100,
                           gamma: float = 2.0,
                           seed: int = 0) -> Dict[str, np.ndarray]:
    """
    Vanilla MLP regression baseline with optional RFF inputs.
    """
    torch = _require_torch()
    from torch import nn
    
    torch.manual_seed(seed)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    y = df['Mkt-RF_std'].values
    features = _prepare_features(df, use_rff, n_rff, gamma, seed).astype(np.float32)
    
    forecasts, realized, positions = [], [], []
    
    iterator = tqdm(
        range(train_window, len(df)),
        total=len(df) - train_window,
        desc="MLP Baseline",
        leave=False
    )
    for t in iterator:
        train_idx = slice(t - train_window, t)
        X_train_std, X_test_std = _standardize_window_features(features, train_idx, t)
        y_train = y[train_idx].astype(np.float32)
        
        input_dim = X_train_std.shape[1]
        layers = []
        prev = input_dim
        for h in hidden_sizes:
            layers.append(nn.Linear(prev, h))
            nn.init.xavier_uniform_(layers[-1].weight)
            layers.append(nn.ReLU())
            prev = h
        layers.append(nn.Linear(prev, 1))
        nn.init.xavier_uniform_(layers[-1].weight)
        model = nn.Sequential(*layers).to(device)
        
        optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        loss_fn = nn.MSELoss()
        
        X_train_tensor = torch.from_numpy(X_train_std.astype(np.float32)).to(device)
        y_train_tensor = torch.from_numpy(y_train).to(device)
        X_test_tensor = torch.from_numpy(X_test_std.astype(np.float32)).to(device)
        
        for _ in range(epochs):
            optimizer.zero_grad()
            preds = model(X_train_tensor).squeeze()
            loss = loss_fn(preds, y_train_tensor)
            loss.backward()
            optimizer.step()
        
        model.eval()
        with torch.no_grad():
            y_pred = model(X_test_tensor).item()
        
        forecasts.append(y_pred)
        realized.append(y[t])
        positions.append(y_pred)
    
    return {
        'forecasts': np.array(forecasts),
        'realized': np.array(realized),
        'positions': np.array(positions),
        'dates': df.index[train_window:].values
    }


def siren_recursive_baseline(df: pd.DataFrame,
                             train_window: int = 12,
                             hidden_dim: int = 64,
                             hidden_layers: int = 3,
                             omega_0: float = 30.0,
                             epochs: int = 300,
                             lr: float = 1e-3,
                             use_rff: bool = False,
                             n_rff: int = 100,
                             gamma: float = 2.0,
                             seed: int = 0) -> Dict[str, np.ndarray]:
    """
    SIREN regression baseline with optional RFF inputs.
    """
    torch = _require_torch()
    from torch import nn
    
    torch.manual_seed(seed)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    y = df['Mkt-RF_std'].values
    features = _prepare_features(df, use_rff, n_rff, gamma, seed).astype(np.float32)
    
    class SineLayer(nn.Module):
        def __init__(self, in_features, out_features, w0=1.0, is_first=False):
            super().__init__()
            self.linear = nn.Linear(in_features, out_features)
            self.w0 = w0
            self.is_first = is_first
            self.init_weights()
        
        def init_weights(self):
            with torch.no_grad():
                if self.is_first:
                    bound = 1 / self.linear.in_features
                else:
                    bound = math.sqrt(6 / self.linear.in_features) / self.w0
                self.linear.weight.uniform_(-bound, bound)
        
        def forward(self, x):
            return torch.sin(self.w0 * self.linear(x))
    
    class SirenNet(nn.Module):
        def __init__(self, input_dim, hidden_dim, hidden_layers, w0):
            super().__init__()
            layers_list = [SineLayer(input_dim, hidden_dim, w0=w0, is_first=True)]
            for _ in range(hidden_layers - 1):
                layers_list.append(SineLayer(hidden_dim, hidden_dim, w0=w0))
            self.layers = nn.ModuleList(layers_list)
            self.final_linear = nn.Linear(hidden_dim, 1)
            with torch.no_grad():
                bound = math.sqrt(6 / hidden_dim) / w0
                self.final_linear.weight.uniform_(-bound, bound)
        
        def forward(self, x):
            for layer in self.layers:
                x = layer(x)
            return self.final_linear(x)
    
    forecasts, realized, positions = [], [], []

    
    iterator = tqdm(
        range(train_window, len(df)),
        total=len(df) - train_window,
        desc="SIREN Baseline",
        leave=False
    )
    for t in iterator:
        train_idx = slice(t - train_window, t)
        X_train_std, X_test_std = _standardize_window_features(features, train_idx, t)
        y_train = y[train_idx].astype(np.float32)
        
        input_dim = X_train_std.shape[1]
        model = SirenNet(input_dim, hidden_dim, hidden_layers, omega_0).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        loss_fn = nn.MSELoss()
        
        X_train_tensor = torch.from_numpy(X_train_std.astype(np.float32)).to(device)
        y_train_tensor = torch.from_numpy(y_train).to(device)
        X_test_tensor = torch.from_numpy(X_test_std.astype(np.float32)).to(device)
        
        for _ in range(epochs):
            optimizer.zero_grad()
            preds = model(X_train_tensor).squeeze()
            loss = loss_fn(preds, y_train_tensor)
            loss.backward()
            optimizer.step()
        
        model.eval()
        with torch.no_grad():
            y_pred = model(X_test_tensor).item()
        
        forecasts.append(y_pred)
        realized.append(y[t])
        positions.append(y_pred)
    
    return {
        'forecasts': np.array(forecasts),
        'realized': np.array(realized),
        'positions': np.array(positions),
        'dates': df.index[train_window:].values
    }


def run_complexity_sweep(df: pd.DataFrame,
                        train_window: int = 12,
                        ridge_alphas: list = None,
                        n_rff_list: list = None,
                        gamma: float = 2.0,
                        n_seeds: int = 10) -> pd.DataFrame:
    """
    Sweep over model complexity (P) and ridge penalty (z) as in paper Figures 7-8.
    
    Args:
        df: DataFrame with standardized data
        train_window: Training window size (T)
        ridge_alphas: List of ridge penalties to try
        n_rff_list: List of RFF counts to try (P/2 values)
        gamma: RFF bandwidth
        n_seeds: Number of random seeds to average over
        
    Returns:
        DataFrame with performance metrics for each (complexity, ridge_alpha) pair
    """
    if ridge_alphas is None:
        # Paper uses log10(z) ∈ {-3, -2, -1, 0, 1, 2, 3}
        ridge_alphas = [10**i for i in range(-3, 4)]
    
    if n_rff_list is None:
        # Paper uses P from 2 to 12,000 (so n_rff from 1 to 6,000)
        # This creates complexity c = P/T up to 1000 for T=12
        n_rff_list = [1, 2, 5, 10, 20, 50, 100, 200, 500]
    
    pred_cols = [col for col in df.columns if col.endswith('_std')]
    n_predictors = len(pred_cols)
    X_raw = df[pred_cols].values
    
    results = []
    phi_cache: Dict[int, Dict[int, np.ndarray]] = {}
    
    total_runs = len(ridge_alphas) * len(n_rff_list)
    run_count = 0
    
    for ridge_alpha in ridge_alphas:
        for n_rff in n_rff_list:
            seed_cache = phi_cache.setdefault(n_rff, {})
            run_count += 1
            print(f"Progress: {run_count}/{total_runs} | "
                  f"T={train_window}, P={2*n_rff}, z={ridge_alpha:.1e}")
            
            # Average over multiple random seeds (as in paper)
            metrics_list = []
            
            for seed in range(n_seeds):
                try:
                    if seed not in seed_cache:
                        rff = RandomFourierFeatures(
                            n_rff, n_predictors, gamma=gamma, seed=seed
                        )
                        seed_cache[seed] = rff.transform(X_raw)
                    phi_matrix = seed_cache[seed]
                    
                    # Run out-of-sample forecast
                    oos_results = recursive_oos_forecast(
                        df, n_rff=n_rff, train_window=train_window,
                        ridge_alpha=ridge_alpha, gamma=gamma, seed=seed,
                        phi_matrix=phi_matrix
                    )
                    
                    # Compute metrics
                    metrics = compute_oos_metrics(
                        oos_results['realized'],
                        oos_results['forecasts'],
                        oos_results['positions']
                    )
                    metrics_list.append(metrics)
                    
                except Exception as e:
                    print(f"\nWarning: Failed for n_rff={n_rff}, alpha={ridge_alpha}, seed={seed}: {e}")
                    continue
            
            if metrics_list:
                # Average metrics across seeds
                avg_metrics = {
                    key: np.mean([m[key] for m in metrics_list])
                    for key in metrics_list[0].keys()
                }
                
                # Store results
                results.append({
                    'train_window': train_window,
                    'n_features': 2 * n_rff,  # P in paper
                    'complexity': (2 * n_rff) / train_window,  # c = P/T
                    'ridge_alpha': ridge_alpha,
                    'log10_alpha': np.log10(ridge_alpha),
                    **avg_metrics
                })
    
    print()  # New line after progress
    return pd.DataFrame(results)


def main():
    """Main reproduction script."""
    
    print("=" * 80)
    print("Kelly, Malamud, Zhou (2024) - Reproduction Script")
    print("The Virtue of Complexity in Return Prediction")
    print("=" * 80)
    print()
    
    # Load prepared data
    data_path = Path('ridge_panel.csv')
    if not data_path.exists():
        print(f"Error: {data_path} not found. Run ridge_data_prep.py first.")
        return
    
    df = pd.read_csv(data_path, index_col=0)
    print(f"Loaded data: {len(df)} months from {df.index[0]} to {df.index[-1]}")
    print()
    
    # ========================================================================
    # Example 1: Single model evaluation (quick test)
    # ========================================================================
    print("Example 1: Single Model Evaluation")
    print("-" * 80)
    
    # Paper's main results use T=12, P=12000, z=10^3
    n_rff = 50  # P = 100 (smaller for quick demo)
    train_window = 12
    ridge_alpha = 1e3
    n_rff_list = [1, 2, 5, 10, 20, 50, 100, 200]  # P = 2 * n_rff
    rff_extremes = [min(n_rff_list), max(n_rff_list)]
    baseline_records = []
    
    def evaluate_baseline(title: str, method_label: str, run_fn):
        print(title)
        print("-" * 80)
        try:
            results = run_fn()
        except ImportError as exc:
            print(f"  Skipped: {exc}")
            print()
            return
        metrics = compute_oos_metrics(
            results['realized'],
            results['forecasts'],
            results['positions']
        )
        print(f"  R²:                {metrics['oos_r2']:.4f}")
        print(f"  Mean Return:       {metrics['mean_return']:.4f}")
        print(f"  Volatility:        {metrics['volatility']:.4f}")
        print(f"  Sharpe Ratio:      {metrics['sharpe_ratio']:.4f}")
        print(f"  Alpha:             {metrics['alpha']:.4f}")
        print(f"  Information Ratio: {metrics['information_ratio']:.4f}")
        print()
        baseline_records.append({
            'method': method_label,
            **metrics
        })
    
    print(f"Configuration: T={train_window}, P={2*n_rff}, z={ridge_alpha:.1e}")
    
    oos_results = recursive_oos_forecast(
        df, n_rff=n_rff, train_window=train_window,
        ridge_alpha=ridge_alpha, gamma=2.0, seed=42
    )
    
    metrics = compute_oos_metrics(
        oos_results['realized'],
        oos_results['forecasts'],
        oos_results['positions']
    )
    
    print("\nOut-of-Sample Performance:")
    print(f"  R²:                {metrics['oos_r2']:.4f}")
    print(f"  Mean Return:       {metrics['mean_return']:.4f}")
    print(f"  Volatility:        {metrics['volatility']:.4f}")
    print(f"  Sharpe Ratio:      {metrics['sharpe_ratio']:.4f}")
    print(f"  Alpha:             {metrics['alpha']:.4f}")
    print(f"  Information Ratio: {metrics['information_ratio']:.4f}")
    print()
    
    # ========================================================================
    # Baseline: Volatility-timed momentum
    # ========================================================================
    evaluate_baseline(
        "Baseline: Volatility-Timed Momentum Strategy",
        "Volatility-Timed Momentum",
        lambda: volatility_timed_momentum(df, momentum_window=12, vol_window=12, target_vol=1.0)
    )
    
    # ========================================================================
    # Baseline: Gaussian kernel regression
    # ========================================================================
    evaluate_baseline(
        "Baseline: Gaussian Kernel Regression",
        "Gaussian Kernel Regression",
        lambda: gaussian_kernel_regression(df, train_window=12, bandwidth=1.0)
    )
    
    evaluate_baseline(
        "Baseline: Vanilla MLP (Raw Predictors)",
        "MLP (Raw)",
        lambda: mlp_recursive_baseline(df, train_window=12, use_rff=False, hidden_sizes=(64, 32),
                                       epochs=150, lr=1e-3, seed=42)
    )
    
    for rff_val in rff_extremes:
        evaluate_baseline(
            f"Baseline: Vanilla MLP (RFF Input, P={2*rff_val})",
            f"MLP (RFF, P={2*rff_val})",
            lambda rff=rff_val: mlp_recursive_baseline(
                df, train_window=12, use_rff=True, n_rff=rff, gamma=2.0,
                hidden_sizes=(128, 64), epochs=150, lr=1e-3, seed=42
            )
        )
    
    evaluate_baseline(
        "Baseline: SIREN (Raw Predictors)",
        "SIREN (Raw)",
        lambda: siren_recursive_baseline(df, train_window=12, use_rff=False,
                                         hidden_dim=64, hidden_layers=3, epochs=200,
                                         lr=1e-3, seed=42)
    )
    
    for rff_val in rff_extremes:
        evaluate_baseline(
            f"Baseline: SIREN (RFF Input, P={2*rff_val})",
            f"SIREN (RFF, P={2*rff_val})",
            lambda rff=rff_val: siren_recursive_baseline(
                df, train_window=12, use_rff=True, n_rff=rff,
                gamma=2.0, hidden_dim=64, hidden_layers=3,
                epochs=200, lr=1e-3, seed=42
            )
        )
    
    if baseline_records:
        baseline_path = Path('ridge_baselines.csv')
        pd.DataFrame(baseline_records).to_csv(baseline_path, index=False)
        print(f"Baseline metrics saved to: {baseline_path}")
        print()
    
    # ========================================================================
    # Example 2: Complexity sweep (reproduce paper Figures 7-8)
    # ========================================================================
    print("\nExample 2: Complexity Sweep (Virtue of Complexity)")
    print("-" * 80)
    print("This reproduces the analysis in paper Section V.C")
    print("Sweeping over model complexity (P) and ridge penalty (z)")
    print()
    
    # Full reproduction matching paper
    # Paper: P ∈ [2, 12000], log10(z) ∈ [-3, 3]
    # This achieves complexity c up to 1000 for T=12
    
    ridge_alphas = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
    
    results_df = run_complexity_sweep(
        df,
        train_window=12,
        ridge_alphas=ridge_alphas,
        n_rff_list=n_rff_list,
        gamma=2.0,
        n_seeds=5  # Paper uses 1000; increase for more stable estimates
    )
    
    # Save results
    output_path = Path('ridge_results.csv')
    results_df.to_csv(output_path, index=False)
    print(f"\nResults saved to: {output_path}")
    print()
    
    # Display summary
    print("Summary by Ridge Penalty:")
    print("-" * 80)
    summary = results_df.groupby('log10_alpha').agg({
        'oos_r2': 'mean',
        'sharpe_ratio': 'mean',
        'information_ratio': 'mean'
    }).round(4)
    print(summary)
    print()
    
    print("Summary by Complexity (c = P/T):")
    print("-" * 80)
    # Bin complexity for readability
    results_df['complexity_bin'] = pd.cut(results_df['complexity'], 
                                          bins=[0, 0.5, 1.0, 2.0, 10.0, 100.0],
                                          labels=['<0.5', '0.5-1', '1-2', '2-10', '>10'])
    summary2 = results_df.groupby('complexity_bin').agg({
        'oos_r2': 'mean',
        'sharpe_ratio': 'mean',
        'information_ratio': 'mean'
    }).round(4)
    print(summary2)
    print()
    
    # Key findings
    print("Key Findings (matching paper):")
    print("-" * 80)
    
    # Find best Sharpe ratio
    best_idx = results_df['sharpe_ratio'].idxmax()
    best = results_df.loc[best_idx]
    print(f"Best Sharpe Ratio: {best['sharpe_ratio']:.4f}")
    print(f"  at P={best['n_features']:.0f}, z={best['ridge_alpha']:.1e}")
    print(f"  (complexity c = {best['complexity']:.2f})")
    print()
    
    # Compare ridgeless vs. high regularization
    ridgeless = results_df[results_df['ridge_alpha'] == 1e-3]
    high_reg = results_df[results_df['ridge_alpha'] == 1e3]
    
    print(f"Ridgeless (z=1e-3) avg Sharpe:     {ridgeless['sharpe_ratio'].mean():.4f}")
    print(f"High Ridge (z=1e3) avg Sharpe:     {high_reg['sharpe_ratio'].mean():.4f}")
    print()
    
    print("=" * 80)
    print("Reproduction complete!")
    print()
    print("To reproduce paper figures, use the saved ridge_results.csv with:")
    print("  - Plot Sharpe ratio vs. complexity (c) for different ridge penalties")
    print("  - Plot R² vs. complexity (expect negative R² but positive Sharpe!)")
    print("  - Compare performance across training windows T ∈ {12, 60, 120}")
    print("=" * 80)


if __name__ == '__main__':
    main()

