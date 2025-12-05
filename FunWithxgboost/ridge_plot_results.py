#!/usr/bin/env python
"""
Plot results from ridge_reproduce.py to match Kelly-Malamud-Zhou (2024) figures.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Iterable, Tuple, Optional


DEFAULT_FORMATS: Tuple[str, ...] = ('png', 'pdf')
BASELINE_STYLES = {
    'Volatility-Timed Momentum': {'color': '#d62728', 'linestyle': '--'},
    'Gaussian Kernel Regression': {'color': '#1f77b4', 'linestyle': '--'},
    'MLP (Raw)': {'color': '#ff7f0e', 'linestyle': '--'},
    'MLP (RFF)': {'color': '#9467bd', 'linestyle': '--'},
    'SIREN (Raw)': {'color': '#2ca02c', 'linestyle': '--'},
    'SIREN (RFF)': {'color': '#8c564b', 'linestyle': '--'}
}


def _save_figure(fig: plt.Figure, output_dir: Path, filename: str,
                 formats: Iterable[str]):
    """Save figure in requested formats."""
    for fmt in formats:
        output_file = output_dir / f"{filename}.{fmt}"
        fig.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_file}")


def _draw_baseline_hlines(ax, metric: str, baseline_df: Optional[pd.DataFrame]):
    if baseline_df is None or baseline_df.empty or metric not in baseline_df.columns:
        return
    labels_added = set()
    for _, row in baseline_df.iterrows():
        value = row.get(metric)
        if pd.isna(value):
            continue
        method = row.get('method', 'Baseline')
        style = BASELINE_STYLES.get(method, {'color': 'black', 'linestyle': '--'})
        label = f"{method} ({value:.3f})"
        if label in labels_added:
            label = None
        ax.axhline(y=value, linewidth=1.5, alpha=0.8, label=label, **style)
        if label:
            labels_added.add(label)


def _draw_baseline_vlines(ax, metric: str, baseline_df: Optional[pd.DataFrame]):
    if baseline_df is None or baseline_df.empty or metric not in baseline_df.columns:
        return
    labels_added = set()
    for _, row in baseline_df.iterrows():
        value = row.get(metric)
        if pd.isna(value):
            continue
        method = row.get('method', 'Baseline')
        style = BASELINE_STYLES.get(method, {'color': 'black', 'linestyle': '--'})
        label = f"{method} ({value:.3f})"
        if label in labels_added:
            label = None
        ax.axvline(x=value, linewidth=1.5, alpha=0.8, label=label, **style)
        if label:
            labels_added.add(label)

def plot_voc_curves(results_df: pd.DataFrame, output_dir: Path = Path('.'),
                    formats: Iterable[str] = DEFAULT_FORMATS,
                    baseline_df: Optional[pd.DataFrame] = None):
    """
    Plot "Virtue of Complexity" curves matching paper Figures 7-8.
    
    Creates plots of:
    - Out-of-sample R²
    - Mean return
    - Volatility  
    - Sharpe ratio
    - Information ratio
    
    All vs. model complexity (c = P/T) for different ridge penalties (z).
    """
    output_dir.mkdir(exist_ok=True)
    
    # Get unique ridge alphas
    alphas = sorted(results_df['ridge_alpha'].unique())
    
    # Color map for different ridge penalties
    colors = plt.cm.viridis(np.linspace(0, 1, len(alphas)))
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Virtue of Complexity: Out-of-Sample Performance vs. Model Complexity\n' +
                 'Kelly, Malamud, Zhou (2024) Reproduction', 
                 fontsize=14, fontweight='bold')
    
    metrics = [
        ('oos_r2', 'Out-of-Sample R²', axes[0, 0]),
        ('mean_return', 'Mean Return', axes[0, 1]),
        ('volatility', 'Volatility', axes[0, 2]),
        ('sharpe_ratio', 'Sharpe Ratio', axes[1, 0]),
        ('information_ratio', 'Information Ratio', axes[1, 1]),
        ('beta_norm', 'Beta Norm', axes[1, 2])
    ]
    
    for metric, title, ax in metrics:
        for alpha, color in zip(alphas, colors):
            subset = results_df[results_df['ridge_alpha'] == alpha].sort_values('complexity')
            
            label = f'z={alpha:.0e}'
            
            # Plot all data
            ax.plot(subset['complexity'], subset[metric], 
                   marker='o', markersize=3, linewidth=1.5, 
                   color=color, label=label, alpha=0.8)
        
        ax.set_xlabel('Complexity (c = P/T)', fontsize=11)
        ax.set_ylabel(title, fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        _draw_baseline_hlines(ax, metric, baseline_df)
        ax.legend(fontsize=8, loc='best', ncol=2)
        
        # Use log scale for x-axis if max complexity > 10
        max_c = results_df['complexity'].max()
        if max_c > 10:
            ax.set_xscale('log')
            ax.set_xlabel('Complexity (c = P/T, log scale)', fontsize=11)
        
        # Add vertical line at interpolation boundary c=1
        ax.axvline(x=1.0, color='red', linestyle='--', linewidth=1.5, 
                  alpha=0.5, label='c=1')
    
    plt.tight_layout()
    
    _save_figure(fig, output_dir, 'voc_curves', formats)
    plt.close(fig)


def plot_sharpe_heatmap(results_df: pd.DataFrame, output_dir: Path = Path('.'),
                        formats: Iterable[str] = DEFAULT_FORMATS,
                        baseline_df: Optional[pd.DataFrame] = None):
    """
    Create heatmap of Sharpe ratio vs. (complexity, ridge_alpha).
    """
    output_dir.mkdir(exist_ok=True)
    
    # Pivot to create heatmap data
    pivot_data = results_df.pivot_table(
        values='sharpe_ratio',
        index='log10_alpha',
        columns='complexity',
        aggfunc='mean'
    )
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    im = ax.imshow(pivot_data.values, aspect='auto', cmap='RdYlGn', 
                   interpolation='nearest', origin='lower')
    
    # Set ticks
    ax.set_xticks(range(len(pivot_data.columns)))
    ax.set_xticklabels([f'{c:.1f}' for c in pivot_data.columns], rotation=45)
    ax.set_yticks(range(len(pivot_data.index)))
    ax.set_yticklabels([f'{a:.0f}' for a in pivot_data.index])
    
    ax.set_xlabel('Complexity (c = P/T)', fontsize=12)
    ax.set_ylabel('log₁₀(Ridge Penalty z)', fontsize=12)
    ax.set_title('Sharpe Ratio Heatmap\nKelly, Malamud, Zhou (2024)', 
                fontsize=14, fontweight='bold')
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Sharpe Ratio', fontsize=11)
    
    plt.tight_layout()
    
    _save_figure(fig, output_dir, 'sharpe_heatmap', formats)
    plt.close(fig)


def plot_r2_vs_sharpe(results_df: pd.DataFrame, output_dir: Path = Path('.'),
                      formats: Iterable[str] = DEFAULT_FORMATS,
                      baseline_df: Optional[pd.DataFrame] = None):
    """
    Scatter plot showing R² vs. Sharpe ratio (key paper insight).
    
    Shows that negative R² can still yield positive Sharpe ratios.
    """
    output_dir.mkdir(exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Color by complexity
    scatter = ax.scatter(results_df['oos_r2'], results_df['sharpe_ratio'],
                        c=np.log10(results_df['complexity'] + 0.1),
                        cmap='plasma', alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
    
    ax.axhline(y=0, color='red', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(x=0, color='red', linestyle='--', linewidth=1, alpha=0.5)
    
    ax.set_xlabel('Out-of-Sample R²', fontsize=12)
    ax.set_ylabel('Sharpe Ratio', fontsize=12)
    ax.set_title('R² vs. Sharpe Ratio: Negative R² Can Yield Positive Sharpe!\n' +
                 'Kelly, Malamud, Zhou (2024)', 
                 fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    _draw_baseline_vlines(ax, 'oos_r2', baseline_df)
    _draw_baseline_hlines(ax, 'sharpe_ratio', baseline_df)
    
    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('log₁₀(Complexity)', fontsize=11)
    
    # Annotate key insight
    ax.text(0.05, 0.95, 
           'Key Insight: Many models with R² < 0\nstill achieve Sharpe > 0',
           transform=ax.transAxes, fontsize=11,
           verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    _save_figure(fig, output_dir, 'r2_vs_sharpe', formats)
    plt.close(fig)


def plot_complexity_comparison(results_df: pd.DataFrame, output_dir: Path = Path('.'),
                               formats: Iterable[str] = DEFAULT_FORMATS,
                               baseline_df: Optional[pd.DataFrame] = None):
    """
    Bar chart comparing low vs. high complexity models.
    """
    output_dir.mkdir(exist_ok=True)
    
    # Define complexity bins (adapt to max complexity)
    max_c = results_df['complexity'].max()
    if max_c > 100:
        bins = [0, 0.5, 1.5, 10, 100, max_c + 1]
        labels = ['Low (c<0.5)', 'Medium (0.5≤c<1.5)', 'High (1.5≤c<10)', 
                  'Very High (10≤c<100)', 'Extreme (c≥100)']
    else:
        bins = [0, 0.5, 1.5, max_c + 1]
        labels = ['Low (c<0.5)', 'Medium (0.5≤c<1.5)', 'High (c≥1.5)']
    
    results_df['complexity_category'] = pd.cut(
        results_df['complexity'],
        bins=bins,
        labels=labels
    )
    
    # Aggregate by category
    summary = results_df.groupby('complexity_category').agg({
        'sharpe_ratio': 'mean',
        'information_ratio': 'mean',
        'oos_r2': 'mean'
    })
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Performance by Model Complexity\nKelly, Malamud, Zhou (2024)', 
                 fontsize=14, fontweight='bold')
    
    metrics = [
        ('sharpe_ratio', 'Sharpe Ratio', axes[0]),
        ('information_ratio', 'Information Ratio', axes[1]),
        ('oos_r2', 'Out-of-Sample R²', axes[2])
    ]
    
    x = np.arange(len(summary.index))
    width = 0.6
    
    for metric, title, ax in metrics:
        # Adapt colors to number of categories
        n_cats = len(summary.index)
        colors_list = plt.cm.viridis(np.linspace(0.2, 0.9, n_cats))
        
        bars = ax.bar(x, summary[metric], width, 
                     color=colors_list, alpha=0.8)
        
        ax.set_xlabel('Complexity Category', fontsize=11)
        ax.set_ylabel(title, fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(summary.index, rotation=15, ha='right')
        ax.grid(True, alpha=0.3, axis='y')
        ax.axhline(y=0, color='red', linestyle='--', linewidth=1, alpha=0.5)
        _draw_baseline_hlines(ax, metric, baseline_df)
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.3f}',
                   ha='center', va='bottom' if height >= 0 else 'top',
                   fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    
    _save_figure(fig, output_dir, 'complexity_comparison', formats)
    plt.close(fig)


def main():
    """Generate all plots."""
    
    print("=" * 80)
    print("Plotting Kelly-Malamud-Zhou (2024) Results")
    print("=" * 80)
    print()
    
    # Load results
    results_path = Path('ridge_results.csv')
    if not results_path.exists():
        print(f"Error: {results_path} not found. Run ridge_reproduce.py first.")
        return
    
    results_df = pd.read_csv(results_path)
    print(f"Loaded {len(results_df)} result rows")
    print()
    
    baseline_path = Path('ridge_baselines.csv')
    if baseline_path.exists():
        baseline_df = pd.read_csv(baseline_path)
        print(f"Loaded baseline metrics from {baseline_path}")
    else:
        baseline_df = None
        print("Warning: Baseline metrics file ridge_baselines.csv not found. "
              "Dashed reference lines will be skipped.")
    print()
    
    # Create output directory
    output_dir = Path('ridge_plots')
    output_dir.mkdir(exist_ok=True)
    
    # Generate plots
    print("Generating plots...")
    print("-" * 80)
    
    plot_voc_curves(results_df, output_dir, baseline_df=baseline_df)
    plot_sharpe_heatmap(results_df, output_dir, baseline_df=baseline_df)
    plot_r2_vs_sharpe(results_df, output_dir, baseline_df=baseline_df)
    plot_complexity_comparison(results_df, output_dir, baseline_df=baseline_df)
    
    print()
    print("=" * 80)
    print(f"All plots saved to: {output_dir}/")
    print("=" * 80)


if __name__ == '__main__':
    main()

