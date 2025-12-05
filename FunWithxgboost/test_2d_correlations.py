#!/usr/bin/env python3
"""
Test script for the new 2D correlation plots with profile overlays.
This script tests the correlation analysis functionality without running the full pipeline.
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# Add current directory to path
sys.path.append('.')

def test_2d_correlation_plot():
    """Test the 2D correlation plotting functionality with synthetic data."""
    print("=== Testing 2D Correlation Plots ===")
    
    # Create synthetic data similar to the physics data
    np.random.seed(42)
    n_samples = 10000
    
    # Create correlated features
    feature1 = np.random.normal(0, 1, n_samples)
    feature2 = np.random.normal(0, 1, n_samples)
    feature3 = np.random.normal(0, 1, n_samples)
    
    # Create isoET with correlations
    isoET = 0.3 * feature1 + 0.2 * feature2 + 0.1 * feature3 + np.random.normal(0, 0.5, n_samples)
    
    # Create DataFrame
    df = pd.DataFrame({
        'feature1': feature1,
        'feature2': feature2, 
        'feature3': feature3,
        'recoisoET': isoET,
        'label': np.random.choice([0, 1], n_samples, p=[0.7, 0.3])  # 70% background, 30% signal
    })
    
    # Filter to background only
    df_background = df[df['label'] == 0]
    print(f"Created synthetic dataset with {len(df_background)} background samples")
    
    # Test the plotting functions
    from main_training import BinnedTrainingPipeline
    
    # Create a minimal config for testing
    config = {
        'data': {
            'iso_column': 'recoisoET',
            'label_column': 'label'
        },
        'output': {
            'save_plots': True,
            'output_dir': 'test_output'
        }
    }
    
    # Create pipeline instance
    pipeline = BinnedTrainingPipeline.__new__(BinnedTrainingPipeline)
    pipeline.config = config
    pipeline.output_dir = 'test_output'
    
    # Create output directory
    os.makedirs('test_output', exist_ok=True)
    
    # Test correlation analysis
    features = ['feature1', 'feature2', 'feature3']
    iso_col = 'recoisoET'
    
    # Calculate correlations
    correlations = {}
    iso_values = df_background[iso_col].values
    
    print(f"\nCorrelations between {iso_col} and features:")
    print("-" * 50)
    
    for feature in features:
        feature_values = df_background[feature].values
        valid_mask = ~(np.isnan(iso_values) | np.isnan(feature_values))
        
        if np.sum(valid_mask) > 10:
            iso_clean = iso_values[valid_mask]
            feature_clean = feature_values[valid_mask]
            r, p = stats.pearsonr(iso_clean, feature_clean)
            correlations[feature] = {'r': r, 'p': p, 'n': len(iso_clean)}
            print(f"{feature:<15} r={r:.3f}, p={p:.2e}, n={len(iso_clean)}")
    
    # Test the plotting functions
    print("\n--- Testing 2D correlation plots ---")
    
    # Test single 2D plot
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    pipeline._plot_2d_correlation_with_profile(ax, df_background, 'feature1', iso_col, correlations['feature1']['r'])
    plt.tight_layout()
    plt.savefig('test_output/test_2d_correlation.png', dpi=150, bbox_inches='tight')
    print("✓ Single 2D correlation plot created: test_output/test_2d_correlation.png")
    plt.close()
    
    # Test 2D correlation grid
    pipeline._create_2d_correlation_grid(df_background, features, iso_col, correlations)
    print("✓ 2D correlation grid created: test_output/2d_correlations_background.png")
    
    # Test full correlation analysis
    print("\n--- Testing full correlation analysis ---")
    pipeline.analyze_feature_correlations(df)
    print("✓ Full correlation analysis completed")
    
    print("\n=== Test completed successfully ===")
    print("Generated files:")
    print("- test_output/test_2d_correlation.png")
    print("- test_output/2d_correlations_background.png") 
    print("- test_output/feature_correlations_background.png")

if __name__ == "__main__":
    test_2d_correlation_plot()















