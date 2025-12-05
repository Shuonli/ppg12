#!/usr/bin/env python3
"""
Test script for the y-axis zoom functionality in 2D correlation plots.
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Add current directory to path
sys.path.append('.')

def test_y_axis_zoom():
    """Test the y-axis zoom functionality with synthetic data."""
    print("=== Testing Y-Axis Zoom for 2D Correlation Plots ===")
    
    # Create synthetic data with a wide range of isoET values
    np.random.seed(42)
    n_samples = 10000
    
    # Create features
    feature1 = np.random.normal(0, 1, n_samples)
    feature2 = np.random.normal(0, 1, n_samples)
    
    # Create isoET with some outliers to test zooming
    isoET = 0.3 * feature1 + 0.2 * feature2 + np.random.normal(0, 0.5, n_samples)
    
    # Add some outliers to test zoom functionality
    outlier_indices = np.random.choice(n_samples, size=100, replace=False)
    isoET[outlier_indices] = np.random.uniform(10, 50, 100)  # High outliers
    
    # Create DataFrame
    df = pd.DataFrame({
        'feature1': feature1,
        'feature2': feature2,
        'recoisoET': isoET,
        'label': np.random.choice([0, 1], n_samples, p=[0.7, 0.3])
    })
    
    # Filter to background only
    df_background = df[df['label'] == 0]
    print(f"Created synthetic dataset with {len(df_background)} background samples")
    print(f"isoET range: [{df_background['recoisoET'].min():.2f}, {df_background['recoisoET'].max():.2f}]")
    
    # Test both zoomed and unzoomed versions
    from main_training import BinnedTrainingPipeline
    
    # Create a minimal config for testing
    config = {
        'data': {
            'iso_column': 'recoisoET',
            'label_column': 'label'
        },
        'output': {
            'save_plots': True,
            'output_dir': 'test_y_zoom_output'
        },
        'plotting': {
            'zoom_y_axis': True
        }
    }
    
    # Create pipeline instance
    pipeline = BinnedTrainingPipeline.__new__(BinnedTrainingPipeline)
    pipeline.config = config
    pipeline.output_dir = 'test_y_zoom_output'
    
    # Create output directory
    os.makedirs('test_y_zoom_output', exist_ok=True)
    
    # Calculate correlations
    from scipy import stats
    correlations = {}
    iso_values = df_background['recoisoET'].values
    
    for feature in ['feature1', 'feature2']:
        feature_values = df_background[feature].values
        valid_mask = ~(np.isnan(iso_values) | np.isnan(feature_values))
        iso_clean = iso_values[valid_mask]
        feature_clean = feature_values[valid_mask]
        r, p = stats.pearsonr(iso_clean, feature_clean)
        correlations[feature] = {'r': r, 'p': p, 'n': len(iso_clean)}
        print(f"{feature}: r={r:.3f}")
    
    # Test zoomed version
    print("\n--- Testing with Y-axis zoom enabled ---")
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    pipeline._plot_2d_correlation_with_profile(ax, df_background, 'feature1', 'recoisoET', 
                                             correlations['feature1']['r'], zoom_y=True)
    plt.tight_layout()
    plt.savefig('test_y_zoom_output/2d_correlation_zoomed.png', dpi=150, bbox_inches='tight')
    print("✓ Zoomed 2D correlation plot created: test_y_zoom_output/2d_correlation_zoomed.png")
    plt.close()
    
    # Test unzoomed version
    print("\n--- Testing with Y-axis zoom disabled ---")
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    pipeline._plot_2d_correlation_with_profile(ax, df_background, 'feature1', 'recoisoET', 
                                             correlations['feature1']['r'], zoom_y=False)
    plt.tight_layout()
    plt.savefig('test_y_zoom_output/2d_correlation_unzoomed.png', dpi=150, bbox_inches='tight')
    print("✓ Unzoomed 2D correlation plot created: test_y_zoom_output/2d_correlation_unzoomed.png")
    plt.close()
    
    print("\n=== Test completed successfully ===")
    print("Compare the two plots to see the difference:")
    print("- 2d_correlation_zoomed.png: Focuses on 0.5-99.5% range of isoET")
    print("- 2d_correlation_unzoomed.png: Shows full range including outliers")

if __name__ == "__main__":
    test_y_axis_zoom()















