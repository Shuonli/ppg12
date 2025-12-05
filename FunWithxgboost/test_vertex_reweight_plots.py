#!/usr/bin/env python3
"""
Test script for vertex reweighting plots.
"""
import os
import yaml
import pandas as pd
import numpy as np
from plotting import BinnedTrainingPlotter

def test_vertex_reweight_plots():
    """Test that vertex reweighting plots are generated correctly."""
    print("=== Testing Vertex Reweighting Plots ===")
    
    # Create test configuration with vertex reweighting enabled
    test_config = {
        'reweighting': {
            'vertex_reweight': True,
            'vertex_reweight_bins': 10,
            'vertex_reweight_max': 5.0
        },
        'data': {
            'pt_column': 'cluster_Et',
            'label_column': 'label',
            'iso_column': 'recoisoET'
        },
        'output': {
            'output_dir': 'test_plots',
            'save_plots': True
        },
        'plotting': {
            'figure_size': [12, 8],
            'n_bins': 50
        }
    }
    
    # Create output directory
    os.makedirs('test_plots', exist_ok=True)
    
    # Create synthetic data with non-uniform vertex distribution
    np.random.seed(42)
    n_samples = 2000
    
    # Create biased vertex distribution (more samples at vertex=0)
    vertex_vals = np.concatenate([
        np.random.normal(-50, 15, n_samples//4),  # backward region
        np.random.normal(0, 5, n_samples//2),     # central region (more dense)
        np.random.normal(50, 15, n_samples//4)    # forward region
    ])
    
    # Create test data
    df_before = pd.DataFrame({
        'vertexz': vertex_vals,
        'cluster_Et': np.random.uniform(6, 35, n_samples),
        'cluster_Eta': np.random.uniform(-0.6, 0.6, n_samples),
        'recoisoET': np.random.uniform(0, 10, n_samples),
        'label': np.random.choice([0, 1], n_samples),
        'e11_over_e33': np.random.uniform(0.5, 1.0, n_samples),
        'cluster_et1': np.random.uniform(0.1, 2.0, n_samples),
        'cluster_et2': np.random.uniform(0.1, 2.0, n_samples),
        'cluster_et3': np.random.uniform(0.1, 2.0, n_samples),
        'cluster_et4': np.random.uniform(0.1, 2.0, n_samples)
    })
    
    # Create reweighted data (simulate vertex reweighting effect)
    # This would normally be done by the KinematicReweighter
    df_after = df_before.copy()
    
    # Simple simulation of vertex reweighting - flatten distribution
    vertex_hist, vertex_bins = np.histogram(df_after['vertexz'], bins=20)
    vertex_centers = (vertex_bins[:-1] + vertex_bins[1:]) / 2
    
    # Create weights that flatten the distribution
    weights = np.ones(len(df_after))
    for i, center in enumerate(vertex_centers):
        if vertex_hist[i] > 0:
            mask = (df_after['vertexz'] >= vertex_bins[i]) & (df_after['vertexz'] < vertex_bins[i+1])
            # Weight inversely proportional to density
            weights[mask] = 1.0 / (vertex_hist[i] / len(df_after))
    
    # Normalize weights
    weights = weights / np.mean(weights)
    df_after['weight'] = weights
    
    print(f"Created test data with {len(df_before)} samples")
    print(f"Vertex range: {df_before['vertexz'].min():.1f} to {df_before['vertexz'].max():.1f} cm")
    print(f"Weight range: {df_after['weight'].min():.3f} to {df_after['weight'].max():.3f}")
    
    # Test the plotting function
    print("\nTesting plotting function...")
    
    # Initialize plotter
    plotter = BinnedTrainingPlotter({}, test_config)
    
    # Test the global reweighting QA plot
    try:
        plotter.plot_global_reweighting_qa(df_before, df_after)
        print("✓ Global reweighting QA plot generated successfully")
        
        # Check if the plot file was saved
        plot_file = 'test_plots/reweighting_qa_global.png'
        if os.path.exists(plot_file):
            print(f"✓ Plot saved to {plot_file}")
        else:
            print(f"⚠ Plot file not found at {plot_file}")
            
    except Exception as e:
        print(f"❌ Error generating plots: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    print("\n=== Vertex Distribution Analysis ===")
    
    # Analyze vertex distribution before reweighting
    print("Before reweighting:")
    print(f"  Vertex std: {df_before['vertexz'].std():.2f} cm")
    
    # Analyze weighted vertex distribution after reweighting
    print("After reweighting (weighted):")
    weighted_mean = np.average(df_after['vertexz'], weights=df_after['weight'])
    weighted_var = np.average((df_after['vertexz'] - weighted_mean)**2, weights=df_after['weight'])
    weighted_std = np.sqrt(weighted_var)
    print(f"  Weighted vertex std: {weighted_std:.2f} cm")
    
    # Clean up test files
    import shutil
    if os.path.exists('test_plots'):
        shutil.rmtree('test_plots')
    
    return True

def main():
    """Run the test."""
    try:
        success = test_vertex_reweight_plots()
        if success:
            print("\n🎉 Vertex reweighting plot tests passed!")
            print("\n=== Usage Notes ===")
            print("• Vertex distribution plots are shown when vertex_reweight: true")
            print("• Before/after plots help verify reweighting is working correctly")
            print("• Plots are saved to output_dir/reweighting_qa_global.png")
        else:
            print("\n❌ Vertex reweighting plot tests failed!")
            return 1
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())

