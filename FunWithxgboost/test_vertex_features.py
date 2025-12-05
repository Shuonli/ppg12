#!/usr/bin/env python3
"""
Test script for vertex binning and reweighting features.
"""
import os
import yaml
import pandas as pd
import numpy as np
from main_training import BinnedTrainingPipeline

def test_vertex_configuration():
    """Test that vertex configuration is parsed correctly."""
    print("=== Testing Configuration Parsing ===")
    
    # Test pT binning mode (default)
    print("Testing pT binning mode...")
    pipeline_pt = BinnedTrainingPipeline("config.yaml")
    edges, labels, mode, column = pipeline_pt._get_binning_config()
    print(f"pT mode - Edges: {edges}, Labels: {labels}, Mode: {mode}, Column: {column}")
    assert mode == "pt"
    assert column == "cluster_Et"
    assert len(edges) == 6
    assert len(labels) == 5
    
    # Test vertex binning mode
    print("Testing vertex binning mode...")
    pipeline_vertex = BinnedTrainingPipeline("config_vertex_test.yaml")
    edges, labels, mode, column = pipeline_vertex._get_binning_config()
    print(f"Vertex mode - Edges: {edges}, Labels: {labels}, Mode: {mode}, Column: {column}")
    assert mode == "vertex"
    assert column == "vertexz"
    assert len(edges) == 5
    assert len(labels) == 4
    
    print("✓ Configuration parsing tests passed")

def test_reweighting_options():
    """Test that reweighting options are configured correctly."""
    print("\n=== Testing Reweighting Options ===")
    
    # Load configurations
    with open("config.yaml", "r") as f:
        config_pt = yaml.safe_load(f)
    
    with open("config_vertex_test.yaml", "r") as f:
        config_vertex = yaml.safe_load(f)
    
    # Check that vertex reweighting is disabled in default config
    assert config_pt["reweighting"]["vertex_reweight"] == False
    print("✓ Default config has vertex reweighting disabled")
    
    # Check that vertex reweighting is enabled in test config
    assert config_vertex["reweighting"]["vertex_reweight"] == True
    print("✓ Test config has vertex reweighting enabled")
    
    print("✓ Reweighting options tests passed")

def test_vertex_reweighting_function():
    """Test vertex reweighting function with synthetic data."""
    print("\n=== Testing Vertex Reweighting Function ===")
    
    from reweighting import KinematicReweighter
    
    # Create test configuration
    test_config = {
        'reweighting': {
            'vertex_reweight': True,
            'vertex_reweight_bins': 10,
            'vertex_reweight_max': 5.0,
            'class_reweight': False,
            'et_reweight': False,
            'eta_reweight': False
        },
        'data': {
            'pt_column': 'cluster_Et',
            'label_column': 'label'
        }
    }
    
    # Create synthetic data with non-uniform vertex distribution
    np.random.seed(42)
    n_samples = 1000
    
    # Create biased vertex distribution (more samples at vertex=0)
    vertex_vals = np.concatenate([
        np.random.normal(-5, 1, n_samples//3),  # backward region
        np.random.normal(0, 0.5, n_samples//2),  # central region (more dense)
        np.random.normal(5, 1, n_samples - n_samples//3 - n_samples//2)  # forward region
    ])
    
    df_test = pd.DataFrame({
        'vertexz': vertex_vals,
        'cluster_Et': np.random.uniform(10, 30, n_samples),
        'label': np.random.choice([0, 1], n_samples)
    })
    
    print(f"Created test data with {len(df_test)} samples")
    print(f"Vertex range: {df_test['vertexz'].min():.2f} to {df_test['vertexz'].max():.2f}")
    
    # Test reweighting
    reweighter = KinematicReweighter(test_config)
    
    # Apply class reweighting first to create weight column
    df_reweighted = reweighter.apply_class_reweighting(df_test)
    
    # Apply vertex reweighting
    df_reweighted = reweighter.apply_vertex_reweighting(df_reweighted, label=0)
    df_reweighted = reweighter.apply_vertex_reweighting(df_reweighted, label=1)
    
    print(f"Weight column created with range: {df_reweighted['weight'].min():.3f} to {df_reweighted['weight'].max():.3f}")
    print(f"Mean weight: {df_reweighted['weight'].mean():.3f}")
    
    # Verify weights exist and are reasonable
    assert 'weight' in df_reweighted.columns
    assert df_reweighted['weight'].min() > 0
    assert df_reweighted['weight'].max() <= test_config['reweighting']['vertex_reweight_max']
    
    print("✓ Vertex reweighting function tests passed")

def main():
    """Run all tests."""
    try:
        test_vertex_configuration()
        test_reweighting_options()
        test_vertex_reweighting_function()
        
        print("\n🎉 All tests passed! Vertex binning and reweighting features are working correctly.")
        
        print("\n=== Usage Examples ===")
        print("1. To enable vertex reweighting (flatten vertex distribution):")
        print("   Set reweighting.vertex_reweight: true in config.yaml")
        
        print("\n2. To use vertex binning instead of pT binning:")
        print("   Set binning.mode: vertex in config.yaml")
        print("   Configure binning.vertex_edges and binning.vertex_labels")
        
        print("\n3. Example vertex training:")
        print("   Use config_vertex_test.yaml for vertex binning with vertex reweighting")
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())

