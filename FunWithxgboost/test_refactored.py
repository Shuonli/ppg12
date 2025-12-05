#!/usr/bin/env python3
"""
Test script for the refactored binned training code.
"""
import os
import sys
import yaml

def test_imports():
    """Test that all modules can be imported."""
    try:
        from data_loader import DataLoader
        from reweighting import KinematicReweighter
        from model_builder import ModelBuilder
        from plotting import BinnedTrainingPlotter
        from main_training import BinnedTrainingPipeline
        print("✓ All modules imported successfully")
        return True
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False

def test_config_loading():
    """Test configuration file loading."""
    try:
        with open('config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        
        required_keys = ['training', 'model', 'binning', 'data', 'reweighting', 'output']
        for key in required_keys:
            if key not in config:
                print(f"✗ Missing config key: {key}")
                return False
        
        print("✓ Configuration loaded successfully")
        return True
    except Exception as e:
        print(f"✗ Config loading failed: {e}")
        return False

def test_class_instantiation():
    """Test that classes can be instantiated."""
    try:
        with open('config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        
        from data_loader import DataLoader
        from reweighting import KinematicReweighter
        from model_builder import ModelBuilder
        
        # Test instantiation
        loader = DataLoader(config)
        reweighter = KinematicReweighter(config)
        model_builder = ModelBuilder(config)
        
        print("✓ All classes instantiated successfully")
        return True
    except Exception as e:
        print(f"✗ Class instantiation failed: {e}")
        return False

def main():
    """Run all tests."""
    print("Testing Refactored Binned Training Code")
    print("=" * 50)
    
    tests = [
        ("Module Imports", test_imports),
        ("Configuration Loading", test_config_loading),
        ("Class Instantiation", test_class_instantiation),
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        if test_func():
            passed += 1
        else:
            print(f"  {test_name} failed")
    
    print("\n" + "=" * 50)
    print(f"Tests passed: {passed}/{total}")
    
    if passed == total:
        print("🎉 All tests passed! The refactored code is ready to use.")
        return 0
    else:
        print("❌ Some tests failed. Please check the errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
