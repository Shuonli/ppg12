# Vertex Binning and Reweighting Features

This document describes the new vertex-based features added to the binned training pipeline.

## Overview

Two new features have been added:

1. **Vertex Distribution Reweighting**: Flatten the vertex z distribution to ensure uniform representation
2. **Vertex Binning**: Train models in vertex z bins instead of pT bins

## Feature 1: Vertex Distribution Reweighting

### Purpose
- Flatten the vertex z distribution to reduce bias from non-uniform vertex distribution
- Ensures equal representation across different vertex z positions
- Applied globally before any binning (pT or vertex)

### Configuration

```yaml
reweighting:
  vertex_reweight: true          # Enable vertex distribution flattening
  vertex_reweight_bins: 20       # Number of bins for reweighting
  vertex_reweight_max: 50.0      # Maximum weight cap to prevent outliers
```

### How it Works
- Creates a histogram of vertex z distribution per class (signal/background)
- Computes weights as inverse of the probability density
- Applies weight cap to prevent extreme outliers
- Normalizes weights to have mean = 1.0

### Quality Assurance Plots
When vertex reweighting is enabled, the global reweighting QA plot automatically includes:
- **Vertex Z Before Reweighting**: Shows the original vertex distribution
- **Vertex Z After Reweighting**: Shows the flattened vertex distribution (weighted)
- Saved to `output_dir/reweighting_qa_global.png`

## Feature 2: Vertex Binning

### Purpose
- Train separate models for different vertex z regions
- Study detector acceptance and efficiency as a function of vertex position
- Alternative to pT binning when vertex position is the primary variable of interest

### Configuration

```yaml
binning:
  mode: vertex                   # Switch from "pt" to "vertex" binning
  vertex_edges: [-10.0, -5.0, 0.0, 5.0, 10.0]  # Bin edges in cm
  vertex_labels: ['backward', 'central', 'forward', 'far_forward']  # Bin names

training:
  train_single_model: false      # Train separate models per vertex bin
```

### Vertex Bin Definitions
- **backward**: -10 to -5 cm (upstream)
- **central**: -5 to 0 cm (central detector)
- **forward**: 0 to 5 cm (downstream)
- **far_forward**: 5 to 10 cm (far downstream)

## Usage Examples

### Example 1: Enable vertex reweighting only
```yaml
reweighting:
  vertex_reweight: true

binning:
  mode: pt  # Keep pT binning
```

### Example 2: Switch to vertex binning without reweighting
```yaml
reweighting:
  vertex_reweight: false

binning:
  mode: vertex
  vertex_edges: [-8.0, -4.0, 0.0, 4.0, 8.0]
  vertex_labels: ['backward', 'central', 'forward']

training:
  train_single_model: false
```

### Example 3: Both vertex reweighting and vertex binning
```yaml
reweighting:
  vertex_reweight: true
  vertex_reweight_bins: 15
  vertex_reweight_max: 30.0

binning:
  mode: vertex
  vertex_edges: [-10.0, -5.0, 0.0, 5.0, 10.0]
  vertex_labels: ['backward', 'central', 'forward', 'far_forward']

training:
  train_single_model: false
```

## Implementation Details

### Processing Order
1. **Data Loading**: Load all data using single_file_set or per_bin_file_sets
2. **Global Reweighting**: Apply all reweighting (class, eta, pT, vertex) to full dataset
3. **Binning**: Split data into bins based on binning mode (pT or vertex)
4. **Training**: Train models per bin or single model for all bins

### Key Implementation Points
- Reweighting always happens before binning
- Vertex binning uses the 'vertexz' column
- Backward compatibility maintained with existing pT binning
- Configuration validation ensures consistent settings

### Files Modified
- `reweighting.py`: Added `apply_vertex_reweighting()` method
- `main_training.py`: Added vertex binning support
- `config.yaml`: Added vertex configuration options
- `data_loader.py`: Updated to support both binning modes

### Output Differences
When using vertex binning:
- Model files named with vertex bin labels (e.g., `backward`, `central`)
- Metrics include `bin_min`/`bin_max` instead of `et_min`/`et_max`
- Plots labeled with vertex positions in cm instead of pT in GeV

### ROOT File Naming
Configure the prefix for ROOT/TMVA model files:
```yaml
output:
  root_file_prefix: "model"  # Creates: model_single_tmva.root, model_backward_tmva.root, etc.
```

Examples:
- `root_file_prefix: "vertex_model"` → `vertex_model_backward_tmva.root`
- `root_file_prefix: "bdt"` → `bdt_single_tmva.root`

## Testing

Use the provided test configuration and script:

```bash
# Test with vertex binning configuration
python main_training.py config_vertex_test.yaml

# Run feature tests
python test_vertex_features.py
```

## Performance Considerations

- Vertex reweighting adds minimal computational overhead
- Vertex binning training time depends on data distribution across bins
- Memory usage similar to pT binning
- Single model training recommended for initial testing

## Tips for Usage

1. **Start with reweighting**: Enable vertex reweighting first to understand its impact
2. **Check data distribution**: Ensure sufficient statistics in each vertex bin
3. **Monitor weights**: Check vertex reweighting weight distributions
4. **Compare results**: Compare vertex vs pT binning performance
5. **Combine wisely**: Use both features when vertex position significantly affects performance
