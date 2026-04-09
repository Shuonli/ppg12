---
description: Write new analysis macros (ROOT C++ or Python) following PPG12 repo conventions
---

You write new analysis code for the PPG12 sPHENIX isolated photon analysis. Follow the established conventions in this repo exactly.

## ROOT Macro Conventions

### Config Loading (yaml-cpp)
Every macro that reads config must use:
```cpp
gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
YAML::Node configYaml = YAML::LoadFile(configname);
```

### Function Signature Pattern
```cpp
void MacroName(const std::string &configname = "config_bdt_nom.yaml", bool isMC = true)
```

### Tree and Node Names
- Tree: `slimtree` — always from config: `configYaml["input"]["tree"].as<string>()`
- Cluster node: from config: `configYaml["input"]["cluster_node_name"].as<string>()`
- Current node: `CLUSTERINFO_CEMC`; legacy: `CLUSTERINFO_CEMC_NO_SPLIT`

### Output Filenames
Always include `var_type` from config to prevent overwrites:
```cpp
string var_type = configYaml["output"]["var_type"].as<string>();
string outname = outdir + "/result_" + var_type + ".root";
```

### Parametric Cuts (never hardcode flat thresholds)
```cpp
// BDT threshold
float bdt_min = configYaml["analysis"]["tight_bdt_min_intercept"].as<float>()
              + configYaml["analysis"]["tight_bdt_min_slope"].as<float>() * cluster_ET;

// Isolation
float iso_max = configYaml["analysis"]["reco_iso_max_b"].as<float>()
              + configYaml["analysis"]["reco_iso_max_s"].as<float>() * cluster_ET;
```

### ABCD Regions
- Tight + Isolated = A (signal)
- Tight + Non-isolated = B
- Non-tight + Isolated = C
- Non-tight + Non-isolated = D (control)

### Plotting Macros
```cpp
#include "plotcommon.h"  // provides: NptBins=12, ptRanges[], frames, legend strings, init_plot()
// Also loads BlairUtils.C and sPhenixStyle.C

void MyPlot() {
    init_plot();
    // ptRanges[NptBins+1] = {8,10,12,14,16,18,20,22,24,26,28,30,35}
    // strleg1 = "#bf{#it{sPHENIX}} Internal"
    // strleg2_1 includes luminosity
}
```

### Config Archival
Follow the `SaveYamlToRoot()` pattern from `RecoEffCalculator.C` — store a copy of the config YAML in the output ROOT file for reproducibility.

## Python Conventions

### Config Handling
```python
from ruamel.yaml import YAML  # preserves comments and ordering
yaml = YAML()
with open(config_path) as f:
    config = yaml.load(f)
```

### ROOT I/O
```python
import uproot
tree = uproot.open(path)["slimtree"]
```

### Reproducibility
```python
seed = config['training']['global_seed']  # always 42
np.random.seed(seed)
```

### CLI Pattern
```python
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--config', default='config.yaml')
parser.add_argument('--results', default='results/')
parser.add_argument('--outdir', default='.')
parser.add_argument('--figdir', default='figures/')
```

### Module Structure
Follow existing patterns:
- `data_loader.py` — data loading and preprocessing
- `reweighting.py` — kinematic and class reweighting
- `model_builder.py` — model configuration and training
- `plotting.py` — visualization utilities

### Feature Set
25 features defined in `FunWithxgboost/config.yaml` under `data.features`. When adding/removing features, update BOTH the config AND any C++ code that reads the model (`apply_BDT.C` feature vectors must match).

### Systematic Variants
When adding a new systematic: add to `VARIANTS` dict in `make_bdt_variations.py`, optionally add to `SYST_TYPES`/`SYST_GROUPS`, then re-run the 3-step pipeline (generate configs, run efficiency, aggregate systematics).

## Testing
Follow existing test patterns (`test_*.py` in `FunWithxgboost/`).
