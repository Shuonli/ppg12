# PPG12 sPHENIX Isolated Photon Analysis

This repository contains the complete analysis pipeline for measuring isolated photon cross-sections in proton-proton collisions using the sPHENIX detector at RHIC. The analysis uses machine learning (XGBoost) for photon identification and produces publication-quality results compared with JETPHOX NLO theory predictions.

## Table of Contents
- [Overview](#overview)
- [Quick Start](#quick-start)
- [Directory Structure](#directory-structure)
- [Analysis Pipeline](#analysis-pipeline)
- [Configuration](#configuration)
- [Key Workflows](#key-workflows)
- [Systematic Uncertainties](#systematic-uncertainties)
- [Output Files](#output-files)
- [Documentation](#documentation)

## Overview

**Goal**: Measure the differential cross-section of isolated photons as a function of transverse momentum (pT) and compare with theoretical predictions.

**Key Steps**:
1. Process raw sPHENIX data into analysis trees
2. Train machine learning models to identify photons vs jet backgrounds
3. Calculate selection efficiencies from simulation
4. Extract photon yields from data using sideband subtraction
5. Apply unfolding corrections and produce final cross-sections

**Physics Motivation**: Isolated photons are a clean probe of QCD processes and can constrain parton distribution functions. They provide direct information about hard scattering processes in proton-proton collisions.

## Quick Start

### Prerequisites
- ROOT 6.x with RooUnfold
- Python 3.x with: numpy, pandas, matplotlib, scikit-learn, xgboost, PyYAML, uproot
- Access to sPHENIX computing environment
- Condor batch system (for parallel processing)

### Basic Workflow

```bash
# 1. Create analysis trees from raw data
cd anatreemaker/macro_maketree
# [Configure and run Fun4All macros]

# 2. Train BDT models
cd ../../FunWithxgboost
root -l -b -q 'BDTinput.C'  # Prepare training data
python main_training.py --config config_nom.yaml  # Train models

# 3. Apply BDT to data
root -l -b -q 'apply_BDT.C'

# 4. Calculate efficiencies and yields
cd ../efficiencytool
./oneforall.sh config_nom.yaml

# 5. Generate final plots
cd ../plotting
root -l -b -q 'plot_final.C'
```

## Directory Structure

### Core Analysis Components

#### `anatreemaker/`
**Purpose**: Converts raw sPHENIX DST files into analysis-ready ROOT trees.

- `source/CaloAna24.cc` - Main analysis module
  - Processes calorimeter clusters (CEMC, IHCAL, OHCAL)
  - Extracts shower shape variables
  - Stores isolation, truth-matching, MBD timing
  - Outputs "slimtree" with ~100+ features per cluster

**Output**: ROOT files with TTree containing cluster properties for ML training and physics analysis.

#### `FunWithxgboost/`
**Purpose**: Machine learning pipeline for photon identification.

**Key Files**:
- `BDTinput.C` - Prepares training samples from anatrees
- `main_training.py` - Orchestrates BDT training pipeline
- `data_loader.py` - Loads data with background flattening
- `reweighting.py` - Kinematic reweighting (eta, pT, vertex)
- `plotting.py` - ROC curves, feature importance, QA plots
- `apply_BDT.C` - Applies trained models to data

**Training Features**:
- Shower shapes: weta, wphi, chi2, e11/e33, e32/e35
- Energy ratios: et1-4
- Isolation: cone energies at R=0.2, 0.3, 0.4

**Training Modes**:
- Single model for all kinematics
- pT-binned: separate models for 5-10, 10-20, 20+ GeV
- Vertex-binned: backward/central/forward regions

**Configuration**: `config_nom.yaml` (see Configuration section)

#### `efficiencytool/`
**Purpose**: Efficiency calculations and photon yield extraction.

**Key Macros**:
- `RecoEffCalculator.C` - Calculates ID efficiencies from simulation
  - Builds truth-reco response matrices for unfolding
  - Handles cross-section weighting for multi-sample MC
  - Categories: tight_iso, tight_noniso, nontight_iso, nontight_noniso

- `CalculatePhotonYield.C` - Extracts yields from data
  - ABCD sideband background subtraction
  - Bayesian unfolding (RooUnfold)
  - Efficiency corrections
  - Final cross-section calculation

- `MergeSim.C` - Combines photon/jet samples with proper weights
- `Cluster_rbr.C` - Run-by-run QA

**Automation Scripts**:
- `oneforall.sh` - Runs full efficiency chain for a given config
- `oneforall.sub` - Condor submission script

#### `plotting/`
**Purpose**: Publication-quality plots and final results.

**Key Macros**:
- `plot_final.C` - Main results: data vs JETPHOX NLO theory
- `plot_efficiency.C` - Efficiency vs pT
- `plot_purity_*.C` - Purity estimation (signal fraction)
- `plot_sideband_*.C` - ABCD method validation
- `CONF_plots.C` - Conference presentation figures

**Automation Scripts**:
- `make_all_selection.sh` - Generates all selection plots
- `make_syst.sh` - Systematic uncertainty plots

#### `NLO/`
**Purpose**: Process JETPHOX NLO theory predictions.

- `MakeJetPHOXhisto.C` - Converts theory output to ROOT histograms

#### `toymcunfold/`
**Purpose**: Unfolding validation studies.

- `macro/jet_pt_unfolding.C` - Tests Bayesian unfolding closure

### Supporting Components

- `dataana/` - Data processing utilities
- `eventskimmer/` - Event pre-selection
- `saturationcheck/` - Calorimeter saturation studies
- `showershapecheck/` - Shower shape validation
- `sideband/` - Sideband method development
- `photontruthisopythia/` - Truth isolation studies
- `simcrosssection/` - Cross-section calculations
- `PPG12-analysis-note/` - LaTeX documentation

## Analysis Pipeline

### Complete Data Flow

```
┌─────────────────────────────────────────────────────────────────┐
│ 1. Raw Data (sPHENIX DST files)                                 │
└──────────────────────┬──────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ 2. anatreemaker/CaloAna24.cc                                     │
│    → Produces "slimtree" ROOT files with cluster properties     │
└──────────────────────┬──────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ 3. FunWithxgboost/BDTinput.C                                     │
│    → Prepares training samples (signal: photons, bkg: jets)     │
└──────────────────────┬──────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ 4. FunWithxgboost/main_training.py                               │
│    → Trains XGBoost models with reweighting                     │
│    → Outputs: bdt_model*.json, bdt_model*.root                  │
└──────────────────────┬──────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ 5. FunWithxgboost/apply_BDT.C                                    │
│    → Scores clusters, adds BDT branches to trees                │
└──────────────────────┬──────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ 6. efficiencytool/MergeSim.C                                     │
│    → Combines photon5/10/20 + jet10/15/20/30/50 samples         │
└──────────────────────┬──────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ 7. efficiencytool/RecoEffCalculator.C                            │
│    → Calculates efficiencies from simulation                    │
│    → Builds response matrices for unfolding                     │
└──────────────────────┬──────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ 8. efficiencytool/CalculatePhotonYield.C                         │
│    → ABCD sideband background subtraction (data)                │
│    → Bayesian unfolding with response matrix (data + MC)        │
│    → Efficiency corrections                                     │
│    → Final cross-section calculation                            │
└──────────────────────┬──────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ 9. plotting/plot_final.C                                         │
│    → Data vs theory comparison with systematics                 │
│    → Publication-ready figures                                  │
└─────────────────────────────────────────────────────────────────┘
```

### Monte Carlo Samples

**Signal Samples** (Pythia8 isolated photons):
- `photon5/` - Generated pT: 5-10 GeV
- `photon10/` - Generated pT: 10-20 GeV
- `photon20/` - Generated pT: 20+ GeV

**Background Samples** (Pythia8 dijets):
- `jet10/`, `jet15/`, `jet20/`, `jet30/`, `jet50/` - Various pT thresholds

Samples are combined with cross-section weights to avoid MC statistical limitations.

## Configuration

### Main Configuration Files

#### 1. BDT Training: `FunWithxgboost/config_nom.yaml`

```yaml
# Data loading
data_loading:
  mode: "single"  # or "binned"
  signal_file: "path/to/photon_sample.root"
  background_file: "path/to/jet_sample.root"

# Binning
binning:
  mode: "single"  # "single", "pt_binned", or "vertex_binned"
  pt_bins: [5, 10, 20, 50]
  vertex_bins: [-10, -3, 3, 10]

# Reweighting
reweighting:
  class_weight: true
  eta_flatten: true
  et_flatten: true
  vertex_flatten: false

# Model parameters
model:
  type: "xgboost"  # or "histgradient"
  max_depth: 6
  n_estimators: 100
  learning_rate: 0.1
```

#### 2. Efficiency Calculation: `efficiencytool/config_nom.yaml`

```yaml
# Cut definitions
cuts:
  eta_range: [-1.1, 1.1]
  vertex_range: [-10, 10]
  iso_threshold: 0.1
  bdt_threshold: 0.5

# Shower shape criteria
tight_cuts:
  chi2: 4.0
  weta: 0.015
  wphi: 0.03
```

### Systematic Variation Configs

Located in `efficiencytool/config_*.yaml`:
- `config_iso05.yaml`, `config_iso20.yaml` - Isolation ±50%
- `config_escale.yaml` - Energy scale ±1%
- `config_eres.yaml` - Energy resolution
- `config_ntf1.yaml`, `config_ntf3.yaml` - Template fit variations
- `config_mbdeffup.yaml`, `config_mbdeffdown.yaml` - MBD efficiency
- Many others for comprehensive systematics

## Key Workflows

### Training New BDT Models

```bash
cd FunWithxgboost

# 1. Prepare training data
root -l -b -q 'BDTinput.C'

# 2. Edit configuration
nano config_nom.yaml

# 3. Train models
python main_training.py --config config_nom.yaml

# 4. Check output
ls -lh bdt_model*.json bdt_model*.root

# 5. Apply to data
root -l -b -q 'apply_BDT.C'
```

### Running Efficiency Calculation

```bash
cd efficiencytool

# Single configuration
./oneforall.sh config_nom.yaml

# All systematic variations (parallel on Condor)
for config in config_*.yaml; do
  condor_submit oneforall.sub -append "arguments=$config"
done
```

### Generating Final Plots

```bash
cd plotting

# Main results
root -l -q 'plot_final.C'

# All selection plots
./make_all_selection.sh

# Systematic uncertainty plots
./make_syst.sh
```

## Systematic Uncertainties

The analysis includes comprehensive systematic uncertainties:

### Experimental Systematics
1. **Energy Scale**: ±1% calorimeter calibration uncertainty
2. **Energy Resolution**: Smearing variations
3. **Isolation Cut**: Variations of cone energy threshold
4. **Shower Shape Cuts**: Looser/tighter selections
5. **MBD Efficiency**: Trigger efficiency uncertainty
6. **Vertex Selection**: Range variations
7. **BDT Threshold**: Different working points

### Methodology Systematics
1. **Template Fit**: Number of templates in sideband method
2. **Unfolding**: Bayesian iteration count, prior choice
3. **Background Subtraction**: ABCD method variations
4. **MC Sample Merging**: Weighting scheme

### Automation
```bash
cd efficiencytool
# Submit all systematics in parallel
./submit_all_systematics.sh

# After jobs complete
cd ../plotting
root -l -q 'calcSyst.C'  # Aggregate systematic uncertainties
```

## Output Files

### BDT Models
- `FunWithxgboost/bdt_model.json` - XGBoost model (JSON format)
- `FunWithxgboost/bdt_model.root` - TMVA-compatible ROOT format
- `FunWithxgboost/bdt_model_vertex_*.json` - Vertex-binned models
- `FunWithxgboost/binned_models/*.pkl` - pT-binned models (pickle)

### Efficiency Files
- `efficiencytool/RecoEffCalculator.root` - Efficiency histograms and response matrices
- `efficiencytool/rbrQA.root` - Run-by-run QA histograms

### Yield Files
- `efficiencytool/Photon_final_*.root` - Final cross-section results

### Plots
- `plotting/figures/*.pdf` - Publication-ready plots
- `plotting/rootFiles/*.root` - Intermediate ROOT files for plotting

## Documentation

### Analysis Note
Complete analysis documentation in LaTeX:
```bash
cd PPG12-analysis-note
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

### Additional Documentation
- `FunWithxgboost/README.md` - BDT training details
- `FunWithxgboost/REFACTORING_SUMMARY.md` - Recent code improvements
- `FunWithxgboost/VERTEX_FEATURES_README.md` - Vertex binning details
- `FunWithxgboost/QUICK_START.md` - BDT quick start guide

## Physics Definitions

### Isolated Photon
A photon with little nearby hadronic activity:
```
Isolation = Σ(E_cluster in cone R < 0.3) - E_photon
Isolated if: Isolation / E_photon < 0.1
```

### Shower Shape Variables
- **weta, wphi**: Width of EM shower in η and φ
- **chi2**: Shower shape χ² from expected profile
- **e11/e33**: Energy in 1×1 core / 3×3 tower region
- **e32/e35**: Energy in 3×2 / 3×5 asymmetric regions

### Selection Categories
- **Tight Isolated**: Pass tight shower cuts + isolation
- **Tight Non-isolated**: Pass tight shower cuts, fail isolation
- **Non-tight Isolated**: Fail tight shower cuts, pass isolation
- **Non-tight Non-isolated**: Fail both (used for background estimation)

### ABCD Sideband Method
```
        │ Isolated │ Non-isolated
────────┼──────────┼──────────────
Tight   │    A     │      B
        │ (signal) │  (bkg region)
────────┼──────────┼──────────────
Non-    │    C     │      D
tight   │(bkg reg.)│  (bkg region)

Background in A = (B × C) / D
```

## Troubleshooting

### Common Issues

**Issue**: BDT training fails with "No signal events"
- Check input file paths in config
- Verify truth matching in BDTinput.C
- Check pT/eta cuts aren't too restrictive

**Issue**: Efficiency calculation crashes
- Verify all MC samples are merged properly
- Check config file has correct paths
- Ensure sufficient statistics in response matrix

**Issue**: Negative yields after subtraction
- Check sideband purity (need sufficient background in D region)
- May need to adjust shower shape or isolation cuts
- Verify BDT threshold isn't too loose

**Issue**: Poor unfolding closure
- Check response matrix for diagonal dominance
- May need to adjust pT binning
- Verify MC has sufficient statistics

## Contributors

This analysis is part of the sPHENIX PPG12 working group studying isolated photon production.

## References

- sPHENIX: https://www.sphenix.bnl.gov/
- JETPHOX: https://lapth.cnrs.fr/PHOX_FAMILY/
- XGBoost: https://xgboost.readthedocs.io/
- RooUnfold: http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html

## License

Internal sPHENIX collaboration code. Check with collaboration before external use.
