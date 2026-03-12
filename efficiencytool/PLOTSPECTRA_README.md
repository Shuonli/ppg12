# PlotSpectra Analysis Scripts

Scripts for analyzing truth and reconstructed cluster ET/pT spectra from sPHENIX simulation samples.

## Overview

The `PlotSpectra.C` macro creates weighted ET/pT spectra with optional event-level filtering by maximum photon or jet pT. These scripts help run the analysis across all simulation samples.

## Files

### Core Macro
- **PlotSpectra.C** - Main ROOT macro for creating spectra

### Execution Scripts
- **submit_plot_spectra.sh** - Main submission script (Condor or local)
- **run_plot_spectra.sh** - Sequential bash script for all samples
- **run_plot_spectra_condor.sh** - Condor worker script
- **run_plot_spectra.sub** - Condor submission file
- **combine_plot_spectra.sh** - Combine output files with hadd

## Quick Start

### Method 1: Condor (Recommended - Parallel Execution)

```bash
# Submit all samples to Condor (runs in parallel)
./submit_plot_spectra.sh condor

# Monitor jobs
condor_q

# Check logs
tail -f logs/plot_spectra_photon5.out

# After jobs complete, combine output
./combine_plot_spectra.sh
```

### Method 2: Local Sequential Execution

```bash
# Run all samples sequentially (slow!)
./submit_plot_spectra.sh local

# Or run specific sample types
./run_plot_spectra.sh photon    # Only photon samples
./run_plot_spectra.sh jet        # Only jet samples
./run_plot_spectra.sh nofilter   # All samples, no pT filters
```

### Method 3: Individual Sample

```bash
# Run single sample directly
root -l 'PlotSpectra.C("photon5", 5, 12)'
```

## Sample-Specific pT Filters

The default filters applied for each sample:

| Sample | Photon pT Filter | Jet pT Filter | Cross-section Weight |
|--------|------------------|---------------|---------------------|
| photon5 | 5 - 12 GeV | - | 1121.99 |
| photon10 | 12 - 25 GeV | - | 53.25 |
| photon20 | 25 - 100 GeV | - | 1.0 |
| jet10 | - | 10 - 15 GeV | 549832 |
| jet15 | - | 15 - 20 GeV | 56042 |
| jet20 | - | 20 - 30 GeV | 8554 |
| jet30 | - | 30 - 50 GeV | 344 |
| jet50 | - | 50 - 100 GeV | 1.0 |

## Output Files

### Individual Sample Files
Format: `PlotSpectra_{sample}_{filters}.root`

Examples:
- `PlotSpectra_photon5_photon5to12.root`
- `PlotSpectra_jet10_jet10to15.root`
- `PlotSpectra_photon10.root` (no filters)

### Combined Files (after running combine script)
- `PlotSpectra_photon_combined.root` - All photon samples
- `PlotSpectra_jet_combined.root` - All jet samples
- `PlotSpectra_all_combined.root` - All samples

## Histograms Created

Each output file contains:

### Main Spectra
- `h_truth_photon_pt` - All truth photon pT
- `h_max_truth_photon_pt` - Maximum photon pT per event
- `h_truth_jet_pt` - All truth jet pT
- `h_max_truth_jet_pt` - Maximum jet pT per event
- `h_reco_cluster_et` - Reconstructed cluster ET

### QA Histograms
- `h_vertex_z` - Vertex z distribution
- `h_nparticles` - Number of truth particles
- `h_nclusters` - Number of reco clusters
- `h_event_counts` - Event selection cutflow

### Metadata (TNamed objects)
- `filetype` - Sample name
- `cross_section_weight` - Weight applied
- `min_photon_pt`, `max_photon_pt` - Photon filter range
- `min_jet_pt`, `max_jet_pt` - Jet filter range

## Usage Examples

### Run specific samples
```bash
# Just photon samples
./run_plot_spectra.sh photon

# Just jet samples
./run_plot_spectra.sh jet
```

### Submit to Condor and monitor
```bash
./submit_plot_spectra.sh condor
condor_q -nobatch
condor_tail -f logs/plot_spectra_photon5.log
```

### Combine and analyze
```bash
# After jobs complete
./combine_plot_spectra.sh

# Inspect combined file
root -l PlotSpectra_photon_combined.root
root [1] .ls
root [2] h_truth_photon_pt->Draw()
```

### Check metadata
```bash
root -l PlotSpectra_photon5_photon5to12.root
root [1] TNamed *meta = (TNamed*)_file0->Get("cross_section_weight")
root [2] meta->GetTitle()
```

## Event Selection

Applied cuts (in order):
1. Vertex cut: |vertexz_truth| < 30 cm
2. Find max photon pT (|eta| < 0.7)
3. Find max jet pT
4. Apply photon pT filter (if specified)
5. Apply jet pT filter (if specified)

Only events passing all cuts have histograms filled.

## Troubleshooting

### No output files
- Check if input files exist: `/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run21/{sample}/condorout_waveform/caloana0220.root`
- Verify ROOT environment is set up

### Condor jobs stuck
```bash
condor_q -better-analyze JOB_ID
condor_rm JOB_ID  # Cancel if needed
```

### Memory issues
- Increase `request_memory` in `run_plot_spectra.sub` (default: 4GB)

### Check logs
```bash
# Output and errors
cat logs/plot_spectra_photon5.out
cat logs/plot_spectra_photon5.err

# Condor system log
cat logs/plot_spectra_photon5.log
```

## Advanced Usage

### Custom filters
```bash
# Run macro directly with custom cuts
root -l 'PlotSpectra.C("photon5", 7, 10, -1, -1)'  # 7-10 GeV photons
root -l 'PlotSpectra.C("jet20", -1, -1, 22, 28)'   # 22-28 GeV jets
```

### No filters (all events)
```bash
./run_plot_spectra.sh nofilter
# Or individually:
root -l 'PlotSpectra.C("photon5", -1, -1, -1, -1)'
```

### Process subset on Condor
Edit `run_plot_spectra.sub` to comment out unwanted samples, then:
```bash
condor_submit run_plot_spectra.sub
```

## Notes

- Cross-section weights are already applied in histograms
- Photon samples normalized to photon20
- Jet samples normalized to jet50
- All samples use `slimtree` from `caloana0220.root`
- Branch naming: `CLUSTERINFO_CEMC_NO_SPLIT`
