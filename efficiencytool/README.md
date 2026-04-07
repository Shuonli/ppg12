# efficiencytool

ROOT/C++ analysis macros for the sPHENIX PPG12 isolated photon cross-section measurement. All macros use `yaml-cpp` for configuration and are compiled/run via the ROOT interpreter.

## Quick Start

### Standard shower shape analysis

```bash
# Run with default config (config_showershape.yaml)
bash run_showershape.sh

# Run with a specific config
bash run_showershape.sh config_showershape_nom.yaml
```

### Double-interaction MC blending (Au+Au pileup)

```bash
bash run_showershape_double.sh [config_showershape.yaml]
```

### Generate shower shape systematic variation configs

```bash
python make_showershape_variations.py config_showershape.yaml [--outdir DIR]
```

---

## Core Analysis Chain

### RecoEffCalculator_TTreeReader.C

Production efficiency calculator, modernized with `TTreeReader`. Loops over photon/jet simulation trees with cross-section weighting, vertex reweighting, optional JES calibration, and MBD pileup rejection. Applies shower shape selection (tight/non-tight) and isolation cuts (iso/non-iso), computes photon reconstruction efficiencies, fills ABCD-method histograms, and builds RooUnfold response matrices. Adds NPB timing/BDT score rejection, cluster timing cuts, and MBD t0 corrections via `MbdPileupHelper.h`.

**Function signature:**
```cpp
void RecoEffCalculator_TTreeReader(const std::string &configname = "config_showershape.yaml", ...)
```

### RecoEffCalculator.C

Earlier version of the efficiency calculator (no `TTreeReader`, no NPB BDT). Retained for reference.

### MergeSim.C

Merges per-sample-type ROOT output files from `RecoEffCalculator` into combined files using `TFileMerger`. Combines photon samples (photon5/10/20) and jet samples (jet10–30).

### CalculatePhotonYield.C

Final-stage macro that extracts the isolated photon yield using the ABCD sideband method. Takes ABCD histograms from data and efficiency/purity from simulation, applies efficiency corrections, performs unfolding via RooUnfold, and computes the corrected photon cross section. Supports data and MC closure modes.

### Closure.C

Full and half-closure tests for the RooUnfold unfolding procedure. Full closure unfolds MC-reco with the MC response matrix and compares to MC truth. Half closure splits the MC sample: one half builds the response, the other is unfolded and compared to truth.

### VertexReweighting.C

Computes data/MC vertex-z reweighting ratios. Loads vertex-z histograms from data, photon MC, and jet MC output files, normalizes to unit area, and saves ratio histograms for use as reweighting factors in other macros.

---

## Shower Shape Studies

### ShowerShapeCheck.C

Primary shower shape study macro for both data and simulation. Applies full event selection with `MbdPileupHelper`, MBD t0 corrections, and NPB rejection (timing and/or BDT score). Fills shower shape distributions for tight/non-tight, iso/non-iso categories in eta and ET bins.

**Function signature:**
```cpp
void ShowerShapeCheck(
    const std::string &configname        = "config_showershape.yaml",
    const std::string  filetype          = "jet12_double",
    bool               doinclusive       = true,
    bool               do_vertex_scan    = false,
    float              mix_weight        = 1.0,
    const std::string  vtxscan_sim_override = ""
)
```

**Parameters:**
| Parameter | Description |
|-----------|-------------|
| `configname` | YAML config file (cuts, input paths, BDT model, etc.) |
| `filetype` | Sample key: `"data"`, `"photon10"`, `"jet12"`, `"photon10_nom"`, `"photon10_double"`, etc. |
| `doinclusive` | If true, includes all cluster kinematics; if false, applies tighter photon selection |
| `do_vertex_scan` | If true, fills only `h_vertexz` (fast Pass 1 for vertex reweighting) |
| `mix_weight` | Multiplicative weight applied to all histogram fills — used to embed single/double interaction blending fractions |
| `vtxscan_sim_override` | If non-empty, overrides the vtxscan sim file path used for vertex reweighting in Pass 2 (allows using the pre-merged combined vtxscan) |

### run_showershape.sh

Standard two-pass pipeline for single-interaction MC samples:

- **Pass 1** (vertex scan): runs all samples in parallel with `do_vertex_scan=true` to fill vertex-z distributions.
- **Pass 2** (full analysis): runs all samples in parallel using Pass 1 output for vertex reweighting.
- After Pass 2: `hadd` merges photon5/10/20 → `MC_efficiencyshower_shape_signal_*.root` and jet5/8/12/20/40 → `MC_efficiencyshower_shape_jet_inclusive_*.root`.

### run_showershape_double.sh — Double-Interaction MC Blending

Used when the Au+Au MC sample is a physical mixture of single-interaction (`nom`) and double-interaction (`double`) events. Current fractions:

```
SINGLE_FRAC = 0.813   (nom)
DOUBLE_FRAC = 0.187   (double)
```

**Pipeline (three stages):**

**Stage 1 — Pass 1 (vertex scan, parallel):**
Each MC sub-sample is run with its physical fraction as `mix_weight`, so that when the vtxscan outputs are added together the result represents the correctly-weighted blended vertex distribution:

```
ShowerShapeCheck(config, "photon10_double", false, true, 0.187)
ShowerShapeCheck(config, "jet12_double",    true,  true, 0.187)
ShowerShapeCheck(config, "photon10_nom",    false, true, 0.813)
ShowerShapeCheck(config, "jet12_nom",       true,  true, 0.813)
ShowerShapeCheck(config, "data",            false, true)
```

**Stage 2 — hadd vtxscan files:**
Merges nom and double vtxscan outputs into one combined vtxscan per sample type:
```
photon10_combined_*_vtxscan.root = photon10_nom_*_vtxscan + photon10_double_*_vtxscan
jet12_combined_inclusive_*_vtxscan.root = jet12_nom_*_vtxscan + jet12_double_*_vtxscan
```
Because `mix_weight` was already embedded, plain `hadd` yields the correctly blended distribution.

**Stage 3 — Pass 2 (full analysis, parallel):**
Both nom and double use the **same combined vtxscan** (via `vtxscan_sim_override`) so that the vertex reweighting ratio is computed against the blended MC distribution. Each sub-sample applies its own `mix_weight` to all histogram fills:

```
ShowerShapeCheck(config, "photon10_double", false, false, 0.187, photon10_combined_vtxscan)
ShowerShapeCheck(config, "jet12_double",    true,  false, 0.187, jet12_combined_vtxscan)
ShowerShapeCheck(config, "photon10_nom",    false, false, 0.813, photon10_combined_vtxscan)
ShowerShapeCheck(config, "jet12_nom",       true,  false, 0.813, jet12_combined_vtxscan)
ShowerShapeCheck(config, "data",            false)
```

**Stage 4 — hadd final outputs:**
```
photon10_combined_*.root = photon10_nom_*.root + photon10_double_*.root
jet12_combined_inclusive_*.root = jet12_nom_inclusive_*.root + jet12_double_inclusive_*.root
```
Since `mix_weight` was already embedded per-event, plain `hadd` is correct.

**Output files** land in `results/` and follow the naming pattern:
`MC_efficiencyshower_shape_{sample}_{suffix}.root`

### make_showershape_variations.py

Generates systematic variation YAML configs from a nominal base config. Each entry in `VARIANTS` produces one `config_showershape_{name}.yaml`. Supported override keys include `vertex_cut`, `run_min/max`, `lumi`, timing cuts (`cluster_mbd_time_min/max`), MBD sigma cut, iso scale/shift, NPB score cut, and others. See the script header for the full list of supported keys and their YAML paths.

### showershape_vertex_check.C

Simulation-only macro evaluating shower shape distributions binned by vertex-z, using both cut-based selection and TMVA RBDT rescoring (NPB score and SS-BDT). Used to study vertex-dependent shower shape systematics.

---

## Additional Macros

### NPB_PurityStudy.C

Studies isolated photon purity with respect to non-physical backgrounds (NPB) using the ABCD method. The NPB dimension uses MBD timing spread or NPB BDT score cuts. Outputs run-by-run purity estimates via `TGraphErrors`.

### DoubleInteractionCheck.C

Studies MBD-based signatures of double (pileup) interactions. Uses TMVA RBDT for NPB rescoring and computes MBD time spread as a pileup proxy. Fills 2D histograms of MBD sigma vs. cluster pT for different selection categories.

### FindETCut.C

Finds the isolation ET cut value that retains a given efficiency fraction (90%, 80%, 70%) as a function of reco cluster pT. Fits the cutoff vs. pT with a linear function to determine the pT-dependent isolation boundary.

### FindTruthETCut.C

Loads 2D truth-level pT vs. isolation ET histograms and plots the fraction of direct and fragmentation photons retained as a function of an isoET cutoff for several pT ranges.

---

## Timing Studies

### plot_cluster_time.C

Comprehensive cluster timing analysis macro. Fills histograms of cluster timing relative to MBD time with optional leading-tower-time corrections. Applies full shower shape selection and fills timing distributions per ET bin for all four ABCD categories.

### plot_cluster_mbd_time_eta.C

Studies cluster timing relative to MBD time as a function of cluster eta. Fills 2D histograms and profiles of (cluster − MBD) time vs. eta in ET bins.

### plot_cluster_jet_time.C

Studies cluster timing vs. recoil jet pT correlation. Applies tight shower shape selection and fills 2D timing histograms split by pT bins.

### merge_cluster_time.C

Merges per-sample-type cluster time analysis ROOT files into combined photon and jet files.

### time_energy_corr.C

Studies the correlation between cluster tower timing and energy within EMCal clusters. Fills 2D histograms of tower energy vs. time for understanding timing-based NPB rejection.

---

## Plotting and QA

### Cluster_rbr.C / Cluster_rbr_TTreeReader.C

Run-by-run QA of EMCal cluster yields. Applies photon selection using YAML-configured shower shape cuts and accumulates cluster counts per run number, divided by per-run luminosity for normalized yield plots. `_TTreeReader` version uses the modernized branch reading API.

### MakeEffPlots.C

Reads pre-computed `TEfficiency` objects from simulation efficiency output and generates plots of reco, ID, and isolation efficiency as functions of truth photon pT, stratified by eta bins.

### MakeSatPlots.C

Reads reco efficiency, saturation probability, and saturation+bad-chi2 probability from single-gamma simulation output and produces plots as functions of truth photon pT by eta bin.

### PlotSpectra.C

Reads BDT-applied simulation output files and fills truth photon pT, truth jet pT, and reco cluster ET spectra with cross-section weighting. Used for verifying pT-hat range population of each sample type.

### draw.C

Utility library providing generic histogram drawing functions (`draw_1D_single_plot`, `draw_1D_multiple_plot`, `draw_1D_multiple_plot_ratio`). Included by other macros for standardized plotting.

### makeHist_showerShapes_v6.C

Shower shape histogram production macro (by Yeonju Go) for studying EMCal cluster shower shapes. Supports streak event rejection, vertex/MBD reweighting, and optional isolation requirements with configurable pT/eta/vertex binning.

---

## Truth-Level Studies

### TruthJetEff.C

Computes photon reconstruction efficiencies and fills ABCD histograms for both photon simulation and jet background samples. Also builds RooUnfold response matrices. An alternative efficiency calculator with a different config structure.

### AnalyzeTruthPhotonTowers.C

Loops over simulation trees to study EMCal tower-level properties (energy and timing) for truth-matched photon clusters. Fills 2D/3D histograms of tower energy vs. time, identifying leading and non-leading towers.

### CheckPhotonPtRange.C

Diagnostic macro that reads tower-level histograms from truth photon tower analysis output and prints statistical summaries (entry counts, energy ranges, mean energies).

---

## Configuration

All macros load a YAML config file (default `config_showershape.yaml`). Key sections:

| Section | Purpose |
|---------|---------|
| `input` | Input file paths, BDT model name, tree/node names |
| `output` | Output file path prefixes and variant type label |
| `analysis` | Eta bins, vertex cut, timing cuts, MBD sigma cut, iso scale/shift, NPB score cut, run range, luminosity |

Use `make_showershape_variations.py` to generate systematic variation configs from a nominal base.
