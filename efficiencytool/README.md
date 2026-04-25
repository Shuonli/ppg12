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

### Double-interaction MC blending (single-pass, condor)

```bash
cd efficiencytool
mkdir -p logs
# 0 mrad
condor_submit list_file=showershape_di_jobs_0rad.list submit_showershape_di.sub
# 1.5 mrad
condor_submit list_file=showershape_di_jobs_1p5rad.list submit_showershape_di.sub
# after all jobs finish
bash hadd_showershape_di.sh config_showershape_0rad.yaml
bash hadd_showershape_di.sh config_showershape_1p5rad.yaml
```

The retired two-pass reco-vertex scheme (`run_showershape_double_reco_legacy.sh`)
is kept on disk for reference only.

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

### Single-pass DI (double-interaction) pipeline — condor

Used when the MC sample is a physical mixture of single-interaction (`nom`) and
double-interaction (`double`) events. Cluster-weighted fractions from
`calc_pileup_range.C` (triple+ folded into double):

| Crossing angle | Run range | `cw_double` | `cw_single` |
|---|---|---|---|
| 0 mrad | 47289–51274 | 0.224 | 0.776 |
| 1.5 mrad | 51274–54000 | 0.079 | 0.921 |

**Design.** Unlike the retired two-pass reco-vertex scheme, vertex reweighting
is now handled inline by `TruthVertexReweightLoader.h`: given a truth vertex
`z_h` (and MB `z_mb` for DI), it returns the data/MC weight from the pre-built
`truth_vertex_reweight/output/{0mrad,1p5mrad}/reweight.root` file. Because the
truth-vertex shape is sample-invariant at first order, the same reweight file
is reused across all MC samples — no refit, no vertex-scan pass, no per-run
`vtxscan_sim_override` wiring. `ShowerShapeCheck.C` already includes and
applies this machinery unconditionally when `truth_vertex_reweight_on: 1` in
the YAML (see `config_showershape_0rad.yaml` / `config_showershape_1p5rad.yaml`).

**Pairs (per crossing angle).** 8 SI/DI sample pairs + 1 data = 17 condor jobs:

| Sample family | `doinclusive` | Pairs |
|---|---|---|
| Photon | `false` | photon{5,10,20}_nom × photon{5,10,20}_double |
| Jet    | `true`  | jet{8,12,20,30,40}_nom × jet{8,12,20,30,40}_double |

The `*_nom` suffix is an alias: `ShowerShapeCheck.C` reads the nominal SI tree
(`photon10`, `jet12`, …) but writes the output under the `*_nom` name so the
SI file does not overwrite the SI-only run from `run_showershape.sh`.
`jet5` is skipped (no DI partner); `jet50_double` is skipped (no analysis use).

**How to run.**
```bash
cd efficiencytool
mkdir -p logs
condor_submit list_file=showershape_di_jobs_0rad.list   submit_showershape_di.sub
condor_submit list_file=showershape_di_jobs_1p5rad.list submit_showershape_di.sub
# after all jobs complete
bash hadd_showershape_di.sh config_showershape_0rad.yaml
bash hadd_showershape_di.sh config_showershape_1p5rad.yaml
```

**What each job does.** `run_showershape_di_job.sh` invokes
```
ShowerShapeCheck.C(config, sample, doinclusive, false, mix_weight)
```
with `do_vertex_scan=false` and no override. `mix_weight` = `SINGLE_FRAC` for
`*_nom` rows, `DOUBLE_FRAC` for `*_double` rows, and `1.0` for data.

**Downstream merge** (`hadd_showershape_di.sh`):
```
signal_combined_{suffix}.root       = hadd(photon{5,10,20}_{nom,double}_{suffix}.root)
jet_inclusive_combined_{suffix}.root = hadd(jet{8,12,20,30,40}_{nom,double}_inclusive_{suffix}.root)
```
Per-event cross-section weights (from `CrossSectionWeights.h`) × `mix_weight`
are already embedded, so plain `hadd` yields the correctly-blended MC.

**Retired legacy script:** `run_showershape_double_reco_legacy.sh` is the old
two-pass reco-vertex scheme (Pass 1 vertex scan → hadd blended vtxscan → Pass 2
full analysis using `vtxscan_sim_override`). It remains on disk for reference
but has been superseded by the single-pass truth-vertex-reweight pipeline above.

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
