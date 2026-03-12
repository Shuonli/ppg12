# efficiencytool

ROOT/C++ analysis macros for the sPHENIX PPG12 isolated photon cross-section measurement. All macros use `yaml-cpp` for configuration and are compiled/run via the ROOT interpreter.

## Core Analysis Chain

### RecoEffCalculator.C
Main simulation efficiency calculator. Loops over photon/jet simulation trees with cross-section weighting, vertex reweighting, optional JES calibration, and MBD pileup rejection. Applies shower shape selection (tight/non-tight) and isolation cuts (iso/non-iso), computes photon reconstruction efficiencies, fills ABCD-method histograms, and builds RooUnfold response matrices. Outputs per-filetype efficiency and response ROOT files.

### RecoEffCalculator_TTreeReader.C
Production version of `RecoEffCalculator.C`, modernized with `TTreeReader`. Adds NPB timing/BDT score rejection, cluster timing cuts, MBD pileup rejection via `MbdPileupHelper.h`, and MBD t0 corrections. This is the primary macro for computing efficiencies and response matrices.

### MergeSim.C
Merges per-sample-type ROOT output files (from `RecoEffCalculator`) into combined files using `TFileMerger`. Combines photon samples (photon5/10/20) and jet samples (jet10-30) into merged efficiency and response matrix files.

### CalculatePhotonYield.C
Final-stage analysis macro that extracts the isolated photon yield using the ABCD sideband method. Takes ABCD histograms from data and efficiency/purity from simulation, applies efficiency corrections, performs unfolding via RooUnfold, and computes the corrected photon cross section. Supports both data and MC closure modes.

### Closure.C
Performs full and half-closure tests for the RooUnfold unfolding procedure. Full closure unfolds MC-reco with the MC response matrix and compares to MC truth. Half closure splits the MC sample: one half builds the response, the other is unfolded and compared to truth.

### VertexReweighting.C
Computes data/MC vertex-z reweighting ratios. Loads vertex-z histograms from data, photon MC, and jet MC output files, normalizes to unit area, and saves the ratio histograms for use as reweighting factors in other macros.

## Shower Shape and BDT Studies

### ShowerShapeCheck.C
Primary shower shape study macro for both data and simulation. Applies full event selection with MBD pileup helper, MBD t0 corrections, and NPB rejection (timing and/or BDT score). Fills shower shape distributions for tight/non-tight, iso/non-iso categories in eta and ET bins.

### showershape_vertex_check.C
Simulation-only macro evaluating shower shape distributions binned by vertex-z, using both cut-based selection and TMVA RBDT rescoring (NPB score and SS-BDT). Used to study vertex-dependent shower shape systematics.

### NPB_PurityStudy.C
Studies isolated photon purity with respect to non-physical backgrounds (NPB) using the ABCD method. The NPB dimension uses MBD timing spread or NPB BDT score cuts. Outputs run-by-run purity estimates via TGraphErrors.

### DoubleInteractionCheck.C
Studies MBD-based signatures of double (pileup) interactions. Uses TMVA RBDT for NPB rescoring and computes MBD time spread as a pileup proxy. Fills 2D histograms of MBD sigma vs. cluster pT for different selection categories.

## Timing Studies

### plot_cluster_time.C
Comprehensive cluster timing analysis macro. Fills histograms of cluster timing relative to MBD time with optional leading-tower-time corrections. Applies full shower shape selection and fills timing distributions per ET bin for all four ABCD categories.

### plot_cluster_mbd_time_eta.C
Studies cluster timing relative to MBD time as a function of cluster eta. Fills 2D histograms and profiles of (cluster - MBD) time vs. eta in ET bins.

### plot_cluster_jet_time.C
Studies cluster timing vs. recoil jet pT correlation. Applies tight shower shape selection and fills 2D timing histograms split by pT bins.

### merge_cluster_time.C
Merges per-sample-type cluster time analysis ROOT files into combined photon and jet files.

### time_energy_corr.C
Studies the correlation between cluster tower timing and energy within EMCal clusters. Fills 2D histograms of tower energy vs. time for understanding timing-based NPB rejection.

## Truth-Level Studies

### TruthJetEff.C
Computes photon reconstruction efficiencies and fills ABCD histograms for both photon simulation and jet background samples. Also builds RooUnfold response matrices. An alternative efficiency calculator with a different config structure.

### AnalyzeTruthPhotonTowers.C
Loops over simulation trees to study EMCal tower-level properties (energy and timing) for truth-matched photon clusters. Fills 2D/3D histograms of tower energy vs. time, identifying leading and non-leading towers.

### FindETCut.C
Finds the isolation ET cut value that retains a given efficiency fraction (90%, 80%, 70%) as a function of reco cluster pT. Fits the cutoff vs. pT with a linear function to determine the pT-dependent isolation boundary.

### FindTruthETCut.C
Loads 2D truth-level pT vs. isolation ET histograms and plots the fraction of direct and fragmentation photons retained as a function of an isoET cutoff for several pT ranges.

### CheckPhotonPtRange.C
Diagnostic macro that reads tower-level histograms from truth photon tower analysis output and prints statistical summaries (entry counts, energy ranges, mean energies).

## Plotting and QA

### MakeEffPlots.C
Reads pre-computed `TEfficiency` objects from simulation efficiency output and generates plots of reco, ID, and isolation efficiency as functions of truth photon pT, stratified by eta bins.

### MakeSatPlots.C
Reads reco efficiency, saturation probability, and saturation+bad-chi2 probability from single-gamma simulation output and produces plots as functions of truth photon pT by eta bin.

### PlotSpectra.C
Reads BDT-applied simulation output files and fills truth photon pT, truth jet pT, and reco cluster ET spectra with cross-section weighting. Used for verifying pT-hat range population of each sample type.

### draw.C
Utility library providing generic histogram drawing functions (`draw_1D_single_plot`, `draw_1D_multiple_plot`, `draw_1D_multiple_plot_ratio`). Included by other macros for standardized plotting.

### Cluster_rbr.C
Run-by-run QA of EMCal cluster yields. Applies photon selection using YAML-configured shower shape cuts and accumulates cluster counts per run number, divided by per-run luminosity for normalized yield plots.

### Cluster_rbr_TTreeReader.C
TTreeReader-based version of `Cluster_rbr.C`. Same run-by-run luminosity-normalized cluster yield QA with modernized branch reading.

## Development/Testing

### TreeLoadingtest.C
Development/debugging macro that sets up all TTree branch addresses (including BDT score branches) and iterates through the tree printing vertex-z values. Used for testing that branches are correctly readable.

### makeHist_showerShapes_v6.C
Shower shape histogram production macro (by Yeonju Go) for studying EMCal cluster shower shapes. Supports streak event rejection, vertex/MBD reweighting, and optional isolation requirements with configurable pT/eta/vertex binning.
