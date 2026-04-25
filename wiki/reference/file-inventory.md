# File Inventory

Every significant file in the repository with a one-line purpose.

## anatreemaker/

| File | Purpose |
|------|---------|
| `source/CaloAna24.cc` | Main Fun4All module: DST to slimtree conversion (2321 lines) |
| `source/CaloAna24.h` | Class declaration and member variables (385 lines) |
| `source/Makefile.am` | Autotools build file |
| `source/autogen.sh` | Build bootstrap script |
| `macro_maketree/data/ana521/Fun4All_runDST.C` | Fun4All macro for data (current production) |
| `macro_maketree/data/ana521/run.sh` | Condor orchestration for data processing |
| `macro_maketree/data/ana521/CondorRun.sh` | Per-job condor wrapper (data) |
| `macro_maketree/sim/run28/{sample}/Fun4All_run_sim.C` | Fun4All macro for simulation |
| `macro_maketree/sim/run28/{sample}/run_condor.sh` | Condor submission for sim |
| `macro_maketree/sim/run28/hadd_combined.sh` | Merge per-job sim outputs into combined.root |
| `macro_maketree/sim/create_file_lists.sh` | Generate DST file lists for all sim samples |
| `macro_maketree/sim/run28/cleanup_and_run.sh` | Remove old output and resubmit all samples |

## FunWithxgboost/

| File | Purpose |
|------|---------|
| `BDTinput.C` | Feature extraction: slimtree to training text files (1040 lines) |
| `main_training.py` | BDT training orchestrator (1295 lines) |
| `data_loader.py` | Data loading and signal/background labeling (338 lines) |
| `model_builder.py` | sklearn Pipeline construction (77 lines) |
| `reweighting.py` | Kinematic reweighting: class, ET, eta, vertex (165 lines) |
| `plotting.py` | Training diagnostic plots (1222 lines) |
| `apply_BDT.C` | Apply TMVA BDT scores to slimtrees — dual cluster-variant (split + nosplit) |
| `train_npb_score.py` | NPB score BDT training (1073 lines) |
| `config.yaml` | BDT training configuration |
| `config_nom.yaml` | Analysis config for BDTinput/apply_BDT |
| `config_npb_training_{split,nosplit}.yaml` | Per-variant NPB score training configs |
| `make_variant_configs.py` | Generate per-model per-variant training configs (`--variant {split,nosplit}`) |
| `train_variant_model.py` | CLI wrapper for training one model from a variant config |
| `submit_variant_training.sub` | Condor submit (default `VARIANT=split`, override via `-append`) |
| `run_variant_training_local.sh` | Non-Condor local driver, `./run_variant_training_local.sh {split\|nosplit}` |
| `export_cluster_arrays_csv.C` | Export tower arrays to CSV (194 lines) |

## efficiencytool/

### Core Analysis Macros

| File | Purpose |
|------|---------|
| `RecoEffCalculator_TTreeReader.C` | Two-pass efficiency calculator with ABCD (2719 lines) |
| `CalculatePhotonYield.C` | ABCD purity, unfolding, cross-section (1170 lines) |
| `MergeSim.C` | Merge per-sample MC outputs (94 lines) |
| `ShowerShapeCheck.C` | Shower shape distributions in ABCD regions (1917 lines) |
| `DoubleInteractionCheck.C` | Toy vertex-shift pileup simulation (1446 lines) |
| `CrossSectionWeights.h` | All MC cross-section constants (119 lines) |
| `MbdPileupHelper.h` | MBD timing pileup metric calculator (207 lines) |

### Study Macros

| File | Purpose |
|------|---------|
| `Closure.C` | MC closure test for pipeline validation |
| `RecoEffCalculator.C` | Older efficiency calculator (pre-TTreeReader version) |
| `FindBDTCut.C` | BDT threshold optimization scan |
| `FindETCut.C` | Isolation ET cut parameter optimization |
| `NPB_PurityStudy.C` | NPB purity across crossing angle periods |
| `PlotROC.C` | ROC curves for BDT and shower shapes |
| `VertexReweighting.C` | Vertex z reweighting histogram generation |
| `SaturationStudy.C` | ADC saturation effects |
| `Cluster_rbr.C` | Run-by-run cluster QA |
| `LumiCalculator.C` | Luminosity from trigger scalers |
| `IsoROC_calculator.C` | Isolation cut ROC curves |
| `AnalyzeTruthPhotonTowers.C` | Tower-level truth photon analysis |
| `TruthJetEff.C` | Truth jet efficiency studies |
| `TruthJetInETWindow.C` | Truth jets within cluster ET windows |
| `BDTScoreVsET.C` | BDT score vs cluster ET visualization |
| `PlotSpectra.C` | MC spectrum overlay and validation |
| `MakeEffPlots.C` | Efficiency plotting utility |
| `MakeSatPlots.C` | Saturation study plotting |
| `PlotTruthJetWindows.C` | Truth jet pT window visualization |
| `CheckPhotonPtRange.C` | Truth photon pT range validation |
| `time_energy_corr.C` | Tower timing vs energy correction (contains 16.67 ns bug) |
| `calc_pileup_range.C` | Pileup rate calculation for crossing angle periods |
| `plot_cluster_time.C` | Cluster timing distributions |
| `investigate_truth_photons.C` | Truth photon debugging |
| `investigate_vertex_shift.C` | Vertex shift studies |
| `investigate_deltaR.C` | DeltaR matching studies |
| `showershape_vertex_check.C` | Vertex-dependent shower shape study |

### Scripts

| File | Purpose |
|------|---------|
| `oneforall.sh` | Orchestrate MergeSim + CalculatePhotonYield |
| `oneforall_tree.sh` | Two-pass RecoEffCalculator orchestration |
| `oneforall.sub` | HTCondor submit file (one job per config) |
| `submit_showershape_di.sub` | HTCondor submit for single-pass showershape DI pipeline (17 jobs / crossing angle) |
| `run_showershape_di_job.sh` | Per-job executable for the showershape DI pipeline |
| `showershape_di_jobs_{0rad,1p5rad}.list` | 17-row job lists (8 SI/DI pairs + 1 data) consumed by `submit_showershape_di.sub` |
| `run_showershape_double_reco_legacy.sh` | Retired two-pass reco-vertex blending (kept for history) |
| `run_double_interaction.sh` | Toy double-interaction simulation |
| `TruthVertexReweightLoader.h` | Loads `h_w_iterative` per-period reweight histogram, applies `w(z)` / `w(z_hard)*w(z_mb)` |
| `truth_vertex_reweight/fit_truth_vertex_reweight.py` | Iterative data-driven fit producing per-period `reweight.root` |
| `truth_vertex_reweight/output/{0mrad,1p5mrad}/reweight.root` | Per-period truth-vertex reweight histograms consumed by `ShowerShapeCheck.C` and `RecoEffCalculator_TTreeReader.C` |

### Python

| File | Purpose |
|------|---------|
| `make_bdt_variations.py` | Generate 30+ systematic variation configs (336 lines) |
| `make_showershape_variations.py` | Shower shape config variants |

### YAML Configs

37 `config_bdt_*.yaml` files (systematic variations), plus `config_showershape_*.yaml` files.

## plotting/

### Shared Utilities

| File | Purpose |
|------|---------|
| `plotcommon.h` | Shared pT bins, frames, legend strings (92 lines) |
| `BlairUtils.C` | Error-band construction, TGraph helpers (543 lines) |

### Final Results

| File | Purpose |
|------|---------|
| `plot_final_selection.C` | Primary final cross-section plot (816 lines) |
| `calc_syst_bdt.py` | Automated systematic aggregation (534 lines) |

### Analysis Plots

| File | Purpose |
|------|---------|
| `plot_efficiency.C` | Reco/iso/ID/total/MBD efficiencies (258 lines) |
| `plot_purity_selection.C` | Purity with fit bands (114 lines) |
| `plot_mc_purity_correction.C` | MC purity correction ratio visualization |
| `plot_sideband_selection.C` | ABCD yield plots (235 lines) |
| `plot_showershapes_variations.C` | Per-config shower shapes with chi2/ndf (670 lines) |
| `plot_showershapes_selections.C` | Shower shapes with selection cuts and DI mixing |
| `plot_unfold_iter.C` | Unfolding convergence (356 lines) |
| `plot_reweight.C` | Unfolding prior comparison |
| `Closure.C` | Full/half closure test driver (in plotting/) |
| `plot_isoET.C` | Isolation ET distributions (310 lines) |
| `CONF_plots.C` | Isolation ET template fits for conference plots |
| `plot_background_recoisoET_overlay.C` | Background reco isoET overlay |
| `plot_double_interaction.C` | Double-interaction visualization (630 lines) |
| `plot_npb_time_bkgsub.C` | NPB timing with background subtraction |
| `plot_response.C` | Response matrix visualization (58 lines) |
| `plot_bdt_variations.C` | Cross-section comparison across BDT model variants |
| `chi2_weight_scan.C` | Double-interaction weight chi2 scan (835 lines) |
| `plot_vertex_check.C` | Vertex distribution studies |
| `plot_SB.C` | Signal/background ratio vs ET |
| `plot_cluster_timing.C` | Cluster timing distributions |
| `plot_cluster_mbd_time.C` | Cluster-MBD timing correlations |
| `plot_tower_timing.C` | Tower-level timing studies (1695 lines) |
| `plot_cluster_rbr_QA.C` | Run-by-run cluster QA |
| `plot_mbd_sigma_efficiency.C` | MBD sigma selection efficiency |

### Report Generators

| File | Purpose |
|------|---------|
| `make_showershape_report.py` | LaTeX shower shape report (447 lines) |
| `make_selection_report.py` | LaTeX selection/efficiency report (413 lines) |
| `make_comparison_report.py` | Side-by-side comparison report (456 lines) |

### Batch Scripts

| File | Purpose |
|------|---------|
| `make_all_bdt_selection.sh` | Run selection plots for all configs |
| `make_selection_plots.sh` | Run 11 plot macros for one config |
| `make_syst.sh` | Run legacy C++ systematic macros |

## PPG12-analysis-note/

| File | Purpose |
|------|---------|
| `main.tex` | Master document (sPH-ppg-2024-012, version 3) |
| `introduction.tex` | Physics motivation |
| `selection.tex` | Data/MC samples, trigger, luminosity |
| `reconstruction.tex` | Cluster reconstruction, isolation, BDT |
| `analysis.tex` | Cross-section extraction, ABCD, unfolding |
| `systematics.tex` | Systematic uncertainty breakdown |
| `results.tex` | Final cross-section, NLO comparison |
| `conclusion.tex` | Summary |
| `appendix.tex` | Additional material |
| `pileup_mixing.tex` | Double-interaction studies |
| `double_interaction_efficiency.tex` | DI efficiency correction |
| `defs.sty` | Custom LaTeX macros |
| `sphenix.cls` | sPHENIX document class |
| `cite.bib` | References |

## Other Directories

| Directory | Purpose |
|-----------|---------|
| `NLO/` | JETPHOX NLO theory predictions (MakeJetPHOXhisto.C) |
| `lumi/` | Per-run luminosity data (60cmLumi_fromJoey.list) |
| `pythia_reweight/` | PYTHIA8 reweighting studies for unfolding prior |
| `simcrosssection/` | MC cross-section validation |
| `eventskimmer/` | DST event skimmer for data quality filtering |
| `calo_code/` | Generic calorimeter processing macros |
| `dataana/` | Legacy ABCD sideband prototype |
| `sideband/` | Legacy sideband studies |
| `showershapecheck/` | Standalone shower shape studies (legacy) |
| `toymcunfold/` | Toy MC unfolding validation |
| `FunWithTMVA/` | Legacy TMVA-based training (obsolete) |
