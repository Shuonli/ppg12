# Completeness Review -- Phase 2.5 Gap Analysis

**Generated**: 2026-04-08

## 1. Reports Reviewed

The following Phase 1 reports were examined:
- `phase1/anatreemaker.md` -- covers `anatreemaker/` (CaloAna24, Fun4All macros, condor infra)
- `phase1/bdt_pipeline.md` -- covers `FunWithxgboost/` (BDTinput.C, training, apply_BDT.C, configs)
- `phase1/efficiency_tool.md` -- covers `efficiencytool/` (RecoEffCalculator, CalculatePhotonYield, ShowerShapeCheck, DoubleInteractionCheck, MergeSim, systematics)
- `phase1/plotting.md` -- covers `plotting/` (all plot macros, calc_syst_bdt.py, report generators)
- `phase1/analysis_note.md` -- **MISSING** (file does not exist; no analysis note report was produced)

---

## 2. Missing Phase 1 Report: Analysis Note

**Importance: CRITICAL**

The file `wiki/_reports/phase1/analysis_note.md` does not exist. The PPG12 analysis note lives at `PPG12-analysis-note/` and contains:

| File | Description |
|------|-------------|
| `main.tex` | Main document (sPH-ppg-2024-012) |
| `introduction.tex` | Physics motivation, RHIC/sPHENIX context |
| `reconstruction.tex` | Cluster reconstruction, isolation |
| `selection.tex` | Photon selection, BDT, ABCD method |
| `analysis.tex` | Cross-section extraction, unfolding |
| `pileup_mixing.tex` | Double-interaction/pileup studies |
| `double_interaction_efficiency.tex` | DI efficiency correction |
| `systematics.tex` | Systematic uncertainty breakdown |
| `results.tex` | Final cross-section, NLO comparison |
| `conclusion.tex` | Summary |
| `appendix.tex` | Additional material |
| `defs.sty` | Custom LaTeX macros (`\etg`, `\pp`, `\isoET`, etc.) |
| `sphenix.cls` | sPHENIX document class |
| `cite.bib` | BibTeX references |

This is the primary physics document for the analysis and should have a dedicated wiki page.

---

## 3. Completely Missed Directories

The following top-level directories were not covered by any Phase 1 report. They are listed with contents and importance ratings.

### 3.1 `NLO/` -- JETPHOX NLO Theory Comparison

**Importance: CRITICAL**

Contains the NLO pQCD theory predictions used for the final cross-section comparison plot (`plot_final_selection.C`).

| File | Description |
|------|-------------|
| `MakeJetPHOXhisto.C` | Reads JETPHOX NLO output trees, rebins into analysis pT bins, outputs `rootFiles/jetPHOX_{scale}.root` |
| `PlotTruthIso.C` | Studies truth isolation dependence of NLO cross-section across isolation cuts and pT bins |
| `run_jetphox_histos.sh` | Runs MakeJetPHOXhisto for 3 scale variations (0.5x, 1x, 2x mu) in parallel |
| `rootFiles/` | Output ROOT files with rebinned NLO spectra |
| `logs/` | Execution logs |

This directory also contains the Werner Vogelsang NLO data files under `plotting/sphenix_nlo/`:
- `photons_newphenix_sc{05,1,2}.dat` -- photon cross-section at 3 scale variations
- `jets_newphenix_sc{05,1,2}.dat` -- jet cross-section at 3 scale variations
- `pi0_newphenix*.dat` -- pi0 cross-section (for context)
- `nlo_plot.C`, `nlo_plot_count.C` -- NLO plotting macros

Without documenting this directory, the provenance of the NLO theory curves on the final plot is unclear.

### 3.2 `lumi/` -- Luminosity Data

**Importance: CRITICAL**

Contains:
- `60cmLumi_fromJoey.list` -- per-run luminosity values (pb^-1) for trigger bits 10, 18, 22, 30 (corrected and uncorrected). This file is the primary luminosity input for `LumiCalculator.C` and is referenced in YAML configs.

This is a single-file directory but critically important: all cross-section normalization depends on it.

### 3.3 `pythia_reweight/` -- PYTHIA8 Spectrum Reweighting Studies

**Importance: USEFUL**

Standalone PYTHIA8 event generation for studying unfolding prior sensitivity.

| File | Description |
|------|-------------|
| `gen_bins.C` | PYTHIA8 event generation in pThat bins (5 bins: 5-50 GeV), with FastJet R=0.4 jet reco. Outputs truth photon/jet/pi0 spectra per bin. |
| `gen_reweighted.C` | Two-pass PYTHIA8: soft QCD (pThat < 20 GeV, no bias) + hard QCD (pThat >= 20 GeV, pT^4 bias with reweighting). |
| `plot_spectra.C` | Plots and compares cross-sections from both generation modes |
| `output_bins.root`, `output_reweighted.root` | Output ROOT files |

### 3.4 `simcrosssection/` -- MC Cross-Section Validation

**Importance: USEFUL**

Validates MC sample cross-section normalization constants.

| File | Description |
|------|-------------|
| `calculate_average_ratio.C` | Computes average ratio of `sim_cross_counting` bin1/bin2 across condor output files -- validates MC generation efficiency |
| `Fun4AllPythia.C` | PYTHIA8 cross-section extraction from Fun4All |
| `PlotCombine.C` | Overlays photon5/10/20 MC spectra to verify cross-section weighting produces a smooth combined spectrum |
| `draw.C` | Plotting utility |

### 3.5 `eventskimmer/` -- DST Event Skimmer Module

**Importance: USEFUL**

A Fun4All SubsysReco module that filters DST events based on a skip-event list.

| File | Description |
|------|-------------|
| `source/JetDSTSkimmer.cc` | Reads a list of (run, event) pairs to skip; returns DISCARDEVENT for matching events |
| `source/JetDSTSkimmer.h` | Class declaration |
| `source/Makefile.am`, `configure.ac`, `autogen.sh` | Build infrastructure |
| `macro/Fun4All_runDST.C` | Fun4All macro using the skimmer (ana468 production) |
| `macro/run.sh` | Condor submission script |
| `macro/streak_analysis_event_list.txt` | List of events identified as streak/noise artifacts |

This module was used for data quality filtering (removing streaks/noise events) in earlier productions.

### 3.6 `calo_code/` -- Calorimeter Processing Macros

**Importance: USEFUL**

Generic calorimeter analysis macros, not specific to PPG12 but used in the production chain.

| File | Description |
|------|-------------|
| `Fun4All_JetSkimmedProductionYear2.C` | Standard jet-skimmed production macro (Year 2) |
| `Fun4All_Year2_Fitting.C` | Tower waveform fitting macro |
| `Fun4All_Prdf_Combiner.C` | PRDF event combination macro |
| `run2ppcalo_physics_new_2024p024_v001.yaml` | Official sPHENIX production DST config |
| `run_eventcombine.sh`, `run_fitting.sh`, `run_jetskimmer.sh` | Execution scripts |

### 3.7 `dataana/` -- Early Sideband Analysis

**Importance: MINOR**

Legacy/early prototype of the ABCD sideband analysis.

| File | Description |
|------|-------------|
| `SideBand.C` | Early ABCD sideband implementation reading slimtree via YAML config (predates RecoEffCalculator) |
| `results/` | Output directory |

### 3.8 `sideband/` -- Legacy Sideband Studies

**Importance: MINOR**

Another legacy sideband directory with early analysis code.

| File | Description |
|------|-------------|
| `slimtreeana.C` | Early slimtree analysis with hardcoded paths and multiple cluster node comparisons |
| `PlotNTsystematic.C` | Non-tight systematic study |
| `PlotPhotonYield.C` | Early photon yield plotting |
| `PlotSubstruction.C` | Sideband subtraction plotting |

### 3.9 `showershapecheck/` -- Standalone Shower Shape Studies

**Importance: MINOR**

Earlier, standalone shower shape analysis code (predates ShowerShapeCheck.C in efficiencytool).

| File | Description |
|------|-------------|
| `MC_showershape_single_saturationtest.C` | Saturation effects on shower shapes |
| `MC_showershape_treeana.C` | Tree-based shower shape analysis |
| `MC_showershape_trainningprep.C` | Training data preparation from shower shapes |
| `PlotIsoCorrelation.C` | Isolation-shower-shape correlation study |
| `Cluster_timing.C` | Cluster timing analysis |
| `*.csv` | Training data exports |

### 3.10 `saturationcheck/` -- Tower Saturation Studies

**Importance: MINOR**

Studies of ADC saturation effects on calorimeter towers.

| File | Description |
|------|-------------|
| `Fun4All_Saturation_Test.C` | Fun4All macro to test waveform saturation handling |
| `Waveformtest.C` | Waveform analysis |
| `MakeSatPlot.C` | Saturation plotting |
| `draw.C`, `Calo_Calib.C`, `Calo_Fitting.C` | Support macros |

### 3.11 `waveformsimbackup/` -- CaloWaveformSim Backup

**Importance: MINOR**

Backup of a custom `CaloWaveformSim` SubsysReco module (author: Shuhang Li). This is a modified version of the standard sPHENIX waveform simulation code.

### 3.12 `toymcunfold/` -- Toy MC Unfolding Studies

**Importance: MINOR**

Simple toy MC for validating unfolding algorithms.

| File | Description |
|------|-------------|
| `macro/jet_pt_unfolding.C` | Toy jet pT unfolding exercise |
| `macro/resolutionsim.C` | Resolution simulation for unfolding tests |

### 3.13 `photontruthisopythia/` -- Truth Isolation Efficiency

**Importance: MINOR**

Studies truth-level photon isolation efficiency using PYTHIA MC.

| File | Description |
|------|-------------|
| `macro/truthisoeff.C` | Plots truth isolation efficiency for direct vs fragmentation photons at different cone sizes |

### 3.14 `FunWithTMVA/` -- TMVA BDT Training (Legacy)

**Importance: MINOR**

Predecessor to the XGBoost pipeline. Uses ROOT's TMVA for BDT training.

| File | Description |
|------|-------------|
| `TMVA_train.C`, `TMVA_train_nt.C`, `TMVA_train_test.C` | TMVA training macros |
| `run_TMVA.sh`, `run_TMVA_nt.sh` | Execution scripts |
| `oneforall.sh`, `oneforall.sub` | Pipeline and condor submission |
| `sh_makelist.sh` | Dataset preparation |

### 3.15 `ResonancePeakIsoAnalysis/` -- Separate Analysis Module

**Importance: MINOR**

An independent Fun4All analysis module for resonance peak / isolation studies (likely pi0/eta). Not part of the main PPG12 pipeline but shares some code patterns.

| File | Description |
|------|-------------|
| `src/caloTreeGen.cc`, `caloTreeGen.h` | Calorimeter tree generator module |
| `src/ClusterIso.cc`, `ClusterIso.h` | Cluster isolation calculator |
| `macros/AnalyzeTriggerGroupings.*` | Trigger analysis macros |
| `macros/Fun4All_CaloTreeGen.C` | Fun4All driver macro |

### 3.16 Other Uncovered Directories

| Directory | Contents | Importance |
|-----------|----------|------------|
| `config/` | Single `config.yaml` (early analysis config with NO_SPLIT, ana462 era) | MINOR |
| `external/` | Contains `coresoftware/` (sPHENIX framework checkout) | MINOR |
| `figures/` | Output PDF figures from sideband comparisons | MINOR |
| `FM_dedx/` | Empty directory | MINOR |
| `tsanalysis/` | Empty directory | MINOR |
| `reports/` | `mc_purity_correction_investigation.tex` -- a standalone LaTeX report on MC purity correction investigation | USEFUL |
| `proposal_example/` | Example thesis proposal PDF | MINOR |

---

## 4. Second LaTeX Documents (Not in Any Report)

### 4.1 `ppg12_conf_note/` -- Conference Note

**Importance: USEFUL**

A separate conference note (`main.tex`) with its own figures directory. Not documented in any report.

### 4.2 `thesis_proposal/` -- Thesis Proposal Document

**Importance: MINOR**

LaTeX thesis proposal with sections on sPHENIX detector, xjg measurement, Au+Au analysis, and schedule.

### 4.3 `FunWithxgboost/sPHENIX-Photon-Identification-TF-Analysis-Note/` -- Photon ID Technical Note

**Importance: USEFUL**

A separate sPHENIX technical note on photon identification:
- `main.tex`, `introduction.tex`, `selection.tex`, `clustering.tex`, `identification.tex`, `appendix.tex`
- Has its own `README.md`, `defs.sty`, `sphenix.cls`, `cite.bib`

The BDT pipeline report mentions `FunWithxgboost/report.tex` but does NOT cover this subdirectory.

### 4.4 `FunWithxgboost/report.tex` -- Auto-Generated Training Report

**Importance: USEFUL**

Auto-generated LaTeX report from the BDT training pipeline. Mentioned but not detailed in the BDT report.

---

## 5. Files Within Covered Directories That Were Not Mentioned

### 5.1 efficiencytool/ -- Files Not in Phase 1 Report

The efficiency_tool report covers the major macros but misses several secondary files:

| File | Description | Importance |
|------|-------------|------------|
| `RecoEffCalculator.C` | Legacy efficiency calculator (no TTreeReader, no NPB) | MINOR |
| `TreeLoadingtest.C` | Test macro for tree loading | MINOR |
| `time_energy_corr.C` | Time-energy correlation analysis | MINOR |
| `Cluster_rbr_TTreeReader.C` | Run-by-run cluster QA (TTreeReader version) | USEFUL |
| `plot_cluster_mbd_time_eta.C` | Cluster-MBD timing vs eta | MINOR |
| `plot_cluster_jet_time.C` | Cluster-jet timing study | MINOR |
| `merge_cluster_time.C` | TFileMerger for cluster timing output files | MINOR |
| `CheckPhotonPtRange.C` | Histogram range validation for truth photon analysis | MINOR |
| `test_vertex_reweight.C` | Unit tests for vertex reweighting logic | USEFUL |
| `test_mbd_cut.C` | Unit test for MBD pileup rejection | USEFUL |
| `BDTScoreVsET.C` | 2D BDT score vs ET distributions (signal vs inclusive) | USEFUL |
| `TruthJetInETWindow.C` | Truth jet pT in reco ET windows (signal vs background) | USEFUL |
| `PlotTruthJetWindows.C` | Plotting companion for TruthJetInETWindow | USEFUL |
| `showershape_vertex_check.C` | Shower shape vertex-z dependence study | USEFUL |
| `investigate_truth_photons.C` | Truth photon content: single vs double-interaction MC | USEFUL |
| `investigate_vertex_shift.C` | Vertex shift effects in double-interaction MC | USEFUL |
| `investigate_deltaR.C` | Truth-reco matching deltaR in single vs double MC | USEFUL |
| `draw.C` | Legacy drawing utility | MINOR |
| `MakeSatPlots.C` | Saturation plotting | MINOR |
| `MakeEffPlots.C` | Efficiency plotting utility | MINOR |
| `makeHist_showerShapes_v6.C` | Legacy shower shape histogram maker | MINOR |
| `FindTruthETCut.C` | Truth ET cut optimization | MINOR |
| `PlotSpectra.C` | Truth/reco spectra with pT-hat window filtering | USEFUL |
| `calc_pileup_range.C` | Computes aggregate pileup rate across run range from trigger scaler files | USEFUL |
| `modify_notebook.py` | Python script to modify Jupyter notebooks programmatically | MINOR |
| `modify_notebook_pdf.py` | Similar notebook modifier | MINOR |
| `eff_plot.ipynb` | Jupyter notebook for efficiency plotting | MINOR |

**Condor submit files not documented:**
- `oneforall_tree.sub`, `oneforall_tree_jet.sub` -- tree-based pipeline condor submission
- `run_plot_spectra.sub` -- condor submission for PlotSpectra
- `run_showershape.sub` -- condor submission for shower shape analysis

**Shell scripts not documented:**
- `run_all.sh` -- runs all samples
- `run_cluster_time.sh`, `run_tower_time.sh` -- timing analysis drivers
- `run_closure_test.sh` -- closure test execution
- `run_plot_spectra.sh`, `submit_plot_spectra.sh`, `combine_plot_spectra.sh`, `run_plot_spectra_condor.sh` -- PlotSpectra pipeline
- `run_vertex_reweighting.sh`, `run_vertex_check.sh`, `run_vertex_check_signal.sh` -- vertex reweighting
- `run_saturation_study.sh` -- saturation study driver
- `hadd_isoroc.sh`, `run_isoroc.sh` -- isolation ROC pipeline
- `oneforoneforall.sh` -- meta-script that runs oneforall for all configs
- `run_truthjet_windows.sh` -- truth jet window analysis driver
- `run_showershape_inclusive.sh`, `run_showershapephoton.sh`, `run_showershapejet.sh` -- sample-specific shower shape drivers
- `run_showershape_double_auto.sh`, `oneforall_tree_double_auto.sh`, `run_double_auto.sh` -- automated double-interaction pipelines

**Documentation not mentioned:**
- `UNFOLDING_CODE_REVIEW.md` -- unfolding code review notes
- `CLOSURE_TEST_README.md` -- documents closure test bugs and fixes
- `CLOSURE_EXTENSION_SUMMARY.md` -- closure test extension summary
- `PLOTSPECTRA_README.md` -- PlotSpectra usage guide
- `embed_test/` -- embedding test subdirectory with its own ShowerShapeCheck.C and MbdPileupHelper.h
- `reports/` -- contains `plot_report_truth.C`, `plot_report_vertex.C`, and `double_interaction_efficiency_report.tex`
- `config_isoroc.yaml`, `jeteff.yaml` -- specialized analysis configs

### 5.2 FunWithxgboost/ -- Files Not in Phase 1 Report

| File | Description | Importance |
|------|-------------|------------|
| `Fun4All_test_photoncluster.C` | Fun4All macro for testing PhotonClusterBuilder (RawClusterLikelihoodProfile) | MINOR |
| `test_apply_npb_tmva_from_txt.C` | Test macro: applies NPB TMVA model to txt-file input | MINOR |
| `test_refactored.py` | Test script for refactored training code | MINOR |
| `test_vertex_features.py` | Test vertex feature handling | MINOR |
| `test_vertex_reweight_plots.py` | Test vertex reweighting visualization | MINOR |
| `test_2d_correlations.py` | Test 2D correlation plotting | MINOR |
| `test_y_zoom.py` | Test y-axis zoom in plots | MINOR |
| `apply_BDT_condor.sub` | Condor submission for BDT application | USEFUL |
| `submit_split_training.sub` | Condor submission for split-model training | USEFUL |
| `config_vertex_test.yaml` | Test config for vertex features | MINOR |
| `config_nom_split.yaml` | Split-model analysis config | MINOR |
| `config_npb_training_split.yaml` | NPB training config (split variant) | MINOR |
| `split_configs/` | 11 per-model-variant training configs (base, base_vr, base_v0, etc.) | USEFUL |
| `npb_models/npb_score_metadata.yaml`, `npb_score_split_metadata.yaml` | NPB model metadata | MINOR |
| `binned_training.ipynb` | Main training notebook | USEFUL |
| `binned_training_original.ipynb`, `binned_training_refactored.ipynb`, `binned_training_refactored_backup.ipynb` | Notebook variants | MINOR |
| `model_training.ipynb`, `makebdt.ipynb` | Earlier training notebooks | MINOR |
| `practice_sklearn.ipynb`, `test_xgb.ipynb`, `Untitled.ipynb` | Scratch/test notebooks | MINOR |
| `SIG_algo.ipynb`, `SIG_test.ipynb` | Signal algorithm notebooks (possibly financial, not physics) | MINOR |
| `TPC_rate_regression.ipynb` | TPC rate regression (not related to PPG12) | MINOR |

### 5.3 plotting/ -- Files Not in Phase 1 Report

The plotting report is quite comprehensive. The `plotting/sphenix_nlo/` subdirectory (Werner Vogelsang NLO data files and plotting macros) is mentioned as input to `plot_final_selection.C` but not separately documented.

Also not explicitly documented:
- `plotting/showershape_report/` -- auto-generated LaTeX report with per-config subdirectories
- `plotting/selection_report/` -- auto-generated selection report
- `plotting/comparison_report/` -- auto-generated comparison report

These are output artifacts of `make_showershape_report.py`, `make_selection_report.py`, `make_comparison_report.py`.

---

## 6. Root-Level Files Not in Any Report

| File | Description | Importance |
|------|-------------|------------|
| `compare_segments.py` | Compares DST segment numbers between two directories to find mismatches | USEFUL |
| `SYST_PIPELINE_README.md` | Documents the 3-step systematic pipeline (referenced in CLAUDE.md) | USEFUL |
| `PLOTTING_UPDATE_SUMMARY.md` | Describes shower shape plotting updates | MINOR |
| `README.md` | Main repository README | USEFUL |
| `CLAUDE.md` | Claude agent instructions and pipeline overview | MINOR (meta) |

---

## 7. Topic Gap Summary

| Topic | Coverage Status | Importance | Key Files |
|-------|----------------|------------|-----------|
| Analysis note (PPG12-analysis-note/) | NOT COVERED (report missing) | CRITICAL | main.tex + 10 section files |
| NLO theory predictions | NOT COVERED | CRITICAL | NLO/MakeJetPHOXhisto.C, PlotTruthIso.C, plotting/sphenix_nlo/ |
| Luminosity calculation | BRIEFLY MENTIONED in efficiency_tool | CRITICAL | lumi/60cmLumi_fromJoey.list, efficiencytool/LumiCalculator.C |
| Pileup rate computation | NOT COVERED | USEFUL | efficiencytool/calc_pileup_range.C |
| Condor job infrastructure | PARTIALLY COVERED (per directory) | USEFUL | *.sub files across efficiencytool/, FunWithxgboost/ |
| PYTHIA reweighting studies | NOT COVERED | USEFUL | pythia_reweight/ |
| MC cross-section validation | NOT COVERED | USEFUL | simcrosssection/ |
| Conference note | NOT COVERED | USEFUL | ppg12_conf_note/ |
| Photon ID technical note | NOT COVERED | USEFUL | FunWithxgboost/sPHENIX-Photon-Identification-TF-Analysis-Note/ |
| Event skimmer module | NOT COVERED | USEFUL | eventskimmer/ |
| Truth jet pT windowing | NOT COVERED | USEFUL | TruthJetInETWindow.C, PlotTruthJetWindows.C |
| BDT score vs ET visualization | NOT COVERED | USEFUL | BDTScoreVsET.C |
| PlotSpectra analysis | NOT COVERED | USEFUL | PlotSpectra.C + 5 driver scripts |
| Isolation ROC studies | NOT COVERED | USEFUL | IsoROC_calculator.C + hadd_isoroc.sh + run_isoroc.sh |
| Unit tests | NOT COVERED | USEFUL | test_vertex_reweight.C, test_mbd_cut.C |
| Double-interaction investigation | NOT COVERED | USEFUL | investigate_truth_photons.C, investigate_vertex_shift.C, investigate_deltaR.C |
| Legacy TMVA pipeline | NOT COVERED | MINOR | FunWithTMVA/ |
| Legacy sideband analysis | NOT COVERED | MINOR | dataana/, sideband/ |
| Saturation check studies | NOT COVERED | MINOR | saturationcheck/, showershapecheck/ |
| Waveform simulation backup | NOT COVERED | MINOR | waveformsimbackup/ |
| Toy MC unfolding | NOT COVERED | MINOR | toymcunfold/ |

---

## 8. Recommendations

### Must-Have Wiki Pages (CRITICAL)

1. **Analysis Note** -- Write the missing `phase1/analysis_note.md` report covering `PPG12-analysis-note/` structure, sections, custom macros, and figure pipeline.

2. **NLO Theory Predictions** -- Document JETPHOX setup, `MakeJetPHOXhisto.C`, Werner Vogelsang data files, and how theory curves appear in the final plot. This is essential for understanding the physics result.

3. **Luminosity** -- Document `lumi/60cmLumi_fromJoey.list` format, `LumiCalculator.C` logic, and how luminosity values flow into configs and cross-section computation.

### Should-Have Wiki Pages (USEFUL)

4. **Supporting Studies** -- A combined page covering:
   - PYTHIA reweighting (`pythia_reweight/`)
   - MC cross-section validation (`simcrosssection/`)
   - Truth jet windowing (`TruthJetInETWindow.C`)
   - BDT score studies (`BDTScoreVsET.C`)
   - Isolation ROC (`IsoROC_calculator.C`)
   - Pileup rate computation (`calc_pileup_range.C`)

5. **Condor Infrastructure** -- A cross-cutting page documenting all `.sub` files and their orchestration scripts across the repository.

6. **Secondary Documents** -- Brief pages for:
   - Conference note (`ppg12_conf_note/`)
   - Photon ID technical note (`FunWithxgboost/sPHENIX-Photon-Identification-TF-Analysis-Note/`)

### Nice-to-Have (MINOR)

7. **Legacy/Archived Code** -- A single page noting the existence and purpose of `FunWithTMVA/`, `dataana/`, `sideband/`, `showershapecheck/`, `saturationcheck/`, `toymcunfold/`, `waveformsimbackup/` without deep documentation.

---

## 9. Existing Documentation Not Reflected in Wiki

These README/documentation files exist in the repo but their contents are not reflected in any Phase 1 report:

| File | Reflected? |
|------|-----------|
| `README.md` (root) | No |
| `efficiencytool/README.md` | Partially (efficiency report covers same content) |
| `plotting/README.md` | Not checked |
| `SYST_PIPELINE_README.md` | Not reflected |
| `PLOTTING_UPDATE_SUMMARY.md` | Not reflected |
| `efficiencytool/UNFOLDING_CODE_REVIEW.md` | Not reflected |
| `efficiencytool/CLOSURE_TEST_README.md` | Not reflected |
| `efficiencytool/CLOSURE_EXTENSION_SUMMARY.md` | Not reflected |
| `efficiencytool/PLOTSPECTRA_README.md` | Not reflected |
| `FunWithxgboost/README.md` | Partially (BDT report covers same content) |
| `FunWithxgboost/REFACTORING_README.md` | Not reflected |
| `FunWithxgboost/REFACTORING_SUMMARY.md` | Not reflected |
| `FunWithxgboost/VERTEX_FEATURES_README.md` | Not reflected |
| `FunWithxgboost/QUICK_START.md` | Not reflected |
| `FunWithxgboost/RIDGE_REPRODUCTION_README.md` | Not reflected (non-physics) |
| `FunWithxgboost/TEST_SET_EVALUATION_UPDATE.md` | Not reflected |
| `FunWithxgboost/OPTUNA_HYPERPARAMETER_TUNING_IMPLEMENTATION.md` | Not reflected |
| `FunWithxgboost/README_split_training.md` | Not reflected |
| `reports/mc_purity_correction_investigation.tex` | Not reflected |
