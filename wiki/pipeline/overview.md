# Pipeline Overview

## End-to-End Flow

```
DST files (data: ana521_2025p007_v001, sim: run28 MC)
  |
  |  [CaloAna24.cc via Fun4All]
  v
slimtree ROOT files (tree: slimtree)
  |
  +-- [BDTinput.C] --> training text files (shapes_*.txt)
  |                        |
  |                        +-- [main_training.py] --> BDT models (binned_models/*.root)
  |                                                       |
  +-- [apply_BDT.C] <------------------------------------+
  |         |
  |         v
  |   scored slimtrees (*_with_bdt_split.root)
  |         |
  |         +-- [RecoEffCalculator_TTreeReader.C, two-pass]
  |                    |
  |                    v
  |              per-sample efficiency + ABCD ROOT files
  |                    |
  |                    +-- [MergeSim.C]
  |                           |
  |                           v
  |                     merged MC files
  |                           |
  |                           +-- [CalculatePhotonYield.C]
  |                                      |
  |                                      v
  |                                Photon_final_*.root
  |                                      |
  |                 +--------------------+-------------------+
  |                 v                                        v
  |    [calc_syst_bdt.py]                     [plot_final_selection.C]
  |           |                                        |
  |           v                                        v
  |    syst_sum.root --------------------------> figures/*.pdf
  |                                                    |
  +---------------------------------------------> analysis note
```

## Stage Summaries

### Stage 1: Tree Making (`anatreemaker/`)

CaloAna24.cc is a Fun4All SubsysReco that reads sPHENIX DST files and produces compact ROOT trees (`slimtree`). It extracts cluster kinematics, shower shapes, isolation energies, truth particles, jets, MBD timing, and trigger info. Data uses production tag `ana521_2025p007_v001`; simulation uses run 28 MC.

Key file: `anatreemaker/source/CaloAna24.cc` (2321 lines)

### Stage 2: BDT Training (`FunWithxgboost/`)

Two sub-stages: (a) `BDTinput.C` extracts features from slimtrees to text files, (b) `main_training.py` trains XGBoost models in ET bins. 11 model variants with different feature subsets are trained. A separate NPB score model is trained via `train_npb_score.py`.

Key files: `BDTinput.C`, `main_training.py`, `config.yaml`

### Stage 3: BDT Application (`FunWithxgboost/apply_BDT.C`)

Reads slimtrees, loads TMVA models, computes BDT scores for all 11 model variants plus the NPB score, and writes them as new branches. Output is a clone of the input tree with score branches appended.

Key file: `apply_BDT.C` (492 lines)

### Stage 4: Efficiency and Yield (`efficiencytool/`)

Three sub-stages orchestrated by `oneforall.sh`:

1. **RecoEffCalculator_TTreeReader.C** -- Two-pass (vertex scan + full analysis). Reads scored slimtrees, applies all cuts, fills ABCD histograms, computes TEfficiency objects, builds response matrices. Runs per MC sample and on data.

2. **MergeSim.C** -- Merges per-sample outputs into combined signal (photon5+10+20), jet (jet5+8+12+20+30+40), and response matrix files. Cross-section weights were pre-applied during RecoEffCalculator.

3. **CalculatePhotonYield.C** -- Extracts purity via ABCD with signal leakage corrections, performs Bayesian unfolding (RooUnfoldBayes), applies efficiency corrections, computes final cross-section.

### Stage 5: Plotting and Systematics (`plotting/`)

`calc_syst_bdt.py` aggregates systematic variations in three tiers: per-type, per-group, total. `plot_final_selection.C` produces the final cross-section plot with JETPHOX NLO comparison and systematic bands.

### Orchestration

The `oneforall.sub` HTCondor submit file queues one job per `config_bdt_*.yaml`, enabling parallel systematic variation processing. The 3-step systematic pipeline is: (1) generate configs with `make_bdt_variations.py`, (2) submit with `condor_submit oneforall.sub`, (3) aggregate with `calc_syst_bdt.py`.

## Auxiliary Pipelines

These are significant analysis workflows that run alongside the main cross-section pipeline.

### Double-Interaction Blending (showershape, single-pass)

Full GEANT double-interaction MC samples blended with single-interaction MC at cluster-weighted fractions (22.4% for 0 mrad, 7.9% for 1.5 mrad crossing angle; computed by `calc_pileup_range.C`). Eight SI/DI pairs are available: photon5/10/20 and jet8/12/20/30/40 (see [Double-Interaction Efficiency](../concepts/double-interaction-efficiency.md) for the sample list).

The pipeline is now **single-pass with truth-vertex reweighting**: `ShowerShapeCheck.C` reads `h_w_iterative` from `efficiencytool/truth_vertex_reweight/output/{0mrad,1p5mrad}/reweight.root` via `TruthVertexReweightLoader.h` and applies `w(z_truth)` (single MC) or `w(z_hard)*w(z_mb)` (double MC) per event, multiplied by `mix_weight` (= cluster-weighted DI fraction or its SI complement). The previous two-pass reco-vertex approach has been retired (`run_showershape_double.sh` → `run_showershape_double_reco_legacy.sh`).

Condor packaging: `submit_showershape_di.sub` with row lists `showershape_di_jobs_{0rad,1p5rad}.list` (17 rows = 8 SI/DI pairs + 1 data per crossing-angle config); per-job driver `run_showershape_di_job.sh`.

Key files: `efficiencytool/submit_showershape_di.sub`, `efficiencytool/run_showershape_di_job.sh`, `efficiencytool/ShowerShapeCheck.C` (1917 lines), `efficiencytool/TruthVertexReweightLoader.h`, `efficiencytool/truth_vertex_reweight/config.yaml`.

### Shower Shape Analysis (`ShowerShapeCheck.C`)

Fills normalized distributions of 11 shower-shape variables (weta_cogx, wphi_cogx, wr_cogx, et1-4, e11/e33, e32/e35, bdt_score, npb_score) broken down by ABCD region, eta bin, and pT bin. Supports mix_weight for double-interaction blending.

Key files: `efficiencytool/ShowerShapeCheck.C`, `efficiencytool/run_showershape.sh`, `plotting/plot_showershapes_variations.C`

### Toy Double-Interaction Simulation (`DoubleInteractionCheck.C`)

Takes single-interaction MC events and artificially shifts the vertex to simulate pileup. Draws random second vertex from data distribution, recomputes cluster kinematics using CEMC geometry (R=93.5 cm), re-evaluates BDT and NPB scores, and re-classifies tight/non-tight and iso/non-iso. Includes Gaussian vertex smearing at 5/10/15/20 cm levels.

Key files: `efficiencytool/DoubleInteractionCheck.C` (1446 lines), `efficiencytool/run_double_interaction.sh`, `plotting/plot_double_interaction.C`

## Key Configurations

| Config | Location | Purpose |
|--------|----------|---------|
| `config.yaml` | `FunWithxgboost/` | BDT training parameters |
| `config_nom.yaml` | `FunWithxgboost/` | BDT application + feature extraction |
| `config_bdt_nom.yaml` | `efficiencytool/` | Nominal efficiency/yield analysis |
| `config_bdt_*.yaml` | `efficiencytool/` | Systematic variation configs (37 total) |
| `config_showershape_*.yaml` | `efficiencytool/` | Shower shape study configs |
| `config_npb_training.yaml` | `FunWithxgboost/` | NPB score training |
