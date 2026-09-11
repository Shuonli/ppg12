# PPG12: isolated prompt photons in p+p at √s = 200 GeV (sPHENIX)

Analysis code for the sPHENIX PPG12 measurement of the differential cross section of isolated prompt photons in p+p collisions at √s = 200 GeV with the Run 24 data (64.37 pb⁻¹). The result is reported for photon E_T of 12–32 GeV within |η| < 0.7, with truth-level isolation E_T^iso < 4 GeV in an R = 0.3 cone, and is compared with PYTHIA 8, NLO (JETPHOX, Vogelsang) and NNLO (NNLOJET) calculations. The analysis note is the `PPG12-analysis-note` submodule.

## Contents
- [Setup](#setup)
- [Repository layout](#repository-layout)
- [Pipeline](#pipeline)
- [Data, simulation and selection](#data-simulation-and-selection)
- [Signal extraction and corrections](#signal-extraction-and-corrections)
- [Configuration](#configuration)
- [Systematic uncertainties](#systematic-uncertainties)
- [Double-interaction studies](#double-interaction-studies)
- [Outputs](#outputs)
- [Documentation](#documentation)

## Setup

**Step 0, environment.** Every interactive shell and every condor job wrapper must first run

```bash
source /sphenix/user/shuhangli/ppg12/env.sh
```

`env.sh` at the repo root is the tracked equivalent of `/sphenix/u/shuhang98/setup.sh`, which 70 tracked `*.sh` scripts still source directly. It sources the sPHENIX `new` release (`sphenix_setup.sh -n new`: gcc 14.2.0, ROOT 6.32.06, cmake 3.31.0), puts `$MYINSTALL` on `LD_LIBRARY_PATH` and `ROOT_INCLUDE_PATH` through `setup_local.sh`, activates the analysis python venv, and exports `PPG12_ROOT`. Export `MYINSTALL`, `PPG12_SPHENIX_RELEASE` or `PPG12_VENV` before sourcing to override the defaults.

The ROOT macros read their YAML configs with yaml-cpp. Build it once into `$MYINSTALL` following [`docs/BUILD_yaml-cpp.md`](docs/BUILD_yaml-cpp.md) (yaml-cpp 0.8.0 at commit `73ef006`, shared library, `lib64/` layout). Caveat: 59 tracked macros, including the whole nominal chain (`RecoEffCalculator_TTreeReader.C`, `MergeSim.C`, `merge_periods.C`, `CalculatePhotonYield.C`, `plot_final_selection.C`, `BDTinput.C`), still call `gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so")` with the literal path. Until those lines are migrated, that exact path must exist, whatever `MYINSTALL` is set to.

Other requirements: RooUnfold for the unfolding, HTCondor for production runs, and a Python 3 environment with numpy, pandas, scipy, matplotlib, seaborn, xgboost, scikit-learn, optuna, joblib, uproot, awkward, PyYAML, ruamel.yaml and PyMuPDF.

## Repository layout

| Path | Contents |
|---|---|
| `anatreemaker/` | `source/CaloAna24.cc` Fun4All module (DST → `slimtree`, cluster node `CLUSTERINFO_CEMC`) and `macro_maketree/` production macros and condor scripts for data and simulation |
| `FunWithxgboost/` | Photon-ID BDT (`BDTinput.C`, `main_training.py`, `config.yaml`), non-physics-background (NPB) BDT (`train_npb_score.py`), scoring (`apply_BDT.C`), TMVA exports in `binned_models/` and `npb_models/` |
| `efficiencytool/` | Efficiency, purity, unfolding and cross section (`RecoEffCalculator_TTreeReader.C`, `MergeSim.C`, `merge_periods.sh`, `CalculatePhotonYield.C`), the `oneforall*` drivers and condor files, `make_bdt_variations.py`, the `config_bdt_*.yaml` configs, `CrossSectionWeights.h`, tower masks, truth-vertex reweighting and double-interaction studies |
| `plotting/` | Final and QA figures (`plot_final_selection.C`), systematic aggregation (`calc_syst_bdt.py`), paper-figure macros (`paper/`), NNLOJET tables (`theory/`), `rootFiles/syst_sum.root`, shared style and binning (`plotcommon.h`) |
| `NLO/` | JETPHOX histogram production and theory corrections (see `NLO/README.md`) |
| `sPhenix-sc{05,1,2}.dat` | NLO tables read by the final-figure macros |
| `lumi/` | Luminosity tables and cross-checks |
| `showershapecheck/`, `simcrosssection/` | Shower-shape validation and generator cross-section studies |
| `scripts/` | Results-site deployment (`deploy_pages.py`, `site_config.yaml`) |
| `wiki/` | Pipeline, concept and physics documentation, starting at `wiki/index.md` |
| `docs/` | Build notes |
| `PPG12-analysis-note/`, `ppg12_conf_note/` | Overleaf submodules for the analysis note and the conference note |

Study reports, condor job lists, notebooks and local assistant files are kept out of git on purpose (see `.gitignore`).

## Pipeline

1. **Trees.** `anatreemaker` turns DSTs into `slimtree` files with cluster, isolation, truth and MBD information.
2. **BDT training.** Photon ID: `root -l -b -q BDTinput.C`, then `python main_training.py --config config.yaml` (25 shower-shape and kinematic features, models per E_T bin, exported to TMVA). NPB: `python train_npb_score.py`.
3. **Scoring.** `root -l -b -q 'apply_BDT.C("config_nom.yaml", "<sample>", "<input file>")'` writes the photon-ID and NPB scores into the trees.
4. **Per-period efficiency and data histograms (condor phase 1).** For each crossing-angle period, `bash oneforall_tree_double_dispatch.sh config_bdt_nom_0rad.yaml` (and `_1p5mrad`) runs `RecoEffCalculator_TTreeReader.C` on the single- and double-interaction MC and on data, then `MergeSim.C`.
5. **All-range result (condor phase 2).** `bash oneforall.sh config_bdt_nom.yaml` merges the two periods' MC (`merge_periods.sh`) and data histograms (`hadd`), then runs `CalculatePhotonYield.C` (purity, unfolding, efficiency and luminosity corrections), writing `efficiencytool/results/Photon_final_bdt_nom.root`.
6. **Systematics.** See [Systematic uncertainties](#systematic-uncertainties).
7. **Figures.** In `plotting/`: `root -l -b -q 'plot_final_selection.C("bdt_nom")'` writes `figures/final_bdt_nom.pdf`. `bash paper/make_paper_figures.sh` writes the paper figures into a separate `PPG12-Paper/` checkout, which is not part of this repository.
8. **Note.** `cd PPG12-analysis-note && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex`

## Data, simulation and selection

**Run periods.** The analysis combines two crossing-angle periods. MC histograms are filled per period with the weight L_i / ΣL and summed.

| Period | Runs | Luminosity (beam-delivered) | Double-interaction fraction (cluster-weighted) |
|---|---|---|---|
| 0 mrad | 47289–51273 | 47.2076 pb⁻¹ | 0.224 |
| 1.5 mrad | 51274–54000 | 17.1642 pb⁻¹ | 0.079 |
| All | 47289–54000 | 64.3718 pb⁻¹ | |

**Simulation.** PYTHIA 8 prompt-photon samples `photon5`, `photon10`, `photon20` (truth p_T 0–14, 14–30, 30+ GeV) and inclusive jet samples `jet8`, `jet12`, `jet20`, `jet30`, `jet40`, each with a full-GEANT double-interaction partner (`*_double`). Samples are combined with the cross-section weights in `efficiencytool/CrossSectionWeights.h`, blended at the period's double-interaction fraction and reweighted to the data vertex distribution per period (`efficiencytool/truth_vertex_reweight/`).

**Selection** (`efficiencytool/config_bdt_nom.yaml`):
- |η| < 0.7 and MBD vertex |z| < 60 cm
- clusters whose centre tower lies in the φ-symmetry tower mask `mask_phisymm_tight` are removed in data and MC
- pre-selection (`analysis.common`): NPB score > 0.5, E11/E33 < 0.98, et1 > 0.6, E32/E35 < 1.0
- tight identification: BDT score > 0.815625 − 0.0015625·E_T, with the E_T-binned `base_v3E` models (bin edges 8, 15, 35 GeV)
- non-tight identification: an E_T-dependent BDT-score window below the tight threshold (`analysis.non_tight`)
- reco isolation: topo-cluster cone R = 0.4, isolated if E_T^iso < 0.490 + 0.037·E_T GeV, with the non-isolated sideband starting 0.8 GeV above the cut
- truth fiducial definition: direct or fragmentation photon with E_T^iso < 4 GeV in R = 0.3

Reco E_T bins are 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36 GeV (the same array as `ptRanges` in `plotting/plotcommon.h`). Truth bins extend from 8 to 45 GeV for the unfolding.

## Signal extraction and corrections

- **Purity.** Data-driven double sideband (ABCD) method with signal-leakage corrections from MC:

  ```
              Isolated     Non-isolated
  Tight       A (signal)   B
  Non-tight   C            D
  ```

  The purity is fit as a function of E_T (nominal [1/1] Padé form) and applied bin by bin.
- **Trigger.** Data clusters are weighted by the inverse L1 trigger efficiency.
- **Unfolding.** D'Agostini Bayesian unfolding with RooUnfold, 2 iterations.
- **Efficiencies.** Reconstruction, identification and isolation efficiencies from MC as a function of truth p_T, plus the MBD trigger and event-selection efficiency.
- **Cross section.** Scaled by the luminosity. The statistical uncertainty of the ABCD purity is propagated with 20,000 toys (Poisson re-throws of the region yields, Gaussian re-throws of the leakage corrections).

## Configuration

- `config_bdt_<variant>.yaml` is the all-range target of a variant. `config_bdt_<variant>_0rad.yaml` and `_1p5mrad.yaml` are its per-period merge feeders, which set the run range, the period luminosity with `lumi_target: 64.3718`, and the truth-vertex reweight file. The double-interaction fraction follows from `run_min`, unless `analysis.double_frac_override` is set.
- `var_type` becomes the output file suffix (`Photon_final_<var_type>.root`) and must be unique.
- `python make_bdt_variations.py config_bdt_nom.yaml` regenerates every variant config from the nominal one.

## Systematic uncertainties

Full procedure in [`SYST_PIPELINE_README.md`](SYST_PIPELINE_README.md). In short:

```bash
cd efficiencytool
python make_bdt_variations.py config_bdt_nom.yaml   # 1. variant configs
condor_submit oneforall_tree_double.sub             # 2a. per-period feeders (all config_bdt*.yaml)
condor_submit oneforall_all.sub                     # 2b. all-range oneforall.sh jobs, after 2a drains
cd ../plotting
python calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures   # 3. aggregate
```

`oneforall_all.sub` queues `all_range_configs.list`, which is not tracked: create it locally with one bare config name per line. Condor only accepts submit files whose event `log` is under `/tmp`. `unfold_iter3/4` come from `efficiencytool/postprocess_unfold_iter_scan.py` and `iso_generator` from `plotting/build_iso_generator_variant.C`, without condor jobs.

| Group | Variations |
|---|---|
| Energy scale | ±1.48% cluster energy scale, linear non-linearity (0 at 10 GeV, +1% at 30 GeV) |
| Energy resolution | larger of no extra smearing and a wider data resolution |
| Purity | tight BDT threshold +0.05, non-tight BDT threshold −0.10, non-isolated sideband start (0.1 / 1.0 GeV above the cut), Erf instead of Padé purity fit, purity-fit 68% band, MC purity-closure correction |
| Efficiency | MC isolation pedestal shift off, HERWIG instead of PYTHIA isolation efficiency |
| Unfolding | unweighted response prior, 3 and 4 iterations |
| Double interaction | blending fraction at the data χ² fit values (0.185 at 0 mrad, 0.000 at 1.5 mrad), both periods in one variation |
| NPB | NPB score cut 0.3 and 0.7 |
| Flat | luminosity +9.13% / −6.75%, MBD vertex efficiency ±6% |

Each variation is symmetrized or combined according to its role, the sources are added in quadrature within a group, and the total adds the groups and the flat terms per bin into `plotting/rootFiles/syst_sum.root`.

## Double-interaction studies

Shower-shape comparisons of single- and double-interaction MC use `efficiencytool/ShowerShapeCheck.C`, submitted with `submit_showershape_di.sub` (one row per sample with its mixing weight) and merged with `hadd_showershape_di.sh` and `merge_periods_showershape.sh`. The data-driven scan of the blending fraction is `plotting/chi2_weight_scan.C`.

## Outputs

| Location | Files |
|---|---|
| `efficiencytool/results/` | `MC_efficiency_*`, `MC_response_*`, `data_histo_<var_type>.root`, `Photon_final_<var_type>.root` (with `_mc` closure companions and the `_0rad` / `_1p5mrad` feeders) |
| `plotting/rootFiles/` | `syst_bdt_<type>.root`, `syst_bdt_group_<group>.root`, `syst_bdt_total.root`, `syst_sum.root` |
| `plotting/figures/` | `final_<var_type>.pdf`, `syst_bdt_rel_<type>.pdf`, `syst_bdt_breakdown.pdf` |

Only `plotting/rootFiles/syst_sum.root` is tracked. The rest is regenerated.

## Documentation

- [`SYST_PIPELINE_README.md`](SYST_PIPELINE_README.md): systematic-variation workflow
- [`wiki/index.md`](wiki/index.md): pipeline stages, concepts, configuration schema and physics background
- [`NLO/README.md`](NLO/README.md): JETPHOX production
- [`FunWithxgboost/README.md`](FunWithxgboost/README.md): BDT training details
- [`docs/BUILD_yaml-cpp.md`](docs/BUILD_yaml-cpp.md): yaml-cpp build
- Results site: https://shuonli.github.io/ppg12/, deployed with `python3 scripts/deploy_pages.py`

## License

Internal sPHENIX collaboration code.
