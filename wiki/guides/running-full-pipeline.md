# Running the Full Pipeline

Step-by-step commands to run the PPG12 analysis from scratch.

## Prerequisites

- sPHENIX software environment sourced
- `libCaloAna24.so` built (see Stage 1)
- yaml-cpp available at `/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so`
- RooUnfold available

## Stage 1: Tree Making

Only needed if you need to reprocess DSTs. Usually the slimtrees already exist.

```bash
# Build the library (out-of-source: autogen.sh runs configure in the
# directory it is called from, so call it from build/)
cd /sphenix/user/shuhangli/ppg12/anatreemaker/source
mkdir -p build && cd build
../autogen.sh --prefix=$MYINSTALL
make install

# Run on data (edit runList.txt first)
cd /sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521
./run.sh

# Run on sim (all samples)
cd /sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28
./cleanup_and_run.sh

# After condor jobs complete, merge sim outputs
./hadd_combined.sh
```

## Stage 2: BDT Training

Only needed if retraining models. Usually the trained models already exist in `binned_models/`.

```bash
cd /sphenix/user/shuhangli/ppg12/FunWithxgboost

# Extract features from slimtrees, one text file per sample. The in-force photon-ID
# model (base_v3E, suffix _split via input.cluster_variants in config_nom.yaml) was
# trained on the run28 set below. jet10/jet15 have no usable run28 trees (their
# condorout/combined.root files are empty stubs) and jet50 is not part of the
# training set. config_nom.yaml (CLUSTERINFO_CEMC) writes shapes_split_<sample>.txt,
# config_nom_nosplit.yaml writes shapes_nosplit_<sample>.txt.
for s in photon5 photon10 photon20 jet5 jet12 jet20 jet30 jet40 data; do
  root -l -b -q "BDTinput.C(\"config_nom.yaml\", \"$s\")"
done
# (condor alternative: condor_submit submit_bdtinput.sub, which queues the same list
#  from bdtinput_joblist.txt for config_nom.yaml and config_nom_nosplit.yaml)

# Generate the per-model training configs variant_configs/config_<model>_split.yaml
# (file set = shapes_split_{photon5,photon10,photon20,jet5,jet12,jet20,jet30,jet40}.txt)
python make_variant_configs.py --variant split

# Train the nominal photon-ID BDT. Writes binned_models/model_base_v3E_split_single_tmva.root,
# the file apply_BDT.C loads for the CLUSTERINFO_CEMC entry of input.cluster_variants
# (model_suffix "_split", the legacy analysis.use_split_bdt flag is ignored when that
# block is present). Condor alternative for all 11 model variants:
# condor_submit submit_variant_training.sub
python main_training.py --config variant_configs/config_base_v3E_split.yaml

# Train NPB score model (writes npb_models/npb_score_split_tmva.root)
python train_npb_score.py --config config_npb_training_split.yaml
```

## Stage 3: Apply BDT Scores

```bash
cd /sphenix/user/shuhangli/ppg12/FunWithxgboost

# Apply to all sim samples
for sample in photon5 photon10 photon20 jet5 jet8 jet12 jet20 jet30 jet40; do
  root -l -b -q "apply_BDT.C(\"config_nom.yaml\", \"$sample\")"
done

# Apply to data (uses glob from config)
root -l -b -q 'apply_BDT.C("config_nom.yaml", "data")'
```

## Stage 4: Efficiency and Yield (Nominal)

The current production pipeline uses **double-interaction (DI) mixing** (cluster-weighted mix_weight per crossing angle: 0.224 for 0mrad, 0.079 for 1.5mrad). The DI orchestrator runs the full per-config pipeline in one shell call.

```bash
cd /sphenix/user/shuhangli/ppg12/efficiencytool

# The nominal is a triplet: config_bdt_nom.yaml (all-range target, lumi 64.3718)
# plus the two merge-feeders config_bdt_nom_0rad.yaml (runs 47289–51273,
# lumi 47.2076) and config_bdt_nom_1p5mrad.yaml (runs 51274–54000, lumi 17.1642),
# both with lumi_target 64.3718 so their MC hadd reproduces the all-range MC.

# Step 1 — per-period tree stage (RecoEffCalculator_TTreeReader on the 8 SI/DI
# sample pairs + data, then MergeSim + CalculatePhotonYield for the period). The
# dispatcher derives DOUBLE_FRAC from analysis.run_min (0.224 for 0 mrad, 0.079
# for 1.5 mrad):
bash oneforall_tree_double_dispatch.sh config_bdt_nom_0rad.yaml
bash oneforall_tree_double_dispatch.sh config_bdt_nom_1p5mrad.yaml

# Equivalent explicit form (DOUBLE_FRAC as second argument):
#   bash oneforall_tree_double.sh config_bdt_nom_0rad.yaml    0.224
#   bash oneforall_tree_double.sh config_bdt_nom_1p5mrad.yaml 0.079
# or derive it from the run range with calc_pileup_range.C (~30 s extra):
#   bash oneforall_tree_double_auto.sh 0mrad   config_bdt_nom_0rad.yaml
#   bash oneforall_tree_double_auto.sh 1.5mrad config_bdt_nom_1p5mrad.yaml

# Step 2 — all-range merge + yield. oneforall.sh sees both feeders on disk and runs
# merge_periods.sh (hadd of the feeder MC outputs), hadd of the two per-period
# data_histo files, CalculatePhotonYield x 2 and the selection plots.
bash oneforall.sh config_bdt_nom.yaml
```

The legacy `oneforall_tree.sh` (single-interaction only, no DI mix) is kept for backward compat but should NOT be used for the production cross-section.

This produces:
- Per-period MC (feeders): `results/MC_efficiency_{photon5,photon10,photon20,jet8,jet12,jet20,jet30,jet40}_bdt_nom_{0rad,1p5mrad}.root`, merged into `results/MC_efficiency_bdt_nom_{0rad,1p5mrad}.root`
- Per-period data and yields: `results/data_histo_bdt_nom_{0rad,1p5mrad}.root`, `results/Photon_final_bdt_nom_{0rad,1p5mrad}.root`
- All-range merged MC: `results/MC_efficiency_bdt_nom.root`, `MC_efficiency_jet_bdt_nom.root`, `MC_response_bdt_nom.root`
- All-range data: `results/data_histo_bdt_nom.root`
- Final: `results/Photon_final_bdt_nom.root` (cross-section), `Photon_final_bdt_nom_mc.root` (MC closure)

### Full-run-range mechanics (0mrad + 1.5mrad combined)

The all-range nominal `config_bdt_nom.yaml` reads a **lumi-weighted hadd of the two per-period merge-feeder MC outputs**. Per-event `lumi_weight = lumi/lumi_target` in `RecoEffCalculator_TTreeReader.C` pre-scales the per-period MC; plain hadd via `merge_periods.sh` then reproduces the all-range expectation. `oneforall.sh` dispatches on the filename: `*_0rad.yaml` / `*_1p5mrad.yaml` → `MergeSim.C` then `CalculatePhotonYield × 2` (no selection plots), bare name with both companions on disk → `merge_periods.sh` + data hadd + `CalculatePhotonYield × 2` + selection plots, bare name without companions → standalone `MergeSim.C` + `CalculatePhotonYield × 2` + selection plots. The same triplet layout applies to every generated variant (`config_bdt_<name>{,_0rad,_1p5mrad}.yaml`).

See `wiki/pipeline/full-run-range.md` for the merge mechanics, `lumi_target` semantics, and `oneforall.sh` dispatch rules.

## Stage 5a: Systematic Variations (full per-sample re-run on HTCondor)

When BDT cuts, isolation parameters, or any per-sample MC content changes, the variants must be **fully re-run** through the TTree-reader stage — `MergeSim`/`oneforall.sub` alone is NOT sufficient.

The production end-to-end run is two condor submissions, each gated on the previous:

### Phase 1 — Per-config TTree → MergeSim → CalculatePhotonYield

```bash
cd /sphenix/user/shuhangli/ppg12/efficiencytool

# Generate variation configs (175 configs from 63 variants as of 2026-09: one
# bare-name all-range config per variant plus _0rad/_1p5mrad merge-feeders for
# the 56 variants that do not pin a run range)
python make_bdt_variations.py config_bdt_nom.yaml

# Submit Phase 1: one job per config_bdt_*.yaml.
# oneforall_tree_double_dispatch.sh derives DOUBLE_FRAC from each config's
# analysis.run_min (0.224 for run_min<51274, 0.079 otherwise, or
# analysis.double_frac_override when set) and invokes oneforall_tree_double.sh.
# Bare-name configs whose _0rad/_1p5mrad companions exist on disk exit cleanly
# via the dispatch — they have no per-sample MC; handled in Phase 2.
# NOTE: the schedd rejects submit files whose event `log` is not under /tmp/.
# oneforall_tree_double.sub still has `log = logs/$(cfg).tree_double.log`, so set
# `log = /tmp/ppg12_tree_double_$(Cluster).log` first (oneforall_tree_double_rerun.sub
# already does this and queues from syst_configs_rerun.list).
condor_submit oneforall_tree_double.sub
```

`oneforall_tree_double.sub` requests **16 GB per job** because the DI pipeline runs 17 parallel ROOT processes per pass (8 SI + 8 DI MC + 1 data). With the older 6 GB request all jobs went on hold with `Job has gone over cgroup memory limit of 6144 megabytes` (peak ~6.1 GB per process group). 16 GB gives comfortable headroom.

**Pass 1 auto-skip**: under `analysis.truth_vertex_reweight_on=1` (current production default), `oneforall_tree_double.sh` auto-skips Pass 1 (vtxscan) + the vtxscan hadd because the reco-vertex reweight branch is force-disabled downstream in `RecoEffCalculator_TTreeReader.C` and the vtxscan output is never read. This halves per-config wall-clock from ~45 min to ~22 min. Legacy configs with `truth_vertex_reweight_on=0` still run both passes.

Wait for completion (~30–60 min per job, run in parallel):
```bash
# Poll until queue empties
while [ "$(condor_q $USER -nobatch | grep '^Total for query:' | awk '{print $4}')" != "0" ]; do
  sleep 60
done
```

Per-config logs go to `efficiencytool/logs/{config}.tree_double.{log,out,err}`. Inspect `.err` files for any held / failed job.

### Phase 2 — Cross-period merge + CalculatePhotonYield for all-range configs

Every bare-name config that has `_0rad` / `_1p5mrad` companions needs this stage (the nominal `config_bdt_nom.yaml` and every variant that does not pin a run range). Each consumes Phase 1's per-period merge-feeder MC (`MC_efficiency_bdt_<name>_{0rad,1p5mrad}.root`) and data (`data_histo_bdt_<name>_{0rad,1p5mrad}.root`).

```bash
condor_submit oneforall_all.sub  # queues the bare-name configs listed in all_range_configs.list
                                 # (same /tmp log caveat as Phase 1: edit the log line first)
```

`oneforall_all.sub` requests **4 GB per job** (CalculatePhotonYield is single-threaded and modest). Each job runs `oneforall.sh` which dispatches via filename: bare name with both companions on disk → `merge_periods.sh` (lumi-weighted plain hadd of the per-period merge-feeder MC outputs) → `hadd` of the two per-period `data_histo` files into the all-range `data_histo_bdt_<name>.root` → `CalculatePhotonYield × 2` → `make_selection_plots.sh`.

Outputs: `Photon_final_bdt_<name>.root` + `Photon_final_bdt_<name>_mc.root` per bare-name config (`bdt_nom` for the nominal).

### Stage 5a-bis — Aggregate systematics

```bash
cd /sphenix/user/shuhangli/ppg12/plotting
python calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures
```

### Notes on which submit file to use

- `oneforall_tree_double.sub` (Phase 1, 16 GB, queue = `config_bdt*.yaml` glob) — full per-config re-run from TTree. Use whenever cuts or per-sample MC content changed.
- `oneforall_tree_double_rerun.sub` (Phase 1, 16 GB, queue = `syst_configs_rerun.list`: 62 entries, nominal, the systematic-budget variants with their feeders, and the `mc_purity_correction_pol2` cross-check triplet) — scoped copy whose `log` already sits under `/tmp/`, so it submits without editing.
- `oneforall_all.sub` (Phase 2, 4 GB, queue = `all_range_configs.list`) — all-range targets only.
- `oneforall.sub` (legacy, 2 GB, one job per `config_bdt*.yaml`) — runs `oneforall.sh` for every config. Only correct when per-sample MC is already current; otherwise produces stale-cut results. Mostly superseded by the Phase 1 + Phase 2 split.
- `oneforall_tree.sub` (legacy, 6 GB) — single-interaction-only TTree reader. Do not use for production.
- Except for `oneforall_tree_double_rerun.sub`, these files still carry `log = logs/...`, which the schedd now rejects. Change it to `log = /tmp/<name>_$(Cluster).log` (keep `output`/`error` under `logs/`) before `condor_submit`.

## Stage 5b: Final Plots

```bash
cd /sphenix/user/shuhangli/ppg12/plotting

# Final cross-section plot
root -l -b -q 'plot_final_selection.C("bdt_nom")'

# All selection plots for nominal
bash make_selection_plots.sh bdt_nom

# All selection plots for all configs
bash make_all_bdt_selection.sh
```

## Stage 5c: Reports

```bash
cd /sphenix/user/shuhangli/ppg12/plotting

# Selection report
python3 make_selection_report.py

# Shower shape report
python3 make_showershape_report.py

# Comparison report
python3 make_comparison_report.py --pair bdt_nom bdt_tightup_p05
```

## Stage 6: Analysis Note

```bash
cd /sphenix/user/shuhangli/ppg12/PPG12-analysis-note
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

## Quick Check: Is My Pipeline Run Complete?

After running the efficiency pipeline, verify these files exist:

```bash
cd /sphenix/user/shuhangli/ppg12/efficiencytool/results
ls Photon_final_bdt_nom.root          # Final result
ls MC_efficiency_bdt_nom.root          # Merged signal MC
ls MC_efficiency_jet_bdt_nom.root      # Merged jet MC
ls data_histo_bdt_nom.root             # Data histograms
ls MC_response_bdt_nom.root            # Response matrix
```

For systematic variations, check all variants:
```bash
ls Photon_final_bdt_*.root | grep -v _mc.root | wc -l   # one per config that has been run
                                                        # (175 generated configs as of 2026-09, plus hand-made ones)
```

For the full-run-range pipeline, also verify the nominal triplet:
```bash
ls Photon_final_bdt_nom_{0rad,1p5mrad}.root   # per-period merge-feeders
ls Photon_final_bdt_nom.root                  # all-range target (beam-delivered lumi 64.3718 pb^-1)
```
