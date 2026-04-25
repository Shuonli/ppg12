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
# Build the library
cd /sphenix/user/shuhangli/ppg12/anatreemaker/source
./autogen.sh --prefix=$MYINSTALL
cd build && make install

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

# Extract features from slimtrees (one per sample)
root -l -b -q 'BDTinput.C("config_nom.yaml", "photon5")'
root -l -b -q 'BDTinput.C("config_nom.yaml", "photon10")'
root -l -b -q 'BDTinput.C("config_nom.yaml", "photon20")'
root -l -b -q 'BDTinput.C("config_nom.yaml", "jet10")'
root -l -b -q 'BDTinput.C("config_nom.yaml", "jet15")'
root -l -b -q 'BDTinput.C("config_nom.yaml", "jet20")'
root -l -b -q 'BDTinput.C("config_nom.yaml", "jet30")'
root -l -b -q 'BDTinput.C("config_nom.yaml", "jet50")'

# Train BDT (produces models in binned_models/)
python main_training.py --config config.yaml

# Train NPB score model
python train_npb_score.py --config config_npb_training.yaml
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

# Full per-config pipeline: 2-pass TTree reader (16 SI/DI samples) + hadd +
# MergeSim + CalculatePhotonYield x 2. Two ways to invoke:

# (a) Auto-derive DOUBLE_FRAC from run range via calc_pileup_range.C (~30s extra)
bash oneforall_tree_double_auto.sh 1.5mrad config_bdt_nom.yaml
bash oneforall_tree_double_auto.sh 0mrad   config_bdt_0rad.yaml

# (b) Pass DOUBLE_FRAC explicitly (faster):
bash oneforall_tree_double.sh config_bdt_nom.yaml  0.079
bash oneforall_tree_double.sh config_bdt_0rad.yaml 0.224
```

The legacy `oneforall_tree.sh` (single-interaction only, no DI mix) is kept for backward compat but should NOT be used for the production cross-section.

This produces (for config_bdt_nom.yaml):
- Per-sample MC: `results/MC_efficiency_{photon5,photon10,photon20,jet8,jet12,jet20,jet30,jet40}_bdt_nom.root`
- Merged: `results/MC_efficiency_bdt_nom.root`, `MC_efficiency_jet_bdt_nom.root`, `MC_response_bdt_nom.root`
- Data: `results/data_histo_bdt_nom.root`
- Final: `results/Photon_final_bdt_nom.root` (cross-section), `Photon_final_bdt_nom_mc.root` (MC closure)

### Full-run-range pipeline (0mrad + 1.5mrad combined)

The all-range nominal `config_bdt_all.yaml` reads a **lumi-weighted hadd of two per-period merge-feeder MC outputs**. Per-event `lumi_weight = lumi/lumi_target` in `RecoEffCalculator_TTreeReader.C` pre-scales the per-period MC; plain hadd via `merge_periods.sh` then reproduces the all-range expectation.

```bash
# 1. Per-period merge-feeders (each is a full Stage 4 run):
bash oneforall_tree_double.sh config_bdt_0rad.yaml 0.224
bash oneforall_tree_double.sh config_bdt_nom.yaml  0.079

# 2. Cross-period merge + CalculatePhotonYield (auto-dispatched by oneforall.sh
#    when var_type ends in _all):
bash oneforall.sh config_bdt_all.yaml
```

See `wiki/pipeline/full-run-range.md` for the merge mechanics, `lumi_target` semantics, and `oneforall.sh` dispatch rules.

## Stage 5a: Systematic Variations (full per-sample re-run on HTCondor)

When BDT cuts, isolation parameters, or any per-sample MC content changes, the variants must be **fully re-run** through the TTree-reader stage — `MergeSim`/`oneforall.sub` alone is NOT sufficient.

The production end-to-end run is two condor submissions, each gated on the previous:

### Phase 1 — Per-config TTree → MergeSim → CalculatePhotonYield

```bash
cd /sphenix/user/shuhangli/ppg12/efficiencytool

# Generate variation configs (87 configs as of April 2026)
python make_bdt_variations.py config_bdt_nom.yaml

# Submit Phase 1: one job per config_bdt_*.yaml.
# oneforall_tree_double_dispatch.sh derives DOUBLE_FRAC from each config's
# analysis.run_min (0.224 for run_min<51274, 0.079 otherwise) and invokes
# oneforall_tree_double.sh. All-range configs (var_type ending in _all) exit
# cleanly via the dispatch — they have no per-sample MC; handled in Phase 2.
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

Only the 2 `*_all.yaml` configs (`config_bdt_all.yaml`, `config_bdt_allz_all.yaml`) need this stage. They consume Phase 1's per-period merge-feeder MC (`MC_efficiency_bdt_{0rad,nom}.root` for `bdt_all`; `MC_efficiency_bdt_allz_{0rad,1p5mrad}.root` for `bdt_allz_all`).

```bash
condor_submit oneforall_all.sub  # queues exactly 2 jobs (the *_all configs)
```

`oneforall_all.sub` requests **4 GB per job** (CalculatePhotonYield is single-threaded and modest). Each job runs `oneforall.sh` which dispatches via filename: `*_all.yaml` → `merge_periods.sh` (lumi-weighted plain hadd of per-period merge-feeder outputs) → `RecoEffCalculator_TTreeReader.C(config, "data")` (produces the all-range `data_histo_{var_type}.root` since `merge_periods.sh` handles only MC) → `CalculatePhotonYield × 2`.

Outputs: `Photon_final_bdt_all.root` + `Photon_final_bdt_all_mc.root` and the analogous `bdt_allz_all` pair.

### Stage 5a-bis — Aggregate systematics

```bash
cd /sphenix/user/shuhangli/ppg12/plotting
python calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures
```

### Notes on which submit file to use

- `oneforall_tree_double.sub` (Phase 1, 16 GB) — full per-config re-run from TTree. Use whenever cuts or per-sample MC content changed.
- `oneforall_all.sub` (Phase 2, 4 GB, 2 jobs) — all-range targets only.
- `oneforall.sub` (legacy, 2 GB, 88 jobs) — runs `oneforall.sh` for every config_bdt_*.yaml. Only correct when per-sample MC is already current; otherwise produces stale-cut results. Mostly superseded by the Phase 1 + Phase 2 split.
- `oneforall_tree.sub` (legacy, 6 GB) — single-interaction-only TTree reader. Do not use for production.

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
python3 make_comparison_report.py --pair bdt_nom bdt_tightbdt50
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
ls Photon_final_bdt_*.root | wc -l    # Should be ~87 (April 2026: 87 VARIANTS in make_bdt_variations.py)
```

For the full-run-range pipeline, also verify:
```bash
ls Photon_final_bdt_{0rad,nom,all}.root           # 60cm fiducial set
ls Photon_final_bdt_allz_{0rad,1p5mrad,all}.root  # all-z lumi cross-check set
```
