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

```bash
cd /sphenix/user/shuhangli/ppg12/efficiencytool

# Two-pass efficiency calculation (all samples)
bash oneforall_tree.sh config_bdt_nom.yaml

# Merge and compute cross-section
bash oneforall.sh config_bdt_nom.yaml
```

This produces:
- `results/MC_efficiency_bdt_nom.root` (merged signal MC)
- `results/MC_efficiency_jet_bdt_nom.root` (merged jet MC)
- `results/data_histo_bdt_nom.root` (data)
- `results/Photon_final_bdt_nom.root` (final cross-section)
- `results/Photon_final_bdt_nom_mc.root` (MC closure)

## Stage 5a: Systematic Variations

```bash
cd /sphenix/user/shuhangli/ppg12/efficiencytool

# Generate variation configs
python make_bdt_variations.py config_bdt_nom.yaml

# Run all variations on HTCondor
condor_submit oneforall.sub

# Wait for all jobs to complete, then aggregate
cd /sphenix/user/shuhangli/ppg12/plotting
python calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures
```

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
ls Photon_final_bdt_*.root | wc -l    # Should be ~37
```
