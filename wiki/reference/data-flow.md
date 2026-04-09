# Data Flow Reference

## File Naming Conventions

### Stage 1: anatreemaker output

| Type | Pattern | Example |
|------|---------|---------|
| Data (per job) | `OUTTREE_DST_Jet_{run}_{segment}.root` | From condor |
| Sim (per job) | `caloana.root` in `OutDir_*/` | From condor |
| Sim (merged) | `condorout/combined.root` | Via hadd_combined.sh |

### Stage 2: BDT training files

| Type | Pattern | Location |
|------|---------|----------|
| MC features | `shapes_split_{sample}.txt` | `FunWithxgboost/` |
| Data NPB | `shapes_split_data_npb.txt` | `FunWithxgboost/` |
| Legacy | `NOSPLIT/shapes_{sample}.txt` | `FunWithxgboost/NOSPLIT/` |

### Stage 3: BDT-scored trees

| Type | Pattern | Location |
|------|---------|----------|
| Sim | `{sample}/bdt_split.root` | `FunWithxgboost/{sample}/` |
| Data | `part_{N}_with_bdt_split.root` | `anatreemaker/.../condorout/` |

### Stage 4: Efficiency tool output

All in `efficiencytool/results/`:

| Type | Pattern |
|------|---------|
| Per-sample MC | `MC_efficiency_{sample}_{var_type}.root` |
| Merged signal MC | `MC_efficiency_{var_type}.root` |
| Merged jet MC | `MC_efficiency_jet_{var_type}.root` |
| Response (per-sample) | `MC_response_{sample}_{var_type}.root` |
| Response (merged) | `MC_response_{var_type}.root` |
| Data | `data_histo_{var_type}.root` |
| Vtxscan (sim) | `MC_efficiency_{sample}_{var_type}_vtxscan.root` |
| Vtxscan (data) | `data_histo_{var_type}_vtxscan.root` |
| Final result | `Photon_final_{var_type}.root` |
| MC closure | `Photon_final_{var_type}_mc.root` |
| Shower shape (MC) | `MC_efficiencyshower_shape_{sample}_{var_type}.root` |
| Shower shape (data) | `data_histoshower_shape_{var_type}.root` |
| Double interaction | `{outfile}_double_interaction_check_{sample}.root` |

### Stage 5: Plotting output

| Type | Pattern | Location |
|------|---------|----------|
| Figures | `{plot_type}_{suffix}.pdf` | `plotting/figures/` |
| Shower shapes | `{prefix}_{var}_eta{i}_pt{j}_{cut}.pdf` | `plotting/figures/{suffix}/` |
| Systematics | `syst_bdt_{type}.root` | `plotting/rootFiles/` |
| Total syst | `syst_sum.root` | `plotting/rootFiles/` |

## Branch Contracts Between Stages

### anatreemaker -> BDTinput.C

Cluster branches with `_{CLUSTERINFO_CEMC}` suffix. All shower shape raw variables (e11-e77, et1-4, weta_cogx, wphi_cogx, w32/w52/w72), plus kinematics (Et, Eta, Phi), truth info (truthtrkID, pid, photonclass), and isolation (iso_02/03/04, iso_03_60_*).

### anatreemaker -> apply_BDT.C

Same branches as BDTinput.C. Additionally reads energy windows for ratio computation (e13, e15, e17, e31, e51, e71, e37, e53).

### apply_BDT.C -> RecoEffCalculator

All original slimtree branches PLUS the new BDT score branches:
- `cluster_bdt_{node}_{model_name}[ncluster]` -- 11 model scores
- `cluster_npb_score_{node}[ncluster]` -- NPB score

The model selected at runtime is determined by `bdt_et_bin_edges` + `bdt_et_bin_models` in the config.

### RecoEffCalculator -> MergeSim

Per-sample ROOT files containing:
- ABCD yield histograms (`h_tight_iso_cluster_0`, etc.)
- Signal-matched histograms (`_signal_0`)
- Background histograms (`_notmatch_0`)
- TEfficiency objects (`eff_reco_eta_0`, etc.)
- Truth pT histograms
- Response matrices (RooUnfoldResponse)

MergeSim adds these via TFileMerger -- cross-section weights are pre-applied.

### MergeSim -> CalculatePhotonYield

Merged files consumed:
- Signal MC: `MC_efficiency_{var_type}.root` (leakage fractions, efficiencies)
- Data: `data_histo_{var_type}.root` (ABCD yields)
- Response: `MC_response_{var_type}.root` (unfolding matrix)
- Jet MC (closure): `MC_efficiency_jet_{var_type}.root`

### CalculatePhotonYield -> Plotting

`Photon_final_{var_type}.root` containing:
- `h_unfold_sub_result` -- THE primary result
- Purity: `gpurity`, `gpurity_leak`, fit functions
- R factors: `h_R`, `h_R_notmatch`
- Leakage: `h_leak_B/C/D`
- All iterations: `h_unfold_sub_{1..10}`

### calc_syst_bdt.py -> plot_final_selection.C

`syst_sum.root` containing: `h_sum_low`, `h_sum_high`, `h_sum_rel_low`, `h_sum_rel_high`

## Config Archival

`SaveYamlToRoot()` pattern: the full YAML config is saved as a TNamed object in output ROOT files for reproducibility. This means any output file can be traced back to its exact configuration.

## Data Paths

| Purpose | Path |
|---------|------|
| Sim BDT-scored trees | `/sphenix/user/shuhangli/ppg12/FunWithxgboost/{sample}/bdt_split.root` |
| Data BDT-scored trees | `.../anatreemaker/macro_maketree/data/ana521/condorout/part_*_with_bdt_split.root` |
| Efficiency results | `/sphenix/user/shuhangli/ppg12/efficiencytool/results/` |
| BDT models | `/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/` |
| NPB model | `/sphenix/user/shuhangli/ppg12/FunWithxgboost/npb_models/npb_score_tmva.root` |
