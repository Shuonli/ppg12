# PPG12 End-to-End Data Flow Report

## 1. Pipeline Overview

```
DST files
  |
  |  [CaloAna24.cc via Fun4All macros]
  v
slimtree ROOT files  (tree: slimtree, per-run or per-job)
  |
  +-- [hadd] --> combined.root per sample
  |
  +-- [BDTinput.C] --> training text files (shapes_*.txt)
  |                        |
  |                        +-- [main_training.py] --> BDT models (binned_models/*.root)
  |                                                       |
  +-- [apply_BDT.C] <------------------------------------+
  |         |
  |         v
  |   scored slimtree ROOT files (*_with_bdt_split.root, bdt_split.root)
  |         |
  |         +-- [RecoEffCalculator_TTreeReader.C, two-pass]
  |                    |
  |                    v
  |              per-sample efficiency/ABCD ROOT files + response matrices
  |                    |
  |                    +-- [MergeSim.C]
  |                           |
  |                           v
  |                     merged MC efficiency + response ROOT files
  |                           |
  |                           +-- [CalculatePhotonYield.C]
  |                                      |
  |                                      v
  |                                Photon_final_*.root (cross-section)
  |                                      |
  |                 +--------------------+-------------------+
  |                 |                                        |
  |                 v                                        v
  |    [calc_syst_bdt.py]                     [plot_final_selection.C]
  |           |                                        |
  |           v                                        v
  |    syst_sum.root ---> [plot_final_selection.C] --> figures/*.pdf
  |                                                       |
  |                                                       v
  +--------------------------------------------> PPG12-analysis-note/Figures/
```

---

## 2. Stage-by-Stage File Handoffs

### Stage 1: DST --> slimtree (anatreemaker)

**Producer:** `anatreemaker/source/CaloAna24.cc` invoked via Fun4All macros

**Data macro:** `anatreemaker/macro_maketree/data/ana521/Fun4All_runDST.C`
- Input: DST files (`DST_Jet`, `DST_JETCALO`), production tag `ana521_2025p007_v001`
- Output: `OUTTREE_DST_Jet_*.root` per condor job, containing TTree `slimtree`

**Sim macro:** `anatreemaker/macro_maketree/sim/run28/{sample}/Fun4All_run_sim.C`
- Input: 4-5 DST streams (G4Hits, DST_CALO_CLUSTER, DST_MBD_EPD, DST_TRUTH_JET)
- Output: `caloana.root` per condor job, containing TTree `slimtree`
- Merged via `hadd_combined.sh` into `condorout/combined.root` per sample

**Samples produced:**
- Signal: `photon5/`, `photon10/`, `photon20/` (truth pT windows 0-14, 14-30, 30+ GeV)
- Background: `jet5/`, `jet8/`, `jet10/`, `jet12/`, `jet15/`, `jet20/`, `jet30/`, `jet40/`, `jet50/`
- Double-interaction: `photon10_double/`, `jet12_double/`
- Minimum bias: `mb/`

**Output tree name:** `slimtree` (hardcoded in CaloAna24.cc)

**Key output histograms (non-tree):**
- `sim_cross_counting` (TH1I): bin 1 = N_Generator_Accepted, bin 2 = N_Processed (MC only)
- `tracking_radiograph` (TH3F): conversion study vertex map (MC only)

---

### Stage 2a: slimtree --> training text files (BDTinput.C)

**Producer:** `FunWithxgboost/BDTinput.C`
**Config:** `FunWithxgboost/config_nom.yaml` (keys: `input.cluster_node_name`, `input.tree`, `input.photon_jet_file_root_dir`)

**Input:** slimtree ROOT files (pre-BDT-scoring, from anatreemaker `combined.root` or data files)
- Sim path: `{photon_jet_file_root_dir}{filetype}/condorout/combined.root`
  - Note: config_nom.yaml points to `anatreemaker/macro_maketree/sim/run28/`
- Data path: from config `input.data_file`

**Branches read from slimtree (all cluster branches suffixed `_{CLUSTERINFO_CEMC}`):**
- Event-level: `mbdnorthhit`, `mbdsouthhit`, `vertexz`, `energy_scale`, `runnumber`, `scaledtrigger[64]`, `livetrigger[64]`
- Truth (MC): `nparticles`, `particle_trkid[]`, `particle_photonclass[]`, `njet_truth`, `jet_truth_Pt[]`
- Per-cluster: `cluster_Et`, `cluster_Eta`, `cluster_Phi`, `cluster_prob`, `cluster_truthtrkID`
- Shower shape: `cluster_et1-4`, `cluster_weta_cogx`, `cluster_wphi_cogx`, `cluster_nsaturated`
- Energy windows: `cluster_e11`, `cluster_e22`, `cluster_e13`, `cluster_e15`, `cluster_e17`, `cluster_e31`, `cluster_e51`, `cluster_e71`, `cluster_e33`, `cluster_e35`, `cluster_e37`, `cluster_e53`, `cluster_e73`, `cluster_e77`
- Widths: `cluster_w32`, `cluster_w52`, `cluster_w72`
- Energy sums: `cluster_e32`
- Isolation: `cluster_iso_02`, `cluster_iso_03`, `cluster_iso_04`, `cluster_iso_03_60_{emcal,hcalin,hcalout}`
- NPB data: `mbd_time`, `njet`, `jet_Pt[]`, `jet_Phi[]`, `cluster_time_array[]`, `cluster_e_array[]`, `cluster_ownership_array[]`

**Derived features computed (not raw branches):**
13 energy ratios: `e11_over_e33`, `e32_over_e35`, `e11_over_e22`, `e11_over_e13`, `e11_over_e15`, `e11_over_e17`, `e11_over_e31`, `e11_over_e51`, `e11_over_e71`, `e22_over_e33`, `e22_over_e35`, `e22_over_e37`, `e22_over_e53`

**Output files:**
- MC: `shapes_split_{filetype}.txt` (e.g., `shapes_split_photon5.txt`)
- Data NPB: `shapes_split_data_npb.txt`
- Legacy: `NOSPLIT/shapes_{filetype}.txt`

**Output format:** Space-separated, 30 columns with header row:
```
cluster_Et cluster_Eta cluster_Phi vertexz
e11_over_e33 e32_over_e35 e11_over_e22 e11_over_e13 e11_over_e15 e11_over_e17
e11_over_e31 e11_over_e51 e11_over_e71 e22_over_e33 e22_over_e35 e22_over_e37 e22_over_e53
cluster_prob cluster_weta_cogx cluster_wphi_cogx
cluster_et1 cluster_et2 cluster_et3 cluster_et4
cluster_w32 cluster_w52 cluster_w72
recoisoET is_tight pid
```

Last 3 columns (`recoisoET`, `is_tight`, `pid`) are metadata -- not used as BDT training features.

---

### Stage 2b: training text files --> BDT models (main_training.py)

**Producer:** `FunWithxgboost/main_training.py` (BinnedTrainingPipeline class)
**Config:** `FunWithxgboost/config.yaml`

**Input files consumed:** Text files listed under `data.single_file_set` in config.yaml.
Currently uses `NOSPLIT/shapes_*.txt` (legacy naming).
- Signal: `NOSPLIT/shapes_photon{5,10,20}.txt`
- Background: `NOSPLIT/shapes_jet{10,15,20,30,50}.txt`

**Signal selection:** `pid in [1, 2]` (direct + fragmentation photon)
**Background selection:** `pid not in [1, 2]`

**25 training features** (from `config.yaml: data.features`, exact order matters for TMVA positional indexing):
```
cluster_Et, cluster_weta_cogx, cluster_wphi_cogx, vertexz, cluster_Eta,
e11_over_e33, cluster_et1, cluster_et2, cluster_et3, cluster_et4,
e32_over_e35, cluster_w32, cluster_w52, cluster_w72,
e11_over_e22, e11_over_e13, e11_over_e15, e11_over_e17,
e11_over_e31, e11_over_e51, e11_over_e71,
e22_over_e33, e22_over_e35, e22_over_e37, e22_over_e53
```

**Output models (11 variants x 2 naming patterns):**
- Standard: `binned_models/model_{variant}_single_tmva.root`
- Split: `binned_models/model_{variant}_split_single_tmva.root`

Variant names: `base`, `base_vr`, `base_v0`, `base_v1`, `base_v2`, `base_v3`, `base_E`, `base_v0E`, `base_v1E`, `base_v2E`, `base_v3E`

Each model uses a DIFFERENT subset of the 25 features (see Section 4 for the exact feature vectors). Feature order is positional (`f0, f1, f2, ...`) in the TMVA file.

**NPB model (separate training):**
- Producer: `FunWithxgboost/train_npb_score.py`
- Config: `FunWithxgboost/config_npb_training.yaml`
- Output: `npb_models/npb_score_tmva.root`
- Signal: jet MC clusters; Background: NPB-tagged data clusters

---

### Stage 3: slimtree + BDT models --> scored slimtree (apply_BDT.C)

**Producer:** `FunWithxgboost/apply_BDT.C`
**Config:** `FunWithxgboost/config_nom.yaml` (uses `input.cluster_node_name`, `analysis.use_split_bdt`)

**Input:** Pre-scoring slimtree ROOT files + trained TMVA models
- Sim: `{sample}/condorout/combined.root` (from anatreemaker)
- Data: individual part files from anatreemaker condor output
- Models: `binned_models/model_{name}[_split]_single_tmva.root` and `npb_models/npb_score_tmva.root`

**Branches read from slimtree (same cluster node suffix pattern):**
- `vertexz`, `ncluster_{node}`
- `cluster_Et`, `cluster_Eta`, `cluster_Phi`, `cluster_et1-4`, `cluster_weta_cogx`, `cluster_wphi_cogx`
- `cluster_e11`, `cluster_e22`, `cluster_e33`, `cluster_e35`, `cluster_e32`
- `cluster_e13`, `cluster_e15`, `cluster_e17`, `cluster_e31`, `cluster_e51`, `cluster_e71`
- `cluster_e37`, `cluster_e53`
- `cluster_prob`, `cluster_w32`, `cluster_w52`, `cluster_w72`
- `cluster_iso_02`, `cluster_iso_03`, `cluster_iso_04`, `cluster_iso_03_60_{emcal,hcalin,hcalout}`

**Derived features recomputed in apply_BDT.C:**
- `e11_over_e33_BDT = (cluster_e11 > 0) ? cluster_e11 / cluster_e33 : 0`
- `e32_over_e35_BDT = (cluster_e32 > 0) ? cluster_e32 / cluster_e35 : 0`

Note: Division guard checks `numerator > 0` (apply_BDT.C) vs `denominator != 0 && isfinite` (BDTinput.C). Functionally equivalent for positive energies but technically different.

**ET threshold for scoring:** `cluster_Et > 7 GeV` (else score = -1)

**Output branches added (per cluster array, for each of 11 models):**
```
cluster_bdt_{clusternodename}_{model_name}[ncluster_{clusternodename}]/F
```
Example: `cluster_bdt_CLUSTERINFO_CEMC_base_v3E`

**NPB score branch:**
```
cluster_npb_score_{clusternodename}[ncluster_{clusternodename}]/F
```

**Output files:**
- Sim: `{filetype}/bdt_split.root` (e.g., `FunWithxgboost/photon10/bdt_split.root`)
- Data: `{inputfilename_without_ext}_with_bdt_split.root`

The output is a full clone of the input slimtree with the BDT branches appended.

---

### Stage 4a: scored slimtree --> efficiency + ABCD histograms (RecoEffCalculator_TTreeReader.C)

**Producer:** `efficiencytool/RecoEffCalculator_TTreeReader.C`
**Config:** `efficiencytool/config_bdt_*.yaml`
**Orchestration:** `efficiencytool/oneforall_tree.sh` (two-pass: vertex scan then full analysis)

**Input files:**
- Sim: `{photon_jet_file_root_dir}{filetype}{photon_jet_file_branch_dir}`
  - With `config_bdt_nom.yaml`: `/sphenix/user/shuhangli/ppg12/FunWithxgboost/{filetype}/bdt_split.root`
  - Note: `photon_jet_file_root_dir` in `config_bdt_nom.yaml` points to `FunWithxgboost/` (the BDT-scored output), NOT to `anatreemaker/`
- Data: `{data_file}` glob pattern (e.g., `part_*_with_bdt_split.root`)

**Key difference from Stage 2a:** RecoEffCalculator reads the BDT-scored trees from FunWithxgboost, not the raw trees from anatreemaker.

**Branches consumed from scored slimtree:**

*Event-level:*
`mbdnorthhit`, `mbdsouthhit`, `vertexz`, `vertexz_truth`, `runnumber`, `scaledtrigger[64]`, `livetrigger[64]`, `energy_scale`, `pythiaid`, `mbd_time`, `mbdnorthtmean`, `mbdsouthtmean`, `mbdnorth{q,t}[64]`, `mbdsouth{q,t}[64]`, `trigger_prescale[64]`

*Truth particles (MC):*
`nparticles`, `particle_{E,Pt,Eta,Phi}[]`, `particle_{pid,trkid,photonclass,converted}[]`, `particle_truth_iso_{02,03,04}[]`

*Per-cluster (all suffixed `_{clusternodename}`):*
- Kinematics: `cluster_{E,Et,Eta,Phi}`
- ID variables: `cluster_prob`, `cluster_CNN_prob`, `cluster_{truthtrkID,pid}` (MC), `cluster_nsaturated`
- Shower shape: `cluster_{et1-4}`, `cluster_{weta,wphi}`, `cluster_{weta_cogx,wphi_cogx}`, `cluster_{ietacent,iphicent}`, `cluster_{detamax,dphimax}`
- Energy windows: `cluster_e{11,22,13,15,17,31,51,71,33,35,37,53,73,55,57,75,77}`, `cluster_e32`, `cluster_{e1-4}` (raw quadrant energies)
- Widths: `cluster_{w32,w52,w72}`, `cluster_{e52,e72}`
- Isolation: `cluster_iso_{02,03,04}`, `cluster_iso_topo_{03,04}`, `cluster_iso_03_{emcal,hcalin,hcalout}`, `cluster_iso_03_{60,70,120}_{emcal,hcalin,hcalout}`, `cluster_iso_{005,01,02}_70_emcal`
- HCal: `cluster_{ihcal,ohcal}_et`, `cluster_{ihcal,ohcal}_et{22,33}`, `cluster_{ihcal,ohcal}_{ieta,iphi}`
- Tower arrays: `cluster_{e,adc,time}_array`, `cluster_{status,ownership}_array`
- **BDT scores** (from apply_BDT.C): `cluster_bdt_{clusternodename}_{model_name}` -- dynamically loaded per config
- **NPB score** (from apply_BDT.C): `cluster_npb_score_{clusternodename}`
- Jets: `njet`, `jet_{E,Pt,Eta,Phi}`
- Truth jets: `njet_truth`, `jet_truth_{E,Pt,Eta,Phi}`

**BDT model selection at runtime:**
The config specifies `bdt_model_name` (fallback) and optionally `bdt_et_bin_edges` + `bdt_et_bin_models` for ET-dependent model selection. The code builds `TTreeReaderArray` objects for each unique model name, then selects per cluster based on ET.

**Double-interaction MC branch aliasing:**
For `filetype == "jet12_double"` or `"photon10_double"`, jet branches use the full container name suffix (e.g., `njet_truth_AntiKt_Truth_r04`, `jet_Pt_AntiKt_unsubtracted_r04_calib`). For standard samples, short names are used (`njet_truth`, `jet_Pt`).

**Two-pass architecture:**
1. Pass 1 (`do_vertex_scan=true`): fills only `h_vertexz`, writes to `{outfile}_vtxscan.root`
2. Pass 2 (`do_vertex_scan=false`): reads vtxscan files for vertex reweighting, runs full analysis

**Output files:**

Per-sample sim:
```
{eff_outfile}_{filetype}_{var_type}.root
```
Example: `results/MC_efficiency_photon10_bdt_nom.root`

Per-sample response matrix:
```
{response_outfile}_{filetype}_{var_type}.root
```
Example: `results/MC_response_photon10_bdt_nom.root`

Data:
```
{data_outfile}_{var_type}.root
```
Example: `results/data_histo_bdt_nom.root`

Vtxscan:
```
{eff_outfile}_{filetype}_{var_type}_vtxscan.root
{data_outfile}_{var_type}_vtxscan.root
```

**Key histogram names in output (per eta bin index):**
- ABCD yields: `h_tight_iso_cluster_0`, `h_tight_noniso_cluster_0`, `h_nontight_iso_cluster_0`, `h_nontight_noniso_cluster_0`
- Signal-matched: `h_tight_iso_cluster_signal_0`, `h_tight_noniso_cluster_signal_0`, etc.
- Background (not-matched): `h_tight_iso_cluster_notmatch_0`, etc.
- Vertex: `h_vertexz`
- Truth pT: `h_truth_pT_0`, `h_truth_pT_vertexcut_0`, `h_truth_pT_vertexcut_mbd_cut_0`
- Common cluster: `h_common_cluster_0`
- Isolation ET: `h_tight_isoET_0_{pTbin}`, `h_nontight_isoET_0_{pTbin}`
- TEfficiency objects: `eff_reco_eta_0`, `eff_iso_eta_0`, `eff_id_eta_0`, `eff_all_eta_0`
- MBD efficiency: `g_mbd_eff`, `g_mbd_eff_north`, `g_mbd_eff_south`
- Response matrix: `h_response_full_0`, `h_response_half_0`, `h_pT_truth_response_0`, `h_pT_reco_response_0`

---

### Stage 4b: per-sample --> merged MC files (MergeSim.C)

**Producer:** `efficiencytool/MergeSim.C`
**Config:** same `config_bdt_*.yaml`

**Input (photon signal, 3 files):**
```
{eff_outfile}_photon5_{var_type}.root
{eff_outfile}_photon10_{var_type}.root
{eff_outfile}_photon20_{var_type}.root
```

**Output (merged signal):**
```
{eff_outfile}_{var_type}.root
```
Example: `results/MC_efficiency_bdt_nom.root`

**Input (jet background, 6 files):**
```
{eff_outfile}_jet{5,8,12,20,30,40}_{var_type}.root
```

**Output (merged jet):**
```
{eff_outfile}_jet_{var_type}.root
```
Example: `results/MC_efficiency_jet_bdt_nom.root`

**Input (response matrix, 3 photon files):**
```
{response_outfile}_photon{5,10,20}_{var_type}.root
```

**Output (merged response):**
```
{response_outfile}_{var_type}.root
```
Example: `results/MC_response_bdt_nom.root`

Cross-section weights are NOT applied in MergeSim -- they are pre-applied in RecoEffCalculator via `PPG12::GetSampleConfig(filetype).weight`, making TFileMerger addition correct.

**Important:** MergeSim expects exactly 6 jet samples (jet5/8/12/20/30/40). It does NOT include jet10, jet15, or jet50. The `oneforall_tree.sh` orchestration matches this exactly.

---

### Stage 4c: merged files --> cross-section (CalculatePhotonYield.C)

**Producer:** `efficiencytool/CalculatePhotonYield.C`
**Config:** same `config_bdt_*.yaml`

**Input files consumed:**
1. Merged signal MC: `{eff_outfile}_{var_type}.root`
2. Data (or MC-as-data for closure): `{data_outfile}_{var_type}.root` (data) or `{eff_outfile}_jet_{var_type}.root` (MC closure)
3. Response matrix: `{response_outfile}_{var_type}.root`

**Key histograms read:**
- From signal MC: `h_tight_iso_cluster_signal_0`, `h_tight_noniso_cluster_signal_0`, `h_nontight_iso_cluster_signal_0`, `h_nontight_noniso_cluster_signal_0` (for leakage), TEfficiency objects `eff_reco_eta_0`, `eff_iso_eta_0`, `eff_id_eta_0`, truth pT histograms
- From data: `h_tight_iso_cluster_0`, `h_tight_noniso_cluster_0`, `h_nontight_iso_cluster_0`, `h_nontight_noniso_cluster_0`
- From response: `response_matrix_full_0` (RooUnfoldResponse)

**Output file:**
```
{final_outfile}_{var_type}.root          (data)
{final_outfile}_{var_type}_mc.root       (MC closure)
```
Example: `results/Photon_final_bdt_nom.root`, `results/Photon_final_bdt_nom_mc.root`

**Key histograms written:**
- `h_unfold_sub_result` -- final unfolded, efficiency-corrected cross-section (THE primary result)
- `h_unfold_sub_result_woeff` -- unfolded without efficiency correction
- `gpurity`, `gpurity_leak` -- purity TGraphErrors
- `g_purity_truth` -- MC truth purity (MC mode only)
- `f_purity_fit`, `f_purity_leak_fit` -- purity fit functions
- `confInt`, `confInt_leak` -- confidence interval graphs
- `h_R`, `h_R_notmatch` -- R-factor histograms
- `h_leak_B`, `h_leak_C`, `h_leak_D` -- leakage fractions
- All 10 unfolding iterations: `h_unfold_sub_{1..10}`, `h_unfold_sub_leak_{1..10}`
- `h_data_sub_leak` -- purity-corrected data yield
- `h_truth_pT_0` -- PYTHIA8 truth spectrum
- `h_common_cluster_0`, `h_tight_iso_cluster_0` -- passed through from input

---

### Stage 5a: systematic aggregation (calc_syst_bdt.py)

**Producer:** `plotting/calc_syst_bdt.py`

**Input:** Multiple `Photon_final_bdt_{variant_name}.root` files from `efficiencytool/results/`
- Reads `h_unfold_sub_result` from each variant and nominal

**Output (to `plotting/rootFiles/`):**
- Per-type: `syst_bdt_{type}.root` (h_dev_low, h_dev_high, h_dev_rel_low, h_dev_rel_high)
- Per-group: `syst_bdt_{group}.root`
- Total: `syst_bdt_total.root`
- **`syst_sum.root`** (h_sum_low, h_sum_high, h_sum_rel_low, h_sum_rel_high) -- consumed by plot_final_selection.C

**Output figures (to `plotting/figures/`):**
- `syst_bdt_rel_{type}.pdf`
- `syst_bdt_breakdown.pdf`

---

### Stage 5b: final plots (plot_final_selection.C)

**Producer:** `plotting/plot_final_selection.C`

**Input files:**
- `efficiencytool/results/Photon_final_{tune}.root` -- reads `h_unfold_sub_result`, `h_truth_pT_0`, `h_common_cluster_0`, `h_tight_iso_cluster_0`, `h_data_sub_leak`, `h_unfold_sub_result_woeff`
- `efficiencytool/results/Photon_final_{tune}_mc.root` -- MC closure
- `plotting/rootFiles/syst_sum.root` -- reads `h_sum_low`, `h_sum_high`
- `NLO/rootFiles/jetPHOX_{10,05,20}.root` -- reads `h_truth_pT`
- `plotting/sphenix_nlo/photons_newphenix_sc{1,05,2}.dat` -- Werner Vogelsang NLO text data
- Luminosity from config YAML `efficiencytool/config_{tune}.yaml: analysis.lumi`

**Output figures (to `plotting/figures/`):**
- `final_{tune}.pdf` -- primary cross-section plot with ratio panel
- `final_common_cluster_{tune}.pdf`
- `final_tight_iso_cluster_{tune}.pdf`
- `final_all_{tune}.pdf` -- pipeline stages
- `final_phenix_{tune}.pdf` -- PHENIX comparison

### Stage 5c: additional plots --> analysis note

Several plotting macros save directly to the analysis note:
- `plotting/plot_showershapes.C` --> `PPG12-analysis-note/Figures/showershapes/`
- `plotting/plot_showershapes_selections.C` --> `PPG12-analysis-note/Figures/showershapes_selections/`
- `plotting/plot_npb_time_bkgsub.C` --> `PPG12-analysis-note/Figures/npb_time_bkgsub/`
- `plotting/plot_vertex_check.C` --> `PPG12-analysis-note/Figures/`

Most other plotting macros save to `plotting/figures/`.

---

## 3. Branch Contract: anatreemaker --> BDT Pipeline

### Branches written by CaloAna24.cc and read by BDTinput.C

All branches match. BDTinput.C reads a strict subset of what CaloAna24 writes. Key branches at this interface:

| Category | Branches | Written by CaloAna24 | Read by BDTinput.C |
|----------|----------|---------------------|--------------------|
| Event | `mbdnorthhit`, `mbdsouthhit`, `vertexz`, `energy_scale`, `runnumber` | Yes | Yes |
| Trigger | `scaledtrigger[64]`, `livetrigger[64]` | Yes | Yes |
| Truth | `nparticles`, `particle_trkid[]`, `particle_photonclass[]` | Yes (MC) | Yes |
| Cluster count | `ncluster_{node}` | Yes | Yes |
| Kinematics | `cluster_{Et,Eta,Phi}_{node}` | Yes | Yes |
| Shower shape | `cluster_{et1-4,weta_cogx,wphi_cogx}_{node}` | Yes | Yes |
| Energy windows | `cluster_e{11,22,13,15,17,31,51,71,33,35,37,53,73,77}_{node}`, `cluster_e32_{node}` | Yes | Yes |
| Widths | `cluster_{w32,w52,w72}_{node}` | Yes | Yes |
| Isolation | `cluster_iso_{02,03,04}_{node}`, `cluster_iso_03_60_{emcal,hcalin,hcalout}_{node}` | Yes | Yes |
| ID | `cluster_prob_{node}`, `cluster_truthtrkID_{node}`, `cluster_nsaturated_{node}` | Yes | Yes |
| NPB data | `mbd_time`, `njet`, `jet_Pt[]`, `jet_Phi[]`, `cluster_time_array[]`, `cluster_e_array[]`, `cluster_ownership_array[]` | Yes | Yes (data NPB only) |

### Branches written by CaloAna24.cc and read by apply_BDT.C

All branches match. apply_BDT.C reads the same core cluster branches as BDTinput.C (kinematics, shower shape, energy windows, widths) plus isolation branches for NPB score computation. No mismatches found.

**One notable difference in derived feature computation:**
- BDTinput.C uses `safe_div` lambda: checks `denominator != 0 && isfinite(result)`
- apply_BDT.C checks `numerator > 0` (e.g., `cluster_e11 > 0`)
- For positive energies this is functionally equivalent, but edge cases with zero numerator and nonzero denominator would produce 0.0 in both, so no actual mismatch in practice.

---

## 4. Branch Contract: apply_BDT.C --> Efficiency Tool

### Branches added by apply_BDT.C and consumed by RecoEffCalculator_TTreeReader.C

| Branch Pattern | Producer | Consumer | Notes |
|----------------|----------|----------|-------|
| `cluster_bdt_{node}_{model_name}[ncluster]` | apply_BDT.C | RecoEffCalculator | Dynamically loaded per config; model names from `bdt_model_name` and `bdt_et_bin_models` |
| `cluster_npb_score_{node}[ncluster]` | apply_BDT.C | RecoEffCalculator | Used for NPB rejection (`npb_cut_on`, `npb_score_cut`) |

**Contract verification:** RecoEffCalculator builds TTreeReaderArray objects for BDT branches based on the model names specified in the YAML config (`bdt_model_name`, `bdt_et_bin_models`). It deduplicates model names before creating readers. The branch name format `cluster_bdt_{clusternodename}_{model_name}` in the consumer exactly matches the `Form("cluster_bdt_%s_%s", ...)` pattern used by the producer.

**Potential failure mode:** If the config used by RecoEffCalculator requests a BDT model name that was not applied by apply_BDT.C, the TTreeReaderArray creation will fail at runtime. There is no compile-time or config-time validation of model name consistency.

---

## 5. File Naming Conventions

### Pattern: `{base}_{qualifier}_{var_type}.root`

The `var_type` field from config is the universal suffix that tags all outputs from a given analysis configuration. It prevents filename collisions between systematic variations.

| Stage | Base | Qualifier | Example |
|-------|------|-----------|---------|
| RecoEffCalculator (sim) | `MC_efficiency` | `{filetype}` | `MC_efficiency_photon10_bdt_nom.root` |
| RecoEffCalculator (data) | `data_histo` | (none) | `data_histo_bdt_nom.root` |
| Response matrix (sim) | `MC_response` | `{filetype}` | `MC_response_photon10_bdt_nom.root` |
| MergeSim (signal) | `MC_efficiency` | (none) | `MC_efficiency_bdt_nom.root` |
| MergeSim (jet) | `MC_efficiency` | `jet` | `MC_efficiency_jet_bdt_nom.root` |
| MergeSim (response) | `MC_response` | (none) | `MC_response_bdt_nom.root` |
| CalculatePhotonYield (data) | `Photon_final` | (none) | `Photon_final_bdt_nom.root` |
| CalculatePhotonYield (MC) | `Photon_final` | (none) | `Photon_final_bdt_nom_mc.root` |
| Vtxscan | `{base}_{qualifier}` | `vtxscan` | `MC_efficiency_photon10_bdt_nom_vtxscan.root` |
| Systematics | `syst_bdt` | `{type/group}` | `syst_bdt_total.root`, `syst_sum.root` |
| Plots | `{plot_type}` | `{tune}` | `final_bdt_nom.pdf`, `eff_reco_bdt_nom.pdf` |

### BDT-scored tree naming:

| Context | Pattern |
|---------|---------|
| Sim (apply_BDT.C output) | `FunWithxgboost/{filetype}/bdt_split.root` |
| Data (apply_BDT.C output) | `{original_name}_with_bdt_split.root` |

---

## 6. Config Propagation

### Two config families

1. **Training config** (`FunWithxgboost/config.yaml`): controls BDT training features, hyperparameters, ET binning, reweighting. Read by `main_training.py`.

2. **Analysis config** (`FunWithxgboost/config_nom.yaml` or `efficiencytool/config_bdt_*.yaml`): controls cuts, paths, sample selection, BDT model choice, output naming. Read by BDTinput.C, apply_BDT.C, RecoEffCalculator, MergeSim, CalculatePhotonYield, ShowerShapeCheck, DoubleInteractionCheck, and plotting macros.

### Config flow through pipeline stages

| Stage | Config File | Key Fields Used |
|-------|-------------|-----------------|
| BDTinput.C | `config_nom.yaml` | `input.{cluster_node_name, tree, photon_jet_file_root_dir, data_file}`, `analysis.{vertex_cut, trigger_used, npb_cut_on, write_outside_common_pass_mc}` |
| main_training.py | `config.yaml` | `data.{features, single_file_set}`, `binning.*`, `model.params.*`, `reweighting.*`, `output.*` |
| apply_BDT.C | `config_nom.yaml` | `input.cluster_node_name`, `analysis.use_split_bdt` |
| RecoEffCalculator | `config_bdt_*.yaml` | `input.*` (paths, BDT model names), `output.*` (output paths, var_type), `analysis.*` (all cuts, efficiencies) |
| MergeSim.C | `config_bdt_*.yaml` | `output.{eff_outfile, response_outfile, var_type}` |
| CalculatePhotonYield.C | `config_bdt_*.yaml` | `output.*`, `analysis.{lumi, unfold.*, pT_bins, pT_bins_truth, fit_option, mc_purity_correction, fittingerror}` |
| ShowerShapeCheck.C | `config_showershape_*.yaml` | Same schema as efficiency configs |
| plot_final_selection.C | reads config for lumi only | `efficiencytool/config_{tune}.yaml: analysis.lumi` |
| calc_syst_bdt.py | Imports VARIANTS from make_bdt_variations.py | Variant names map to `Photon_final_bdt_{name}.root` |

### Config path difference between BDT and efficiency stages

A critical detail: `config_nom.yaml` (used by BDTinput.C, apply_BDT.C) and `config_bdt_nom.yaml` (used by efficiency tool) have DIFFERENT `photon_jet_file_root_dir` values:

- `config_nom.yaml`: `photon_jet_file_root_dir: /sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/`
  - + `photon_jet_file_branch_dir: /condorout/combined.root`
  - Resolves to: `anatreemaker/macro_maketree/sim/run28/{sample}/condorout/combined.root` (pre-BDT trees)
  
- `config_bdt_nom.yaml`: `photon_jet_file_root_dir: /sphenix/user/shuhangli/ppg12/FunWithxgboost/`
  - + `photon_jet_file_branch_dir: /bdt_split.root`
  - Resolves to: `FunWithxgboost/{sample}/bdt_split.root` (BDT-scored trees)

This is correct: BDTinput.C and apply_BDT.C read/write from the anatreemaker output, while the efficiency tool reads the BDT-scored output from FunWithxgboost.

---

## 7. Mismatches and Gaps

### 7.1 Division guard inconsistency in energy ratios

- **BDTinput.C** `safe_div`: checks `denominator != 0 && isfinite(result)`, returns 0.0 on failure
- **apply_BDT.C**: checks `numerator > 0` (e.g., `cluster_e11 > 0`), returns 0 on failure

These produce identical results when both numerator and denominator are positive (the normal case), but edge cases differ: a cluster with e11=0 and e33>0 gives 0.0 in both, but a cluster with e11>0 and e33=0 gives 0.0 in BDTinput.C (div-by-zero caught) and inf in apply_BDT.C (not caught). In practice, e33=0 with e11>0 is physically impossible (e11 is a subset of e33), so this is cosmetic.

### 7.2 Training text vs training config file set divergence

The training config `config.yaml` references `NOSPLIT/shapes_*.txt` files (legacy naming) while BDTinput.C currently outputs `shapes_split_{filetype}.txt` files (current naming). The `NOSPLIT/` files appear to be from an earlier run of BDTinput.C with a different cluster node. If the training were re-run from scratch, the config.yaml `single_file_set` paths would need updating.

### 7.3 Jet sample set evolution

The current pipeline uses jet5/8/12/20/30/40 (6 samples), while BDTinput.C and the training config reference jet10/15/20/30/50 (5 samples). These are different jet sample sets:
- Training (BDTinput.C + config.yaml): `jet10, jet15, jet20, jet30, jet50`
- Efficiency (oneforall_tree.sh + MergeSim.C): `jet5, jet8, jet12, jet20, jet30, jet40`

This is NOT a bug -- the training and efficiency pipelines intentionally use different jet MC samples. The training needs broad background representation; the efficiency tool needs finer binning matched to CrossSectionWeights.h. However, the jet samples used for training (jet10/15/50) do not have weight entries in `CrossSectionWeights.h` for the newer pipeline samples (jet5/8/12/40), and vice versa. Both sets have their own weight constants and are self-consistent within their respective pipelines.

### 7.4 Dead code in oneforall.sh

`oneforall.sh` has commented-out `RecoEffCalculator.C` (old non-TTreeReader version) calls for jet10/15/20/30 (the old jet sample set). These are legacy references. The active pipeline uses `oneforall_tree.sh` which calls `RecoEffCalculator_TTreeReader.C` for the current jet sample set.

### 7.5 Potential model name mismatch between apply_BDT and efficiency config

If `config_bdt_nom.yaml` specifies `bdt_et_bin_models: ["base_v1E", "base_v1E"]` but apply_BDT.C was run with a config that did not produce the `cluster_bdt_CLUSTERINFO_CEMC_base_v1E` branch (e.g., because a different config or older version of apply_BDT.C was used), RecoEffCalculator will crash at TTreeReaderArray creation. There is no validation layer between these stages.

All 11 model variants ARE written by the current apply_BDT.C regardless of config, so this is not currently a problem. But if apply_BDT.C is modified to write only selected models, consistency would need to be enforced.

### 7.6 pT bin mismatch between config_nom.yaml and plotcommon.h

`config_nom.yaml` uses `pT_bins: [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35]` (10 bins), while `plotcommon.h` uses `ptRanges[13] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35}` (12 bins). The `config_bdt_nom.yaml` uses `pT_bins: [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]` (12 bins, but different upper edges: 32/36 vs 30/35).

This means the pT bins used in the efficiency tool do NOT exactly match the plotting bins in plotcommon.h. The efficiency tool uses config-driven bins, and the plotting macros typically read the histograms as-is without rebinning, so any mismatch would manifest as incorrect bin labels or axis ranges in plots.

### 7.7 config_nom.yaml vs config_bdt_nom.yaml isolation parameters differ

- `config_nom.yaml`: `reco_iso_max_b: 1.08128`, `reco_iso_max_s: 0.0299107`
- `config_bdt_nom.yaml`: `reco_iso_max_b: 0.502095`, `reco_iso_max_s: 0.0433036`

These are intentionally different configurations (nom vs bdt_nom), but it means BDTinput.C (which uses config_nom.yaml) applies different isolation criteria than the efficiency tool (which uses config_bdt_nom.yaml). Since the BDT training is isolation-agnostic (recoisoET is written but not used as a training feature), this is not a problem.

### 7.8 Produced but not consumed: several branches

CaloAna24.cc writes many branches that no downstream consumer reads:
- `cluster_CNN_prob` (always -1, ONNX disabled)
- `cluster_merged_prob`
- `cluster_ecore`
- `cluster_embed_id`
- `cluster_iso_03_sub1_*`, `cluster_iso_04_sub1_*` (background-subtracted isolation)
- `cluster_iso_topo_soft_*` (soft topo isolation)
- `cluster_iso_04_*` (R=0.4 non-topo, non-sub1) -- only topo_04 or standard iso_04 is used depending on config
- `cluster_{e55,e57,e75}` (only RecoEffCalculator reads them, but no analysis cut or histogram uses them)
- `cluster_{e52,e72}` (read by RecoEffCalculator but not used in any cut or fill besides being available)
- `cluster_{e1-4}` (raw quadrant energies, as opposed to et1-4 fractions)
- `cent`, `Psi2`, `eventnumber` (event-level)
- `daughter_*` branches (single-particle MC only)
- `totalEMCal_energy`, `totalIHCal_energy`, `totalOHCal_energy`
- `currentscaler_*[64]`
- Truth jet container names with `_AntiKt_Truth_r04` suffix (only used for double-interaction MC)
- `tracking_radiograph` histogram

These are not bugs -- they represent exploratory or future-use branches. The tree is intentionally "wide" to support multiple analysis paths.

---

## 8. Summary of Key Handoff Points

| Handoff | Producer File | Consumer File | Data Format | Key Names at Boundary |
|---------|--------------|---------------|-------------|----------------------|
| DST --> slimtree | CaloAna24.cc | BDTinput.C, apply_BDT.C | ROOT TTree `slimtree` | ~100 branches per cluster node |
| slimtree --> training txt | BDTinput.C | main_training.py / data_loader.py | Space-separated text | 30 columns, `pid` for labeling |
| training txt --> models | main_training.py | apply_BDT.C | TMVA ROOT | `model_{name}_single_tmva.root`, positional features |
| slimtree + models --> scored slimtree | apply_BDT.C | RecoEffCalculator | ROOT TTree `slimtree` | `cluster_bdt_{node}_{model}`, `cluster_npb_score_{node}` |
| scored slimtree --> efficiency ROOT | RecoEffCalculator | MergeSim.C | ROOT histograms | `h_tight_iso_cluster_*`, TEfficiency, response matrices |
| per-sample --> merged MC | MergeSim.C | CalculatePhotonYield.C | ROOT histograms | Same histogram names, additively merged |
| merged MC + data --> cross-section | CalculatePhotonYield.C | plot_final_selection.C, calc_syst_bdt.py | ROOT histograms | `h_unfold_sub_result`, `gpurity`, `h_leak_*` |
| cross-section variants --> systematics | calc_syst_bdt.py | plot_final_selection.C | ROOT histograms | `syst_sum.root: h_sum_low/high` |
| cross-section + systematics --> figures | plot_final_selection.C | Analysis note LaTeX | PDF figures | `figures/final_{tune}.pdf` |
