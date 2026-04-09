# anatreemaker -- DST to Slim ROOT Tree Converter

## 1. Purpose

The `anatreemaker/` directory contains the first stage of the PPG12 analysis pipeline. It converts sPHENIX DST (Data Summary Tape) files into compact ROOT trees (`slimtree`) that contain per-event cluster-level, particle-level, jet-level, and MBD-level information. This "slim tree" is the common input for all downstream steps: BDT training (`BDTinput.C`), BDT scoring (`apply_BDT.C`), efficiency/yield calculations (`CalculatePhotonYield.C`), and shower-shape studies.

The core code is a Fun4All `SubsysReco` module called `CaloAna24`, compiled as a shared library (`libCaloAna24.so`). It is invoked within Fun4All macro scripts that set up the full reconstruction chain (clustering, isolation, jet finding) before calling CaloAna24 to extract variables and write the output tree.

## 2. Key Files

| File | Lines | Description |
|------|-------|-------------|
| `source/CaloAna24.cc` | 2321 | Main analysis module implementation -- event processing, cluster loops, shower-shape computation, isolation ET, truth matching, jet extraction |
| `source/CaloAna24.h` | 385 | Header -- class declaration, member variables, array dimensions, helper function signatures |
| `source/Makefile.am` | 54 | Autotools build file -- links against calo_io, g4eval, globalvertex, mbd_io, sph_onnx |
| `source/configure.ac` | 17 | Autoconf configuration |
| `source/autogen.sh` | 8 | Build bootstrap script (aclocal, libtoolize, automake, autoconf) |
| `macro_maketree/data/ana521/Fun4All_runDST.C` | 393 | Fun4All macro for **data** processing (latest production tag `ana521_2025p007_v001`) |
| `macro_maketree/data/ana521/run.sh` | 150 | Condor orchestration for data -- splits DST lists per run, creates per-job sub-files, submits to HTCondor |
| `macro_maketree/data/ana521/CondorRun.sh` | 20 | Condor job wrapper for data -- sources environment, runs Fun4All macro |
| `macro_maketree/sim/run28/photon10/Fun4All_run_sim.C` | 461 | Fun4All macro for **simulation** processing (run 28, MC) |
| `macro_maketree/sim/run28/photon10/run_condor.sh` | 57 | Condor submission for sim -- splits DST file lists, creates per-job directories |
| `macro_maketree/sim/run28/photon10/CondorRunSim_new.sh` | 19 | Condor job wrapper for sim |
| `macro_maketree/sim/create_file_lists.sh` | 53 | Helper: runs `CreateFileList.pl` for all sim samples across run periods |
| `macro_maketree/sim/run28/hadd_combined.sh` | 70 | Merges per-job `caloana.root` files into `combined.root` via `hadd` |
| `macro_maketree/sim/run28/cleanup_and_run.sh` | 57 | Removes old condor output and resubmits all sim jobs |
| `macro_maketree/sim/run28/cleanup.sh` | 41 | Removes old OutDir* condor output directories only |

## 3. Directory Structure

```
anatreemaker/
  source/
    CaloAna24.cc          -- main module
    CaloAna24.h           -- header
    Makefile.am            -- build system
    configure.ac
    autogen.sh
    build/                 -- compiled objects (RHEL9)
    build_SL7/             -- compiled objects (SL7, legacy)
  macro_maketree/
    data/
      ana450/              -- older production (legacy)
      ana462/              -- older production
      ana462_play/         -- test/playground
      ana468/              -- older production
      ana468newcalib/       -- recalibration variant
      ana484/
      ana502/
      ana509/
      ana518/
      ana521/              -- CURRENT production (tag ana521_2025p007_v001)
        Fun4All_runDST.C   -- data macro
        run.sh             -- condor orchestration
        CondorRun.sh        -- per-job wrapper
        runList.txt         -- list of run numbers to process
        dstLists/           -- per-run DST file lists
        condorout/          -- condor job output directories
      auau_data_new/       -- Au+Au data processing
      auau_embed/          -- Au+Au embedding
      DSTlists/            -- shared DST file lists
    sim/
      create_file_lists.sh -- generates DST file lists for sim
      run28/               -- MC production at run 28 conditions
        photon5/           -- signal: 0-14 GeV truth pT
        photon10/          -- signal: 14-30 GeV truth pT
        photon20/          -- signal: 30+ GeV truth pT
        photon10_double/   -- double-interaction signal MC
        jet5/              -- bkg: finer binning for DI study
        jet8/              -- bkg: finer binning for DI study
        jet10/             -- bkg: 10-14 GeV
        jet12/             -- bkg: 14-21 GeV
        jet12_double/      -- double-interaction jet MC
        jet15/             -- bkg: 15-20 GeV
        jet20/             -- bkg: 20-30 GeV
        jet30/             -- bkg: 30-50 GeV
        jet40/             -- bkg: finer binning for DI study
        jet50/             -- bkg: 50+ GeV
        jet70/             -- bkg: high pT
        mb/                -- minimum bias
        cleanup.sh
        cleanup_and_run.sh
        hadd_combined.sh
      run29/               -- MC production at run 29 conditions
      test/                -- test directory
    hadd_tmp/              -- scratch space for hadd merging
```

Each sim sample subdirectory contains:
- `Fun4All_run_sim.C` -- macro (sample-specific only in file list names)
- `run_condor.sh` -- condor submission (configures DST list paths, lines per job)
- `CondorRunSim_new.sh` -- per-job wrapper
- `g4hits.list`, `dst_calo_cluster.list`, `dst_mbd_epd.list`, `dst_truth_jet.list` -- DST file lists
- `condorout/` -- output ROOT files in `OutDir*/caloana.root`
- `condorout/combined.root` -- merged output

## 4. Entry Points

### Building the library

```bash
cd anatreemaker/source
./autogen.sh --prefix=$MYINSTALL
cd build && make install
```

This produces `libCaloAna24.so` which is loaded by Fun4All macros.

### Running on data (current production)

```bash
cd anatreemaker/macro_maketree/data/ana521
# Edit runList.txt with desired run numbers
./run.sh
```

`run.sh` iterates over run numbers in `runList.txt`, fetches DST file lists (from local cache or via `CreateDstList.pl` with tag `ana521_2025p007_v001`), splits them into ~10 files per condor job, and submits via `condor_submit`. Each job runs `CondorRun.sh` which calls:

```
root -b "Fun4All_runDST.C(\"inputdata.txt\",\"inputdatacalo.txt\",0)"
```

Output: `OUTTREE_DST_Jet_*.root` files per job, later merged by the user.

### Running on simulation

```bash
cd anatreemaker/macro_maketree/sim/run28/photon10
# Ensure DST lists exist (or run create_file_lists.sh)
./run_condor.sh
```

Each job runs `CondorRunSim_new.sh` which calls:
```
root "Fun4All_run_sim.C"
```

Output: `caloana.root` per job. Merge with:
```bash
cd anatreemaker/macro_maketree/sim/run28
./hadd_combined.sh
```

### Batch cleanup and resubmission

```bash
cd sim/run28
./cleanup_and_run.sh    # removes OutDir*, resubmits all samples
./cleanup.sh            # removes OutDir* only
```

## 5. Input -- DST Nodes Read

### Data DSTs (Fun4All_runDST.C)

Two input DST streams:
- `DST_Jet` (via `fname0`) -- contains jet-related containers, trigger info
- `DST_JETCALO` (via `fname1`) -- contains calibrated calo tower containers

Production tag: `ana521_2025p007_v001`, CDB global tag: `ProdA_2024`

### Simulation DSTs (Fun4All_run_sim.C)

Five input DST streams:
- `G4Hits` (`inputFile0`) -- GEANT4 hit information
- `DST_CALO_CLUSTER` (`inputFile1`) -- calo cluster containers
- `DST_MBD_EPD` (`inputFile3`) -- MBD and EPD containers
- `DST_TRUTH_JET` (`inputFile4`) -- truth-level jet containers
- `DST_TRUTH` (`inputFile2`) -- truth particle info (currently commented out)

CDB global tag: `MDC2`, RUNNUMBER set to 28

### DST nodes accessed by CaloAna24 (in process_event):

| Node Name | Type | Usage |
|-----------|------|-------|
| `14001` (Gl1Packet) | Gl1Packet | Trigger bits (scaled, live), scaler values (data only) |
| `TriggerRunInfo` | TriggerRunInfo | Trigger prescale values (data only) |
| `EventHeader` | EventHeader | Event sequence number, event validity |
| `PHHepMCGenEventMap` | PHHepMCGenEventMap | HepMC truth record, Pythia process ID, energy scale (MC only) |
| `G4TruthInfo` | PHG4TruthInfoContainer | Truth particles, vertices, showers (MC only) |
| `MbdPmtContainer` | MbdPmtContainer | MBD PMT hits (128 channels), charge, time |
| `MbdOut` | MbdOut | MBD reconstructed times (t0, north, south) |
| `MbdVertexMap` | MbdVertexMap | MBD vertex z (data only) |
| `GlobalVertexMap` | GlobalVertexMap | MBD vertex z via GlobalVertex (MC only) |
| `CentralityInfo` | CentralityInfo | Centrality percentile (mbd_NS) |
| `EventplaneinfoMap` | EventplaneinfoMap | Second-order event plane (sEPDNS, Psi2) |
| `TOWERGEOM_CEMC` | RawTowerGeomContainer | EMCal tower geometry |
| `TOWERGEOM_HCALIN` | RawTowerGeomContainer | IHCal tower geometry |
| `TOWERGEOM_HCALOUT` | RawTowerGeomContainer | OHCal tower geometry |
| `TOWERINFO_CALIB_CEMC` | TowerInfoContainer | Calibrated EMCal towers |
| `TOWERS_CEMC` | TowerInfoContainer | Raw (uncalibrated) EMCal towers (for ADC values) |
| `TOWERINFO_CALIB_HCALIN` | TowerInfoContainer | Calibrated IHCal towers |
| `TOWERS_HCALIN` | TowerInfoContainer | Raw IHCal towers |
| `TOWERINFO_CALIB_HCALOUT` | TowerInfoContainer | Calibrated OHCal towers |
| `TOWERS_HCALOUT` | TowerInfoContainer | Raw OHCal towers |
| `TOWERINFO_CALIB_CEMC_RETOWER_SUB1` | TowerInfoContainer | Background-subtracted EMCal towers |
| `TOWERINFO_CALIB_HCALIN_SUB1` | TowerInfoContainer | Background-subtracted IHCal towers |
| `TOWERINFO_CALIB_HCALOUT_SUB1` | TowerInfoContainer | Background-subtracted OHCal towers |
| `TOWERINFO_CALIB_CEMC_RETOWER` | TowerInfoContainer | Retowered EMCal (for jet finding) |
| `CLUSTERINFO_CEMC_NO_SPLIT` | RawClusterContainer | EMCal clusters without subcluster splitting |
| `CLUSTERINFO_CEMC` | RawClusterContainer | EMCal clusters with subcluster splitting (default) |
| `TOPOCLUSTER_ALLCALO` | RawClusterContainer | Topological clusters (all calo, standard thresholds) |
| `TOPOCLUSTER_ALLCALO_SOFT` | RawClusterContainer | Topological clusters (all calo, softer thresholds) |
| `AntiKt_unsubtracted_r04` | JetContainer | Anti-kT R=0.4 unsubtracted jets (default) |
| `AntiKt_Truth_r04` | JetContainer | Truth-level anti-kT R=0.4 jets (MC only) |
| `PHGenIntegral` | PHGenIntegral | MC cross-section normalization info |
| `CEMC` (via CaloEvalStack) | CaloEvalStack | Cluster-to-truth matching (MC only) |

### Reconstruction modules registered in Fun4All macros

Both data and sim macros register a full reconstruction chain before CaloAna24:

1. `Process_Calo_Calib()` -- standard sPHENIX calibration (towers, status)
2. `RawClusterBuilderTemplate` -- template-based EMCal clustering (threshold 0.070 GeV), registered twice:
   - `CLUSTERINFO_CEMC_NO_SPLIT` (splitting disabled)
   - `CLUSTERINFO_CEMC` (splitting enabled, default from DST)
3. `ClusterIso` -- cone isolation at R=0.2, 0.3, 0.4 for each cluster container
4. `RawClusterBuilderTopo` -- topological clustering for all calorimeters:
   - `TOPOCLUSTER_ALLCALO` (significance 4.0/2.0/1.0)
   - `TOPOCLUSTER_ALLCALO_SOFT` (significance 3.0/2.0/0.0)
5. `RetowerCEMC` -- retowers EMCal to HCal granularity for jet finding
6. `JetReco` + `FastJetAlgoSub` -- anti-kT R=0.4 jet reconstruction

**Sim-only additional modules:**
- `MbdReco` -- MBD reconstruction
- `GlobalVertexReco` -- global vertex reconstruction
- `JetCalib` -- jet energy calibration (eta + zvtx dependent), output `AntiKt_unsubtracted_r04_calib`
- `TimerStats` -- timing statistics

**Sim-only disabled block (if false):** Full waveform simulation chain (CaloWaveformSim, CaloTowerBuilder, CaloTowerStatus, CaloTowerCalib for all three calorimeters). This block is wrapped in `if(false)` and is not currently active.

## 6. Output -- `slimtree` Branch Catalog

The output ROOT file contains a TTree named `slimtree` and two histograms. All cluster-level branches are duplicated for each cluster container (indexed by `{node}` = `CLUSTERINFO_CEMC_NO_SPLIT` or `CLUSTERINFO_CEMC`). Jet branches are indexed by jet container name (default: `AntiKt_unsubtracted_r04`).

### 6.1 Event-Level Branches

| Branch | Type | Description |
|--------|------|-------------|
| `runnumber` | `Int_t` | Run number |
| `eventnumber` | `Int_t` | Event sequence number |
| `mbd_time` | `Float_t` | MBD t0 |
| `mbd_north_time` | `Float_t` | MBD north time |
| `mbd_south_time` | `Float_t` | MBD south time |
| `mbdnorthhit` | `Int_t` | Number of MBD north hits (charge > 0.4) |
| `mbdsouthhit` | `Int_t` | Number of MBD south hits (charge > 0.4) |
| `mbdnorthq[64]` | `Float_t[64]` | MBD north PMT charges |
| `mbdsouthq[64]` | `Float_t[64]` | MBD south PMT charges |
| `mbdnortht[64]` | `Float_t[64]` | MBD north PMT times (ns) |
| `mbdsoutht[64]` | `Float_t[64]` | MBD south PMT times (ns) |
| `mbdnorthqsum` | `Float_t` | Sum of MBD north charges |
| `mbdsouthqsum` | `Float_t` | Sum of MBD south charges |
| `mbdnorthtmean` | `Float_t` | Mean time of north hits (charge > 0.4) |
| `mbdsouthtmean` | `Float_t` | Mean time of south hits (charge > 0.4) |
| `vertexz` | `Float_t` | Reconstructed vertex z (MBD vertex for data, GlobalVertex/MBD for MC, truth vertex for single particle) |
| `vertexz_truth` | `Float_t` | Truth vertex z (MC only, -9999 for data) |
| `cent` | `Float_t` | Centrality percentile (mbd_NS, -999.99 if unavailable) |
| `Psi2` | `Float_t` | Second-order event plane angle (sEPDNS shifted) |
| `pythiaid` | `Int_t` | Pythia process ID (MC only, -9999 for data) |
| `energy_scale` | `Float_t` | Event energy scale / Q^2 (MC only, -1 for data) |
| `scaledtrigger[64]` | `Bool_t[64]` | Scaled trigger bit decisions (data only) |
| `livetrigger[64]` | `Bool_t[64]` | Live trigger bit decisions (data only) |
| `currentscaler_raw[64]` | `Long64_t[64]` | Raw scaler counts |
| `currentscaler_live[64]` | `Long64_t[64]` | Live scaler counts |
| `currentscaler_scaled[64]` | `Long64_t[64]` | Scaled scaler counts |
| `trigger_prescale[64]` | `Float_t[64]` | Trigger prescale factors (first 32 bits filled) |
| `totalEMCal_energy` | `Float_t` | Total energy in good EMCal towers |
| `totalIHCal_energy` | `Float_t` | Total energy in good IHCal towers |
| `totalOHCal_energy` | `Float_t` | Total energy in good OHCal towers |

### 6.2 Truth Particle Branches (MC only)

Array dimension: `nparticles` (max 10000)

| Branch | Type | Description |
|--------|------|-------------|
| `nparticles` | `Int_t` | Number of stored truth particles |
| `particle_E[nparticles]` | `Float_t` | Particle energy |
| `particle_Pt[nparticles]` | `Float_t` | Particle pT |
| `particle_Eta[nparticles]` | `Float_t` | Particle pseudorapidity |
| `particle_Phi[nparticles]` | `Float_t` | Particle azimuthal angle |
| `particle_pid[nparticles]` | `Int_t` | PDG particle ID (22=photon, 111=pi0-decay photon, 221=eta-decay photon) |
| `particle_trkid[nparticles]` | `Int_t` | GEANT track ID |
| `particle_photonclass[nparticles]` | `Int_t` | Photon classification: 1=direct, 2=fragmentation, 3=decay, 0=unidentified, -1=not photon |
| `particle_photon_mother_pid[nparticles]` | `Int_t` | PID of incoming particle at production vertex |
| `particle_truth_iso_02[nparticles]` | `Float_t` | Truth isolation ET in R<0.2 cone |
| `particle_truth_iso_03[nparticles]` | `Float_t` | Truth isolation ET in R<0.3 cone |
| `particle_truth_iso_04[nparticles]` | `Float_t` | Truth isolation ET in R<0.4 cone |
| `particle_converted[nparticles]` | `Int_t` | Conversion flag: 0=not converted, 1=pair-converted (e+e-), 2=other interaction before CEMC |

### 6.3 Daughter Branches (single-particle MC only)

Array dimension: `ndaughter` (max 100)

| Branch | Type | Description |
|--------|------|-------------|
| `ndaughter` | `Int_t` | Number of secondary daughters |
| `daughter_pid[ndaughter]` | `Int_t` | Daughter PDG ID |
| `daughter_E[ndaughter]` | `Float_t` | Daughter energy |
| `daughter_Pt[ndaughter]` | `Float_t` | Daughter pT |
| `daughter_Eta[ndaughter]` | `Float_t` | Daughter pseudorapidity |
| `daughter_Phi[ndaughter]` | `Float_t` | Daughter azimuthal angle |
| `daughter_vtx_x[ndaughter]` | `Float_t` | Daughter production vertex x |
| `daughter_vtx_y[ndaughter]` | `Float_t` | Daughter production vertex y |
| `daughter_vtx_z[ndaughter]` | `Float_t` | Daughter production vertex z |

### 6.4 Cluster Branches (per cluster container {node})

Array dimension: `ncluster_{node}` (max 10000). Two containers: `CLUSTERINFO_CEMC_NO_SPLIT` and `CLUSTERINFO_CEMC`.

#### Basic kinematics and identification

| Branch | Type | Description |
|--------|------|-------------|
| `ncluster_{node}` | `Int_t` | Number of clusters in this container |
| `cluster_E_{node}[n]` | `Float_t` | Cluster energy (GeV) |
| `cluster_ecore_{node}[n]` | `Float_t` | Cluster core energy (GeV) |
| `cluster_Et_{node}[n]` | `Float_t` | Cluster transverse energy = E / cosh(eta) |
| `cluster_Eta_{node}[n]` | `Float_t` | Cluster pseudorapidity (vertex-corrected) |
| `cluster_Phi_{node}[n]` | `Float_t` | Cluster azimuthal angle (vertex-corrected) |
| `cluster_prob_{node}[n]` | `Float_t` | Cluster chi2 probability |
| `cluster_merged_prob_{node}[n]` | `Float_t` | Merged cluster probability |
| `cluster_CNN_prob_{node}[n]` | `Float_t` | CNN classification probability (currently always -1, ONNX disabled) |
| `cluster_truthtrkID_{node}[n]` | `Int_t` | Truth particle track ID matched to cluster (MC only) |
| `cluster_pid_{node}[n]` | `Int_t` | PID of matched truth particle (MC only; 111/221 for pi0/eta decay photons) |
| `cluster_embed_id_{node}[n]` | `Int_t` | Embedding ID of matched truth particle (MC only) |

#### Isolation energy (tower-level, multiple cone sizes and thresholds)

All isolation branches subtract the cluster's own ET for EMCal cones. Format: `cluster_iso_{R}_{threshold}_{calo}_{node}`.

| Branch | Type | Description |
|--------|------|-------------|
| `cluster_iso_02_{node}[n]` | `Float_t` | RawCluster isolation R=0.2 (from ClusterIso module) |
| `cluster_iso_03_{node}[n]` | `Float_t` | RawCluster isolation R=0.3 |
| `cluster_iso_04_{node}[n]` | `Float_t` | RawCluster isolation R=0.4 |
| `cluster_iso_03_emcal_{node}[n]` | `Float_t` | Tower-level EMCal ISO ET, R=0.3, no tower E cut, minus cluster ET |
| `cluster_iso_03_hcalin_{node}[n]` | `Float_t` | Tower-level IHCal ISO ET, R=0.3, no tower E cut |
| `cluster_iso_03_hcalout_{node}[n]` | `Float_t` | Tower-level OHCal ISO ET, R=0.3, no tower E cut |
| `cluster_iso_03_60_emcal_{node}[n]` | `Float_t` | EMCal ISO, R=0.3, tower E > 0.06 GeV, minus cluster ET |
| `cluster_iso_03_60_hcalin_{node}[n]` | `Float_t` | IHCal ISO, R=0.3, tower E > 0.06 GeV |
| `cluster_iso_03_60_hcalout_{node}[n]` | `Float_t` | OHCal ISO, R=0.3, tower E > 0.06 GeV |
| `cluster_iso_03_70_emcal_{node}[n]` | `Float_t` | EMCal ISO, R=0.3, tower E > 0.07 GeV, minus cluster ET |
| `cluster_iso_03_70_hcalin_{node}[n]` | `Float_t` | IHCal ISO, R=0.3, tower E > 0.07 GeV |
| `cluster_iso_03_70_hcalout_{node}[n]` | `Float_t` | OHCal ISO, R=0.3, tower E > 0.07 GeV |
| `cluster_iso_005_70_emcal_{node}[n]` | `Float_t` | EMCal ISO, R=0.05, tower E > 0.07 GeV, minus cluster ET |
| `cluster_iso_01_70_emcal_{node}[n]` | `Float_t` | EMCal ISO, R=0.1, tower E > 0.07 GeV, minus cluster ET |
| `cluster_iso_02_70_emcal_{node}[n]` | `Float_t` | EMCal ISO, R=0.2, tower E > 0.07 GeV, minus cluster ET |
| `cluster_iso_03_120_emcal_{node}[n]` | `Float_t` | EMCal ISO, R=0.3, tower E > 0.12 GeV, minus cluster ET |
| `cluster_iso_03_120_hcalin_{node}[n]` | `Float_t` | IHCal ISO, R=0.3, tower E > 0.12 GeV |
| `cluster_iso_03_120_hcalout_{node}[n]` | `Float_t` | OHCal ISO, R=0.3, tower E > 0.12 GeV |
| `cluster_iso_04_emcal_{node}[n]` | `Float_t` | EMCal ISO, R=0.4, no tower E cut, minus cluster ET |
| `cluster_iso_04_hcalin_{node}[n]` | `Float_t` | IHCal ISO, R=0.4 |
| `cluster_iso_04_hcalout_{node}[n]` | `Float_t` | OHCal ISO, R=0.4 |
| `cluster_iso_03_sub1_emcal_{node}[n]` | `Float_t` | Background-subtracted EMCal ISO, R=0.3, minus cluster ET |
| `cluster_iso_03_sub1_hcalin_{node}[n]` | `Float_t` | Background-subtracted IHCal ISO, R=0.3 |
| `cluster_iso_03_sub1_hcalout_{node}[n]` | `Float_t` | Background-subtracted OHCal ISO, R=0.3 |
| `cluster_iso_04_sub1_emcal_{node}[n]` | `Float_t` | Background-subtracted EMCal ISO, R=0.4, minus cluster ET |
| `cluster_iso_04_sub1_hcalin_{node}[n]` | `Float_t` | Background-subtracted IHCal ISO, R=0.4 |
| `cluster_iso_04_sub1_hcalout_{node}[n]` | `Float_t` | Background-subtracted OHCal ISO, R=0.4 |
| `cluster_iso_topo_03_{node}[n]` | `Float_t` | Topo-cluster-based ISO, R=0.3, minus cluster ET |
| `cluster_iso_topo_04_{node}[n]` | `Float_t` | Topo-cluster-based ISO, R=0.4, minus cluster ET |
| `cluster_iso_topo_soft_03_{node}[n]` | `Float_t` | Soft topo-cluster ISO, R=0.3, minus cluster ET |
| `cluster_iso_topo_soft_04_{node}[n]` | `Float_t` | Soft topo-cluster ISO, R=0.4, minus cluster ET |

#### Shower shape variables

| Branch | Type | Description |
|--------|------|-------------|
| `cluster_e1_{node}[n]` | `Float_t` | `showershape[8]` -- energy in quadrant 1 of 2x2 around center |
| `cluster_e2_{node}[n]` | `Float_t` | `showershape[9]` -- energy in quadrant 2 |
| `cluster_e3_{node}[n]` | `Float_t` | `showershape[10]` -- energy in quadrant 3 |
| `cluster_e4_{node}[n]` | `Float_t` | `showershape[11]` -- energy in quadrant 4 |
| `cluster_et1_{node}[n]` | `Float_t` | `showershape[0]` -- ET fraction in quadrant 1 |
| `cluster_et2_{node}[n]` | `Float_t` | `showershape[1]` -- ET fraction in quadrant 2 |
| `cluster_et3_{node}[n]` | `Float_t` | `showershape[2]` -- ET fraction in quadrant 3 |
| `cluster_et4_{node}[n]` | `Float_t` | `showershape[3]` -- ET fraction in quadrant 4 |
| `cluster_ietacent_{node}[n]` | `Float_t` | `showershape[4]` -- continuous eta index of cluster center |
| `cluster_iphicent_{node}[n]` | `Float_t` | `showershape[5]` -- continuous phi index of cluster center |
| `cluster_weta_{node}[n]` | `Float_t` | Second moment in eta (integer distance from center tower, owned towers only) |
| `cluster_wphi_{node}[n]` | `Float_t` | Second moment in phi (integer distance) |
| `cluster_weta_cog_{node}[n]` | `Float_t` | Second moment in eta with CoG-corrected center |
| `cluster_wphi_cog_{node}[n]` | `Float_t` | Second moment in phi with CoG-corrected center |
| `cluster_weta_cogx_{node}[n]` | `Float_t` | Second moment in eta with CoG, excluding seed tower |
| `cluster_wphi_cogx_{node}[n]` | `Float_t` | Second moment in phi with CoG, excluding seed tower |
| `cluster_detamax_{node}[n]` | `Int_t` | Max delta-eta extent of cluster (in tower units from lead tower) |
| `cluster_dphimax_{node}[n]` | `Int_t` | Max delta-phi extent of cluster |
| `cluster_nsaturated_{node}[n]` | `Int_t` | Number of saturated towers in cluster |

#### Energy sums in NxM tower windows (centered on cluster CoG)

| Branch | Type | Description |
|--------|------|-------------|
| `cluster_e11_{node}[n]` | `Float_t` | 1x1 (seed tower energy from E77 matrix) |
| `cluster_e22_{node}[n]` | `Float_t` | 2x2 (from showershape: sum of e1+e2+e3+e4) |
| `cluster_e13_{node}[n]` | `Float_t` | 1x3 (1 eta x 3 phi towers) |
| `cluster_e15_{node}[n]` | `Float_t` | 1x5 |
| `cluster_e17_{node}[n]` | `Float_t` | 1x7 |
| `cluster_e31_{node}[n]` | `Float_t` | 3x1 (3 eta x 1 phi) |
| `cluster_e51_{node}[n]` | `Float_t` | 5x1 |
| `cluster_e71_{node}[n]` | `Float_t` | 7x1 |
| `cluster_e33_{node}[n]` | `Float_t` | 3x3 |
| `cluster_e35_{node}[n]` | `Float_t` | 3x5 |
| `cluster_e37_{node}[n]` | `Float_t` | 3x7 |
| `cluster_e53_{node}[n]` | `Float_t` | 5x3 |
| `cluster_e73_{node}[n]` | `Float_t` | 7x3 |
| `cluster_e55_{node}[n]` | `Float_t` | 5x5 |
| `cluster_e57_{node}[n]` | `Float_t` | 5x7 |
| `cluster_e75_{node}[n]` | `Float_t` | 7x5 |
| `cluster_e77_{node}[n]` | `Float_t` | 7x7 (full matrix sum) |
| `cluster_w32_{node}[n]` | `Float_t` | Second moment in eta for 3x2 strip (3 eta x 2 phi) |
| `cluster_e32_{node}[n]` | `Float_t` | Energy sum in 3x2 strip |
| `cluster_w52_{node}[n]` | `Float_t` | Second moment in eta for 5x2 strip |
| `cluster_e52_{node}[n]` | `Float_t` | Energy sum in 5x2 strip |
| `cluster_w72_{node}[n]` | `Float_t` | Second moment in eta for 7x2 strip |
| `cluster_e72_{node}[n]` | `Float_t` | Energy sum in 7x2 strip |

#### 7x7 Tower Arrays

| Branch | Type | Description |
|--------|------|-------------|
| `cluster_e_array_{node}[n][49]` | `Float_t[49]` | Calibrated tower energies in 7x7 grid (row-major, eta outer) |
| `cluster_adc_array_{node}[n][49]` | `Float_t[49]` | Raw ADC tower energies in 7x7 grid |
| `cluster_time_array_{node}[n][49]` | `Float_t[49]` | Tower times in 7x7 grid |
| `cluster_e_array_idx_{node}[n][49]` | `Int_t[49]` | Tower info keys for each 7x7 position |
| `cluster_status_array_{node}[n][49]` | `Int_t[49]` | Tower status flags |
| `cluster_ownership_array_{node}[n][49]` | `Int_t[49]` | Whether tower belongs to this cluster (1) or not (0) |

#### HCal energy behind cluster

| Branch | Type | Description |
|--------|------|-------------|
| `cluster_ihcal_et_{node}[n]` | `Float_t` | IHCal ET in closest tower (1x1) |
| `cluster_ohcal_et_{node}[n]` | `Float_t` | OHCal ET in closest tower (1x1) |
| `cluster_ihcal_et22_{node}[n]` | `Float_t` | IHCal ET in 2x2 window behind cluster |
| `cluster_ohcal_et22_{node}[n]` | `Float_t` | OHCal ET in 2x2 window behind cluster |
| `cluster_ihcal_et33_{node}[n]` | `Float_t` | IHCal ET in 3x3 window behind cluster |
| `cluster_ohcal_et33_{node}[n]` | `Float_t` | OHCal ET in 3x3 window behind cluster |
| `cluster_ihcal_ieta_{node}[n]` | `Int_t` | IHCal matched tower ieta |
| `cluster_ihcal_iphi_{node}[n]` | `Int_t` | IHCal matched tower iphi |
| `cluster_ohcal_ieta_{node}[n]` | `Int_t` | OHCal matched tower ieta |
| `cluster_ohcal_iphi_{node}[n]` | `Int_t` | OHCal matched tower iphi |

### 6.5 Jet Branches (per jet container {jetnode})

Default: `AntiKt_unsubtracted_r04`. Array dimension: `njet_{jetnode}` (max 1000).

| Branch | Type | Description |
|--------|------|-------------|
| `njet_{jetnode}` | `Int_t` | Number of jets |
| `jet_E_{jetnode}[n]` | `Float_t` | Jet energy |
| `jet_Et_{jetnode}[n]` | `Float_t` | Jet transverse energy |
| `jet_Pt_{jetnode}[n]` | `Float_t` | Jet transverse momentum |
| `jet_Eta_{jetnode}[n]` | `Float_t` | Jet pseudorapidity |
| `jet_Phi_{jetnode}[n]` | `Float_t` | Jet azimuthal angle |
| `jet_time_{jetnode}[n]` | `Float_t` | Energy-weighted mean tower time (-999 if no towers > 0.1 GeV) |
| `jet_emcal_calo_E_{jetnode}[n]` | `Float_t` | EMCal energy contribution to jet |
| `jet_ihcal_calo_E_{jetnode}[n]` | `Float_t` | IHCal energy contribution to jet |
| `jet_ohcal_calo_E_{jetnode}[n]` | `Float_t` | OHCal energy contribution to jet |

### 6.6 Truth Jet Branches (MC only, per container {truthjetnode})

Default: `AntiKt_Truth_r04`. Array dimension: `njet_truth_{truthjetnode}` (max 1000).

| Branch | Type | Description |
|--------|------|-------------|
| `njet_truth_{truthjetnode}` | `Int_t` | Number of truth jets |
| `jet_truth_E_{truthjetnode}[n]` | `Float_t` | Truth jet energy |
| `jet_truth_Et_{truthjetnode}[n]` | `Float_t` | Truth jet ET |
| `jet_truth_Pt_{truthjetnode}[n]` | `Float_t` | Truth jet pT |
| `jet_truth_Eta_{truthjetnode}[n]` | `Float_t` | Truth jet eta |
| `jet_truth_Phi_{truthjetnode}[n]` | `Float_t` | Truth jet phi |

### 6.7 Histograms (not in tree)

| Object | Type | Description |
|--------|------|-------------|
| `sim_cross_counting` | `TH1I` | 2-bin histogram: bin 1 = N_Generator_Accepted_Event, bin 2 = N_Processed_Event (MC only, for cross-section normalization) |
| `tracking_radiograph` | `TH3F` | 3D histogram (x,y,z) of shower daughter vertices, weighted by momentum (MC only, for conversion studies) |

## 7. Cuts and Filters

### Event-level cuts

1. **Specific bad events skipped** (hardcoded):
   - Run 51576, event 2261368
   - Run 51909, event 11450906

2. **Vertex cut (data only)**: `|vertexz| < 200 cm` (configurable via `m_vertex_cut`, default 200.0)

3. **Vertex NaN check**: Events with NaN vertex are skipped (data).

4. **No MBD vertex**: Event skipped if no MBD vertex found (data).

5. **Event saved only if**: at least one cluster passes cuts in any container, OR at least one truth particle is stored (`nparticles > 0`).

### Cluster-level cuts

1. **Minimum cluster ET**: `ET > clusterpTmin` (default 5 GeV, configurable via `set_clusterpTmin()`)

2. **Edge cut on ieta**: Clusters with lead tower `ieta < 3` or `ieta > 92` are skipped (prevents 7x7 matrix from going out of bounds; EMCal has 96 eta bins 0-95)

3. **CNN classification edge cut**: Clusters where lead tower ieta is within 2 of edges (ieta-2 < 0 or ieta+2 >= 96) are skipped for CNN input.

4. **Empty shower shape**: Clusters with `showershape.size() == 0` are skipped.

5. **MC-only**: Clusters with no matched truth primary particle (`maxPrimary == nullptr`) are skipped via `continue`.

### Truth particle cuts

1. **Embedding ID filter**: Only particles with `embed_id >= 1` are kept (filters out background particles in embedding).

2. **Eta cut**: `|eta| < 1.5` for photon conversion check.

3. **ET cut for conversion check**: `ET > 5 GeV`.

4. **Minimum particle pT**: `pT > particlepTmin` (default 1 GeV) for storing in particle arrays.

5. **Eta cut for storing**: `|eta| < 1.5`.

### Jet cuts

1. **Reco jets**: `pT > 5 GeV` (hardcoded).

2. **Truth jets**: `|eta| < 2`, `5 < pT < 100 GeV`.

### Trigger (data only)

Default trigger bits used: `{12, 13, 22, 24, 25, 26, 27, 30, 31, 36, 37, 38}`. These are stored but **not used as event-level cuts** in CaloAna24 itself -- trigger selection is applied downstream.

## 8. Cluster Variable Computation Details

### Shower shapes (7x7 matrix)

The code constructs a 7x7 tower energy matrix (`E77[7][7]`) centered on the cluster center-of-gravity (CoG), NOT the lead tower:

```
avg_eta = showershape[4] + 0.5   // continuous eta index
avg_phi = showershape[5] + 0.5   // continuous phi index
maxieta = floor(avg_eta)          // integer eta bin for matrix center
maxiphi = floor(avg_phi)          // integer phi bin
```

Towers below the minimum energy threshold (`m_shower_shape_min_tower_E = 0.07 GeV`) are zeroed in the matrix. Only towers flagged `isGood()` contribute.

### Width variables

- **weta, wphi**: Energy-weighted second moment using integer distance from center tower (tower index 3), considering only towers owned by the cluster:
  ```
  weta = sum(E[i][j] * (i-3)^2) / sum(E[i][j])   for owned towers
  wphi = sum(E[i][j] * (j-3)^2) / sum(E[i][j])   for owned towers
  ```

- **weta_cog, wphi_cog**: Same but using floating-point distance from the CoG position:
  ```
  cog_eta = 3 + (avg_eta - floor(avg_eta) - 0.5)
  di_float = i - cog_eta
  weta_cog = sum(E * di_float^2) / sum(E)
  ```

- **weta_cogx, wphi_cogx**: Same as cog but **excluding the seed tower** (i=3, j=3). These are the primary BDT features.

### Strip width variables (w32, w52, w72)

These compute the second moment in eta for a 2-column-wide strip (in phi) centered on the cluster:

```
w32: 3x2 strip -- |deta| <= 1, dphi == 0 or j == (3 + signphi)
w52: 5x2 strip -- |deta| <= 2, same phi logic
w72: 7x2 strip -- |deta| <= 3, same phi logic

signphi = +1 if fractional phi > 0.5, else -1
wN2 = sum(E * (ieta - 3)^2) / sum(E)
```

The phi column selection uses `signphi` to pick the adjacent phi column closest to the CoG.

### Energy ratios (eNM)

Energy sums in NxM tower windows centered on the cluster CoG:
- N = number of towers in eta (total width), M = number of towers in phi
- Computed via: `|di| <= (N-1)/2 && |dj| <= (M-1)/2`
- Example: e33 = 3x3 = towers within |deta|<=1 and |dphi|<=1

### Isolation ET calculation (`calculateET`)

Tower-level isolation energy computation:
```cpp
float calculateET(eta, phi, dR, layer, min_E, use_subtracted=false)
```
- Loops over ALL towers in the specified calorimeter layer
- For each good tower within `deltaR < dR`: if `energy > min_E`, adds `energy / cosh(tower_eta)` to total
- Layer: 0=EMCal, 1=IHCal, 2=OHCal
- Tower eta is corrected for vertex position using `getTowerEta()`
- When `use_subtracted=true`, uses the SUB1 (background-subtracted) tower containers

For EMCal isolation branches, the cluster's own ET is subtracted: `iso_emcal = emcalET - cluster_ET`.

### Topo-cluster isolation (`calculateET_topo`)

Sums ET of topological clusters within `deltaR < dR` of the cluster position. Uses vertex-corrected positions.

### HCal tower matching

`find_closest_hcal_tower()` loops over ALL HCal towers to find the one closest in deltaR to the EMCal cluster. Returns matched (ieta, iphi) plus signs indicating which direction the nearest neighbor lies (for 2x2 selection). The 2x2 window uses these signs to select the adjacent tower quadrant.

## 9. Configuration Dependencies

CaloAna24 itself has minimal external configuration -- most parameters are set via setter methods in the Fun4All macro:

| Parameter | Setter | Default | Description |
|-----------|--------|---------|-------------|
| `isMC` | `set_isMC(bool)` | `true` | MC vs data mode |
| `isSingleParticle` | `set_isSingleParticle(bool)` | `false` | Single-particle gun mode (sets isMC=true, uses truth vertex) |
| `clusterpTmin` | `set_clusterpTmin(float)` | 5.0 GeV | Minimum cluster ET to store |
| `particlepTmin` | `set_particlepTmin(float)` | 1.0 GeV | Minimum truth particle pT to store |
| `using_trigger_bits` | `set_using_trigger_bits(vector)` | {12,13,22,24,25,26,27,30,31,36,37,38} | Trigger bits to record |
| jet node name | `set_jet_node_name(string)` | `AntiKt_unsubtracted_r04` | Reco jet container |
| truth jet node name | `set_truth_jet_node_name(string)` | `AntiKt_Truth_r04` | Truth jet container |
| output file name | constructor arg | `caloana.root` | Output ROOT file |

The Fun4All macros also configure:
- CDB global tag: `ProdA_2024` (data) or `MDC2` (sim)
- RUNNUMBER/TIMESTAMP: from DST filename (data) or 28 (sim)
- Cluster builder threshold: `0.070 GeV`
- ClusterIso min pT: 5 GeV, cone radii: 2, 3, 4 (in units of 0.1)
- Topo cluster noise thresholds: EMCal 0.0053, IHCal 0.0351, OHCal 0.0684 (3-sigma pedestal)
- Jet reconstruction: anti-kT R=0.4, using retowered EMCal + IHCal + OHCal

## 10. Hardcoded Constants

| Constant | Value | Location | Description |
|----------|-------|----------|-------------|
| `nclustermax` | 10000 | CaloAna24.h | Max clusters per container per event |
| `nparticlesmax` | 10000 | CaloAna24.h | Max truth particles per event |
| `ndaughtermax` | 100 | CaloAna24.h | Max daughter particles (single-particle mode) |
| `nclustercontainer` | 2 | CaloAna24.h | Number of cluster containers (NO_SPLIT + CEMC) |
| `arrayntower` | 49 | CaloAna24.h | 7x7 tower array size |
| `njetmax` | 1000 | CaloAna24.h | Max jets per container |
| `njet_truthmax` | 1000 | CaloAna24.h | Max truth jets per container |
| `njetcontainermax` | 8 | CaloAna24.h | Max jet container slots |
| `nRadii` | 3 | CaloAna24.h | Number of isolation radii (R=0.2, 0.3, 0.4) |
| `m_vertex_cut` | 200.0 cm | CaloAna24.h | Vertex z cut (data only) |
| `m_shower_shape_min_tower_E` | 0.07 GeV | CaloAna24.h | Min tower energy for shower shape computation |
| `truthisocut` | 4 GeV | CaloAna24.h | Truth isolation cut (declared but not used in CaloAna24 itself) |
| `clusterpTmin` | 5.0 GeV | CaloAna24.h | Default min cluster ET |
| `particlepTmin` | 1.0 GeV | CaloAna24.h | Default min truth particle pT |
| EMCal eta bins | 96 (0-95) | process_event | Used for edge checks, wrapping |
| EMCal phi bins | 256 (0-255) | process_event | Used for phi wrapping |
| IHCal eta bins | 24 | process_event | Used for shift_tower_index |
| IHCal phi bins | 64 | process_event | Used for phi wrapping |
| OHCal eta/phi bins | 24 / 64 | process_event | Same as IHCal |
| MBD channels | 128 | process_event | 64 north + 64 south |
| MBD charge threshold | 0.4 | process_event | Minimum charge for hit counting |
| Conversion vertex radius | 93 cm | process_event | Max radius for pre-CEMC interaction check |
| Conversion momentum fraction | 0.4 | process_event | Daughter momentum > 40% of parent = "bad" photon |
| pi0 mass | 0.1349799931 GeV | process_event | For invariant mass matching (currently in commented-out block) |
| eta mass | 0.5478500128 GeV | process_event | For invariant mass matching (commented out) |
| merger cone | 0.001 | process_event | Cone for self-ET subtraction in truth isolation |
| CNN input dim | 5x5 | process_event | CNN classifier input (disabled) |
| Bad event: run 51576, evt 2261368 | -- | process_event | Hardcoded bad event skip |
| Bad event: run 51909, evt 11450906 | -- | process_event | Hardcoded bad event skip |
| Jet min pT (reco) | 5 GeV | process_event | Hardcoded jet pT cut |
| Truth jet eta range | [-2, 2] | process_event | Hardcoded truth jet eta cut |
| Truth jet pT range | [5, 100] GeV | process_event | Hardcoded truth jet pT cut |
| ONNX model path | `/sphenix/u/shuhang98/core_patch/...` | CaloAna24.h | CNN model path (unused -- inference disabled) |

## 11. MC vs Data Differences

| Aspect | Data | MC |
|--------|------|-----|
| Trigger processing | GL1 packet decoded, scaler/prescale recorded | Skipped entirely |
| Vertex source | MbdVertexMap (data-specific code path) | GlobalVertexMap with MBD vertex type; single-particle mode uses truth vertex |
| Vertex cut | `|z| < 200 cm` enforced | No vertex cut applied |
| MC cross-section | -- | `PHGenIntegral` node read for N_accepted and N_processed (fills `sim_cross_counting`) |
| HepMC record | -- | Pythia process ID, energy scale extracted; photon classification via `photon_type()` |
| Truth matching | -- | `CaloEvalStack` / `CaloRawClusterEval` used to match clusters to truth primaries |
| Truth particles | `nparticles` = 0 | Primary particles stored with isolation, conversion status |
| Daughters | `ndaughter` = 0 | Stored only in single-particle mode |
| Truth jets | Not stored | `AntiKt_Truth_r04` jets stored |
| Jet container | `AntiKt_unsubtracted_r04` | `AntiKt_unsubtracted_r04_calib` (calibrated, sim only) |
| Cluster skip | -- | Clusters with no truth match (`maxPrimary == nullptr`) are skipped |
| `isMC` flag | `false` | `true` |
| CDB tag | `ProdA_2024` | `MDC2` |
| Run number | From DST filename | Hardcoded to 28 |
| Pedestal overlay | -- | Random pedestals from run 54256 overlaid (sim Fun4All macro) |
| Waveform sim | -- | Full waveform chain exists but is disabled (`if(false)` block) |
| Cluster builder | Runs from DST; both CEMC nodes built | Same, but with `setMinTowerEnergy(0.070)` in ClusterIso |
| MBD reco | Not run (assumed in DST) | `MbdReco` explicitly registered |
| GlobalVertexReco | Not run | Explicitly registered |
| JetCalib | Not run | Applied to produce calibrated jet output |

## 12. Gotchas and Non-Obvious Details

1. **Two cluster containers always stored**: Every cluster branch exists twice -- once for `CLUSTERINFO_CEMC_NO_SPLIT` and once for `CLUSTERINFO_CEMC`. Downstream code must specify which node to use. The current pipeline uses `CLUSTERINFO_CEMC`.

2. **CNN probability always -1**: The ONNX inference is fully disabled (both the session creation and the inference call are commented out). `cluster_CNN_prob` is hardcoded to -1.

3. **7x7 matrix centered on CoG, not lead tower**: The shower-shape matrix uses `floor(showershape[4] + 0.5)` and `floor(showershape[5] + 0.5)` as center, which is the cluster CoG index, NOT the lead tower index used earlier for CNN input and detamax/dphimax computation. These are different in general.

4. **MC clusters without truth match are silently skipped**: The `if (!maxPrimary) continue;` on line 1213 means MC clusters that cannot be matched to a truth primary are not stored at all. This is intentional but means MC cluster counts may be lower than data for the same events.

5. **Trigger prescale only fills first 32 bits**: The loop `for (int i = 0; i < 32; i++)` fills `trigger_prescale[0..31]`, but the array is 64 elements. Elements 32-63 remain at default (-1).

6. **Vertex handling differs between data and MC**: Data uses `MbdVertexMap` directly. MC uses `GlobalVertexMap` and iterates through MBD-typed vertices. This means the vertex value can differ even for the same underlying data if the vertex maps are populated differently.

7. **sim_cross_counting only filled at event 1**: The PHGenIntegral node is read only when `m_eventnumber == 1`. This means only the first segment's cross-section info is captured. For multi-segment jobs, only the first segment contributes.

8. **Two hardcoded bad events**: Run 51576 event 2261368 and run 51909 event 11450906 are unconditionally skipped. The reason is not documented in the code.

9. **e22 computed differently from other eNM**: While e33, e55, e77 etc. are computed from the E77 matrix, e22 is computed from `showershape[8]+[9]+[10]+[11]` (the four quadrant energies from `RawCluster::get_shower_shapes()`). This is consistent but uses a different code path.

10. **Shower shape minimum tower energy is hardcoded to 0.07 GeV**: Both `RawClusterBuilderTemplate` threshold and `m_shower_shape_min_tower_E` use 0.07 GeV, but they are set independently. Changing one without the other would create inconsistency.

11. **Isolation ET subtracts cluster ET for EMCal only**: The HCal isolation branches do NOT subtract the cluster ET. EMCal branches do: `cluster_iso_03_emcal = emcalET_03 - ET`. This is because the cluster is in the EMCal only.

12. **calculateET loops over ALL towers**: The isolation function iterates every tower in the calorimeter for every cluster. This is O(N_towers * N_clusters) per event.

13. **find_closest_hcal_tower loops over ALL HCal towers**: Same brute-force approach for HCal matching. No spatial indexing.

14. **Phi wrapping in 7x7 matrix vs tower map**: The `shift_tower_index()` function handles phi wrapping (mod 256 for EMCal, mod 64 for HCal) but sets eta to -1 for out-of-range values. When eta becomes -1, the tower key encode may produce unexpected results, but the subsequent `ieta < 0` check catches this.

15. **Jet component type IDs are magic numbers**: Jet tower components are identified by integer IDs (14/29/25/28 = EMCal, 15/30/26 = IHCal, 16/31/27 = OHCal). These correspond to `Jet::SRC` enum values in the sPHENIX codebase but are hardcoded as integers.

16. **Waveform simulation block disabled**: The full tower-level simulation (CaloWaveformSim + CaloTowerBuilder) in the sim macro is wrapped in `if(false)`. The DSTs are assumed to already contain calibrated tower info. The MC recalibration (`CEMC_MC_RECALIB`) is also inside this disabled block.

17. **Pedestal overlay in sim**: The sim macro adds real pedestal data from run 54256 via `Fun4AllNoSyncDstInputManager` with `Repeat()`. This adds realistic noise to the simulation.

18. **Single-particle mode overrides vertex**: When `isSingleParticle = true`, `vertexz` is overwritten with `vertexz_truth`. This ensures single-particle efficiency studies use the truth vertex.

19. **photon_type classification follows Pythia shower history**: The `photon_type()` function traces back through HepMC vertices, skipping "pass-through" vertices (single photon in, single photon out). Classification: 2->2 with quarks/gluons = direct (class 1), 1->2 with quark split = fragmentation (class 2), incoming hadron (pid > 37) = decay (class 3).

20. **Conversion tagging radius = 93 cm**: Daughters of truth photons that interact within r < 93 cm (EMCal inner radius is ~93.5 cm) and carry > 40% of the photon energy are flagged. Electron/positron pairs get `converted=1`, other species get `converted=2`.

21. **Data production versions**: The directory tree shows evolution through ana450, ana462, ana468, ana484, ana502, ana509, ana518, to ana521 (current). Each represents a different sPHENIX production tag with updated calibrations.

22. **sim create_file_lists.sh type mapping**: The MC production uses specific type IDs: jet5=36, jet12=39, jet20=21, jet30=11, jet40=19, jet50=34, photon5=27, photon10=28, photon20=29. These are sPHENIX simulation production type codes.
