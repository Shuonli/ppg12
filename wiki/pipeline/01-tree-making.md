# Stage 1: Tree Making (anatreemaker)

## Key Files

| File | Purpose |
|------|---------|
| `anatreemaker/source/CaloAna24.cc` | Main Fun4All SubsysReco module (2321 lines) |
| `anatreemaker/source/CaloAna24.h` | Class declaration (385 lines) |
| `anatreemaker/macro_maketree/data/ana521/Fun4All_runDST.C` | Data processing macro (current production) |
| `anatreemaker/macro_maketree/sim/run28/{sample}/Fun4All_run_sim.C` | Sim processing macro |
| `anatreemaker/macro_maketree/data/ana521/run.sh` | Condor orchestration for data |
| `anatreemaker/macro_maketree/sim/run28/hadd_combined.sh` | Merge per-job sim outputs |

## Building

```bash
cd anatreemaker/source
./autogen.sh --prefix=$MYINSTALL
cd build && make install  # produces libCaloAna24.so
```

## Running

**Data:**
```bash
cd anatreemaker/macro_maketree/data/ana521
# Edit runList.txt with desired run numbers
./run.sh  # submits condor jobs
```
Output: `OUTTREE_DST_Jet_*.root` per job. Production tag: `ana521_2025p007_v001`.

**Simulation:**
```bash
cd anatreemaker/macro_maketree/sim/run28/photon10
./run_condor.sh  # submits condor jobs
# Merge outputs:
cd .. && ./hadd_combined.sh  # produces condorout/combined.root per sample
```

## Output Tree: `slimtree`

The tree name `slimtree` is hardcoded. All cluster branches are duplicated for two containers: `CLUSTERINFO_CEMC` (with splitting, the default) and `CLUSTERINFO_CEMC_NO_SPLIT` (without splitting, legacy).

### Key Event-Level Branches

| Branch | Type | Notes |
|--------|------|-------|
| `runnumber` | Int_t | |
| `vertexz` | Float_t | MBD vertex (data), GlobalVertex (MC) |
| `vertexz_truth` | Float_t | MC only, -9999 for data |
| `scaledtrigger[64]` | Bool_t[] | Trigger decisions (data only) |
| `mbdnorthhit`, `mbdsouthhit` | Int_t | MBD hit multiplicities |
| `mbdnorth{q,t}[64]`, `mbdsouth{q,t}[64]` | Float_t[] | MBD PMT charges and times |
| `mbd_time` | Float_t | MBD t0 |
| `pythiaid` | Int_t | Pythia process ID (MC only) |
| `energy_scale` | Float_t | Event Q^2 (MC only) |

### Key Cluster Branches (per `{node}`)

| Branch | Type | Notes |
|--------|------|-------|
| `ncluster_{node}` | Int_t | Max 10000 |
| `cluster_Et_{node}[n]` | Float_t | E / cosh(eta) |
| `cluster_Eta_{node}[n]` | Float_t | Vertex-corrected |
| `cluster_et1-4_{node}[n]` | Float_t | Ring energy fractions |
| `cluster_weta_cogx_{node}[n]` | Float_t | Width in eta, CoG, excl. seed |
| `cluster_wphi_cogx_{node}[n]` | Float_t | Width in phi, CoG, excl. seed |
| `cluster_e{NM}_{node}[n]` | Float_t | NxM tower window energies (e11 through e77) |
| `cluster_w{32,52,72}_{node}[n]` | Float_t | Strip width variables |
| `cluster_iso_03_{node}[n]` | Float_t | RawCluster isolation R=0.3 |
| `cluster_iso_topo_{03,04}_{node}[n]` | Float_t | Topo-cluster isolation |
| `cluster_iso_03_60_{emcal,hcalin,hcalout}_{node}[n]` | Float_t | Threshold isolation components |
| `cluster_truthtrkID_{node}[n]` | Int_t | Truth matching (MC only) |
| `cluster_e_array_{node}[n*49]` | Float_t | 7x7 tower energy grid |

### Truth Particle Branches (MC only)

| Branch | Type | Notes |
|--------|------|-------|
| `nparticles` | Int_t | Max 10000 |
| `particle_photonclass[n]` | Int_t | 1=direct, 2=fragmentation, 3=decay, 0=unid, -1=not photon |
| `particle_truth_iso_{02,03,04}[n]` | Float_t | Truth isolation in cone R=0.2/0.3/0.4 |
| `particle_converted[n]` | Int_t | 0=no, 1=pair-converted, 2=other interaction before CEMC |

### Non-Tree Output

| Object | Type | Notes |
|--------|------|-------|
| `sim_cross_counting` | TH1I | Bin 1 = N_Generator_Accepted, bin 2 = N_Processed (MC only) |

## Reconstruction Chain (registered before CaloAna24)

1. `Process_Calo_Calib()` -- standard tower calibration
2. `RawClusterBuilderTemplate` -- EMCal clustering (threshold 0.070 GeV), two containers
3. `ClusterIso` -- cone isolation R=0.2, 0.3, 0.4
4. `RawClusterBuilderTopo` -- topological clustering (standard + soft thresholds)
5. `RetowerCEMC` -- retower EMCal to HCal granularity
6. `JetReco` + `FastJetAlgoSub` -- anti-kT R=0.4 jet reconstruction

Sim-only additions: `MbdReco`, `GlobalVertexReco`, `JetCalib`.

## Shower Shape Computation

The 7x7 tower matrix (`E77[7][7]`) is centered on the cluster center-of-gravity (CoG), NOT the lead tower. Towers below `m_shower_shape_min_tower_E = 0.07 GeV` are zeroed. The `weta_cogx` and `wphi_cogx` variables use the CoG-corrected center and **exclude the seed tower** -- these are the primary BDT features.

## Cuts Applied in CaloAna24

- Minimum cluster ET: `> clusterpTmin` (default 5 GeV)
- Edge cut: lead tower ieta in [3, 92] (prevents 7x7 overflow; EMCal has 96 eta bins)
- Vertex: `|vertexz| < 200 cm` (data only)
- Bad events: run 51576/event 2261368, run 51909/event 11450906

## Gotchas

- Cluster branches appear for BOTH `CLUSTERINFO_CEMC` and `CLUSTERINFO_CEMC_NO_SPLIT` containers. Downstream code must select the correct one via config `input.cluster_node_name`.
- CNN classification (`cluster_CNN_prob`) is always -1 (ONNX runtime disabled).
- The `hadd_combined.sh` must be run manually after condor jobs complete.
- Production tags matter: `ana521` is current. Older directories (`ana450`, `ana462`, etc.) are legacy.
