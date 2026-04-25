# Double Interaction / Pileup Pipeline (PPG12)

## Physics Motivation

In pp collisions at RHIC, two separate collisions can occur in the same bunch crossing (double interaction / pileup). This shifts the reconstructed vertex (average of two collision vertices) and distorts cluster kinematics (ET, eta) and shower-shape BDT scores. Two complementary approaches quantify the effect: full GEANT double-interaction MC samples blended with single-interaction MC, and a toy vertex-shift simulation.

## Two Approaches

### 1. Full GEANT Double-Interaction MC Blending (condor single-pass truth-vertex)

Uses actual double-interaction MC samples produced through full detector simulation, blended with single-interaction MC at proper physics fractions.

**Double interaction fractions** (cluster-weighted, from `calc_pileup_range.C`; triple+ folded into double):
- **0 mrad crossing angle**: 22.4% cluster-weighted (`DOUBLE_FRAC=0.224`, event-level `f_double`=11.1%)
- **1.5 mrad crossing angle**: 7.9% cluster-weighted (`DOUBLE_FRAC=0.079`, event-level `f_double`=3.9%)

**MC samples** (in `anatreemaker/macro_maketree/sim/run28/`) — **8 SI/DI pairs**:

| channel | SI | DI | truth pT |
|---|---|---|---|
| signal | `photon5` | `photon5_double` | 0-14 GeV |
| signal | `photon10` | `photon10_double` | 14-30 GeV |
| signal | `photon20` | `photon20_double` | 30+ GeV |
| BG | `jet8` | `jet8_double` | 9-14 GeV |
| BG | `jet12` | `jet12_double` | 14-21 GeV |
| BG | `jet20` | `jet20_double` | 21-32 GeV |
| BG | `jet30` | `jet30_double` | 32-42 GeV |
| BG | `jet40` | `jet40_double` | 42+ GeV |

SI aliased to `{sample}_nom` in analysis outputs to distinguish from pure-SI runs. **Skipped**: `jet5` (no DI partner exists), `jet50_double` (no SI partner in mc_samples.list).

**Single-pass truth-vertex reweighting pipeline** (primary; replaces prior two-pass reco-vertex scheme):

1. `ShowerShapeCheck.C` reads per-event truth vertex (hard-scatter `vertexz_truth`, and for DI also the min-bias `vertexz_truth_mb`) and applies `TruthVertexWeight(h, z_hard, z_mb)` factorized reweight from `TruthVertexReweightLoader.h`. The reweight histogram `h_w_iterative` is precomputed in `efficiencytool/truth_vertex_reweight/output/{0mrad,1p5mrad}/reweight.root`.
2. Per-event weight = `cross_section × mix_weight × TruthVertexWeight(...)`. `mix_weight = DOUBLE_FRAC` for `_double` samples and `SINGLE_FRAC = 1 - DOUBLE_FRAC` for `_nom` samples.
3. Config flag: `analysis.truth_vertex_reweight_on: 1` + `analysis.truth_vertex_reweight_file: ...reweight.root`. Mutually exclusive with the legacy `vertex_reweight_on` (reco-vertex), which is silently forced off when the truth path is on.
4. `ShowerShapeCheck.C` uses `filetype.compare(size-7, 7, "_double") == 0` to detect DI samples for the R=0.4 truth-jet branch selection. `_nom` aliases resolve to the SI `combined.root` at `anatreemaker/macro_maketree/sim/run28/{SI}/condorout/combined.root`.

**Condor packaging** (17 jobs per crossing-angle config: 8 `_nom` × 1 SI + 8 `_double` × 1 DI + 1 data):

```bash
cd efficiencytool
condor_submit list_file=showershape_di_jobs_0rad.list   submit_showershape_di.sub  # 0 mrad
condor_submit list_file=showershape_di_jobs_1p5rad.list submit_showershape_di.sub  # 1.5 mrad
bash hadd_showershape_di.sh   # final per-channel merge (signal + inclusive BG)
```

**ShowerShapeCheck.C signature** (still takes `do_vertex_scan` and `vtxscan_sim_override` for legacy compat, but the DI pipeline passes `do_vertex_scan=false` and empty override; the reweight path is driven entirely by config):
```cpp
void ShowerShapeCheck(configname, filetype, doinclusive, do_vertex_scan, mix_weight, vtxscan_sim_override)
```

**Legacy two-pass reco-vertex pipeline** (`run_showershape_double_reco_legacy.sh`, renamed from `run_showershape_double.sh`): preserved for reference/cross-check. Filled `h_vertexz` in Pass 1, used the blended MC reco-vertex distribution as the reweight reference in Pass 2. Retired in favour of the truth-vertex single-pass flow above.

### 2. Toy Double-Interaction Simulation (`DoubleInteractionCheck.C`)

Takes existing single-interaction MC events and artificially shifts the vertex to simulate pileup effects. Used for studying purity degradation and shower-shape distortion.

**How to run**:
```bash
cd efficiencytool
bash run_double_interaction.sh config_showershape.yaml
cd ../plotting && root -l -b -q 'plot_double_interaction.C("config_showershape.yaml")'
```

**Toy algorithm** (sim only):
1. Draw random second vertex from unrestricted data vertex distribution (`vertex_scan_data_file`)
2. Compute double vertex: `double_vtx = (original_vtx + random_vtx) / 2.0`
3. Recalculate cluster kinematics using CEMC cylinder geometry (radius = 93.5 cm):
   - `new_z/r = old_z/r + (old_vtx - new_vtx) / 93.5`
   - `new_eta = asinh(new_z/r)`, `new_ET = E / cosh(new_eta)`
4. Re-evaluate NPB and shower-shape BDT scores from scratch at new vertex/ET/eta
5. Re-classify tight/non-tight and iso/non-iso under shifted kinematics
6. Additional vertex smearing: Gaussian 5, 10, 15, 20 cm on top of `double_vtx`

**MC samples** (finer jet binning than main pipeline):
- Signal: photon5 (0-14 GeV), photon10 (14-30), photon20 (30+)
- Background: jet5 (7-9), jet8 (9-14), jet12 (14-21), jet20 (21-32), jet30 (32-42), jet40 (42+)

**4 analysis tasks**:
1. ABCD yields as 2D (ET vs NPB score) — pileup contamination
2. MBD avg sigma_t vs run number (data only) — pileup rate evolution by crossing angle
3. ABCD yields as 2D (ET vs MBD avg sigma) — data pileup metric vs ABCD
4. Single vs double interaction ABCD (sim only) — truth purity comparison

**Shower shape distributions**: Recorded per (interaction_type, variable, eta_bin, pT_bin) for 11 variables (weta_cogx, wphi_cogx, wr_cogx, et1-4, e11/e33, e32/e35, bdt_score, npb_score) and 6 interaction types (single, double, double_smear5/10/15/20). Data: split by data_0mrad, data_1p5mrad.

## Key Files

| File | Purpose |
|------|---------|
| `efficiencytool/submit_showershape_di.sub` | Condor submit for single-pass truth-vertex DI blending (primary) |
| `efficiencytool/run_showershape_di_job.sh` | Per-job wrapper invoked by condor (4 args: config, sample, doinc, mix_weight) |
| `efficiencytool/showershape_di_jobs_{0rad,1p5rad}.list` | 17-row job list per crossing-angle config |
| `efficiencytool/hadd_showershape_di.sh` | Per-channel merge (signal_combined + jet_inc_combined) |
| `efficiencytool/run_showershape_double_reco_legacy.sh` | Legacy two-pass reco-vertex pipeline (retired reference) |
| `efficiencytool/ShowerShapeCheck.C` | Shower shape analysis (single-pass truth-vertex when `truth_vertex_reweight_on: 1`) |
| `efficiencytool/TruthVertexReweightLoader.h` | Factorized reweight loader `TruthVertexWeight(h, z_hard, z_mb)` |
| `efficiencytool/truth_vertex_reweight/` | Iterative reweight fit (`fit_truth_vertex_reweight.py`, precomputed `reweight.root` per period) |
| `efficiencytool/DoubleInteractionCheck.C` | Toy vertex-shift simulation (~1446 lines) |
| `efficiencytool/run_double_interaction.sh` | Toy simulation orchestration |
| `efficiencytool/MbdPileupHelper.h` | MBD pileup metric calculator |
| `efficiencytool/NPB_PurityStudy.C` | NPB purity study (0mrad/1.5mrad split) |
| `efficiencytool/CrossSectionWeights.h` | MC cross-section weights (all 8 SI/DI pairs + full jet binning) |
| `plotting/plot_double_interaction.C` | Plotting macro for toy sim results |
| `anatreemaker/macro_maketree/sim/run28/photon{5,10,20}_double/` | Full GEANT DI signal MC (3 pT bins) |
| `anatreemaker/macro_maketree/sim/run28/jet{8,12,20,30,40}_double/` | Full GEANT DI BG MC (5 pT bins) |

## MBD Pileup Metrics (MbdPileupHelper.h)

Computes timing spread from 128 MBD PMTs (64 north + 64 south):
- `avgsigma`: average of north/south timing RMS (primary metric, used for pileup rejection via `mbd_avg_sigma_max`)
- `prodsigma`, `maxsigma`, `proddelta`, `avgdelta`, `maxdelta`: additional metrics
- PMT selection: `|time| < 25 ns`, `charge > 0.4`, minimum 2 hits per side

## Config Dependencies

- Uses showershape configs (`config_showershape*.yaml`), same YAML schema as main efficiency configs
- `analysis.vertex_scan_data_file` — unrestricted vertex distribution for toy double interaction
- `analysis.mbd_t0_correction_file` — MBD timing corrections per run (default: `MbdOut.corr`)
- `analysis.npb_tmva_file` — NPB TMVA model (re-scored at shifted vertex in toy sim)
- `input.bdt_et_bin_edges` / `bdt_et_bin_models` — ET-binned BDT models

## Output Files

**Full GEANT blending** (single-pass via `submit_showershape_di.sub`):
- `results/MC_efficiencyshower_shape_{sample}_{suffix}.root` (one per sample in `{photon,jet}{pT}_{nom,double}`)
- After `hadd_showershape_di.sh`: `results/MC_efficiencyshower_shape_signal_combined_{suffix}.root` (6 files merged) and `results/MC_efficiencyshower_shape_jet_inclusive_combined_{suffix}.root` (10 files merged). Cross-section × mix_weight baked into per-event weight, so plain hadd is correct.

**Legacy two-pass output** (from `run_showershape_double_reco_legacy.sh`, retired):
- `results/MC_efficiencyshower_shape_{sample}_{suffix}_vtxscan.root` (Pass 1)
- `results/MC_efficiencyshower_shape_{sample}_combined[_inclusive]_{suffix}.root` (Pass 2, merged)

**Toy simulation** (from `run_double_interaction.sh`):
- `results/{eff_outfile}_double_interaction_check_{filetype}[_inclusive].root`
- `results/{data_outfile}_double_interaction_check.root`
- Histogram naming: `h_{tight/nontight}_{iso/noniso}_cluster_{single/double}[_smearN][_signal/_incsig]_{etabin}`
