# Double Interaction / Pileup Pipeline (PPG12)

## Physics Motivation

In pp collisions at RHIC, two separate collisions can occur in the same bunch crossing (double interaction / pileup). This shifts the reconstructed vertex (average of two collision vertices) and distorts cluster kinematics (ET, eta) and shower-shape BDT scores. Two complementary approaches quantify the effect: full GEANT double-interaction MC samples blended with single-interaction MC, and a toy vertex-shift simulation.

## Two Approaches

### 1. Full GEANT Double-Interaction MC Blending (`run_showershape_double.sh`)

Uses actual double-interaction MC samples produced through full detector simulation, blended with single-interaction MC at proper physics fractions.

**Double interaction fractions** (cluster-weighted, from `calc_pileup_range.C`; triple+ folded into double):
- **0 mrad crossing angle**: 22.4% cluster-weighted (`DOUBLE_FRAC=0.224`, event-level `f_double`=11.1%)
- **1.5 mrad crossing angle**: 7.9% cluster-weighted (`DOUBLE_FRAC=0.079`, event-level `f_double`=3.9%)

**MC samples** (in `anatreemaker/macro_maketree/sim/run28/`):
- `photon10_double/` — double-interaction photon (14-30 GeV truth pT)
- `jet12_double/` — double-interaction jet (14-21 GeV pT)
- Corresponding nominal: `photon10/` (aka `photon10_nom`), `jet12/`

**Two-pass pipeline** (`run_showershape_double.sh`):

```
Pass 1 (vertex scan): ShowerShapeCheck(config, sample, ..., do_vertex_scan=true, mix_weight=frac)
  → Fills only h_vertexz weighted by physics fraction
  → 5 parallel jobs: photon10_double(0.224), jet12_double(0.224), photon10_nom(0.776), jet12_nom(0.776), data

hadd: nom_vtxscan + double_vtxscan → combined_vtxscan (blended MC vertex distribution)

Pass 2 (full analysis): ShowerShapeCheck(config, sample, ..., mix_weight=frac, vtxscan_sim_override=combined)
  → Full ABCD/efficiency/purity analysis
  → All MC samples use the SAME blended vtxscan for vertex reweighting
  → mix_weight applied to all histogram fills

hadd: nom + double → combined outputs
```

**ShowerShapeCheck.C signature**:
```cpp
void ShowerShapeCheck(configname, filetype, doinclusive, do_vertex_scan, mix_weight, vtxscan_sim_override)
```
- `mix_weight` — multiplicative fraction for all histogram fills (0.224 or 0.776 for 0 mrad); `cross_weight *= mix_weight`
- `vtxscan_sim_override` — path to blended vtxscan ROOT file for consistent vertex reweighting
- `do_vertex_scan` — if true, only fills `h_vertexz` (fast Pass 1)

**How to run**:
```bash
cd efficiencytool
bash run_showershape_double.sh config_showershape.yaml        # 0 mrad (22.4%)
bash run_showershape_double.sh config_showershape.yaml 0.079  # 1.5 mrad (7.9%)
```

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
| `efficiencytool/run_showershape_double.sh` | Full GEANT blending pipeline (2-pass) |
| `efficiencytool/ShowerShapeCheck.C` | Shower shape analysis (supports mix_weight + vtxscan override) |
| `efficiencytool/DoubleInteractionCheck.C` | Toy vertex-shift simulation (~1535 lines) |
| `efficiencytool/run_double_interaction.sh` | Toy simulation orchestration |
| `efficiencytool/MbdPileupHelper.h` | MBD pileup metric calculator |
| `efficiencytool/NPB_PurityStudy.C` | NPB purity study (0mrad/1.5mrad split) |
| `efficiencytool/CrossSectionWeights.h` | MC cross-section weights (includes jet5/8/12/40) |
| `plotting/plot_double_interaction.C` | Plotting macro for toy sim results |
| `anatreemaker/macro_maketree/sim/run28/photon10_double/` | Full GEANT double-interaction MC production |
| `anatreemaker/macro_maketree/sim/run28/jet12_double/` | Full GEANT double-interaction MC production |

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

**Full GEANT blending** (from `run_showershape_double.sh`):
- `results/MC_efficiencyshower_shape_{sample}_{suffix}_vtxscan.root` (Pass 1)
- `results/MC_efficiencyshower_shape_{sample}_combined[_inclusive]_{suffix}.root` (Pass 2, merged)

**Toy simulation** (from `run_double_interaction.sh`):
- `results/{eff_outfile}_double_interaction_check_{filetype}[_inclusive].root`
- `results/{data_outfile}_double_interaction_check.root`
- Histogram naming: `h_{tight/nontight}_{iso/noniso}_cluster_{single/double}[_smearN][_signal/_incsig]_{etabin}`
