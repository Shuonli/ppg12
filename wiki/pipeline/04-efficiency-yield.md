# Stage 4: Efficiency and Yield (efficiencytool)

## Key Files

| File | Lines | Purpose |
|------|-------|---------|
| `efficiencytool/RecoEffCalculator_TTreeReader.C` | 2719 | Main workhorse: cuts, ABCD histograms, efficiencies, response matrices |
| `efficiencytool/MergeSim.C` | 94 | Merges per-sample MC outputs via TFileMerger |
| `efficiencytool/CalculatePhotonYield.C` | 1170 | ABCD purity, unfolding, efficiency correction, cross-section |
| `efficiencytool/CrossSectionWeights.h` | 119 | All MC cross-section constants |
| `efficiencytool/oneforall.sh` | 28 | Orchestrates MergeSim + CalculatePhotonYield |
| `efficiencytool/oneforall_tree.sh` | 39 | Two-pass orchestration for RecoEffCalculator |

## RecoEffCalculator_TTreeReader.C

### Signature

```cpp
void RecoEffCalculator_TTreeReader(
    const std::string &configname,
    const std::string filetype = "jet40",
    bool do_vertex_scan = false,
    float mix_weight = 1.0,
    const std::string vtxscan_sim_override = "")
```

### Two-Pass Architecture

**Pass 1** (`do_vertex_scan=true`): Only fills `h_vertexz`. Output: `{outfile}_vtxscan.root`. Used by Pass 2 to compute `vertex_weight = h_vertexz->Interpolate(vertexz)` for the legacy reco-vertex reweight.

**Pass 2** (`do_vertex_scan=false`): Reads vtxscan files for vertex reweighting, runs full analysis with all cuts and histogram filling.

**Pass 1 auto-skip**: under `analysis.truth_vertex_reweight_on=1` (current production), the reco-vertex reweight branch is force-disabled in `RecoEffCalculator_TTreeReader.C` and the vtxscan output is never read. `oneforall_tree_double.sh` auto-detects this and skips Pass 1 + the vtxscan hadd, halving the per-config wall-clock (~45 min → ~22 min). Pass 2 then receives an empty `vtxscan_sim_override` argument. Legacy configs with `truth_vertex_reweight_on=0` still run both passes.

### Cross-Section Weighting

Uses `PPG12::GetSampleConfig(filetype)` from `CrossSectionWeights.h`:
```cpp
float cross_weight = GetSampleConfig(filetype).weight * mix_weight;
```
Applied to ALL histogram fills as `weight = cross_weight * vertex_weight`.

### Energy Scale/Resolution Smearing (MC only)

```cpp
cluster_Et *= clusterescale;          // analysis.cluster_escale (default 1.0)
cluster_Et *= rand->Gaus(1, clustereres);  // analysis.cluster_eres (default 0.0)
```

### MC Isolation Adjustment

```cpp
recoisoET = recoisoET * mc_iso_scale;  // code default: 1.0, nominal config: 1.2
recoisoET += mc_iso_shift;             // code default: 0.0, nominal config: 0.1
```

### Three Factorized Efficiencies

All computed as TEfficiency objects per eta bin:

1. **`eff_reco_eta_{ieta}`** -- denominator: truth photons (pid=22, photonclass<3, truth_iso < cut); numerator: those with reco cluster match within dR < 0.1
2. **`eff_iso_eta_{ieta}`** -- denominator: reco-matched; numerator: passing reco isolation
3. **`eff_id_eta_{ieta}`** -- denominator: reco + iso; numerator: passing tight BDT/shower-shape cuts

### BDT Model Selection at Runtime

Config specifies `bdt_et_bin_edges` + `bdt_et_bin_models` for ET-dependent model selection:
```yaml
bdt_et_bin_edges: [8, 15, 35]    # N+1 edges
bdt_et_bin_models: [base_v1E, base_v1E]  # N models
```

### Output Files

Per-sample sim: `{eff_outfile}_{filetype}_{var_type}.root`
Data: `{data_outfile}_{var_type}.root`
Response: `{response_outfile}_{filetype}_{var_type}.root`

### Key Output Histograms

ABCD yields: `h_tight_iso_cluster_0`, `h_tight_noniso_cluster_0`, `h_nontight_iso_cluster_0`, `h_nontight_noniso_cluster_0`

Signal-matched: `h_tight_iso_cluster_signal_0`, etc.

Background: `h_tight_iso_cluster_notmatch_0`, etc.

## MergeSim.C

Merges per-sample outputs using TFileMerger (pure histogram addition -- weights were pre-applied).

**Signal merge** (3 -> 1): photon5 + photon10 + photon20 -> `{eff_outfile}_{var_type}.root`

**Jet merge** (6 -> 1): jet5 + jet8 + jet12 + jet20 + jet30 + jet40 -> `{eff_outfile}_jet_{var_type}.root`

**Response merge** (3 -> 1): photon response files

Note: MergeSim uses exactly 6 jet samples (jet5/8/12/20/30/40). It does NOT include jet10, jet15, or jet50.

## CalculatePhotonYield.C

### Signature

```cpp
void CalculatePhotonYield(const std::string &configname, bool isMC = false)
```

### Key Constants

```cpp
float mbdcorr = 25.2/42 / 0.57;              // MBD efficiency correction = 1.053
float solid_angle = 2 * M_PI * 0.7 * 2;       // 2pi * delta_eta
float nsimevents = 1E7;                        // generated events per sample
float simluminosity = nsimevents / photon20cross;  // sim luminosity
float jetevents = 0.3555 * 1E7;
float jetluminosity = jetevents / jet50cross;
```

### ABCD Purity Extraction

See [ABCD Method](../concepts/abcd-method.md) for full details.

1. Compute signal leakage fractions (cB, cC, cD) from MC
2. Run 20,000 Poisson toy MC experiments to solve the self-consistent ABCD equation
3. Fit purity vs ET with smooth function (see Purity Fitting below)

### Purity Fitting

- `fit_option=0` (default): Error function `[0]*Erf((x-[1])/[2])` over [8, 32] GeV
- `fit_option=1`: Pade rational function `([0]+[1]*x)/(1+[2]*x)`

68.3% CL confidence intervals via `TVirtualFitter::GetConfidenceIntervals()`.

### Bayesian Unfolding

Uses RooUnfoldBayes with response matrix `response_matrix_full_0`. Iterations 1-10 stored; default selection: `resultit` (usually 2). Leakage correction toggle: `resultleak` (0 or 1).

### Efficiency Correction

```cpp
float total_eff = eff_reco_val * eff_iso_val * eff_id_val * vertex_eff_val;
```

### Cross-Section Formula

```
d^3 sigma / (dpT deta dphi) = N_signal / (Lumi * dpT * deta * dphi * epsilon_total)
```

Luminosity from config `analysis.lumi` (pb^-1).

### Primary Output

`{final_outfile}_{var_type}.root` containing:
- `h_unfold_sub_result` -- THE final cross-section histogram
- `gpurity`, `gpurity_leak` -- purity TGraphErrors
- `h_R`, `h_R_notmatch` -- R-factor histograms
- `h_leak_B/C/D` -- leakage fractions
- All 10 unfolding iterations

## Orchestration

```bash
cd efficiencytool

# Two-pass for all samples
bash oneforall_tree.sh config_bdt_nom.yaml

# Then merge and compute yield
bash oneforall.sh config_bdt_nom.yaml
```

`oneforall_tree.sh` runs Pass 1 (vertex scan for all 10 samples in parallel), then Pass 2 (full analysis for all 10 samples in parallel).

`oneforall.sh` runs MergeSim.C then CalculatePhotonYield.C (isMC=true, then isMC=false).

## Gotchas

- The two-pass architecture is essential: Pass 1 must complete before Pass 2 starts
- mc_iso_scale/mc_iso_shift code defaults are 1.0/0.0 -- the values 1.2/0.1 are config overrides, not hardcoded
- MergeSim expects exactly 6 jet samples -- adding jet10/15/50 requires modifying MergeSim.C
- The response matrix is built only for tight + isolated + signal-matched clusters
- `jetevents = 0.3555 * 1E7` and `mbdcorr = 25.2/42/0.57` are hardcoded physics constants not documented elsewhere
