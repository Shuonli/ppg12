# Efficiency Tool -- Pipeline Stage Report

**Directory**: `/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/`
**Generated**: 2026-04-08

---

## 1. Purpose

The `efficiencytool/` directory implements the core physics analysis chain that transforms BDT-scored ROOT trees (from the upstream `FunWithxgboost/apply_BDT.C` stage) into a measured isolated photon cross section. It performs:

1. **Reconstruction efficiency calculation** -- truth-matching signal MC to reco clusters; computing reco, isolation, and identification efficiencies as TEfficiency objects.
2. **ABCD background subtraction** -- data-driven estimation of background contamination using a 2D sideband method (tight/non-tight vs isolated/non-isolated), with signal leakage corrections from simulation.
3. **Bayesian unfolding** -- correcting the purity-subtracted yield from reco-level ET bins to truth-level pT bins using a RooUnfold response matrix.
4. **Cross-section extraction** -- dividing the unfolded yield by luminosity, bin width, efficiency, and solid angle.
5. **Shower shape studies** -- comparing data and MC shower shape variables in ABCD regions.
6. **Double interaction / pileup studies** -- both full GEANT MC blending and toy vertex-shift simulation.
7. **Systematic variation generation** -- programmatic creation of 30+ analysis configs from a nominal base.

---

## 2. Key Files

### C++ Macros

| File | Lines | Description |
|------|-------|-------------|
| `RecoEffCalculator_TTreeReader.C` | 2719 | Main workhorse: reads slimtree via TTreeReader, applies all cuts, fills ABCD histograms for data/sim, computes TEfficiency objects for signal MC, builds RooUnfoldResponse matrices. Two-pass architecture (vertex scan then full analysis). |
| `CalculatePhotonYield.C` | 1170 | ABCD purity extraction, signal leakage correction, purity fitting (erf or Pade), Bayesian unfolding (RooUnfoldBayes, 1-10 iterations), efficiency correction, cross-section scaling. |
| `MergeSim.C` | 94 | TFileMerger wrapper: merges per-sample photon MC (photon5/10/20), jet MC (jet5/8/12/20/30/40), and response matrix files into combined outputs. |
| `ShowerShapeCheck.C` | 1917 | Shower shape distribution analysis: fills histograms of 11 shower-shape variables per ABCD region, per (eta_bin, pT_bin). Supports mix_weight and vtxscan_sim_override for double-interaction blending. |
| `DoubleInteractionCheck.C` | 1446 | Toy vertex-shift simulation: draws random second vertex from data, computes shifted cluster kinematics, re-evaluates BDT and NPB scores, fills ABCD yields under single/double/smeared scenarios. |
| `CrossSectionWeights.h` | 119 | Header-only: all MC cross-section constants (9 jet + 3 photon samples), `SampleConfig` struct with kinematic windows, `GetSampleConfig()` lookup function. |
| `MbdPileupHelper.h` | 207 | Pileup detection from MBD PMT timing: `calculateMbdPileupMetrics()`, `isMbdPileup()`, `checkMbdPileup()`. Computes avgsigma, prodsigma, maxsigma, deltas from 128 PMTs. |

### Other Notable C++ Macros

| File | Lines | Description |
|------|-------|-------------|
| `Closure.C` | ~480 | Closure test: runs CalculatePhotonYield on MC-as-data to verify pipeline roundtrip. |
| `FindBDTCut.C` | ~260 | Scans BDT threshold values to optimize signal/background separation. |
| `FindETCut.C` | ~210 | Scans isolation ET cut parameters (intercept, slope) for optimization. |
| `NPB_PurityStudy.C` | ~890 | Studies non-physics-background purity across 0mrad/1.5mrad crossing angle periods. |
| `PlotROC.C` | ~510 | ROC curve generation for BDT and shower-shape discriminants. |
| `VertexReweighting.C` | ~280 | Generates vertex-z reweighting histograms (data/MC ratio). |
| `SaturationStudy.C` | ~480 | Studies ADC saturation effects on cluster energy and shower shapes. |
| `Cluster_rbr.C` | ~1270 | Run-by-run QA: cluster rate, energy scale, shower shapes per run number. |
| `LumiCalculator.C` | ~190 | Luminosity calculation from trigger scalers and run lists. |
| `IsoROC_calculator.C` | ~570 | Isolation cut receiver-operating-characteristic curves. |
| `AnalyzeTruthPhotonTowers.C` | ~1060 | Tower-level analysis of truth-matched photon clusters. |
| `TruthJetEff.C` | ~1950 | Truth jet efficiency studies. |

### Shell Scripts

| File | Lines | Description |
|------|-------|-------------|
| `oneforall.sh` | 28 | Orchestrates `MergeSim.C` + `CalculatePhotonYield.C` (MC=true, MC=false). Called by condor. |
| `oneforall_tree.sh` | 39 | Two-pass orchestration: Pass 1 (vertex scan for all 10 samples in parallel), then Pass 2 (full analysis for all 10 samples in parallel). |
| `oneforall_tree_double.sh` | ~240 | Full pipeline including double-interaction blending for RecoEffCalculator. |
| `run_showershape_double.sh` | 61 | Two-pass double-interaction blending pipeline for ShowerShapeCheck. |
| `run_double_interaction.sh` | 41 | Runs DoubleInteractionCheck for all samples (photon5/10/20, jet5/8/12/20/30/40, data) then merges. |
| `run_showershape.sh` | ~70 | Runs ShowerShapeCheck for multiple samples and configs. |
| `oneforall.sub` | 28 | HTCondor submit file: queues one job per `config_bdt_*.yaml` file. |

### Python

| File | Lines | Description |
|------|-------|-------------|
| `make_bdt_variations.py` | 336 | Generates 30+ systematic variation configs from a nominal base. Defines VARIANTS, SYST_TYPES, SYST_GROUPS, OVERRIDE_MAP. Uses ruamel.yaml for comment-preserving YAML manipulation. |
| `make_showershape_variations.py` | ~210 | Generates showershape-specific config variants. |

### YAML Configs

37 `config_bdt_*.yaml` files in the main directory, plus backups in `config_backup/`. Additionally, `config_showershape_*.yaml` files for shower shape studies.

---

## 3. MergeSim.C

### Purpose

Merges per-sample efficiency and response matrix ROOT files into combined files using `TFileMerger`. This is pure histogram addition -- the cross-section weights were already applied during `RecoEffCalculator_TTreeReader.C` (each sample's histograms are pre-weighted by `cross_weight`).

### What It Merges

**Signal photon merge** (3 files into 1):
- Input: `{eff_outfile}_photon5_{var_type}.root`, `{eff_outfile}_photon10_{var_type}.root`, `{eff_outfile}_photon20_{var_type}.root`
- Output: `{eff_outfile}_{var_type}.root`

**Background jet merge** (6 files into 1):
- Input: `{eff_outfile}_jet5_{var_type}.root` through `{eff_outfile}_jet40_{var_type}.root`
- Output: `{eff_outfile}_jet_{var_type}.root`

**Response matrix merge** (3 photon files into 1):
- Input: `{response_outfile}_photon5_{var_type}.root`, `{response_outfile}_photon10_{var_type}.root`, `{response_outfile}_photon20_{var_type}.root`
- Output: `{response_outfile}_{var_type}.root`

### How Cross-Section Weights Are Applied

The weights are NOT applied in MergeSim.C. They are applied upstream in `RecoEffCalculator_TTreeReader.C` at histogram fill time via the `weight` variable (which equals `cross_weight = GetSampleConfig(filetype).weight * mix_weight`). The sample weights are computed as ratios relative to a reference sample:
- Photon samples: normalized to `photon20cross` (e.g., `photon5cross / photon20cross`)
- Jet samples: normalized to `jet50cross` (e.g., `jet10cross / jet50cross`)

Since histograms are filled with these relative weights, `hadd`-style merging (TFileMerger) correctly produces the combined weighted spectrum.

---

## 4. CalculatePhotonYield.C -- Full ABCD Implementation

### Function Signature

```cpp
void CalculatePhotonYield(const std::string &configname, bool isMC = false)
```

### Key Constants

```cpp
float mbdcorr = 25.2/42 / 0.57;              // MBD efficiency correction factor
float solid_angle = 2 * M_PI * 0.7 * 2;       // 2pi * delta_eta (eta range [-0.7, 0.7])
float nsimevents = 1E7;                        // generated events per sample
float simluminosity = nsimevents / photon20cross;  // sim luminosity (pb^-1)
float jetevents = 0.3555 * 1E7;               // effective jet events (trigger fraction)
float jetluminosity = jetevents / jet50cross;   // jet sim luminosity
```

### Input Files

- **Sim signal**: `{eff_outfile}_{var_type}.root` (merged photon MC)
- **Data/MC-as-data**: `{data_outfile}_{var_type}.root` (data) or `{eff_outfile}_jet_{var_type}.root` (MC closure)
- **Unfolding response**: `{response_outfile}_{var_type}.root`

### ABCD Region Definition

The ABCD method uses two orthogonal discriminators:
- **Tight vs Non-tight**: BDT score threshold (parametric: `bdt_min_intercept + bdt_min_slope * ET`)
- **Isolated vs Non-isolated**: Isolation ET within parametric cone cut

```
            Isolated (iso)       Non-isolated (noniso)
Tight       A = signal region    B = bkg template
Non-tight   C = bkg control      D = bkg normalization
```

Histogram names in the ROOT files:
- `h_tight_iso_cluster_0` (A), `h_tight_noniso_cluster_0` (B)
- `h_nontight_iso_cluster_0` (C), `h_nontight_noniso_cluster_0` (D)
- Signal-matched: `h_tight_iso_cluster_signal_0`, etc.
- Not-matched (background): `h_tight_iso_cluster_notmatch_0`, etc.

### Signal Leakage Corrections (cB, cC, cD)

Signal leakage fractions quantify how much true signal leaks into the non-signal ABCD regions:

```cpp
h_leak_B = h_tight_noniso_cluster_signal / h_tight_iso_cluster_signal;   // cB = B_sig / A_sig
h_leak_C = h_nontight_iso_cluster_signal / h_tight_iso_cluster_signal;   // cC = C_sig / A_sig
h_leak_D = h_nontight_noniso_cluster_signal / h_tight_iso_cluster_signal; // cD = D_sig / A_sig
```

These are computed bin-by-bin from the merged signal MC.

### ABCD Equation (with leakage)

The function `myfunc()` encodes the self-consistent ABCD equation:

```
NsigA = NA - R * (NB - cB * NsigA) * (NC - cC * NsigA) / (ND - cD * NsigA)
```

This is solved numerically by finding the root of:
```
f(NsigA) = NsigA - (NA - R * (NB - cB*NsigA) * (NC - cC*NsigA) / (ND - cD*NsigA))
```

where `R = 1` (the double ratio is already embedded in the formula).

### R Factor Calculation

Two R-factor calculations are performed:

1. **Data-driven R**: Background yields are estimated as `data - sim_signal_leak` for each ABCD region:
   ```
   R = (A_bkg / B_bkg) * (D_bkg / C_bkg)
   ```

2. **MC truth R**: Uses `_notmatch` histograms (non-truth-matched clusters directly from simulation):
   ```
   R_notmatch = (A_notmatch / B_notmatch) * (D_notmatch / C_notmatch)
   ```

### Purity Extraction via Toy MC

For each pT bin, 20,000 Poisson toy samples are drawn:
1. Compute effective number of events: `A_eff = A^2 / sum(w^2)`
2. Draw Poisson: `NA_eff = Poisson(A_eff)`, rescale: `NA = A * NA_eff / A_eff`
3. Solve ABCD equation (without leakage) -> `h_result`; compute purity = `NsigA / NA`
4. Repeat with Gaussian-smeared leakage corrections (cB, cC, cD) -> `h_result_leak`
5. Fit each purity distribution with a Gaussian -> central value and width

### Purity Fitting

The per-bin purity values are fit with a smooth function:
- **Default (fit_option=0)**: Error function: `[0]*Erf((x - [1])/[2])` over [8, 32] GeV
- **Pade (fit_option=1)**: Rational function: `([0] + [1]*x) / (1 + [2]*x)`

Confidence intervals are computed at 68.3% CL via `TVirtualFitter::GetConfidenceIntervals()`.

The `fittingerror` config parameter shifts the purity by +/- 1 confidence interval for systematic studies.

### MC Purity Correction

When `mc_purity_correction = 1`, the code:
1. Runs MC closure (isMC=true) to get truth purity `g_purity_truth`
2. Fits the MC truth purity with the same functional form
3. Computes ratio: `mc_corr = truth_purity / fitted_purity` per bin
4. Applies this ratio to the data purity: `purity_data_corrected = fitted_purity * mc_corr`

### Unfolding

Uses RooUnfoldBayes (iterative Bayesian unfolding):
- Response matrix: `response_matrix_full_0` from merged photon response files
- Number of iterations: 1 to 10, stored as separate histograms
- Selected iteration: `configYaml["analysis"]["unfold"]["resultit"]` (default: 2)
- Leakage correction toggle: `resultleak` (0 = no leakage, 1 = with leakage)
- Response reweighting: optional (controlled by `reweight` config, currently forced to 0 in code)

### Efficiency Correction

After unfolding, the spectrum is divided by the total efficiency per bin:

```cpp
float total_eff = eff_reco_val * eff_iso_val * eff_id_val * vertex_eff_val;
```

where:
- `eff_reco_val` = reconstruction efficiency (TEfficiency `eff_reco_eta_0`)
- `eff_iso_val` = isolation efficiency (TEfficiency `eff_iso_eta_0`)
- `eff_id_val` = identification/tight-cut efficiency (TEfficiency `eff_id_eta_0`)
- `vertex_eff_val` = vertex + MBD efficiency: `h_truth_pT_vertexcut_mbd_cut / h_truth_pT_vertexcut + mbd_eff_scale`

### Cross-Section Formula

The final differential cross-section is:

```
d^3 sigma / (dpT deta dphi) = N_signal / (Lumi * dpT * deta * dphi * epsilon_total)
```

Implemented via `scale_histogram(h, luminosity)`:
```cpp
scale = 1.0 / binwidth / luminosity;
h->SetBinContent(ibin, h->GetBinContent(ibin) * scale);
```

The luminosity is read from config (`analysis.lumi`, in pb^-1). For MC closure, `jetluminosity` is used instead.

### Output File

`{final_outfile}_{var_type}[_mc].root` containing:
- `h_unfold_sub_result` -- the final cross-section histogram (selected iteration, with or without leakage)
- `h_unfold_sub_result_woeff` -- same without efficiency correction
- `gpurity`, `gpurity_leak` -- purity TGraphErrors (with/without leakage)
- `g_purity_truth` -- MC truth purity (MC only)
- `f_purity_fit`, `f_purity_leak_fit` -- purity fit functions
- `confInt`, `confInt_leak` -- confidence interval graphs
- `h_R`, `h_R_notmatch` -- R-factor histograms
- `h_leak_B`, `h_leak_C`, `h_leak_D` -- leakage fraction histograms
- All 10 unfolding iterations: `h_unfold_sub_{1..10}`, `h_unfold_sub_leak_{1..10}`
- Purity distributions per bin: `h_purity_{1..N}`, `h_purity_leak_{1..N}`

---

## 5. RecoEffCalculator_TTreeReader.C

### Function Signature

```cpp
void RecoEffCalculator_TTreeReader(
    const std::string &configname,
    const std::string filetype = "jet40",
    bool do_vertex_scan = false,
    float mix_weight = 1.0,
    const std::string vtxscan_sim_override = "")
```

### Two-Pass Architecture

**Pass 1 (vertex scan)**: `do_vertex_scan = true`
- Only fills `h_vertexz` histogram (vertex-z distribution)
- Writes to `{outfile}_vtxscan.root`
- Used to build on-the-fly vertex reweighting ratio (data/MC)

**Pass 2 (full analysis)**: `do_vertex_scan = false`
- Reads vtxscan files from Pass 1 (or `vtxscan_sim_override` for blended MC)
- Computes data/MC vertex-z ratio for reweighting
- Performs full event loop with all cuts and histogram filling

### Input File Resolution

- Simulation: `{photon_jet_file_root_dir}{filetype}{photon_jet_file_branch_dir}`
  - Example: `/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root`
- Data: `{data_file}` (supports wildcards)
- `_nom` aliases: `photon10_nom` reads from `photon10` path; `jet12_nom` reads from `jet12` path

### Cross-Section Weighting

Uses `PPG12::GetSampleConfig(filetype)` from `CrossSectionWeights.h`:
```cpp
float cross_weight = GetSampleConfig(filetype).weight * mix_weight;
```
Applied to all histogram fills as `weight = cross_weight * vertex_weight`.

### Vertex Reweighting

For simulation, vertex-z distributions are reweighted to match data:
1. Primary method: on-the-fly ratio from Pass 1 vtxscan files (`h_vertexz` data / sim, both normalized to unit integral)
2. Fallback: pre-made `vertex_reweight_file` containing `h_vertexz_ratio_data_over_mccombined`

The weight per event: `vertex_weight = h_vertex_reweight->GetBinContent(h_vertex_reweight->FindBin(*vertexz))`

### Energy Scale/Resolution Smearing

For simulation only:
```cpp
cluster_Et[icluster] = cluster_Et[icluster] * clusterescale;
if (clustereres > 0)
    cluster_Et[icluster] = cluster_Et[icluster] * rand->Gaus(1, clustereres);
```
Controlled by `analysis.cluster_escale` (default 1.0) and `analysis.cluster_eres` (default 0.0).

### MC Isolation Fudge

Simulation isolation ET is adjusted to match data:
```cpp
recoisoET = recoisoET * mc_iso_scale;
recoisoET += mc_iso_shift;
```
Default: `mc_iso_scale = 1.2`, `mc_iso_shift = 0.1`.

### Efficiencies Computed

Three factorized efficiencies, computed as TEfficiency objects per eta bin:

1. **`eff_reco_eta_{ieta}`**: Reconstruction efficiency
   - Denominator: all truth photons (pid=22, photonclass<3, truth_iso < truthisocut, in eta bin, truth vertex within cut)
   - Numerator: truth photons with a reco cluster match within `dR < eff_dR` (default 0.1)
   
2. **`eff_iso_eta_{ieta}`**: Isolation efficiency
   - Denominator: reconstructed truth photons (passed reco match)
   - Numerator: reconstructed truth photons passing reco isolation cut
   
3. **`eff_id_eta_{ieta}`**: Identification (tight-cut) efficiency
   - Denominator: reconstructed + isolated truth photons
   - Numerator: above that also pass the tight shower-shape/BDT cuts

Additional efficiencies: `eff_converts_eta` (conversion probability), `eff_all_eta` (total combined efficiency).

### Truth Matching Logic

1. Build `particle_trkidmap`: maps `particle_trkid[i]` -> particle index `i`
2. Build `photon_reco` map: for each truth particle that is a photon (pid=22) with photonclass<3 (direct or fragmentation), truth isolation < cut, and within eta/pT range
3. For each reco cluster, look up `cluster_truthtrkID[icluster]` in `particle_trkidmap`
4. If found and the truth particle is in `photon_reco`, and `dR < eff_dR`: mark as matched
5. Unmatched clusters fill `_notmatch` histograms (background for R-factor)

### Response Matrix

Built only for tight + isolated + signal-matched clusters:
```cpp
responses_full[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
```
Optional response reweighting via rational function `f_reweight`.

Half-sample response matrices (`responses_half`) are also built for closure tests (first half of events only).

### Trigger Selection

Data events must fire at least one trigger from `trigger_used` list:
```cpp
// trigger_used can be a scalar int or a YAML sequence
// Each trigger is checked: scaledtrigger[trigger_id] must be true
```
Default: trigger 30 (photon trigger). Multi-trigger support: `[26, 29, 30, 31, 36, 37, 38]`.

### Run Range Filtering

Data events are filtered by run number:
```cpp
if (*runnumber < run_min || *runnumber > run_max) continue;
```
Also supports `run_list_file` for explicit run lists. Default: `run_min: 51274`, `run_max: 54000` (1.5 mrad period).

### MBD Pileup Cut

Event-level pileup rejection via MBD timing RMS:
```cpp
mbd_avg_sigma = 0.5 * (mbd_std_north + mbd_std_south);
if (mbd_avg_sigma < mbd_avg_sigma_min || mbd_avg_sigma > mbd_avg_sigma_max) continue;
```
Default: no cut (min = -inf, max = +inf).

### Output Files

Per-sample:
- `{eff_outfile}_{filetype}_{var_type}.root` -- efficiency histograms, ABCD yields
- `{response_outfile}_{filetype}_{var_type}.root` -- response matrix

Data:
- `{data_outfile}_{var_type}.root`

---

## 6. ShowerShapeCheck.C

### Function Signature

```cpp
void ShowerShapeCheck(
    const std::string &configname,
    const std::string filetype,
    bool doinclusive = true,
    bool do_vertex_scan = false,
    float mix_weight = 1.0,
    const std::string vtxscan_sim_override = "")
```

### Purpose

Fills normalized shower-shape distributions for 11 variables, broken down by:
- ABCD region (tight_iso, tight_noniso, nontight_iso, nontight_noniso)
- eta bin
- pT bin
- Signal vs inclusive (truth-matched vs all clusters)

The 11 shower-shape variables: `weta_cogx`, `wphi_cogx`, `wr_cogx`, `et1`, `et2`, `et3`, `et4`, `e11/e33`, `e32/e35`, `bdt_score`, `npb_score`.

### mix_weight Mechanism

The `mix_weight` parameter is a multiplicative fraction applied to `cross_weight`:
```cpp
cross_weight *= mix_weight;
```
All histogram fills use this combined weight. For double-interaction blending:
- Nominal MC samples: `mix_weight = 0.813` (= 1 - 0.187)
- Double-interaction MC samples: `mix_weight = 0.187`

### Vertex Scan Passes

Same two-pass architecture as `RecoEffCalculator_TTreeReader.C`:
- **Pass 1** (`do_vertex_scan=true`): only fills `h_vertexz` with mix_weight-scaled entries
- **Pass 2** (`do_vertex_scan=false`): full analysis using blended vtxscan for vertex reweighting

### Output

`{eff_outfile}shower_shape_{filetype}[_inclusive]_{suffix}.root` containing:
- Shower-shape distributions: `h_{var}_{region}_{etabin}_{ptbin}`
- Vertex distributions: `h_vertexz`
- ABCD yield histograms (same naming as RecoEffCalculator)

---

## 7. DoubleInteractionCheck.C

### Function Signature

```cpp
void DoubleInteractionCheck(
    const std::string &configname,
    const std::string filetype,
    bool doinclusive = true)
```

### Toy Vertex-Shift Algorithm (Simulation Only)

For each MC event:

1. Draw a random second vertex from the unrestricted data vertex distribution (`vertex_scan_data_file`).
2. Compute double vertex: `double_vtx = (original_vtx + random_vtx) / 2.0`
3. Recalculate cluster kinematics using CEMC cylinder geometry:
   ```cpp
   constexpr float cemc_radius_cm = 93.5f;
   new_z_over_r = sinh(old_eta) + (old_vtxz - new_vtxz) / 93.5;
   new_eta = asinh(new_z_over_r);
   cluster_energy = old_et * cosh(old_eta);
   new_et = cluster_energy / cosh(new_eta);
   ```
4. Re-evaluate NPB score from scratch at the shifted vertex/ET/eta (using `buildNpbFeatureVector()` and TMVA::RBDT)
5. Re-evaluate shower-shape BDT score from scratch (using `buildSsBdtFeatureVector()` with model-specific feature vectors for base_v3E, base_v2E, base_v1E, etc.)
6. Re-classify tight/non-tight and iso/non-iso under shifted kinematics

### Additional Vertex Smearing

On top of the double vertex, Gaussian smearing is applied at 4 levels:
```cpp
std::vector<float> smear_sigmas = {5.0, 10.0, 15.0, 20.0};  // in cm
smeared_vtx = double_vtx + Gaus(0, sigma);
```

### Crossing Angle Classification (Data Only)

```cpp
bool is_0mrad(int run) { return (run >= 47000 && run <= 51274); }
bool is_1p5mrad(int run) { return (run >= 51274 && run <= 54000); }
```

MBD avg sigma is computed per event and histogrammed vs run number for monitoring pileup rate evolution.

### 4 Analysis Tasks

1. ABCD yields as 2D histograms (ET vs NPB score) for pileup contamination studies
2. MBD avg sigma_t vs run number (data only) for pileup rate tracking
3. ABCD yields as 2D (ET vs MBD avg sigma) for data pileup metric
4. Single vs double interaction ABCD comparison (sim only) for truth purity

### Output Histograms

Naming pattern: `h_{tight/nontight}_{iso/noniso}_cluster_{single/double}[_smearN][_signal/_incsig]_{etabin}`

6 interaction types: `single`, `double`, `double_smear5`, `double_smear10`, `double_smear15`, `double_smear20`.

Shower shape distributions: per (interaction_type, variable, eta_bin, pT_bin) for all 11 variables.

---

## 8. CrossSectionWeights.h -- All Constants

### Namespace: `PPG12`

**Photon-jet samples** (cross sections in pb):
| Variable | Value | Truth pT Window |
|----------|-------|------------------|
| `photon5cross` | 146,359.3 | 0--14 GeV |
| `photon10cross` | 6,944.675 | 14--30 GeV |
| `photon20cross` | 130.4461 | 30+ GeV |

**QCD jet samples** (cross sections in pb):
| Variable | Value | Truth Jet pT Window | Cluster ET Upper |
|----------|-------|---------------------|-----------------|
| `jet5cross` | 1.3878e+08 | 7--9 GeV | 10 GeV |
| `jet8cross` | 1.15e+07 | 9--14 GeV | 15 GeV |
| `jet10cross` | 3.997e+06 | 10--15 GeV | 18 GeV |
| `jet12cross` | 1.4903e+06 | 14--21 GeV | 23 GeV |
| `jet15cross` | 4.073e+05 | 15--21 GeV | 23 GeV |
| `jet20cross` | 6.2623e+04 | 21--32 GeV | 35 GeV |
| `jet30cross` | 2.5298e+03 | 32--42 GeV | 45 GeV |
| `jet40cross` | 1.3553e+02 | 42--100 GeV | 100 GeV |
| `jet50cross` | 7.3113 | 52--100 GeV | 100 GeV |

**Default generated events**: `default_nsimevents = 1E7`

**Weight normalization**:
- Photon samples: `weight = photon{N}cross / photon20cross`
- Jet samples: `weight = jet{N}cross / jet50cross`

### SampleConfig Struct

```cpp
struct SampleConfig {
    float photon_pt_lower, photon_pt_upper;
    float jet_pt_lower, jet_pt_upper;
    float cluster_ET_upper;
    float weight;
    bool isbackground;
    bool valid;
};
```

`GetSampleConfig(filetype)` returns the correct config for any sample name string including `_double` and `_nom` variants.

---

## 9. MbdPileupHelper.h -- API

### `calculateMbdPileupMetrics(northt, southt, northq, southq, hitcut, time_cut, charge_cut)`

Template function accepting any array type. Returns `MbdPileupResult`:
- Iterates over 64 north + 64 south PMT channels
- Selection: `|time| < 25 ns`, `charge > 0.4`, minimum `hitcut` (default 2) per side
- Computes per-side: mean time, RMS time, max-min time spread
- Metrics: `avgsigma` (primary), `prodsigma`, `maxsigma`, `avgdelta`, `proddelta`, `maxdelta`

### `isMbdPileup(result, strength, comfort_cut, strict_cut, draconian_cut)`

Three modes:
- `COMFORT`: `prodsigma >= 1.0`
- `STRICT`: `avgsigma >= 2.0`
- `DRACONIAN`: `avgsigma >= 1.5`

### `checkMbdPileup(...)` -- convenience one-call wrapper

---

## 10. Config Schema -- `config_bdt_nom.yaml`

### `input` Section

| Field | Value | Description |
|-------|-------|-------------|
| `photon_jet_file_root_dir` | `/sphenix/user/.../FunWithxgboost/` | Root directory for sim trees |
| `photon_jet_file_branch_dir` | `/bdt_split.root` | Suffix appended to sample name |
| `data_file` | `/.../part_*_with_bdt_split.root` | Data file glob pattern |
| `tree` | `slimtree` | TTree name |
| `cluster_node_name` | `CLUSTERINFO_CEMC` | Cluster branch prefix |
| `bdt_model_name` | `base_E` | Fallback BDT model name |
| `bdt_et_bin_edges` | `[8, 15, 35]` | ET bin boundaries for binned BDT (N+1 values) |
| `bdt_et_bin_models` | `[base_v1E, base_v1E]` | BDT model per ET bin (N values) |

### `output` Section

| Field | Value | Description |
|-------|-------|-------------|
| `eff_outfile` | `.../results/MC_efficiency` | Base path for efficiency output |
| `response_outfile` | `.../results/MC_response` | Base path for response matrix |
| `data_outfile` | `.../results/data_histo` | Base path for data histograms |
| `var_type` | `bdt_nom` | Unique suffix for this config |
| `final_outfile` | `.../results/Photon_final` | Base path for final cross-section |

### `analysis` Section -- Top Level

| Field | Value | Description |
|-------|-------|-------------|
| `eta_bins` | `[-0.7, 0.7]` | Pseudorapidity bins |
| `vertex_cut` | `60.0` | Maximum |vertex_z| (cm) |
| `truth_iso_max` | `4.0` | Truth isolation ET cut (GeV) |
| `reco_iso_min` | `-20.0` | Reco isolation ET lower bound |
| `reco_iso_max_b` | `0.502095` | Isolation cut intercept |
| `reco_iso_max_s` | `0.0433036` | Isolation cut slope (per GeV) |
| `reco_noniso_min_shift` | `0.8` | Gap between iso max and noniso min (GeV) |
| `reco_noniso_max` | `20.0` | Non-iso upper bound |
| `reco_min_ET` | `5.0` | Minimum reco cluster ET (GeV) |
| `cone_size` | `3` | Isolation cone R=0.3 |
| `eff_dR` | `0.1` | Truth-reco matching dR |
| `trigger_used` | `30` | Trigger bit index |
| `pT_bins` | `[8,10,12,14,16,18,20,22,24,26,28,32,36]` | Reco pT bin edges (12 bins) |
| `pT_bins_truth` | `[7,8,10,12,...,36,45]` | Truth pT bin edges (extended for overflow) |
| `use_topo_iso` | `2` | 0=standard cone, 1=topo R=0.3, 2=topo R=0.4 |
| `iso_threshold` | `1` | Use threshold-based isolation (120 MeV) |
| `mc_iso_shift` | `0.1` | MC isolation shift (GeV) |
| `mc_iso_scale` | `1.2` | MC isolation scale factor |
| `vertex_reweight_on` | `1` | Enable vertex reweighting |
| `fit_option` | `1` | Purity fit: 0=erf, 1=Pade |
| `cluster_escale` | `1.0` | Energy scale factor |
| `cluster_eres` | `0.0` | Energy resolution smearing |
| `mc_purity_correction` | `0` | Apply MC purity correction |
| `n_nt_fail` | `0` | Minimum number of tight cuts failed for non-tight |
| `bdt_on` | `1` | Include BDT in tight/non-tight classification |
| `run_min` | `51274` | Data run range lower bound |
| `run_max` | `54000` | Data run range upper bound |
| `lumi` | `16.2735` | Integrated luminosity (pb^-1) |

### `analysis.tight` Section

| Field | Value | Description |
|-------|-------|-------------|
| `bdt_min` | `0.7` | BDT lower bound (flat, overridden by parametric) |
| `bdt_max` | `1.0` | BDT upper bound |
| `bdt_min_slope` | `-0.015` | BDT threshold slope |
| `bdt_min_intercept` | `0.80` | BDT threshold intercept |
| `weta_cogx_max` | `1.0` | Max weta_cogx (disabled; set to 1.0) |
| `et1_min` | `0.5` | Min et1 fraction |
| `e32_over_e35_min` | `0.8` | Min e32/e35 ratio |
| (many others) | `[0, 1]` | Effectively disabled (full range) |

Parametric BDT cut: `tight_bdt_min = 0.80 - 0.015 * ET`

### `analysis.non_tight` Section

| Field | Value | Description |
|-------|-------|-------------|
| `bdt_max` | `0.5` | BDT upper bound for non-tight |
| `bdt_min` | `0.02` | BDT lower bound for non-tight |
| `bdt_min_slope` | `-0.015` | Non-tight BDT max slope |
| `bdt_min_intercept` | `0.80` | Non-tight BDT max intercept |

Non-tight parametric BDT upper cut: `non_tight_bdt_max = 0.80 - 0.015 * ET` (matches tight threshold)

### `analysis.common` Section

| Field | Value | Description |
|-------|-------|-------------|
| `wr_cogx_bound` | `0.0` | Min wphi/weta ratio |
| `cluster_weta_cogx_bound` | `2.0` | Max weta_cogx |
| `npb_cut_on` | `1` | Enable NPB score cut |
| `npb_score_cut` | `0.5` | NPB score threshold |
| `e11_over_e33_max` | `0.98` | Max e11/e33 (removes bad clusters) |

### `analysis.unfold` Section

| Field | Value | Description |
|-------|-------|-------------|
| `resultit` | `2` | Selected unfolding iteration |
| `resultleak` | `1` | 1 = use leakage-corrected purity |
| `reweight` | `1` | Response matrix reweighting |

---

## 11. Shell Scripts -- Orchestration

### `oneforall.sh`

Called per-config by HTCondor (or manually). Steps:
1. `MergeSim.C(config)` -- merge per-sample files (assumes RecoEffCalculator already run)
2. `CalculatePhotonYield.C(config, true)` -- MC closure (isMC=true)
3. `CalculatePhotonYield.C(config, false)` -- Data analysis (isMC=false)

### `oneforall_tree.sh`

Full two-pass pipeline. Steps:
1. **Pass 1**: 10 parallel `RecoEffCalculator_TTreeReader.C` calls with `do_vertex_scan=true` (photon5/10/20, jet5/8/12/20/30/40, data)
2. Wait for all to complete
3. **Pass 2**: 10 parallel `RecoEffCalculator_TTreeReader.C` calls with full analysis

### `run_showershape_double.sh`

Double-interaction blending pipeline for shower shapes:

```
Pass 1: vertex scans (5 parallel jobs)
  photon10_double (mix_weight=0.187)
  jet12_double    (mix_weight=0.187)
  photon10_nom    (mix_weight=0.813)
  jet12_nom       (mix_weight=0.813)
  data            (mix_weight=1.0)

hadd: nom + double vtxscans -> combined vtxscan per sample type

Pass 2: full analysis (5 parallel jobs)
  All MC samples use SAME combined vtxscan for vertex reweighting
  Each sample applies its own mix_weight to all fills

hadd: nom + double full outputs -> combined outputs
```

Default fraction: 18.7% double (0 mrad). Accepts second argument for 7.2% (1.5 mrad).

### `run_double_interaction.sh`

Runs `DoubleInteractionCheck.C` for all sample types (10 parallel jobs), then merges:
- Signal: photon5+10+20 -> `_signal.root`
- Jet: jet5+8+12+20+30+40 -> `_jet.root`
- Combined: signal+jet -> `_inclusive.root`

### `oneforall.sub`

HTCondor submit file:
```
executable = /bin/bash
arguments = "oneforall.sh $(cfg)"
queue cfg matching (config_bdt*.yaml)
```
Automatically queues one job per matching YAML config file.

---

## 12. Systematics -- `make_bdt_variations.py`

### Architecture

1. Define VARIANTS list of dicts with override parameters
2. Define SYST_TYPES mapping each systematic source to a mode and group
3. Define SYST_GROUPS for quadrature combination
4. Deep-copy nominal config, apply overrides, set `var_type`, write to file

### OVERRIDE_MAP

Maps flat parameter names to YAML paths:
```python
"tight_bdt_min_intercept" -> analysis.tight.bdt_min_intercept
"nt_bdt_max" -> analysis.non_tight.bdt_max
"npb_score_cut" -> analysis.common.npb_score_cut
"clusterescale" -> analysis.cluster_escale
"lumi" -> analysis.lumi
```
(22 override keys total)

### Variant List (30+ configs)

**Cross-checks (no syst_type)**:
- `nom` -- identical to base
- `0rad` -- run range 47289-51274, lumi=32.6574 pb^-1
- `all` -- run range 47289-54000, lumi=49.562 pb^-1
- `mbdrms01/03/05` -- MBD avg sigma cuts
- `tightbdtflat06` -- flat BDT cut at 0.6
- `mciso_noscaleshift`, `mciso_no_shift`, `mciso_no_scale` -- MC isolation variants
- `common_wr03/05` -- common wr_cogx cut variations
- `split` -- CLUSTERINFO_CEMC split node
- `etbin_E_E`, `etbin_v1E_v3E` -- alternate ET-binned BDT models
- `b2bjet` -- back-to-back jet cut
- `timingcut_2/5` -- cluster-MBD timing window

**Tight BDT cut** (syst_type=`tight_bdt`, group=`eff`):
- `tightbdt50` (intercept=0.70, down), `tightbdt70` (intercept=0.90, up)

**Non-tight BDT boundary** (syst_type=`nt_bdt`, group=`purity`):
- `ntbdtmin05` (min=0.05), `ntbdtmin10` (min=0.10), both one-sided

**Non-iso sideband** (syst_type=`noniso`, group=`purity`):
- `noniso04` (shift=0.1, down), `noniso10` (shift=1.0, up)

**NPB cut** (syst_type=`npb_cut`, group=`eff`):
- `npb03` (cut=0.3, down), `npb07` (cut=0.7, up)

**Purity fit** (syst_type=`purity_fit`, group=`purity`):
- `purity_pade` (fit_option=0, one-sided)

**MC purity correction** (syst_type=`mc_purity_correction`, group=`purity`):
- `mc_purity_correction` (mc_purity_correction=1, one-sided)

**Vertex reweighting** (syst_type=`vtx_reweight`, group=`eff`):
- `vtxreweight0` (vertex_reweight_on=0, one-sided)

**Unfolding reweighting** (syst_type=`reweight`, group=`unfold`):
- `no_unfolding_reweighting` (reweight=0, one-sided)

**BDT model** (syst_type=`bdt_model`, group=`eff`):
- `etbin_v3E_v3E` (both bins use base_v3E, max-envelope)

**Energy scale** (syst_type=`escale`, group=`escale`):
- `energyscale26up` (scale=1.026, down), `energyscale26down` (scale=0.974, up)

**Energy resolution** (syst_type=`eres`, group=`eres`):
- `energyresolution5` (eres=0.05, max-envelope)
- `energyresolution7/8` (cross-checks)

### Systematic Modes

- `two_sided`: (up, down) pair -- take half-difference as uncertainty
- `one_sided`: single variation -- take full difference from nominal
- `max`: take maximum deviation among variants (envelope method)

### Quadrature Groups

```python
SYST_GROUPS = {
    "purity": ["noniso", "nt_bdt", "purity_fit", "mc_purity_correction"],
    "eff":    ["tight_bdt", "npb_cut"],
    "escale": ["escale"],
    "eres":   ["eres"],
    "mbd":    ["mbd"],
    "unfolding": ["reweight"],
}
FINAL_SYSTS = ["purity", "eff", "escale", "eres", "mbd", "unfolding"]
```

### Luminosity Uncertainty

```python
LUMI_SYST = {"down": (25.2 - 23.5) / 25.2, "up": (27.5 - 25.2) / 25.2}
# down: 6.7%, up: 9.1%
```

---

## 13. Cross-Section Constants -- Where They Appear

The cross-section weights are centralized in `CrossSectionWeights.h` (the "single source of truth"):

| Constant | Value (pb) |
|----------|-----------|
| `photon5cross` | 146,359.3 |
| `photon10cross` | 6,944.675 |
| `photon20cross` | 130.4461 |
| `jet5cross` | 1.3878e+08 |
| `jet8cross` | 1.15e+07 |
| `jet10cross` | 3.997e+06 |
| `jet12cross` | 1.4903e+06 |
| `jet15cross` | 4.073e+05 |
| `jet20cross` | 6.2623e+04 |
| `jet30cross` | 2.5298e+03 |
| `jet40cross` | 1.3553e+02 |
| `jet50cross` | 7.3113 |

Files that `#include "CrossSectionWeights.h"`:
- `RecoEffCalculator_TTreeReader.C` (line 23)
- `CalculatePhotonYield.C` (line 5)
- `ShowerShapeCheck.C` (line 2)
- `DoubleInteractionCheck.C` (line 2)

All use `using namespace PPG12;` to access the constants directly.

---

## 14. Input/Output

### Input Files (from BDT Stage)

**Simulation** (per sample, BDT-scored):
- Pattern: `{photon_jet_file_root_dir}{sample}{photon_jet_file_branch_dir}`
- Example: `/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root`
- Samples: photon5, photon10, photon20, jet5, jet8, jet10, jet12, jet15, jet20, jet30, jet40, jet50
- Double-interaction: photon10_double, jet12_double

**Data** (BDT-scored):
- `/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/condorout/part_*_with_bdt_split.root`

### Output Files (in `results/`)

**Per-sample efficiency files**:
- `MC_efficiency_{sample}_{var_type}.root` (one per photon/jet sample)
- `MC_efficiency_{sample}_{var_type}_vtxscan.root` (vertex scan, Pass 1)

**Merged files** (from MergeSim):
- `MC_efficiency_{var_type}.root` (merged photon5+10+20)
- `MC_efficiency_jet_{var_type}.root` (merged jets)
- `MC_response_{var_type}.root` (merged response matrices)

**Data histograms**:
- `data_histo_{var_type}.root`
- `data_histo_{var_type}_vtxscan.root`

**Final cross-section**:
- `Photon_final_{var_type}.root` (data)
- `Photon_final_{var_type}_mc.root` (MC closure)

**Shower shape outputs**:
- `MC_efficiencyshower_shape_{sample}[_inclusive]_{suffix}.root`
- `MC_efficiencyshower_shape_{sample}_combined[_inclusive]_{suffix}.root` (blended)

**Double interaction outputs**:
- `MC_efficiency_double_interaction_check_{sample}[_inclusive].root`
- `data_histo_double_interaction_check.root`

---

## 15. ABCD Details -- Parametric Thresholds

### Isolation (iso / noniso)

```
recoiso_max = reco_iso_max_b + reco_iso_max_s * clusterET
            = 0.502095 + 0.0433036 * ET  (GeV)

recononiso_min = recoiso_max + reco_noniso_min_shift
               = recoiso_max + 0.8  (GeV gap)

Isolated:     reco_iso_min (-20) < recoisoET < recoiso_max
Non-isolated: recononiso_min < recoisoET < reco_noniso_max (20)
```

For MC: `recoisoET = recoisoET * mc_iso_scale + mc_iso_shift` (default: `* 1.2 + 0.1`).

### Tight vs Non-tight

**Common pre-selection** (applied to both tight and non-tight):
- `cluster_prob` in [0, 1]
- `e11/e33` in [0, 0.98]
- `wr_cogx > 0.0` (wphi/weta ratio)
- `cluster_weta_cogx < 2.0`
- NPB score > 0.5 (when `npb_cut_on = 1`)
- Optional back-to-back jet cut

**Tight** (all of these must pass):
- `weta_cogx` in [min, max] (parametric: `max = max_b + max_s * ET`)
- `wphi_cogx` in [min, max] (parametric)
- `et1` in [min, max] (parametric: `min = min_b + min_s * ET`)
- `et2`, `et3` in ranges
- `e11/e33` in range
- `e32/e35` in range
- `et4` in range
- `prob` in range
- BDT score in [`bdt_min_intercept + bdt_min_slope * ET`, `bdt_max`]
  - Nominal: `[0.80 - 0.015*ET, 1.0]`

**Non-tight** (cluster must be within non-tight bounds AND fail enough tight cuts):
1. All non-tight ranges must be satisfied (wider than tight)
2. BDT score must be in [`non_tight_bdt_min`, `non_tight_bdt_max_intercept + non_tight_bdt_max_slope * ET`]
   - Nominal: `[0.02, 0.80 - 0.015*ET]`
3. Count number of failed tight sub-cuts (`nfail`), using `{var}_on` flags to select which cuts count
4. Must fail more than `n_nt_fail` (default 0, so must fail at least 1)
5. Optional: specific `{var}_fail` flags can require certain cuts to be failed

In the nominal config (BDT-based analysis), `bdt_on=1` and most other `*_on` flags are 0, so the BDT score is the primary tight/non-tight discriminant.

---

## 16. Gotchas and Synchronization Requirements

### Constants That Must Stay in Sync

1. **Cross-section weights**: All files include `CrossSectionWeights.h` -- this is now the single source of truth. Any change to cross sections must be made there only.

2. **pT bins**: `analysis.pT_bins` in the YAML config must match `plotcommon.h::ptRanges` in the plotting directory. Reco bins: `[8,10,12,14,16,18,20,22,24,26,28,32,36]` (12 bins).

3. **BDT feature vectors**: The feature order in `DoubleInteractionCheck.C::buildSsBdtFeatureVector()` must exactly match the training feature order in `apply_BDT.C` for each model variant. Model-specific branches exist (base_v3E, base_v2E, base_v1E, etc.).

4. **`var_type` uniqueness**: Each config must have a unique `output.var_type`. Duplicate var_types silently overwrite output files.

5. **Run range and luminosity**: `analysis.run_min/run_max` must correspond to `analysis.lumi`. Key values:
   - 0 mrad: 47289-51274, lumi=32.6574
   - 1.5 mrad: 51274-54000, lumi=16.2735 (or 16.8588 in CLAUDE.md)
   - All: 47289-54000, lumi=49.562

6. **Trigger index**: `trigger_used: 30` is the photon trigger. Using wrong trigger will silently select wrong events.

### Common Failure Modes

1. **Missing vtxscan files**: Pass 2 fails if Pass 1 vtxscan files don't exist. The code falls back to `vertex_reweight_file` but warns.

2. **Zero-content histograms**: ABCD equation root-finding (`TF1::GetX`) can return NaN/Inf for bins with zero counts. Guarded by `root == root` check.

3. **Purity > 1 or < 0**: The self-consistent ABCD equation can produce unphysical purities in low-statistics bins. These are passed through without clamping.

4. **mc_purity_correction depends on MC run first**: Data run with `mc_purity_correction=1` requires the MC output file (`_mc.root`) to already exist with `g_mc_purity_fit_ratio`.

5. **jet sample overlap**: jet10 (10-15 GeV) and jet12 (14-21 GeV) have overlapping truth pT windows (14-15 GeV). The `cluster_ET_upper` cap (18 vs 23 GeV) prevents double-counting at reco level.

6. **Hardcoded paths**: JES calibration file (`/sphenix/user/hanpuj/.../output_forshuhang.root`), MBD t0 correction (`MbdOut.corr`), and yaml-cpp library path are hardcoded.

7. **`reweight` forced to 0**: In `CalculatePhotonYield.C` line 836, `reweight = 0;` overrides the config value. The comment says "we move the reweight handling upstream".

8. **jetevents factor**: In `CalculatePhotonYield.C`, `float jetevents = 0.3555 * 1E7` applies a trigger efficiency factor of 0.3555 to the total jet events. This factor must match the actual trigger efficiency for the selected trigger.
