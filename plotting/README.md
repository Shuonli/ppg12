# Plotting Macros

ROOT/C++ plotting macros for the sPHENIX direct photon cross section analysis (PPG12). All macros include `plotcommon.h` and are run via the ROOT interpreter.

## Common Utilities

| File | Description |
|------|-------------|
| `plotcommon.h` | Common header: defines pT bins (8--35 GeV, 10 bins), frame histograms, legend strings, `init_plot()`, and `calcDelta()` helper. |
| `BlairUtils.C` | Plot styling utilities: `TGraphErrors` division with error propagation, error-band construction, `TH1`-to-`TGraph` conversion, text/marker drawing helpers. |
| `draw.C` | Generic drawing functions: `draw_1D_single_plot`, `draw_1D_multiple_plot_ratio` with configurable axes, legends, and ratio pads. |
| `sPhenixStyle.C` | sPHENIX style settings (included via `plotcommon.h`). |

## Final Results

| File | What It Plots | Input Files |
|------|---------------|-------------|
| `plot_final.C` | Differential cross section (d sigma / dET) vs ET. Data compared to JETPHOX NLO (with scale-variation band) and PYTHIA8 truth. Ratio panels: data/NLO and NLO/data. | `Photon_final_bdt_base.root`, `syst_sum.root`, `jetPHOX_{05,10,20}.root`, `Photon_final_nom_mc.root` |
| `plot_final_backup.C` | Backup/earlier version of the final cross section plot. | Same as `plot_final.C` |
| `plot_final_backup250327_v1.C` | Another archived version of `plot_final.C`. | Same as `plot_final.C` |

## Efficiency

| File | What It Plots | Input Files |
|------|---------------|-------------|
| `plot_efficiency.C` | Reconstruction, isolation, identification, total, MBD trigger, and combined (reco x ID x iso) efficiencies vs truth ET. | `MC_efficiency_bdt_base_v3E_2.root`, `Photon_final_nom.root` |
| `plot_photonjeteff.C` | Efficiency of finding a back-to-back jet for signal photon clusters, as a function of photon ET for different jet pT thresholds. | `MC_efficiency_bdt_none.root` |

## Purity

| File | What It Plots | Input Files |
|------|---------------|-------------|
| `plot_purity.C` | Photon purity vs ET with and without signal-leakage correction, compared to MC truth purity. | `Photon_final_nomtestv3.root`, `Photon_final_nomtestv3_mc.root` |
| `plot_purity_sim.C` | MC-only purity evaluation using PYTHIA signal and inclusive jet samples. | `MC_efficiency_nom.root`, `MC_efficiency_jet_nom.root` |
| `plot_purity_sim_selection.C` | Same as `plot_purity_sim.C` but parameterized by selection suffix to compare different cut variations. | `MC_efficiency_{suffix}.root`, `MC_efficiency_jet_{suffix}.root` |
| `plot_purity_selection.C` | Data purity under different selection criteria (parameterized by suffix). | `Photon_final_{suffix}.root` |
| `plot_pid.C` | Particle species composition (pion, electron, photon, eta, etc.) in tight-iso clusters from jet MC. | `MC_efficiency_jet_nom.root` |

## ABCD Sideband Method

| File | What It Plots | Input Files |
|------|---------------|-------------|
| `plot_sideband.C` | Full ABCD sideband analysis: ET spectra in A/B/C/D regions, NPB-weighted spectra, B/A C/A D/A ratios (data vs signal MC vs jet MC), signal leakage fractions, isolation ET distributions per pT bin. | `data_histo_468.root`, `MC_efficiency_468.root`, `MC_efficiency_jet_nom.root` |
| `plot_sideband_selection.C` | Same ABCD analysis parameterized by a selection suffix for systematic variation studies. | `data_histo_{suffix}.root`, `MC_efficiency_{suffix}.root`, `MC_efficiency_jet_{suffix}.root` |
| `plot_sideband_sim.C` | ABCD sideband ratios using only MC (jet MC as pseudo-data, signal MC for leakage). | `MC_efficiency_jet_nom.root`, `MC_efficiency_nom.root` |
| `plot_sideband_sim_selection.C` | Same as `plot_sideband_sim.C` parameterized by suffix. | `MC_efficiency_jet_{suffix}.root`, `MC_efficiency_{suffix}.root` |
| `plot_sideband_compare_left_right.C` | Compares ABCD histograms between two different analysis configurations (right vs left) to check sideband stability. | Two sets of `Photon_final_{suffix}.root` |
| `plot_SB.C` | Signal-to-background ratio (S/B) vs ET using signal and jet MC shower-shape files, scaled by cross sections. | `MC_efficiencyshower_shape_signal.root`, `MC_efficiencyshower_shape_jet.root` |

## Unfolding

| File | What It Plots | Input Files |
|------|---------------|-------------|
| `plot_response.C` | Response matrix (reco ET vs truth ET) as 2D color plot, with and without reweighting. | `Photon_final_nom.root`, `Photon_final_nr.root` |
| `plot_unfold_iter.C` | Bayesian unfolding convergence: relative change in unfolded spectrum vs iteration number. Uses RooUnfold. | `Photon_final_nom.root` (response matrix) |
| `Closure.C` | Unfolding closure test: full-sample and half-sample closure using RooUnfoldBayes. Compares truth vs reco vs unfolded spectra. Reads config from YAML. | `MC_response_photon10.root`, `config.yaml` |

## Shower Shapes & BDT

| File | What It Plots | Input Files |
|------|---------------|-------------|
| `plot_showershapes.C` | Shower shape variable distributions (probability, CNN prob, energy ratios e17/e77, e37/e77, HCAL contributions) comparing data vs signal MC vs jet MC, per pT bin. | `data_histoshower_shape_.root`, MC shower shape files, `config_nom.yaml` |
| `plot_showershapes_selections.C` | Extended shower shape comparison across different selection variations. | Multiple shower shape ROOT files |
| `plot_vertex_check.C` | Signal (NPB) and background (BDT) efficiency stability vs vertex z shifts. 2D correlation plots and efficiency curves for different vertex positions. | `MC_vertex_check_photon.root`, `MC_vertex_check_jet.root` |

## Cluster & Event Timing

| File | What It Plots | Input Files |
|------|---------------|-------------|
| `plot_cluster_timing.C` | Cluster-MBD time difference (delta-T) distributions for different selections (common, tight-iso, tight-noniso, NPB-scored) projected per pT bin. | `cluster_time_analysis_data.root` |
| `plot_cluster_timing_compare_leading.C` | Timing analysis comparing leading vs sub-leading clusters. | Cluster timing ROOT files |
| `plot_cluster_mbd_time.C` | Cluster-MBD time correlation distributions. | Cluster timing ROOT files |
| `plot_clustertime.C` | Cluster timing quality distributions. | Cluster timing ROOT files |
| `plot_tower_timing.C` | Detailed tower-level timing analysis and performance. | Tower timing ROOT files |
| `plot_leading_tower_mean_time_2d3d.C` | 2D and 3D distributions of leading tower mean time, comparing data/MC and different selections. | Cluster timing ROOT files |
| `plot_npb_time_bkgsub.C` | NPB score vs cluster-MBD time with background subtraction per pT bin. | Cluster timing ROOT files with NPB filtering |
| `plot_time_energy_corr.C` | Time-energy correlation in the calorimeter. | `time_energy_corr.root` |
| `plot_time_energy_corr_fixtime.C` | Time-energy correlation with fixed time window. | `time_energy_corr.root` |

## Special Studies

| File | What It Plots | Input Files |
|------|---------------|-------------|
| `plot_reweight.C` | Purity-corrected data spectrum vs PYTHIA signal spectrum (both shape-normalized), and their ratio used as the reweighting function. | `Photon_final_nom.root` |
| `plot_combine.C` | Ratio of leading photon pT spectra from different pT-hat samples (photon20/photon10, photon10/photon5) to validate cross-section weighting. | `MC_efficiency_photon{5,10,20}_mbdeffup.root` |
| `plot_double_interaction.C` | Double interaction / pileup analysis: ABCD yields with NPB < 0.5, fraction of events with NPB < 0.5 vs ET. | Double interaction check ROOT files (data, signal, jet) |
| `plot_rbrQA.C` | Run-by-run QA: cluster counts (common, tight-iso, tight-noniso, nontight-iso, nontight-noniso) per run with 2.5th/95th/97.5th percentile thresholds. | `rbrQA.root`, `rbrQA_nc.root` |
| `plot_vertex.C` | Vertex z distributions for data and MC (shape-normalized), and data/MC ratio fitted with a Gaussian. | `MC_efficiency_nom.root`, `data_histo_nom.root` |
| `plot_xjg.C` | Photon-jet momentum imbalance (x_Jgamma) 2D and 1D distributions for signal MC, inclusive jet MC, and data. | `MC_efficiency_nom.root`, `MC_efficiency_jet_nom.root`, `data_histo_nom.root` |
| `plot_truth_iso.C` | Truth- and reco-level isolation energy distributions per pT bin for signal and background MC. | `MC_efficiencyshower_shape_signal.root`, `MC_efficiencyshower_shape_jet.root` |
| `plot_satdiff.C` | Impact of saturation correction: ratios of ABCD region counts and common cluster ET between old and new calibrations. | `data_histo_nom468.root`, `data_histo_468newcalib.root` |
| `plot_mbd_sigma_efficiency.C` | MBD average sigma distribution and selection efficiency as a function of the MBD sigma cut. | `data_histoshower_shape_.root` |
| `CONF_plots.C` | Conference-style combined plots: isolation ET spectra for tight/non-tight clusters with Crystal Ball fits, per pT bin range. | `data_histo_{tune}.root`, `MC_efficiency_{tune}.root` |

## Systematic Uncertainties

All `syst_*.C` macros compute systematic deviations by comparing a nominal analysis to one or more variations. Each produces a ROOT file in `rootFiles/` containing `h_dev_low`, `h_dev_high`, `h_dev_rel_low`, `h_dev_rel_high` histograms.

| File | Systematic Source | Variations Compared |
|------|-------------------|---------------------|
| `calcSyst.C` | **Combined**: sums all sources in quadrature, adds luminosity uncertainty (16.6 pb^-1). Reads individual `syst_*.root` files and outputs `syst_sum.root`. | All below |
| `syst_purity.C` | **Purity** (envelope): combines isolation, non-tight fraction, and fit variations. | `syst_iso.root`, `syst_ntf.root`, `syst_fit.root` |
| `syst_iso.C` | **Isolation threshold**: effect of varying the isolation cut boundary. | nom vs iso20 vs iso05 |
| `syst_nt.C` | **Non-tight definition**: effect of changing the tight/non-tight discrimination. | nom vs id1 |
| `syst_ntf.C` | **Non-tight fraction**: ABCD purity extraction method variation. | nom vs ntf1 vs ntf3 |
| `syst_fit.C` | **Fit method**: ABCD fit procedure variations. | nom vs fit variations |
| `syst_escale.C` | **Energy scale**: calorimeter calibration uncertainty. | nom vs escale (symmetric) |
| `syst_eres.C` | **Energy resolution**: smearing uncertainty. | nom vs eres |
| `syst_eff.C` | **Efficiency**: tight selection and fudge-factor variations. | nom vs tight, fudge |
| `syst_tight.C` | **Tight selection**: shower shape cut variation. | nom vs tight variation |
| `syst_fudge.C` | **Fudge factor**: MC correction factor variation. | nom vs fudge |
| `syst_mbd.C` | **MBD trigger efficiency**: up/down variations. | nom vs mbdeffup vs mbdeffdown |
| `syst_nor.C` | **Normalization / unfolding**: unfolding procedure variation. | nom vs nor |
| `syst_vertex.C` | **Vertex**: vertex reconstruction variation. | nom vs vertex |
| `syst_CNNvar.C` | **CNN classifier**: CNN-based selection variation. | nom vs CNN variation |
| `examplecalcSyst.C` | Example/template for systematic calculations. | -- |

## NLO Comparison (subdirectory `sphenix_nlo/`)

| File | Description |
|------|-------------|
| `sphenix_nlo/nlo_plot.C` | Plots NLO theory predictions (JETPHOX) for comparison. |
| `sphenix_nlo/nlo_plot_count.C` | NLO prediction bin counts / validation. |

## Output

Plots are saved to `figures/` (PDF/PNG) or `../PPG12-analysis-note/Figures/` for the analysis note.
