# Donut iso variants: inner-R truth-purity scan

**Date**: 2026-04-19
**Scope**: cross-check only (all 7 variants have `syst_type=None`, `syst_role=None`)

## Motivation

Quantify the truth-matched photon purity (`g_purity_truth` from `CalculatePhotonYield.C`) under a new isolation definition: an **annular donut** built from tower ETs with configurable inner radius. Two outer-cone families, both at R_outer=0.4 with 120 MeV tower threshold (EMCal+HCalIn+HCalOut):

- **donutFull** (`use_topo_iso=4`): raw cone `cluster_iso_04` minus raw inner `cluster_iso_{005,0075,02}` — cluster ET subtracted as scalar, the two scalar subtractions cancel in the difference.
- **donutExcl** (`use_topo_iso=5`): exclusion cone `cluster_iso_excl_04` minus exclusion inner `cluster_iso_excl_{005,0075,01,02}` — cluster's own EMCal towers removed by ownership at tower-level, not scalar.

The `donutFull_01` variant was deliberately skipped because `cluster_iso_01` (raw family at R=0.10) does not exist in the current slimtree. Logged in `reports/anatreemaker_reprocess_TODO.md` (full copy at `/sphenix/u/shuhang98/ppg12_anatreemaker_reprocess_TODO.md`).

## Method

Two-pass pipeline per variant:

1. **Pass 1** — run full `RecoEffCalculator_TTreeReader + MergeSim + CalculatePhotonYield` with loose iso cut (`reco_iso_max_b=100, s=0`, effectively no cut). Produces `MC_efficiency_bdt_donut*.root` with the unbiased 2D `h_singal_reco_isoET_0`.
2. **Cut derivation** — `FindETCut.C` now takes `(var_type, target_eff)` and emits `DERIVED_ISO_CUT var_type=... intercept=... slope=...`. Driver `derive_iso_cuts.py` loops over the 7 variants, parses, writes `derived_iso_cuts.yaml`.
3. **Pass 2** — regenerate configs via `make_bdt_variations.py` with the derived (b, s), re-run pipeline via `oneforall_donut.sub` on condor.

Both passes on condor cluster 2701653 (Pass 1) and 2701654 (Pass 2), ~16-20 min each for all 7 variants in parallel.

## Derived 80% iso cuts

| Variant | b (intercept) | s (slope) |
|---|---|---|
| donutFull_005 | 0.1089 | 0.0224 |
| donutFull_0075 | 0.2642 | 0.0061 |
| donutFull_02 | 0.2305 | -0.0012 |
| donutExcl_005 | 0.2457 | 0.0111 |
| donutExcl_0075 | 0.2774 | 0.0053 |
| donutExcl_01 | 0.2921 | 0.0017 |
| donutExcl_02 | 0.2305 | -0.0012 |

Note: `donutFull_02` and `donutExcl_02` landed on identical fit bins. Expected: at R_inner=0.2 the EMCal ownership-exclusion and the scalar ET subtraction produce nearly identical annulus integrals (the cluster footprint is fully contained within R<0.2).

## Truth purity results

Values from `g_purity_truth` in `Photon_final_bdt_<var>_mc.root` (MC-as-data; the non-`_mc` file returns zero because `_signal` histograms are absent from data).

| Variant | y(9 GeV) | y(21 GeV) | y(34 GeV) |
|---|---|---|---|
| nom (reference) | 0.581 | 0.743 | 0.797 |
| donutFull_005 | 0.684 | 0.786 | 0.849 |
| donutFull_0075 | 0.608 | 0.805 | 0.873 |
| donutFull_02 | 0.496 | 0.728 | 0.805 |
| donutExcl_005 | 0.648 | 0.803 | 0.858 |
| donutExcl_0075 | 0.606 | 0.804 | 0.874 |
| donutExcl_01 | 0.577 | 0.793 | 0.869 |
| donutExcl_02 | 0.496 | 0.728 | 0.805 |

### Observations

- **All donut variants with R_inner <= 0.10 recover more signal purity than nominal** at all pT bins. donutFull_005 at 34 GeV reaches 0.85 vs nominal 0.80.
- **Trend with inner R** (Excl family): y(9 GeV) decreases smoothly 0.648 -> 0.606 -> 0.577 -> 0.496 as R_inner tightens from 0.05 -> 0.20. Physical: larger R_inner admits less annulus area, looser iso cut in absolute terms -> more background contamination.
- **donutFull vs donutExcl at matched R_inner**: donutFull variants have marginally higher purity at low pT (smaller annulus, since the -ET scalar subtracts a bit more than the cluster footprint); convergence at high pT.
- **donutFull_02 == donutExcl_02** to 4 decimals across all 12 bins -> iso cut equivalence is exact, not just approximate.

## Plots

- `plotting/figures/purity_sim_donutFull_scan.pdf` — 3 curves (R=0.05, 0.075, 0.20)
- `plotting/figures/purity_sim_donutExcl_scan.pdf` — 4 curves (R=0.05, 0.075, 0.10, 0.20)

Both show rising purity from ~0.5-0.7 at 8 GeV to ~0.80-0.87 at 30+ GeV. Cosmetics review: PASS with two WARNs: (a) x-axis extends to 40 GeV per shared `frame_et_rec` convention (8-40), not 8-36; (b) the "Donut iso" header text overlaps curves near pT=20 GeV. Both are cosmetic and non-blocking.

## Integration with `make_selection_report.py`

Added a top-level "Donut inner-R truth-purity scan" section before per-config sections, containing the two overlay PDFs as subfloats. Each of the 7 donut configs also appears in its own per-config section via existing auto-discovery (`config_bdt_donut*.yaml` glob).

End-to-end test: `python3 plotting/make_selection_report.py --output-dir /tmp/donut_selection_report` produces 66 sections, 1483 figures, both donut overlay PDFs copied.

## Caveats and follow-ups

1. **All-calo outer cone includes HCal** — the donut is not EMCal-only. A pure-EMCal donut requires new EMCal-only inner-ring branches (`cluster_iso_005_emcal`, `_0075_emcal`, `_01_emcal`, `_02_emcal`) at matching threshold (120 MeV or no threshold). Logged in the anatreemaker reprocess TODO.
2. **`donutFull_01` skipped** pending upstream `cluster_iso_01` branch (all-calo, 120 MeV, R<0.1).
3. **Degraded mirrors** (`Cluster_rbr.C`, `EtaMigrationStudy.C`, `compare_energy_response.C`) silently use the full outer cone without inner subtraction when run with a donut variant. Do not run donut variants through these utilities expecting donut semantics.
4. **MC iso scale/shift** (`mc_iso_scale`/`shift`) are set to 1.0/0.0 for all donut variants because the nominal scale/shift was tuned on the full-cone topo04 iso. Retuning for donut would be needed if these ever become systematic variants.
5. **Disk quota on `/sphenix/user/shuhangli` (~5 TB)** was hit during this session. Large writes under `reports/` failed with EDQUOT. TODO: prune stale `efficiencytool/results/` entries before next deploy.

## Reviewer pass matrix

| Artifact | 2a physics | 2b cosmetics | 2c re-derivation | 3.5 fix-validator |
|---|---|---|---|---|
| use_topo_iso=4,5 switch | PASS (2 WARNs on utility mirrors) | N/A | N/A | PASS (pipeline produced sane MC_efficiency and Photon_final outputs) |
| donutFull_scan.pdf | N/A | PASS + 2 WARNs (conv x-range, header overlap) | PASS | N/A |
| donutExcl_scan.pdf | N/A | PASS + 2 WARNs | PASS | N/A |
| 80% iso cuts | N/A | N/A | PASS (derived cuts physically reasonable, lower magnitude than 70 MeV innerRtower as expected from 120 MeV threshold) | PASS |
| g_purity_truth (7 variants) | N/A | N/A | PASS (y in [0.49, 0.87], smooth trend with R_inner) | PASS |

## Files touched

- `efficiencytool/RecoEffCalculator_TTreeReader.C` — added case 4/5 + 7 branch decls
- `efficiencytool/ShowerShapeCheck.C` — mirror
- `efficiencytool/compare_efficiency_deltaR.C` — mirror
- `efficiencytool/Cluster_rbr.C` — degraded mirror
- `efficiencytool/EtaMigrationStudy.C` — degraded mirror
- `efficiencytool/compare_energy_response.C` — degraded mirror
- `efficiencytool/make_bdt_variations.py` — 7 new VARIANTS + derived cuts
- `efficiencytool/FindETCut.C` — parameterized entry point
- `efficiencytool/derive_iso_cuts.py` — new driver (96 lines)
- `efficiencytool/derived_iso_cuts.yaml` — output of derive driver
- `efficiencytool/oneforall_donut.sub` — new condor submit
- `plotting/plot_purity_sim_donut_scan.C` — new overlay macro (131 lines)
- `plotting/make_selection_report.py` — added donut section + all_src_paths entries
- `plotting/figures/purity_sim_donutFull_scan.pdf` — output
- `plotting/figures/purity_sim_donutExcl_scan.pdf` — output
- 7 `efficiencytool/config_bdt_donut*.yaml` — generated configs (auto-regen from make_bdt_variations.py)
