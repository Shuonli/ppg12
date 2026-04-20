# MC Samples

## Signal Samples (Prompt Photon MC)

| Sample | Truth pT Window | Cross-section (pb) | Weight (/ photon20) | Path |
|--------|-----------------|--------------------|--------------------|------|
| photon5 | 0--14 GeV | 146,359.3 | 1122.0 | `sim/run28/photon5/` |
| photon10 | 14--22 GeV | 6,944.675 | 53.24 | `sim/run28/photon10/` |
| photon20 | 22+ GeV | 130.4461 | 1.0 (ref) | `sim/run28/photon20/` |

Note: before 2026-04-20 the photon10/photon20 boundary was at 30 GeV. Lowered to 22 GeV after empirically verifying (a) photon20 is fully populated above its pT-hat turn-on by pT ≥ 22 GeV, (b) photon10/photon20 weighted spectra agree to within 2-5% over [22, 30] GeV, and (c) photon20 has ~7× more events/bin than photon10 in that range — so using photon20 for pT ≥ 22 delivers ~7× better MC template statistics. The 22 GeV edge also aligns with an existing `plotcommon.h:ptRanges` bin edge (old 30 GeV cut split bin [28, 32] across two samples).

All signal samples combined via MergeSim.C (photon5 + photon10 + photon20).

## Background Samples (QCD Jet MC)

### Main Pipeline (MergeSim uses 6)

| Sample | Truth Jet pT | Cross-section (pb) | Weight (/ jet50) | Cluster ET Upper | Path |
|--------|-------------|--------------------|--------------------|-----------------|------|
| jet5 | 7--9 GeV | 1.3878e+08 | 1.898e+07 | 10 GeV | `sim/run28/jet5/` |
| jet8 | 9--14 GeV | 1.15e+07 | 1.573e+06 | 15 GeV | `sim/run28/jet8/` |
| jet12 | 14--21 GeV | 1.4903e+06 | 2.038e+05 | 23 GeV | `sim/run28/jet12/` |
| jet20 | 21--32 GeV | 6.2623e+04 | 8,565 | 35 GeV | `sim/run28/jet20/` |
| jet30 | 32--42 GeV | 2.5298e+03 | 346.0 | 45 GeV | `sim/run28/jet30/` |
| jet40 | 42--100 GeV | 1.3553e+02 | 18.54 | 100 GeV | `sim/run28/jet40/` |

### BDT Training (BDTinput.C uses 5)

| Sample | Used in training |
|--------|-----------------|
| jet10 | Yes (training only) |
| jet15 | Yes (training only) |
| jet20 | Yes (training + pipeline) |
| jet30 | Yes (training + pipeline) |
| jet50 | Yes (training only) |

BDTinput.C uses jet10/15/20/30/50. MergeSim uses jet5/8/12/20/30/40. The overlap is jet20 and jet30.

### Additional Samples

| Sample | pT Window | Purpose |
|--------|-----------|---------|
| jet10 | 10--15 GeV | BDT training |
| jet15 | 15--20 GeV | BDT training |
| jet50 | 52--100 GeV | BDT training |
| jet70 | high pT | Exploratory |
| mb | minimum bias | QA |

## Double-Interaction Samples

Eight SI/DI pairs are now produced (Apr 2026 reprocess). DI samples inherit their SI partner's cross-section and truth pT window; `efficiencytool/CrossSectionWeights.h:62-134` treats each `{sample}_double` / `{sample}_nom` alias as identical to the base SI sample (with two exceptions noted below).

| DI sample | Truth pT / jet pT window | Cross-section (pb) | Path |
|-----------|--------------------------|--------------------|------|
| photon5_double | 0--14 GeV photon | 146,359.3 | `sim/run28/photon5_double/` |
| photon10_double | 14--22 GeV photon | 6,944.675 | `sim/run28/photon10_double/` |
| photon20_double | 22+ GeV photon | 130.4461 | `sim/run28/photon20_double/` |
| jet8_double | 9--14 GeV jet | 1.15e+07 | `sim/run28/jet8_double/` |
| jet12_double | 14--21 GeV jet | 1.4903e+06 | `sim/run28/jet12_double/` |
| jet20_double | 21--32 GeV jet | 6.2623e+04 | `sim/run28/jet20_double/` |
| jet30_double | 32--42 GeV jet | 2.5298e+03 | `sim/run28/jet30_double/` |
| jet40_double | 42--100 GeV jet | 1.3553e+02 | `sim/run28/jet40_double/` |

All DI aliases now inherit their SI partner's exact truth-pT window. Prior to 2026-04-20, `photon10_double`/`photon10_nom` carried `[10, 100]` and `jet12_double`/`jet12_nom` carried `[10, 100]` (legacy behaviour from before the DI expansion), causing double-counting against adjacent samples in the DI-blended truth spectrum; fixed in `CrossSectionWeights.h:62-65, 99-102` (see `reports/xsec_weight_audit.md`).

Naming convention: `{sample}_double` = DI variant; `{sample}_nom` = alias for the SI variant used in DI blending contexts (equivalent to `{sample}` in `GetSampleConfig`).

Cluster-weighted DI mixing fractions (applied as `mix_weight` in `ShowerShapeCheck.C` and `RecoEffCalculator_TTreeReader.C`): 22.4% at 0 mrad, 7.9% at 1.5 mrad. See [Pileup Fractions](pileup-fractions.md) for the calculation.

**Skipped samples:**
- `jet5_double` — not produced. `jet5` is kept in the SI pipeline without a DI partner.
- `jet50_double` — produced but unused: `jet50` is a BDT-training-only sample, not in `FunWithxgboost/mc_samples.list` or `MergeSim.C`, so there is no SI partner in the main pipeline.

Used in the single-pass truth-vertex-reweight showershape DI pipeline (`submit_showershape_di.sub`; see [Double-Interaction Efficiency](../concepts/double-interaction-efficiency.md)).

## File Locations

Base directory: `/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/`

Each sample directory contains:
- `Fun4All_run_sim.C` -- macro
- `run_condor.sh` -- condor submission
- `condorout/combined.root` -- merged output (the input to downstream stages)

BDT-scored versions: `/sphenix/user/shuhangli/ppg12/FunWithxgboost/{sample}/bdt_split.root`

## Generated Events

Default: `nsimevents = 1E7` per sample (from CrossSectionWeights.h).

## Weight Application

Weights are applied in `RecoEffCalculator_TTreeReader.C` at histogram fill time:
```cpp
cross_weight = PPG12::GetSampleConfig(filetype).weight * mix_weight;
weight = cross_weight * vertex_weight;
```

Photon weights are relative to photon20cross. Jet weights are relative to jet50cross. After weighting, TFileMerger addition in MergeSim.C correctly produces the combined weighted spectrum.

## Cross-Section Source

All cross-section values are from PYTHIA8 generation. They are centralized in `efficiencytool/CrossSectionWeights.h` and accessed via `PPG12::GetSampleConfig(filetype)`.

## Run Conditions

MC is produced at run 28 conditions. CDB global tag: `MDC2`. This corresponds to the detector configuration during early Run 24 data taking.

## Sim-Specific DST Streams

Each sim sample reads 4 DST streams:
1. `G4Hits` -- GEANT4 hit information
2. `DST_CALO_CLUSTER` -- calo cluster containers
3. `DST_MBD_EPD` -- MBD and EPD containers
4. `DST_TRUTH_JET` -- truth-level jet containers
