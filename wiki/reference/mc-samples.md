# MC Samples

## Signal Samples (Prompt Photon MC)

| Sample | Truth pT Window | Cross-section (pb) | Weight (/ photon20) | Path |
|--------|-----------------|--------------------|--------------------|------|
| photon5 | 0--14 GeV | 146,359.3 | 1122.0 | `sim/run28/photon5/` |
| photon10 | 14--30 GeV | 6,944.675 | 53.24 | `sim/run28/photon10/` |
| photon20 | 30+ GeV | 130.4461 | 1.0 (ref) | `sim/run28/photon20/` |

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

| Sample | Purpose | Mixing Fraction |
|--------|---------|-----------------|
| photon10_double | Full GEANT double-interaction signal | 18.7% (0 mrad), 7.2% (1.5 mrad) |
| jet12_double | Full GEANT double-interaction background | Same fractions |

Used in the two-pass blending pipeline (`run_showershape_double.sh`).

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
