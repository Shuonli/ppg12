# Full Run Range Cross-Section Pipeline

The cross-section can be measured per-period (0 mrad only, 1.5 mrad only) or for the
**full run range** (47289-54000, ~48.9 pb⁻¹). The full-range workflow combines
per-period MC outputs into an all-range MC via lumi-weighted aggregation.

## Architecture

Per-period MC carries different physics conditions (DI fraction, truth-vertex
distribution) that must stay correct per crossing angle. The all-range MC is
produced by:

1. Running the standard per-period pipeline (RecoEffCalculator + MergeSim) for
   each crossing angle, with a **per-event lumi weight** baked into the MC fills.
2. Combining the pre-scaled per-period MC outputs via plain `TFileMerger` (hadd).

The per-event lumi weight is `lumi_weight = lumi / lumi_target`, applied in
`RecoEffCalculator_TTreeReader.C` to MC events only (gated by `issim`). When
`lumi_target = lumi` (default), the weight is 1 and per-period analysis is
unchanged. When `lumi_target = sum(L_periods)`, each per-period MC histogram is
pre-scaled to its contribution to the all-range expectation, and a plain hadd
across periods reproduces the all-range MC.

## Configuration

### Plan B restructure (April 2026): nominal is all-range, variants auto-expand

The nominal cross-section IS the all-range result. `config_bdt_nom.yaml` is the
all-range base + nominal analysis target (BDT cuts are the AxX anchor pattern:
tight 0.80→0.70, nt_hi 0.80→0.70, nt_lo 0.70→0.40). Each physics variant
generates 3 on-disk configs via `make_bdt_variations.py`:

| Filename | Role | Source |
|----------|------|--------|
| `config_bdt_X.yaml` | All-range analysis target (variant X). Reads hadd of the two merge-feeders. | Variant's own overrides applied to the all-range base. |
| `config_bdt_X_0rad.yaml` | 0 mrad merge-feeder. Per-event `lumi_weight = lumi/lumi_target` pre-scales MC. | Variant overrides + `PER_PERIOD_OVERRIDES["0rad"]`. |
| `config_bdt_X_1p5mrad.yaml` | 1.5 mrad merge-feeder. Same. | Variant overrides + `PER_PERIOD_OVERRIDES["1p5mrad"]`. |

### Canonical configs

| File | lumi (pb⁻¹) | lumi_target (pb⁻¹) | Role |
|------|-------------|--------------------|------|
| `config_bdt_nom.yaml` | 48.9309 | 48.9309 | All-range nominal analysis target (base for `make_bdt_variations.py`). truth_vertex_reweight_on=0 (merged MC is pre-reweighted). |
| `config_bdt_nom_0rad.yaml` | 32.6574 | 48.9309 | 0 mrad merge-feeder for nom |
| `config_bdt_nom_1p5mrad.yaml` | 16.2735 | 48.9309 | 1.5 mrad merge-feeder for nom |
| `config_bdt_allz.yaml` | 64.3718 | 64.3718 | All-range all-z cross-check (beam-delivered lumi, vertex_cut_truth=9999) |
| `config_bdt_allz_0rad.yaml` | 47.2076 | 64.3718 | 0 mrad allz merge-feeder |
| `config_bdt_allz_1p5mrad.yaml` | 17.1642 | 64.3718 | 1.5 mrad allz merge-feeder |

`truth_vertex_reweight_on=0` on the all-range configs because the merged MC
carries the per-period reweight baked into per-event weights; applying again
would double-count. Merge-feeders have `truth_vertex_reweight_on=1` with the
period-specific reweight file.

### Systematic variants

`make_bdt_variations.py:apply_overrides` auto-defaults `lumi_target := lumi` for
every variant unless the variant entry explicitly overrides. This prevents the
base config's `lumi_target=48.9309` from silently propagating to systematic
variants, which would otherwise produce wrong per-period cross-sections (MC
scaled to all-range but compared against per-period data).

## Workflow

```bash
cd efficiencytool

# 1. Per-period RecoEffCalculator for each MC sample × period (long, on condor).
#    With lumi_target set in the per-period configs, MC is pre-scaled.
condor_submit submit_nosplit_tree.sub  # or your equivalent submit infrastructure

# 2. Per-period MergeSim (one job per period config).
root -l -b -q 'MergeSim.C("config_bdt_0rad.yaml")'
root -l -b -q 'MergeSim.C("config_bdt_nom.yaml")'

# 3. Cross-period merge (TFileMerger / hadd). Replaces the lumi-weighted
#    histogram-level scale that lived in the old merge_periods.C.
bash merge_periods.sh config_bdt_0rad.yaml config_bdt_nom.yaml config_bdt_all.yaml

# 4. Data-side RecoEffCalculator for the all-range (one job).
root -l -b -q 'RecoEffCalculator_TTreeReader.C("config_bdt_all.yaml", "data")'

# 5. CalculatePhotonYield × 2 (MC closure + data cross-section).
root -l -b -q 'CalculatePhotonYield.C("config_bdt_all.yaml", true)'
root -l -b -q 'CalculatePhotonYield.C("config_bdt_all.yaml", false)'
```

Steps 3–5 are wrapped by `oneforall.sh` (which is queued by `oneforall.sub`).

## `oneforall.sh` dispatch

`oneforall.sh` inspects the config filename:
- **`*_all.yaml`** → invokes `merge_periods.sh` with auto-derived per-period
  configs (`_all.yaml → _0rad.yaml + _1p5mrad.yaml`, with the bare
  `config_bdt_all.yaml` mapping to `config_bdt_nom.yaml` for the historical
  `_nom` partner), then runs `RecoEffCalculator_TTreeReader.C(config, "data")`
  to produce the all-range `data_histo_{var_type}.root` (since `merge_periods.sh`
  handles only MC), then `CalculatePhotonYield × 2`.
- **anything else** → invokes `MergeSim.C` (the per-period merge of per-sample
  MC files), then CalculatePhotonYield × 2.

`set -eo pipefail` at the top of `oneforall.sh` aborts on merge failure so
CalculatePhotonYield never runs against stale or empty data. `set -u` is
enabled AFTER sourcing `setup.sh` because `sphenix_setup.sh` references
unbound vars (e.g. `PGHOST`) which would trip strict mode.

## Production HTCondor submission (two-phase)

End-to-end production run uses two submit files, each gated on the previous:

### Phase 1: `oneforall_tree_double.sub` (per-config TTree → merge → yield)

Queues one job per `config_bdt_*.yaml`. Each job runs
`oneforall_tree_double_dispatch.sh` → derives DOUBLE_FRAC from
`analysis.run_min` (0.224 for run_min<51274, else 0.079) → invokes
`oneforall_tree_double.sh`. All-range configs (`var_type` ending in `_all`)
exit cleanly via the dispatch — they have no per-sample MC.

Memory: **16 GB** per job. The DI pipeline launches 17 parallel ROOT
processes per pass (8 SI + 8 DI MC + 1 data), each with TTreeReader buffers,
peaking ~8–12 GB. The historical 6 GB request was too tight (jobs went on
hold at ~6.1 GB).

```bash
condor_submit oneforall_tree_double.sub
```

### Phase 2: `oneforall_all.sub` (all-range merge + yield, 2 jobs)

Queues exactly the 2 `*_all.yaml` configs. Each invokes `oneforall.sh`
which dispatches to `merge_periods.sh` + `CalculatePhotonYield × 2`.

Memory: **4 GB** per job (CalculatePhotonYield is single-threaded with a
modest footprint).

```bash
# Wait for Phase 1 to complete, then:
condor_submit oneforall_all.sub
```

### Other submit files (legacy / specialized)

- `oneforall.sub` (2 GB, 88 jobs) — runs `oneforall.sh` for every config. Only
  correct when per-sample MC is already current; otherwise produces stale-cut
  results. Mostly superseded by Phase 1 + Phase 2.
- `oneforall_tree.sub` (6 GB) — single-interaction-only TTree reader (no DI
  mix). Legacy; do not use for production.

## Order dependency

`oneforall.sub` queues all 99 configs with no DAG dependency. The all-range
configs (`config_bdt_all.yaml`, `config_bdt_allz_all.yaml`) need their
per-period merge-feeder MC outputs to exist before `merge_periods.sh` runs. In
practice: submit per-period configs first, wait for completion, then submit
all-range configs. Or use DAGMan if you re-run frequently.

## Per-period standalone analysis

The merge-feeder configs (`config_bdt_0rad.yaml`, `config_bdt_nom.yaml`, etc.)
have `lumi_target ≠ lumi`, so their MC is pre-scaled to the all-range
expectation. Running them as standalone per-period analyses gives the WRONG
cross-section because data is normalized by the per-period lumi but MC is
scaled to a different target.

For a per-period standalone cross-check, create a one-off config that omits
`lumi_target` (or sets `lumi_target=lumi` explicitly) so `lumi_weight=1` and the
MC is unscaled. This is a deliberate trade-off documented in the wave-pipeline
discussion: "per-period MC is pre-scaled by default; per-period standalone is a
separate one-off."

## Allz cross-check

Tests the |z|<60 cm fiducial chain by widening the MBD-eff truth denominator to
all z while keeping the analysis-fiducial |z_reco|<60 cut on data and the
MBD-eff numerator. Paired with the beam-delivered (allz) lumi.

The new `vertex_cut_truth` config field controls the truth-vertex cut on the
MBD-eff denominator (line ~1785 of RecoEffCalculator_TTreeReader.C) independently
from the analysis-wide `vertex_cut`. Setting `vertex_cut_truth: 9999.0` removes
the truth cut while leaving everything else at 60 cm.

The per-period configs `allz_{0rad,1p5mrad}` are merge-feeders just like
`{0rad,nom}`. The combined allz MC is hadded into `MC_efficiency_bdt_allz_all.root`
and analyzed with `config_bdt_allz_all.yaml`.

## Files

| File | Purpose |
|------|---------|
| `efficiencytool/RecoEffCalculator_TTreeReader.C:~360` | Reads `lumi`, `lumi_target`; bakes `cross_weight *= lumi/lumi_target` for MC |
| `efficiencytool/merge_periods.C` | TFileMerger of two per-period MC outputs into combined; writes `merge_lumi_*` provenance TNamed |
| `efficiencytool/merge_periods.sh` | Wrapper invoking the ROOT macro |
| `efficiencytool/oneforall.sh` | Dispatch by `*_all.yaml` glob; runs merge_periods.sh OR MergeSim.C + CalculatePhotonYield × 2 |
| `efficiencytool/make_bdt_variations.py:apply_overrides` | Auto-defaults `lumi_target := lumi` for variants without explicit override |

## See also

- [Configuration Schema](../reference/config-schema.md) — `lumi_target`, `vertex_cut_truth` field definitions
- [Constants Sync](../reference/constants-sync.md) — luminosity sources, BDT cut canonical values
- [Stage 4: Efficiency and Yield](04-efficiency-yield.md) — per-period efficiency calculation details
