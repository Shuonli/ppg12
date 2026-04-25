# Pileup Fractions

## Overview

Double-interaction (pileup) fractions used for MC blending are computed from per-run luminosity data by `efficiencytool/calc_pileup_range.C`. The calculation derives Poisson interaction statistics from the corrected MBD trigger rate and converts event-level fractions to cluster-weighted fractions suitable for MC sample blending.

## Calculation Method

### Input

Per-run ROOT files from `analyze.C` containing:
- `nmbdc` — Poisson-corrected MBD count (full pile-up correction)
- `ttseg` — live time in seconds

Beam crossing rate: `rB = 111 bunches x 78200 Hz = 8.682 MHz` (Run-3 standard).

### Step 1: Mean collisions per crossing

Per-run Poisson-corrected mean:
```
mu_corr(run) = nmbdc / (rB * ttseg)
```

### Step 2: Triggered interaction-mix fractions

Conditioned on at least one collision (i.e., a triggered event), the probability of exactly `n` collisions:
```
P(n >= 1) = 1 - exp(-mu)
f_single  = P(n=1 | n>=1) = mu * exp(-mu) / P(n>=1)
f_double  = P(n=2 | n>=1) = mu^2/2 * exp(-mu) / P(n>=1)
f_triple+ = 1 - f_single - f_double
```

These are computed per-run and then luminosity-weighted (weight = `nmbdc` per run), because `f(mu)` is nonlinear and `<f(mu)> != f(<mu>)` (Jensen's inequality).

### Step 3: Cluster-weighted mixing fractions

A double-interaction event produces ~2x as many clusters as single, triple ~3x, etc. To correctly weight MC samples at the cluster level:
```
cw_single = 1 * f_single / (1*f_single + 2*f_double + 3*f_triple+)
cw_double = (2*f_double + 3*f_triple+) / (1*f_single + 2*f_double + 3*f_triple+)
```

Triple+ events are folded into the double fraction because only double-interaction MC samples are available.

## Results

Computed 2026-04-10 from the full good run list (`fullgoodrunlist.list`).

### 0 mrad crossing angle (runs 47289--51274)

| Quantity | Value |
|----------|-------|
| Runs processed | 607 |
| `mu_corr` (aggregate) | 0.236 |
| P(pileup) = P(>=2 \| >=1) | 12.1% |
| `f_single` (event-level) | 87.9% |
| `f_double` (event-level) | 11.1% |
| `f_triple+` (event-level) | 1.1% |
| **`cw_single`** (cluster-weighted) | **77.6%** |
| **`cw_double`** (cluster-weighted) | **22.4%** |
| Integrated luminosity | 67.7 pb^-1 |

### 1.5 mrad crossing angle (runs 51274--54000)

| Quantity | Value |
|----------|-------|
| Runs processed | 258 |
| `mu_corr` (aggregate) | 0.078 |
| P(pileup) = P(>=2 \| >=1) | 4.0% |
| `f_single` (event-level) | 96.0% |
| `f_double` (event-level) | 3.9% |
| `f_triple+` (event-level) | 0.1% |
| **`cw_single`** (cluster-weighted) | **92.1%** |
| **`cw_double`** (cluster-weighted) | **7.9%** |
| Integrated luminosity | 12.5 pb^-1 |

## Usage in Pipeline

The cluster-weighted fractions (`cw_double`) are used as `mix_weight` (DI jobs) and `1 - cw_double` (SI jobs) in `ShowerShapeCheck.C` and `RecoEffCalculator_TTreeReader.C`. For the single-pass showershape DI pipeline the values are hardcoded into the condor row lists:

```bash
# Showershape DI (single-pass, truth-vertex reweight):
cd efficiencytool
condor_submit list_file=showershape_di_jobs_0rad.list   submit_showershape_di.sub  # 0 mrad (0.224 / 0.776)
condor_submit list_file=showershape_di_jobs_1p5rad.list submit_showershape_di.sub  # 1.5 mrad (0.079 / 0.921)

# Legacy / retired (reco-vertex two-pass):
bash run_showershape_double_reco_legacy.sh config.yaml 0.224

# Auto (main cross-section pipeline; computes fraction from data):
bash run_double_auto.sh --mode crosssection 1.5mrad
```

The auto scripts run `calc_pileup_range.C` for the appropriate run range and extract `cw_double` automatically. The showershape DI row lists currently pin the fractions; regenerate them if the good-run list changes.

## How to Recompute

```bash
cd efficiencytool
root -b -q -l 'calc_pileup_range.C(47289, 51274)'   # 0 mrad
root -b -q -l 'calc_pileup_range.C(51274, 54000)'   # 1.5 mrad
```

Recompute if the good run list changes or if run range boundaries are revised.

## Files That Use These Fractions

| File | How it enters |
|------|---------------|
| `efficiencytool/showershape_di_jobs_{0rad,1p5rad}.list` | Per-row `mix_weight` column (0.224/0.776 or 0.079/0.921) consumed by `submit_showershape_di.sub` |
| `efficiencytool/run_showershape_double_reco_legacy.sh` | Legacy two-pass `DOUBLE_FRAC` default (0.224) — retired |
| `efficiencytool/oneforall_tree_double.sh` | `DOUBLE_FRAC` default (0.224) |
| `efficiencytool/run_double_auto.sh` | Computed dynamically from `calc_pileup_range.C` |
| `efficiencytool/oneforall_tree_double_auto.sh` | Computed dynamically from `calc_pileup_range.C` |
| `efficiencytool/run_showershape_double_auto.sh` | Computed dynamically from `calc_pileup_range.C` (still dispatches to legacy script) |
| `efficiencytool/compare_efficiency_deltaR.C` | Hardcoded `DOUBLE_FRAC = 0.224f` |
| `efficiencytool/ShowerShapeCheck.C` | Receives via `mix_weight` argument |
| `efficiencytool/RecoEffCalculator_TTreeReader.C` | Receives via `mix_weight` argument |

## See Also

- [Double-Interaction Efficiency](../concepts/double-interaction-efficiency.md) -- Physics study of pileup effects on efficiency
- [Constants Sync](constants-sync.md) -- Other constants that must stay in sync across files
