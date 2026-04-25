# L1 Photon-4-GeV Trigger (Bit 30) Turn-On vs Cluster ET

**Date**: 2026-04-21
**Author**: PPG12 ask-shuonli
**Question**: Using the data tree, for events where the MBD scaled bit is on, what fraction also have bit 30 (Photon-4-GeV) on, as a function of cluster ET starting at 7 GeV?

## TL;DR

- **MBD reference bit = 10** (`MBD_N&S_geq_1`). Live-on fraction ~99.93% per event; scaled-down prescale ~1660 (too aggressive for reference — used *live* decision).
- **Bit 30 = `Photon_4_GeV_plus_MBD_NS`**. Analysis trigger.
- **Integrated efficiency at ET ≥ 8 GeV = 0.9958 ± 0.0001** (per-cluster, 355,656 / 357,166 clusters).
- **Fine-binned turn-on near threshold**:
  - 7.0–7.5 GeV: 0.988
  - 7.5–8.0 GeV: 0.993
  - 8.0–8.5 GeV: 0.995
  - 8.5–9.0 GeV: 0.996
  - plateau 1.000 by ~27 GeV
- **Recommendation**: The PPG12 pipeline currently applies **no L1 plateau correction** (verified in `RecoEffCalculator_TTreeReader.C`, `CalculatePhotonYield.C`, `MergeSim.C`, `LumiCalculator.C`). The ~0.4% low bias is roughly flat across the analysis pT range and should either be (a) applied as a per-ET correction `1/ε_L1(ET)` or (b) quoted as a known uncorrected low-side bias. Not large relative to the ~3–10% purity systematics, but it is a systematic one-sided bias, not a statistical fluctuation.

## Measurement

**Sample**: 10 of 77 `chunk_aa/part_*_with_bdt_split.root` files, 6,619,944 events.
**Tree**: `slimtree`.
**Branches**: `scaledtrigger[64]`, `livetrigger[64]`, `ET_core_calib`, fiducial vars.

**Fills** (both num & den):
- Fiducial cut: `|cluster_eta| < 0.7`, `|vertexz| < 60 cm`.
- Per-cluster: one entry per cluster at its `ET_core_calib`.
- Cross-check: per leading-ET cluster per event.

**Denominator**: `livetrigger[10] == 1` (MBD NandS coincidence, ~100% live).
**Numerator**: `livetrigger[10] == 1 AND livetrigger[30] == 1`.
**Binning**: 0.5 GeV from 7 to 20, 1 GeV from 20 to 40.

**Why `livetrigger` not `scaledtrigger`**: the bit-10 scaled trigger has prescale ~1660 (≈0.08% scaled-on rate) — insufficient statistics. Prescale is a per-bit random decimation independent of cluster ET, so it cancels exactly in the ratio. Verified by the physics reviewer.

## Per-analysis-bin efficiency (per cluster)

| pT bin [GeV] | N_den | N_num | ε_L1(bit30\|MBD) | bias if ε=1 |
|---|---|---|---|---|
| 8–10 | 241,685 | 240,598 | 0.9955 ± 0.0001 | +0.45% |
| 10–12 | 72,328 | 72,063 | 0.9963 ± 0.0002 | +0.37% |
| 12–14 | 24,904 | 24,808 | 0.9961 ± 0.0004 | +0.39% |
| 14–16 | 9,704 | 9,668 | 0.9963 ± 0.0006 | +0.37% |
| 16–18 | 4,217 | 4,207 | 0.9976 ± 0.0007 | +0.24% |
| 18–20 | 2,008 | 2,002 | 0.9970 ± 0.0012 | +0.30% |

Higher-pT bins (20–36 GeV) each sit at ≥0.997 with stats-dominated error.

## Per-cluster vs per-leading-cluster

Agreement within 0.0001 everywhere. Cluster multiplicity > 1 at ET > 8 GeV is very rare in this dataset, so per-cluster fill does not bias vs the "intrinsic" leading-cluster turn-on.

## Caveats & follow-ups

1. **Conditional eff, not absolute**. Bit 30 name `Photon_4_GeV_plus_MBD_NS` and bit 10 `MBD_N&S_geq_1` both require the same MBD-coincidence primitive, so what is measured here is ε(photon-primitive | MBD-NS). This is the correct multiplicative correction for the analysis (which already cuts on bit 30, inheriting the MBD-NS requirement). The absolute MBD-NS efficiency is covered separately by the MBD efficiency study.

2. **Not split by eta**. The physics reviewer flagged an eta-dependence cross-check (L1 photon primitive is a 4×4 tower sum; edge/crack clusters may have a softer turn-on). Worth checking as a follow-up before quoting a single per-bin correction.

3. **Memory flag "bit 26 vs bit 30"**. Only the legacy `makeHist_showerShapes_v6.C:282` uses bit 26 (OR of 24/25/26/27). The active pipeline (`config_bdt_nom*.yaml`, `RecoEffCalculator_TTreeReader.C:1519–1530`) cuts on bit 30 — consistent with this measurement.

## Files

| Path | Contents |
|---|---|
| `efficiencytool/results/trigger_bit30_turnon.root` | `h_num_{percluster,leading}`, `h_den_{percluster,leading}`, `h_eff_{percluster,leading}_binomial` |
| `plotting/figures/trigger_bit30_turnon.pdf` | Turn-on plot, per-cluster + leading-cluster overlays |
| `reports/trigger_bit30_turnon.md` | This report |

## Reviewer pass matrix

| Artifact | 2a physics | 2b cosmetics | 2c re-derivation |
|---|---|---|---|
| Method + code path | PASS (Q1–Q5) | N/A | N/A |
| Turn-on plot PDF | N/A | PASS (minor nit: add lumi tag to legend) | N/A |
| Integrated & per-bin numbers | N/A | N/A | PASS (15/15 values match ROOT to <1e-4) |

## Recommendation (concrete)

Either:
- **Apply correction**: multiply the data-side yield by `1/ε_L1(ET)` per pT bin. Values above. Assign ~0.05–0.10% per-bin systematic on the fit parameter (statistical error on ε_L1 at current sample size). Implementation point: `RecoEffCalculator_TTreeReader.C` data loop where the trigger cut is applied (line 1519–1530).
- **Document only**: add to the analysis note "no L1 plateau correction applied; one-sided low bias ~0.4% roughly flat across ET ≥ 8 GeV" and leave as an open item. Matches existing memory-flagged concern from the cut-efficiency audit.
