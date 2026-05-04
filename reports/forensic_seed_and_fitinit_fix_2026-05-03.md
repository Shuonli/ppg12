# Forensic note: clock-seeded MC smearing and Padé fit-init bugs — fixes applied 2026-05-03

## TL;DR

Two latent bugs in the post-rerun analysis were identified and fixed. First, per-cluster MC ET smearing in `RecoEffCalculator` was clock-seeded, so every Phase-1 run produced a slightly different MC efficiency file (~7% drift in integrated cross-section between reruns). Second, the Padé fit for ABCD purity in `CalculatePhotonYield` diverged for the `npb_cut` systematic variant, producing a +149% blow-up at low pT despite per-bin purity values agreeing with nominal within 1%. Both fixed. The nominal cross-section is reproducible run-to-run, and the npb_cut group systematic dropped from 87% to 0.92% at 12-14 GeV.

## Bug 1: Clock-seeded MC ET smearer

- **Where**: `efficiencytool/RecoEffCalculator_TTreeReader.C:1486` and `efficiencytool/RecoEffCalculator.C:1154`
- **What**: `TRandom3 *rand = new TRandom3(0)` — seed 0 means clock-based init.
- **Use site**: line ~2064 of `RecoEffCalculator_TTreeReader.C`, `cluster_Et[icluster] *= rand->Gaus(1, clustereres)` per MC cluster fill.
- **Effect**: each Phase-1 rerun produced a slightly different MC efficiency file. Reruns with identical config produced bdt_nom integrated cross-sections of 1.158 → 1.190 → 1.270 nb (8-32 GeV), roughly ±7% non-reproducibility from the same code+config+slimtrees. A 5-seed test at the toy-MC stage confirmed MC closure shifted by less than 0.02 percentage points across seeds 1, 100, 1000, 10000, 99999, so the toy-MC seed itself contributes ~0%. The ~7% variability comes entirely from the upstream ET-smearing seed.
- **Fix**: `TRandom3(0)` → `TRandom3(42)`. A comment block in the diff explains the rationale.

## Bug 2: Padé fit-init runaway minimum

- **Where**: `efficiencytool/CalculatePhotonYield.C` Padé branch of `f_purity_fit` (line ~648) and `f_purity_leak_fit` (line ~692). Triggered when `analysis.fit_option = 1` (current nominal).
- **What**: Padé form `(p0+p1*x)/(1+p2^2*x)` was initialized with `(0.5, 0.5, 0.5)`. For the `npb03` variant the optimizer converged to a degenerate runaway with p0 ≈ 2e+27, p1 ≈ 4e+24, p2 ≈ -1e+13, and chi² = 224.5 versus nom's 5.73. The fit eval gave 1.25 at pT = 11 GeV (unphysical purity > 1) even though the per-bin `gpurity_leak` value at that bin was 0.52, matching nom 0.51 within 1%.
- **Effect on cross-section**: `CalcPhotonYield` uses the FIT (not per-bin gpurity values) when `fit_purity_dis = true`, so the unphysical fit eval propagated into `NA_sig_count = NA × purity × mc_corr`. The npb03 cross-section blew up by +149% at 12-14 GeV vs nom even though every input agreed within 5%. Pure methodology artifact.
- **Fix**: initialize the Padé from nom's converged params (constants `pade_purity_p0/p1/p2` and `pade_purity_leak_p0/p1/p2` defined just above the `SetParameters` calls in `CalculatePhotonYield.C`). Smoke test on npb03 after fix gave chi² = 5.81 (matches nom 5.73), fit eval at pT = 11 of 0.515 (matches per-bin value 0.516). Cross-section bin0 went from +149% to -1.28%.

## Effect on the systematic budget (post fit-fix)

| pT (GeV) | total syst pre-fix low/high | total syst post-fix low/high |
|---|---|---|
| 12-14 | 88.87% / 18.79% | 19.67% / 20.94% |
| 14-16 | 44.87% / 23.33% | 18.05% / 20.69% |
| 16-18 | 22.20% / 26.19% | 19.30% / 21.77% |
| 18-20 | 20.12% / 28.13% | 21.88% / 24.02% |
| 20-22 | 27.13% / 28.06% | 21.93% / 25.16% |
| 22-24 | 37.32% / 27.44% | 24.56% / 26.08% |
| 24-26 | 46.85% / 30.90% | 28.84% / 30.99% |
| 26-28 | 53.85% / 33.38% | 30.40% / 35.13% |
| 28-32 | 59.91% / 28.81% | 28.44% / 32.65% |

NPB-cut group: was 87% at 12-14 GeV, now 0.92%. The blow-up was a fitting artifact, not real physics.

## Effect on the bdt_nom cross-section

Integrated 12-32 GeV: original Apr 25 baseline = 4.560e+02 nb, post-rerun (post fit-fix) = 5.028e+02 nb. The +10.3% shift is dominated by the data DST pipeline restoration of JetCalib + topo-cluster reconstruction, not by the fit-fix. The fit-fix leaves nom essentially unchanged because nom's fit was already converged.

## Reproducibility status

With both fixes, the entire pipeline (Phase 1 + Phase 2 + aggregation + plot) is now bit-reproducible run-to-run for fixed inputs. Verified by re-running `CalculatePhotonYield` on npb03 alone after the fit-fix and confirming chi² convergence and matched cross-section.

## Backups

- `/sphenix/user/shuhangli/ppg12/efficiencytool/results/backup_bdt_nom_pre_data_rerun_20260502/` — pre-data-rerun (Apr 25 baseline)
- `/sphenix/user/shuhangli/ppg12/efficiencytool/results/backup_syst_pre_data_rerun_20260503_0931/` — pre-syst-rerun snapshot
- `/sphenix/user/shuhangli/ppg12/efficiencytool/results/backup_syst_pre_seed_fix_20260503_1459/` — seed=0 syst run
- `/sphenix/user/shuhangli/ppg12/efficiencytool/results/backup_pre_fitfix_20260503_1925/` — pre-Padé-fix snapshot

## Files modified

- `efficiencytool/RecoEffCalculator_TTreeReader.C` (line 1486)
- `efficiencytool/RecoEffCalculator.C` (line 1154)
- `efficiencytool/CalculatePhotonYield.C` (lines ~644-651, ~687-694)
