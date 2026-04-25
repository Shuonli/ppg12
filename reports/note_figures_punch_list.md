# Analysis-Note Figure Regeneration — Macro-by-Macro Punch List

**Date:** 2026-04-23
**Decisions locked:** default tune = `bdt_nom`; for §5, update note `\includegraphics` + prose (option b).

## Summary table

| # | Macro | Figures (in note) | Stale? | Change required |
|---|-------|-------------------|--------|-----------------|
| 1 | `plotting/plot_final_selection.C` | 11 `Figures/results/final*.pdf` | YES | Macro OK. Note: 11 `\includegraphics` renames. |
| 2 | `plotting/calc_syst_bdt.py` | 18 `Figures/systematic/syst_*.pdf` | YES | Macro OK (taxonomy current). Note: rewrite §5 for new taxonomy. |
| 3 | `plotting/plot_efficiency.C` | 6 `Figures/analysis/eff_*.pdf` | YES | Macro: update axis range for new pT. Note: 6 renames. |
| 4 | `plotting/plot_xjg.C` | `eff_iso.pdf` (duplicate of #3), xjg_*.pdf | YES | Macro: input filename change (3 lines). |
| 5 | `plotting/plot_purity_selection.C` | `Figures/analysis/purity_nom.pdf` | YES | Macro OK. Note: 1 rename. |
| 6 | `plotting/plot_purity_sim_selection.C` | `Figures/analysis/purity_sim_nom.pdf` | YES | Macro OK. Note: 1 rename. |
| 7 | `plotting/plot_sideband.C` | `et_sbs.pdf`, `et_sbs_ratio.pdf`, `leakage_fraction_et.pdf`, `SB.pdf`, `sideband_diagram.png`, `iso_ET_pt*.pdf` | LIKELY | To inspect — macro-side changes TBD, same naming pattern expected |
| 8 | `plotting/plot_reweight.C` | `response_reweight.pdf` | LIKELY | To inspect |
| 9 | `plotting/plot_response.C` | `response.pdf` | LIKELY | To inspect |
| 10 | `plotting/plot_unfold_iter.C` | `unfold_iter.pdf`, `closure_full.pdf`, `closure_half.pdf` | LIKELY | To inspect |
| 11 | (several scripts, §3) | ~25 `Figures/reconstruction/*.pdf` | YES | Deferred — §3 is large, lower priority |
| 12 | (trigger macro, §2) | `Figures/etc/*.pdf` | YES | Task #7 (Phase 3a) |
| 13 | External / manual | `Figures/introduction/*.png` | NO | Feynman diagrams, external images |

---

## Details per macro

### 1. `plotting/plot_final_selection.C` — §6 final cross-section

**Inputs (current — today 15:07):**
- `efficiencytool/results/Photon_final_bdt_nom.root`
- `efficiencytool/results/Photon_final_bdt_nom_mc.root`
- lumi read dynamically from `efficiencytool/config_bdt_nom.yaml` (→ 64.3718)

**Output:** `plotting/figures/final_bdt_nom.pdf`, `final_phenix_bdt_nom.pdf`, `final_all_bdt_nom.pdf`, `final_common_cluster_bdt_nom.pdf`, `final_tight_iso_cluster_bdt_nom.pdf`

**Macro changes:** NONE. Stale `datalumi = 49.562` fallback is overridden by YAML read, safe.

**Note changes (results.tex):** rename 11 `\includegraphics`:
| Old path (stale) | New path (from bdt_nom regen) |
|---|---|
| `Figures/results/final.pdf` | `Figures/results/final_bdt_nom.pdf` |
| `Figures/results/final_phenix.pdf` | `Figures/results/final_phenix_bdt_nom.pdf` |
| `Figures/results/final_all.pdf` | `Figures/results/final_all_bdt_nom.pdf` |
| `Figures/results/final_common_cluster.pdf` | `Figures/results/final_common_cluster_bdt_nom.pdf` |
| `Figures/results/final_tight_iso_cluster.pdf` | `Figures/results/final_tight_iso_cluster_bdt_nom.pdf` |
| `Figures/results/final_JETPHOX.pdf`, `final_JETPHOX_rewe.pdf` | likely dropped or replaced by `final_phenix_bdt_nom.pdf`; to confirm vs current prose |
| `Figures/results/final_Pythia.pdf`, `final_Pythia_rewe.pdf` | same as above |
| `Figures/results/final_sphenix.pdf`, `final_sphenix_.pdf` | same |

**Regen command:**
```bash
cd plotting && root -l -b -q 'plot_final_selection.C("bdt_nom")'
cp figures/final_bdt_nom.pdf ../PPG12-analysis-note/Figures/results/
cp figures/final_phenix_bdt_nom.pdf ../PPG12-analysis-note/Figures/results/
cp figures/final_all_bdt_nom.pdf ../PPG12-analysis-note/Figures/results/
cp figures/final_common_cluster_bdt_nom.pdf ../PPG12-analysis-note/Figures/results/
cp figures/final_tight_iso_cluster_bdt_nom.pdf ../PPG12-analysis-note/Figures/results/
```

---

### 2. `plotting/calc_syst_bdt.py` — §5 systematics (TAXONOMY CHANGE)

**Inputs:** `efficiencytool/results/Photon_final_<var_type>.root` for every systematic variation defined in `make_bdt_variations.py:VARIANTS`.

**Output:** `plotting/figures/syst_bdt_rel_<type>.pdf` per category + `syst_bdt_breakdown.pdf` (total).

**Macro changes:** NONE for the plotting logic. But must verify `LUMI_SYST` value in `make_bdt_variations.py` matches the allz-setup uncertainty (up / down percentages) — spot-check needed.

**Note changes — §5 systematics.tex — MAJOR:**

The old 9-rel + 8-spec + 1-sum set (18 figures) maps onto a new 13-channel taxonomy that doesn't line up 1-to-1. Mapping table:

| Old filename (note) | New filename (pipeline) | Semantic meaning |
|---|---|---|
| `syst_rel_tight.pdf` | `syst_bdt_rel_tight_bdt.pdf` | tight-BDT-threshold variation |
| `syst_rel_ntf.pdf` | `syst_bdt_rel_nt_bdt.pdf` | non-tight-BDT-range |
| `syst_rel_nt.pdf` | part of above | (merged into `nt_bdt`) |
| `syst_rel_iso.pdf` | `syst_bdt_rel_noniso.pdf` | iso-sideband shift |
| `syst_rel_fit.pdf` | `syst_bdt_rel_purity_fit.pdf` | purity fit form |
| `syst_rel_nor.pdf` | `syst_bdt_rel_reweight.pdf` | unfolding reweight |
| `syst_rel_escale.pdf` | `syst_bdt_rel_escale.pdf` | same name |
| `syst_rel_eres.pdf` | `syst_bdt_rel_eres.pdf` | same name |
| `syst_rel_mbd.pdf` | `syst_bdt_rel_timing.pdf`? | MBD vs timing — check mapping |
| `syst_sum_rel.pdf` | `syst_bdt_breakdown.pdf` | total |
| (NEW — not in note) | `syst_bdt_rel_mc_purity_correction.pdf` | MC purity closure correction — maps to new §5 subsection from Phase 2 |
| (NEW) | `syst_bdt_rel_vtx_reweight.pdf` | truth-vertex reweight closure — Phase 2 new subsection |
| (NEW) | `syst_bdt_rel_bdt_model.pdf` | BDT-model variant |
| (NEW) | `syst_bdt_rel_b2bjet.pdf` | back-to-back jet |
| (NEW) | `syst_bdt_rel_npb_cut.pdf` | NPB cut variation |
| `syst_spec_*.pdf` (8) | (not produced by new pipeline) | spectrum-level plots dropped |

**Note changes:**
1. Replace 9 `syst_rel_*` `\includegraphics` in `systematics.tex` with `syst_bdt_rel_*`.
2. Remove 8 `syst_spec_*` `\includegraphics` (new pipeline doesn't produce them).
3. Replace `syst_sum_rel.pdf` → `syst_bdt_breakdown.pdf`.
4. Add 5 new `\includegraphics` for the new channels (mc_purity_correction, vtx_reweight, bdt_model, b2bjet, npb_cut).
5. Rewrite §5 prose: discuss 13 channels + lumi instead of the old 9 rel + 1 sum. The Phase 2 §5 subsections (mc_purity_closure, tower_acc, l1_plateau, truth_vtx) already partially aligned with the new taxonomy — revise to flow naturally.

**Regen command:**
```bash
cd plotting && python3 calc_syst_bdt.py \
  --results ../efficiencytool/results \
  --outdir rootFiles \
  --figdir figures
cp figures/syst_bdt_rel_*.pdf figures/syst_bdt_breakdown.pdf ../PPG12-analysis-note/Figures/systematic/
```

---

### 3. `plotting/plot_efficiency.C` — §4 efficiency (6 figures)

**Inputs (current, correct naming):**
- `efficiencytool/results/MC_efficiency_bdt_nom.root`
- `efficiencytool/results/Photon_final_bdt_nom.root`

**Output:** `eff_reco_bdt_nom.pdf`, `eff_iso_bdt_nom.pdf`, `eff_id_bdt_nom.pdf`, `eff_total_bdt_nom.pdf`, `eff_mbd_bdt_nom.pdf`, `eff_mbd_frac_bdt_nom.pdf`, `eff_photon_bdt_nom.pdf`

**Macro changes:**
- Line ~29: `frame_et_truth->GetXaxis()->SetRangeUser(10, 26)` → **`10, 32`** (new reporting range per Phase 1 binning decision)
- Similar `(10, 26)` occurrences in other canvases — sweep the file and widen to `(10, 32)`

**Note changes (analysis.tex):** rename 6 `\includegraphics`:
| Old | New |
|---|---|
| `Figures/analysis/eff_reco.pdf` | `Figures/analysis/eff_reco_bdt_nom.pdf` |
| `Figures/analysis/eff_id.pdf` | `Figures/analysis/eff_id_bdt_nom.pdf` |
| `Figures/analysis/eff_iso.pdf` | `Figures/analysis/eff_iso_bdt_nom.pdf` |
| `Figures/analysis/eff_photon.pdf` | `Figures/analysis/eff_photon_bdt_nom.pdf` |
| `Figures/analysis/eff_mbd.pdf` | `Figures/analysis/eff_mbd_bdt_nom.pdf` |
| `Figures/analysis/eff_total.pdf` | `Figures/analysis/eff_total_bdt_nom.pdf` (if present in note — check) |

**Regen:**
```bash
cd plotting && root -l -b -q 'plot_efficiency.C("bdt_nom")'
cp figures/eff_*_bdt_nom.pdf ../PPG12-analysis-note/Figures/analysis/
```

---

### 4. `plotting/plot_xjg.C` — §4 x_Jγ (mostly diagnostic)

**Inputs (currently STALE names hardcoded):**
- `MC_efficiency_nom.root` (should be `MC_efficiency_bdt_nom.root`)
- `MC_efficiency_jet_nom.root` (should be `MC_efficiency_jet_bdt_nom.root`)
- `data_histo_nom.root` (should be `Photon_final_bdt_nom.root` — the data-histogram structure moved)

**Output:** `xjg_signalmc.pdf`, `xjg_incmc.pdf`, `xjg_data.pdf`, `xjg_data_signalmc_incmc.pdf`, `eff_iso.pdf` (DUPLICATE — already produced by `plot_efficiency.C`)

**Macro changes (3–4 lines at plot_xjg.C:17–19):**
```c++
// Before:
TFile *fmc = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_nom.root");
TFile *fmc_inclusive = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_nom.root");
TFile *fdata = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_nom.root");

// After:
TFile *fmc = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root");
TFile *fmc_inclusive = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_bdt_nom.root");
TFile *fdata = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root");
```

**Note changes:** `eff_iso.pdf` is the only note-referenced output and already covered by `plot_efficiency.C` (#3). If `plot_xjg.C` is rerun, its `eff_iso.pdf` should NOT overwrite the `plot_efficiency.C` one — simplest to comment out line 142 (`SaveAs eff_iso.pdf`) so there's no conflict, or rename its output.

---

### 5. `plotting/plot_purity_selection.C` — §4 data purity

**Inputs (current, correct):** `Photon_final_bdt_nom.root` + `Photon_final_bdt_nom_mc.root`

**Output:** `plotting/figures/purity_bdt_nom.pdf`

**Macro changes:** NONE.

**Note change (analysis.tex):**
- `\includegraphics{Figures/analysis/purity_nom.pdf}` → `\includegraphics{Figures/analysis/purity_bdt_nom.pdf}`

**Regen:**
```bash
cd plotting && root -l -b -q 'plot_purity_selection.C("bdt_nom")'
cp figures/purity_bdt_nom.pdf ../PPG12-analysis-note/Figures/analysis/
```

---

### 6. `plotting/plot_purity_sim_selection.C` — §4 simulation purity

**Inputs (current, correct):** `Photon_final_bdt_nom_mc.root`

**Output:** `plotting/figures/purity_sim_bdt_nom.pdf`

**Macro changes:** NONE.

**Note change (analysis.tex):**
- `\includegraphics{Figures/analysis/purity_sim_nom.pdf}` → `\includegraphics{Figures/analysis/purity_sim_bdt_nom.pdf}`

**Regen:**
```bash
cd plotting && root -l -b -q 'plot_purity_sim_selection.C("bdt_nom")'
cp figures/purity_sim_bdt_nom.pdf ../PPG12-analysis-note/Figures/analysis/
```

---

### 7–10. To-inspect (Priority-3 deeper dive)

| Macro | Note figures | Expected changes |
|-------|--------------|------------------|
| `plotting/plot_sideband.C` | `et_sbs.pdf`, `et_sbs_ratio.pdf`, `leakage_fraction_et.pdf`, `SB.pdf`, `iso_ET_pt*.pdf`, `iso_ET_sim_pt*.pdf`, `iso_ET_pt*_comb.pdf` | Likely: input Photon_final_bdt_nom.root (correct naming); output gets `_bdt_nom` suffix. Note needs renames. |
| `plotting/plot_reweight.C` | `response_reweight.pdf` | Same pattern |
| `plotting/plot_response.C` | `response.pdf` | Same pattern |
| `plotting/plot_unfold_iter.C` | `unfold_iter.pdf`, `closure_full.pdf`, `closure_half.pdf` | Same pattern |

**Policy for these four:** inspect the macro header, apply the same "change input suffix → `bdt_nom` if needed; rename note `\includegraphics`" recipe.

---

### 11. §3 reconstruction figures (~25 figures, deferred)

Produced by a mix of:
- `efficiencytool/ShowerShapeCheck.C` (direct ROOT output)
- `efficiencytool/FindETCut.C` (direct)
- plotting macros like `plot_showershapes_selections.C`, `plot_iso_et_cut.C` (TBD)

**Recommendation:** handle as a Phase 4 deeper sweep once Priority-1/2/3 land. Lower impact because Phase 2 already inserted the main §3 methodology updates (topo-iso R=0.4 description, truth-vertex reweight pointer); the figures are supporting material.

---

### 12. §2 trigger efficiency figure — already queued

Task #7 (Phase 3a): regenerate `Figures/etc/Photon_4_GeV_ConstantFit.pdf` + maxClusterEnergy overlay against bit-30 selection. Covered separately.

### 13. §1 Feynman diagrams — no regen

`Figures/introduction/compton_scattering.png`, `annihilation.png`, `fragmentation.png`, `pT_Spectrum_IsoET_4GeV_Cut.pdf` — all exist on disk, LOW-priority per earlier plan. No script regens them.

---

## Proposed execution plan (1–2 hours total)

**Priority-1 (25 min, no macro changes):**
1. Run `plot_final_selection.C` + copy 5 PDFs (§6).
2. Run `plot_purity_selection.C` + `plot_purity_sim_selection.C` + copy 2 PDFs (§4).
3. Update `results.tex` + `analysis.tex` `\includegraphics` calls to new `_bdt_nom` filenames. Total rename edits: ~13 in `results.tex`, ~2 in `analysis.tex`.

**Priority-2 (30 min, one macro change + note prose rewrite):**
1. Widen `plot_efficiency.C` X-axis range 10–26 → 10–32 (5 occurrences expected).
2. Run `plot_efficiency.C` + copy 6 PDFs.
3. Update `analysis.tex` `\includegraphics` (6 renames).
4. Run `calc_syst_bdt.py` + copy `syst_bdt_rel_*` + `syst_bdt_breakdown.pdf` to `Figures/systematic/`.
5. Rewrite `systematics.tex` §5 for the new 13-channel taxonomy (retain the Phase 2 subsections since they already align).

**Priority-3 (30 min, inspect-then-regen):**
1. Open `plot_sideband.C`, `plot_reweight.C`, `plot_response.C`, `plot_unfold_iter.C` — confirm each reads `bdt_nom`; run each; copy outputs; rename note refs.

**Priority-4 (deferred):** §3 reconstruction figures.

---

## Risks

1. **§5 taxonomy rewrite is nontrivial** — Phase 2 already wrote 4 new subsections (mc_purity_closure, tower_acc, l1_plateau, truth_vtx) that partially anticipate the new taxonomy. Need to reconcile those with the `calc_syst_bdt.py` output names (e.g., tower_acc subsection → no direct `syst_bdt_rel_tower_acc.pdf` exists; this systematic was added to the note but isn't in the automated pipeline).
2. **MBD systematic mapping ambiguity** — old `syst_rel_mbd.pdf` vs new `syst_bdt_rel_timing.pdf` may not be the same thing. Verify against `make_bdt_variations.py:SYST_TYPES` before renaming.
3. **Missing `syst_spec_*` family** — the new pipeline doesn't produce spectrum-level syst plots. Note §5 currently shows 8 of these; removing them requires prose that references them to be removed/rewritten.
4. **`plot_xjg.C` produces `eff_iso.pdf`** with no suffix, conflicting with `plot_efficiency.C`'s `eff_iso_bdt_nom.pdf`. Run order matters, or one save should be commented out.
