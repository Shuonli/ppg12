# PPG12 — Isolated Photon Cross-Section in pp at sqrt(s) = 200 GeV, sPHENIX

## Wiki

For detailed pipeline documentation, see `wiki/index.md`. Consult it before exploring unfamiliar parts of the codebase. The wiki covers pipeline stages, physics concepts, config schema, constants sync issues, and step-by-step guides.

## Pipeline

1. **anatreemaker/source/CaloAna24.cc** — DST files → slimtree ROOT files (tree: `slimtree`, node: `CLUSTERINFO_CEMC`)
2. **FunWithxgboost/BDTinput.C** — slimtree → training text files (`NOSPLIT/shapes_*.txt`)
3. **FunWithxgboost/main_training.py** — XGBoost training → models in `binned_models/`
4. **FunWithxgboost/apply_BDT.C** — scores clusters, writes BDT branches to slimtree
5. **efficiencytool/oneforall.sh** — runs `MergeSim.C` then `CalculatePhotonYield.C` (ABCD + unfolding + cross-section)
6. **plotting/** — final plots, systematic aggregation (`calc_syst_bdt.py`), cross-section vs JETPHOX NLO

## Key Commands

```bash
# BDT training
cd FunWithxgboost && python main_training.py --config config.yaml

# Apply BDT to data/sim
root -l -b -q 'apply_BDT.C("config_nom.yaml", "photon5")'

# Run efficiency + yield for one config (per-period dispatches to MergeSim;
# *_all.yaml dispatches to merge_periods.sh — see wiki/pipeline/full-run-range.md)
cd efficiencytool && bash oneforall.sh config_bdt_nom.yaml

# Lumi-weighted cross-period MC merge (TFileMerger of pre-scaled per-period
# outputs into combined all-range MC). Per-event lumi_weight = lumi/lumi_target
# in RecoEffCalculator does the scaling upstream; this is now a plain hadd.
# Naming convention: config_bdt_X.yaml (all-range) + config_bdt_X_{0rad,1p5mrad}.yaml (merge-feeders)
bash merge_periods.sh config_bdt_nom_0rad.yaml config_bdt_nom_1p5mrad.yaml config_bdt_nom.yaml

# Generate 23+ systematic variation configs
python make_bdt_variations.py config_bdt_nom.yaml

# Run all variations in parallel on HTCondor
condor_submit oneforall.sub

# Aggregate systematics and plot
cd plotting && python calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures

# Final cross-section plot
root -l -b -q 'plot_final.C'

# Compile analysis note
cd PPG12-analysis-note && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

## Physics Parameters

- **Eta**: `[-0.7, 0.7]` (barrel EMCal)
- **pT bins (reco)**: `[8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]` GeV (12 bins, from `plotcommon.h`)
- **pT bins (truth)**: extended on both sides for unfolding overflow
- **Reco isolation**: topo-cluster R=0.4 (nominal, `use_topo_iso: 2`), parametric cut `reco_iso_max = reco_iso_max_b + reco_iso_max_s * ET`
- **Truth isolation**: cone R=0.3, `iso_ET_truth < 4 GeV` (fiducial definition for cross-section)
- **BDT threshold**: parametric `tight_bdt_min = intercept + slope * ET` (ET-dependent, not flat)
- **Luminosity**: 16.6 pb^-1 (Run 24)
- **Cluster node**: `CLUSTERINFO_CEMC` (current), `CLUSTERINFO_CEMC_NO_SPLIT` (legacy)
- **Global seed**: 42

## Configuration Files

- **FunWithxgboost/config.yaml** — training: 25 features, XGB hyperparams, ET binning, reweighting
- **efficiencytool/config_bdt_*.yaml** — analysis: cuts, paths, `var_type` suffix
- `var_type` is appended to output filenames (`Photon_final_{var_type}.root`) — must be unique per variation

## MC Samples

- Signal: photon5 (0-14 GeV), photon10 (14-30), photon20 (30+)
- Background: jet10/15/20/30/50
- Cross-section weights normalized per sample; constants must match across MergeSim.C, RecoEffCalculator.C, CalculatePhotonYield.C

## ABCD Background Subtraction

```
            Isolated     Non-isolated
Tight       A (signal)   B (bkg template)
Non-tight   C (bkg)      D (bkg control)

Background in A = R * (B * C) / D  (with signal leakage corrections cB, cC, cD)
```

## Plotting Conventions

- `#include "plotcommon.h"` — shared pT bins (NptBins=12), frames, legend strings
- `SetsPhenixStyle()` — sPHENIX publication style
- `BlairUtils.C` — error-band construction, TGraph operations
- Legend: `"#bf{#it{sPHENIX}} Internal"`, `"#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, 16.6 pb^{-1}"`
- Figures saved to `plotting/figures/` as PDF

## Systematics Workflow

See `SYST_PIPELINE_README.md` for full details. Three steps:
1. `make_bdt_variations.py` generates 23+ configs from nominal
2. `condor_submit oneforall.sub` runs efficiency tool on all configs
3. `calc_syst_bdt.py` aggregates into groups (purity, eff, escale, eres) via quadrature + luminosity

## Custom Agents

- `physics-reviewer` — read-only physics code review
- `data-explorer` — explore ROOT trees, configs, result consistency
- `report-writer` — write/edit LaTeX analysis note and reports
- `code-writer` — write new analysis macros following repo conventions

## Analysis Note

Located at `PPG12-analysis-note/` — document class `sphenix`, tag `sPH-ppg-2024-012`.
Custom macros in `defs.sty`: `\etg`, `\pp`, `\isoET`, `\purity`, `\eff`, `\lumi`, `\wetacogx`, etc.
Several plotting macros save figures directly to `PPG12-analysis-note/Figures/`.

## Data Paths

- Sim input: `/sphenix/user/shuhangli/ppg12/FunWithxgboost/` (per-sample BDT-scored trees)
- Data input: `/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/condorout/part_*_with_bdt_split.root`
- Results: `/sphenix/user/shuhangli/ppg12/efficiencytool/results/`
- BDT models: `/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/`

## GitHub Pages

Analysis results published at `https://shuonli.github.io/ppg12/` (enable in repo Settings → Pages → gh-pages branch).

```bash
# Deploy reports + selected figures to the site
python3 scripts/deploy_pages.py

# Preview without pushing
python3 scripts/deploy_pages.py --dry-run
```

- Config: `scripts/site_config.yaml` — controls which reports and figures are published
- Reports in `reports/` (`.pdf`, `.md`) are auto-discovered
- Figures require explicit allowlist in `site_config.yaml` under `figures.include`
- Site content lives on the orphan `gh-pages` branch (managed via `_gh-pages-worktree/`)
- Skill: `/deploy-pages` wraps the deploy script

## References

- `README.md` — full pipeline overview and physics motivation
- `SYST_PIPELINE_README.md` — systematic variation workflow (3-step pipeline)
- `PLOTTING_UPDATE_SUMMARY.md` — recent plotting code changes
