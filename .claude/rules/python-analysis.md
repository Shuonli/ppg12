---
globs: "**/*.py"
---

# Python Analysis Conventions (PPG12)

- Training config: `FunWithxgboost/config.yaml` — 25 features under `data.features`, XGB hyperparams, ET binning
- Efficiency configs: `efficiencytool/config_bdt_*.yaml` — analysis cuts, paths, `var_type`
- 25 features in `config.yaml` must stay in sync with `apply_BDT.C` feature vectors — update both when adding/removing features
- `global_seed: 42` — use for all random state initialization (numpy, sklearn, xgboost)
- `ruamel.yaml` (not PyYAML) for config generation — preserves comments and field ordering (see `make_bdt_variations.py`)
- `uproot` for ROOT file I/O; tree name is always `slimtree`
- `SYST_TYPES` and `SYST_GROUPS` dicts in `make_bdt_variations.py` define the systematics taxonomy
- When adding a new systematic variant: add to `VARIANTS` dict, optionally add to `SYST_TYPES`/`SYST_GROUPS`, re-run the 3-step pipeline
- CLI pattern: `argparse` with `--config`, `--results`, `--outdir`, `--figdir` (see existing scripts)
- Report generation scripts (`make_showershape_report.py`, `make_selection_report.py`, `make_comparison_report.py`) discover configs, collect figures, generate LaTeX
