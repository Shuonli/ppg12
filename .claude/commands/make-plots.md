---
description: Generate analysis plots (selection, systematics, final cross-section)
---

Generate analysis plots.

Arguments: $ARGUMENTS
- `selection <suffix>` → run `make_selection_plots.sh` for one variant (e.g., `selection bdt_nom`)
- `selection-all` → run `make_all_bdt_selection.sh` for ALL variants
- `syst-bdt` → run `calc_syst_bdt.py` (Python-based systematics aggregation + plots)
- `final` → run `plot_final.C`
- `all <suffix>` → selection + syst-bdt + final

Working directory: `/sphenix/user/shuhangli/ppg12/plotting`

## Steps

1. Parse the first word of `$ARGUMENTS` to determine plot type.

2. Source environment and run the appropriate command:
   - `selection <suffix>`: `bash make_selection_plots.sh <suffix>`
   - `selection-all`: `bash make_all_bdt_selection.sh`
   - `syst-bdt`: `python calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures`
   - `final`: `root -l -b -q 'plot_final.C'`
   - `all <suffix>`: run selection, then syst-bdt, then final sequentially

3. After running, list new/updated PDF files in `figures/` from the last 5 minutes.

4. Report any ROOT errors found in the terminal output.
