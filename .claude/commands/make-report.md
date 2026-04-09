---
description: Generate automated LaTeX reports (selection, showershape, or comparison)
---

Generate an automated LaTeX report from analysis output.

Arguments: $ARGUMENTS (first word is report type, rest are options)
- `selection [--compile]` → run `make_selection_report.py`
- `showershape [--suffix X] [--compile] [--mixed]` → run `make_showershape_report.py`
- `comparison [--compile]` → run `make_comparison_report.py`

Working directory: `/sphenix/user/shuhangli/ppg12/plotting`

## Steps

1. Parse report type from first word of `$ARGUMENTS`.

2. Source environment and run the corresponding Python script, passing remaining arguments directly:
   ```
   source /sphenix/u/shuhang98/setup.sh && python make_{type}_report.py {remaining_args}
   ```

3. Report which `.tex` files were generated and their locations.

4. If `--compile` was included, report whether pdflatex succeeded and the output PDF path.
