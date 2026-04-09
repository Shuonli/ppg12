---
description: Compile the PPG12 analysis note (pdflatex + bibtex)
---

Compile the PPG12 analysis note.

Arguments: $ARGUMENTS
- (empty) → full build: pdflatex → bibtex → pdflatex x2
- `quick` → single pdflatex pass (fast iteration, no bibtex)
- `clean` → remove aux files first, then full build

Working directory: `/sphenix/user/shuhangli/ppg12/PPG12-analysis-note`

## Steps

1. If `$ARGUMENTS` contains `clean`:
   ```
   rm -f main.aux main.bbl main.blg main.log main.out main.toc main.lof main.lot
   ```

2. If `$ARGUMENTS` contains `quick`:
   ```
   pdflatex main.tex
   ```
   Otherwise (full build):
   ```
   pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
   ```

3. After compilation:
   - Check if `main.pdf` was produced and report its size
   - Grep `main.log` for warnings about undefined references or missing figures
   - Report any LaTeX errors (lines starting with `!`)
