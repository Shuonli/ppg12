---
description: Deploy updated reports and figures to GitHub Pages
---

Deploy the PPG12 analysis results site to GitHub Pages.

Arguments: $ARGUMENTS (optional flags passed to deploy_pages.py)
- `--dry-run` — show what would be deployed without pushing
- `--message "msg"` — custom commit message

## Steps

1. Run the deploy script from the repo root:
   ```
   cd /sphenix/user/shuhangli/ppg12 && python3 scripts/deploy_pages.py $ARGUMENTS
   ```

2. The script auto-compiles any `.tex` reports whose PDF is missing or stale (runs `pdflatex` twice). Check the output for compilation warnings.

3. Report which files were compiled, deployed, the total size, and the site URL.

3. If this is the first deployment, remind the user to enable GitHub Pages:
   - Go to https://github.com/Shuonli/ppg12/settings/pages
   - Source: Deploy from a branch → `gh-pages` / `/ (root)`
   - Save, wait 1-2 minutes

## Adding figures to the site

Edit `scripts/site_config.yaml` and add filenames to `figures.include`:
```yaml
figures:
  include:
    - "final_all_bdt_nom.pdf"
    - "syst_bdt_breakdown.pdf"
```
Then re-run `/deploy-pages`.
