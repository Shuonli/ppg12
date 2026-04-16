---
description: Deploy updated reports and figures to GitHub Pages
---

Deploy the PPG12 analysis results site to GitHub Pages.

Arguments: $ARGUMENTS (optional flags passed to deploy_pages.py)

## Common flags

- `--dry-run` — show what would be deployed without pushing
- `--message "msg"` — custom commit message (default now includes time-of-day: `Update site: YYYY-MM-DD HH:MM TZ`)
- `--tag NAME` — after a successful deploy, tag the gh-pages commit as `site/NAME` and main HEAD as `src/NAME`; push both. Use for milestones (e.g. `--tag pwg-2026-03-15`)
- `--no-tag-main` — when `--tag` is given, skip the `src/NAME` tag on main

## Time-travel / past snapshots

Every deploy is already a full snapshot on `gh-pages`. Helpers:

- `--snapshot-at REF` — materialize a past gh-pages state in `_gh-pages-snapshot/` (a second worktree). `REF` can be a commit hash, a tag (`site/2026-04-15`), or a date (`2026-04-15` → last deploy on or before that day)
- `--close-snapshot` — remove the snapshot worktree
- `--backfill-tags` — create `site/YYYY-MM-DD[-N]` tags for every existing gh-pages commit and push to origin. One-shot.

Examples:
```
# Deploy + name this state
python3 scripts/deploy_pages.py --tag DAC-review-2026-Q2

# Look at the site as it was two weeks ago
python3 scripts/deploy_pages.py --snapshot-at 2026-04-02
# ...browse _gh-pages-snapshot/index.html ...
python3 scripts/deploy_pages.py --close-snapshot

# Grab one file from a past deploy without a worktree
git show site/2026-04-13:reports/showershape_report.pdf > /tmp/shower_old.pdf
```

## Steps

1. Run the deploy script from the repo root:
   ```
   cd /sphenix/user/shuhangli/ppg12 && python3 scripts/deploy_pages.py $ARGUMENTS
   ```

2. The script auto-compiles any `.tex` reports whose PDF is missing or stale (runs `pdflatex` twice). Check the output for compilation warnings.

3. The deploy now diff-copies into the worktree (only changed files touch disk). Per-directory summary is printed (`reports: 2 copied, 0 removed (of 17)`).

4. Report which files were compiled, deployed, the total size, the site URL, and any tags created.

5. If this is the first deployment, remind the user to enable GitHub Pages:
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
