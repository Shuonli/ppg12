#!/usr/bin/env python3
"""Deploy PPG12 analysis results to GitHub Pages.

Collects reports and figures from the repository, generates an index.html,
and pushes everything to the gh-pages branch via a git worktree.

Usage:
    python3 scripts/deploy_pages.py [--dry-run] [--message MSG] [--config PATH]
"""

import argparse
import fnmatch
import html
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


# ---------------------------------------------------------------------------
# YAML parsing — prefer ruamel.yaml (repo dependency), fall back to PyYAML,
# then to a minimal built-in parser for the flat structure we need.
# ---------------------------------------------------------------------------

def _load_yaml(path: Path) -> dict:
    """Load a YAML file, trying available libraries in order."""
    text = path.read_text()

    # Try ruamel.yaml first (already used in make_bdt_variations.py)
    try:
        from ruamel.yaml import YAML
        yaml = YAML()
        from io import StringIO
        return yaml.load(StringIO(text))
    except ImportError:
        pass

    # Try PyYAML
    try:
        import yaml
        return yaml.safe_load(text)
    except ImportError:
        pass

    # Minimal fallback for the simple structure in site_config.yaml.
    # Handles top-level keys with scalar or list values nested one level.
    return _parse_simple_yaml(text)


def _parse_simple_yaml(text: str) -> dict:
    """Parse a very simple YAML subset (one level of nesting, scalars and
    lists).  Good enough for site_config.yaml but not a general parser."""
    result: dict = {}
    current_section = None
    current_key = None
    for raw_line in text.splitlines():
        stripped = raw_line.split("#")[0].rstrip()  # remove comments
        if not stripped:
            continue
        indent = len(raw_line) - len(raw_line.lstrip())
        content = stripped.strip()
        if indent == 0 and content.endswith(":"):
            current_section = content[:-1]
            result[current_section] = {}
            current_key = None
        elif indent > 0 and current_section is not None:
            if content.startswith("- "):
                val = content[2:].strip().strip('"').strip("'")
                if current_key and isinstance(result[current_section].get(current_key), list):
                    result[current_section][current_key].append(val)
            elif ":" in content:
                k, _, v = content.partition(":")
                k = k.strip()
                v = v.strip().strip('"').strip("'")
                if v == "[]":
                    result[current_section][k] = []
                elif v == "":
                    result[current_section][k] = []
                    current_key = k
                    continue
                else:
                    try:
                        v = int(v)
                    except ValueError:
                        try:
                            v = float(v)
                        except ValueError:
                            pass
                    result[current_section][k] = v
                current_key = k
    return result


# ---------------------------------------------------------------------------
# Git helpers
# ---------------------------------------------------------------------------

def git(*args: str, cwd: Path | None = None, check: bool = True) -> subprocess.CompletedProcess:
    """Run a git command, returning the CompletedProcess."""
    return subprocess.run(
        ["git", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=check,
    )


def find_repo_root() -> Path:
    """Detect the repository root via git."""
    proc = git("rev-parse", "--show-toplevel")
    return Path(proc.stdout.strip())


# ---------------------------------------------------------------------------
# Worktree management
# ---------------------------------------------------------------------------

def ensure_worktree(repo_root: Path, worktree_dir: str, branch: str) -> Path:
    """Ensure the gh-pages worktree exists and is on the right branch.

    Returns the absolute path to the worktree directory.
    """
    wt_path = repo_root / worktree_dir
    if wt_path.exists():
        # Verify it is on the expected branch.
        proc = git("rev-parse", "--abbrev-ref", "HEAD", cwd=wt_path, check=False)
        current_branch = proc.stdout.strip()
        if current_branch == branch:
            # Pull latest to avoid conflicts.
            git("pull", "--ff-only", "origin", branch, cwd=wt_path, check=False)
            return wt_path
        # Branch mismatch — remove stale worktree and re-create.
        print(f"Worktree branch mismatch (have {current_branch}, want {branch}). Re-creating.")
        git("worktree", "remove", "--force", str(wt_path), cwd=repo_root, check=False)

    # Create fresh worktree.
    print(f"Creating worktree at {wt_path} on branch {branch} ...")
    git("worktree", "add", "--detach", str(wt_path), cwd=repo_root)
    git("checkout", branch, cwd=wt_path)
    return wt_path


# ---------------------------------------------------------------------------
# Report / figure collection
# ---------------------------------------------------------------------------

def _matches_any(name: str, patterns: list[str]) -> bool:
    return any(fnmatch.fnmatch(name, p) for p in patterns)


def _prettify_filename(name: str) -> str:
    """Turn a filename stem into a human-readable title."""
    stem = Path(name).stem
    return stem.replace("_", " ").replace("-", " ").title()


def _extract_md_metadata(path: Path) -> tuple[str, str]:
    """Extract title (first # heading) and description (first paragraph) from
    a Markdown file.  Returns (title, description)."""
    title = _prettify_filename(path.name)
    description = ""
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return title, description
    lines = text.splitlines()
    found_heading = False
    para_lines: list[str] = []
    for line in lines:
        stripped = line.strip()
        if not found_heading and stripped.startswith("#"):
            title = stripped.lstrip("#").strip()
            found_heading = True
            continue
        if found_heading:
            if stripped == "":
                if para_lines:
                    break  # end of first paragraph
                continue
            para_lines.append(stripped)
    description = " ".join(para_lines)[:200]
    return title, description


def _extract_tex_title(tex_path: Path) -> str | None:
    """Try to find a meaningful title in a .tex file.

    Handles multi-line \\title{...} blocks common in sPHENIX reports, which
    often contain \\vspace, \\large, \\textbf wrappers around the actual title.
    Falls back to the first \\section{} heading if no \\title is found.
    """
    try:
        text = tex_path.read_text(errors="replace")
    except OSError:
        return None

    def _clean_latex(raw: str) -> str:
        """Strip LaTeX markup from a title string."""
        t = raw.replace("\\\\", " ").replace("\n", " ")
        t = re.sub(r"\\[a-zA-Z]+\*?\{?", "", t)
        t = t.replace("$", "").replace("{", "").replace("}", "")
        t = re.sub(r"\s+", " ", t)
        return t.strip()

    # First try: find the \LARGE \textbf{...} inside a \title block (sPHENIX style).
    m = re.search(r"\\LARGE\s*\\textbf\{([^}]+)\}", text)
    if m:
        title = _clean_latex(m.group(1))
        return title if title else None

    # Second try: simple single-line \title{...}.
    m = re.search(r"\\title\{([^%{}\n]+)\}", text)
    if m:
        title = _clean_latex(m.group(1))
        return title if title else None

    # Third try: first \section{...}.
    m = re.search(r"\\section\{([^}]+)\}", text)
    if m:
        return m.group(1).strip()

    return None


def _get_report_dirs(repo_root: Path, cfg: dict) -> list[Path]:
    """Return the list of report source directories from config."""
    rcfg = cfg.get("reports", {})
    # Support both source_dirs (list) and source_dir (single string).
    dirs = rcfg.get("source_dirs", None)
    if dirs is None:
        dirs = [rcfg.get("source_dir", "reports")]
    return [repo_root / d for d in dirs]


def compile_tex_reports(repo_root: Path, cfg: dict) -> list[str]:
    """Compile .tex reports across all report directories that are newer than
    their PDF.

    Returns a list of compiled filenames (stems).
    """
    rcfg = cfg.get("reports", {})
    exclude = rcfg.get("exclude", [])
    compiled: list[str] = []

    for source in _get_report_dirs(repo_root, cfg):
        if not source.is_dir():
            continue
        for tex in sorted(source.glob("*.tex")):
            # Skip excluded files (e.g., main.tex, helper .tex).
            if _matches_any(tex.name, exclude):
                continue
            # Skip empty files.
            if tex.stat().st_size == 0:
                continue
            pdf = tex.with_suffix(".pdf")
            # Compile if PDF missing or older than .tex source or any
            # included figure (symlinked PDFs in the same directory).
            needs_compile = not pdf.exists()
            if not needs_compile and pdf.stat().st_mtime < tex.stat().st_mtime:
                needs_compile = True
            if not needs_compile:
                # Check if any symlinked PDF figure is newer than the report PDF
                pdf_mtime = pdf.stat().st_mtime
                for fig in source.glob("*.pdf"):
                    if fig == pdf:
                        continue
                    # Resolve symlinks to get the real file's mtime
                    real = fig.resolve()
                    if real.exists() and real.stat().st_mtime > pdf_mtime:
                        needs_compile = True
                        break
            if not needs_compile:
                continue
            print(f"Compiling {source.name}/{tex.name} ...")
            # Run pdflatex twice for cross-references.
            for pass_num in range(2):
                proc = subprocess.run(
                    ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex.name],
                    cwd=source,
                    capture_output=True,
                    text=True,
                )
                if proc.returncode != 0 and pass_num == 0:
                    print(f"  WARNING: pdflatex failed for {tex.name}:")
                    for line in proc.stdout.splitlines()[-15:]:
                        print(f"    {line}")
                    break
            if pdf.exists():
                compiled.append(f"{source.name}/{tex.stem}")
                print(f"  -> {pdf.name} ({pdf.stat().st_size / 1024:.0f} KB)")
            else:
                print(f"  WARNING: {pdf.name} not produced")
    return compiled


def collect_reports(repo_root: Path, cfg: dict) -> list[dict]:
    """Collect report files from all configured source directories."""
    rcfg = cfg.get("reports", {})
    include = rcfg.get("include", ["*.pdf", "*.md"])
    exclude = rcfg.get("exclude", ["*.aux", "*.log", "*.out", "*.py"])

    reports: list[dict] = []
    seen_names: set[str] = set()

    for source in _get_report_dirs(repo_root, cfg):
        if not source.is_dir():
            continue

        # Label for provenance (e.g., "reports", "efficiencytool/reports").
        rel_source = source.relative_to(repo_root)

        for f in sorted(source.iterdir()):
            if not f.is_file():
                continue
            # Skip symlinks (e.g., symlinked figure PDFs in efficiencytool/reports).
            if f.is_symlink():
                continue
            if not _matches_any(f.name, include):
                continue
            if _matches_any(f.name, exclude):
                continue

            # Deduplicate by name (first dir wins).
            if f.name in seen_names:
                continue
            seen_names.add(f.name)

            mtime = datetime.fromtimestamp(f.stat().st_mtime, tz=timezone.utc)
            size_mb = f.stat().st_size / (1024 * 1024)

            # Extract metadata.
            title = _prettify_filename(f.name)
            description = f"From {rel_source}"

            if f.suffix == ".pdf":
                # Check for a companion .tex with the same stem.
                tex = f.with_suffix(".tex")
                tex_title = _extract_tex_title(tex) if tex.exists() else None
                if tex_title:
                    title = tex_title
            elif f.suffix == ".md":
                title, description = _extract_md_metadata(f)

            reports.append({
                "path": f,
                "name": f.name,
                "title": title,
                "description": description,
                "date": mtime,
                "size_mb": size_mb,
            })

    # Sort newest first.
    reports.sort(key=lambda r: r["date"], reverse=True)
    return reports


def collect_figures(repo_root: Path, cfg: dict) -> list[dict]:
    """Collect figure files according to the explicit allowlist in config."""
    fcfg = cfg.get("figures", {})
    source = repo_root / fcfg.get("source_dir", "plotting/figures")
    patterns = fcfg.get("include", [])

    if not source.is_dir() or not patterns:
        return []

    figures: list[dict] = []
    seen: set[str] = set()
    for pattern in patterns:
        for f in sorted(source.glob(pattern)):
            if not f.is_file() or f.name in seen:
                continue
            seen.add(f.name)
            figures.append({
                "path": f,
                "name": f.name,
                "size_mb": f.stat().st_size / (1024 * 1024),
            })

    figures.sort(key=lambda fig: fig["name"])
    return figures


# ---------------------------------------------------------------------------
# Size check
# ---------------------------------------------------------------------------

def check_total_size(reports: list[dict], figures: list[dict], max_mb: float) -> None:
    """Abort if the total deployment size exceeds the limit."""
    all_items = reports + figures
    total = sum(item["size_mb"] for item in all_items)
    if total > max_mb:
        print(f"ERROR: Total deployment size {total:.1f} MB exceeds limit of {max_mb:.0f} MB.")
        print("Largest files:")
        for item in sorted(all_items, key=lambda x: x["size_mb"], reverse=True)[:10]:
            print(f"  {item['size_mb']:6.2f} MB  {item['name']}")
        sys.exit(1)


# ---------------------------------------------------------------------------
# HTML generation
# ---------------------------------------------------------------------------

INDEX_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <link rel="stylesheet" href="style.css">
</head>
<body>
  <header class="site-header">
    <div class="container">
      <h1><span class="accent">sPHENIX</span> {title}</h1>
      <div class="subtitle">{subtitle}</div>
      <div class="subtitle">Last updated: {timestamp}</div>
    </div>
  </header>

  <main class="container">
    <section id="reports">
      <h2 class="section-title">Reports</h2>
{reports_block}
    </section>

    <section id="figures">
      <h2 class="section-title">Figures</h2>
{figures_block}
    </section>
  </main>

  <footer class="site-footer">
    <div class="container">
      Auto-generated by <code>deploy_pages.py</code> &middot;
      PPG12 &middot;
      <a href="{repo_url}">sPHENIX Collaboration</a>
    </div>
  </footer>
</body>
</html>
"""


def _build_reports_table(reports: list[dict]) -> str:
    if not reports:
        return "      <p>No reports available.</p>"
    rows: list[str] = []
    rows.append('      <table class="reports-table">')
    rows.append("        <thead>")
    rows.append("          <tr><th>Date</th><th>Title</th><th>Description</th><th>Download</th></tr>")
    rows.append("        </thead>")
    rows.append("        <tbody>")
    for r in reports:
        date_str = r["date"].strftime("%Y-%m-%d")
        title_esc = html.escape(r["title"])
        desc_esc = html.escape(r["description"])
        link_target = f"reports/{html.escape(r['name'])}"
        target_attr = ' target="_blank" rel="noopener"' if r["name"].endswith(".pdf") else ""
        size_str = f"{r['size_mb']:.1f} MB" if r["size_mb"] >= 0.1 else f"{r['size_mb'] * 1024:.0f} KB"
        rows.append(f'          <tr>')
        rows.append(f'            <td>{date_str}</td>')
        rows.append(f'            <td><a href="{link_target}"{target_attr}>{title_esc}</a></td>')
        rows.append(f'            <td>{desc_esc}</td>')
        rows.append(f'            <td>{size_str}</td>')
        rows.append(f'          </tr>')
    rows.append("        </tbody>")
    rows.append("      </table>")
    return "\n".join(rows)


def _build_figures_grid(figures: list[dict]) -> str:
    if not figures:
        return "      <p>No figures selected for publication. Edit <code>site_config.yaml</code> to add figures.</p>"
    cards: list[str] = []
    cards.append('      <div class="figures-grid">')
    for fig in figures:
        name_esc = html.escape(fig["name"])
        link = f"figures/{html.escape(fig['name'])}"
        # For PDF figures, link opens in new tab; for images, show inline.
        if fig["name"].lower().endswith(".pdf"):
            cards.append(f'        <div class="figure-card">')
            cards.append(f'          <a href="{link}" target="_blank" rel="noopener">')
            cards.append(f'            <div class="figure-placeholder">PDF</div>')
            cards.append(f'          </a>')
            cards.append(f'          <p class="caption">{name_esc}</p>')
            cards.append(f'        </div>')
        else:
            cards.append(f'        <div class="figure-card">')
            cards.append(f'          <a href="{link}" target="_blank" rel="noopener">')
            cards.append(f'            <img src="{link}" alt="{name_esc}" loading="lazy">')
            cards.append(f'          </a>')
            cards.append(f'          <p class="caption">{name_esc}</p>')
            cards.append(f'        </div>')
    cards.append('      </div>')
    return "\n".join(cards)


def generate_index_html(reports: list[dict], figures: list[dict], cfg: dict) -> str:
    """Render the index.html from collected metadata."""
    site = cfg.get("site", {})
    title = site.get("title", "PPG12 Analysis Results")
    subtitle = site.get("subtitle", "Isolated Photon Cross-Section, pp \u221As = 200 GeV \u2014 sPHENIX")
    repo_url = site.get("repo_url", "https://github.com/Shuonli/ppg12")
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    return INDEX_HTML_TEMPLATE.format(
        title=html.escape(title),
        subtitle=html.escape(subtitle),
        timestamp=timestamp,
        repo_url=html.escape(repo_url),
        reports_block=_build_reports_table(reports),
        figures_block=_build_figures_grid(figures),
    )


# ---------------------------------------------------------------------------
# Deploy
# ---------------------------------------------------------------------------

def copy_to_worktree(
    wt_path: Path,
    reports: list[dict],
    figures: list[dict],
) -> list[str]:
    """Copy collected files into the worktree, returning a list of deployed
    relative paths."""
    deployed: list[str] = []

    # Clean and re-populate reports/.
    reports_dest = wt_path / "reports"
    if reports_dest.exists():
        shutil.rmtree(reports_dest)
    reports_dest.mkdir(parents=True, exist_ok=True)
    for r in reports:
        dst = reports_dest / r["name"]
        shutil.copy2(r["path"], dst)
        deployed.append(f"reports/{r['name']}")

    # Clean and re-populate figures/.
    figures_dest = wt_path / "figures"
    if figures_dest.exists():
        shutil.rmtree(figures_dest)
    figures_dest.mkdir(parents=True, exist_ok=True)
    for fig in figures:
        dst = figures_dest / fig["name"]
        shutil.copy2(fig["path"], dst)
        deployed.append(f"figures/{fig['name']}")

    return deployed


def deploy(wt_path: Path, message: str) -> bool:
    """Stage, commit, and push.  Returns True if changes were pushed."""
    git("add", "-A", cwd=wt_path)

    # Check for changes.
    proc = git("diff", "--cached", "--quiet", cwd=wt_path, check=False)
    if proc.returncode == 0:
        print("No changes to deploy.")
        return False

    git("commit", "-m", message, cwd=wt_path)
    git("push", "origin", "gh-pages", cwd=wt_path)
    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Deploy PPG12 analysis results to GitHub Pages.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be deployed without touching the worktree.",
    )
    parser.add_argument(
        "--message",
        default=None,
        help="Custom commit message (default: 'Update site: YYYY-MM-DD').",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Path to site_config.yaml (default: scripts/site_config.yaml relative to repo root).",
    )
    args = parser.parse_args()

    # ----- Locate repo and config -----------------------------------------
    repo_root = find_repo_root()
    config_path = Path(args.config) if args.config else repo_root / "scripts" / "site_config.yaml"

    if config_path.exists():
        print(f"Loading config from {config_path}")
        cfg = _load_yaml(config_path)
    else:
        print(f"Config not found at {config_path}, using defaults.")
        cfg = {}

    # ----- Compile .tex reports if needed -----------------------------------
    compiled = compile_tex_reports(repo_root, cfg)
    if compiled:
        print(f"Compiled {len(compiled)} report(s): {', '.join(compiled)}")

    # ----- Collect files --------------------------------------------------
    reports = collect_reports(repo_root, cfg)
    figures = collect_figures(repo_root, cfg)

    max_mb = cfg.get("figures", {}).get("max_total_mb", 50)
    check_total_size(reports, figures, max_mb)

    total_mb = sum(r["size_mb"] for r in reports) + sum(f["size_mb"] for f in figures)

    # ----- Dry-run output -------------------------------------------------
    if args.dry_run:
        print(f"\n{'=' * 60}")
        print(f"DRY RUN — nothing will be written")
        print(f"{'=' * 60}")
        print(f"\nReports ({len(reports)}):")
        for r in reports:
            print(f"  {r['size_mb']:6.2f} MB  {r['name']}  [{r['title']}]")
        print(f"\nFigures ({len(figures)}):")
        for fig in figures:
            print(f"  {fig['size_mb']:6.2f} MB  {fig['name']}")
        print(f"\nTotal: {total_mb:.2f} MB / {max_mb:.0f} MB limit")

        # Generate and show the index path that would be written.
        index_html = generate_index_html(reports, figures, cfg)
        preview_path = repo_root / "_dry_run_index.html"
        preview_path.write_text(index_html)
        print(f"\nGenerated index.html preview: {preview_path}")
        return

    # ----- Ensure worktree ------------------------------------------------
    deploy_cfg = cfg.get("deploy", {})
    wt_dir = deploy_cfg.get("worktree_dir", "_gh-pages-worktree")
    branch = deploy_cfg.get("branch", "gh-pages")
    wt_path = ensure_worktree(repo_root, wt_dir, branch)

    # ----- Copy files and generate index ----------------------------------
    deployed = copy_to_worktree(wt_path, reports, figures)

    index_html = generate_index_html(reports, figures, cfg)
    (wt_path / "index.html").write_text(index_html)
    deployed.append("index.html")

    print(f"\nPrepared {len(deployed)} files ({total_mb:.2f} MB):")
    for d in deployed:
        print(f"  {d}")

    # ----- Commit and push ------------------------------------------------
    message = args.message or f"Update site: {datetime.now(timezone.utc).strftime('%Y-%m-%d')}"
    pushed = deploy(wt_path, message)

    if pushed:
        print(f"\nDeployed to https://shuonli.github.io/ppg12/")
    print("Done.")


if __name__ == "__main__":
    main()
