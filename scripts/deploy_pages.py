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
from zoneinfo import ZoneInfo

_EASTERN = ZoneInfo("America/New_York")
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


def _clean_latex(raw: str) -> str:
    """Strip LaTeX markup from a snippet so it renders cleanly as plain HTML."""
    t = raw.replace("\n", " ")

    # \\ followed (after optional [Npt] dim) by \large/\Large acts as a
    # subtitle separator.
    def _sep(m: "re.Match[str]") -> str:
        prev = t[m.start() - 1] if m.start() > 0 else ""
        return " " if prev in ":-" else ": "
    t = re.sub(r"\\\\\s*(?:\[-?[\d.]+pt\])?\s*\\(?:large|Large|LARGE)\b\s*",
               _sep, t)
    t = t.replace("\\\\", " ")
    t = re.sub(r"\[-?[\d.]+pt\]", "", t)

    t = t.replace("\\ ", " ").replace("\\,", " ")
    # Thin-space / negative-space / kern shortcuts: \!  \;  \:  \>  \<  \@.
    t = re.sub(r"\\[!;:<>@]", "", t)
    t = t.replace("\\%", "%").replace("\\&", "&")
    t = t.replace("---", "\u2014").replace("--", "\u2013")

    # PPG12 macros — expand BEFORE the generic command strip.
    macro_map = {
        r"\\etg(?![a-zA-Z])":       "E_T^\u03b3",
        r"\\etgtruth(?![a-zA-Z])":  "E_T^(\u03b3,truth)",
        r"\\isoET(?![a-zA-Z])":     "E_T^iso",
        r"\\isoETtruth(?![a-zA-Z])":"E_T^(iso,truth)",
        r"\\DR(?![a-zA-Z])":        "\u0394R",
        r"\\Delta(?![a-zA-Z])":     "\u0394",
        r"\\gamma(?![a-zA-Z])":     "\u03b3",
        r"\\pp(?![a-zA-Z])":        "pp",
        r"\\sqs(?![a-zA-Z])":       "\u221As",
        r"\\purity(?![a-zA-Z])":    "P",
        r"\\eff(?![a-zA-Z])":       "E",
        r"\\lumi(?![a-zA-Z])":      "L",
        r"\\ET(?![a-zA-Z])":        "E_T",
    }
    for pat, rep in macro_map.items():
        t = re.sub(pat, rep, t)

    t = re.sub(r"\\text(?:bf|it|rm|sf|tt|sc|md|up|sl)\s*\{([^{}]*)\}", r"\1", t)
    t = re.sub(r"\\(?:large|Large|LARGE|small|Huge|huge|normalsize|bf|it|rm|sf)\b",
               "", t)

    # Strip environment begin/end wrappers (e.g. \begin{abstract}, \end{abstract}).
    t = re.sub(r"\\(?:begin|end)\{[^}]*\}", "", t)

    t = re.sub(r"\\[a-zA-Z]+\*?\{?", "", t)
    t = t.replace("$", "").replace("{", "").replace("}", "")
    t = re.sub(r"\s+", " ", t)
    t = re.sub(r"\s+([,.;:])", r"\1", t)
    t = re.sub(r"([\u0394\u03b3])\s+([a-zA-Z])", r"\1\2", t)
    return t.strip()


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

    # 2b: multi-line \title{...} with possible nested braces.
    m_start = re.search(r"\\title\{", text)
    if m_start:
        depth, i = 1, m_start.end()
        while i < len(text) and depth > 0:
            if text[i] == "{":
                depth += 1
            elif text[i] == "}":
                depth -= 1
            i += 1
        if depth == 0:
            raw = text[m_start.end() : i - 1]
            title = _clean_latex(raw)
            if title:
                return title

    # Third try: first \section{...}.
    m = re.search(r"\\section\{([^}]+)\}", text)
    if m:
        return m.group(1).strip()

    return None


def _extract_tex_tags(tex_path: Path) -> list[str]:
    """Return explicit tags declared via a `%% tags: a, b, c` line at the top
    of a .tex file.  Empty list if no such directive."""
    try:
        text = tex_path.read_text(errors="replace")
    except OSError:
        return []
    m = re.search(r"^\s*%%\s*tags\s*:\s*(.+)$", text, re.MULTILINE | re.IGNORECASE)
    if not m:
        return []
    return [t.strip().lower() for t in m.group(1).split(",") if t.strip()]


# Display order + keyword inference for grouping reports.
# Each entry: (tag_id, human_label, filename/title substrings that trigger it).
TAG_DEFINITIONS: list[tuple[str, str, list[str]]] = [
    ("cross-section",      "Final Cross-Section",          ["final_cross", "cross_section", "cross-section"]),
    ("double-interaction", "Double Interaction / Pileup",  ["double_interaction", "double-interaction", "pileup"]),
    ("bdt",                "BDT Training & Systematics",   ["bdt"]),
    ("shower-shape",       "Shower Shape",                 ["showershape", "shower_shape", "shower-shape", "cluster_size", "converted_photon"]),
    ("efficiency",         "Efficiency & Energy Response", ["efficiency", "eta_migration", "deltar", "delta_r", "energy_response"]),
    ("isolation",          "Isolation",                    ["isolation", "iso_et"]),
    ("purity",             "Purity",                       ["purity"]),
    ("vertex",             "Vertex",                       ["vertex"]),
    ("selection",          "Selection & Config Comparisons", ["selection", "comparison"]),
    ("updates",            "Status Updates & Summaries",   ["post_preliminary", "updates"]),
]
TAG_LABELS: dict[str, str] = {tid: lbl for tid, lbl, _ in TAG_DEFINITIONS}
TAG_ORDER: list[str] = [tid for tid, _, _ in TAG_DEFINITIONS] + ["other"]
TAG_LABELS["other"] = "Other"


def _infer_tag(filename: str, title: str) -> str:
    """Guess a report's primary tag from filename/title substrings."""
    key = filename.lower()
    for tid, _lbl, substrings in TAG_DEFINITIONS:
        for s in substrings:
            if s in key:
                return tid
    # Fall back to checking the title (looser match on a space-separated form).
    tlower = title.lower()
    for tid, _lbl, substrings in TAG_DEFINITIONS:
        for s in substrings:
            if s.replace("_", " ").replace("-", " ") in tlower:
                return tid
    return "other"


def _extract_tex_description(tex_path: Path, max_len: int = 200) -> str | None:
    """Extract a short description from a .tex report.

    Preference: \\begin{abstract}...\\end{abstract}. Falls back to the first
    substantive paragraph after \\begin{document}, skipping \\maketitle,
    \\tableofcontents, comments, and lone macro invocations.
    """
    try:
        text = tex_path.read_text(errors="replace")
    except OSError:
        return None

    m = re.search(r"\\begin\{abstract\}(.+?)\\end\{abstract\}", text, re.DOTALL)
    if m:
        cleaned = _clean_latex(m.group(1))
        if cleaned:
            return (cleaned[: max_len - 1].rstrip() + "\u2026") if len(cleaned) > max_len else cleaned

    doc_start = re.search(r"\\begin\{document\}", text)
    if doc_start:
        body = text[doc_start.end() :]
        # Drop comment-only lines, then split into paragraphs.
        body = "\n".join(ln for ln in body.splitlines() if not ln.lstrip().startswith("%"))
        for para in re.split(r"\n\s*\n", body):
            para = para.strip()
            if not para:
                continue
            if re.match(r"^\\(?:maketitle|tableofcontents|clearpage|newpage|thispagestyle|pagestyle|noindent|vspace|hspace)\b",
                        para):
                continue
            cleaned = _clean_latex(para)
            if len(cleaned) >= 40:
                return (cleaned[: max_len - 1].rstrip() + "\u2026") if len(cleaned) > max_len else cleaned
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

            mtime = datetime.fromtimestamp(f.stat().st_mtime, tz=_EASTERN)
            size_mb = f.stat().st_size / (1024 * 1024)

            # Extract metadata.
            title = _prettify_filename(f.name)
            description = f"From {rel_source}"

            explicit_tags: list[str] = []
            if f.suffix == ".pdf":
                # Check for a companion .tex with the same stem.
                tex = f.with_suffix(".tex")
                if tex.exists():
                    tex_title = _extract_tex_title(tex)
                    if tex_title:
                        title = tex_title
                    tex_desc = _extract_tex_description(tex)
                    if tex_desc:
                        description = tex_desc
                    explicit_tags = _extract_tex_tags(tex)
            elif f.suffix == ".md":
                title, description = _extract_md_metadata(f)

            # Primary tag: first explicit tag, else infer from name/title.
            primary_tag = explicit_tags[0] if explicit_tags else _infer_tag(f.name, title)
            if primary_tag not in TAG_LABELS:
                primary_tag = "other"

            reports.append({
                "path": f,
                "name": f.name,
                "title": title,
                "description": description,
                "date": mtime,
                "size_mb": size_mb,
                "tag": primary_tag,
            })

    # Sort newest first.
    reports.sort(key=lambda r: r["date"], reverse=True)
    return reports


def collect_past_deploys(repo_root: Path) -> list[dict]:
    """List every site/* tag with its commit date, short hash, and subject.

    Lightweight tags delegate to the pointed-to commit, so committerdate and
    contents:subject resolve correctly for backfilled tags.
    """
    proc = git(
        "for-each-ref",
        "--format=%(refname:short)%09%(objectname:short)%09%(committerdate:iso-strict)%09%(contents:subject)",
        "refs/tags/site/*",
        cwd=repo_root,
        check=False,
    )
    entries: list[dict] = []
    for line in proc.stdout.splitlines():
        parts = line.split("\t", 3)
        if len(parts) != 4:
            continue
        tag, sha, iso_date, subject = parts
        entries.append({
            "tag": tag,
            "sha": sha,
            "date": iso_date,
            "message": subject.strip(),
        })
    entries.sort(key=lambda e: e["date"], reverse=True)
    # Annotate each entry with the tag of the chronologically prior deploy
    # so the UI can render a GitHub compare link.
    for i, e in enumerate(entries):
        e["prev_tag"] = entries[i + 1]["tag"] if i + 1 < len(entries) else None
    return entries


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

# ---------------------------------------------------------------------------
# PDF thumbnails (page 1 preview for each report)
# ---------------------------------------------------------------------------

def generate_thumbnails(reports: list[dict], out_dir: Path, dpi: int = 72) -> int:
    """Render page 1 of each PDF to a PNG thumbnail in ``out_dir``.

    Uses PyMuPDF (``fitz``).  Skips when the existing thumbnail is newer than
    the source PDF.  Each report dict gains a ``thumb`` key pointing at the
    thumbnail filename (relative to out_dir's parent, i.e. ``_thumbs/x.png``)
    if rendering succeeded; otherwise ``thumb`` stays unset.
    """
    try:
        import fitz
    except ImportError:
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    zoom = dpi / 72.0
    matrix = fitz.Matrix(zoom, zoom)
    generated = 0
    for r in reports:
        if not r["name"].lower().endswith(".pdf"):
            continue
        thumb_path = out_dir / f"{Path(r['name']).stem}.png"
        r["thumb"] = f"_thumbs/{thumb_path.name}"
        # Skip if thumb exists and is newer than the PDF.
        if thumb_path.exists() and thumb_path.stat().st_mtime >= r["path"].stat().st_mtime:
            continue
        try:
            doc = fitz.open(r["path"])
            if doc.page_count == 0:
                doc.close()
                continue
            page = doc[0]
            pix = page.get_pixmap(matrix=matrix, alpha=False)
            pix.save(thumb_path)
            doc.close()
            generated += 1
        except Exception as e:  # noqa: BLE001
            print(f"  WARNING: thumbnail failed for {r['name']}: {e}")
            r.pop("thumb", None)
    if generated:
        print(f"  thumbnails: {generated} regenerated, {sum(1 for r in reports if 'thumb' in r) - generated} reused")
    return generated


# ---------------------------------------------------------------------------
# Per-report detail pages
# ---------------------------------------------------------------------------

DETAIL_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title} — PPG12</title>
  <link rel="stylesheet" href="../style.css">
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
</head>
<body>
  <header class="site-header">
    <div class="container">
      <h1><a href="../index.html" style="color:#B9D9EB;text-decoration:none">&larr; sPHENIX PPG12</a></h1>
      <div class="subtitle">{tag_label}</div>
    </div>
  </header>

  <main class="container detail-main">
    <h2 class="detail-title">{title}</h2>
    <div class="detail-meta">
      <span class="detail-date">{date}</span>
      <span class="detail-size">{size}</span>
      <a class="detail-download pdf-link" href="{pdf_name}" target="_blank" rel="noopener">Download PDF</a>
{source_link}
    </div>

    <p class="detail-abstract">{abstract}</p>

    <iframe class="detail-iframe" src="{pdf_name}" title="{title}"></iframe>

{versions_block}
  </main>

  <footer class="site-footer">
    <div class="container">
      <a href="../index.html">Back to the index</a> &middot; PPG12 &middot;
      <a href="{repo_url}">sPHENIX Collaboration</a>
    </div>
  </footer>
</body>
</html>
"""


def _build_versions_block(report_name: str, deploys: list[dict], repo_url: str) -> str:
    """List the deploys (newest first) that contain this report name."""
    if not deploys:
        return ""
    rows: list[str] = []
    rows.append('    <section class="detail-versions">')
    rows.append('      <h3>Past versions of this report</h3>')
    rows.append('      <ul class="version-list">')
    for d in deploys[:15]:
        tag_esc = html.escape(d["tag"])
        date_esc = html.escape(d["date"][:10])
        url = f"{repo_url}/raw/{tag_esc}/reports/{report_name}"
        rows.append(f'        <li><a href="{url}" target="_blank" rel="noopener">'
                    f'{date_esc} · {tag_esc}</a></li>')
    if len(deploys) > 15:
        rows.append(f'        <li class="muted">&hellip; {len(deploys) - 15} older deploys on '
                    f'<a href="{repo_url}/commits/gh-pages" target="_blank" rel="noopener">GitHub</a></li>')
    rows.append('      </ul>')
    rows.append('    </section>')
    return "\n".join(rows)


def _deploys_containing_report(
    repo_root: Path,
    report_name: str,
    deploys: list[dict],
    limit: int = 30,
) -> list[dict]:
    """Check which deploys (tags) actually contain ``reports/<report_name>``.

    Uses ``git cat-file -e`` which is cheap.  Returns deploys sorted newest first.
    """
    present: list[dict] = []
    for d in deploys:
        ref = f"{d['tag']}:reports/{report_name}"
        proc = git("cat-file", "-e", ref, cwd=repo_root, check=False)
        if proc.returncode == 0:
            present.append(d)
            if len(present) >= limit:
                break
    return present


def generate_detail_pages(
    reports: list[dict],
    deploys: list[dict],
    wt_path: Path,
    repo_root: Path,
    repo_url: str,
) -> int:
    """Write a detail HTML page alongside each PDF in the worktree."""
    reports_dir = wt_path / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    for r in reports:
        if not r["name"].lower().endswith(".pdf"):
            continue
        stem = Path(r["name"]).stem
        out_path = reports_dir / f"{stem}.html"

        size_str = (f"{r['size_mb']:.1f} MB" if r["size_mb"] >= 0.1
                    else f"{r['size_mb'] * 1024:.0f} KB")
        tag_id = r.get("tag", "other")
        tag_label = TAG_LABELS.get(tag_id, "Report")
        date_str = r["date"].strftime("%Y-%m-%d")

        # Link back to source .tex on GitHub if we can find it on main HEAD.
        source_link = ""
        source_rel = _find_source_tex(repo_root, stem)
        if source_rel:
            source_link = (f'      <a class="detail-source" '
                           f'href="{repo_url}/blob/main/{source_rel}" '
                           f'target="_blank" rel="noopener">View .tex on GitHub</a>')

        present = _deploys_containing_report(repo_root, r["name"], deploys)
        versions_block = _build_versions_block(r["name"], present, repo_url)

        out_path.write_text(DETAIL_HTML_TEMPLATE.format(
            title=html.escape(r["title"]),
            tag_label=html.escape(tag_label),
            date=html.escape(date_str),
            size=html.escape(size_str),
            pdf_name=html.escape(r["name"]),
            source_link=source_link,
            abstract=html.escape(r["description"]),
            versions_block=versions_block,
            repo_url=html.escape(repo_url),
        ))
        count += 1
    return count


def _find_source_tex(repo_root: Path, stem: str) -> str | None:
    """Return the repo-relative path to ``<stem>.tex`` if it exists anywhere."""
    # git ls-files is fast and scoped to tracked files.
    proc = git("ls-files", f"*{stem}.tex", cwd=repo_root, check=False)
    for line in proc.stdout.splitlines():
        if Path(line).stem == stem:
            return line
    return None


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
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
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
{recent_rail}
{reports_section}
{figures_section}
{deploys_section}
  </main>

  <footer class="site-footer">
    <div class="container">
      Auto-generated by <code>deploy_pages.py</code>{build_info} &middot;
      PPG12 &middot;
      <a href="{repo_url}">sPHENIX Collaboration</a>
    </div>
  </footer>

  <script>
  // Zero-dependency filter: as the user types in the reports-filter input,
  // hide rows whose text content does not match. Empty groups collapse.
  (function () {{
    var input = document.getElementById('reports-filter');
    if (!input) return;
    input.addEventListener('input', function () {{
      var q = input.value.trim().toLowerCase();
      document.querySelectorAll('#reports .reports-group').forEach(function (grp) {{
        var visibleRows = 0;
        grp.querySelectorAll('tbody tr').forEach(function (tr) {{
          var hit = q === '' || tr.textContent.toLowerCase().indexOf(q) !== -1;
          tr.style.display = hit ? '' : 'none';
          if (hit) visibleRows += 1;
        }});
        grp.style.display = (visibleRows === 0 && q !== '') ? 'none' : '';
      }});
      // Ungrouped flat table case.
      document.querySelectorAll('#reports > table tbody tr').forEach(function (tr) {{
        var hit = q === '' || tr.textContent.toLowerCase().indexOf(q) !== -1;
        tr.style.display = hit ? '' : 'none';
      }});
    }});
  }})();
  </script>
</body>
</html>
"""


def _render_reports_table(reports: list[dict]) -> list[str]:
    """Inner table body for a single group of reports (no wrapping section)."""
    rows: list[str] = []
    rows.append('      <table class="reports-table">')
    rows.append("        <thead>")
    rows.append("          <tr><th>Preview</th><th>Date</th><th>Title</th><th>Description</th><th>Download</th></tr>")
    rows.append("        </thead>")
    rows.append("        <tbody>")
    for r in reports:
        date_str = r["date"].strftime("%Y-%m-%d")
        title_esc = html.escape(r["title"])
        desc_esc = html.escape(r["description"])
        name_esc = html.escape(r["name"])
        is_pdf = r["name"].endswith(".pdf")
        stem = Path(r["name"]).stem
        detail_href = f"reports/{html.escape(stem)}.html" if is_pdf else f"reports/{name_esc}"
        pdf_href = f"reports/{name_esc}"
        size_str = f"{r['size_mb']:.1f} MB" if r["size_mb"] >= 0.1 else f"{r['size_mb'] * 1024:.0f} KB"

        thumb_rel = r.get("thumb")
        if thumb_rel:
            thumb_href = f"reports/{html.escape(thumb_rel)}"
            preview = (f'<a href="{detail_href}"><img class="report-thumb" '
                       f'src="{thumb_href}" alt="" loading="lazy"></a>')
        else:
            preview = '<span class="report-thumb-placeholder">PDF</span>' if is_pdf else ""

        rows.append('          <tr>')
        rows.append(f'            <td data-label="Preview" class="preview-cell">{preview}</td>')
        rows.append(f'            <td data-label="Date">{date_str}</td>')
        rows.append(f'            <td data-label="Title"><a href="{detail_href}">{title_esc}</a></td>')
        rows.append(f'            <td data-label="Description">{desc_esc}</td>')
        rows.append(f'            <td data-label="Size"><a href="{pdf_href}" target="_blank" rel="noopener" class="pdf-link">PDF&nbsp;·&nbsp;{size_str}</a></td>')
        rows.append('          </tr>')
    rows.append("        </tbody>")
    rows.append("      </table>")
    return rows


def _build_reports_section(reports: list[dict]) -> str:
    if not reports:
        return ""
    buckets: dict[str, list[dict]] = {}
    for r in reports:
        buckets.setdefault(r.get("tag", "other"), []).append(r)

    rows: list[str] = ['    <section id="reports">']
    rows.append('      <div class="reports-header">')
    rows.append('        <h2 class="section-title">Reports</h2>')
    rows.append('        <input id="reports-filter" type="search" '
                'placeholder="Filter by title, description, tag, date…" '
                'aria-label="Filter reports">')
    rows.append('      </div>')

    if len(buckets) <= 1:
        rows.extend(_render_reports_table(reports))
    else:
        for tag_id in TAG_ORDER:
            group = buckets.get(tag_id)
            if not group:
                continue
            label = TAG_LABELS.get(tag_id, tag_id.title())
            rows.append(f'      <div class="reports-group" id="group-{html.escape(tag_id)}">')
            rows.append(f'        <h3 class="group-title">{html.escape(label)}'
                        f' <span class="group-count">({len(group)})</span></h3>')
            rows.extend(_render_reports_table(group))
            rows.append('      </div>')
    rows.append("    </section>")
    return "\n".join(rows)


def _build_recent_rail(reports: list[dict], n: int = 3) -> str:
    """Horizontal strip of the ``n`` most recent reports (thumbnail + meta)."""
    pdfs = [r for r in reports if r["name"].lower().endswith(".pdf")]
    if not pdfs:
        return ""
    recent = pdfs[:n]
    rows: list[str] = ['    <section id="recent" class="recent-rail">']
    rows.append('      <h2 class="section-title">Recently updated</h2>')
    rows.append('      <div class="recent-grid">')
    for r in recent:
        stem = Path(r["name"]).stem
        detail_href = f"reports/{html.escape(stem)}.html"
        thumb_rel = r.get("thumb")
        thumb_html = (f'<img class="recent-thumb" src="reports/{html.escape(thumb_rel)}" alt="" loading="lazy">'
                      if thumb_rel else '<div class="recent-thumb-placeholder">PDF</div>')
        date_str = r["date"].strftime("%Y-%m-%d")
        title_esc = html.escape(r["title"])
        tag_label = html.escape(TAG_LABELS.get(r.get("tag", "other"), ""))
        rows.append('        <a class="recent-card" href="' + detail_href + '">')
        rows.append(f'          {thumb_html}')
        rows.append(f'          <div class="recent-body">')
        rows.append(f'            <div class="recent-tag">{tag_label}</div>')
        rows.append(f'            <div class="recent-title">{title_esc}</div>')
        rows.append(f'            <div class="recent-date">{date_str}</div>')
        rows.append('          </div>')
        rows.append('        </a>')
    rows.append('      </div>')
    rows.append('    </section>')
    return "\n".join(rows)


def _build_figures_section(figures: list[dict]) -> str:
    if not figures:
        return ""
    cards: list[str] = []
    cards.append('    <section id="figures">')
    cards.append('      <h2 class="section-title">Figures</h2>')
    cards.append('      <div class="figures-grid">')
    for fig in figures:
        name_esc = html.escape(fig["name"])
        link = f"figures/{html.escape(fig['name'])}"
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
    cards.append('    </section>')
    return "\n".join(cards)


_MONTH_LABELS = ["", "January", "February", "March", "April", "May", "June",
                 "July", "August", "September", "October", "November", "December"]


def _month_label(ym: str) -> str:
    """Turn 'YYYY-MM' into 'Month YYYY'."""
    try:
        y, m = ym.split("-")
        return f"{_MONTH_LABELS[int(m)]} {y}"
    except (ValueError, IndexError):
        return ym


def _render_deploys_rows(deploys: list[dict], repo_url: str) -> list[str]:
    """Inner <tr> rows for a single month bucket (no wrapping <table>)."""
    rows: list[str] = []
    for d in deploys:
        tag_esc = html.escape(d["tag"])
        date_esc = html.escape(d["date"][:16].replace("T", " "))
        sha_esc = html.escape(d["sha"])
        msg_esc = html.escape(d["message"])
        tree_url = f"{repo_url}/tree/{tag_esc}"
        commit_url = f"{repo_url}/commit/{sha_esc}"
        prev = d.get("prev_tag")
        if prev:
            cmp_url = f"{repo_url}/compare/{html.escape(prev)}...{tag_esc}"
            diff_cell = (f'<a href="{cmp_url}" target="_blank" rel="noopener" '
                         f'title="Compare with {html.escape(prev)}">diff</a>')
        else:
            diff_cell = '<span class="muted">\u2014</span>'
        rows.append('            <tr>')
        rows.append(f'              <td data-label="Date">{date_esc}</td>')
        rows.append(f'              <td data-label="Tag"><a href="{tree_url}" target="_blank" rel="noopener">{tag_esc}</a></td>')
        rows.append(f'              <td data-label="Commit"><a href="{commit_url}" target="_blank" rel="noopener">{sha_esc}</a></td>')
        rows.append(f'              <td data-label="&Delta;">{diff_cell}</td>')
        rows.append(f'              <td data-label="Message">{msg_esc}</td>')
        rows.append('            </tr>')
    return rows


def _build_deploys_section(deploys: list[dict], repo_url: str) -> str:
    if not deploys:
        return ""
    # Bucket by YYYY-MM, preserving newest-first order.
    months: dict[str, list[dict]] = {}
    for d in deploys:
        key = d["date"][:7]  # YYYY-MM
        months.setdefault(key, []).append(d)
    month_keys = list(months.keys())  # already in newest-first order thanks to caller sort

    rows: list[str] = []
    rows.append('    <section id="past-deploys">')
    rows.append('      <h2 class="section-title">Past deploys</h2>')
    rows.append('      <details class="deploys-details">')
    rows.append(f'        <summary>Browse {len(deploys)} past deploy(s) across {len(months)} month(s)</summary>')

    for i, key in enumerate(month_keys):
        bucket = months[key]
        open_attr = " open" if i == 0 else ""
        label = _month_label(key)
        rows.append(f'        <details class="deploys-month"{open_attr}>')
        rows.append(f'          <summary>{html.escape(label)} '
                    f'<span class="group-count">({len(bucket)})</span></summary>')
        rows.append('          <table class="deploys-table">')
        rows.append('            <thead>')
        rows.append('              <tr><th>Date</th><th>Tag</th><th>Commit</th>'
                    '<th title="Compare with previous deploy">&Delta;</th><th>Message</th></tr>')
        rows.append('            </thead>')
        rows.append('            <tbody>')
        rows.extend(_render_deploys_rows(bucket, repo_url))
        rows.append('            </tbody>')
        rows.append('          </table>')
        rows.append('        </details>')

    rows.append('        <p class="deploys-footnote">')
    rows.append('          Tag links open the full site snapshot; &ldquo;diff&rdquo; shows what changed vs. the previous deploy. '
                f'See also <a href="{repo_url}/tags" target="_blank" rel="noopener">all tags</a> '
                f'and the <a href="{repo_url}/commits/gh-pages" target="_blank" rel="noopener">gh-pages commit history</a>.')
    rows.append('        </p>')
    rows.append('      </details>')
    rows.append('    </section>')
    return "\n".join(rows)


def _build_build_info(repo_root: Path, repo_url: str) -> str:
    """Return a " · Built from main@<sha>" snippet for the footer, or ""."""
    proc = git("rev-parse", "--short", "HEAD", cwd=repo_root, check=False)
    sha = proc.stdout.strip()
    if not sha or proc.returncode != 0:
        return ""
    branch_proc = git("rev-parse", "--abbrev-ref", "HEAD", cwd=repo_root, check=False)
    branch = branch_proc.stdout.strip() or "main"
    commit_url = f"{repo_url}/commit/{sha}"
    return (f' &middot; Built from <code>{html.escape(branch)}@'
            f'<a href="{html.escape(commit_url)}" target="_blank" rel="noopener">'
            f'{html.escape(sha)}</a></code>')


def generate_index_html(
    reports: list[dict],
    figures: list[dict],
    deploys: list[dict],
    cfg: dict,
    repo_root: Path,
) -> str:
    """Render the index.html from collected metadata."""
    site = cfg.get("site", {})
    title = site.get("title", "PPG12 Analysis Results")
    subtitle = site.get("subtitle", "Isolated Photon Cross-Section, pp \u221As = 200 GeV \u2014 sPHENIX")
    repo_url = site.get("repo_url", "https://github.com/Shuonli/ppg12")
    timestamp = datetime.now(_EASTERN).strftime("%Y-%m-%d %H:%M %Z")

    return INDEX_HTML_TEMPLATE.format(
        title=html.escape(title),
        subtitle=html.escape(subtitle),
        timestamp=timestamp,
        repo_url=html.escape(repo_url),
        recent_rail=_build_recent_rail(reports),
        reports_section=_build_reports_section(reports),
        figures_section=_build_figures_section(figures),
        deploys_section=_build_deploys_section(deploys, repo_url),
        build_info=_build_build_info(repo_root, repo_url),
    )


# ---------------------------------------------------------------------------
# Deploy
# ---------------------------------------------------------------------------

def _needs_copy(src: Path, dst: Path) -> bool:
    """True if dst is missing or its size/mtime differs from src."""
    if not dst.exists():
        return True
    ss, ds = src.stat(), dst.stat()
    if ss.st_size != ds.st_size:
        return True
    return abs(ss.st_mtime - ds.st_mtime) > 1.0


def _sync_dir(dst_dir: Path, items: list[dict], label: str) -> tuple[int, int]:
    """Sync items into dst_dir: copy new/changed, remove stale.

    Only prunes top-level files whose extension matches one of the expected
    items.  Sub-directories (``_thumbs/``) and ``.html`` detail pages are
    preserved — they're managed by separate passes.
    """
    dst_dir.mkdir(parents=True, exist_ok=True)
    expected = {item["name"]: item for item in items}
    expected_exts = {Path(name).suffix.lower() for name in expected}

    removed = 0
    for existing in list(dst_dir.iterdir()):
        if not existing.is_file():
            continue
        if existing.name in expected:
            continue
        if existing.suffix.lower() not in expected_exts:
            continue  # leave unrelated artifacts alone (html, thumbs, …)
        existing.unlink()
        removed += 1

    copied = 0
    for name, item in expected.items():
        dst = dst_dir / name
        if _needs_copy(item["path"], dst):
            shutil.copy2(item["path"], dst)
            copied += 1

    if copied or removed:
        print(f"  {label}: {copied} copied, {removed} removed (of {len(items)})")
    else:
        print(f"  {label}: unchanged ({len(items)} files)")
    return copied, removed


def copy_to_worktree(
    wt_path: Path,
    reports: list[dict],
    figures: list[dict],
) -> list[str]:
    """Diff-copy reports/figures into the worktree.

    Only changed files are copied; files no longer in the include set are
    removed.  Git still sees only content deltas, but we avoid churning the
    filesystem on every deploy.
    """
    _sync_dir(wt_path / "reports", reports, "reports")
    _sync_dir(wt_path / "figures", figures, "figures")
    return [f"reports/{r['name']}" for r in reports] + [f"figures/{f['name']}" for f in figures]


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
# Tags, snapshots, backfill
# ---------------------------------------------------------------------------

def apply_tags(repo_root: Path, wt_path: Path, tag_name: str, tag_main: bool) -> None:
    """Tag the current gh-pages HEAD as site/<name> and (optionally) main HEAD
    as src/<name>.  Tags are pushed to origin.  Refuses to overwrite existing
    tags — delete and re-run if you need to move a tag."""
    site_tag = f"site/{tag_name}"
    src_tag = f"src/{tag_name}"

    # gh-pages side.
    proc = git("tag", site_tag, cwd=wt_path, check=False)
    if proc.returncode != 0:
        print(f"  WARNING: could not create {site_tag}: {proc.stderr.strip()}")
        print(f"    (delete with `git tag -d {site_tag}` or pick a different name)")
    else:
        print(f"  Tagged gh-pages as {site_tag}")
        git("push", "origin", site_tag, cwd=wt_path)

    # main side.
    if tag_main:
        proc = git("tag", src_tag, cwd=repo_root, check=False)
        if proc.returncode != 0:
            print(f"  WARNING: could not create {src_tag}: {proc.stderr.strip()}")
        else:
            head = git("rev-parse", "--abbrev-ref", "HEAD", cwd=repo_root).stdout.strip()
            print(f"  Tagged {head} HEAD as {src_tag}")
            git("push", "origin", src_tag, cwd=repo_root)


_DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


def _resolve_ref(repo_root: Path, ref: str, branch: str) -> str:
    """Resolve a user-supplied ref to a commit hash.

    If ref looks like YYYY-MM-DD, return the most recent branch commit on or
    before that date.  Otherwise return ref unchanged (git will validate)."""
    if _DATE_RE.match(ref):
        proc = git(
            "rev-list", "-1",
            f"--before={ref} 23:59:59",
            branch,
            cwd=repo_root,
            check=False,
        )
        resolved = proc.stdout.strip()
        if not resolved:
            print(f"ERROR: no {branch} commit found on or before {ref}.")
            sys.exit(1)
        print(f"Resolved {ref} -> {resolved[:8]}")
        return resolved
    return ref


def snapshot_at(repo_root: Path, ref: str, deploy_cfg: dict) -> Path:
    """Create a detached worktree at a past gh-pages state for local browsing."""
    snap_dir = deploy_cfg.get("snapshot_dir", "_gh-pages-snapshot")
    branch = deploy_cfg.get("branch", "gh-pages")
    snap_path = repo_root / snap_dir

    resolved = _resolve_ref(repo_root, ref, branch)

    # Remove any existing snapshot worktree at that path.
    existing = git("worktree", "list", "--porcelain", cwd=repo_root).stdout
    if str(snap_path) in existing:
        print(f"Removing existing snapshot at {snap_path} ...")
        git("worktree", "remove", "--force", str(snap_path), cwd=repo_root, check=False)
    elif snap_path.exists():
        print(f"ERROR: {snap_path} exists but is not a git worktree. Remove it manually.")
        sys.exit(1)

    print(f"Creating snapshot worktree at {snap_path} (ref: {resolved[:12]}) ...")
    git("worktree", "add", "--detach", str(snap_path), resolved, cwd=repo_root)

    proc = git("log", "-1", "--format=%h %ad %s", "--date=short", cwd=snap_path)
    print(f"Snapshot commit: {proc.stdout.strip()}")
    print(f"\nBrowse at: {snap_path}")
    print(f"Open: file://{snap_path}/index.html")
    print(f"Close later: python3 scripts/deploy_pages.py --close-snapshot")
    return snap_path


def close_snapshot(repo_root: Path, deploy_cfg: dict) -> None:
    snap_dir = deploy_cfg.get("snapshot_dir", "_gh-pages-snapshot")
    snap_path = repo_root / snap_dir
    if not snap_path.exists():
        print(f"No snapshot worktree at {snap_path}.")
        return
    print(f"Removing snapshot worktree at {snap_path} ...")
    git("worktree", "remove", "--force", str(snap_path), cwd=repo_root)
    print("Done.")


def backfill_tags(repo_root: Path, deploy_cfg: dict, push: bool = True) -> None:
    """Create site/YYYY-MM-DD[-N] tags for every existing gh-pages commit.

    For dates with a single deploy: site/YYYY-MM-DD.
    For dates with multiple deploys: site/YYYY-MM-DD-1, ...-2, ... in
    chronological order (oldest first).  Tags that already exist are skipped.
    """
    branch = deploy_cfg.get("branch", "gh-pages")
    proc = git(
        "log", branch,
        "--reverse",
        "--format=%H %ad",
        "--date=format:%Y-%m-%d",
        cwd=repo_root,
    )
    lines = [ln for ln in proc.stdout.strip().splitlines() if ln]
    if not lines:
        print(f"No commits on {branch}.")
        return

    by_date: dict[str, list[str]] = {}
    for line in lines:
        h, date = line.split(None, 1)
        by_date.setdefault(date, []).append(h)

    tags: list[tuple[str, str]] = []
    for date, hashes in by_date.items():
        if len(hashes) == 1:
            tags.append((f"site/{date}", hashes[0]))
        else:
            for i, h in enumerate(hashes, start=1):
                tags.append((f"site/{date}-{i}", h))

    existing = set(
        git("tag", "-l", "site/*", cwd=repo_root).stdout.strip().splitlines()
    )

    created: list[str] = []
    for tag_name, commit in tags:
        if tag_name in existing:
            continue
        proc = git("tag", tag_name, commit, cwd=repo_root, check=False)
        if proc.returncode == 0:
            created.append(tag_name)
        else:
            print(f"  WARNING: could not create {tag_name}: {proc.stderr.strip()}")

    print(f"Backfill: {len(created)} new tags created, "
          f"{len(tags) - len(created)} already existed.")

    if created and push:
        print(f"Pushing {len(created)} tag(s) to origin ...")
        # Push in batches to avoid huge command lines.
        batch_size = 50
        for i in range(0, len(created), batch_size):
            batch = created[i : i + batch_size]
            git("push", "origin", *batch, cwd=repo_root)
        print("Done.")


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
        help="Custom commit message (default: 'Update site: YYYY-MM-DD HH:MM TZ').",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Path to site_config.yaml (default: scripts/site_config.yaml relative to repo root).",
    )
    parser.add_argument(
        "--tag",
        metavar="NAME",
        default=None,
        help="After deploy, tag gh-pages as site/NAME and main HEAD as src/NAME; push both tags.",
    )
    parser.add_argument(
        "--no-tag-main",
        action="store_true",
        help="When --tag is given, skip the src/NAME tag on main.",
    )
    parser.add_argument(
        "--snapshot-at",
        metavar="REF",
        default=None,
        help="Materialize a past gh-pages snapshot in _gh-pages-snapshot/ and exit. "
             "REF may be a commit hash, tag (e.g. site/2026-04-15), or YYYY-MM-DD.",
    )
    parser.add_argument(
        "--close-snapshot",
        action="store_true",
        help="Remove the snapshot worktree and exit.",
    )
    parser.add_argument(
        "--backfill-tags",
        action="store_true",
        help="Create site/YYYY-MM-DD[-N] tags for every existing gh-pages commit, push, and exit.",
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

    deploy_cfg = cfg.get("deploy", {})

    # ----- Standalone subcommands (exit early) ----------------------------
    if args.backfill_tags:
        backfill_tags(repo_root, deploy_cfg)
        return
    if args.close_snapshot:
        close_snapshot(repo_root, deploy_cfg)
        return
    if args.snapshot_at:
        snapshot_at(repo_root, args.snapshot_at, deploy_cfg)
        return

    # ----- Compile .tex reports if needed -----------------------------------
    compiled = compile_tex_reports(repo_root, cfg)
    if compiled:
        print(f"Compiled {len(compiled)} report(s): {', '.join(compiled)}")

    # ----- Collect files --------------------------------------------------
    reports = collect_reports(repo_root, cfg)
    figures = collect_figures(repo_root, cfg)
    past_deploys = collect_past_deploys(repo_root)

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
        print(f"\nPast deploys found: {len(past_deploys)} site/* tag(s)")

        # Generate and show the index path that would be written.
        index_html = generate_index_html(reports, figures, past_deploys, cfg, repo_root)
        preview_path = repo_root / "_dry_run_index.html"
        preview_path.write_text(index_html)
        print(f"\nGenerated index.html preview: {preview_path}")
        return

    # ----- Ensure worktree ------------------------------------------------
    wt_dir = deploy_cfg.get("worktree_dir", "_gh-pages-worktree")
    branch = deploy_cfg.get("branch", "gh-pages")
    wt_path = ensure_worktree(repo_root, wt_dir, branch)

    # ----- Copy files -----------------------------------------------------
    deployed = copy_to_worktree(wt_path, reports, figures)

    # ----- Render PDF thumbnails (page 1 preview per report) --------------
    thumbs_dir = wt_path / "reports" / "_thumbs"
    generate_thumbnails(reports, thumbs_dir)

    # ----- Per-report detail pages ----------------------------------------
    site = cfg.get("site", {})
    repo_url = site.get("repo_url", "https://github.com/Shuonli/ppg12")
    n_details = generate_detail_pages(reports, past_deploys, wt_path, repo_root, repo_url)
    print(f"  detail pages: {n_details} written")

    # ----- Generate index.html --------------------------------------------
    index_html = generate_index_html(reports, figures, past_deploys, cfg, repo_root)
    (wt_path / "index.html").write_text(index_html)
    deployed.append("index.html")

    print(f"\nPrepared {len(deployed)} files ({total_mb:.2f} MB):")
    for d in deployed:
        print(f"  {d}")

    # ----- Commit and push ------------------------------------------------
    message = args.message or f"Update site: {datetime.now(_EASTERN).strftime('%Y-%m-%d %H:%M %Z')}"
    pushed = deploy(wt_path, message)

    # ----- Tagging (after deploy so tag points at the freshly pushed commit)
    if args.tag:
        apply_tags(repo_root, wt_path, args.tag, tag_main=not args.no_tag_main)

    if pushed:
        print(f"\nDeployed to https://shuonli.github.io/ppg12/")
    print("Done.")


if __name__ == "__main__":
    main()
