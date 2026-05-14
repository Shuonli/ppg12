#!/usr/bin/env python3
"""
Audit MC merge files in efficiencytool/results/ for TFileMerger silent partial-merge bug.
For each merged file with both _0rad and _1p5mrad siblings, compare integral of every
TH1/TH2 key against sum(0rad + 1p5mrad). Broken: rel diff > 5%.
"""
import sys
import os
import glob
from datetime import date

import uproot
import numpy as np

RESULTS = sys.argv[1] if len(sys.argv) > 1 else "/sphenix/user/shuhangli/ppg12/efficiencytool/results"
TOL = 0.05

PATTERNS = ["MC_efficiency_bdt_*.root", "MC_efficiency_jet_bdt_*.root", "MC_response_bdt_*.root"]
EXCLUDE_SUFFIXES = ("_0rad.root", "_1p5mrad.root", "_0mrad.root", ".broken_TFileMerger", ".pre_hadd_backup")


def is_excluded(p):
    name = os.path.basename(p)
    return any(name.endswith(suf) for suf in EXCLUDE_SUFFIXES)


def integrate(arr):
    return float(np.sum(np.asarray(arr)))


def audit_pair(merged_path, p0_path, p1_path):
    """Return list of (key, i0, i1, im, rel) for broken keys.

    Detects three failure modes:
      1. Empty merged file (no keys, e.g. ~700 B header-only stub from a
         silent hadd failure). Flagged as a single synthetic 'EMPTY' entry.
      2. Stale merged file (mtime < max per-period mtime − 60 s). Flagged
         as a single synthetic 'STALE' entry — the merged file is older
         than its inputs, so even if the integrals happen to match the
         current per-period files they aren't from this run.
      3. Per-key relative-integral mismatch > TOL (the classic
         TFileMerger silent partial-merge bug).
    """
    broken = []

    # ---- Mode 1: empty/header-only merged file
    merged_size = os.path.getsize(merged_path)
    if merged_size < 2048:
        with uproot.open(merged_path) as fm:
            n_keys = len(fm.keys())
        if n_keys == 0:
            broken.append(("__EMPTY__", 0.0, 0.0, 0.0,
                           float(merged_size)))
            return broken

    # ---- Mode 2: stale merged file (older than per-period siblings)
    m_ts = os.path.getmtime(merged_path)
    p0_ts = os.path.getmtime(p0_path)
    p1_ts = os.path.getmtime(p1_path)
    max_period = max(p0_ts, p1_ts)
    if m_ts < max_period - 60:
        broken.append(("__STALE__",
                       m_ts, max_period, m_ts, max_period - m_ts))
        return broken

    # ---- Mode 3: per-key integral comparison (full union of keys so a key
    # missing from the merged file is visible, not silently dropped).
    with uproot.open(merged_path) as fm, uproot.open(p0_path) as f0, uproot.open(p1_path) as f1:
        m_keys = set(fm.keys())
        s_keys = set(f0.keys()) & set(f1.keys())
        union = m_keys | s_keys
        for k in union:
            try:
                if k in f0 and k in f1:
                    v0 = f0[k].values()
                    v1 = f1[k].values()
                    i0 = integrate(v0)
                    i1 = integrate(v1)
                else:
                    continue
                if k in fm:
                    om = fm[k]
                    if not hasattr(om, "values"):
                        continue
                    im = integrate(om.values())
                else:
                    # Key in both periods but missing from merged → break
                    im = 0.0
            except Exception:
                continue
            ssum = i0 + i1
            if ssum == 0:
                continue
            rel = abs(im - ssum) / abs(ssum)
            if rel > TOL:
                broken.append((k, i0, i1, im, rel))
    return broken


def main():
    candidates = []
    for pat in PATTERNS:
        for p in glob.glob(os.path.join(RESULTS, pat)):
            if is_excluded(p):
                continue
            candidates.append(p)
    candidates.sort()

    print(f"Auditing {len(candidates)} candidate merged files in {RESULTS}", flush=True)
    n_ok = 0
    n_skip = 0
    n_broken = 0
    broken_variants = {}

    for i, merged in enumerate(candidates):
        base = merged[:-len(".root")]
        p0 = base + "_0rad.root"
        p1 = base + "_1p5mrad.root"
        if not (os.path.exists(p0) and os.path.exists(p1)):
            n_skip += 1
            continue
        try:
            br = audit_pair(merged, p0, p1)
        except Exception as e:
            print(f"  ERROR on {os.path.basename(merged)}: {e}", flush=True)
            continue
        if br:
            n_broken += 1
            broken_variants[merged] = br
            print(f"  [{i+1}/{len(candidates)}] BROKEN: {os.path.basename(merged)}: {len(br)} keys", flush=True)
        else:
            n_ok += 1
        if (i+1) % 20 == 0:
            print(f"  ... {i+1}/{len(candidates)} processed (ok={n_ok}, broken={n_broken}, skip={n_skip})", flush=True)

    distinct_keys = set()
    for items in broken_variants.values():
        for k, *_ in items:
            distinct_keys.add(k)

    today = date.today().isoformat()
    report_path = f"/sphenix/user/shuhangli/ppg12/reports/merge_audit_{today}.md"
    os.makedirs(os.path.dirname(report_path), exist_ok=True)

    with open(report_path, "w") as r:
        r.write(f"# Merge audit -- {today}\n\n")
        r.write(f"Audited **{n_ok + n_broken}** variants ({n_skip} skipped: standalone). "
                f"**{n_broken}** broken ({len(distinct_keys)} distinct keys affected).\n\n")
        if n_broken == 0:
            r.write("**No breakages found.**\n")
        else:
            r.write("## Broken variants\n\n")
            for merged, items in sorted(broken_variants.items()):
                r.write(f"### `{os.path.basename(merged)}`\n\n")
                r.write("| key | 0rad | 1p5mrad | merged | sum | rel diff |\n")
                r.write("|---|---:|---:|---:|---:|---:|\n")
                for k, i0, i1, im, rel in sorted(items, key=lambda x: -x[4]):
                    ssum = i0 + i1
                    r.write(f"| `{k}` | {i0:.4e} | {i1:.4e} | {im:.4e} | {ssum:.4e} | {rel*100:.2f}% |\n")
                r.write("\n")
            r.write("## Remediation\n\n```bash\ncd /sphenix/user/shuhangli/ppg12/efficiencytool\n")
            for merged in sorted(broken_variants.keys()):
                base = merged[:-len(".root")]
                rel_merged = os.path.relpath(merged, "/sphenix/user/shuhangli/ppg12/efficiencytool")
                rel_p0 = os.path.relpath(base + "_0rad.root", "/sphenix/user/shuhangli/ppg12/efficiencytool")
                rel_p1 = os.path.relpath(base + "_1p5mrad.root", "/sphenix/user/shuhangli/ppg12/efficiencytool")
                r.write(f"hadd -f {rel_merged} {rel_p0} {rel_p1}\n")
            r.write("```\n")

    print(f"\nAudit complete: {n_ok} ok, {n_broken} broken, {n_skip} standalone-skip.", flush=True)
    print(f"Report: {report_path}", flush=True)
    if n_broken > 0:
        print("\nProposed hadd commands:", flush=True)
        for merged in sorted(broken_variants.keys()):
            base = merged[:-len(".root")]
            rel_merged = os.path.relpath(merged, "/sphenix/user/shuhangli/ppg12/efficiencytool")
            rel_p0 = os.path.relpath(base + "_0rad.root", "/sphenix/user/shuhangli/ppg12/efficiencytool")
            rel_p1 = os.path.relpath(base + "_1p5mrad.root", "/sphenix/user/shuhangli/ppg12/efficiencytool")
            print(f"  hadd -f {rel_merged} {rel_p0} {rel_p1}", flush=True)

    prev_audits = sorted(glob.glob("/sphenix/user/shuhangli/ppg12/reports/merge_audit_*.md"))
    prev_audits = [p for p in prev_audits if os.path.realpath(p) != os.path.realpath(report_path)]
    if prev_audits:
        prev = prev_audits[-1]
        print(f"\nPrior audit: {prev}", flush=True)
        with open(prev) as f:
            prev_text = f.read()
        prev_broken = set()
        for ln in prev_text.split("\n"):
            if ln.startswith("### "):
                prev_broken.add(ln.strip("# `"))
        new_broken = set(os.path.basename(p) for p in broken_variants.keys())
        added = new_broken - prev_broken
        resolved = prev_broken - new_broken
        if added:
            print(f"  newly BROKEN: {sorted(added)}", flush=True)
        if resolved:
            print(f"  RESOLVED:     {sorted(resolved)}", flush=True)
        if not added and not resolved:
            print("  no diff vs prior audit", flush=True)
    else:
        print("\nNo prior audit -- diff skipped.", flush=True)

    return n_broken


if __name__ == "__main__":
    sys.exit(0 if main() == 0 else 1)
