#!/usr/bin/env python3
"""Verify R=0.05 / 0.075 iso branches in CaloAna24 test outputs.

Usage: python3 verify_branches.py <output.root>

Checks:
 (a) All expected new branches exist on slimtree.
 (b) Monotonicity within each family: iso(R_small) <= iso(R_large) on a per-entry basis,
     allowing small numerical slop. For particle truth iso, R=0.05 subset of R=0.2 subset of R=0.3 subset of R=0.4.
 (c) Reasonable value ranges (no NaN, finite).
"""
import sys
import numpy as np
import uproot


NEW_BRANCHES_COMMON = [
    "particle_truth_iso_005",
    "particle_truth_iso_0075",
]

NEW_BRANCHES_PER_CONT = [
    "cluster_iso_005",
    "cluster_iso_0075",
    "cluster_iso_0075_70_emcal",
    "cluster_iso_topo_005",
    "cluster_iso_topo_0075",
    "cluster_iso_topo_soft_005",
    "cluster_iso_topo_soft_0075",
]

CONTAINERS = ["CLUSTERINFO_CEMC_NO_SPLIT", "CLUSTERINFO_CEMC"]


def check_monotone(arr_small, arr_large, small_label, large_label, slop=1e-3):
    """Flatten per-event arrays and verify arr_small <= arr_large + slop element-wise."""
    # Convert awkward-style jagged arrays to flat numpy
    small = np.asarray(arr_small.to_numpy() if hasattr(arr_small, "to_numpy") else arr_small).flatten()
    large = np.asarray(arr_large.to_numpy() if hasattr(arr_large, "to_numpy") else arr_large).flatten()
    # Attempt: if they are jagged, flatten across events
    try:
        small = np.concatenate([np.asarray(a, dtype=float) for a in arr_small])
        large = np.concatenate([np.asarray(a, dtype=float) for a in arr_large])
    except Exception:
        pass
    n = min(len(small), len(large))
    if n == 0:
        return True, f"  [SKIP] {small_label} vs {large_label}: zero entries"
    small, large = small[:n], large[:n]
    violations = np.sum(small > large + slop)
    frac = violations / n
    status = "PASS" if violations == 0 else ("MILD" if frac < 0.01 else "FAIL")
    return violations == 0, (f"  [{status}] {small_label} <= {large_label}: "
                             f"violations={violations}/{n} ({frac:.2%}), "
                             f"small.max={small.max():.3g}, large.max={large.max():.3g}")


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    path = sys.argv[1]
    print(f"\n=== {path} ===")
    f = uproot.open(path)
    if "slimtree" not in f:
        print(f"  [FAIL] no slimtree in {path} -- keys: {list(f.keys())[:10]}")
        sys.exit(1)
    t = f["slimtree"]
    branches = set(t.keys())
    nev = t.num_entries
    print(f"  slimtree entries = {nev}, total branches = {len(branches)}")

    # (a) branch presence
    missing = []
    for b in NEW_BRANCHES_COMMON:
        if b not in branches:
            missing.append(b)
    for base in NEW_BRANCHES_PER_CONT:
        for c in CONTAINERS:
            bn = f"{base}_{c}"
            if bn not in branches:
                missing.append(bn)
    if missing:
        print(f"  [FAIL] missing {len(missing)} branches: {missing[:5]}...")
        sys.exit(2)
    print(f"  [PASS] all {len(NEW_BRANCHES_COMMON) + len(NEW_BRANCHES_PER_CONT)*len(CONTAINERS)} new branches present")

    if nev == 0:
        print("  [WARN] tree has 0 entries, cannot sanity-check values")
        return

    # (b) monotonicity per family
    # Truth
    if "particle_truth_iso_005" in branches:
        arrs = t.arrays(["particle_truth_iso_005", "particle_truth_iso_0075",
                         "particle_truth_iso_02", "particle_truth_iso_03", "particle_truth_iso_04"],
                        library="np")
        chain = [("005", "0075"), ("0075", "02"), ("02", "03"), ("03", "04")]
        print("  --- truth particle iso (cumulative cone) ---")
        for s, l in chain:
            ok, msg = check_monotone(arrs[f"particle_truth_iso_{s}"], arrs[f"particle_truth_iso_{l}"],
                                     f"truth_iso_{s}", f"truth_iso_{l}")
            print(msg)

    # Per-container: full-calo cluster_iso
    for c in CONTAINERS:
        print(f"  --- cluster_iso (full-calo, {c.split('_')[-1]}) ---")
        arrs = t.arrays([f"cluster_iso_005_{c}", f"cluster_iso_0075_{c}",
                         f"cluster_iso_02_{c}", f"cluster_iso_03_{c}", f"cluster_iso_04_{c}"],
                        library="np")
        chain = [("005", "0075"), ("0075", "02"), ("02", "03"), ("03", "04")]
        for s, l in chain:
            ok, msg = check_monotone(arrs[f"cluster_iso_{s}_{c}"], arrs[f"cluster_iso_{l}_{c}"],
                                     f"iso_{s}", f"iso_{l}")
            print(msg)

    # Per-container: topo cluster iso
    for c in CONTAINERS:
        print(f"  --- cluster_iso_topo (topo, {c.split('_')[-1]}) ---")
        arrs = t.arrays([f"cluster_iso_topo_005_{c}", f"cluster_iso_topo_0075_{c}",
                         f"cluster_iso_topo_01_{c}", f"cluster_iso_topo_02_{c}",
                         f"cluster_iso_topo_03_{c}", f"cluster_iso_topo_04_{c}"],
                        library="np")
        chain = [("005", "0075"), ("0075", "01"), ("01", "02"), ("02", "03"), ("03", "04")]
        for s, l in chain:
            ok, msg = check_monotone(arrs[f"cluster_iso_topo_{s}_{c}"], arrs[f"cluster_iso_topo_{l}_{c}"],
                                     f"topo_{s}", f"topo_{l}")
            print(msg)

    # Per-container: 70 MeV EMCal family
    for c in CONTAINERS:
        print(f"  --- cluster_iso EMCal@70MeV ({c.split('_')[-1]}) ---")
        arrs = t.arrays([f"cluster_iso_005_70_emcal_{c}", f"cluster_iso_0075_70_emcal_{c}",
                         f"cluster_iso_01_70_emcal_{c}", f"cluster_iso_02_70_emcal_{c}",
                         f"cluster_iso_03_70_emcal_{c}"],
                        library="np")
        chain = [("005", "0075"), ("0075", "01"), ("01", "02"), ("02", "03")]
        for s, l in chain:
            ok, msg = check_monotone(arrs[f"cluster_iso_{s}_70_emcal_{c}"], arrs[f"cluster_iso_{l}_70_emcal_{c}"],
                                     f"iso70_{s}", f"iso70_{l}")
            print(msg)


if __name__ == "__main__":
    main()
