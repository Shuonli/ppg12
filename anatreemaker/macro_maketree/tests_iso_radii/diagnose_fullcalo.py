#!/usr/bin/env python3
"""Diagnose the 1/7 full-calo monotonicity violation.

cluster_iso_005 / cluster_iso_0075 : manual sum of calculateET over 3 layers, minus cluster ET
cluster_iso_02 / cluster_iso_03 / cluster_iso_04 : RawCluster::get_et_iso builtin

These use different definitions of "cluster self-energy subtraction" and different
tower-level logic, so a small mismatch at R=0.2 between the two families is expected.
"""
import sys
import numpy as np
import uproot

path = sys.argv[1] if len(sys.argv) > 1 else "sim_photon10/caloana.root"
t = uproot.open(path)["slimtree"]

for c in ["CLUSTERINFO_CEMC_NO_SPLIT", "CLUSTERINFO_CEMC"]:
    print(f"\n=== {c} ===")
    arrs = t.arrays([f"cluster_Et_{c}",
                     f"cluster_iso_005_{c}", f"cluster_iso_0075_{c}",
                     f"cluster_iso_02_{c}", f"cluster_iso_03_{c}", f"cluster_iso_04_{c}"],
                    library="np")
    # Flatten across events
    Et = np.concatenate([np.asarray(a, dtype=float) for a in arrs[f"cluster_Et_{c}"]])
    i005 = np.concatenate([np.asarray(a, dtype=float) for a in arrs[f"cluster_iso_005_{c}"]])
    i0075 = np.concatenate([np.asarray(a, dtype=float) for a in arrs[f"cluster_iso_0075_{c}"]])
    i02 = np.concatenate([np.asarray(a, dtype=float) for a in arrs[f"cluster_iso_02_{c}"]])
    i03 = np.concatenate([np.asarray(a, dtype=float) for a in arrs[f"cluster_iso_03_{c}"]])
    i04 = np.concatenate([np.asarray(a, dtype=float) for a in arrs[f"cluster_iso_04_{c}"]])

    print(f"  nclusters = {len(Et)}")
    for i in range(len(Et)):
        print(f"  cluster[{i}] ET={Et[i]:.3f}  iso_005={i005[i]:.3f}  iso_0075={i0075[i]:.3f}"
              f"  iso_02={i02[i]:.3f}  iso_03={i03[i]:.3f}  iso_04={i04[i]:.3f}"
              + (" <-- 0075>02 mismatch" if i0075[i] > i02[i] + 1e-3 else ""))
