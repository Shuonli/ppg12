#!/usr/bin/env python3
"""Quick diagnostic: read first ~500k events of each photon sample,
compute the leading-truth-photon-pT distribution with cross-section
weight, and report photon10/photon5 and photon20/photon10 ratios.
If the cross-sections are correctly tuned, the ratios should be ~1.0
in regions where BOTH samples are fully efficient.
"""
import os, sys, numpy as np, uproot

XSEC = {"photon5": 146359.3, "photon10": 6944.675, "photon20": 130.4461}
SLIMTREE = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/{sample}/condorout/combined.root"
BIN_EDGES = np.arange(0.0, 50.001, 1.0)
ETA_MAX = 0.7
TARGET_EVENTS = 500_000  # subsample for speed

def fill_sample(sample):
    path = SLIMTREE.format(sample=sample)
    f = uproot.open(path)
    t = f["slimtree"]
    n_total = t.num_entries
    print(f"[{sample}] {n_total:,} events in slimtree, sampling first {TARGET_EVENTS:,}")
    H = np.zeros(len(BIN_EDGES) - 1, dtype=np.float64)
    Hsq = np.zeros_like(H)
    # Per-event weight uses TOTAL n_total (so weight is xsec/N_processed for the full sample)
    w_evt = XSEC[sample] / n_total
    n_seen = 0
    n_with_photon = 0
    for arr in t.iterate(["particle_pid","particle_Pt","particle_Eta","nparticles"],
                          step_size="200 MB", library="np", entry_stop=TARGET_EVENTS):
        pids = arr["particle_pid"]
        pts = arr["particle_Pt"]
        etas = arr["particle_Eta"]
        nev = len(arr["nparticles"])
        for i in range(nev):
            mask = (pids[i] == 22) & (np.abs(etas[i]) < ETA_MAX)
            if not mask.any(): continue
            n_with_photon += 1
            mx = float(pts[i][mask].max())
            bi = int(np.searchsorted(BIN_EDGES, mx, side="right") - 1)
            if 0 <= bi < len(H):
                H[bi] += w_evt
                Hsq[bi] += w_evt * w_evt
        n_seen += nev
        if n_seen >= TARGET_EVENTS: break
    print(f"[{sample}] processed {n_seen:,}, {n_with_photon:,} with truth photon in |eta|<0.7")
    return H, Hsq, n_seen

H = {}; HSQ = {}; N = {}
for s in ["photon5", "photon10", "photon20"]:
    H[s], HSQ[s], N[s] = fill_sample(s)

# Scale up to match if all events were processed (extrapolate from sample fraction)
for s in ["photon5", "photon10", "photon20"]:
    f_used = N[s] / uproot.open(SLIMTREE.format(sample=s))["slimtree"].num_entries
    print(f"[{s}] sample fraction used: {f_used:.4f}")

print()
print(f"{'pT':>5}  {'photon5':>11}  {'photon10':>11}  {'photon20':>11}  "
      f"{'10/5':>7}  {'20/10':>7}")
print("-" * 72)
centers = 0.5 * (BIN_EDGES[1:] + BIN_EDGES[:-1])
for i, c in enumerate(centers):
    if 8 <= c <= 36:
        h5 = H["photon5"][i]
        h10 = H["photon10"][i]
        h20 = H["photon20"][i]
        r10 = h10/h5 if h5>0 else 0
        r20 = h20/h10 if h10>0 else 0
        print(f"{c:5.1f}  {h5:11.3e}  {h10:11.3e}  {h20:11.3e}  {r10:7.3f}  {r20:7.3f}")
