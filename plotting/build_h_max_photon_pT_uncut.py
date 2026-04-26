#!/usr/bin/env python3
"""
Build the unrestricted leading-truth-photon-pT histogram for each
photon sample (photon5, photon10, photon20) by reading the slimtree
input directly with vectorized awkward operations.

Bypasses the per-sample event cut applied in
RecoEffCalculator_TTreeReader.C so the resulting distributions can
be ratio-divided to demonstrate where each sample becomes fully
efficient (analysis-note Fig 3) and stitched cleanly for the
Hagedorn fit (Fig 4).

Output: plotting/photon_max_pT_uncut.root
"""
import os, sys, time
import numpy as np
import awkward as ak
import uproot

XSEC = {"photon5": 146359.3, "photon10": 6944.675, "photon20": 130.4461}
SLIMTREE = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/{sample}/condorout/combined.root"
OUT = "/sphenix/user/shuhangli/ppg12/plotting/photon_max_pT_uncut.root"
BIN_EDGES = np.arange(0.0, 50.001, 0.5)  # 0.5 GeV bins for cleaner ratios
ETA_MAX = 0.7
N_PER_SAMPLE = 100_000_000  # use all events (~10M per file)

def fill_one(sample):
    path = SLIMTREE.format(sample=sample)
    f = uproot.open(path)
    t = f["slimtree"]
    n_total = t.num_entries
    n_use = min(n_total, N_PER_SAMPLE)
    print(f"[{sample}] {n_total:,} total, processing {n_use:,}", flush=True)
    H = np.zeros(len(BIN_EDGES) - 1, dtype=np.float64)
    Hsq = np.zeros_like(H)
    w_evt = XSEC[sample] / n_total
    n_seen = 0
    t0 = time.time()
    for arr in t.iterate(["particle_pid","particle_Pt","particle_Eta"],
                          step_size="500 MB", entry_stop=n_use):
        pids, pts, etas = arr["particle_pid"], arr["particle_Pt"], arr["particle_Eta"]
        # Vectorized leading photon: mask photons in |eta|<0.7, then take max per event
        mask = (pids == 22) & (abs(etas) < ETA_MAX)
        sel_pt = pts[mask]
        # events with at least one good photon
        n_good = ak.num(sel_pt)
        has_photon = n_good > 0
        max_pt = ak.to_numpy(ak.fill_none(ak.max(sel_pt, axis=1), 0.0))
        max_pt_good = max_pt[ak.to_numpy(has_photon)]
        # bin
        H_chunk, _ = np.histogram(max_pt_good, bins=BIN_EDGES)
        H += H_chunk * w_evt
        Hsq += H_chunk * (w_evt ** 2)
        n_seen += len(pids)
        if n_seen >= n_use:
            break
        elapsed = time.time() - t0
        rate = n_seen / elapsed
        eta_sec = (n_use - n_seen) / rate if rate > 0 else 0
        print(f"[{sample}]   {n_seen:,} ({100*n_seen/n_use:.0f}%)  "
              f"rate={rate:.0f} ev/s  ETA {eta_sec:.0f}s", flush=True)
    print(f"[{sample}] done in {time.time()-t0:.0f}s, integral={H.sum():.3e}", flush=True)
    return H, Hsq

def main():
    out = uproot.recreate(OUT)
    H, HSQ = {}, {}
    for s in ("photon5","photon10","photon20"):
        H[s], HSQ[s] = fill_one(s)
        out[f"h_max_photon_pT_{s}"] = (H[s], BIN_EDGES)
        # Also save sum-of-weights-squared so the consumer macro can
        # set the proper bin errors (sqrt(Hsq)) instead of the default
        # sqrt(values) that uproot writes when given a (vals, edges)
        # tuple.
        out[f"h_max_photon_pT_{s}_sumw2"] = (HSQ[s], BIN_EDGES)
    centers = 0.5*(BIN_EDGES[1:]+BIN_EDGES[:-1])
    print(f"\n{'pT':>5}  {'photon5':>11}  {'photon10':>11}  {'photon20':>11}  "
          f"{'10/5':>7}  {'20/10':>7}")
    for i,c in enumerate(centers):
        if 9 <= c <= 36 and (H['photon5'][i]>0 or H['photon10'][i]>0 or H['photon20'][i]>0):
            r10 = H['photon10'][i]/H['photon5'][i] if H['photon5'][i]>0 else 0
            r20 = H['photon20'][i]/H['photon10'][i] if H['photon10'][i]>0 else 0
            print(f"{c:5.1f}  {H['photon5'][i]:11.3e}  {H['photon10'][i]:11.3e}  "
                  f"{H['photon20'][i]:11.3e}  {r10:7.3f}  {r20:7.3f}")
    print(f"\nWrote {OUT}")

if __name__ == "__main__":
    sys.exit(main())
