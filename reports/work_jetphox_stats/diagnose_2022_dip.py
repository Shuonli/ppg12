#!/usr/bin/env python3
"""
Diagnose the apparent dip at pT [20, 22] GeV in JETPHOX NLO spectra.

Tests:
  1. Fine-binned pT histogram (1 GeV bins) from 8 to 40 GeV for each of
     the 6 CT14nlo pawres files (direct + frag, 3 scales). See whether the
     dip is in the raw ntuple data or only after binning into PPG12 2-GeV bins.
  2. Break out direct-only vs frag-only contributions in PPG12 bins.
  3. Test hypothesis: is it localized to one scale, or present in all 3?
"""
import os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAWRES = "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres"
OUT_DIR = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/work_jetphox_stats"
NSEG = 100
DETA = 1.4

# PPG12 truth binning
PT_BINS_PPG12 = [7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45]

# Fine binning
PT_FINE_EDGES = np.arange(7, 41, 1.0)   # 1 GeV bins

def read_and_fill(path, edges):
    """Read a JETPHOX pawres file, fill a TH1 with photon pT (index 0) under |y|<0.7.
       Return (values, errors) arrays in pb/GeV/unit_y (divide by deta here).
    """
    with uproot.open(path) as f:
        t2 = f["t2"]
        # UserInfo — Python uproot cannot read TList[TVectorF] easily, pull via hex of header
        # Hardcode from our earlier verification
        u = t2["energy"].array(library="np", entry_stop=None)  # trigger preload? no — stream
    # actually stream with uproot
    with uproot.open(path) as f:
        t2 = f["t2"]
        step = 2_000_000
        hist = np.zeros(len(edges) - 1)
        sumw2 = np.zeros(len(edges) - 1)
        # Use awkward arrays to handle variable-length photon/jet record.
        import awkward as ak
        for arrs in t2.iterate(
            ["energy", "px", "py", "pz", "pdf_weight"],
            step_size=step,
            library="ak",
        ):
            # First element per event (= photon)
            px = ak.to_numpy(arrs["px"][:, 0])
            py = ak.to_numpy(arrs["py"][:, 0])
            pz = ak.to_numpy(arrs["pz"][:, 0])
            e  = ak.to_numpy(arrs["energy"][:, 0])
            w  = ak.to_numpy(arrs["pdf_weight"][:, 0]).astype(np.float64)
            pt = np.sqrt(px * px + py * py)
            safe = (e + pz > 0) & (e - pz > 0)
            y = np.where(safe, 0.5 * np.log((e + pz) / np.clip(e - pz, 1e-9, None)), 999.0)
            mask = (np.abs(y) < 0.7)
            h, _ = np.histogram(pt[mask], bins=edges, weights=w[mask])
            h2, _ = np.histogram(pt[mask], bins=edges, weights=w[mask] ** 2)
            hist += h
            sumw2 += h2
    return hist, np.sqrt(sumw2)

def normalize(hist, err, edges, xsec, nb_evt=1_000_000, nseg=NSEG):
    """Apply norma = xsec/nb_evt/nseg, divide by bin width, divide by deta."""
    norma = xsec / nb_evt / nseg
    bw = np.diff(edges)
    return hist * norma / bw / DETA, err * norma / bw / DETA

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Hardcode xsec from UserInfo (verified in /tmp/verify_pawres output)
    # CT14nlo
    samples = {
        "d_10": (f"{PAWRES}/ggdrhic_nlo_10.root", 78979.5),
        "o_10": (f"{PAWRES}/ggorhic_nlo_10.root", 21562.8),
        "d_05": (f"{PAWRES}/ggdrhic_nlo_05.root", 122029.0),
        "o_05": (f"{PAWRES}/ggorhic_nlo_05.root", 35526.2),
        "d_20": (f"{PAWRES}/ggdrhic_nlo_20.root", 54139.3),
        "o_20": (f"{PAWRES}/ggorhic_nlo_20.root", 14332.1),
    }

    # Fine binning
    print("=== Fine-bin (1 GeV) direct-only + fragmentation-only, mu=1.0 ===")
    fine_d, _ = read_and_fill(samples["d_10"][0], PT_FINE_EDGES)
    fine_o, _ = read_and_fill(samples["o_10"][0], PT_FINE_EDGES)
    fine_d_norm, _ = normalize(fine_d, np.zeros_like(fine_d), PT_FINE_EDGES, samples["d_10"][1])
    fine_o_norm, _ = normalize(fine_o, np.zeros_like(fine_o), PT_FINE_EDGES, samples["o_10"][1])
    centers = 0.5 * (PT_FINE_EDGES[:-1] + PT_FINE_EDGES[1:])
    print(f"{'pT':>6} {'direct':>12} {'frag':>12} {'total':>12}")
    for i, c in enumerate(centers):
        print(f"{c:6.1f} {fine_d_norm[i]:12.4g} {fine_o_norm[i]:12.4g} {fine_d_norm[i]+fine_o_norm[i]:12.4g}")

    # PPG12 bins direct-only + frag-only
    print("\n=== PPG12 bins direct-only + frag-only, 3 scales ===")
    bins = np.array(PT_BINS_PPG12, dtype=float)
    for scale in ["10", "05", "20"]:
        hd, _ = read_and_fill(samples[f"d_{scale}"][0], bins)
        ho, _ = read_and_fill(samples[f"o_{scale}"][0], bins)
        nd, _ = normalize(hd, np.zeros_like(hd), bins, samples[f"d_{scale}"][1])
        no_, _ = normalize(ho, np.zeros_like(ho), bins, samples[f"o_{scale}"][1])
        print(f"\n--- scale {scale} ---")
        print(f"{'pT bin':>12} {'direct':>10} {'frag':>10} {'total':>10}")
        for i in range(len(bins) - 1):
            print(f"{bins[i]:4.0f}-{bins[i+1]:4.0f} {nd[i]:10.4g} {no_[i]:10.4g} {nd[i]+no_[i]:10.4g}")

    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                                   gridspec_kw={"height_ratios": [3, 2]})
    # Fine binning (mu=1 direct+frag total)
    fine_total = fine_d_norm + fine_o_norm
    ax1.step(centers, fine_total, where="mid", color="C0", label="fine 1-GeV bins (d+o)")
    ax1.step(centers, fine_d_norm, where="mid", color="C1", ls="--", label="direct")
    ax1.step(centers, fine_o_norm, where="mid", color="C2", ls="--", label="frag")

    # PPG12 bins (same scale, direct+frag sum)
    ppg_h, _ = read_and_fill(samples["d_10"][0], bins)
    ppg_ho, _ = read_and_fill(samples["o_10"][0], bins)
    n_d, _ = normalize(ppg_h, np.zeros_like(ppg_h), bins, samples["d_10"][1])
    n_o, _ = normalize(ppg_ho, np.zeros_like(ppg_ho), bins, samples["o_10"][1])
    ppg_centers = 0.5 * (bins[:-1] + bins[1:])
    ax1.plot(ppg_centers, n_d + n_o, "o", color="red", label="PPG12 2-GeV bins (d+o)")
    ax1.set_yscale("log")
    ax1.set_ylabel(r"d$^2\sigma$/d$p_T$d$y$ (pb/GeV)")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_title("CT14nlo, mu=1.0, |y|<0.7 — fine vs PPG12 binning")

    # Residual from smooth power law
    # Fit a power law to the fine-binned total in the range 10 < pT < 35 GeV
    mask_fit = (centers > 10) & (centers < 35) & (fine_total > 0)
    slope, intercept = np.polyfit(np.log(centers[mask_fit]), np.log(fine_total[mask_fit]), 1)
    smooth = np.exp(intercept) * centers ** slope
    ax2.plot(centers, fine_total / smooth, "C0-", label="fine / power-law")
    ax2.plot(ppg_centers, (n_d + n_o) / (np.exp(intercept) * ppg_centers ** slope), "ro", label="PPG12 / power-law")
    ax2.axhline(1.0, color="k", lw=0.5)
    ax2.axvspan(20, 22, color="red", alpha=0.2, label="[20,22] bin")
    ax2.set_xlabel(r"$p_T$ (GeV)")
    ax2.set_ylabel("data / pT^{%.2f} power law" % slope)
    ax2.set_ylim(0.5, 1.8)
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    path = f"{OUT_DIR}/diagnose_2022_dip.pdf"
    plt.savefig(path)
    plt.savefig(path.replace(".pdf", ".png"), dpi=150)
    print(f"\nPlot: {path}")

if __name__ == "__main__":
    main()
