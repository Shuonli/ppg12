#!/usr/bin/env python3
"""
lumi_convention_validation.py — validates the post-fix L(noVtx) × eff_new
cross-check against the old L(60cm) × eff_old for both crossing-angle
periods (1.5 mrad and 0 mrad).

Uses existing first-pass vtxscan files to apply the fixed reweight
(rebin + smooth + zero fallback) and propagate through the efficiency
calculation, without re-running MergeSim.

Output: results/lumi_convention_validation.csv
        results/lumi_convention_validation.npz
"""

import os
import numpy as np
import uproot

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PT_EDGES = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
N_PT_BINS = len(PT_EDGES) - 1

# Luminosity totals from the two lumi files (summed per period)
L_60cm = {"1p5mrad": 16.8790, "0mrad": 32.7034}
L_noVtx = {"1p5mrad": 17.1642, "0mrad": 47.2284}

# Per-period vtxscan files + production MC file
PERIODS = {
    "1p5mrad": {
        "label": "1.5 mrad",
        "data_vtxscan": "results/data_histo_bdt_nom_vtxscan.root",
        "sim_vtxscan": "results/MC_efficiency_photon10_bdt_nom_vtxscan.root",
        "mc_prod": "results/MC_efficiency_nom.root",
    },
    "0mrad": {
        "label": "0 mrad",
        "data_vtxscan": "results/data_histo_bdt_0rad_vtxscan.root",
        "sim_vtxscan": "results/MC_efficiency_photon10_bdt_0rad_vtxscan.root",
        "mc_prod": "results/MC_efficiency_bdt_0rad.root",
    },
}


# ---------------------------------------------------------------------------
# Reweight helpers
# ---------------------------------------------------------------------------
def rebin_1d(arr, factor):
    n = len(arr)
    new_n = n // factor
    return arr[: new_n * factor].reshape(new_n, factor).sum(axis=1)


def smooth_3bin(arr, n_passes=1):
    for _ in range(n_passes):
        new = np.copy(arr)
        for i in range(1, len(arr) - 1):
            new[i] = (arr[i - 1] + 2.0 * arr[i] + arr[i + 1]) / 4.0
        arr = new
    return arr


def compute_new_reweight(h_data, h_sim, rebin=5, smooth=1):
    """Apply the fixed reweight logic: rebin + normalize + divide + smooth."""
    h_d = rebin_1d(h_data, rebin)
    h_s = rebin_1d(h_sim, rebin)
    h_d = h_d / h_d.sum() if h_d.sum() > 0 else h_d
    h_s = h_s / h_s.sum() if h_s.sum() > 0 else h_s
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = h_d / h_s
    ratio[~np.isfinite(ratio)] = 0.0
    ratio[ratio < 0] = 0.0
    ratio = smooth_3bin(ratio, smooth)
    return ratio


def compute_old_reweight(h_data, h_sim):
    """Old reweight: 1cm bins, fallback=1.0 for non-finite/zero bins."""
    h_d = h_data / h_data.sum() if h_data.sum() > 0 else h_data
    h_s = h_sim / h_sim.sum() if h_sim.sum() > 0 else h_sim
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = h_d / h_s
    ratio[~np.isfinite(ratio)] = 1.0
    ratio[ratio <= 0] = 1.0
    return ratio


def apply_reweight_1cm(h_sim_1cm, reweight_rb, edges_rb, centers_1cm):
    """Linearly interpolate rebinned reweight to 1 cm bins and apply."""
    centers_rb = 0.5 * (edges_rb[:-1] + edges_rb[1:])
    weights = np.zeros(len(centers_1cm))
    for i, z in enumerate(centers_1cm):
        if z < edges_rb[0] or z > edges_rb[-1]:
            weights[i] = 0.0
            continue
        idx = np.searchsorted(centers_rb, z)
        if idx == 0:
            weights[i] = reweight_rb[0]
        elif idx >= len(reweight_rb):
            weights[i] = reweight_rb[-1]
        else:
            x0, x1 = centers_rb[idx - 1], centers_rb[idx]
            y0, y1 = reweight_rb[idx - 1], reweight_rb[idx]
            weights[i] = y0 + (y1 - y0) * (z - x0) / (x1 - x0)
    return h_sim_1cm * weights, weights


# ---------------------------------------------------------------------------
# Per-period validation
# ---------------------------------------------------------------------------
def validate_period(period_key, cfg):
    print("\n" + "=" * 72)
    print(f"VALIDATION: {cfg['label']}")
    print("=" * 72)

    # Load vtxscan files
    fd = uproot.open(cfg["data_vtxscan"])
    fs = uproot.open(cfg["sim_vtxscan"])
    h_data = fd["h_vertexz"].values()
    h_sim = fs["h_vertexz"].values()
    edges = fd["h_vertexz"].axis().edges()
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Load production MC histograms
    try:
        fp = uproot.open(cfg["mc_prod"])
        h_vtxcut = fp["h_truth_pT_vertexcut_0"].values()
        h_mbd_vtx = fp["h_truth_pT_vertexcut_mbd_cut_0"].values()
        pt_edges_mc = fp["h_truth_pT_vertexcut_0"].axis().edges()
    except Exception as e:
        print(f"  WARNING: cannot load {cfg['mc_prod']}: {e}")
        return None

    print(f"\nInputs:")
    print(f"  Data vtxscan: {cfg['data_vtxscan']}")
    print(f"  Sim  vtxscan: {cfg['sim_vtxscan']}")
    print(f"  MC prod:      {cfg['mc_prod']}")
    print(f"  Data vtxscan integral:  {h_data.sum():.0f}, RMS = "
          f"{np.sqrt(np.average(centers**2, weights=h_data)):.1f} cm")
    print(f"  Sim vtxscan integral:   {h_sim.sum():.0f}, RMS = "
          f"{np.sqrt(np.average(centers**2, weights=h_sim)):.1f} cm")

    # Compute reweights
    old_ratio = compute_old_reweight(h_data, h_sim)
    new_ratio_rb = compute_new_reweight(h_data, h_sim, rebin=5, smooth=1)
    new_edges_rb = edges[::5][: len(new_ratio_rb) + 1]

    # Apply to MC
    h_sim_rw_new, per_bin_w_new = apply_reweight_1cm(h_sim, new_ratio_rb, new_edges_rb, centers)
    h_sim_rw_old = h_sim * old_ratio

    # Fraction within |z|<60
    mask_60 = np.abs(centers) < 60
    data_frac_60 = h_data[mask_60].sum() / h_data.sum()
    old_frac_60 = h_sim_rw_old[mask_60].sum() / h_sim_rw_old.sum()
    new_frac_60 = h_sim_rw_new[mask_60].sum() / h_sim_rw_new.sum()
    raw_frac_60 = h_sim[mask_60].sum() / h_sim.sum()

    print(f"\n=== Reweighted MC fraction within |z|<60 ===")
    print(f"  Data (target):    {data_frac_60:.4f}")
    print(f"  Raw MC:           {raw_frac_60:.4f}")
    print(f"  OLD reweighted:   {old_frac_60:.4f}")
    print(f"  NEW reweighted:   {new_frac_60:.4f}")

    # Simulate post-fix h_truth_pT by scaling h_truth_pT_vertexcut by 1/data_frac_60.
    # This assumes the fixed reweight brings the MC truth vertex distribution
    # to match data so that h_vtxcut / h_truth_new = data_frac_60.
    h_truth_new = h_vtxcut / data_frac_60

    # Per-pT cross-check
    print(f"\n=== Per-pT L × eff comparison ===")
    L60 = L_60cm[period_key]
    Lnv = L_noVtx[period_key]
    print(f"  L(60cm)={L60:.4f} pb^-1, L(noVtx)={Lnv:.4f} pb^-1, ratio={Lnv/L60:.4f}")
    print(f"\n  {'pT bin':<14} {'eff_old':>10} {'eff_new':>10} "
          f"{'L60*eOld':>12} {'Lnov*eNew':>12} {'Ratio':>8}")
    print("  " + "-" * 68)

    pt_lo_list, pt_hi_list = [], []
    eff_old_list, eff_new_list = [], []
    Lold_list, Lnew_list, ratio_list = [], [], []

    for i in range(N_PT_BINS):
        pt_lo, pt_hi = PT_EDGES[i], PT_EDGES[i + 1]
        bin_mask = (pt_edges_mc[:-1] >= pt_lo - 0.01) & (pt_edges_mc[:-1] < pt_hi - 0.01)
        n_vtxcut = h_vtxcut[bin_mask].sum()
        n_mbd_vtx = h_mbd_vtx[bin_mask].sum()
        n_truth_new = h_truth_new[bin_mask].sum()

        if n_vtxcut == 0 or n_truth_new == 0:
            continue

        eff_old = n_mbd_vtx / n_vtxcut
        eff_new = n_mbd_vtx / n_truth_new
        L_old_product = L60 * eff_old
        L_new_product = Lnv * eff_new
        ratio = L_new_product / L_old_product if L_old_product > 0 else 0

        print(f"  {pt_lo:4.0f}-{pt_hi:4.0f} GeV  "
              f"{eff_old:10.4f} {eff_new:10.4f} "
              f"{L_old_product:12.4f} {L_new_product:12.4f} {ratio:8.4f}")

        pt_lo_list.append(pt_lo)
        pt_hi_list.append(pt_hi)
        eff_old_list.append(eff_old)
        eff_new_list.append(eff_new)
        Lold_list.append(L_old_product)
        Lnew_list.append(L_new_product)
        ratio_list.append(ratio)

    # Integrated
    n_t_new = h_truth_new.sum()
    n_v = h_vtxcut.sum()
    n_m = h_mbd_vtx.sum()
    int_ratio = 0.0
    int_eff_old = 0.0
    int_eff_new = 0.0
    if n_v > 0 and n_t_new > 0:
        int_eff_old = n_m / n_v
        int_eff_new = n_m / n_t_new
        prod_old = L60 * int_eff_old
        prod_new = Lnv * int_eff_new
        int_ratio = prod_new / prod_old
        print(f"\n  {'Integrated':<14} {int_eff_old:10.4f} {int_eff_new:10.4f} "
              f"{prod_old:12.4f} {prod_new:12.4f} {int_ratio:8.4f}")
        print(f"  Cross-section shift: {(1/int_ratio - 1) * 100:+.2f}%")

    return {
        "label": cfg["label"],
        "z_centers": centers,
        "h_data_vtx": h_data,
        "h_sim_vtx": h_sim,
        "h_sim_rw_old": h_sim_rw_old,
        "h_sim_rw_new": h_sim_rw_new,
        "old_ratio": old_ratio,
        "new_ratio_1cm": per_bin_w_new,
        "new_ratio_rb": new_ratio_rb,
        "new_edges_rb": new_edges_rb,
        "data_frac_60": data_frac_60,
        "raw_frac_60": raw_frac_60,
        "old_frac_60": old_frac_60,
        "new_frac_60": new_frac_60,
        "pt_lo": np.array(pt_lo_list),
        "pt_hi": np.array(pt_hi_list),
        "eff_old": np.array(eff_old_list),
        "eff_new": np.array(eff_new_list),
        "L_old_product": np.array(Lold_list),
        "L_new_product": np.array(Lnew_list),
        "ratio": np.array(ratio_list),
        "integrated_ratio": int_ratio,
        "integrated_eff_old": int_eff_old,
        "integrated_eff_new": int_eff_new,
    }


def main():
    os.makedirs("results", exist_ok=True)

    all_results = {}
    for period_key, cfg in PERIODS.items():
        all_results[period_key] = validate_period(period_key, cfg)

    # Save CSV
    csv_path = "results/lumi_convention_validation.csv"
    with open(csv_path, "w") as fh:
        fh.write("period,pt_lo,pt_hi,eff_old,eff_new,L60cm_x_eff_old,"
                 "LnoVtx_x_eff_new,ratio\n")
        for period_key, r in all_results.items():
            if r is None:
                continue
            for i in range(len(r["pt_lo"])):
                fh.write(f"{period_key},{r['pt_lo'][i]:.0f},{r['pt_hi'][i]:.0f},"
                         f"{r['eff_old'][i]:.6f},{r['eff_new'][i]:.6f},"
                         f"{r['L_old_product'][i]:.6f},{r['L_new_product'][i]:.6f},"
                         f"{r['ratio'][i]:.6f}\n")
    print(f"\nSaved: {csv_path}")

    # Save NPZ for plotting
    npz_data = {}
    for period_key, r in all_results.items():
        if r is None:
            continue
        prefix = f"{period_key}_"
        for key in ["z_centers", "old_ratio", "new_ratio_1cm", "new_ratio_rb",
                    "h_data_vtx", "h_sim_vtx", "h_sim_rw_old", "h_sim_rw_new",
                    "new_edges_rb", "pt_lo", "pt_hi", "eff_old", "eff_new",
                    "L_old_product", "L_new_product", "ratio"]:
            npz_data[prefix + key] = np.asarray(r[key])
        for key in ["data_frac_60", "raw_frac_60", "old_frac_60", "new_frac_60",
                    "integrated_ratio", "integrated_eff_old", "integrated_eff_new"]:
            npz_data[prefix + key] = np.asarray(r[key])
    npz_data["pt_edges"] = PT_EDGES
    np.savez("results/lumi_convention_validation.npz", **npz_data)
    print("Saved: results/lumi_convention_validation.npz")


if __name__ == "__main__":
    main()
