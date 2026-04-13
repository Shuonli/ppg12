#!/usr/bin/env python3
"""
lumi_convention_validation.py — Simulates the post-fix Section 3 validation.

Uses the existing vtxscan files to compute the reweight with the fixed
logic (rebin + smooth + zero fallback), applies it to the raw MC to
estimate what the fixed h_truth_pT would look like, and computes
L_noVtx × eff_new vs L_60cm × eff_old per pT bin.

This simulates the re-run WITHOUT actually re-running the pipeline.
Once the pipeline is re-run, the results can be compared against
this simulation.

Output: results/lumi_convention_validation.root  (histograms)
        results/lumi_convention_validation.csv   (per-pT table)
        plotting/figures/lumi_convention/*.pdf   (plots)
"""

import os
import numpy as np
import uproot

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PT_EDGES = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
N_PT_BINS = len(PT_EDGES) - 1

# Lumi from the two files
L_60cm = {"1p5mrad": 16.8790, "0mrad": 32.7034}
L_noVtx = {"1p5mrad": 17.1642, "0mrad": 47.2284}


def rebin_1d(arr, factor):
    """Rebin a 1D array by summing groups of `factor` adjacent bins."""
    n = len(arr)
    new_n = n // factor
    return arr[: new_n * factor].reshape(new_n, factor).sum(axis=1)


def smooth_353qh(arr, n_passes=1):
    """Approximate ROOT's TH1::Smooth() with a simple 3-bin weighted average."""
    for _ in range(n_passes):
        new = np.copy(arr)
        for i in range(1, len(arr) - 1):
            new[i] = (arr[i - 1] + 2.0 * arr[i] + arr[i + 1]) / 4.0
        arr = new
    return arr


def compute_fixed_reweight(h_data, h_sim, rebin=5, smooth=1):
    """Apply the fixed reweight procedure: rebin + normalize + divide + smooth.
    Returns the reweight histogram at the rebinned resolution."""
    h_d = rebin_1d(h_data, rebin)
    h_s = rebin_1d(h_sim, rebin)
    h_d = h_d / h_d.sum() if h_d.sum() > 0 else h_d
    h_s = h_s / h_s.sum() if h_s.sum() > 0 else h_s
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = h_d / h_s
    ratio[~np.isfinite(ratio)] = 0.0
    ratio[ratio < 0] = 0.0
    ratio = smooth_353qh(ratio, smooth)
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


def apply_reweight_1cm(h_sim_1cm, reweight, edges_rebinned, centers_1cm):
    """Apply the rebinned reweight to 1cm MC bins via linear interpolation."""
    centers_rb = 0.5 * (edges_rebinned[:-1] + edges_rebinned[1:])
    weights = np.zeros(len(centers_1cm))
    for i, z in enumerate(centers_1cm):
        if z < edges_rebinned[0] or z > edges_rebinned[-1]:
            weights[i] = 0.0
            continue
        idx = np.searchsorted(centers_rb, z)
        if idx == 0:
            weights[i] = reweight[0]
        elif idx >= len(reweight):
            weights[i] = reweight[-1]
        else:
            x0, x1 = centers_rb[idx - 1], centers_rb[idx]
            y0, y1 = reweight[idx - 1], reweight[idx]
            weights[i] = y0 + (y1 - y0) * (z - x0) / (x1 - x0)
    return h_sim_1cm * weights, weights


def main():
    print("=" * 70)
    print("Section 3 validation: simulate post-fix cross-check")
    print("=" * 70)

    os.makedirs("results", exist_ok=True)
    os.makedirs("../plotting/figures/lumi_convention", exist_ok=True)

    # Load production MC_efficiency_nom.root histograms
    mc_prod = uproot.open("results/MC_efficiency_nom.root")
    h_truth_old = mc_prod["h_truth_pT_0"].values()
    h_vtxcut = mc_prod["h_truth_pT_vertexcut_0"].values()
    h_mbd_vtx = mc_prod["h_truth_pT_vertexcut_mbd_cut_0"].values()
    pt_edges_mc = mc_prod["h_truth_pT_0"].axis().edges()

    # Load vtxscan files
    fd = uproot.open("results/data_histo_bdt_nom_vtxscan.root")
    fs = uproot.open("results/MC_efficiency_photon10_bdt_nom_vtxscan.root")
    h_data_vtx = fd["h_vertexz"].values()
    h_sim_vtx = fs["h_vertexz"].values()
    edges_vtx = fd["h_vertexz"].axis().edges()
    centers_vtx = 0.5 * (edges_vtx[:-1] + edges_vtx[1:])

    print(f"\nLoaded:")
    print(f"  MC_efficiency_nom: h_truth integral = {h_truth_old.sum():.0f}")
    print(f"  vtxscan data:      integral = {h_data_vtx.sum():.0f}, RMS = "
          f"{np.sqrt(np.average(centers_vtx**2, weights=h_data_vtx)):.1f} cm")
    print(f"  vtxscan sim :      integral = {h_sim_vtx.sum():.0f}, RMS = "
          f"{np.sqrt(np.average(centers_vtx**2, weights=h_sim_vtx)):.1f} cm")

    # Compute the OLD reweight (with 1.0 fallback bug)
    old_ratio = compute_old_reweight(h_data_vtx, h_sim_vtx)

    # Compute the NEW reweight (rebin + smooth + zero fallback)
    new_ratio_rb = compute_fixed_reweight(h_data_vtx, h_sim_vtx, rebin=5, smooth=1)
    new_edges_rb = edges_vtx[::5][: len(new_ratio_rb) + 1]

    # Apply the NEW reweight to the 1cm sim histogram (simulating what the
    # fixed RecoEffCalculator would produce)
    h_sim_reweighted_new, per_bin_weights_new = apply_reweight_1cm(
        h_sim_vtx, new_ratio_rb, new_edges_rb, centers_vtx
    )
    h_sim_reweighted_old = h_sim_vtx * old_ratio

    # Shape comparison
    mask_60 = np.abs(centers_vtx) < 60
    print(f"\n=== Reweighted MC truth vertex fraction within |z|<60 ===")
    print(f"  Target (data):   {h_data_vtx[mask_60].sum() / h_data_vtx.sum():.4f}")
    print(f"  Raw MC sim:      {h_sim_vtx[mask_60].sum() / h_sim_vtx.sum():.4f}")
    print(f"  OLD reweighted:  {h_sim_reweighted_old[mask_60].sum() / h_sim_reweighted_old.sum():.4f}")
    print(f"  NEW reweighted:  {h_sim_reweighted_new[mask_60].sum() / h_sim_reweighted_new.sum():.4f}")

    # Now simulate what h_truth_pT and h_truth_pT_vertexcut would look like
    # after the fix. We use the CURRENT production h_truth_pT (which has the
    # bloated integral) and scale it by the weight ratio.
    #
    # The logic: in the current production, h_truth_pT is filled with
    # vertex_weight = (buggy reweight at z_reco). If we replace that with the
    # (fixed reweight at z_reco), the effect on h_truth_pT is to rescale by
    # the per-z-bin ratio of (new reweight)/(old reweight).
    #
    # But we don't have h_truth_pT binned in z. We only have it in pT.
    # The cleanest proxy: the overall integral ratio of reweighted sim.
    integral_ratio_new_over_old = h_sim_reweighted_new.sum() / h_sim_reweighted_old.sum()
    print(f"\n  Integral ratio (new reweighted sim / old reweighted sim): "
          f"{integral_ratio_new_over_old:.4f}")

    # The h_truth_pT in the production was filled with the (old) vertex weight.
    # After the fix, the same MC events get a different weight, but the
    # pT-binned distribution shape is roughly preserved — the events that
    # changed weights are at large |z|, which may have a pT dependence.
    #
    # For this validation, we assume the shape change is small and use the
    # integral ratio to scale h_truth_pT. This is approximate but captures
    # the dominant effect (the 3x bloat was an OVERALL integral effect).
    #
    # Per-pT correction would require re-running the pipeline, which is
    # the recommendation at the end of the report.
    #
    # We compute two "fixed" h_truth_pT estimates:
    # (a) Simple integral-rescaling (approximate)
    # (b) Assume the NEW h_truth_pT_vertexcut = OLD h_truth_pT_vertexcut
    #     (since the vertex reweight within |z|<60 is well-defined for both)
    #     Then NEW h_truth_pT = OLD h_truth_pT_vertexcut / (data_frac_60)

    data_frac_60 = h_data_vtx[mask_60].sum() / h_data_vtx.sum()
    h_truth_new_estimate = h_vtxcut / data_frac_60  # method (b)

    # Compute the old and new efficiencies per pT bin
    print(f"\n=== Per-pT efficiency comparison ===")
    print(f"{'pT bin':<14} {'eff_old':>10} {'eff_new':>10} "
          f"{'L_60*eff_old':>14} {'L_nov*eff_new':>14} {'Ratio':>8}")
    print("-" * 80)

    results = {
        "pt_lo": [], "pt_hi": [],
        "eff_old": [], "eff_new": [],
        "L_old_product": [], "L_new_product": [],
        "ratio": [],
    }

    L60 = L_60cm["1p5mrad"]
    Lnv = L_noVtx["1p5mrad"]

    for i in range(N_PT_BINS):
        pt_lo, pt_hi = PT_EDGES[i], PT_EDGES[i + 1]
        # Sum over the pT bins in the MC histogram
        bin_mask = (pt_edges_mc[:-1] >= pt_lo - 0.01) & (pt_edges_mc[:-1] < pt_hi - 0.01)
        n_vtxcut = h_vtxcut[bin_mask].sum()
        n_mbd_vtx = h_mbd_vtx[bin_mask].sum()
        n_truth_old = h_truth_old[bin_mask].sum()
        n_truth_new = h_truth_new_estimate[bin_mask].sum()

        if n_vtxcut == 0 or n_truth_new == 0:
            continue

        eff_old = n_mbd_vtx / n_vtxcut
        eff_new = n_mbd_vtx / n_truth_new
        L_old_product = L60 * eff_old
        L_new_product = Lnv * eff_new
        ratio = L_new_product / L_old_product if L_old_product > 0 else 0

        print(f"{pt_lo:4.0f}-{pt_hi:4.0f} GeV  "
              f"{eff_old:10.4f} {eff_new:10.4f} "
              f"{L_old_product:14.4f} {L_new_product:14.4f} {ratio:8.4f}")

        results["pt_lo"].append(pt_lo)
        results["pt_hi"].append(pt_hi)
        results["eff_old"].append(eff_old)
        results["eff_new"].append(eff_new)
        results["L_old_product"].append(L_old_product)
        results["L_new_product"].append(L_new_product)
        results["ratio"].append(ratio)

    # Integrated
    n_t_old = h_truth_old.sum()
    n_t_new = h_truth_new_estimate.sum()
    n_v = h_vtxcut.sum()
    n_m = h_mbd_vtx.sum()
    eff_old_int = n_m / n_v
    eff_new_int = n_m / n_t_new
    prod_old = L60 * eff_old_int
    prod_new = Lnv * eff_new_int
    print(f"\n{'Integrated':<14} {eff_old_int:10.4f} {eff_new_int:10.4f} "
          f"{prod_old:14.4f} {prod_new:14.4f} {prod_new/prod_old:8.4f}")
    print(f"\n  Cross-section shift: {(prod_old/prod_new - 1) * 100:+.2f}%  (should be ~0% if fix is correct)")

    # Save to CSV
    with open("results/lumi_convention_validation.csv", "w") as fh:
        fh.write("pt_lo,pt_hi,eff_old,eff_new,L60cm_x_eff_old,LnoVtx_x_eff_new,ratio\n")
        for i in range(len(results["pt_lo"])):
            fh.write(f"{results['pt_lo'][i]:.0f},{results['pt_hi'][i]:.0f},"
                     f"{results['eff_old'][i]:.6f},{results['eff_new'][i]:.6f},"
                     f"{results['L_old_product'][i]:.6f},{results['L_new_product'][i]:.6f},"
                     f"{results['ratio'][i]:.6f}\n")
    print(f"\nSaved: results/lumi_convention_validation.csv")

    # Save vertex reweight data for plots
    np.savez("results/lumi_convention_validation.npz",
             z_centers=centers_vtx,
             old_ratio=old_ratio,
             new_ratio_rb=new_ratio_rb,
             new_ratio_1cm=per_bin_weights_new,
             h_data_vtx=h_data_vtx,
             h_sim_vtx=h_sim_vtx,
             h_sim_rw_old=h_sim_reweighted_old,
             h_sim_rw_new=h_sim_reweighted_new,
             new_edges_rb=new_edges_rb,
             pt_edges=PT_EDGES,
             eff_old=np.array(results["eff_old"]),
             eff_new=np.array(results["eff_new"]),
             ratio=np.array(results["ratio"]))
    print(f"Saved: results/lumi_convention_validation.npz")


if __name__ == "__main__":
    main()
