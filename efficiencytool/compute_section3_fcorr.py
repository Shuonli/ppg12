#!/usr/bin/env python3
"""
compute_section3_fcorr.py — Steps B3-B5 of the Section 3 luminosity correction.

Convolves the MC-derived MBD response function (from extract_mbd_response.py)
with the data z-vertex distribution (from extract_data_vertex.py) to compute
the per-pT, per-period f_corr correction factor.

Then adds the Poisson-weighted pileup correction using the double-interaction
response function (Pass₂).

This implements the formula from the sPHENIX luminosity note Section 3:
  f_corr = [1/(1-Prob0)] × [Prob1 × Pass1 + Prob2 × Pass2 + ...]

where:
  Pass_k(pT) = ∫ P_k(MBD+vtx | z) × f_data(z) dz
  P_1 = single-interaction response function (from MC)
  P_2 = double-interaction response function (from MC)

Finally, computes L_eff per period and validates against the old approach.

Output: results/section3_fcorr.root, results/section3_fcorr.csv

Usage:
    python3 compute_section3_fcorr.py
"""

import argparse
import os

import numpy as np
import uproot

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PT_EDGES = np.array([8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36], dtype=float)
N_PT_BINS = len(PT_EDGES) - 1
XING_RUN = 51274

# Event-level double-interaction fractions from calc_pileup_range.C
# P(n>=2|n>=1), lumi-weighted, NOT cluster-weighted
# mu is derived from these fractions by solving P(>=2|>=1) = 1 - mu*exp(-mu)/(1-exp(-mu))
EVENT_FRAC_DOUBLE = {
    "0mrad": 0.111,    # from .claude/rules/double-interaction.md
    "1p5mrad": 0.039,  # from .claude/rules/double-interaction.md
}


def mu_from_event_frac(f_double, tol=1e-8):
    """Solve P(>=2|>=1) = f_double for mu using bisection (robust).
    P(1|>=1) = mu*exp(-mu)/(1-exp(-mu)), so P(>=2|>=1) = 1 - P(1|>=1)."""
    from math import exp

    def p_ge2(mu):
        if mu < 1e-15:
            return 0.0
        e = exp(-mu)
        return 1.0 - mu * e / (1.0 - e)

    # Bisection on [0, 2]
    lo, hi = 0.001, 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if p_ge2(mid) < f_double:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def poisson_fracs(mu, n_max=5):
    """Compute P(n|n>=1) for n=1..n_max."""
    from math import exp, factorial
    p_ge1 = 1 - exp(-mu)
    fracs = []
    for n in range(1, n_max + 1):
        pn = (mu**n * exp(-mu) / factorial(n)) / p_ge1
        fracs.append(pn)
    return np.array(fracs)


def load_lumi_file(path):
    """Parse lumi file → {run: bit30corr}."""
    lumi = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 9:
                continue
            try:
                rn = int(parts[0])
            except ValueError:
                continue
            lumi[rn] = float(parts[7])  # Bit30Corr
    return lumi


def main():
    parser = argparse.ArgumentParser(description="Compute Section 3 f_corr")
    parser.add_argument("--results", default="results")
    parser.add_argument("--lumi-60cm",
                        default="/sphenix/user/shuhangli/ppg12/lumi/60cmLumi_fromJoey.list")
    parser.add_argument("--lumi-novtx",
                        default="/sphenix/user/shuhangli/ppg12/lumi/allzLumi_fromJoey.list")
    args = parser.parse_args()

    # --- Load MBD response functions (from B1) ---
    resp_file = os.path.join(args.results, "mbd_response_function.root")
    print(f"Loading MBD response: {resp_file}")
    f_resp = uproot.open(resp_file)

    h_resp_s = f_resp["h_response_single_2d"]
    response_single = h_resp_s.values()  # [z_bin, pt_bin]
    response_double = f_resp["h_response_double_2d"].values()
    z_edges = h_resp_s.axis(0).edges()  # z-axis edges from histogram
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    dz = z_edges[1] - z_edges[0]

    print(f"  Response shape: {response_single.shape} (z × pT)")
    print(f"  z range: [{z_edges[0]:.0f}, {z_edges[-1]:.0f}] cm, dz = {dz:.0f} cm")

    # --- Load data vertex distributions (from B2) ---
    vtx_file = os.path.join(args.results, "data_vertex_distributions.root")
    print(f"Loading data vertex: {vtx_file}")
    f_vtx = uproot.open(vtx_file)

    data_vtxz = {
        "0mrad": f_vtx["h_vtxz_period_0mrad"].values(),
        "1p5mrad": f_vtx["h_vtxz_period_1p5mrad"].values(),
        "all": f_vtx["h_vtxz_all"].values(),
    }

    # Normalize to probability distributions
    data_vtxz_norm = {}
    for key, h in data_vtxz.items():
        total = h.sum()
        data_vtxz_norm[key] = h / total if total > 0 else h
        print(f"  {key}: {total:.0f} events, RMS = {np.sqrt(np.average(z_centers**2, weights=h)):.1f} cm")

    # --- Load lumi files ---
    lumi_60cm = load_lumi_file(args.lumi_60cm)
    lumi_novtx = load_lumi_file(args.lumi_novtx)

    # Compute totals per period
    L = {}
    for period, (rmin, rmax) in [("0mrad", (47289, 51274)), ("1p5mrad", (51274, 54000))]:
        L[f"60cm_{period}"] = sum(v for r, v in lumi_60cm.items() if rmin <= r <= rmax and v > 0)
        L[f"novtx_{period}"] = sum(v for r, v in lumi_novtx.items() if rmin <= r <= rmax and v > 0)

    print(f"\nLuminosity totals:")
    for k, v in L.items():
        print(f"  {k}: {v:.4f} pb^-1")

    # --- Compute f_corr per pT bin by convolving response × data vtx ---
    print(f"\n{'='*80}")
    print("STEP B3: f_corr(pT) = ∫ P(MBD+vtx | z) × f_data(z) dz")
    print(f"{'='*80}")

    # For single interaction (Pass1): convolve response_single × data_vtxz
    # For double interaction (Pass2): convolve response_double × data_vtxz
    results = {}

    for period in ["0mrad", "1p5mrad"]:
        f_data = data_vtxz_norm[period]

        pass1 = np.zeros(N_PT_BINS)
        pass2 = np.zeros(N_PT_BINS)

        for ipt in range(N_PT_BINS):
            # Convolve: sum over z bins of [response(z, pT) × f_data(z)]
            pass1[ipt] = np.sum(response_single[:, ipt] * f_data)
            pass2[ipt] = np.sum(response_double[:, ipt] * f_data)

        results[period] = {"pass1": pass1, "pass2": pass2}

        print(f"\n--- {period} ---")
        print(f"{'pT bin':<14} {'Pass1(single)':>14} {'Pass2(double)':>14} {'Ratio P2/P1':>12}")
        print("-" * 56)
        for ipt in range(N_PT_BINS):
            ratio = pass2[ipt] / pass1[ipt] if pass1[ipt] > 0 else 0
            print(f"{PT_EDGES[ipt]:4.0f}-{PT_EDGES[ipt+1]:4.0f} GeV  "
                  f"{pass1[ipt]:14.6f} {pass2[ipt]:14.6f} {ratio:12.4f}")

    # --- Compute Poisson-weighted f_corr (Step B5) ---
    print(f"\n{'='*80}")
    print("STEP B5: Poisson-weighted f_corr = [1/(1-Prob0)] × Σ Prob_k × Pass_k")
    print(f"{'='*80}")

    for period in ["0mrad", "1p5mrad"]:
        f_d = EVENT_FRAC_DOUBLE[period]
        mu = mu_from_event_frac(f_d)
        print(f"  μ derived from f_double={f_d}: μ = {mu:.4f}")
        fracs = poisson_fracs(mu)  # P(n=1|n>=1), P(n=2|n>=1), ...
        pass1 = results[period]["pass1"]
        pass2 = results[period]["pass2"]

        # f_corr = P(1|>=1) × Pass1 + P(2|>=1) × Pass2 + P(>=3|>=1) × Pass2 (truncation)
        f_corr_pileup = fracs[0] * pass1 + (fracs[1] + sum(fracs[2:])) * pass2
        f_corr_single = pass1  # no pileup (just Pass1)

        results[period]["f_corr_single"] = f_corr_single
        results[period]["f_corr_pileup"] = f_corr_pileup
        results[period]["pileup_correction"] = np.where(
            f_corr_single > 0, f_corr_pileup / f_corr_single, 1.0)

        print(f"\n--- {period} (μ = {mu}) ---")
        print(f"  P(1|≥1) = {fracs[0]:.4f}, P(2|≥1) = {fracs[1]:.4f}, P(≥3|≥1) = {sum(fracs[2:]):.6f}")
        print(f"\n{'pT bin':<14} {'f_corr(no PU)':>14} {'f_corr(PU)':>14} {'PU correction':>14}")
        print("-" * 58)
        for ipt in range(N_PT_BINS):
            print(f"{PT_EDGES[ipt]:4.0f}-{PT_EDGES[ipt+1]:4.0f} GeV  "
                  f"{f_corr_single[ipt]:14.6f} {f_corr_pileup[ipt]:14.6f} "
                  f"{results[period]['pileup_correction'][ipt]:13.4f} "
                  f"({(results[period]['pileup_correction'][ipt]-1)*100:+.1f}%)")

    # --- Step B4: Compute L_eff and validate ---
    print(f"\n{'='*80}")
    print("STEP B4: VALIDATION — L(noVtx) × f_corr vs L(60cm) × vertex_eff")
    print(f"{'='*80}")

    # Load vertex_eff from production (MC_efficiency_nom.root)
    mc_file = os.path.join(args.results, "MC_efficiency_nom.root")
    try:
        f_mc = uproot.open(mc_file)
        h_vtxcut = f_mc["h_truth_pT_vertexcut_0"].values()
        h_mbd_vtx = f_mc["h_truth_pT_vertexcut_mbd_cut_0"].values()
        mc_edges = f_mc["h_truth_pT_vertexcut_0"].axis().edges()

        # Extract vertex_eff per analysis pT bin
        vtx_eff_old = np.zeros(N_PT_BINS)
        for ipt in range(N_PT_BINS):
            mask = (mc_edges[:-1] >= PT_EDGES[ipt] - 0.01) & (mc_edges[:-1] < PT_EDGES[ipt+1] - 0.01)
            d = h_vtxcut[mask].sum()
            n = h_mbd_vtx[mask].sum()
            vtx_eff_old[ipt] = n / d if d > 0 else 0
        has_old = True
    except Exception as e:
        print(f"  WARNING: Cannot load old vertex_eff from {mc_file}: {e}")
        vtx_eff_old = np.full(N_PT_BINS, 0.68)
        has_old = False

    for period in ["0mrad", "1p5mrad"]:
        l60 = L[f"60cm_{period}"]
        lnv = L[f"novtx_{period}"]
        f_corr_s = results[period]["f_corr_single"]
        f_corr_p = results[period]["f_corr_pileup"]

        L_eff_single = lnv * f_corr_s
        L_eff_pileup = lnv * f_corr_p
        L_old = l60 * vtx_eff_old

        results[period]["L_eff_single"] = L_eff_single
        results[period]["L_eff_pileup"] = L_eff_pileup
        results[period]["L_old"] = L_old

        label = "0 mrad" if period == "0mrad" else "1.5 mrad"
        print(f"\n--- {label} ---")
        print(f"  L(60cm) = {l60:.4f},  L(noVtx) = {lnv:.4f}")
        print(f"\n{'pT bin':<14} {'L×vtx_eff OLD':>14} {'L×f(no PU)NEW':>14} {'Ratio NEW/OLD':>14} {'L×f(PU) NEW':>14} {'Ratio PU/OLD':>14}")
        print("-" * 86)
        for ipt in range(N_PT_BINS):
            ratio_s = L_eff_single[ipt] / L_old[ipt] if L_old[ipt] > 0 else 0
            ratio_p = L_eff_pileup[ipt] / L_old[ipt] if L_old[ipt] > 0 else 0
            print(f"{PT_EDGES[ipt]:4.0f}-{PT_EDGES[ipt+1]:4.0f} GeV  "
                  f"{L_old[ipt]:14.4f} {L_eff_single[ipt]:14.4f} {ratio_s:14.4f} "
                  f"{L_eff_pileup[ipt]:14.4f} {ratio_p:14.4f}")

        # Integrated
        L_old_int = np.sum(L_old * (PT_EDGES[1:] - PT_EDGES[:-1])) / np.sum(PT_EDGES[1:] - PT_EDGES[:-1])
        L_eff_s_int = np.sum(L_eff_single * (PT_EDGES[1:] - PT_EDGES[:-1])) / np.sum(PT_EDGES[1:] - PT_EDGES[:-1])
        L_eff_p_int = np.sum(L_eff_pileup * (PT_EDGES[1:] - PT_EDGES[:-1])) / np.sum(PT_EDGES[1:] - PT_EDGES[:-1])
        print(f"\n  Integrated (pT-weighted avg):")
        print(f"    OLD:             {L_old_int:.4f}")
        print(f"    NEW (no pileup): {L_eff_s_int:.4f}  ratio = {L_eff_s_int/L_old_int:.4f}")
        print(f"    NEW (pileup):    {L_eff_p_int:.4f}  ratio = {L_eff_p_int/L_old_int:.4f}")
        print(f"    → Cross-section change (single): {(L_old_int/L_eff_s_int - 1)*100:+.2f}%")
        print(f"    → Cross-section change (pileup): {(L_old_int/L_eff_p_int - 1)*100:+.2f}%")

    # --- Save results ---
    outpath_root = os.path.join(args.results, "section3_fcorr.root")
    with uproot.recreate(outpath_root) as fout:
        for period in ["0mrad", "1p5mrad"]:
            for key in ["pass1", "pass2", "f_corr_single", "f_corr_pileup",
                        "pileup_correction", "L_eff_single", "L_eff_pileup", "L_old"]:
                fout[f"h_{key}_{period}"] = (results[period][key], PT_EDGES)
        fout["h_vtx_eff_old"] = (vtx_eff_old, PT_EDGES)
        # pt_edges embedded in histogram axes — no separate save

    # CSV
    outpath_csv = os.path.join(args.results, "section3_fcorr.csv")
    with open(outpath_csv, "w") as fh:
        header = "pt_lo,pt_hi"
        for period in ["0mrad", "1p5mrad"]:
            header += f",pass1_{period},pass2_{period}"
            header += f",fcorr_single_{period},fcorr_pileup_{period}"
            header += f",pileup_corr_{period}"
            header += f",Leff_single_{period},Leff_pileup_{period},Lold_{period}"
        header += ",vtx_eff_old"
        fh.write(header + "\n")

        for ipt in range(N_PT_BINS):
            row = f"{PT_EDGES[ipt]:.0f},{PT_EDGES[ipt+1]:.0f}"
            for period in ["0mrad", "1p5mrad"]:
                r = results[period]
                row += f",{r['pass1'][ipt]:.8f},{r['pass2'][ipt]:.8f}"
                row += f",{r['f_corr_single'][ipt]:.8f},{r['f_corr_pileup'][ipt]:.8f}"
                row += f",{r['pileup_correction'][ipt]:.6f}"
                row += f",{r['L_eff_single'][ipt]:.6f},{r['L_eff_pileup'][ipt]:.6f},{r['L_old'][ipt]:.6f}"
            row += f",{vtx_eff_old[ipt]:.6f}"
            fh.write(row + "\n")

    print(f"\nSaved: {outpath_root}")
    print(f"Saved: {outpath_csv}")


if __name__ == "__main__":
    main()
