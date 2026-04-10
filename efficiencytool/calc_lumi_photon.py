#!/usr/bin/env python3
"""
calc_lumi_photon.py — Per-run photon luminosity from slimtree scalers.

Follows the sPHENIX luminosity note:
  L_unc = N_live_mbd / (P_trig × sigma_MBD)
  L_corr = L_unc × F_corr

Trigger 30 only (active in 0 mrad runs, 47289-51274).
MBD N&S>=1 trigger = bit 10.
sigma_MBD = 25.2 mb (Vernier scan).

Outputs:
  results/lumi_photon_perrun.csv    — per-run luminosity table
  results/lumi_photon_qa.root       — TGraphErrors for run-by-run QA
  plotting/figures/lumi_photon_qa.pdf — QA plots
"""

import argparse
import glob
import os
import sys
from collections import defaultdict

import numpy as np
import yaml

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SIGMA_MBD_PB = 25.2e9          # 25.2 mb in pb
MBD_BIT = 10                    # MBD N&S>=1 trigger bit
TRIG_BIT = 30                   # analysis trigger
XING_RUN = 51274                # 0 mrad / 1.5 mrad boundary

# MBD single-collision firing probabilities (GEANT4, from lumi note Table 1)
P_BOTH = 0.562
P_N_ONLY = 0.175
P_S_ONLY = 0.175
P_NEITHER = 0.088
# Derived: P(not fire N) = P_S_ONLY + P_NEITHER = 0.263
# Derived: P(not fire S) = P_N_ONLY + P_NEITHER = 0.263
P_NOT_N = P_S_ONLY + P_NEITHER  # 0.263
P_NOT_S = P_N_ONLY + P_NEITHER  # 0.263

# ---------------------------------------------------------------------------
# Pileup correction
# ---------------------------------------------------------------------------

def p_mbd_fires(mu):
    """P(MBD N&S>=1 fires | mean interactions = mu).

    P(MBD|n) = 1 - 2*P_NOT_N^n + P_NEITHER^n  (inclusion-exclusion)
    P(fires) = sum_{n>=1} P(MBD|n) * Poisson(n, mu)
             = 1 - 2*exp(-mu*(1-P_NOT_N)) + exp(-mu*(1-P_NEITHER)) - (1-2+1)*exp(-mu)
             = 1 - 2*exp(-0.737*mu) + exp(-0.912*mu)
    """
    return 1.0 - 2.0 * np.exp(-0.737 * mu) + np.exp(-0.912 * mu)


def mu_from_mbd_fraction(f_obs, mu_max=5.0):
    """Invert p_mbd_fires(mu) = f_obs to find mu (Newton's method)."""
    if f_obs <= 0:
        return 0.0
    # initial guess: f ~ mu * P_BOTH for small mu
    mu = f_obs / P_BOTH
    for _ in range(50):
        val = p_mbd_fires(mu) - f_obs
        # derivative: dp/dmu = 2*0.737*exp(-0.737*mu) - 0.912*exp(-0.912*mu)
        deriv = 2.0 * 0.737 * np.exp(-0.737 * mu) - 0.912 * np.exp(-0.912 * mu)
        if abs(deriv) < 1e-15:
            break
        mu_new = mu - val / deriv
        mu_new = max(0, min(mu_new, mu_max))
        if abs(mu_new - mu) < 1e-10:
            mu = mu_new
            break
        mu = mu_new
    return mu


def f_corr(mu):
    """Pileup correction factor: <n | n>=1> = mu / (1 - exp(-mu))."""
    if mu < 1e-10:
        return 1.0
    return mu / (1.0 - np.exp(-mu))


# ---------------------------------------------------------------------------
# Load official luminosity file for comparison
# ---------------------------------------------------------------------------

def load_official_lumi(path):
    """Parse 60cmLumi_fromJoey.list → {run: (bit30corr, bit30uc)}"""
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
            lumi[rn] = (float(parts[7]), float(parts[8]))  # Bit30Corr, Bit30UC
    return lumi


# ---------------------------------------------------------------------------
# Pass 1: per-run scaler accumulation
# ---------------------------------------------------------------------------

def pass_scalers(data_glob, treename):
    """Read all data files and accumulate per-run scaler max/min."""
    import uproot

    files = sorted(glob.glob(data_glob))
    if not files:
        print(f"ERROR: no files matching {data_glob}")
        sys.exit(1)
    print(f"Processing {len(files)} data files for scaler extraction...")

    # per-run accumulators: max and min of scaler arrays
    run_live_max = defaultdict(lambda: np.full(64, -1, dtype=np.int64))
    run_live_min = defaultdict(lambda: np.full(64, np.iinfo(np.int64).max, dtype=np.int64))
    run_scaled_max = defaultdict(lambda: np.full(64, -1, dtype=np.int64))
    run_scaled_min = defaultdict(lambda: np.full(64, np.iinfo(np.int64).max, dtype=np.int64))
    run_nevents = defaultdict(int)

    branches = ["runnumber", "currentscaler_live", "currentscaler_scaled"]

    for fi, fname in enumerate(files):
        try:
            f = uproot.open(fname)
            tree = f[treename]
        except Exception as e:
            print(f"  SKIP {fname}: {e}")
            continue

        for batch in tree.iterate(branches, step_size=200000, library="np"):
            rns = batch["runnumber"]
            live = batch["currentscaler_live"]
            scaled = batch["currentscaler_scaled"]

            for rn in np.unique(rns):
                mask = rns == rn
                n = int(np.sum(mask))
                run_nevents[rn] += n

                live_sub = live[mask]
                scaled_sub = scaled[mask]

                mx = live_sub.max(axis=0)
                mn = live_sub.min(axis=0)
                run_live_max[rn] = np.maximum(run_live_max[rn], mx)
                run_live_min[rn] = np.minimum(run_live_min[rn], mn)

                mx2 = scaled_sub.max(axis=0)
                mn2 = scaled_sub.min(axis=0)
                run_scaled_max[rn] = np.maximum(run_scaled_max[rn], mx2)
                run_scaled_min[rn] = np.minimum(run_scaled_min[rn], mn2)

        if (fi + 1) % 10 == 0 or fi == len(files) - 1:
            print(f"  [{fi+1}/{len(files)}] {len(run_nevents)} runs so far")

    # compute ranges
    results = {}
    for rn in sorted(run_nevents.keys()):
        live_range = run_live_max[rn] - run_live_min[rn]
        scaled_range = run_scaled_max[rn] - run_scaled_min[rn]
        results[rn] = {
            "nevents": run_nevents[rn],
            "live_range": live_range,
            "scaled_range": scaled_range,
        }
    return results


# ---------------------------------------------------------------------------
# Luminosity computation
# ---------------------------------------------------------------------------

def compute_luminosity(scaler_results, official_lumi):
    """Compute per-run luminosity from scaler data."""
    lumi_table = {}

    for rn in sorted(scaler_results.keys()):
        info = scaler_results[rn]
        lr = info["live_range"]
        sr = info["scaled_range"]

        n_live_mbd = int(lr[MBD_BIT])
        n_live_trig = int(lr[TRIG_BIT])
        n_scaled_trig = int(sr[TRIG_BIT])

        # skip runs with no MBD counts
        if n_live_mbd <= 0:
            continue

        # prescale for trigger 30
        if n_scaled_trig > 0:
            prescale = n_live_trig / n_scaled_trig
            trig_active = True
        else:
            prescale = -1.0
            trig_active = False

        # uncorrected luminosity
        if trig_active and prescale > 0:
            l_unc = n_live_mbd / (prescale * SIGMA_MBD_PB)
        else:
            l_unc = 0.0

        # pileup correction using official file ratio (if available)
        f_corr_val = 1.0
        if rn in official_lumi:
            oc, ou = official_lumi[rn]
            if ou > 0:
                f_corr_val = oc / ou

        l_corr = l_unc * f_corr_val

        # official values for comparison
        off_corr, off_unc = official_lumi.get(rn, (0.0, 0.0))

        is_0mrad = rn < XING_RUN

        lumi_table[rn] = {
            "nevents": info["nevents"],
            "n_live_mbd": n_live_mbd,
            "n_live_trig": n_live_trig,
            "n_scaled_trig": n_scaled_trig,
            "prescale": prescale,
            "trig_active": trig_active,
            "l_unc": l_unc,
            "f_corr": f_corr_val,
            "l_corr": l_corr,
            "off_corr": off_corr,
            "off_unc": off_unc,
            "is_0mrad": is_0mrad,
        }

    return lumi_table


# ---------------------------------------------------------------------------
# Pass 2: per-run photon counting
# ---------------------------------------------------------------------------

def pass_photons(data_glob, treename, config, lumi_table):
    """Count photon candidates per run (trigger 30 only).
    Uses awkward arrays for vectorized jagged-array operations."""
    import uproot
    import awkward as ak

    # cuts from config
    eta_lo, eta_hi = config["analysis"]["eta_bins"]
    vtx_cut = config["analysis"]["vertex_cut"]
    node = config["input"]["cluster_node_name"]
    et_min = 7.0   # matching Cluster_rbr.C
    et_max = 40.0

    files = sorted(glob.glob(data_glob))
    print(f"\nProcessing {len(files)} data files for photon counting...")

    # branch names include cluster node suffix
    br_et = f"cluster_Et_{node}"
    br_eta = f"cluster_Eta_{node}"

    branches = ["runnumber", "scaledtrigger", "vertexz", br_et, br_eta]

    photon_counts = defaultdict(int)
    event_counts = defaultdict(int)

    for fi, fname in enumerate(files):
        try:
            f = uproot.open(fname)
            tree = f[treename]
        except Exception:
            continue

        for batch in tree.iterate(branches, step_size=200000):
            rns = ak.to_numpy(batch["runnumber"])
            trig = ak.to_numpy(batch["scaledtrigger"])
            vtxz = ak.to_numpy(batch["vertexz"])
            cet = batch[br_et]    # awkward jagged array
            ceta = batch[br_eta]  # awkward jagged array

            # event-level mask: trigger 30 + vertex cut
            evt_mask = trig[:, TRIG_BIT].astype(bool) & (np.abs(vtxz) < vtx_cut)

            # cluster-level mask (vectorized over jagged arrays)
            cl_mask = (cet > et_min) & (cet < et_max) & (abs(ceta) < abs(eta_hi))
            # per-event photon count (sum over clusters per event)
            n_phot_per_evt = ak.to_numpy(ak.sum(cl_mask, axis=1))

            # accumulate per-run, applying event mask
            for rn in np.unique(rns):
                rm = (rns == rn) & evt_mask
                event_counts[rn] += int(np.sum(rm))
                photon_counts[rn] += int(n_phot_per_evt[rm].sum())

        if (fi + 1) % 10 == 0 or fi == len(files) - 1:
            print(f"  [{fi+1}/{len(files)}] photon counting...")

    return dict(photon_counts), dict(event_counts)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def save_csv(lumi_table, photon_counts, event_counts, outpath):
    """Save per-run results to CSV."""
    with open(outpath, "w") as fh:
        fh.write("run,period,nevents,nevents_trig30,n_photon,"
                 "n_live_mbd,prescale_30,trig_active,"
                 "L_unc_calc,F_corr,L_corr_calc,L_corr_official,L_unc_official,"
                 "ratio_corr,nphoton_per_lumi\n")
        for rn in sorted(lumi_table.keys()):
            d = lumi_table[rn]
            n_evt_trig = event_counts.get(rn, 0)
            n_phot = photon_counts.get(rn, 0)
            ratio = d["l_corr"] / d["off_corr"] if d["off_corr"] > 0 else 0.0
            npl = n_phot / d["l_corr"] if d["l_corr"] > 0 else 0.0
            period = "0mrad" if d["is_0mrad"] else "1.5mrad"
            fh.write(f"{rn},{period},{d['nevents']},{n_evt_trig},{n_phot},"
                     f"{d['n_live_mbd']},{d['prescale']:.4f},{d['trig_active']},"
                     f"{d['l_unc']:.10f},{d['f_corr']:.6f},{d['l_corr']:.10f},"
                     f"{d['off_corr']:.10f},{d['off_unc']:.10f},"
                     f"{ratio:.6f},{npl:.2f}\n")
    print(f"Saved {outpath}")


def make_plots(lumi_table, photon_counts, event_counts, figdir):
    """Generate QA plots matching the lumi note style."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(figdir, exist_ok=True)

    # Collect arrays for 0 mrad runs where trigger 30 is active
    runs_0 = []
    lumi_calc = []
    lumi_off = []
    ratio_lumi = []
    nphoton_per_lumi = []
    fcorr_vals = []

    for rn in sorted(lumi_table.keys()):
        d = lumi_table[rn]
        if not d["trig_active"]:
            continue
        if d["l_corr"] <= 0 or d["off_corr"] <= 0:
            continue
        n_phot = photon_counts.get(rn, 0)
        if n_phot < 10:  # skip very low-stats runs
            continue

        runs_0.append(rn)
        lumi_calc.append(d["l_corr"])
        lumi_off.append(d["off_corr"])
        ratio_lumi.append(d["l_corr"] / d["off_corr"])
        npl = n_phot / d["l_corr"]
        nphoton_per_lumi.append(npl)
        fcorr_vals.append(d["f_corr"])

    runs_0 = np.array(runs_0)
    lumi_calc = np.array(lumi_calc)
    lumi_off = np.array(lumi_off)
    ratio_lumi = np.array(ratio_lumi)
    nphoton_per_lumi = np.array(nphoton_per_lumi)
    fcorr_vals = np.array(fcorr_vals)

    n_phot_arr = np.array([photon_counts.get(rn, 0) for rn in runs_0])
    npl_err = np.sqrt(n_phot_arr) / lumi_calc

    print(f"\n=== QA Summary (trigger 30, 0 mrad runs) ===")
    print(f"Runs with trigger 30 active: {len(runs_0)}")
    print(f"Run range: {runs_0.min()} - {runs_0.max()}")
    print(f"Total calculated lumi: {lumi_calc.sum():.4f} pb^-1")
    print(f"Total official lumi:   {lumi_off.sum():.4f} pb^-1")
    print(f"Mean L_calc/L_off ratio: {ratio_lumi.mean():.4f} ± {ratio_lumi.std():.4f}")
    print(f"Mean F_corr: {fcorr_vals.mean():.4f}")
    print(f"Mean Nphoton/Lint: {nphoton_per_lumi.mean():.1f}")

    # --- Figure 1: Nphoton / Lint vs run number ---
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.errorbar(runs_0, nphoton_per_lumi, yerr=npl_err, fmt=".", ms=3, lw=0.5, color="black")
    mean_npl = nphoton_per_lumi.mean()
    ax.axhline(mean_npl, color="red", lw=1.5, label=f"mean = {mean_npl:.0f}")
    ax.axvline(XING_RUN, color="blue", ls="--", lw=0.8, alpha=0.5, label="0→1.5 mrad boundary")
    ax.set_xlabel("Run Number")
    ax.set_ylabel(r"$N_{\gamma\mathrm{-cand}} / \mathcal{L}_{\mathrm{int}}$ [pb]")
    ax.set_title(r"Photon QA: $N_{\gamma}/\mathcal{L}$ vs Run (trigger 30, 0 mrad)")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "lumi_photon_nphot_vs_run.pdf"))
    print(f"Saved {figdir}/lumi_photon_nphot_vs_run.pdf")
    plt.close()

    # --- Figure 2: Calculated vs Official luminosity ratio ---
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(runs_0, ratio_lumi, ".", ms=3, color="black")
    ax.axhline(1.0, color="red", lw=1, ls="--")
    ax.axhline(ratio_lumi.mean(), color="blue", lw=1.5,
               label=f"mean = {ratio_lumi.mean():.4f}")
    ax.set_xlabel("Run Number")
    ax.set_ylabel(r"$\mathcal{L}_{\mathrm{calc}} / \mathcal{L}_{\mathrm{official}}$")
    ax.set_title("Luminosity cross-check: calculated (scaler) vs official (Bit30Corr)")
    ax.set_ylim(ratio_lumi.mean() - 0.15, ratio_lumi.mean() + 0.15)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "lumi_photon_ratio_vs_run.pdf"))
    print(f"Saved {figdir}/lumi_photon_ratio_vs_run.pdf")
    plt.close()

    # --- Figure 3: F_corr vs run number ---
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(runs_0, fcorr_vals, ".", ms=3, color="black")
    ax.set_xlabel("Run Number")
    ax.set_ylabel(r"$F_{\mathrm{corr}}$ (pileup correction)")
    ax.set_title("Per-run pileup correction factor (Bit30Corr / Bit30UC)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "lumi_photon_fcorr_vs_run.pdf"))
    print(f"Saved {figdir}/lumi_photon_fcorr_vs_run.pdf")
    plt.close()

    # --- Figure 4: 1D distribution of Nphoton/Lint ---
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(nphoton_per_lumi, bins=40, color="gray", edgecolor="black", lw=0.5,
            density=True, label=f"All 0 mrad (mean={mean_npl:.0f})")
    ax.axvline(mean_npl, color="red", lw=1.5)
    ax.set_xlabel(r"$N_{\gamma\mathrm{-cand}} / \mathcal{L}_{\mathrm{int}}$ [pb]")
    ax.set_ylabel("Normalized counts")
    ax.set_title(r"Distribution of $N_{\gamma}/\mathcal{L}$ across runs")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "lumi_photon_nphot_dist.pdf"))
    print(f"Saved {figdir}/lumi_photon_nphot_dist.pdf")
    plt.close()

    # --- Figure 5: Per-run luminosity ---
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(runs_0, lumi_calc * 1e3, ".", ms=3, color="black", label="Calculated")
    ax.plot(runs_0, lumi_off * 1e3, ".", ms=2, color="red", alpha=0.5, label="Official")
    ax.set_xlabel("Run Number")
    ax.set_ylabel(r"$\mathcal{L}_{\mathrm{int}}$ [nb$^{-1}$]")
    ax.set_title("Per-run integrated luminosity (trigger 30)")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "lumi_photon_perrun.pdf"))
    print(f"Saved {figdir}/lumi_photon_perrun.pdf")
    plt.close()


def save_root(lumi_table, photon_counts, event_counts, outpath):
    """Save TGraphErrors to ROOT file."""
    try:
        import ROOT
    except ImportError:
        print("ROOT not available; skipping ROOT output")
        return

    runs_0 = []
    npl, npl_err = [], []
    lc, lo, ratio = [], [], []

    for rn in sorted(lumi_table.keys()):
        d = lumi_table[rn]
        if not d["trig_active"] or d["l_corr"] <= 0 or d["off_corr"] <= 0:
            continue
        n_phot = photon_counts.get(rn, 0)
        if n_phot < 10:
            continue
        runs_0.append(float(rn))
        v = n_phot / d["l_corr"]
        e = np.sqrt(n_phot) / d["l_corr"]
        npl.append(v)
        npl_err.append(e)
        lc.append(d["l_corr"])
        lo.append(d["off_corr"])
        ratio.append(d["l_corr"] / d["off_corr"])

    n = len(runs_0)
    if n == 0:
        print("No valid runs for ROOT output")
        return

    x = np.array(runs_0)
    xe = np.zeros(n)
    fout = ROOT.TFile(outpath, "RECREATE")

    gr_npl = ROOT.TGraphErrors(n, x, np.array(npl), xe, np.array(npl_err))
    gr_npl.SetName("gr_nphoton_per_lumi")
    gr_npl.SetTitle("N_{#gamma-cand} / L_{int} vs Run;Run Number;N_{#gamma}/L [pb]")
    gr_npl.Write()

    gr_ratio = ROOT.TGraphErrors(n, x, np.array(ratio), xe, xe)
    gr_ratio.SetName("gr_lumi_ratio")
    gr_ratio.SetTitle("L_{calc}/L_{official} vs Run;Run Number;Ratio")
    gr_ratio.Write()

    gr_lcalc = ROOT.TGraphErrors(n, x, np.array(lc), xe, xe)
    gr_lcalc.SetName("gr_lumi_calc")
    gr_lcalc.SetTitle("Calculated luminosity;Run Number;L [pb^{-1}]")
    gr_lcalc.Write()

    gr_loff = ROOT.TGraphErrors(n, x, np.array(lo), xe, xe)
    gr_loff.SetName("gr_lumi_official")
    gr_loff.SetTitle("Official luminosity;Run Number;L [pb^{-1}]")
    gr_loff.Write()

    fout.Close()
    print(f"Saved {outpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Per-run photon luminosity calculator")
    parser.add_argument("--config", default="config_bdt_nom.yaml")
    parser.add_argument("--lumi-file",
                        default="/sphenix/user/shuhangli/ppg12/lumi/60cmLumi_fromJoey.list")
    parser.add_argument("--results", default="results")
    parser.add_argument("--figdir", default="../plotting/figures")
    parser.add_argument("--skip-photons", action="store_true",
                        help="Skip photon counting pass (scalers only)")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    data_glob = config["input"]["data_file"]
    treename = config["input"]["tree"]

    os.makedirs(args.results, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    # Load official lumi
    official_lumi = load_official_lumi(args.lumi_file)
    print(f"Official lumi file: {len(official_lumi)} runs")

    # Pass 1: scalers
    scaler_results = pass_scalers(data_glob, treename)
    print(f"\nScaler pass complete: {len(scaler_results)} runs")

    # Compute luminosity
    lumi_table = compute_luminosity(scaler_results, official_lumi)

    # Count trigger 30 active runs
    n_trig_active = sum(1 for d in lumi_table.values() if d["trig_active"])
    n_0mrad = sum(1 for d in lumi_table.values() if d["is_0mrad"])
    print(f"Runs with trigger 30 active: {n_trig_active}")
    print(f"0 mrad runs: {n_0mrad}")

    # Pass 2: photon counting
    photon_counts, event_counts = {}, {}
    if not args.skip_photons:
        photon_counts, event_counts = pass_photons(data_glob, treename, config, lumi_table)
        print(f"Photon counting complete: {len(photon_counts)} runs with photons")

    # Save CSV
    csv_path = os.path.join(args.results, "lumi_photon_perrun.csv")
    save_csv(lumi_table, photon_counts, event_counts, csv_path)

    # Make plots
    make_plots(lumi_table, photon_counts, event_counts, args.figdir)

    # Save ROOT
    root_path = os.path.join(args.results, "lumi_photon_qa.root")
    save_root(lumi_table, photon_counts, event_counts, root_path)

    # Summary of trigger 30 coverage
    trig30_off = sum(1 for d in lumi_table.values()
                     if not d["trig_active"] and d["n_live_mbd"] > 0)
    print(f"\n=== IMPORTANT ===")
    print(f"Trigger 30 is DISABLED in {trig30_off} runs (1.5 mrad period).")
    print(f"QA plots cover only 0 mrad runs ({n_trig_active} runs).")
    print(f"For 1.5 mrad coverage, a different trigger (e.g., 36) is needed.")


if __name__ == "__main__":
    main()
