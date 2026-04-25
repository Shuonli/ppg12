#!/usr/bin/env python3
"""
compare_single_double_showershapes.py
Quantitative comparison of shower shape distributions between single and
double interaction events from both analysis pipelines:
  1. Toy simulation  (DoubleInteractionCheck.C)  — vertex-shift method
  2. Full GEANT blending (run_showershape_double.sh) — weighted hadd of
     photon10_nom / photon10_double MC at physics fractions

Outputs:
  - Summary tables (mean, RMS, KS test, chi2/ndf) as CSV
  - Per-variable mean-shift vs pT plots
  - Overall summary figure: fractional mean shift heatmap
  - LaTeX-ready table fragments
"""
import argparse
import os
import sys
import numpy as np
import uproot

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SS_VARS = [
    "weta_cogx", "wphi_cogx", "wr_cogx", "et1", "et2", "et3", "et4",
    "e11_to_e33", "e32_to_e35", "bdt_score", "npb_score",
]

INTERACTION_TYPES = ["single", "double", "double_smear5", "double_smear10",
                     "double_smear15", "double_smear20"]

NICE_NAMES = {
    "weta_cogx":   r"$w_{\eta}^{\mathrm{CoG}}$",
    "wphi_cogx":   r"$w_{\phi}^{\mathrm{CoG}}$",
    "wr_cogx":     r"$w_{r}^{\mathrm{CoG}}$",
    "et1":         r"$f_{1}$",
    "et2":         r"$f_{2}$",
    "et3":         r"$f_{3}$",
    "et4":         r"$f_{4}$",
    "e11_to_e33":  r"$E_{1\times1}/E_{3\times3}$",
    "e32_to_e35":  r"$E_{3\times2}/E_{3\times5}$",
    "bdt_score":   r"BDT score",
    "npb_score":   r"NPB score",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def hist_stats(h):
    """Return (mean, rms, integral) from a uproot TH1."""
    vals = h.values()
    edges = h.axis().edges()
    centers = 0.5 * (edges[:-1] + edges[1:])
    total = vals.sum()
    if total <= 0:
        return np.nan, np.nan, 0.0
    mean = np.sum(centers * vals) / total
    rms = np.sqrt(np.sum((centers - mean)**2 * vals) / total)
    return mean, rms, float(total)


def ks_test(h1, h2):
    """KS test statistic between two uproot TH1 histograms (binned)."""
    v1 = h1.values().astype(float)
    v2 = h2.values().astype(float)
    s1, s2 = v1.sum(), v2.sum()
    if s1 <= 0 or s2 <= 0:
        return np.nan
    c1 = np.cumsum(v1) / s1
    c2 = np.cumsum(v2) / s2
    return float(np.max(np.abs(c1 - c2)))


def chi2_ndf(h1, h2):
    """Chi2/ndf between two normalized histograms."""
    v1 = h1.values().astype(float)
    v2 = h2.values().astype(float)
    s1, s2 = v1.sum(), v2.sum()
    if s1 <= 0 or s2 <= 0:
        return np.nan, 0
    n1 = v1 / s1
    n2 = v2 / s2
    # Poisson-like errors on normalized histograms
    e1 = np.sqrt(v1) / s1
    e2 = np.sqrt(v2) / s2
    e2_tot = e1**2 + e2**2
    mask = e2_tot > 0
    if mask.sum() == 0:
        return np.nan, 0
    chi2 = np.sum((n1[mask] - n2[mask])**2 / e2_tot[mask])
    ndf = int(mask.sum()) - 1
    return float(chi2), max(ndf, 1)


def get_hist_safe(f, name):
    """Return histogram or None if missing."""
    try:
        return f[name]
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Analysis 1: Toy simulation (DoubleInteractionCheck.C output)
# ---------------------------------------------------------------------------
def analyze_toy_simulation(results_dir, outdir, pt_bins):
    """Compare single vs double shower shapes from toy vertex-shift simulation."""
    sig_file = os.path.join(results_dir, "MC_efficiency_double_interaction_check_signal.root")
    jet_file = os.path.join(results_dir, "MC_efficiency_double_interaction_check_jet.root")

    n_pt = len(pt_bins) - 1
    rows = []

    for label, fpath in [("Signal MC", sig_file), ("Inclusive MC", jet_file)]:
        if not os.path.exists(fpath):
            print(f"[WARN] Missing {fpath}, skipping {label}")
            continue
        f = uproot.open(fpath)

        for var in SS_VARS:
            for ipt in range(n_pt):
                h_single = get_hist_safe(f, f"h_ss_single_{var}_eta0_pt{ipt}")
                h_double = get_hist_safe(f, f"h_ss_double_{var}_eta0_pt{ipt}")
                if h_single is None or h_double is None:
                    continue

                m_s, r_s, n_s = hist_stats(h_single)
                m_d, r_d, n_d = hist_stats(h_double)
                ks = ks_test(h_single, h_double)
                c2, ndf = chi2_ndf(h_single, h_double)

                frac_shift = (m_d - m_s) / m_s if abs(m_s) > 1e-12 else np.nan

                rows.append({
                    "sample": label,
                    "variable": var,
                    "pt_bin": ipt,
                    "pt_lo": pt_bins[ipt],
                    "pt_hi": pt_bins[ipt + 1],
                    "mean_single": m_s,
                    "mean_double": m_d,
                    "mean_shift": m_d - m_s,
                    "frac_shift": frac_shift,
                    "rms_single": r_s,
                    "rms_double": r_d,
                    "rms_shift": r_d - r_s,
                    "ks_stat": ks,
                    "chi2_ndf": c2 / ndf if ndf > 0 else np.nan,
                    "entries_single": n_s,
                    "entries_double": n_d,
                })

        # Also analyze smeared variants
        for smear in [5, 10, 15, 20]:
            for var in SS_VARS:
                for ipt in range(n_pt):
                    h_single = get_hist_safe(f, f"h_ss_single_{var}_eta0_pt{ipt}")
                    h_smear = get_hist_safe(f, f"h_ss_double_smear{smear}_{var}_eta0_pt{ipt}")
                    if h_single is None or h_smear is None:
                        continue

                    m_s, r_s, n_s = hist_stats(h_single)
                    m_sm, r_sm, n_sm = hist_stats(h_smear)
                    ks = ks_test(h_single, h_smear)
                    c2, ndf = chi2_ndf(h_single, h_smear)
                    frac_shift = (m_sm - m_s) / m_s if abs(m_s) > 1e-12 else np.nan

                    rows.append({
                        "sample": f"{label} smear{smear}cm",
                        "variable": var,
                        "pt_bin": ipt,
                        "pt_lo": pt_bins[ipt],
                        "pt_hi": pt_bins[ipt + 1],
                        "mean_single": m_s,
                        "mean_double": m_sm,
                        "mean_shift": m_sm - m_s,
                        "frac_shift": frac_shift,
                        "rms_single": r_s,
                        "rms_double": r_sm,
                        "rms_shift": r_sm - r_s,
                        "ks_stat": ks,
                        "chi2_ndf": c2 / ndf if ndf > 0 else np.nan,
                        "entries_single": n_s,
                        "entries_double": n_sm,
                    })

    return rows


# ---------------------------------------------------------------------------
# Analysis 2: Full GEANT blending (ShowerShapeCheck.C combined output)
# ---------------------------------------------------------------------------
def analyze_geant_blending(results_dir, outdir, pt_bins, config_suffix):
    """Compare nominal (single-only) vs blended (combined) MC shower shapes."""
    # File paths following the naming convention
    nom_sig = os.path.join(results_dir, f"MC_efficiencyshower_shape_signal_{config_suffix}.root")
    nom_jet = os.path.join(results_dir, f"MC_efficiencyshower_shape_jet_inclusive_{config_suffix}.root")
    com_sig = os.path.join(results_dir, f"MC_efficiencyshower_shape_signal_combined_{config_suffix}.root")
    com_jet = os.path.join(results_dir, f"MC_efficiencyshower_shape_jet_inclusive_combined_{config_suffix}.root")
    data_f  = os.path.join(results_dir, f"data_histoshower_shape_{config_suffix}.root")

    # h2d variables (ShowerShapeCheck naming)
    h2d_vars = [
        "h2d_weta_cogx", "h2d_wphi_cogx", "h2d_et1", "h2d_et2", "h2d_et3", "h2d_et4",
        "h2d_e11_to_e33", "h2d_e32_to_e35", "h2d_bdt", "h2d_npb_score",
    ]
    var_map = {
        "h2d_weta_cogx": "weta_cogx", "h2d_wphi_cogx": "wphi_cogx",
        "h2d_et1": "et1", "h2d_et2": "et2", "h2d_et3": "et3", "h2d_et4": "et4",
        "h2d_e11_to_e33": "e11_to_e33", "h2d_e32_to_e35": "e32_to_e35",
        "h2d_bdt": "bdt_score", "h2d_npb_score": "npb_score",
    }

    n_pt = len(pt_bins) - 1
    rows = []

    for label, nom_path, com_path in [
        ("Signal MC", nom_sig, com_sig),
        ("Inclusive MC", nom_jet, com_jet),
    ]:
        if not os.path.exists(nom_path) or not os.path.exists(com_path):
            print(f"[WARN] Missing {nom_path} or {com_path}, skipping {label}")
            continue
        f_nom = uproot.open(nom_path)
        f_com = uproot.open(com_path)

        for h2d_name in h2d_vars:
            var = var_map[h2d_name]
            for ipt in range(n_pt):
                for icut, cut_label in [(0, "cut0"), (1, "cut1")]:
                    hname = f"{h2d_name}_eta0_pt{ipt}_cut{icut}"
                    h2_nom = get_hist_safe(f_nom, hname)
                    h2_com = get_hist_safe(f_com, hname)
                    if h2_nom is None or h2_com is None:
                        continue

                    # Project X (shower shape axis) from TH2F
                    v_nom = h2_nom.values().sum(axis=1)  # sum over Y bins
                    v_com = h2_com.values().sum(axis=1)

                    edges = h2_nom.axis(0).edges()
                    centers = 0.5 * (edges[:-1] + edges[1:])

                    s_nom, s_com = v_nom.sum(), v_com.sum()
                    if s_nom <= 0 or s_com <= 0:
                        continue

                    m_nom = np.sum(centers * v_nom) / s_nom
                    m_com = np.sum(centers * v_com) / s_com
                    r_nom = np.sqrt(np.sum((centers - m_nom)**2 * v_nom) / s_nom)
                    r_com = np.sqrt(np.sum((centers - m_com)**2 * v_com) / s_com)

                    # KS (binned)
                    c_nom = np.cumsum(v_nom) / s_nom
                    c_com = np.cumsum(v_com) / s_com
                    ks = float(np.max(np.abs(c_nom - c_com)))

                    # Chi2/ndf
                    n_nom = v_nom / s_nom
                    n_com = v_com / s_com
                    e2 = (np.sqrt(v_nom) / s_nom)**2 + (np.sqrt(v_com) / s_com)**2
                    mask = e2 > 0
                    ndf = int(mask.sum()) - 1
                    c2 = float(np.sum((n_nom[mask] - n_com[mask])**2 / e2[mask])) if ndf > 0 else np.nan

                    frac_shift = (m_com - m_nom) / m_nom if abs(m_nom) > 1e-12 else np.nan

                    rows.append({
                        "sample": label,
                        "variable": var,
                        "cut": cut_label,
                        "pt_bin": ipt,
                        "pt_lo": pt_bins[ipt],
                        "pt_hi": pt_bins[ipt + 1],
                        "mean_nominal": m_nom,
                        "mean_combined": m_com,
                        "mean_shift": m_com - m_nom,
                        "frac_shift": frac_shift,
                        "rms_nominal": r_nom,
                        "rms_combined": r_com,
                        "rms_shift": r_com - r_nom,
                        "ks_stat": ks,
                        "chi2_ndf": c2 / ndf if ndf > 0 else np.nan,
                        "entries_nominal": s_nom,
                        "entries_combined": s_com,
                    })

    # Also compare data against both nominal and combined MC (chi2/ndf)
    if os.path.exists(data_f):
        f_data = uproot.open(data_f)
        for mc_label, mc_path, mc_tag in [
            ("Data vs Nominal MC", nom_jet, "nominal"),
            ("Data vs Combined MC", com_jet, "combined"),
        ]:
            if not os.path.exists(mc_path):
                continue
            f_mc = uproot.open(mc_path)
            for h2d_name in h2d_vars:
                var = var_map[h2d_name]
                for ipt in range(n_pt):
                    for icut, cut_label in [(1, "cut1")]:
                        hname = f"{h2d_name}_eta0_pt{ipt}_cut{icut}"
                        h2_data = get_hist_safe(f_data, hname)
                        h2_mc = get_hist_safe(f_mc, hname)
                        if h2_data is None or h2_mc is None:
                            continue

                        v_data = h2_data.values().sum(axis=1)
                        v_mc = h2_mc.values().sum(axis=1)

                        s_data, s_mc = v_data.sum(), v_mc.sum()
                        if s_data <= 0 or s_mc <= 0:
                            continue

                        n_data = v_data / s_data
                        n_mc = v_mc / s_mc
                        e2 = (np.sqrt(v_data) / s_data)**2 + (np.sqrt(v_mc) / s_mc)**2
                        mask = e2 > 0
                        ndf = int(mask.sum()) - 1
                        c2 = float(np.sum((n_data[mask] - n_mc[mask])**2 / e2[mask])) if ndf > 0 else np.nan

                        c_data = np.cumsum(v_data) / s_data
                        c_mc = np.cumsum(v_mc) / s_mc
                        ks = float(np.max(np.abs(c_data - c_mc)))

                        rows.append({
                            "sample": mc_label,
                            "variable": var,
                            "cut": cut_label,
                            "pt_bin": ipt,
                            "pt_lo": pt_bins[ipt],
                            "pt_hi": pt_bins[ipt + 1],
                            "mean_nominal": np.nan,
                            "mean_combined": np.nan,
                            "mean_shift": np.nan,
                            "frac_shift": np.nan,
                            "rms_nominal": np.nan,
                            "rms_combined": np.nan,
                            "rms_shift": np.nan,
                            "ks_stat": ks,
                            "chi2_ndf": c2 / ndf if ndf > 0 else np.nan,
                            "entries_nominal": s_data,
                            "entries_combined": s_mc,
                        })

    return rows


# ---------------------------------------------------------------------------
# Output: summary tables and LaTeX fragments
# ---------------------------------------------------------------------------
def write_csv(rows, path):
    import csv
    if not rows:
        return
    with open(path, "w", newline="") as fout:
        w = csv.DictWriter(fout, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"  Wrote {path} ({len(rows)} rows)")


def write_latex_summary(toy_rows, geant_rows, outdir, pt_bins):
    """Write a LaTeX table fragment summarizing the key findings."""
    # Aggregate toy simulation: mean fractional shift per variable (averaged over pT)
    tex_path = os.path.join(outdir, "di_shower_shape_summary.tex")
    with open(tex_path, "w") as f:
        f.write("% Auto-generated by compare_single_double_showershapes.py\n")
        f.write("% Toy simulation: mean fractional shift (double - single) / single, signal MC\n")
        f.write("\\begin{tabular}{l" + "c" * (len(pt_bins) - 1) + "}\n")
        f.write("\\toprule\n")
        pt_header = " & ".join([f"${pt_bins[i]}$--${pt_bins[i+1]}$ GeV" for i in range(len(pt_bins) - 1)])
        f.write(f"Variable & {pt_header} \\\\\n")
        f.write("\\midrule\n")

        for var in SS_VARS:
            nice = NICE_NAMES.get(var, var)
            vals = []
            for ipt in range(len(pt_bins) - 1):
                matches = [r for r in toy_rows if r["sample"] == "Signal MC"
                           and r["variable"] == var and r["pt_bin"] == ipt]
                if matches and not np.isnan(matches[0]["frac_shift"]):
                    vals.append(f"{matches[0]['frac_shift']*100:.1f}\\%")
                else:
                    vals.append("---")
            f.write(f"{nice} & " + " & ".join(vals) + " \\\\\n")

        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")

        # GEANT blending summary (cut1 only, inclusive MC)
        f.write("\n\n% Full GEANT blending: chi2/ndf nominal vs combined at cut1, inclusive MC\n")
        f.write("\\begin{tabular}{l" + "c" * (len(pt_bins) - 1) + "}\n")
        f.write("\\toprule\n")
        f.write(f"Variable & {pt_header} \\\\\n")
        f.write("\\midrule\n")

        geant_inc = [r for r in geant_rows if r["sample"] == "Inclusive MC"
                     and r.get("cut") == "cut1"]
        for var in ["weta_cogx", "wphi_cogx", "et1", "et4", "e11_to_e33",
                     "e32_to_e35", "bdt_score", "npb_score"]:
            nice = NICE_NAMES.get(var, var)
            vals = []
            for ipt in range(len(pt_bins) - 1):
                matches = [r for r in geant_inc if r["variable"] == var and r["pt_bin"] == ipt]
                if matches and not np.isnan(matches[0]["chi2_ndf"]):
                    vals.append(f"{matches[0]['chi2_ndf']:.1f}")
                else:
                    vals.append("---")
            f.write(f"{nice} & " + " & ".join(vals) + " \\\\\n")

        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")

    print(f"  Wrote {tex_path}")


def make_summary_plots(toy_rows, geant_rows, outdir, pt_bins):
    """Generate matplotlib summary plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] matplotlib not available, skipping plots")
        return

    n_pt = len(pt_bins) - 1
    pt_centers = [0.5 * (pt_bins[i] + pt_bins[i+1]) for i in range(n_pt)]

    # --- Plot 1: Fractional mean shift vs pT (toy sim, signal MC) ---
    fig, axes = plt.subplots(3, 4, figsize=(16, 10), sharex=True)
    axes = axes.flatten()
    for iv, var in enumerate(SS_VARS):
        ax = axes[iv]
        # Signal MC: single vs double
        shifts = []
        for ipt in range(n_pt):
            matches = [r for r in toy_rows if r["sample"] == "Signal MC"
                       and r["variable"] == var and r["pt_bin"] == ipt]
            shifts.append(matches[0]["frac_shift"] * 100 if matches and not np.isnan(matches[0]["frac_shift"]) else np.nan)
        ax.plot(pt_centers, shifts, "ko-", label="Signal: double", ms=5)

        # Jet MC: single vs double
        shifts_j = []
        for ipt in range(n_pt):
            matches = [r for r in toy_rows if r["sample"] == "Inclusive MC"
                       and r["variable"] == var and r["pt_bin"] == ipt]
            shifts_j.append(matches[0]["frac_shift"] * 100 if matches and not np.isnan(matches[0]["frac_shift"]) else np.nan)
        ax.plot(pt_centers, shifts_j, "rs-", label="Incl. MC: double", ms=5)

        # Smeared variants (signal only, 5 and 20 cm)
        for smear, marker, color in [(5, "^", "blue"), (20, "v", "green")]:
            shifts_sm = []
            for ipt in range(n_pt):
                matches = [r for r in toy_rows if r["sample"] == f"Signal MC smear{smear}cm"
                           and r["variable"] == var and r["pt_bin"] == ipt]
                shifts_sm.append(matches[0]["frac_shift"] * 100 if matches and not np.isnan(matches[0]["frac_shift"]) else np.nan)
            ax.plot(pt_centers, shifts_sm, marker=marker, color=color, ls="--", label=f"Sig: +{smear}cm smear", ms=4)

        ax.axhline(0, color="gray", ls="--", lw=0.5)
        ax.set_title(NICE_NAMES.get(var, var), fontsize=10)
        ax.set_ylabel("Mean shift [%]", fontsize=8)
        if iv >= 8:
            ax.set_xlabel(r"$E_T$ [GeV]", fontsize=9)
        if iv == 0:
            ax.legend(fontsize=6, loc="best")

    # Hide unused subplot
    if len(SS_VARS) < len(axes):
        axes[-1].set_visible(False)

    fig.suptitle("Toy Simulation: Fractional Mean Shift (double vs single)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out1 = os.path.join(outdir, "toy_sim_mean_shift_vs_pt.pdf")
    fig.savefig(out1)
    plt.close(fig)
    print(f"  Wrote {out1}")

    # --- Plot 2: KS statistic heatmap (toy sim, signal MC) ---
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ks_matrix = np.full((len(SS_VARS), n_pt), np.nan)
    for iv, var in enumerate(SS_VARS):
        for ipt in range(n_pt):
            matches = [r for r in toy_rows if r["sample"] == "Signal MC"
                       and r["variable"] == var and r["pt_bin"] == ipt]
            if matches and not np.isnan(matches[0]["ks_stat"]):
                ks_matrix[iv, ipt] = matches[0]["ks_stat"]

    im = ax2.imshow(ks_matrix, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.15)
    ax2.set_xticks(range(n_pt))
    ax2.set_xticklabels([f"{pt_bins[i]}-{pt_bins[i+1]}" for i in range(n_pt)], fontsize=8)
    ax2.set_yticks(range(len(SS_VARS)))
    ax2.set_yticklabels([NICE_NAMES.get(v, v) for v in SS_VARS], fontsize=8)
    ax2.set_xlabel(r"$E_T$ bin [GeV]")
    ax2.set_title("Toy Simulation: KS Statistic (single vs double), Signal MC")
    plt.colorbar(im, ax=ax2, label="KS statistic")
    # Annotate cells
    for iv in range(len(SS_VARS)):
        for ipt in range(n_pt):
            val = ks_matrix[iv, ipt]
            if not np.isnan(val):
                ax2.text(ipt, iv, f"{val:.3f}", ha="center", va="center", fontsize=7,
                         color="white" if val > 0.08 else "black")
    fig2.tight_layout()
    out2 = os.path.join(outdir, "toy_sim_ks_heatmap.pdf")
    fig2.savefig(out2)
    plt.close(fig2)
    print(f"  Wrote {out2}")

    # --- Plot 3: GEANT blending chi2/ndf comparison (nom vs combined vs data) ---
    geant_cut1 = [r for r in geant_rows if r.get("cut") == "cut1"]
    if geant_cut1:
        fig3, ax3 = plt.subplots(figsize=(8, 5))
        key_vars = ["weta_cogx", "wphi_cogx", "et1", "et4", "e11_to_e33",
                     "e32_to_e35", "bdt_score", "npb_score"]
        chi2_nom_vs_com = np.full((len(key_vars), n_pt), np.nan)
        for iv, var in enumerate(key_vars):
            for ipt in range(n_pt):
                matches = [r for r in geant_cut1 if r["sample"] == "Inclusive MC"
                           and r["variable"] == var and r["pt_bin"] == ipt]
                if matches and not np.isnan(matches[0]["chi2_ndf"]):
                    chi2_nom_vs_com[iv, ipt] = matches[0]["chi2_ndf"]

        im3 = ax3.imshow(chi2_nom_vs_com, aspect="auto", cmap="Blues", vmin=0,
                          vmax=min(np.nanmax(chi2_nom_vs_com) * 1.1, 20) if np.any(np.isfinite(chi2_nom_vs_com)) else 5)
        ax3.set_xticks(range(n_pt))
        ax3.set_xticklabels([f"{pt_bins[i]}-{pt_bins[i+1]}" for i in range(n_pt)], fontsize=8)
        ax3.set_yticks(range(len(key_vars)))
        ax3.set_yticklabels([NICE_NAMES.get(v, v) for v in key_vars], fontsize=8)
        ax3.set_xlabel(r"$E_T$ bin [GeV]")
        ax3.set_title(r"GEANT Blending: $\chi^2$/ndf (nominal vs combined), Incl. MC, cut1")
        plt.colorbar(im3, ax=ax3, label=r"$\chi^2$/ndf")
        for iv in range(len(key_vars)):
            for ipt in range(n_pt):
                val = chi2_nom_vs_com[iv, ipt]
                if not np.isnan(val):
                    ax3.text(ipt, iv, f"{val:.1f}", ha="center", va="center", fontsize=7,
                             color="white" if val > 5 else "black")
        fig3.tight_layout()
        out3 = os.path.join(outdir, "geant_blending_chi2_heatmap.pdf")
        fig3.savefig(out3)
        plt.close(fig3)
        print(f"  Wrote {out3}")

    # --- Plot 4: Data-MC chi2 improvement: nominal vs combined ---
    data_nom = [r for r in geant_rows if r["sample"] == "Data vs Nominal MC"]
    data_com = [r for r in geant_rows if r["sample"] == "Data vs Combined MC"]
    if data_nom and data_com:
        fig4, ax4 = plt.subplots(figsize=(8, 5))
        key_vars = ["weta_cogx", "wphi_cogx", "et1", "et4", "e11_to_e33",
                     "e32_to_e35", "bdt_score", "npb_score"]
        improvement = np.full((len(key_vars), n_pt), np.nan)
        for iv, var in enumerate(key_vars):
            for ipt in range(n_pt):
                r_nom = [r for r in data_nom if r["variable"] == var and r["pt_bin"] == ipt]
                r_com = [r for r in data_com if r["variable"] == var and r["pt_bin"] == ipt]
                if r_nom and r_com:
                    c_nom = r_nom[0]["chi2_ndf"]
                    c_com = r_com[0]["chi2_ndf"]
                    if not np.isnan(c_nom) and not np.isnan(c_com) and c_nom > 0:
                        improvement[iv, ipt] = (c_nom - c_com) / c_nom * 100

        im4 = ax4.imshow(improvement, aspect="auto", cmap="RdYlGn", vmin=-20, vmax=20)
        ax4.set_xticks(range(n_pt))
        ax4.set_xticklabels([f"{pt_bins[i]}-{pt_bins[i+1]}" for i in range(n_pt)], fontsize=8)
        ax4.set_yticks(range(len(key_vars)))
        ax4.set_yticklabels([NICE_NAMES.get(v, v) for v in key_vars], fontsize=8)
        ax4.set_xlabel(r"$E_T$ bin [GeV]")
        ax4.set_title(r"Data-MC $\chi^2$/ndf improvement: nominal $\to$ combined [%]")
        plt.colorbar(im4, ax=ax4, label="Improvement [%]")
        for iv in range(len(key_vars)):
            for ipt in range(n_pt):
                val = improvement[iv, ipt]
                if not np.isnan(val):
                    ax4.text(ipt, iv, f"{val:.1f}", ha="center", va="center", fontsize=7)
        fig4.tight_layout()
        out4 = os.path.join(outdir, "data_mc_chi2_improvement.pdf")
        fig4.savefig(out4)
        plt.close(fig4)
        print(f"  Wrote {out4}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Compare single vs double interaction shower shapes")
    parser.add_argument("--results", default="/sphenix/user/shuhangli/ppg12/efficiencytool/results",
                        help="Results directory")
    parser.add_argument("--outdir", default="/sphenix/user/shuhangli/ppg12/plotting/figures/di_comparison",
                        help="Output directory for plots and tables")
    parser.add_argument("--config-suffix", default="showershape",
                        help="Config suffix for GEANT blending files "
                             "(default 'showershape' = all-range bare; "
                             "per-angle feeders: 'showershape_0rad', 'showershape_1p5mrad')")
    parser.add_argument("--pt-bins", nargs="+", type=float, default=[10, 14, 18, 22, 28, 30],
                        help="pT bin edges")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    pt_bins = args.pt_bins

    print("=== Analysis 1: Toy simulation (DoubleInteractionCheck.C) ===")
    toy_rows = analyze_toy_simulation(args.results, args.outdir, pt_bins)
    write_csv(toy_rows, os.path.join(args.outdir, "toy_simulation_comparison.csv"))

    print(f"\n=== Analysis 2: Full GEANT blending ({args.config_suffix}) ===")
    geant_rows = analyze_geant_blending(args.results, args.outdir, pt_bins, args.config_suffix)
    write_csv(geant_rows, os.path.join(args.outdir, f"geant_blending_comparison_{args.config_suffix}.csv"))

    # Also run per-angle feeders if they exist alongside the requested suffix.
    # For the bare all-range suffix "showershape", check both _0rad and _1p5mrad.
    # For a per-angle suffix, check the other period.
    geant_rows_1p5 = []
    companion_suffixes = []
    base = args.config_suffix
    if base.endswith("_0rad"):
        companion_suffixes.append(base[:-len("_0rad")] + "_1p5mrad")
    elif base.endswith("_1p5mrad"):
        companion_suffixes.append(base[:-len("_1p5mrad")] + "_0rad")
    elif base.endswith("_1p5rad"):  # legacy
        companion_suffixes.append(base[:-len("_1p5rad")] + "_0rad")
    else:
        # bare (e.g. "showershape") — show both per-angle feeders if present
        companion_suffixes.extend([base + "_0rad", base + "_1p5mrad"])

    for suffix_alt in companion_suffixes:
        if suffix_alt == args.config_suffix:
            continue
        com_check = os.path.join(args.results, f"MC_efficiencyshower_shape_signal_combined_{suffix_alt}.root")
        if os.path.exists(com_check):
            print(f"\n=== Analysis 2b: Full GEANT blending ({suffix_alt}) ===")
            geant_rows_alt = analyze_geant_blending(args.results, args.outdir, pt_bins, suffix_alt)
            write_csv(geant_rows_alt, os.path.join(args.outdir, f"geant_blending_comparison_{suffix_alt}.csv"))
            geant_rows_1p5 = geant_rows_alt  # keep name for downstream summary function

    print("\n=== Writing summary tables and plots ===")
    write_latex_summary(toy_rows, geant_rows, args.outdir, pt_bins)
    make_summary_plots(toy_rows, geant_rows, args.outdir, pt_bins)

    # Print key findings to stdout
    print("\n" + "=" * 70)
    print("KEY FINDINGS")
    print("=" * 70)

    # Largest shifts in toy simulation
    toy_signal = [r for r in toy_rows if r["sample"] == "Signal MC"]
    if toy_signal:
        toy_signal.sort(key=lambda r: abs(r.get("frac_shift", 0) or 0), reverse=True)
        print("\nToy simulation — largest fractional mean shifts (signal MC, double vs single):")
        for r in toy_signal[:10]:
            print(f"  {r['variable']:15s}  pT {r['pt_lo']}-{r['pt_hi']} GeV:  "
                  f"shift = {r['frac_shift']*100:+.2f}%  KS = {r['ks_stat']:.4f}")

    # GEANT blending — largest chi2/ndf between nominal and combined
    geant_mc = [r for r in geant_rows if r["sample"] == "Inclusive MC" and r.get("cut") == "cut1"]
    if geant_mc:
        geant_mc.sort(key=lambda r: abs(r.get("chi2_ndf", 0) or 0), reverse=True)
        print(f"\nGEANT blending ({args.config_suffix}) — largest chi2/ndf (nominal vs combined), cut1, incl. MC:")
        for r in geant_mc[:10]:
            print(f"  {r['variable']:15s}  pT {r['pt_lo']}-{r['pt_hi']} GeV:  "
                  f"chi2/ndf = {r['chi2_ndf']:.2f}  frac_shift = {r['frac_shift']*100:+.3f}%")

    print(f"\nAll outputs in: {args.outdir}/")


if __name__ == "__main__":
    main()
