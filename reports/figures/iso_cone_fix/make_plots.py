#!/usr/bin/env python3
"""Diagnostic plots for the isoET cone-centering fix.

Compares BEFORE (cluster-COG axis, no threshold, -cluster_ET scalar subtraction
path) vs AFTER (CoG-tower axis, 120 MeV threshold, manual path) isolation
computations in CaloAna24.cc.

Inputs are ROOT files produced by the slimtree maker test pipeline at
anatreemaker/macro_maketree/tests_iso_radii/{data, sim_photon10}/.
"""

from __future__ import annotations
import sys, os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BEFORE = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/tests_iso_radii/data/OUTTREE_before_fix.root"
AFTER  = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/tests_iso_radii/data/OUTTREE_DST_Jet_run2pp_ana521_2025p007_v001-00047289-00000.root"
OUTDIR = os.path.dirname(os.path.abspath(__file__))

C = "CLUSTERINFO_CEMC"


def flatten(arr_of_arr):
    return np.concatenate([np.asarray(x, dtype=float) for x in arr_of_arr])


def load(path, extra=None):
    t = uproot.open(path)["slimtree"]
    br = [f"cluster_Et_{C}", f"cluster_Eta_{C}",
          f"cluster_iso_005_{C}", f"cluster_iso_0075_{C}",
          f"cluster_iso_02_{C}", f"cluster_iso_03_{C}", f"cluster_iso_04_{C}",
          f"cluster_iso_topo_005_{C}", "vertexz"]
    if extra:
        br += extra
    a = t.arrays(br, library="np")
    ET  = flatten(a[f"cluster_Et_{C}"])
    eta = flatten(a[f"cluster_Eta_{C}"])
    i005 = flatten(a[f"cluster_iso_005_{C}"])
    i0075 = flatten(a[f"cluster_iso_0075_{C}"])
    i02 = flatten(a[f"cluster_iso_02_{C}"])
    i03 = flatten(a[f"cluster_iso_03_{C}"])
    i04 = flatten(a[f"cluster_iso_04_{C}"])
    it005 = flatten(a[f"cluster_iso_topo_005_{C}"])
    vz = np.concatenate([np.full(len(x), v) for v, x in zip(a["vertexz"], a[f"cluster_Et_{C}"])])
    out = {
        "ET": ET, "eta": eta, "vz": vz,
        "i005": i005, "i0075": i0075, "i02": i02, "i03": i03, "i04": i04,
        "it005": it005,
    }
    if extra:
        for b in extra:
            short = b.replace(f"_{C}", "")
            out[short] = flatten(a[b])
    return out


def mask_et(d, et_min=1.0):
    m = d["ET"] > et_min
    return {k: v[m] for k, v in d.items()}


def legend_counts(b, a):
    return f"before ($n={len(b['ET'])}$)", f"after ($n={len(a['ET'])}$)"


def plot_1d_iso_over_et(b, a):
    fig, ax = plt.subplots(figsize=(6, 4.2))
    bins = np.linspace(-1.2, 1.0, 56)
    lb, la = legend_counts(b, a)
    ax.hist(b["i005"] / b["ET"], bins=bins, histtype="stepfilled", color="#d62728",
            alpha=0.45, label=lb, edgecolor="#d62728", linewidth=1.4)
    ax.hist(a["i005"] / a["ET"], bins=bins, histtype="step", color="#1f77b4",
            linewidth=2.0, label=la)
    ax.axvline(-1, color="gray", linestyle=":", alpha=0.5)
    ax.annotate("unphysical\n(cone misses cluster)",
                xy=(-0.97, ax.get_ylim()[1]*0.7), xytext=(-0.75, ax.get_ylim()[1]*0.85),
                fontsize=8, ha="left", color="#d62728",
                arrowprops=dict(arrowstyle="->", color="#d62728", lw=0.8))
    ax.set_xlabel(r"cluster_iso_005 / cluster_$E_T$")
    ax.set_ylabel("clusters / bin")
    ax.set_title(r"$R=0.05$ iso ratio: before vs after fix (data, $E_T>1$ GeV)")
    ax.legend(loc="upper right", frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "iso005_ratio_1d.pdf"))
    plt.close(fig)


def plot_scatter_vs_vz(b, a):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, d, title, col in zip(axes, (b, a), ("BEFORE fix", "AFTER fix"),
                                 ("#d62728", "#1f77b4")):
        ax.scatter(np.abs(d["vz"]), d["i005"] / d["ET"], s=14, alpha=0.7,
                   color=col, edgecolor="none")
        ax.axhline(0, color="black", linewidth=0.5)
        ax.axhline(-1, color="gray", linestyle=":", alpha=0.5)
        ax.set_xlim(0, 150)
        ax.set_ylim(-1.2, 1.0)
        ax.set_xlabel(r"$|v_z|$ [cm]")
        ax.set_title(f"{title} ($n={len(d['ET'])}$)")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r"cluster_iso_005 / cluster_$E_T$")
    fig.suptitle(r"$R=0.05$ iso ratio vs $|v_z|$ (data, $E_T>1$ GeV)", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "iso005_vs_vz.pdf"))
    plt.close(fig)


def plot_pathological_by_vz(b, a):
    edges = np.array([0, 20, 40, 60, 80, 100, 150])
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = edges[1:] - edges[:-1]

    def frac(d):
        rat = d["i005"] / d["ET"]
        out = np.zeros(len(edges) - 1)
        for i in range(len(edges) - 1):
            m = (np.abs(d["vz"]) >= edges[i]) & (np.abs(d["vz"]) < edges[i + 1])
            out[i] = np.mean(((rat > -1.1) & (rat < -0.9))[m]) if m.sum() else 0.0
        return out

    fb = frac(b) * 100
    fa = frac(a) * 100
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    w = 0.4
    ax.bar(centers - w * 5, fb, width=w * 10, alpha=0.6, color="#d62728", label="before fix")
    ax.bar(centers + w * 5, fa, width=w * 10, alpha=0.6, color="#1f77b4", label="after fix")
    ax.set_xlabel(r"$|v_z|$ [cm]")
    ax.set_ylabel(r"cluster_iso_005 $\approx -E_T$ (pathological) fraction [%]")
    ax.set_title("Pathological-cluster rate vs $|v_z|$ (data)")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25, axis="y")
    ax.set_ylim(0, 110)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "pathological_vs_vz.pdf"))
    plt.close(fig)


def plot_iso_cascade(b, a):
    """Cumulative iso(R) for individual clusters; before-fix violates monotonicity
    badly at small R, after-fix is monotone."""
    Rs = np.array([0.05, 0.075, 0.2, 0.3, 0.4])

    def stack(d):
        return np.stack([d["i005"] / d["ET"], d["i0075"] / d["ET"],
                         d["i02"] / d["ET"], d["i03"] / d["ET"],
                         d["i04"] / d["ET"]], axis=1)

    sb = stack(b)
    sa = stack(a)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, s, title, col in zip(axes, (sb, sa), ("BEFORE fix", "AFTER fix"),
                                 ("#d62728", "#1f77b4")):
        for row in s:
            ax.plot(Rs, row, color=col, alpha=0.35, lw=0.8)
        ax.plot(Rs, np.median(s, axis=0), color="black", lw=2.2, label="median")
        ax.axhline(0, color="gray", lw=0.5)
        ax.set_xlabel(r"cone radius $R$")
        ax.set_title(f"{title} ($n={len(s)}$)")
        ax.set_xticks(Rs)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r"cluster_iso_R / cluster_$E_T$")
    axes[0].legend(frameon=False, loc="lower right")
    axes[1].legend(frameon=False, loc="lower right")
    fig.suptitle("Per-cluster iso ratio vs cone radius (monotonicity)", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "iso_cascade.pdf"))
    plt.close(fig)


def plot_excl(a):
    fig, ax = plt.subplots(figsize=(6, 4.2))
    bins = np.linspace(-0.2, 2.0, 45)
    ax.hist(a["cluster_iso_excl_005"], bins=bins, histtype="stepfilled", color="#2ca02c",
            alpha=0.45, label=r"$R=0.05$", edgecolor="#2ca02c", lw=1.3)
    ax.hist(a["cluster_iso_excl_02"], bins=bins, histtype="step", color="#9467bd",
            lw=1.8, label=r"$R=0.2$")
    ax.hist(a["cluster_iso_excl_04"], bins=bins, histtype="step", color="#ff7f0e",
            lw=1.8, label=r"$R=0.4$")
    ax.set_xlabel(r"cluster_iso_excl_R [GeV] (ownership-excluded, no $-E_T$)")
    ax.set_ylabel("clusters / bin")
    ax.set_title(r"New iso_excl distributions, data ($E_T>1$ GeV, threshold 120 MeV)")
    ax.axvline(0, color="black", lw=0.5)
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "iso_excl_distributions.pdf"))
    plt.close(fig)


def main():
    b = mask_et(load(BEFORE))
    a = mask_et(load(AFTER, extra=[f"cluster_iso_excl_005_{C}",
                                   f"cluster_iso_excl_02_{C}",
                                   f"cluster_iso_excl_04_{C}"]))
    # Key numbers — used in the LaTeX report.
    def pathol(d, tag):
        rat = d["i005"] / d["ET"]
        return np.mean((rat > -1.1) & (rat < -0.9))
    print(f"BEFORE: n={len(b['ET'])}, pathological iso_005 frac = {pathol(b,'before')*100:.1f}%")
    print(f"AFTER:  n={len(a['ET'])}, pathological iso_005 frac = {pathol(a,'after' )*100:.1f}%")
    for lo, hi in [(0, 40), (40, 80), (80, 200)]:
        mb = (np.abs(b["vz"]) >= lo) & (np.abs(b["vz"]) < hi)
        ma = (np.abs(a["vz"]) >= lo) & (np.abs(a["vz"]) < hi)
        pb = np.mean(((b["i005"][mb] / b["ET"][mb] > -1.1) &
                      (b["i005"][mb] / b["ET"][mb] < -0.9))) if mb.sum() else 0
        pa = np.mean(((a["i005"][ma] / a["ET"][ma] > -1.1) &
                      (a["i005"][ma] / a["ET"][ma] < -0.9))) if ma.sum() else 0
        print(f"  |vz| {lo:3d}-{hi:3d}cm: before n={mb.sum():3d}, pathol={pb*100:.1f}%   "
              f"after n={ma.sum():3d}, pathol={pa*100:.1f}%")

    plot_1d_iso_over_et(b, a)
    plot_scatter_vs_vz(b, a)
    plot_pathological_by_vz(b, a)
    plot_iso_cascade(b, a)
    plot_excl(a)
    print(f"figures written to {OUTDIR}")


if __name__ == "__main__":
    main()
