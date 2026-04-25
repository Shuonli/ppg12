#!/usr/bin/env python3
"""Plot three-stage efficiency (reco, ID, iso) for converted vs unconverted photons.

Reads rootFiles/efficiencies_summary.json and produces:
  figures/efficiency_vs_pt.pdf           2x2 grid: eps_reco, eps_id, eps_iso, eps_all
  figures/efficiency_cond_vs_pt.pdf      eps_id|reco and eps_iso|reco
  figures/efficiency_ratio_vs_pt.pdf     conv/unconv ratio per stage
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
JSON_IN = HERE / "rootFiles" / "efficiencies_summary.json"
FIG_DIR = HERE / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

STYLE = {
    "unconv": dict(label="Unconverted", color="black", marker="o",
                   fillstyle="full", linestyle="-"),
    "conv": dict(label="Converted", color="red", marker="s",
                 fillstyle="none", linestyle="--"),
}


def sphenix_label(ax, x=0.03, y=0.97, extra=None):
    ax.text(x, y, r"$\bf{\it{sPHENIX}}$ Simulation", transform=ax.transAxes,
            ha="left", va="top", fontsize=10)
    ax.text(x, y - 0.06, r"photon10, $|\eta^{\gamma}|<0.7$",
            transform=ax.transAxes, ha="left", va="top", fontsize=8, color="gray")
    if extra:
        ax.text(x, y - 0.11, extra, transform=ax.transAxes, ha="left", va="top",
                fontsize=8, color="gray")


def errorbar_cat(ax, pt_centers, y, yerr, cat):
    s = STYLE[cat]
    mask = np.isfinite(y) & (np.asarray(y) >= 0)
    ax.errorbar(np.asarray(pt_centers)[mask], np.asarray(y)[mask],
                yerr=np.asarray(yerr)[mask],
                color=s["color"], marker=s["marker"],
                mfc=s["color"] if s["fillstyle"] == "full" else "white",
                linestyle=s["linestyle"], markersize=6, capsize=2,
                linewidth=1.3, label=s["label"])


def load():
    with open(JSON_IN) as f:
        return json.load(f)


def plot_absolute(summary):
    pt_edges = np.array(summary["pt_edges"])
    pt_centers = 0.5 * (pt_edges[:-1] + pt_edges[1:])
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    axes = axes.ravel()

    stages = [
        ("eps_reco", r"$\varepsilon_{\rm reco}$", "Reco efficiency"),
        ("eps_id", r"$\varepsilon_{\rm ID}$", "ID efficiency (reco $\\wedge$ BDT)"),
        ("eps_iso", r"$\varepsilon_{\rm iso}$", "Iso efficiency (reco $\\wedge$ iso)"),
        ("eps_all", r"$\varepsilon_{\rm total}$",
         "Total efficiency (reco $\\wedge$ BDT $\\wedge$ iso)"),
    ]

    for ax, (key, ylab, title) in zip(axes, stages):
        for cat in ("unconv", "conv"):
            y = summary["per_pt"][cat][key]
            yerr = summary["per_pt"][cat][f"{key}_err"]
            errorbar_cat(ax, pt_centers, y, yerr, cat)
        ax.set_xlabel(r"$p_T^{\rm truth}$ [GeV]")
        ax.set_ylabel(ylab)
        ax.set_ylim(0, 1.0)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9, loc="lower right")
        ax.set_title(title, fontsize=11)
        sphenix_label(ax)

    fig.suptitle("Selection efficiency vs truth $p_T$ (absolute, per truth photon)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = FIG_DIR / "efficiency_vs_pt.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


def plot_conditional(summary):
    pt_edges = np.array(summary["pt_edges"])
    pt_centers = 0.5 * (pt_edges[:-1] + pt_edges[1:])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    stages = [
        ("eps_id_cond", r"$\varepsilon_{\rm ID \mid reco}$",
         "ID efficiency given reco"),
        ("eps_iso_cond", r"$\varepsilon_{\rm iso \mid reco}$",
         "Iso efficiency given reco"),
    ]

    for ax, (key, ylab, title) in zip(axes, stages):
        for cat in ("unconv", "conv"):
            y = summary["per_pt"][cat][key]
            yerr = summary["per_pt"][cat][f"{key}_err"]
            errorbar_cat(ax, pt_centers, y, yerr, cat)
        ax.set_xlabel(r"$p_T^{\rm truth}$ [GeV]")
        ax.set_ylabel(ylab)
        ax.set_ylim(0, 1.0)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9, loc="lower right")
        ax.set_title(title, fontsize=11)
        sphenix_label(ax)

    fig.tight_layout()
    out = FIG_DIR / "efficiency_cond_vs_pt.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


def plot_ratio(summary):
    pt_edges = np.array(summary["pt_edges"])
    pt_centers = 0.5 * (pt_edges[:-1] + pt_edges[1:])
    fig, ax = plt.subplots(figsize=(7.0, 5.0))

    stages = [
        ("eps_reco", r"reco", "tab:blue", "o", "-"),
        ("eps_id", r"ID",    "tab:green", "s", "--"),
        ("eps_iso", r"iso",  "tab:orange", "^", "-."),
        ("eps_all", r"total", "k", "D", ":"),
    ]

    for key, label, color, marker, ls in stages:
        u = np.asarray(summary["per_pt"]["unconv"][key], dtype=float)
        c = np.asarray(summary["per_pt"]["conv"][key], dtype=float)
        u_err = np.asarray(summary["per_pt"]["unconv"][f"{key}_err"], dtype=float)
        c_err = np.asarray(summary["per_pt"]["conv"][f"{key}_err"], dtype=float)
        mask = (u > 0) & np.isfinite(u) & np.isfinite(c)
        ratio = np.where(mask, c / np.where(u > 0, u, 1), np.nan)
        # Error propagation for independent binomials (approximate)
        with np.errstate(invalid="ignore", divide="ignore"):
            rel_err = np.sqrt((u_err / np.where(u > 0, u, 1)) ** 2
                               + (c_err / np.where(c > 0, c, 1)) ** 2)
        err = ratio * rel_err
        ax.errorbar(pt_centers[mask], ratio[mask], yerr=err[mask],
                    color=color, marker=marker, linestyle=ls,
                    markersize=6, capsize=2, linewidth=1.3,
                    label=label, mfc=color)

    ax.axhline(1.0, color="gray", linestyle=":", lw=0.8)
    ax.set_xlabel(r"$p_T^{\rm truth}$ [GeV]")
    ax.set_ylabel(r"$\varepsilon_{\rm conv} / \varepsilon_{\rm unconv}$")
    ax.set_ylim(0.0, 1.5)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc="lower right", ncol=2)
    sphenix_label(ax)
    ax.set_title("Converted / unconverted efficiency ratio per stage", fontsize=11)
    fig.tight_layout()
    out = FIG_DIR / "efficiency_ratio_vs_pt.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    summary = load()
    plot_absolute(summary)
    plot_conditional(summary)
    plot_ratio(summary)


if __name__ == "__main__":
    main()
