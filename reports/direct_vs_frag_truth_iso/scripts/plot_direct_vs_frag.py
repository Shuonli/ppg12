#!/usr/bin/env python3
"""
Plot direct vs fragmentation truth-iso distributions and f_direct(iso_cut).

Reads ``direct_vs_frag_truth.root`` produced by
``extract_direct_vs_frag.py`` and writes PDFs under ``../figures/``.

Figure 1: iso_ET_truth stacked (direct + frag) per pT bin, 3 panels, log-y.
Figure 2: f_direct vs iso_cut, 3 curves (one per pT bin), linear-y.
"""

from __future__ import annotations

import argparse
import math
import os

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

# Must match extract_direct_vs_frag.py
PT_BINS = [(8.0, 14.0), (14.0, 22.0), (22.0, 36.0)]
ISO_CUTS = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, math.inf]

# sPHENIX-ish palette: black headline, blue direct, red frag.
COLOR_DIRECT = "#1f3d9e"
COLOR_FRAG = "#c62d2d"

LEGEND_SUBTITLE = (r"$\bf{\mathit{sPHENIX}}$ Internal"
                   "\n"
                   r"$p$+$p$  $\sqrt{s}=200$ GeV, Pythia truth")


def read_hists(path):
    """Return dict[(cls, pt_idx)] -> (values, edges)."""
    f = uproot.open(path)
    out = {}
    for key in f.keys():
        if not key.startswith("h_iso_"):
            continue
        # name = h_iso_{direct|frag}_pt{lo}_{hi}
        # Skip per-sample histograms (they have a trailing _photon*)
        # e.g. h_iso_direct_pt8_14_photon5
        base = key.split(";")[0]
        parts = base.split("_")
        if len(parts) != 5:
            # Per-sample hists have 6 tokens; skip
            continue
        cls_name = parts[2]
        if cls_name not in ("direct", "frag"):
            continue
        pt_lo = int(parts[3][2:])  # "pt8" -> 8
        pt_hi = int(parts[4])
        # Find pT bin index
        pt_idx = None
        for k, (lo, hi) in enumerate(PT_BINS):
            if int(lo) == pt_lo and int(hi) == pt_hi:
                pt_idx = k
                break
        if pt_idx is None:
            continue
        hist = f[key]
        vals, edges = hist.to_numpy()
        out[(cls_name, pt_idx)] = (vals, edges)
    return out


def read_scan(path):
    """Return dict[(pt_idx, iso_cut)] -> (N_direct, N_frag)."""
    f = uproot.open(path)
    t = f["scan_tree"]
    arr = t.arrays(library="np")
    out = {}
    for i in range(len(arr["pt_bin_low"])):
        lo = float(arr["pt_bin_low"][i])
        hi = float(arr["pt_bin_high"][i])
        iso = float(arr["iso_cut"][i])
        pt_idx = None
        for k, (plo, phi) in enumerate(PT_BINS):
            if abs(plo - lo) < 1e-3 and abs(phi - hi) < 1e-3:
                pt_idx = k
                break
        if pt_idx is None:
            continue
        key = (pt_idx, iso)
        out[key] = (float(arr["N_direct"][i]), float(arr["N_frag"][i]))
    return out


def fig1_stacked_iso(hists, outpath):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5), sharey=True)
    for k, (lo, hi) in enumerate(PT_BINS):
        ax = axes[k]
        d_vals, edges = hists.get(("direct", k), (None, None))
        f_vals, _ = hists.get(("frag", k), (None, None))
        if d_vals is None or f_vals is None:
            ax.text(0.5, 0.5, "missing", transform=ax.transAxes, ha="center")
            continue
        centers = 0.5 * (edges[:-1] + edges[1:])
        widths = edges[1:] - edges[:-1]
        # Stack: frag below, direct on top
        ax.bar(centers, f_vals, width=widths, color=COLOR_FRAG, alpha=0.85,
               label="Fragmentation", edgecolor="none")
        ax.bar(centers, d_vals, width=widths, bottom=f_vals,
               color=COLOR_DIRECT, alpha=0.85, label="Direct",
               edgecolor="none")
        ax.set_yscale("log")
        ax.set_xlabel(r"$E_{T}^{\mathrm{iso,truth}}$ (R=0.3) [GeV]")
        if k == 0:
            ax.set_ylabel("Weighted counts / 0.25 GeV")
        ax.set_xlim(0, 15)
        # Draw vertical line at the nominal iso cut = 4 GeV
        ax.axvline(4.0, color="k", linestyle="--", lw=1.0, alpha=0.6)
        ax.text(4.15, ax.get_ylim()[1] * 0.2 if False else 0.95,
                r"PPG12 iso$<$4", transform=ax.get_xaxis_transform(),
                rotation=90, fontsize=8, ha="left", va="top", color="k",
                alpha=0.7)
        ax.set_title(
            rf"$p_T^{{\gamma,\mathrm{{truth}}}} \in [{int(lo)},{int(hi)})$ GeV"
        )
        if k == 0:
            ax.legend(loc="upper right", frameon=False, fontsize=9,
                      title=LEGEND_SUBTITLE, title_fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close(fig)
    print(f"[plot] wrote {outpath}")


def fig2_fdirect_vs_iso(scan, outpath):
    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    colors = ["#000000", "#1f3d9e", "#c62d2d"]
    markers = ["o", "s", "^"]
    x_for_inf = 14.0  # display "inf" at x=14 GeV

    for k, (lo, hi) in enumerate(PT_BINS):
        xs, ys = [], []
        for iso_cut in ISO_CUTS:
            iso_val = 1e9 if math.isinf(iso_cut) else float(iso_cut)
            entry = scan.get((k, iso_val))
            if entry is None:
                continue
            nd, nf = entry
            tot = nd + nf
            if tot <= 0:
                continue
            fd = nd / tot
            x = x_for_inf if math.isinf(iso_cut) else iso_cut
            xs.append(x); ys.append(fd)
        ax.plot(xs, ys,
                color=colors[k], marker=markers[k], lw=1.4, ms=7,
                label=rf"$p_T^{{\gamma,\mathrm{{truth}}}} \in [{int(lo)},{int(hi)})$ GeV")

    ax.set_xlabel(r"$E_{T}^{\mathrm{iso,truth}}$ cut (R=0.3) [GeV]")
    ax.set_ylabel(r"$f_{\mathrm{direct}} = N_{\gamma,\mathrm{direct}} / "
                  r"(N_{\mathrm{direct}} + N_{\mathrm{frag}})$")
    ax.set_xlim(0, 15)
    ax.set_ylim(0.0, 1.02)
    # Axis tick labels: replace 14 with "inf"
    xticks = [1, 2, 3, 4, 5, 6, 8, 10, x_for_inf]
    xtick_lbls = ["1", "2", "3", "4", "5", "6", "8", "10", r"$\infty$"]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_lbls)
    ax.axvline(4.0, color="gray", linestyle="--", lw=1.0, alpha=0.6)
    ax.text(4.15, 0.05, "PPG12 iso<4", rotation=90, color="gray",
            fontsize=9, ha="left", va="bottom", alpha=0.8)
    ax.axhline(1.0, color="k", linestyle=":", lw=0.8, alpha=0.3)
    ax.legend(loc="lower left", frameon=False, fontsize=10,
              title=LEGEND_SUBTITLE, title_fontsize=9)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close(fig)
    print(f"[plot] wrote {outpath}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",
                   default=os.path.join(os.path.dirname(__file__), "..",
                                        "direct_vs_frag_truth.root"))
    p.add_argument("--figdir",
                   default=os.path.join(os.path.dirname(__file__), "..",
                                        "figures"))
    args = p.parse_args()

    figdir = os.path.abspath(args.figdir)
    os.makedirs(figdir, exist_ok=True)

    inp = os.path.abspath(args.input)
    if not os.path.exists(inp):
        raise FileNotFoundError(f"Input ROOT not found: {inp}")

    hists = read_hists(inp)
    scan = read_scan(inp)

    fig1_stacked_iso(hists, os.path.join(figdir, "iso_ET_truth_stacked.pdf"))
    fig2_fdirect_vs_iso(scan, os.path.join(figdir, "f_direct_vs_iso_cut.pdf"))


if __name__ == "__main__":
    main()
