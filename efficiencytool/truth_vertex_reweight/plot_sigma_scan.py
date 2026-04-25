#!/usr/bin/env python3
"""
Chi2-vs-sigma_tgt scan plot for the truth vertex reweight report.

Shows how the parametric Gaussian reweight chi2 varies as a function of
the target width sigma_tgt, and where the iterative method lands
relative to the best parametric result. One panel per crossing-angle
period (1.5 mrad, 0 mrad).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sphenix_style import set_sphenix_rcparams, sphenix_label, sphenix_period_text

set_sphenix_rcparams()


def load_fit_result(period_dir: Path) -> dict:
    yaml_path = period_dir / "fit_result.yaml"
    with open(yaml_path) as f:
        return yaml.safe_load(f)


def main():
    here = Path(__file__).resolve().parent
    out_dir = here / "output"

    periods = [
        ("1p5mrad", "1.5 mrad"),
        ("0mrad", "0 mrad"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    for ax, (period_key, period_label) in zip(axes, periods):
        doc = load_fit_result(out_dir / period_key)
        s1 = doc["stage1_parametric"]
        it = doc["iterative"]

        grid_sigma = np.asarray(s1["grid_sigma"])
        grid_chi2 = np.asarray(s1["grid_chi2"])
        ndf_param = int(s1["ndf"])
        best_sigma = float(s1["sigma_tgt"])
        best_chi2 = float(s1["chi2"])

        iter_chi2 = float(it["chi2_in_window"])
        iter_ndf = int(it["ndf_in_window"])

        # Normalize to chi2/ndf
        grid_chi2_ndf = grid_chi2 / ndf_param
        best_chi2_ndf = best_chi2 / ndf_param
        iter_chi2_ndf = iter_chi2 / iter_ndf

        # Improvement factor: best parametric chi2/ndf over iterative chi2/ndf
        improvement = best_chi2_ndf / iter_chi2_ndf

        # --- Scan curve ---
        ax.plot(grid_sigma, grid_chi2_ndf, "o-", color="C0", markersize=3,
                lw=1.3, label=r"parametric scan", zorder=3)

        # --- Best-fit sigma vertical line ---
        ax.axvline(best_sigma, color="C0", ls="--", lw=1.2, alpha=0.8)

        # Place the best-sigma annotation in a fixed axes-fraction
        # position so it stays clear of the scan curve regardless of
        # where the minimum falls.  Use the center of the panel
        # horizontally, at 45% up from the bottom (well below the
        # scan curve, well above the iterative line).
        ax.annotate(
            fr"$\sigma_{{\mathrm{{tgt}}}}^{{\mathrm{{best}}}} = {best_sigma:.1f}$ cm",
            xy=(best_sigma, best_chi2_ndf),
            xytext=(0.50, 0.45), textcoords="axes fraction",
            fontsize=9.5, color="C0", ha="center",
            arrowprops=dict(arrowstyle="->", color="C0", alpha=0.6,
                            connectionstyle="arc3,rad=-0.15"),
        )

        # --- Iterative chi2/ndf horizontal line ---
        ax.axhline(iter_chi2_ndf, color="C3", ls="--", lw=1.4, alpha=0.85,
                    label=fr"iterative ($\chi^2/\mathrm{{ndf}} = {iter_chi2_ndf:.1f}$)",
                    zorder=2)

        # --- Reference chi2/ndf = 1 ---
        ax.axhline(1.0, color="gray", ls=":", lw=1.0, alpha=0.6,
                    label=r"$\chi^2/\mathrm{ndf} = 1$", zorder=1)

        # --- Improvement annotation ---
        # Place in the lower-center of the panel, away from the scan
        # curve and the sPHENIX label.
        ax.text(
            0.50, 0.15,
            f"parametric best / iterative\n"
            fr"$= {best_chi2_ndf:.1f} \,/\, {iter_chi2_ndf:.1f} = {improvement:.0f}\times$",
            transform=ax.transAxes,
            ha="center", va="bottom", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray",
                      alpha=0.85),
        )

        # --- Axes formatting ---
        ax.set_yscale("log")
        ax.set_xlabel(r"$\sigma_{\mathrm{tgt}}$ [cm]")
        ax.set_ylabel(r"$\chi^2 / \mathrm{ndf}$")
        ax.set_title(f"{period_label}", fontsize=12)
        ax.legend(fontsize=8.5, loc="upper left")
        ax.grid(alpha=0.25, which="both")

        sphenix_label(ax, period=period_key, loc="upper right", size=9.5)

    fig.suptitle(
        r"Stage 1 parametric $\chi^2$ scan over $\sigma_{\mathrm{tgt}}$"
        r" — $p$+$p$ truth vertex reweight",
        fontsize=11, y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    # Save to reports/figures/
    out_path = here.parents[1] / "reports" / "figures" / "vertex_reweight_sigma_scan.pdf"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    print(f"wrote {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
