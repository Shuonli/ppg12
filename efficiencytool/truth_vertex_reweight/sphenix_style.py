"""
Minimal sPHENIX matplotlib style for the truth-vertex-reweight plots.

Matches the conventions of the main PPG12 ROOT plotting pipeline
(plotcommon.h / sPhenixStyle.C): Helvetica-like sans-serif font, inward
ticks on all four axes, bold axis labels, no error bar caps, and the
canonical two-line annotation

    sPHENIX Internal
    1.5 mrad   (or 0 mrad)

placed as a TLatex-style overlay in plot coordinates. The crossing
angle replaces the beam/luminosity line in these MC-centric plots
(no data lumi is relevant).
"""
from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt


def set_sphenix_rcparams() -> None:
    """Apply sPHENIX-like global matplotlib rcParams. Safe to call
    multiple times."""
    matplotlib.rcParams.update({
        # Font: sans-serif Helvetica-alike
        "font.family":        "sans-serif",
        "font.sans-serif":    ["DejaVu Sans", "Helvetica", "Arial", "sans-serif"],
        "mathtext.fontset":   "dejavusans",
        "font.size":          11,
        # Axis labels: a little larger than tick labels
        "axes.labelsize":     12,
        "axes.titlesize":     12,
        "axes.labelweight":   "normal",
        "axes.linewidth":     1.2,
        # Tick style — inward ticks, on all four sides, like sPhenixStyle
        "xtick.direction":    "in",
        "ytick.direction":    "in",
        "xtick.top":          True,
        "ytick.right":        True,
        "xtick.major.size":   6,
        "xtick.minor.size":   3,
        "ytick.major.size":   6,
        "ytick.minor.size":   3,
        "xtick.major.width":  1.2,
        "ytick.major.width":  1.2,
        "xtick.minor.width":  1.0,
        "ytick.minor.width":  1.0,
        "xtick.labelsize":    10,
        "ytick.labelsize":    10,
        # Legend
        "legend.frameon":     False,
        "legend.fontsize":    9,
        # Figure
        "figure.dpi":         100,
        "savefig.dpi":        150,
        "savefig.bbox":       "tight",
    })


def sphenix_label(
    ax,
    period: str | None = None,
    loc: str = "upper left",
    x: float | None = None,
    y: float | None = None,
    size: float = 10.0,
) -> None:
    """Draw the sPHENIX Internal label on an axes in axes coordinates.

    If ``period`` is "1p5mrad" or "0mrad", a second line with the
    crossing angle ("1.5 mrad" / "0 mrad") is drawn below. Otherwise
    only "sPHENIX Internal" is drawn (used for plots that overlay both
    periods).
    """
    if x is None or y is None:
        positions = {
            "upper left":  (0.035, 0.96),
            "upper right": (0.96, 0.96),
            "lower left":  (0.035, 0.10),
            "lower right": (0.96, 0.10),
        }
        if loc not in positions:
            loc = "upper left"
        x, y = positions[loc]

    ha = "left" if "left" in loc else "right"
    va = "top" if "upper" in loc else "bottom"
    dy = -0.055 if va == "top" else +0.055

    ax.text(
        x, y, r"$\mathbf{\mathit{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha=ha, va=va, fontsize=size, color="k",
    )
    angle_text = None
    if period == "1p5mrad":
        angle_text = "1.5 mrad"
    elif period == "0mrad":
        angle_text = "0 mrad"
    if angle_text is not None:
        ax.text(
            x, y + dy, angle_text,
            transform=ax.transAxes,
            ha=ha, va=va, fontsize=size - 1, color="k",
        )


def sphenix_period_text(period: str) -> str:
    """Human-readable period tag, e.g. '1.5 mrad' / '0 mrad'."""
    if period == "1p5mrad":
        return "1.5 mrad"
    if period == "0mrad":
        return "0 mrad"
    return period
